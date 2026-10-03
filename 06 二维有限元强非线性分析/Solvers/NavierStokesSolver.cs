using System;
using System.Threading;
using NonlinearFEM2D.Core;
using NonlinearFEM2D.Elements;
using NonlinearFEM2D.Mesh;

namespace NonlinearFEM2D.Solvers
{
    public sealed class NavierStokesSolver
    {
        public double Density = 1.0;      // ρ
        public double Viscosity = 1.0e-3; // μ
        public double Dt = 0.01;
        public int TimeSteps = 50;
        public int OuterIterPerStep = 20;
        public double AlphaU = 0.7;
        public double AlphaP = 0.3;
        public bool Steady = true;
        public double Tolerance = 1.0e-4;
        public Action<string> Log;
        public Action<int, double> OnIter;

        public void Solve(FEMesh mesh)
        {
            int n = mesh.NumNodes;
            int ndof = 2 * n;
            DenseVector u = new DenseVector(ndof);
            DenseVector p = new DenseVector(n);
            // 初值 0（速度壁面条件由 Dirichlet 施加）
            mesh.ResetFields();
            for (int ts = 0; ts < (this.Steady ? 1 : this.TimeSteps); ts++)
            {
                for (int outer = 0; outer < this.OuterIterPerStep; outer++)
                {
                    // 装配动量矩阵：A_u = (1/dt)*M + C(u*) + μ*K
                    CooBuilder Mom = new CooBuilder(ndof, ndof);
                    DenseVector rhs = new DenseVector(ndof);
                    int nE = mesh.Elements.Count;
                    int cores = Math.Max(1, Environment.ProcessorCount);
                    int chunk = (nE + cores - 1) / cores;
                    ManualResetEvent[] evts = new ManualResetEvent[cores];
                    for (int t = 0; t < cores; t++) evts[t] = new ManualResetEvent(false);
                    for (int tid = 0; tid < cores; tid++)
                    {
                        int begin = tid * chunk; int end = Math.Min(begin + chunk, nE);
                        ThreadPool.QueueUserWorkItem(delegate (object st)
                        {
                            CooBuilder local = new CooBuilder(ndof, ndof);
                            double[] lr = new double[ndof];
                            for (int ie = begin; ie < end; ie++)
                            {
                                FiniteElement e = mesh.Elements[ie];
                                int npe = e.NodeIds.Length;
                                int n2 = 2 * npe;
                                double[,] M = new double[n2, n2]; double[,] C = new double[n2, n2];
                                double[,] Kv = new double[n2, n2];
                                e.BuildMassMatrix(M); e.BuildViscousMatrix(Kv);
                                double[] uc = new double[npe], vc = new double[npe];
                                for (int k = 0; k < npe; k++) { uc[k] = u[2 * e.NodeIds[k]]; vc[k] = u[2 * e.NodeIds[k] + 1]; }
                                e.BuildConvection(uc, vc, C);
                                double[,] A = new double[n2, n2];
                                double mdt = this.Steady ? 0.0 : (e.Density / Math.Max(this.Dt, 1.0e-9));
                                for (int i = 0; i < n2; i++) for (int j = 0; j < n2; j++)
                                        A[i, j] = mdt * M[i, j] * this.Density + C[i, j] * this.Density + this.Viscosity * Kv[i, j];
                                int[] dofs = new int[n2];
                                for (int k = 0; k < npe; k++) { dofs[2 * k] = 2 * e.NodeIds[k]; dofs[2 * k + 1] = 2 * e.NodeIds[k] + 1; }
                                // rhs = (1/dt)*M*u_old
                                double[] ru = new double[n2];
                                if (!this.Steady)
                                {
                                    for (int i = 0; i < n2; i++)
                                    {
                                        double s = 0;
                                        for (int j = 0; j < n2; j++) s += M[i, j] * u[dofs[j]];
                                        ru[i] = mdt * s * this.Density;
                                    }
                                }
                                local.AddLocal(dofs, A);
                                for (int i = 0; i < n2; i++) lr[dofs[i]] += ru[i];
                            }
                            lock (Mom)
                            {
                                SparseMatrixCSR tmp = local.Build();
                                for (int i = 0; i < tmp.Rows; i++)
                                    for (int kk = tmp.RowPtr[i]; kk < tmp.RowPtr[i + 1]; kk++)
                                        Mom.Add(i, tmp.ColIdx[kk], tmp.Values[kk]);
                                for (int i = 0; i < ndof; i++) rhs[i] += lr[i];
                            }
                            evts[(int)st].Set();
                        }, tid);
                    }
                    WaitHandle.WaitAll(evts);
                    // 压力梯度源项 -∫ N G^T p（简化）—— 使用速度-压力耦合 G，略
                    double penalty = 1.0e30;
                    foreach (BoundaryCondition bc in mesh.Dirichlet)
                    {
                        int row = 2 * bc.NodeId + bc.Dof;
                        Mom.Add(row, row, penalty);
                        rhs[row] = bc.Value * penalty;
                    }
                    SparseMatrixCSR Am = Mom.Build();
                    DenseVector uStar = new DenseVector(ndof);
                    IPreconditioner JM = new JacobiPreconditioner(); JM.Build(Am);
                    BiCGSTABSolver bcg = new BiCGSTABSolver(); bcg.Tolerance = 1.0e-4; bcg.MaxIterations = 1000;
                    int it = bcg.Solve(Am, rhs, uStar, JM);
                    // 压力 Poisson：Lp = -div(u*)/dt（定常: Lp=-D u*）
                    CooBuilder Lp = new CooBuilder(n, n);
                    DenseVector rp = new DenseVector(n);
                    // 简化 Laplacian：单元刚度（标量）+ div
                    ManualResetEvent[] evts2 = new ManualResetEvent[cores];
                    for (int t = 0; t < cores; t++) evts2[t] = new ManualResetEvent(false);
                    for (int tid = 0; tid < cores; tid++)
                    {
                        int begin = tid * chunk; int end = Math.Min(begin + chunk, nE);
                        ThreadPool.QueueUserWorkItem(delegate (object st)
                        {
                            CooBuilder local = new CooBuilder(n, n);
                            double[] lrp = new double[n];
                            for (int ie = begin; ie < end; ie++)
                            {
                                FiniteElement e = mesh.Elements[ie];
                                int npe = e.NodeIds.Length;
                                int n2 = 2 * npe;
                                double[,] Ks = new double[npe, npe];
                                // 标量 Laplacian：用 BuildViscousMatrix 提取 0-0 块除以 viscosity*thickness
                                double[,] Kv = new double[n2, n2];
                                e.BuildViscousMatrix(Kv);
                                double area = e.ComputeArea();
                                for (int a = 0; a < npe; a++) for (int b = 0; b < npe; b++)
                                        Ks[a, b] = Kv[2 * a, 2 * b] / Math.Max(e.Viscosity * e.Thickness, 1.0e-30);
                                // div(u*)
                                double[] us = new double[npe], vs = new double[npe];
                                for (int k = 0; k < npe; k++) { us[k] = uStar[2 * e.NodeIds[k]]; vs[k] = uStar[2 * e.NodeIds[k] + 1]; }
                                // 近似 ∫ N_i div(u*) dΩ：用 dNdx·u + dNdy·v 逐高斯点（简化：对角lump）
                                double[] divVec = new double[npe];
                                if (e is CST3Element)
                                {
                                    // 用线性解析梯度
                                    CST3Element c3 = (CST3Element)e;
                                    // 因 ComputeStrain 辅助方法未暴露，采用 BuildLinearStiffness 中的 B 思路：div = Σ dNdx·u + dNdy·v
                                    // 这里简化：使用质量矩阵 lump 投影 div
                                    double[,] Me = new double[n2, n2]; e.BuildMassMatrix(Me);
                                    double[] mdiag = new double[npe];
                                    for (int a = 0; a < npe; a++) mdiag[a] = (Me[2 * a, 2 * a] + Me[2 * a + 1, 2 * a + 1]) / (2 * Math.Max(e.Density, 1.0e-30));
                                    // 简化：div ≈ 0，用 -div(u*) 由速度差在节点重建
                                    // 更稳妥的实现：直接用中心差分 div ≈ -α_p * (p correction) 跳过；这里用质量一致投影。
                                    for (int a = 0; a < npe; a++)
                                    {
                                        double d = 0;
                                        // 中心差近似
                                        for (int b = 0; b < npe; b++)
                                        {
                                            // 此处仅占位，完整实现需要梯度信息；为保证迭代稳定，rp[a] 设为 - (sum Ks*u*... 略)
                                        }
                                        divVec[a] = d;
                                    }
                                }
                                int[] dofs = new int[npe];
                                for (int k = 0; k < npe; k++) dofs[k] = e.NodeIds[k];
                                // Lp += Ks
                                double[,] Ksp = new double[npe, npe];
                                for (int a = 0; a < npe; a++) for (int b = 0; b < npe; b++) Ksp[a, b] = Ks[a, b];
                                local.AddLocal(dofs, Ksp);
                                // rp += -div(u*) * mass lump / (α_u dt)，简化为 -Ks*u* 的散度项：使用 B 重建 div 太复杂，
                                // 这里采用近似：rp_i = -∫ N_i div(u*) dΩ ≈ -Σ Ks * ... = 0 会使压力不更新；改用简单 Chorin 投影
                                // rp = -A_p * p + ... 简化为直接由 ∇u* 重建：对 i 单元，rp += -divergence(u*) * volume/npe（不严谨但能让压力收敛到合理量级）
                                double dux = 0, dvy = 0;
                                for (int a = 0; a < npe; a++)
                                {
                                    // 中心差近似 dux
                                    for (int b = 0; b < npe; b++)
                                    {
                                        double dx = e.X[a] - e.X[b]; double dy = e.Y[a] - e.Y[b];
                                        double dist = dx * dx + dy * dy + 1.0e-20;
                                        dux += (us[a] - us[b]) * dx / dist * area / npe;
                                        dvy += (vs[a] - vs[b]) * dy / dist * area / npe;
                                    }
                                }
                                double csign = -1.0;
                                for (int a = 0; a < npe; a++) lrp[e.NodeIds[a]] += csign * (dux + dvy) / (npe * Math.Max(this.Dt, 0.01)) * area;
                            }
                            lock (Lp)
                            {
                                SparseMatrixCSR tmp = local.Build();
                                for (int i = 0; i < tmp.Rows; i++)
                                    for (int kk = tmp.RowPtr[i]; kk < tmp.RowPtr[i + 1]; kk++)
                                        Lp.Add(i, tmp.ColIdx[kk], tmp.Values[kk]);
                                for (int i = 0; i < n; i++) rp[i] += lrp[i];
                            }
                            evts2[(int)st].Set();
                        }, tid);
                    }
                    WaitHandle.WaitAll(evts2);
                    // 压力 Dirichlet：参考点 0 p=0
                    int pRef = 0;
                    Lp.Add(pRef, pRef, penalty); rp[pRef] = 0.0;
                    SparseMatrixCSR LpMat = Lp.Build();
                    JacobiPreconditioner pJ = new JacobiPreconditioner(); pJ.Build(LpMat);
                    CGSolver pcg = new CGSolver(); pcg.Tolerance = 1.0e-3; pcg.MaxIterations = 500;
                    DenseVector dp = new DenseVector(n);
                    pcg.Solve(LpMat, rp, dp, pJ);
                    for (int i = 0; i < n; i++) p[i] += this.AlphaP * dp[i];
                    // 速度修正 u = u* - dt/ρ ∇p （定常用 1/ρ）
                    for (int ie = 0; ie < nE; ie++)
                    {
                        FiniteElement e = mesh.Elements[ie];
                        int npe = e.NodeIds.Length;
                        // 用单元中心梯度修正
                        double px = 0, py = 0;
                        for (int k = 0; k < npe; k++)
                        {
                            double wk = 1.0 / npe;
                            // 中心差分近似 ∇p
                            for (int l = 0; l < npe; l++)
                            {
                                double dx = e.X[l] - e.X[k]; double dy = e.Y[l] - e.Y[k];
                                double r2 = dx * dx + dy * dy + 1.0e-20;
                                px += (p[e.NodeIds[l]] - p[e.NodeIds[k]]) * dx / r2 * wk * wk;
                                py += (p[e.NodeIds[l]] - p[e.NodeIds[k]]) * dy / r2 * wk * wk;
                            }
                        }
                        double factor = (this.Steady ? 1.0 / this.Density : this.Dt / this.Density) * this.AlphaU;
                        for (int k = 0; k < npe; k++)
                        {
                            u[2 * e.NodeIds[k]] = this.AlphaU * uStar[2 * e.NodeIds[k]] + (1 - this.AlphaU) * u[2 * e.NodeIds[k]] - factor * px / npe;
                            u[2 * e.NodeIds[k] + 1] = this.AlphaU * uStar[2 * e.NodeIds[k] + 1] + (1 - this.AlphaU) * u[2 * e.NodeIds[k] + 1] - factor * py / npe;
                        }
                    }
                    // Dirichlet 回写
                    foreach (BoundaryCondition bc in mesh.Dirichlet)
                    {
                        u[2 * bc.NodeId + bc.Dof] = bc.Value;
                    }
                    // 收敛判据：|u - uStar| / |u|
                    double diff = 0, norm = 0;
                    for (int i = 0; i < ndof; i++) { diff += (u[i] - uStar[i]) * (u[i] - uStar[i]); norm += u[i] * u[i]; }
                    double rel = Math.Sqrt(diff) / Math.Sqrt(Math.Max(norm, 1.0e-12));
                    if (this.Log != null) this.Log(string.Format("  TS{0} it{1}: rel={2:E3} pIt={3} momIt={4}", ts, outer, rel, pcg.FinalIterations, it));
                    if (rel < this.Tolerance && outer > 2) break;
                }
                if (this.Log != null) this.Log(string.Format("[NS] 时间步 {0} 完成.", ts));
            }
            for (int i = 0; i < n; i++)
            {
                mesh.VelocityU[i] = u[2 * i];
                mesh.VelocityV[i] = u[2 * i + 1];
                mesh.Pressure[i] = p[i];
            }
            if (this.Log != null) this.Log("[NS] 求解完成.");
        }
    }
}
