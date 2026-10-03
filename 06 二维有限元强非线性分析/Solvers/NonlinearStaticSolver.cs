using System;
using System.Threading;
using NonlinearFEM2D.Core;
using NonlinearFEM2D.Elements;
using NonlinearFEM2D.Mesh;

namespace NonlinearFEM2D.Solvers
{
    public enum LoadControlType { Force, Displacement }

    public sealed class NonlinearStaticSolver
    {
        public KrylovKind Krylov = KrylovKind.BiCGSTAB;
        public PrecondKind Precond = PrecondKind.AMG;
        public double Tolerance = 1.0e-6;
        public int MaxIterations = 50;
        public int LoadSteps = 30;
        public bool UseGeometricNonlinearity = false;
        public LoadControlType LoadControl = LoadControlType.Force;
        public double DisplacementTarget = 0.05;
        public int DisplacementNodeId = -1;
        public int DisplacementDof = 1;
        public Action<string> Log;
        public Action<int, double> OnIter;

        public void Solve(FEMesh mesh)
        {
            int ndof = 2 * mesh.NumNodes;
            mesh.ResetFields();
            // 重置每个单元的历史变量
            foreach (FiniteElement e in mesh.Elements) e.Initialize();
            DenseVector U = new DenseVector(ndof);
            DenseVector Fext = new DenseVector(ndof);
            for (int i = 0; i < mesh.NumNodes; i++) { Fext[2 * i] = mesh.NodeForceX[i]; Fext[2 * i + 1] = mesh.NodeForceY[i]; }
            double penalty = 1.0e30;
            int dispRow = -1; double dispTargetAbs = Math.Abs(this.DisplacementTarget);
            double totalLoad = 0;
            int nSteps = this.LoadSteps;
            for (int step = 1; step <= nSteps; step++)
            {
                double lambda = (double)step / nSteps;
                // 构造当前步外荷载 f
                DenseVector f = new DenseVector(ndof);
                if (this.LoadControl == LoadControlType.Force)
                {
                    for (int i = 0; i < ndof; i++) f[i] = lambda * Fext[i];
                    // Dirichlet
                    foreach (BoundaryCondition bc in mesh.Dirichlet)
                    {
                        int row = 2 * bc.NodeId + bc.Dof;
                        f[row] = bc.Value * penalty;
                    }
                }
                else
                {
                    // 位移控：目标节点某个方向逐渐加到 DisplacementTarget，其余节点力=0
                    for (int i = 0; i < ndof; i++) f[i] = 0;
                    foreach (BoundaryCondition bc in mesh.Dirichlet)
                    {
                        int row = 2 * bc.NodeId + bc.Dof;
                        f[row] = bc.Value * penalty;
                    }
                    if (this.DisplacementNodeId < 0) this.DisplacementNodeId = mesh.NumNodes - 1;
                    int drow = 2 * this.DisplacementNodeId + this.DisplacementDof;
                    dispRow = drow;
                    f[drow] = lambda * this.DisplacementTarget * penalty;
                }
                // N-R 迭代
                int it;
                double res0 = 0;
                for (it = 1; it <= this.MaxIterations; it++)
                {
                    // 装配 Kt & Rint
                    int nElem = mesh.Elements.Count;
                    CooBuilder KtCoo = new CooBuilder(ndof, ndof);
                    DenseVector Rint = new DenseVector(ndof);
                    int cores = Math.Max(1, Environment.ProcessorCount);
                    int chunk = (nElem + cores - 1) / cores;
                    ManualResetEvent[] evts = new ManualResetEvent[cores];
                    for (int t = 0; t < cores; t++) evts[t] = new ManualResetEvent(false);
                    for (int tid = 0; tid < cores; tid++)
                    {
                        int begin = tid * chunk; int end = Math.Min(begin + chunk, nElem);
                        ThreadPool.QueueUserWorkItem(delegate (object st)
                        {
                            CooBuilder local = new CooBuilder(ndof, ndof);
                            double[] localR = new double[ndof];
                            for (int ie = begin; ie < end; ie++)
                            {
                                FiniteElement e = mesh.Elements[ie];
                                int npe = e.NodeIds.Length;
                                double[] eu = new double[npe], ev = new double[npe];
                                for (int k = 0; k < npe; k++)
                                {
                                    eu[k] = U[2 * e.NodeIds[k]];
                                    ev[k] = U[2 * e.NodeIds[k] + 1];
                                }
                                double[,] Ke = new double[2 * npe, 2 * npe];
                                double[] fi = new double[2 * npe];
                                e.BuildTangentAndInternal(eu, ev, Ke, fi, this.UseGeometricNonlinearity);
                                int[] dofs = new int[2 * npe];
                                for (int k = 0; k < npe; k++) { dofs[2 * k] = 2 * e.NodeIds[k]; dofs[2 * k + 1] = 2 * e.NodeIds[k] + 1; }
                                local.AddLocal(dofs, Ke);
                                for (int k = 0; k < 2 * npe; k++) localR[dofs[k]] += fi[k];
                            }
                            lock (KtCoo)
                            {
                                SparseMatrixCSR tmp = local.Build();
                                for (int i = 0; i < tmp.Rows; i++)
                                    for (int kk = tmp.RowPtr[i]; kk < tmp.RowPtr[i + 1]; kk++)
                                        KtCoo.Add(i, tmp.ColIdx[kk], tmp.Values[kk]);
                                for (int i = 0; i < ndof; i++) Rint[i] += localR[i];
                            }
                            evts[(int)st].Set();
                        }, tid);
                    }
                    WaitHandle.WaitAll(evts);
                    // 残差 g = f - Rint；Dirichlet 行强置
                    DenseVector g = new DenseVector(ndof);
                    if (this.LoadControl == LoadControlType.Force)
                    {
                        foreach (BoundaryCondition bc in mesh.Dirichlet)
                        {
                            int row = 2 * bc.NodeId + bc.Dof;
                            KtCoo.Add(row, row, penalty);
                        }
                    }
                    else
                    {
                        foreach (BoundaryCondition bc in mesh.Dirichlet)
                        {
                            int row = 2 * bc.NodeId + bc.Dof;
                            KtCoo.Add(row, row, penalty);
                        }
                        KtCoo.Add(dispRow, dispRow, penalty);
                    }
                    SparseMatrixCSR Kt = KtCoo.Build();
                    for (int i = 0; i < ndof; i++) g[i] = f[i] - Rint[i];
                    foreach (BoundaryCondition bc in mesh.Dirichlet)
                    {
                        int row = 2 * bc.NodeId + bc.Dof;
                        g[row] = f[row] - U[row] * penalty; // 将该自由度推向目标位移
                    }
                    if (this.LoadControl == LoadControlType.Displacement)
                        g[dispRow] = f[dispRow] - U[dispRow] * penalty;
                    double res = g.Norm2();
                    if (it == 1) res0 = Math.Max(res, 1.0e-12);
                    double rel = res / res0;
                    if (this.Log != null) this.Log(string.Format("  Step {0}/{1} It {2}: |g|={3:E3} rel={4:E3}",
                        step, nSteps, it, res, rel));
                    if (rel < this.Tolerance || res < 1.0e-8)
                    {
                        // 统计屈服点
                        int yielded = 0; int totalGP = 0; double maxEqp = 0; double maxVM = 0;
                        foreach (FiniteElement e in mesh.Elements)
                        {
                            foreach (Elements.GPHistory gh in e.Histories)
                            {
                                totalGP++;
                                if (gh.Yielded) yielded++;
                                if (gh.EqPlasticStrain > maxEqp) maxEqp = gh.EqPlasticStrain;
                                if (gh.VonMises > maxVM) maxVM = gh.VonMises;
                            }
                        }
                        totalLoad = lambda;
                        if (this.Log != null) this.Log(string.Format("  -> Converged. Yielded GP: {0}/{1} maxEqp={2:E3} maxVM={3:E3}",
                            yielded, totalGP, maxEqp, maxVM));
                        break;
                    }
                    // 求解 Kt dU = g
                    IPreconditioner M = StructuralSolver.BuildPrecond(this.Precond, Kt);
                    ILinearSolver solver = StructuralSolver.BuildSolver(this.Krylov);
                    solver.Tolerance = this.Tolerance * 0.1;
                    solver.MaxIterations = 2000;
                    if (this.OnIter != null) solver.OnIteration = this.OnIter;
                    DenseVector dU = new DenseVector(ndof);
                    int it2 = solver.Solve(Kt, g, dU, M);
                    if (this.Log != null) this.Log(string.Format("    KSP {0}+{1} it={2} res={3:E3}", this.Krylov, this.Precond, it2, solver.FinalResidual));
                    U.Axpy(1.0, dU);
                    // 将当前 U 写回单元历史（应力）——我们已在装配 Kt 时就地更新了 histories（因为 BuildTangentAndInternal 直接写 h）
                    // 这里不用额外操作。
                }
            }
            // 写回位移场
            for (int i = 0; i < mesh.NumNodes; i++)
            {
                mesh.DisplacementU[i] = U[2 * i];
                mesh.DisplacementV[i] = U[2 * i + 1];
            }
            // 节点应力恢复
            int[] cnt = new int[mesh.NumNodes];
            Array.Clear(mesh.StressXX, 0, mesh.NumNodes);
            Array.Clear(mesh.StressYY, 0, mesh.NumNodes);
            Array.Clear(mesh.StressXY, 0, mesh.NumNodes);
            Array.Clear(mesh.VonMises, 0, mesh.NumNodes);
            Array.Clear(mesh.EqPlasticStrain, 0, mesh.NumNodes);
            foreach (FiniteElement e in mesh.Elements)
            {
                int npe = e.NodeIds.Length;
                double sx = 0, sy = 0, sxy = 0, vm = 0, eqp = 0; int ng = e.Histories.Length;
                foreach (Elements.GPHistory gh in e.Histories)
                {
                    sx += gh.Stress[0]; sy += gh.Stress[1]; sxy += gh.Stress[3]; vm += gh.VonMises; eqp += gh.EqPlasticStrain;
                }
                sx /= ng; sy /= ng; sxy /= ng; vm /= ng; eqp /= ng;
                for (int k = 0; k < npe; k++)
                {
                    int nid = e.NodeIds[k];
                    mesh.StressXX[nid] += sx; mesh.StressYY[nid] += sy; mesh.StressXY[nid] += sxy;
                    mesh.VonMises[nid] += vm; mesh.EqPlasticStrain[nid] += eqp;
                    cnt[nid]++;
                }
            }
            for (int i = 0; i < mesh.NumNodes; i++) if (cnt[i] > 0)
                {
                    mesh.StressXX[i] /= cnt[i]; mesh.StressYY[i] /= cnt[i]; mesh.StressXY[i] /= cnt[i];
                    mesh.VonMises[i] /= cnt[i]; mesh.EqPlasticStrain[i] /= cnt[i];
                }
            if (this.Log != null) this.Log(string.Format("[强非线性] 完成. 最终载荷系数={0:F3}", totalLoad));
        }
    }
}
