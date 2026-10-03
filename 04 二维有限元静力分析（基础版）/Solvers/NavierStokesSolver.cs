using System;
using System.Threading.Tasks;
using FEM2D.Core;
using FEM2D.Elements;
using FEM2D.Mesh;

namespace FEM2D.Solvers
{
    public sealed class NSResult
    {
        public DenseVector Velocity; // 3*N: 0=u,1=v,2=p
        public double[] U;
        public double[] V;
        public double[] P;
        public double[] Speed;
        public int Iterations;
        public double FinalResidual;
    }

    /// <summary>
    /// 不可压缩 Navier-Stokes SIMPLE 求解器
    /// 等阶插值（u,v,p 同单元），采用压力稳定（通过压力泊松方程近似，Mass lumped pressure-Laplacian + 欠松弛）。
    /// 非稳态采用向后欧拉：(u^{n+1} - u^n)/dt = -C(u^*)u^* - grad(p) + nu Laplacian(u)
    /// </summary>
    public sealed class NavierStokesSolver
    {
        public double Density { get; set; } = 1.0;
        public double Viscosity { get; set; } = 0.001;
        public double AlphaU { get; set; } = 0.7;
        public double AlphaP { get; set; } = 0.3;
        public int MaxOuterIter { get; set; } = 200;
        public double ConvergenceTol { get; set; } = 1e-5;
        public double TimeStep { get; set; } = 0.0; // 0=稳态
        public int TimeSteps { get; set; } = 1;
        public bool Unsteady => TimeStep > 1e-12;

        public NSResult Solve(FEMesh mesh, Action<string> log = null)
        {
            int nNodes = mesh.NumNodes;
            int nElem = mesh.NumElements;
            int ndof = 3 * nNodes; // u,v,p
            int vdof = 2 * nNodes;
            int pdof = nNodes;

            FiniteElement[] elems = new FiniteElement[nElem];
            for (int e = 0; e < nElem; e++) elems[e] = ElementFactory.Create(mesh.ElementTypes[e]);

            // 预先构造粘性扩散矩阵 + 质量矩阵（速度质量矩阵，集总）
            SparseMatrixCSR.CooBuilder viscBuilder = new SparseMatrixCSR.CooBuilder(vdof);
            SparseMatrixCSR.CooBuilder massBuilder = new SparseMatrixCSR.CooBuilder(vdof);
            SparseMatrixCSR.CooBuilder pMassBuilder = new SparseMatrixCSR.CooBuilder(pdof);
            SparseMatrixCSR.CooBuilder pLapBuilder = new SparseMatrixCSR.CooBuilder(pdof);
            double[] lumpMass = new double[vdof];
            double[] lumpPMass = new double[pdof];

            Parallel.For(0, nElem, e =>
            {
                FiniteElement fe = elems[e];
                int npe = mesh.Elements[e].Length;
                int edof = 2 * npe;
                double[,] visc = new double[edof, edof];
                double[,] me = new double[edof, edof];
                double[,] B = new double[edof, npe];

                fe.ComputeViscous(mesh, e, Viscosity, 1.0, visc);
                fe.ComputeMass(mesh, e, Density, 1.0, me);
                fe.ComputePressureCoupling(mesh, e, 1.0, B);

                // 压力质量矩阵（pMe 必须在 lock 块外声明，以便下面集总循环访问）
                double[,] pMe = new double[npe, npe];
                for (int i = 0; i < npe; i++)
                {
                    for (int j = 0; j < npe; j++)
                    {
                        pMe[i, j] = me[2 * i, 2 * j] / Density;
                    }
                }

                int[] conn = mesh.Elements[e];
                int[] vMap = new int[edof];
                int[] pMap = new int[npe];
                for (int k = 0; k < npe; k++)
                {
                    vMap[2 * k] = 2 * conn[k];
                    vMap[2 * k + 1] = 2 * conn[k] + 1;
                    pMap[k] = conn[k];
                }
                lock (viscBuilder)
                {
                    viscBuilder.AddLocal(vMap, visc);
                    massBuilder.AddLocal(vMap, me);
                    pMassBuilder.AddLocal(pMap, pMe);
                }

                // 集总质量（行和）
                for (int k = 0; k < npe; k++)
                {
                    double sum = 0;
                    for (int j = 0; j < npe; j++) sum += me[2 * k, 2 * j];
                    lumpMass[2 * conn[k]] += sum;
                    lumpMass[2 * conn[k] + 1] += sum;
                    double psum = 0;
                    for (int j = 0; j < npe; j++) psum += pMe[k, j];
                    lumpPMass[conn[k]] += psum;
                }
            });

            // 单独构建压力 Laplacian：∫∇N·∇N，按单元重新装配
            Parallel.For(0, nElem, e =>
            {
                FiniteElement fe = elems[e];
                int npe = mesh.Elements[e].Length;
                double[,] dummy = new double[2 * npe, 2 * npe];
                // 用粘性矩阵（mu=1）中的 u-u 子块作为标量 Laplacian：
                fe.ComputeViscous(mesh, e, 1.0, 1.0, dummy);
                double[,] pLap = new double[npe, npe];
                for (int i = 0; i < npe; i++)
                {
                    for (int j = 0; j < npe; j++)
                    {
                        pLap[i, j] = dummy[2 * i, 2 * j];
                    }
                }
                int[] conn = mesh.Elements[e];
                int[] pMap = new int[npe];
                for (int k = 0; k < npe; k++) pMap[k] = conn[k];
                lock (pLapBuilder) pLapBuilder.AddLocal(pMap, pLap);
            });

            SparseMatrixCSR ViscMat = viscBuilder.Build(true);
            SparseMatrixCSR PLapMat = pLapBuilder.Build(true);

            // 速度解与压力解
            DenseVector u = new DenseVector(ndof);
            DenseVector uStar = new DenseVector(ndof);
            DenseVector uOld = new DenseVector(ndof);
            // 初场：满足 Dirichlet 速度
            foreach ((int node, int dof, double val) in mesh.DirichletBCs)
            {
                u[3 * node + dof] = val;
                uOld[3 * node + dof] = val;
            }

            GMRESSolver gmres = new GMRESSolver
            {
                Tolerance = 1e-6,
                KrylovDim = 50,
                MaxRestart = 100
            };

            int totalSteps = Unsteady ? TimeSteps : 1;
            int totalIter = 0;
            double finalRes = 0;

            for (int ts = 0; ts < totalSteps; ts++)
            {
                if (log != null) log(Unsteady ? $"时间步 {ts + 1}/{totalSteps}, dt = {TimeStep}" : "开始定常 SIMPLE 迭代…");
                double res0 = 1.0;
                for (int it = 0; it < MaxOuterIter; it++)
                {
                    totalIter++;
                    // 1. 装配动量矩阵 A = (1/dt)*M + Visc + C(u) ，添加欠松弛（A/alpha_u）
                    SparseMatrixCSR.CooBuilder momBuilder = new SparseMatrixCSR.CooBuilder(vdof);
                    // 粘性
                    // 直接从已装配的 ViscMat 拷贝三元组
                    for (int r = 0; r < vdof; r++)
                    {
                        for (int k = ViscMat.RowPtr[r]; k < ViscMat.RowPtr[r + 1]; k++)
                        {
                            int c = ViscMat.ColIdx[k];
                            double v = ViscMat.Values[k];
                            momBuilder.Add(r, c, v);
                        }
                    }
                    // 对流（Picard 线性化）
                    Parallel.For(0, nElem, e =>
                    {
                        FiniteElement fe = elems[e];
                        int npe = mesh.Elements[e].Length;
                        int edof = 2 * npe;
                        double[] ue = new double[npe];
                        double[] ve = new double[npe];
                        int[] conn = mesh.Elements[e];
                        for (int k = 0; k < npe; k++)
                        {
                            ue[k] = u[3 * conn[k]];
                            ve[k] = u[3 * conn[k] + 1];
                        }
                        double[,] ce = new double[edof, edof];
                        fe.ComputeConvection(mesh, e, ue, ve, 1.0, ce);
                        int[] vMap = new int[edof];
                        for (int k = 0; k < npe; k++) { vMap[2 * k] = 3 * conn[k]; vMap[2 * k + 1] = 3 * conn[k] + 1; }
                        // 注意：速度自由度位于全局 0,1（跳过 p 间隔），重新映射
                        int[] vMapU = new int[edof];
                        for (int k = 0; k < npe; k++) { vMapU[2 * k] = 2 * conn[k]; vMapU[2 * k + 1] = 2 * conn[k] + 1; }
                        lock (momBuilder)
                        {
                            momBuilder.AddLocal(vMapU, ce);
                        }
                    });
                    // 时间项
                    if (Unsteady)
                    {
                        double dtInv = Density / TimeStep;
                        for (int i = 0; i < vdof; i++) momBuilder.Add(i, i, dtInv * lumpMass[i]);
                    }
                    // 欠松弛：A/AlphaU - (1/AlphaU - 1)*A_diag? 实际标准欠松弛：A_ua = A/alpha_u, rhs += (1/alpha_u - 1) A_diag u^old
                    // 简化：直接对 A 对角加强并调整右端（下面处理）。
                    SparseMatrixCSR Amat = momBuilder.Build(true);
                    double[] AdiagInv = Amat.DiagonalInverse();

                    // 2. 构造 rhs：-G*p + (时间项的 M*uOld/dt)
                    DenseVector rhs = new DenseVector(vdof);
                    // 压力梯度：-∫∇p·v 对应 -B * p_vec（B 对应 ∂N/∂x N_j 等，符号在单元中为 -∫∂N_i/∂x N_j）
                    Parallel.For(0, nElem, e =>
                    {
                        FiniteElement fe = elems[e];
                        int npe = mesh.Elements[e].Length;
                        int edof = 2 * npe;
                        double[,] B = new double[edof, npe];
                        fe.ComputePressureCoupling(mesh, e, 1.0, B);
                        int[] conn = mesh.Elements[e];
                        double[] pe = new double[npe];
                        for (int k = 0; k < npe; k++) pe[k] = u[3 * conn[k] + 2];
                        double[] re = new double[edof];
                        for (int i = 0; i < edof; i++)
                        {
                            double s = 0;
                            for (int j = 0; j < npe; j++) s += B[i, j] * pe[j];
                            // -∫ ∂v/∂x * p - ∫ ∂v/∂y * p → 动量方程 -grad p
                            re[i] = -s;
                        }
                        lock (rhs)
                        {
                            for (int k = 0; k < npe; k++)
                            {
                                rhs[2 * conn[k]] += re[2 * k];
                                rhs[2 * conn[k] + 1] += re[2 * k + 1];
                            }
                        }
                    });
                    // 时间项 rhs 加 M*uOld/dt
                    if (Unsteady)
                    {
                        double coef = Density / TimeStep;
                        for (int i = 0; i < vdof; i++)
                        {
                            rhs[i] += coef * lumpMass[i] * uOld[i];
                        }
                    }

                    // 速度 Dirichlet 边界
                    DenseVector uGuess = new DenseVector(vdof);
                    for (int i = 0; i < vdof; i++)
                    {
                        uGuess[i] = (i % 2 == 0) ? u[2 * (i / 2)] : u[2 * (i / 2) + 1];
                    }
                    // 欠松弛： rhs += (1/alpha_u - 1) diag(A) u_prev
                    for (int i = 0; i < vdof; i++)
                    {
                        double aii = 1.0 / AdiagInv[i];
                        rhs[i] += (1.0 / AlphaU - 1.0) * aii * uGuess[i];
                    }
                    // 施加边界条件到矩阵和右端
                    SparseMatrixCSR.CooBuilder momBC = new SparseMatrixCSR.CooBuilder(vdof);
                    for (int r = 0; r < vdof; r++)
                    {
                        for (int k = Amat.RowPtr[r]; k < Amat.RowPtr[r + 1]; k++)
                        {
                            momBC.Add(r, Amat.ColIdx[k], Amat.Values[k] / AlphaU);
                        }
                    }
                    foreach ((int node, int dof, double val) in mesh.DirichletBCs)
                    {
                        if (dof > 1) continue;
                        int row = 2 * node + dof;
                        double kv = (1.0 / AdiagInv[row]) * 1e14;
                        momBC.Add(row, row, kv);
                        rhs[row] = kv * val;
                    }
                    SparseMatrixCSR AmatBC = momBC.Build(true);
                    double[] AdiagInvBC = AmatBC.DiagonalInverse();

                    // 3. 求解动量方程
                    DenseVector uS_v = new DenseVector(vdof);
                    gmres.Solve(AmatBC, rhs, uS_v, AdiagInvBC);

                    // 4. 压力修正方程：L_p p' = -D * u^*
                    DenseVector divU = new DenseVector(pdof);
                    Parallel.For(0, nElem, e =>
                    {
                        FiniteElement fe = elems[e];
                        int npe = mesh.Elements[e].Length;
                        int edof = 2 * npe;
                        double[,] B = new double[edof, npe];
                        fe.ComputePressureCoupling(mesh, e, 1.0, B);
                        // 散度 D 对应 B^T（B 是 u 到 p 的梯度耦合，散度矩阵 D = -B^T）
                        int[] conn = mesh.Elements[e];
                        double[] us = new double[edof];
                        for (int k = 0; k < npe; k++)
                        {
                            us[2 * k] = uS_v[2 * conn[k]];
                            us[2 * k + 1] = uS_v[2 * conn[k] + 1];
                        }
                        double[] divE = new double[npe];
                        for (int j = 0; j < npe; j++)
                        {
                            double s = 0;
                            for (int i = 0; i < edof; i++)
                            {
                                s += B[i, j] * us[i];
                            }
                            // D u^* = -∫ N_j div u^* dΩ = ∫ ∇N_j · u 分部（略去边界）
                            divE[j] = -s;
                        }
                        lock (divU)
                        {
                            for (int j = 0; j < npe; j++)
                            {
                                divU[conn[j]] += divE[j];
                            }
                        }
                    });

                    // 压力矩阵：Ap = B^T diag(A)^{-1} B；SIMPLE 用 Ap ≈ PLap （标量 Laplacian）× 松弛
                    // 加上压力质量稳定化 τ * Mp 以消除伪压
                    SparseMatrixCSR.CooBuilder pBuilder = new SparseMatrixCSR.CooBuilder(pdof);
                    for (int r = 0; r < pdof; r++)
                    {
                        for (int k = PLapMat.RowPtr[r]; k < PLapMat.RowPtr[r + 1]; k++)
                        {
                            pBuilder.Add(r, PLapMat.ColIdx[k], PLapMat.Values[k]);
                        }
                    }

                    // PSPG 稳定化（简化）：加 τ_M * ∫ ∇p·∇q + τ_C * ∫∇p·∇q 用网格尺度
                    // Dirichlet 压力参考点
                    SparseMatrixCSR Ap = pBuilder.Build(true);
                    double[] pDiagInv = Ap.DiagonalInverse();
                    SparseMatrixCSR.CooBuilder pBC = new SparseMatrixCSR.CooBuilder(pdof);
                    for (int r = 0; r < pdof; r++)
                    {
                        for (int k = Ap.RowPtr[r]; k < Ap.RowPtr[r + 1]; k++)
                        {
                            pBC.Add(r, Ap.ColIdx[k], Ap.Values[k]);
                        }
                    }
                    foreach ((int node, int dof, double val) in mesh.DirichletBCs)
                    {
                        if (dof != 2) continue;
                        int row = node;
                        double kv = (1.0 / pDiagInv[row]) * 1e14;
                        pBC.Add(row, row, kv);
                        divU[row] = kv * val;
                    }
                    SparseMatrixCSR ApBC = pBC.Build(true);
                    double[] pDiagInvBC = ApBC.DiagonalInverse();
                    DenseVector pPrime = new DenseVector(pdof);
                    gmres.Tolerance = 1e-5;
                    gmres.Solve(ApBC, divU, pPrime, pDiagInvBC);

                    // 5. 修正速度和压力
                    double resN = 0.0;
                    // 压力修正
                    for (int i = 0; i < nNodes; i++) u[3 * i + 2] += AlphaP * pPrime[i];
                    // 速度修正：u = u^* - diag(A)^{-1} G p'
                    Parallel.For(0, nElem, e =>
                    {
                        FiniteElement fe = elems[e];
                        int npe = mesh.Elements[e].Length;
                        int edof = 2 * npe;
                        double[,] B = new double[edof, npe];
                        fe.ComputePressureCoupling(mesh, e, 1.0, B);
                        int[] conn = mesh.Elements[e];
                        double[] ppe = new double[npe];
                        for (int j = 0; j < npe; j++) ppe[j] = pPrime[conn[j]];
                        for (int i = 0; i < npe; i++)
                        {
                            double du = 0, dv = 0;
                            for (int j = 0; j < npe; j++)
                            {
                                du += B[2 * i, j] * ppe[j];
                                dv += B[2 * i + 1, j] * ppe[j];
                            }
                            double dIu = AdiagInvBC[2 * conn[i]];
                            double dIv = AdiagInvBC[2 * conn[i] + 1];
                            // uS_v - dI * G p'
                            double nu = uS_v[2 * conn[i]] - dIu * (-du); // G = -B，所以 -G p' = B p'
                            double nv = uS_v[2 * conn[i] + 1] - dIv * (-dv);
                            uS_v[2 * conn[i]] = nu;
                            uS_v[2 * conn[i] + 1] = nv;
                        }
                    });
                    // 写回全局速度（显式边界覆盖）
                    for (int i = 0; i < nNodes; i++)
                    {
                        u[3 * i] = uS_v[2 * i];
                        u[3 * i + 1] = uS_v[2 * i + 1];
                    }
                    foreach ((int node, int dof, double val) in mesh.DirichletBCs)
                    {
                        if (dof <= 1) u[3 * node + dof] = val;
                        else u[3 * node + dof] = val;
                    }
                    // 残差：速度修正量 L2
                    double res = 0.0;
                    for (int i = 0; i < vdof; i++) res += (uS_v[i] - uGuess[i]) * (uS_v[i] - uGuess[i]);
                    res = Math.Sqrt(res) / Math.Max(1e-12, uS_v.Norm2());
                    if (it == 0) res0 = Math.Max(res, 1e-12);
                    double relres = res / res0;
                    finalRes = relres;
                    if (log != null && it % 5 == 0) log($"SIMPLE 外迭代 {it + 1}, 速度相对残差 = {relres:E4}");
                    if (relres < ConvergenceTol && it > 2) break;
                }
                if (Unsteady) uOld.CopyFrom(u);
            }

            NSResult res2 = new NSResult();
            res2.Velocity = u;
            res2.U = new double[nNodes];
            res2.V = new double[nNodes];
            res2.P = new double[nNodes];
            res2.Speed = new double[nNodes];
            double maxS = 0;
            for (int i = 0; i < nNodes; i++)
            {
                res2.U[i] = u[3 * i];
                res2.V[i] = u[3 * i + 1];
                res2.P[i] = u[3 * i + 2];
                double s = Math.Sqrt(res2.U[i] * res2.U[i] + res2.V[i] * res2.V[i]);
                res2.Speed[i] = s;
                if (s > maxS) maxS = s;
            }
            res2.Iterations = totalIter;
            res2.FinalResidual = finalRes;
            return res2;
        }
    }
}
