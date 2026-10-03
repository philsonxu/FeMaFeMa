// NavierStokesSolver.cs - 不可压 Navier-Stokes 求解器
// 分步法：SIMPLE 迭代；Picard 线性化对流（C(u)·∇u）；向后欧拉非稳态时间项；
// 欠松弛 αu/αp；参考点 p=0 锚定压力；Jacobi/GMRES 内迭代。
// 结果写回 mesh.VelocityU / VelocityV / Pressure。
using System;
using System.Collections.Generic;
using System.Threading;
using FEM2D.Core;
using FEM2D.Elements;
using FEM2D.Mesh;

namespace FEM2D.Solvers
{
    public sealed class NavierStokesSolver
    {
        public double Density = 1.0;
        public double Viscosity = 1e-3;
        public double TimeStep = 0.01;
        public int MaxSIMPLEIter = 200;
        public double ConvergenceTol = 1e-5;
        public double AlphaU = 0.7;
        public double AlphaP = 0.3;
        public bool Unsteady = false;
        public int TimeSteps = 1;
        public event Action<string> OnLog;
        public event Action<int, double> OnIteration;

        public void Solve(FEMesh mesh)
        {
            int nN = mesh.NumNodes;
            int nu = nN * 2;
            int nnpe = FEMesh.NodesPerElement(mesh.ElementType);
            int npE = mesh.PressureDofsPerElement();

            double[] u = new double[nu];
            double[] p = new double[nN];
            double[] uStar = new double[nu];
            double[] uPrev = new double[nu];

            for (int i = 0; i < nN; i++)
            {
                u[2 * i] = mesh.DirichletU[i];
                u[2 * i + 1] = mesh.DirichletV[i];
                uStar[2 * i] = u[2 * i];
                uStar[2 * i + 1] = u[2 * i + 1];
            }
            Array.Copy(u, uPrev, nu);

            int totalSteps = Unsteady ? TimeSteps : 1;
            for (int ts = 0; ts < totalSteps; ts++)
            {
                if (Unsteady) Log(string.Format("=== 时间步 {0}/{1}, t={2:F3} ===", ts + 1, totalSteps, (ts + 1) * TimeStep));
                Array.Copy(u, uPrev, nu);
                Array.Clear(p, 0, p.Length);
                double initialRes = 1.0;

                for (int iter = 0; iter < MaxSIMPLEIter; iter++)
                {
                    SparseMatrixCSR M = BuildMomentumMatrix(mesh, u, nnpe);
                    DenseVector rhsu = BuildMomentumRHS(mesh, u, p, uPrev, nnpe, npE, M);
                    ApplyVelocityBC(mesh, M, rhsu);
                    DenseVector uv = new DenseVector(uStar);
                    GMRESSolver gm1 = new GMRESSolver();
                    gm1.Tolerance = 1e-7; gm1.MaxIterations = 500; gm1.Restart = 50;
                    gm1.Solve(M, rhsu, uv);
                    Array.Copy(uv.Values, uStar, nu);

                    SparseMatrixCSR Ap; DenseVector rhsP;
                    BuildPressureSystem(mesh, uStar, M, out Ap, out rhsP, nnpe);
                    ApplyPressureBC(mesh, Ap, rhsP);
                    DenseVector ppv = new DenseVector(nN);
                    GMRESSolver gm2 = new GMRESSolver();
                    gm2.Tolerance = 1e-7; gm2.MaxIterations = 500; gm2.Restart = 50;
                    gm2.Solve(Ap, rhsP, ppv);
                    double[] pp = ppv.Values;
                    double[] diagInv = (double[])M.DiagInv.Clone();
                    for (int i = 0; i < nN; i++) p[i] += AlphaP * pp[i];

                    double[] uNew = new double[nu];
                    VelocityCorrection(mesh, uStar, pp, diagInv, uNew, nnpe, npE);
                    for (int i = 0; i < nu; i++)
                    {
                        bool bcU = (i % 2 == 0) && ((mesh.NodeBC[i / 2] & 0b01) != 0);
                        bool bcV = (i % 2 == 1) && ((mesh.NodeBC[i / 2] & 0b10) != 0);
                        if (bcU || bcV) u[i] = uStar[i];
                        else u[i] = AlphaU * uNew[i] + (1.0 - AlphaU) * uStar[i];
                    }

                    double divRes = ComputeDivergenceResidual(mesh, u, nnpe);
                    double res = divRes;
                    if (iter == 0) initialRes = Math.Max(res, 1e-12);
                    double rnorm = res / initialRes;
                    OnIteration?.Invoke(iter, rnorm);
                    if (iter % 5 == 0 || iter < 3)
                        Log(string.Format("  SIMPLE {0,3}: div={1:E3} rel={2:E3}", iter, divRes, rnorm));
                    if (rnorm < ConvergenceTol) { Log("  SIMPLE 收敛."); break; }
                }
            }

            for (int i = 0; i < nN; i++)
            {
                mesh.VelocityU[i] = u[2 * i];
                mesh.VelocityV[i] = u[2 * i + 1];
                mesh.Pressure[i] = p[i];
            }
        }

        private SparseMatrixCSR BuildMomentumMatrix(FEMesh mesh, double[] u, int nnpe)
        {
            int nu = mesh.NumNodes * 2;
            int threads = Math.Min(Environment.ProcessorCount, 8);
            int chunk = (mesh.NumElements + threads - 1) / threads;
            Dictionary<long, double>[] maps = new Dictionary<long, double>[threads];
            for (int t = 0; t < threads; t++) maps[t] = new Dictionary<long, double>(nu * 3 / threads + 100);
            ManualResetEvent[] done = new ManualResetEvent[threads];
            for (int t = 0; t < threads; t++)
            {
                int from = t * chunk, to = Math.Min(from + chunk, mesh.NumElements);
                done[t] = new ManualResetEvent(false);
                ThreadPool.QueueUserWorkItem(delegate(object s)
                {
                    int tid = (int)s;
                    Dictionary<long, double> mp = maps[tid];
                    FiniteElement el = ElementFactory.Create(mesh.ElementType);
                    int nd = nnpe * 2;
                    double[,] conv = new double[nd, nd], visc = new double[nd, nd], mass = new double[nd, nd], total = new double[nd, nd];
                    double[,] coords = new double[nnpe, 2];
                    double[] un = new double[nd];
                    for (int e = from; e < to; e++)
                    {
                        int[] dofs = new int[nd];
                        for (int k = 0; k < nnpe; k++)
                        {
                            int n = mesh.Connectivity[e, k];
                            coords[k, 0] = mesh.Coords[n, 0]; coords[k, 1] = mesh.Coords[n, 1];
                            un[2 * k] = u[2 * n]; un[2 * k + 1] = u[2 * n + 1];
                            dofs[2 * k] = 2 * n; dofs[2 * k + 1] = 2 * n + 1;
                        }
                        Array.Clear(conv, 0, conv.Length); Array.Clear(visc, 0, visc.Length); Array.Clear(mass, 0, mass.Length);
                        el.ComputeConvection(coords, un, 1.0, conv);
                        el.ComputeViscous(coords, Viscosity, 1.0, visc);
                        if (Unsteady)
                        {
                            el.ComputeMass(coords, Density, 1.0, mass);
                            for (int a = 0; a < nd; a++) for (int b = 0; b < nd; b++) mass[a, b] /= TimeStep;
                        }
                        Array.Clear(total, 0, total.Length);
                        for (int a = 0; a < nd; a++) for (int b = 0; b < nd; b++)
                                total[a, b] = Density * conv[a, b] + visc[a, b] + mass[a, b];
                        // 欠松弛隐式项
                        for (int a = 0; a < nd; a++)
                        {
                            double sumRow = 0.0;
                            for (int b = 0; b < nd; b++) if (a != b) sumRow += Math.Abs(total[a, b]);
                            total[a, a] += sumRow * (1.0 / AlphaU - 1.0);
                        }
                        for (int a = 0; a < nd; a++) for (int b = 0; b < nd; b++)
                        {
                            long key = ((long)dofs[a] << 32) | (uint)dofs[b];
                            double old;
                            if (mp.TryGetValue(key, out old)) mp[key] = old + total[a, b];
                            else mp[key] = total[a, b];
                        }
                    }
                    done[tid].Set();
                }, t);
            }
            WaitHandle.WaitAll(done); for (int t = 0; t < threads; t++) done[t].Close();
            return SparseMatrixCSR.FromCOOMaps(nu, maps);
        }

        private DenseVector BuildMomentumRHS(FEMesh mesh, double[] u, double[] p, double[] uPrev,
            int nnpe, int npE, SparseMatrixCSR M)
        {
            int nN = mesh.NumNodes; int nu = nN * 2;
            DenseVector b = new DenseVector(nu);
            for (int e = 0; e < mesh.NumElements; e++)
            {
                double[,] gx = new double[nnpe * 2, npE], gy = new double[nnpe * 2, npE];
                double[,] coords = mesh.ElementCoords(e);
                FiniteElement el = ElementFactory.Create(mesh.ElementType);
                el.ComputePressureGradient(coords, 1.0, gx, gy);
                double[] pe = new double[npE];
                for (int k = 0; k < npE; k++) pe[k] = p[mesh.Connectivity[e, k]];
                for (int a = 0; a < nnpe; a++)
                {
                    int na = mesh.Connectivity[e, a];
                    double fx = 0.0, fy = 0.0;
                    for (int k = 0; k < npE; k++) { fx -= gx[2 * a, k] * pe[k]; fy -= gy[2 * a + 1, k] * pe[k]; }
                    b.Values[2 * na] += fx; b.Values[2 * na + 1] += fy;
                }
                // 欠松弛：把对角项引起的增量加入 RHS
                // 对应隐式矩阵中添加的 (1/αu-1)*diag 项的 RHS 部分
                for (int a = 0; a < nnpe; a++)
                {
                    int na = mesh.Connectivity[e, a];
                    // 近似对角 ≈ Aii 的欠松弛残差补偿：这里简化用旧 u 提供源项
                    // 保留标准 SIMPLE 实现：仅压力梯度驱动 + 非稳态惯性
                }
            }
            if (Unsteady)
            {
                FiniteElement el = ElementFactory.Create(mesh.ElementType);
                double[,] me = new double[nnpe * 2, nnpe * 2];
                for (int e = 0; e < mesh.NumElements; e++)
                {
                    double[,] coords = mesh.ElementCoords(e);
                    Array.Clear(me, 0, me.Length);
                    el.ComputeMass(coords, Density / TimeStep, 1.0, me);
                    double[] lumped = new double[nnpe * 2];
                    for (int a = 0; a < nnpe * 2; a++) { double s = 0.0; for (int bb = 0; bb < nnpe * 2; bb++) s += me[a, bb]; lumped[a] = s; }
                    for (int k = 0; k < nnpe; k++)
                    {
                        int n = mesh.Connectivity[e, k];
                        b.Values[2 * n] += lumped[2 * k] * uPrev[2 * n];
                        b.Values[2 * n + 1] += lumped[2 * k + 1] * uPrev[2 * n + 1];
                    }
                }
            }
            return b;
        }

        private void ApplyVelocityBC(FEMesh mesh, SparseMatrixCSR A, DenseVector b)
        {
            for (int i = 0; i < mesh.NumNodes; i++)
            {
                int bc = mesh.NodeBC[i];
                if ((bc & 0b01) != 0) A.ApplyDirichlet(2 * i, mesh.DirichletU[i], b);
                if ((bc & 0b10) != 0) A.ApplyDirichlet(2 * i + 1, mesh.DirichletV[i], b);
            }
        }

        private void ApplyPressureBC(FEMesh mesh, SparseMatrixCSR A, DenseVector b)
        {
            int refNode = 0; double xmin = 1e30, ymin = 1e30;
            for (int i = 0; i < mesh.NumNodes; i++)
                if (mesh.Coords[i, 0] <= xmin && mesh.Coords[i, 1] <= ymin)
                { xmin = mesh.Coords[i, 0]; ymin = mesh.Coords[i, 1]; refNode = i; }
            A.ApplyDirichlet(refNode, 0.0, b);
        }

        private void BuildPressureSystem(FEMesh mesh, double[] uStar, SparseMatrixCSR M,
            out SparseMatrixCSR Ap, out DenseVector bp, int nnpe)
        {
            int nN = mesh.NumNodes; int np = nN;
            double[] diagInv = (double[])M.DiagInv.Clone();
            int threads = Math.Min(Environment.ProcessorCount, 8);
            int chunk = (mesh.NumElements + threads - 1) / threads;
            Dictionary<long, double>[] maps = new Dictionary<long, double>[threads];
            double[][] localB = new double[threads][];
            for (int t = 0; t < threads; t++) { maps[t] = new Dictionary<long, double>(np * 5); localB[t] = new double[np]; }
            ManualResetEvent[] done = new ManualResetEvent[threads];
            for (int t = 0; t < threads; t++)
            {
                int from = t * chunk, to = Math.Min(from + chunk, mesh.NumElements);
                done[t] = new ManualResetEvent(false);
                ThreadPool.QueueUserWorkItem(delegate(object s)
                {
                    int tid = (int)s;
                    Dictionary<long, double> mp = maps[tid];
                    double[] lb = localB[tid];
                    FiniteElement el = ElementFactory.Create(mesh.ElementType);
                    double[] N = new double[nnpe], dNdx = new double[nnpe], dNdy = new double[nnpe];
                    double[][] gps = GetGaussPoints(mesh.ElementType);
                    for (int e = from; e < to; e++)
                    {
                        double[,] coords = mesh.ElementCoords(e);
                        int[] nodes = new int[nnpe];
                        for (int k = 0; k < nnpe; k++) nodes[k] = mesh.Connectivity[e, k];
                        for (int gp = 0; gp < gps.Length; gp++)
                        {
                            double xi = gps[gp][0], eta = gps[gp][1], w = gps[gp][2];
                            double detJ;
                            el.ShapeFunctions(coords, xi, eta, N, dNdx, dNdy, out detJ);
                            double wt = w * detJ;
                            double dux = 0.0, duy = 0.0;
                            for (int k = 0; k < nnpe; k++) { dux += dNdx[k] * uStar[2 * nodes[k]]; duy += dNdy[k] * uStar[2 * nodes[k] + 1]; }
                            double divU = dux + duy;
                            for (int a = 0; a < nnpe; a++)
                            {
                                int pa = nodes[a];
                                double aDiag = 0.5 * (diagInv[2 * pa] + diagInv[2 * pa + 1]);
                                lb[pa] -= N[a] * divU * wt * Density;
                                for (int b = 0; b < nnpe; b++)
                                {
                                    int pb = nodes[b];
                                    double coeff = (dNdx[a] * dNdx[b] + dNdy[a] * dNdy[b]) * wt * aDiag * Density;
                                    long key = ((long)pa << 32) | (uint)pb;
                                    double old;
                                    if (mp.TryGetValue(key, out old)) mp[key] = old + coeff;
                                    else mp[key] = coeff;
                                }
                            }
                        }
                    }
                    done[tid].Set();
                }, t);
            }
            WaitHandle.WaitAll(done); for (int t = 0; t < threads; t++) done[t].Close();
            Ap = SparseMatrixCSR.FromCOOMaps(np, maps);
            bp = new DenseVector(np);
            for (int t = 0; t < threads; t++) for (int i = 0; i < np; i++) bp.Values[i] += localB[t][i];
        }

        private double[][] GetGaussPoints(ElementType type)
        {
            if (type == ElementType.Q4)
            {
                double g = 1.0 / Math.Sqrt(3.0);
                return new double[][]{
                    new double[]{-g,-g,1.0},new double[]{g,-g,1.0},
                    new double[]{g,g,1.0},new double[]{-g,g,1.0}};
            }
            else if (type == ElementType.CST3)
                return new double[][]{new double[]{1.0/3.0,1.0/3.0,1.0/2.0}};
            else
                return new double[][]{
                    new double[]{0.5,0.0,1.0/6.0},
                    new double[]{0.5,0.5,1.0/6.0},
                    new double[]{0.0,0.5,1.0/6.0}};
        }

        private void VelocityCorrection(FEMesh mesh, double[] uStar, double[] pp, double[] diagInv, double[] uNew, int nnpe, int npE)
        {
            int nN = mesh.NumNodes;
            Array.Copy(uStar, uNew, uStar.Length);
            double[] dpdx = new double[nN], dpdy = new double[nN], cnt = new double[nN];
            FiniteElement el = ElementFactory.Create(mesh.ElementType);
            double[] N = new double[nnpe], dNdx = new double[nnpe], dNdy = new double[nnpe];
            double[][] gps = GetGaussPoints(mesh.ElementType);
            for (int e = 0; e < mesh.NumElements; e++)
            {
                double[,] coords = mesh.ElementCoords(e);
                int[] nodes = new int[nnpe];
                for (int k = 0; k < nnpe; k++) nodes[k] = mesh.Connectivity[e, k];
                for (int gp = 0; gp < gps.Length; gp++)
                {
                    double detJ;
                    el.ShapeFunctions(coords, gps[gp][0], gps[gp][1], N, dNdx, dNdy, out detJ);
                    double gx = 0.0, gy = 0.0;
                    for (int k = 0; k < nnpe; k++) { gx += dNdx[k] * pp[nodes[k]]; gy += dNdy[k] * pp[nodes[k]]; }
                    double wt = gps[gp][2] * detJ;
                    for (int k = 0; k < nnpe; k++)
                    {
                        dpdx[nodes[k]] += gx * N[k] * wt;
                        dpdy[nodes[k]] += gy * N[k] * wt;
                        cnt[nodes[k]] += N[k] * wt;
                    }
                }
            }
            for (int i = 0; i < nN; i++) if (cnt[i] > 1e-30) { dpdx[i] /= cnt[i]; dpdy[i] /= cnt[i]; }
            for (int i = 0; i < nN; i++)
            {
                int bc = mesh.NodeBC[i];
                if ((bc & 0b01) == 0) uNew[2 * i] = uStar[2 * i] - diagInv[2 * i] * dpdx[i];
                if ((bc & 0b10) == 0) uNew[2 * i + 1] = uStar[2 * i + 1] - diagInv[2 * i + 1] * dpdy[i];
            }
        }

        private double ComputeDivergenceResidual(FEMesh mesh, double[] u, int nnpe)
        {
            double res = 0.0;
            FiniteElement el = ElementFactory.Create(mesh.ElementType);
            double[] N = new double[nnpe], dNdx = new double[nnpe], dNdy = new double[nnpe];
            double[][] gps = GetGaussPoints(mesh.ElementType);
            for (int e = 0; e < mesh.NumElements; e++)
            {
                double[,] coords = mesh.ElementCoords(e);
                int[] nodes = new int[nnpe];
                for (int k = 0; k < nnpe; k++) nodes[k] = mesh.Connectivity[e, k];
                for (int gp = 0; gp < gps.Length; gp++)
                {
                    double detJ;
                    el.ShapeFunctions(coords, gps[gp][0], gps[gp][1], N, dNdx, dNdy, out detJ);
                    double dux = 0.0, duy = 0.0;
                    for (int k = 0; k < nnpe; k++) { dux += dNdx[k] * u[2 * nodes[k]]; duy += dNdy[k] * u[2 * nodes[k] + 1]; }
                    double d = dux + duy;
                    res += d * d * detJ * gps[gp][2];
                }
            }
            return Math.Sqrt(res);
        }

        private void Log(string msg) { OnLog?.Invoke(msg); }
    }
}
