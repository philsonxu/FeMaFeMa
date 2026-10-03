using System;
using System.Collections;
using System.Collections.Generic;

namespace ModalFEM2D.Solvers
{
    using System;
    using System.Threading;
    using ModalFEM2D.Core;
    using ModalFEM2D.Elements;
    using ModalFEM2D.Mesh;

    public sealed class NSResult
    {
        public int OuterIterations;
        public double FinalResidual;
        public double AssemblyMs;
        public double SolveMs;
        public bool Unsteady;
        public int TimeSteps;
        public double TimeStep;
    }

    /// <summary>
    /// Incompressible Navier-Stokes on collocated velocity-pressure nodes using SIMPLE.
    ///   (nu/rho) Lap(u) - (u·∇)u - ∇p/ρ = f;  div u = 0.
    /// Convection is Picard-linearised; transient uses backward Euler (lumped mass);
    /// equal-order interpolation stabilised with pressure Laplacian.
    /// </summary>
    public static class NavierStokesSolver
    {
        public static NSResult SolveSteady(FEMesh mesh, int maxOuter, double tol,
                                            double alphaU, double alphaP, double rho, double nu,
                                            Action<string> log)
        {
            return Run(mesh, maxOuter, tol, alphaU, alphaP, rho, nu, 0.0, 0, log);
        }

        public static NSResult SolveUnsteady(FEMesh mesh, int maxOuter, double tol,
                                              double alphaU, double alphaP, double rho, double nu,
                                              double dt, int nTimeSteps, Action<string> log)
        {
            return Run(mesh, maxOuter, tol, alphaU, alphaP, rho, nu, dt, nTimeSteps, log);
        }

        private static NSResult Run(FEMesh mesh, int maxOuter, double tol, double alphaU, double alphaP,
                                     double rho, double nu, double dt, int nTimeSteps, Action<string> log)
        {
            System.Diagnostics.Stopwatch total = System.Diagnostics.Stopwatch.StartNew();
            int n = mesh.NumNodes; int ndof = 2 * n; int ne = mesh.NumElements;
            mesh.ResetResults();
            for (int i = 0; i < n; i++)
            {
                mesh.VelocityU[i] = mesh.DirichletU[i] >= 0 ? mesh.DirichletUValue[i] : 0.0;
                mesh.VelocityV[i] = mesh.DirichletV[i] >= 0 ? mesh.DirichletVValue[i] : 0.0;
                mesh.Pressure[i] = 0.0;
            }
            double[] uOld = new double[ndof];
            NSResult res = new NSResult { Unsteady = dt > 0, TimeSteps = nTimeSteps, TimeStep = dt };
            int nSteps = dt > 0 ? Math.Max(1, nTimeSteps) : 1;
            double asmMs = 0.0, solMs = 0.0;
            int outerFinal = 0; double resFinal = 0.0;
            int threads = Environment.ProcessorCount;

            for (int ts = 0; ts < nSteps; ts++)
            {
                if (dt > 0)
                {
                    for (int i = 0; i < n; i++) { uOld[2 * i] = mesh.VelocityU[i]; uOld[2 * i + 1] = mesh.VelocityV[i]; }
                    if (log != null) log(string.Format("[NS] 时间步 {0}/{1}, dt={2:G4}", ts + 1, nSteps, dt));
                }
                int outer = 0; double maxDiv = 1e10;
                for (outer = 0; outer < maxOuter; outer++)
                {
                    System.Diagnostics.Stopwatch sw = System.Diagnostics.Stopwatch.StartNew();

                    // ============ Momentum matrix ============
                    CooAssembler MomAsm = new CooAssembler(ndof, ndof);
                    // Parallel assembly
                    CooAssembler[] locals = new CooAssembler[threads];
                    for (int t = 0; t < threads; t++) locals[t] = new CooAssembler(ndof, ndof);
                    WaitHandle[] hs = new WaitHandle[threads];
                    for (int t = 0; t < threads; t++)
                    {
                        int t0 = t; int s = t * ne / threads; int e = (t + 1) * ne / threads;
                        ManualResetEvent mre = new ManualResetEvent(false); hs[t] = mre;
                        ThreadPool.QueueUserWorkItem(_ =>
                        {
                            for (int ei = s; ei < e; ei++)
                            {
                                FiniteElement el = mesh.CreateElement(ei);
                                el.Viscosity = nu; el.Density = rho;
                                double[] eu = new double[el.NodesPerElement];
                                double[] ev = new double[el.NodesPerElement];
                                for (int k = 0; k < el.NodesPerElement; k++)
                                { eu[k] = mesh.VelocityU[el.NodeIds[k]]; ev[k] = mesh.VelocityV[el.NodeIds[k]]; }
                                double[,] visc = el.BuildViscousMatrix();
                                double[,] conv = el.BuildConvectionMatrix(eu, ev);
                                int nd = 2 * el.NodesPerElement;
                                int[] uMap = new int[nd];
                                for (int k = 0; k < el.NodesPerElement; k++)
                                { uMap[2 * k] = 2 * el.NodeIds[k]; uMap[2 * k + 1] = 2 * el.NodeIds[k] + 1; }
                                double[,] Ke = new double[nd, nd];
                                for (int i = 0; i < nd; i++) for (int j = 0; j < nd; j++) Ke[i, j] = (visc[i, j] + conv[i, j]) / rho;
                                if (dt > 0.0)
                                {
                                    double[] lm = el.BuildLumpedMass();
                                    for (int i = 0; i < nd; i++) Ke[i, i] += lm[i] / dt;
                                }
                                locals[t0].AddSymmLocal(uMap, Ke);
                            }
                            mre.Set();
                        });
                    }
                    WaitHandle.WaitAll(hs);
                    for (int t = 0; t < threads; t++) { MomAsm.Merge(locals[t]); hs[t].Dispose(); }

                    // Assemble rhs = 0 for steady; for transient: rhs = M_lump/dt * uOld
                    DenseVector rhsM = new DenseVector(ndof);
                    if (dt > 0)
                    {
                        for (int ei = 0; ei < ne; ei++)
                        {
                            FiniteElement el = mesh.CreateElement(ei);
                            double[] lm = el.BuildLumpedMass();
                            for (int k = 0; k < el.NodesPerElement; k++)
                            {
                                int gid = el.NodeIds[k];
                                rhsM.Values[2 * gid] += lm[2 * k] / dt * uOld[2 * gid];
                                rhsM.Values[2 * gid + 1] += lm[2 * k + 1] / dt * uOld[2 * gid + 1];
                            }
                        }
                    }
                    // apply velocity Dirichlet (penalty method for velocity)
                    double penalty = 1e18;
                    DenseVector uCur = new DenseVector(ndof);
                    for (int i = 0; i < n; i++)
                    {
                        uCur.Values[2 * i] = mesh.VelocityU[i];
                        uCur.Values[2 * i + 1] = mesh.VelocityV[i];
                    }
                    // For implicit under-relaxation we solve: (A/diagFactor) u_new = rhs + (1-αu)/αu A_diag u_old
                    // Easier approach: modify matrix by adding (1/αu - 1) * diag(A) to diag and scale rhs/old term.
                    // We'll implement it by adding the extra diagonal after building CSR.
                    for (int i = 0; i < n; i++)
                    {
                        if (mesh.DirichletU[i] >= 0)
                        { MomAsm.Add(2 * i, 2 * i, penalty); rhsM.Values[2 * i] = penalty * mesh.DirichletUValue[i]; }
                        if (mesh.DirichletV[i] >= 0)
                        { MomAsm.Add(2 * i + 1, 2 * i + 1, penalty); rhsM.Values[2 * i + 1] = penalty * mesh.DirichletVValue[i]; }
                    }
                    SparseMatrixCSR A_base = MomAsm.Build();
                    // Apply under-relaxation: build A_ur = A_base / alpha_u with (1/αu - 1) * diag augmented
                    List<Tuple<int, int, double>> triples = new List<Tuple<int, int, double>>();
                    double[] extraDiag = new double[ndof];
                    for (int i = 0; i < ndof; i++)
                    {
                        double d = 0.0;
                        for (int kk = A_base.RowPtr[i]; kk < A_base.RowPtr[i + 1]; kk++)
                            if (A_base.ColIdx[kk] == i) { d = A_base.Values[kk]; break; }
                        extraDiag[i] = (1.0 / alphaU - 1.0) * d;
                    }
                    for (int kk = 0; kk < A_base.Values.Length; kk++)
                    {
                        // determine row
                        triples.Add(Tuple.Create(0, A_base.ColIdx[kk], 0.0)); // placeholder; we'll rebuild row-by-row below
                    }
                    triples.Clear();
                    for (int i = 0; i < ndof; i++)
                    {
                        for (int kk = A_base.RowPtr[i]; kk < A_base.RowPtr[i + 1]; kk++)
                        {
                            int c = A_base.ColIdx[kk]; double v = A_base.Values[kk] / alphaU;
                            if (c == i) v += extraDiag[i];
                            triples.Add(Tuple.Create(i, c, v));
                        }
                        rhsM.Values[i] += extraDiag[i] * uCur.Values[i];
                    }
                    SparseMatrixCSR A_ur = SparseMatrixCSR.FromCoo(ndof, ndof, triples);
                    sw.Stop(); asmMs += sw.Elapsed.TotalMilliseconds; sw.Restart();
                    IPreconditioner jp1 = new JacobiPreconditioner(A_ur);
                    DenseVector uStar = uCur.Clone();
                    int mIt = KrylovSolvers.BiCGSTAB(A_ur, rhsM, uStar, jp1, 500, 1e-6, null);
                    sw.Stop(); solMs += sw.Elapsed.TotalMilliseconds; sw.Reset();

                    // ============ Pressure Poisson equation ============
                    // Ap = -D * (A_base.diag)^{-1} * G ; rhs = -D u* + stab.
                    sw.Start();
                    double[] invDiagBase = A_base.DiagInv;
                    CooAssembler Pasm = new CooAssembler(n, n);
                    DenseVector prhs = new DenseVector(n);
                    // We need also lumped-mass / area weighting for diagInv better. Simple diag is typical in textbook SIMPLE.
                    double stabCoef = 0.3;
                    for (int ei = 0; ei < ne; ei++)
                    {
                        FiniteElement el = mesh.CreateElement(ei);
                        double[,] G = el.BuildPressureGradient();
                        double[,] Dv = el.BuildDivergence();
                        double[,] Lp = el.BuildPressureLaplacian();
                        int pn = el.NodesPerElement;
                        int[] pMap = new int[pn];
                        double[] us = new double[2 * pn];
                        double[] invd = new double[2 * pn];
                        for (int k = 0; k < pn; k++)
                        {
                            int gid = el.NodeIds[k]; pMap[k] = gid;
                            us[2 * k] = uStar.Values[2 * gid];
                            us[2 * k + 1] = uStar.Values[2 * gid + 1];
                            invd[2 * k] = Math.Abs(invDiagBase[2 * gid]) > 1e-30 ? invDiagBase[2 * gid] : 1.0;
                            invd[2 * k + 1] = Math.Abs(invDiagBase[2 * gid + 1]) > 1e-30 ? invDiagBase[2 * gid + 1] : 1.0;
                        }
                        double[,] Ap = new double[pn, pn];
                        double[] rhsE = new double[pn];
                        for (int i = 0; i < pn; i++)
                        {
                            for (int jj = 0; jj < 2 * pn; jj++) rhsE[i] -= Dv[i, jj] * us[jj];
                            for (int j = 0; j < pn; j++)
                            {
                                double s = 0.0;
                                for (int jj = 0; jj < 2 * pn; jj++) s += Dv[i, jj] * invd[jj] * G[jj, j];
                                Ap[i, j] = -s + stabCoef * Lp[i, j];
                            }
                        }
                        Pasm.AddLocal(pMap, pMap, Ap);
                        for (int i = 0; i < pn; i++) prhs.Values[pMap[i]] += rhsE[i];
                    }
                    foreach (BCEntry bc in mesh.BCs) if (bc.Kind == BCKind.Pressure)
                    {
                        Pasm.Add(bc.Node, bc.Node, penalty); prhs.Values[bc.Node] = penalty * bc.Value;
                    }
                    SparseMatrixCSR ApMat = Pasm.Build();
                    sw.Stop(); asmMs += sw.Elapsed.TotalMilliseconds; sw.Restart();
                    IPreconditioner jp2 = new JacobiPreconditioner(ApMat);
                    DenseVector pNew = new DenseVector(n);
                    for (int i = 0; i < n; i++) pNew.Values[i] = mesh.Pressure[i];
                    int pIt = KrylovSolvers.BiCGSTAB(ApMat, prhs, pNew, jp2, 800, 1e-5, null);
                    sw.Stop(); solMs += sw.Elapsed.TotalMilliseconds; sw.Reset();

                    // ============ Velocity / pressure correction ============
                    sw.Start();
                    for (int ei = 0; ei < ne; ei++)
                    {
                        FiniteElement el = mesh.CreateElement(ei);
                        double[,] G = el.BuildPressureGradient();
                        int pn = el.NodesPerElement;
                        double[] pL = new double[pn];
                        int[] gids = new int[pn];
                        double[] invd = new double[2 * pn];
                        for (int k = 0; k < pn; k++)
                        {
                            gids[k] = el.NodeIds[k];
                            pL[k] = pNew.Values[gids[k]] - mesh.Pressure[gids[k]]; // dp correction
                            invd[2 * k] = Math.Abs(invDiagBase[2 * gids[k]]) > 1e-30 ? invDiagBase[2 * gids[k]] : 1.0;
                            invd[2 * k + 1] = Math.Abs(invDiagBase[2 * gids[k] + 1]) > 1e-30 ? invDiagBase[2 * gids[k] + 1] : 1.0;
                        }
                        double[] dp = new double[2 * pn];
                        for (int jj = 0; jj < 2 * pn; jj++) for (int j = 0; j < pn; j++) dp[jj] += G[jj, j] * pL[j];
                        for (int k = 0; k < pn; k++)
                        {
                            int gid = gids[k];
                            if (mesh.DirichletU[gid] < 0) mesh.VelocityU[gid] = uStar.Values[2 * gid] - alphaU * invd[2 * k] * dp[2 * k];
                            else mesh.VelocityU[gid] = mesh.DirichletUValue[gid];
                            if (mesh.DirichletV[gid] < 0) mesh.VelocityV[gid] = uStar.Values[2 * gid + 1] - alphaU * invd[2 * k + 1] * dp[2 * k + 1];
                            else mesh.VelocityV[gid] = mesh.DirichletVValue[gid];
                            mesh.Pressure[gid] += alphaP * pL[k];
                        }
                    }
                    // residual
                    maxDiv = 0.0;
                    for (int ei = 0; ei < ne; ei++)
                    {
                        FiniteElement el = mesh.CreateElement(ei);
                        double[,] Dv = el.BuildDivergence();
                        int pn = el.NodesPerElement;
                        double[] us = new double[2 * pn];
                        for (int k = 0; k < pn; k++)
                        { int gid = el.NodeIds[k]; us[2 * k] = mesh.VelocityU[gid]; us[2 * k + 1] = mesh.VelocityV[gid]; }
                        for (int i = 0; i < pn; i++)
                        {
                            double dv = 0.0;
                            for (int jj = 0; jj < 2 * pn; jj++) dv += Dv[i, jj] * us[jj];
                            double a = Math.Abs(dv); if (a > maxDiv) maxDiv = a;
                        }
                    }
                    sw.Stop(); asmMs += sw.Elapsed.TotalMilliseconds;
                    outerFinal = outer + 1; resFinal = maxDiv;
                    if (log != null && ((outer + 1) % 5 == 0 || outer == 0))
                        log(string.Format("[SIMPLE] outer={0} u_it={1} p_it={2} max|div u|={3:E3}", outer + 1, mIt, pIt, maxDiv));
                    if (maxDiv < tol) break;
                }
                if (log != null) log(string.Format("[SIMPLE] 本次{0}收敛: {1} 次迭代, |div|={2:E3}", dt > 0 ? "时间步" : "步", outerFinal, resFinal));
                if (dt == 0.0) break;
            }
            total.Stop();
            res.OuterIterations = outerFinal; res.FinalResidual = resFinal;
            res.AssemblyMs = asmMs; res.SolveMs = solMs;
            if (log != null) log(string.Format("[NS] 全部完成, 总耗时 {0:F1} ms (asm {1:F1}, sol {2:F1})", total.Elapsed.TotalMilliseconds, asmMs, solMs));
            return res;
        }
    }
}
