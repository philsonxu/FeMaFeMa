namespace ModalFEM2D.Solvers
{
    using System;
    using System.Collections.Generic;
    using System.Threading;
    using ModalFEM2D.Core;
    using ModalFEM2D.Elements;
    using ModalFEM2D.Mesh;

    public enum KrylovKind { CG, GMRES, BiCGSTAB }
    public enum PrecondKind { None, Jacobi, AMG }

    public sealed class StructuralResult
    {
        public int Iterations;
        public double FinalResidual;
        public string SolverUsed;
        public double AssemblyTimeMs;
        public double SolveTimeMs;
    }

    public static class StructuralSolver
    {
        public static StructuralResult Solve(FEMesh mesh, KrylovKind krylov, PrecondKind precond,
                                              Action<string> log)
        {
            System.Diagnostics.Stopwatch sw = System.Diagnostics.Stopwatch.StartNew();
            int n = mesh.NumNodes * 2;
            int ne = mesh.NumElements;
            CooAssembler Kasm = new CooAssembler(n, n);
            DenseVector f = new DenseVector(n);
            // build loads
            for (int i = 0; i < mesh.NumNodes; i++)
            {
                f.Values[2 * i] = mesh.ForceU[i];
                f.Values[2 * i + 1] = mesh.ForceV[i];
            }
            double penalty = 1e20;
            // parallel assembly of K
            int threads = Environment.ProcessorCount;
            CooAssembler[] locals = new CooAssembler[threads];
            for (int t = 0; t < threads; t++) locals[t] = new CooAssembler(n, n);
            WaitHandle[] hs = new WaitHandle[threads];
            for (int t = 0; t < threads; t++)
            {
                int t0 = t;
                int s = t * ne / threads; int e = (t + 1) * ne / threads;
                ManualResetEvent mre = new ManualResetEvent(false); hs[t] = mre;
                ThreadPool.QueueUserWorkItem(_ =>
                {
                    for (int ei = s; ei < e; ei++)
                    {
                        FiniteElement el = mesh.CreateElement(ei);
                        double[,] Ke = el.BuildStiffnessMatrix();
                        int nd = 2 * el.NodesPerElement;
                        int[] dofMap = new int[nd];
                        for (int k = 0; k < el.NodesPerElement; k++)
                        { dofMap[2 * k] = 2 * el.NodeIds[k]; dofMap[2 * k + 1] = 2 * el.NodeIds[k] + 1; }
                        locals[t0].AddSymmLocal(dofMap, Ke);
                    }
                    mre.Set();
                });
            }
            WaitHandle.WaitAll(hs);
            for (int t = 0; t < threads; t++) { Kasm.Merge(locals[t]); hs[t].Dispose(); }
            // apply Dirichlet
            for (int i = 0; i < mesh.NumNodes; i++)
            {
                if (mesh.DirichletU[i] >= 0)
                {
                    int dof = 2 * i;
                    AddPenalty(Kasm, dof, dof, penalty);
                    f.Values[dof] = penalty * mesh.DirichletUValue[i];
                }
                if (mesh.DirichletV[i] >= 0)
                {
                    int dof = 2 * i + 1;
                    AddPenalty(Kasm, dof, dof, penalty);
                    f.Values[dof] = penalty * mesh.DirichletVValue[i];
                }
            }
            SparseMatrixCSR K = Kasm.Build();
            sw.Stop();
            double asmMs = sw.Elapsed.TotalMilliseconds;
            if (log != null) log(string.Format("[Structural] 装配完成, NNZ={0}, 耗时 {1:F1} ms", K.Values.Length, asmMs));

            // build preconditioner
            IPreconditioner M;
            if (precond == PrecondKind.AMG) { sw.Restart(); M = new AMGPreconditioner(K); sw.Stop(); if (log != null) log(string.Format("[Structural] AMG 构造完成, 耗时 {0:F1} ms", sw.Elapsed.TotalMilliseconds)); }
            else if (precond == PrecondKind.Jacobi) M = new JacobiPreconditioner(K);
            else M = new IdentityPreconditioner();

            DenseVector u = new DenseVector(n);
            IterationCallback cb = (it, rr) =>
            {
                if (log != null && it % 20 == 0) log(string.Format("[{0}] iter {1}, relRes = {2:E3}", krylov, it, rr));
            };
            sw.Restart();
            int iters = 0; string solverName = krylov.ToString();
            if (krylov == KrylovKind.CG) iters = KrylovSolvers.CG(K, f, u, M, 20000, 1e-8, cb);
            else if (krylov == KrylovKind.GMRES) iters = KrylovSolvers.GMRES(K, f, u, M, 60, 200, 1e-8, cb);
            else iters = KrylovSolvers.BiCGSTAB(K, f, u, M, 20000, 1e-8, cb);
            sw.Stop();
            double solMs = sw.Elapsed.TotalMilliseconds;
            if (log != null) log(string.Format("[{0}] 求解完成: {1} 次迭代, 用时 {2:F1} ms", solverName, iters, solMs));
            for (int i = 0; i < mesh.NumNodes; i++) { mesh.DisplacementU[i] = u.Values[2 * i]; mesh.DisplacementV[i] = u.Values[2 * i + 1]; }

            // recover element stresses -> average to nodes
            double[] sxxA = new double[mesh.NumNodes]; double[] syyA = new double[mesh.NumNodes];
            double[] sxyA = new double[mesh.NumNodes]; double[] vmA = new double[mesh.NumNodes];
            int[] cnt = new int[mesh.NumNodes];
            for (int ei = 0; ei < ne; ei++)
            {
                FiniteElement el = mesh.CreateElement(ei);
                int nd = 2 * el.NodesPerElement;
                double[] ue = new double[nd];
                for (int k = 0; k < el.NodesPerElement; k++)
                { ue[2 * k] = mesh.DisplacementU[el.NodeIds[k]]; ue[2 * k + 1] = mesh.DisplacementV[el.NodeIds[k]]; }
                el.ComputeStress(ue);
                for (int k = 0; k < el.NodesPerElement; k++)
                {
                    int nn = el.NodeIds[k];
                    sxxA[nn] += el.StressXX;
                    syyA[nn] += el.StressYY;
                    sxyA[nn] += el.StressXY;
                    vmA[nn] += el.VonMises;
                    cnt[nn]++;
                }
            }
            for (int i = 0; i < mesh.NumNodes; i++)
            {
                if (cnt[i] == 0) continue;
                mesh.StressXX[i] = sxxA[i] / cnt[i]; mesh.StressYY[i] = syyA[i] / cnt[i];
                mesh.StressXY[i] = sxyA[i] / cnt[i]; mesh.VonMises[i] = vmA[i] / cnt[i];
            }
            return new StructuralResult { Iterations = iters, FinalResidual = 1e-8, SolverUsed = solverName + "+" + precond, AssemblyTimeMs = asmMs, SolveTimeMs = solMs };
        }

        private static void AddPenalty(CooAssembler a, int r, int c, double v) { a.Add(r, c, v); }
    }
}
