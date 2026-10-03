namespace ModalFEM2D.Solvers
{
    using System;
    using System.Collections.Generic;
    using System.Threading;
    using ModalFEM2D.Core;
    using ModalFEM2D.Elements;
    using ModalFEM2D.Mesh;

    public sealed class ModalResult
    {
        public int TargetModes;
        public double[] FrequenciesHz;
        public double[] OmegaRad;
        public int SubspaceIterations;
        public double AssemblyMs;
        public double SolveMs;
        public string SolverUsed = "Subspace Iteration + CG/AMG";
    }

    /// <summary>
    /// Modal analysis via subspace iteration (inverse iteration on a block of p+4 vectors,
    /// Rayleigh-Ritz projection, M-orthonormal Gram-Schmidt). Inner solve uses CG+AMG.
    /// We impose Dirichlet constraints by penalising the stiffness matrix (same as structural),
    /// so fixed DOFs produce high-frequency modes that are discarded after extraction.
    /// </summary>
    public static class ModalSolver
    {
        public static ModalResult Solve(FEMesh mesh, int numModes, Action<string> log)
        {
            System.Diagnostics.Stopwatch sw = System.Diagnostics.Stopwatch.StartNew();
            int n = mesh.NumNodes * 2;
            int ne = mesh.NumElements;
            CooAssembler Kasm = new CooAssembler(n, n);
            CooAssembler Massm = new CooAssembler(n, n);
            int threads = Environment.ProcessorCount;
            CooAssembler[] Ks = new CooAssembler[threads];
            CooAssembler[] Ms = new CooAssembler[threads];
            for (int t = 0; t < threads; t++) { Ks[t] = new CooAssembler(n, n); Ms[t] = new CooAssembler(n, n); }
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
                        double[,] Me = el.BuildMassMatrix();
                        int nd = 2 * el.NodesPerElement;
                        int[] dofMap = new int[nd];
                        for (int k = 0; k < el.NodesPerElement; k++)
                        { dofMap[2 * k] = 2 * el.NodeIds[k]; dofMap[2 * k + 1] = 2 * el.NodeIds[k] + 1; }
                        Ks[t0].AddSymmLocal(dofMap, Ke);
                        Ms[t0].AddSymmLocal(dofMap, Me);
                    }
                    mre.Set();
                });
            }
            WaitHandle.WaitAll(hs);
            for (int t = 0; t < threads; t++) { Kasm.Merge(Ks[t]); Massm.Merge(Ms[t]); hs[t].Dispose(); }
            // Apply Dirichlet: big penalty on K (same as structural), zero mass rows (or small diag)
            double penalty = 1e18;
            for (int i = 0; i < mesh.NumNodes; i++)
            {
                if (mesh.DirichletU[i] >= 0) { int d = 2 * i; Kasm.Add(d, d, penalty); Massm.Add(d, d, 1.0); }
                if (mesh.DirichletV[i] >= 0) { int d = 2 * i + 1; Kasm.Add(d, d, penalty); Massm.Add(d, d, 1.0); }
            }
            SparseMatrixCSR K = Kasm.Build();
            SparseMatrixCSR M = Massm.Build();
            sw.Stop();
            double asmMs = sw.Elapsed.TotalMilliseconds;
            if (log != null) log(string.Format("[Modal] K/M 装配完成, K.NNZ={0}, M.NNZ={1}, {2:F1} ms", K.Values.Length, M.Values.Length, asmMs));

            // AMG preconditioner for the inner solves (on K with penalties; high-frequency modes handled by subspace projection)
            sw.Restart();
            IPreconditioner prec = new AMGPreconditioner(K);
            sw.Stop();
            if (log != null) log(string.Format("[Modal] AMG 预条件构造完成 {0:F1} ms", sw.Elapsed.TotalMilliseconds));

            // subspace size
            int p = Math.Min(numModes, n / 2);
            int m = Math.Min(2 * p + 4, n);
            if (log != null) log(string.Format("[Modal] 子空间维数 m={0}, 目标阶次 p={1}", m, p));

            // initial guess X: random
            Random rnd = new Random(1234);
            DenseVector[] X = new DenseVector[m];
            for (int j = 0; j < m; j++)
            {
                X[j] = new DenseVector(n);
                for (int i = 0; i < n; i++) X[j].Values[i] = rnd.NextDouble() - 0.5;
            }
            MOrthonormalize(X, M);

            int maxOut = 60;
            double[] lam = new double[m];
            DenseVector[] Phi = new DenseVector[m];
            int outIt = 0;
            sw.Restart();
            for (outIt = 0; outIt < maxOut; outIt++)
            {
                // Y = K^{-1} M X (block inverse iteration). Each column: solve K y = M x.
                DenseVector[] Y = new DenseVector[m];
                for (int j = 0; j < m; j++)
                {
                    DenseVector rhs = new DenseVector(n);
                    M.MatVec(X[j], rhs);
                    Y[j] = new DenseVector(n);
                    // warm-start from X[j] for faster convergence
                    Y[j].CopyFrom(X[j]);
                    KrylovSolvers.CG(K, rhs, Y[j], prec, 300, 1e-10, null);
                }
                // project K and M onto span(Y): Kp = Y^T K Y, Mp = Y^T M Y
                // K Y gives K*Yj cheaply via matvec (we already have it since K Yj = M Xj by construction, exactly! So K Y = M X = rhs.)
                // Use that identity to avoid an extra K*Y multiply.
                DenseVector[] KY = new DenseVector[m];
                for (int j = 0; j < m; j++) { KY[j] = new DenseVector(n); M.MatVec(X[j], KY[j]); }
                DenseVector[] MY = new DenseVector[m];
                for (int j = 0; j < m; j++) { MY[j] = new DenseVector(n); M.MatVec(Y[j], MY[j]); }

                double[,] Kp = new double[m, m];
                double[,] Mp = new double[m, m];
                for (int i = 0; i < m; i++)
                    for (int j = 0; j < m; j++)
                    { Kp[i, j] = Y[i].Dot(KY[j]); Mp[i, j] = Y[i].Dot(MY[j]); }

                // Solve small generalized eigenvalue problem Kp q = λ Mp q (symmetric, Mp SPD).
                // Use Cholesky of Mp -> standard symmetric eigen -> Jacobi rotations.
                double[,] L = Cholesky(Mp, m);
                double[,] Linv = InvertLowerTri(L, m);
                double[,] As = new double[m, m];
                // A = L^{-1} Kp L^{-T}
                double[,] Ks2 = new double[m, m];
                for (int i = 0; i < m; i++) for (int j = 0; j < m; j++)
                    for (int k = 0; k < m; k++) Ks2[i, j] += Linv[i, k] * Kp[k, j];
                for (int i = 0; i < m; i++) for (int j = 0; j < m; j++)
                    for (int k = 0; k < m; k++) As[i, j] += Ks2[i, k] * Linv[j, k];
                // symmetrize
                for (int i = 0; i < m; i++) for (int j = i + 1; j < m; j++) { As[j, i] = (As[i, j] + As[j, i]) * 0.5; As[i, j] = As[j, i]; }
                double[,] Q = JacobiEigen(As, m, 200, 1e-10, lam);
                // eigenvectors of original problem: Z = L^{-T} Q
                double[,] Z = new double[m, m];
                for (int i = 0; i < m; i++)
                    for (int j = 0; j < m; j++)
                        for (int k = 0; k < m; k++) Z[i, j] += Linv[k, i] * Q[k, j];

                // Ritz vectors: X_new = Y Z
                DenseVector[] Xnew = new DenseVector[m];
                for (int j = 0; j < m; j++) Xnew[j] = new DenseVector(n);
                for (int j = 0; j < m; j++)
                    for (int i = 0; i < m; i++)
                        Xnew[j].Axpy(Z[i, j], Y[i]);
                // M-orthonormalize again for numerical safety
                MOrthonormalize(Xnew, M);
                X = Xnew;
                // check convergence: ratio change of eigenvalues
                Array.Sort(lam);
                // ignore ~penalty modes (very large eigenvalues)
                int keep = 0;
                for (int j = 0; j < m; j++) if (lam[j] < penalty * 1e-6) keep++;
                keep = Math.Min(keep, p);
                if (log != null && (outIt % 5 == 0 || outIt == maxOut - 1))
                {
                    string firstFreq = "";
                    for (int j = 0; j < Math.Min(keep, 4); j++)
                    {
                        double om = Math.Sqrt(Math.Max(lam[j], 0.0));
                        firstFreq += string.Format("f{0}={1:F2}Hz ", j + 1, om / (2 * Math.PI));
                    }
                    log(string.Format("[Modal] 迭代 {0}, 保留 {1} 阶, {2}", outIt + 1, keep, firstFreq));
                }
                if (keep >= p && outIt > 5) break;
            }
            sw.Stop();
            double solMs = sw.Elapsed.TotalMilliseconds;
            // Extract eigenvalues/vectors in ascending order; skip penalty modes
            Array.Sort(lam);
            List<double> keepList = new List<double>();
            for (int j = 0; j < m; j++) if (lam[j] < 1e15) keepList.Add(lam[j]);
            int take = Math.Min(numModes, keepList.Count);
            double[] omegas = new double[take];
            double[] freqs = new double[take];
            for (int j = 0; j < take; j++)
            {
                double lamj = keepList[j];
                omegas[j] = Math.Sqrt(Math.Max(lamj, 0.0));
                freqs[j] = omegas[j] / (2.0 * Math.PI);
            }
            // For shapes we need the actual eigenvectors corresponding to sorted eigenvalues.
            // Re-run Rayleigh-Ritz once more to get eigenvectors aligned with sorted lam, cheap.
            // (We discarded Z above; redo projection on the final X.)
            double[,] Kp2 = new double[m, m];
            double[,] Mp2 = new double[m, m];
            DenseVector[] KX = new DenseVector[m];
            DenseVector[] MX = new DenseVector[m];
            for (int j = 0; j < m; j++) { KX[j] = new DenseVector(n); K.MatVec(X[j], KX[j]); MX[j] = new DenseVector(n); M.MatVec(X[j], MX[j]); }
            for (int i = 0; i < m; i++) for (int j = 0; j < m; j++)
            { Kp2[i, j] = X[i].Dot(KX[j]); Mp2[i, j] = X[i].Dot(MX[j]); }
            double[,] L2 = Cholesky(Mp2, m);
            double[,] Li2 = InvertLowerTri(L2, m);
            double[,] A2 = new double[m, m];
            double[,] tmp = new double[m, m];
            for (int i = 0; i < m; i++) for (int j = 0; j < m; j++) for (int k = 0; k < m; k++) tmp[i, j] += Li2[i, k] * Kp2[k, j];
            for (int i = 0; i < m; i++) for (int j = 0; j < m; j++) for (int k = 0; k < m; k++) A2[i, j] += tmp[i, k] * Li2[j, k];
            for (int i = 0; i < m; i++) for (int j = i + 1; j < m; j++) { A2[j, i] = (A2[i, j] + A2[j, i]) * 0.5; A2[i, j] = A2[j, i]; }
            double[,] Q2 = JacobiEigen(A2, m, 200, 1e-10, lam);
            double[,] Z2 = new double[m, m];
            for (int i = 0; i < m; i++) for (int j = 0; j < m; j++) for (int k = 0; k < m; k++) Z2[i, j] += Li2[k, i] * Q2[k, j];
            // sort (lam, col Z2) ascending by lam
            int[] idx = new int[m]; for (int i = 0; i < m; i++) idx[i] = i;
            Array.Sort(lam, idx);
            mesh.ModalFreqHz = new double[take];
            mesh.ModalShapes = new double[take, n];
            int stored = 0;
            for (int j = 0; j < m && stored < take; j++)
            {
                int col = idx[j];
                if (lam[j] >= 1e15) continue;
                mesh.ModalFreqHz[stored] = Math.Sqrt(Math.Max(lam[j], 0.0)) / (2.0 * Math.PI);
                DenseVector vec = new DenseVector(n);
                for (int i = 0; i < m; i++) vec.Axpy(Z2[i, col], X[i]);
                // normalise to max displacement = 1 for visualisation
                double mx = 0.0;
                for (int i = 0; i < n; i++) if (Math.Abs(vec.Values[i]) > mx) mx = Math.Abs(vec.Values[i]);
                if (mx > 1e-12) vec.Scale(1.0 / mx);
                for (int i = 0; i < n; i++) mesh.ModalShapes[stored, i] = vec.Values[i];
                stored++;
            }
            // store first mode into displacement for immediate display
            if (take > 0)
            {
                for (int i = 0; i < mesh.NumNodes; i++)
                {
                    mesh.DisplacementU[i] = mesh.ModalShapes[0, 2 * i];
                    mesh.DisplacementV[i] = mesh.ModalShapes[0, 2 * i + 1];
                }
            }
            return new ModalResult
            {
                TargetModes = take,
                FrequenciesHz = mesh.ModalFreqHz,
                OmegaRad = (double[])mesh.ModalFreqHz.Clone(),
                SubspaceIterations = outIt,
                AssemblyMs = asmMs,
                SolveMs = solMs,
                SolverUsed = "Subspace Inv+Rayleigh-Ritz+CG+AMG"
            };
        }

        private static void MOrthonormalize(DenseVector[] X, SparseMatrixCSR M)
        {
            int m = X.Length; int n = M.NumRows;
            DenseVector Mx = new DenseVector(n);
            for (int j = 0; j < m; j++)
            {
                // M-orthogonalise against previous vectors
                M.MatVec(X[j], Mx);
                for (int i = 0; i < j; i++)
                {
                    DenseVector Mxi = new DenseVector(n);
                    M.MatVec(X[i], Mxi);
                    double c = X[i].Dot(Mx) / Math.Max(X[i].Dot(Mxi), 1e-30);
                    X[j].Axpy(-c, X[i]);
                    M.MatVec(X[j], Mx);
                }
                double nrm = Math.Sqrt(Math.Max(X[j].Dot(Mx), 1e-30));
                X[j].Scale(1.0 / nrm);
            }
        }

        private static double[,] Cholesky(double[,] A, int n)
        {
            double[,] L = new double[n, n];
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j <= i; j++)
                {
                    double s = A[i, j];
                    for (int k = 0; k < j; k++) s -= L[i, k] * L[j, k];
                    if (i == j) L[i, j] = Math.Sqrt(Math.Max(s, 1e-30));
                    else L[i, j] = s / L[j, j];
                }
            }
            return L;
        }

        private static double[,] InvertLowerTri(double[,] L, int n)
        {
            double[,] Li = new double[n, n];
            for (int i = 0; i < n; i++)
            {
                Li[i, i] = 1.0 / L[i, i];
                for (int j = 0; j < i; j++)
                {
                    double s = 0.0;
                    for (int k = j; k < i; k++) s -= L[i, k] * Li[k, j];
                    Li[i, j] = s / L[i, i];
                }
            }
            return Li;
        }

        /// <summary>Classic Jacobi eigenvalue algorithm for symmetric matrix A. Returns Q with eigenvectors in columns, eigenvalues in lam (ascending).</summary>
        private static double[,] JacobiEigen(double[,] A, int n, int maxSweeps, double tol, double[] lam)
        {
            double[,] D = (double[,])A.Clone();
            double[,] Q = new double[n, n];
            for (int i = 0; i < n; i++) Q[i, i] = 1.0;
            for (int sweep = 0; sweep < maxSweeps; sweep++)
            {
                double off = 0.0;
                for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) off += D[i, j] * D[i, j];
                if (Math.Sqrt(off) < tol) break;
                for (int p = 0; p < n; p++) for (int q = p + 1; q < n; q++)
                {
                    double apq = D[p, q]; if (Math.Abs(apq) < 1e-15) continue;
                    double app = D[p, p], aqq = D[q, q];
                    double theta = (aqq - app) / (2.0 * apq);
                    double t = (theta >= 0 ? 1.0 : -1.0) / (Math.Abs(theta) + Math.Sqrt(theta * theta + 1.0));
                    double c = 1.0 / Math.Sqrt(1.0 + t * t);
                    double s = t * c;
                    D[p, p] = app - t * apq;
                    D[q, q] = aqq + t * apq;
                    D[p, q] = 0.0; D[q, p] = 0.0;
                    for (int k = 0; k < n; k++)
                    {
                        if (k == p || k == q) continue;
                        double dip = D[k, p], diq = D[k, q];
                        D[k, p] = c * dip - s * diq; D[p, k] = D[k, p];
                        D[k, q] = s * dip + c * diq; D[q, k] = D[k, q];
                    }
                    for (int k = 0; k < n; k++)
                    {
                        double gip = Q[k, p], giq = Q[k, q];
                        Q[k, p] = c * gip - s * giq;
                        Q[k, q] = s * gip + c * giq;
                    }
                }
            }
            for (int i = 0; i < n; i++) lam[i] = D[i, i];
            // sort ascending
            int[] idx = new int[n]; for (int i = 0; i < n; i++) idx[i] = i;
            Array.Sort(lam, idx);
            double[,] Qs = new double[n, n];
            for (int j = 0; j < n; j++) for (int i = 0; i < n; i++) Qs[i, j] = Q[i, idx[j]];
            return Qs;
        }
    }
}
