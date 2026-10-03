using System;

namespace FEM2D.Core
{
    public sealed class GMRESSolver
    {
        public int MaxIterations = 2000;
        public int Restart = 80;
        public double Tolerance = 1e-10;
        public IPreconditioner Preconditioner;
        public int LastIterations;
        public double LastResidual;
        public event Action<int, double> OnProgress;

        public bool Solve(SparseMatrixCSR A, DenseVector b, DenseVector x)
        {
            int n = A.N;
            int m = Math.Min(Restart, n);
            IPreconditioner M = Preconditioner ?? new JacobiPreconditioner(A);
            double bnorm = b.Norm2();
            if (bnorm < 1e-30) { x.SetZero(); LastIterations = 0; LastResidual = 0.0; return true; }
            DenseVector r = new DenseVector(n);
            DenseVector w = new DenseVector(n);
            DenseVector z = new DenseVector(n);
            DenseVector[] V = new DenseVector[m + 1];
            for (int k = 0; k <= m; k++) V[k] = new DenseVector(n);
            double[,] H = new double[m + 1, m];
            double[] cs = new double[m];
            double[] sn = new double[m];
            double[] g = new double[m + 1];

            int totalIter = 0; double res = 1.0;
            while (totalIter < MaxIterations)
            {
                A.Multiply(x, w);
                for (int i = 0; i < n; i++) r.Values[i] = b.Values[i] - w.Values[i];
                M.Apply(r, z);
                double beta = z.Norm2();
                res = beta / bnorm;
                if (res < Tolerance) { LastIterations = totalIter; LastResidual = res; return true; }
                for (int i = 0; i < n; i++) V[0].Values[i] = z.Values[i] / beta;
                Array.Clear(g, 0, g.Length); g[0] = beta;
                int j = 0;
                for (; j < m && totalIter < MaxIterations; j++, totalIter++)
                {
                    A.Multiply(V[j], w);
                    M.Apply(w, z);
                    Array.Copy(z.Values, w.Values, n);
                    for (int i = 0; i <= j; i++)
                    { H[i, j] = w.Dot(V[i]); w.Axpy(-H[i, j], V[i]); }
                    H[j + 1, j] = w.Norm2();
                    if (H[j + 1, j] < 1e-30) V[j + 1].SetZero();
                    else for (int i = 0; i < n; i++) V[j + 1].Values[i] = w.Values[i] / H[j + 1, j];
                    for (int i = 0; i < j; i++)
                    {
                        double v1 = H[i, j], v2 = H[i + 1, j];
                        H[i, j] = cs[i] * v1 + sn[i] * v2;
                        H[i + 1, j] = -sn[i] * v1 + cs[i] * v2;
                    }
                    double h1 = H[j, j], h2 = H[j + 1, j];
                    double den = Math.Sqrt(h1 * h1 + h2 * h2);
                    if (den < 1e-30) { cs[j] = 1.0; sn[j] = 0.0; }
                    else { cs[j] = h1 / den; sn[j] = h2 / den; }
                    H[j, j] = cs[j] * h1 + sn[j] * h2; H[j + 1, j] = 0.0;
                    double g0 = g[j];
                    g[j] = cs[j] * g0 + sn[j] * g[j + 1];
                    g[j + 1] = -sn[j] * g0 + cs[j] * g[j + 1];
                    res = Math.Abs(g[j + 1]) / bnorm;
                    OnProgress?.Invoke(totalIter, res);
                    if (res < Tolerance)
                    {
                        BacksubAndUpdate(j, H, g, V, n, M, x, z, w);
                        LastIterations = totalIter + 1; LastResidual = res; return true;
                    }
                }
                int jj = Math.Min(j, m - 1);
                if (jj <= 0) break;
                BacksubAndUpdate(jj, H, g, V, n, M, x, z, w);
            }
            LastIterations = totalIter; LastResidual = res;
            return res < Tolerance;
        }

        private static void BacksubAndUpdate(int j, double[,] H, double[] g, DenseVector[] V, int n, IPreconditioner M, DenseVector x, DenseVector z, DenseVector w)
        {
            // 左预条件 GMRES：Krylov 基 V 已在 M^{-1}A 的作用下构造，
            // 回代得到 y 后直接 x += V*y，不再额外应用预条件。
            double[] y = new double[j + 1];
            for (int ii = j; ii >= 0; ii--)
            {
                double s = g[ii];
                for (int kk = ii + 1; kk <= j; kk++) s -= H[ii, kk] * y[kk];
                y[ii] = s / H[ii, ii];
            }
            for (int kk = 0; kk <= j; kk++) x.Axpy(y[kk], V[kk]);
        }
    }
}
