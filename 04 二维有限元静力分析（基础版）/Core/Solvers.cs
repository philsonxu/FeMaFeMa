using System;

namespace FEM2D.Core
{
    /// <summary>
    /// 重启 GMRES(m) 求解非对称线性系统 A x = b
    /// 使用 Givens 旋转高效最小二乘求解，支持对角预条件。
    /// </summary>
    public sealed class GMRESSolver
    {
        public int MaxRestart { get; set; } = 200;
        public int KrylovDim { get; set; } = 50;
        public double Tolerance { get; set; } = 1e-8;
        public int PrintLevel { get; set; } = 0;
        public int FinalIterations { get; private set; }
        public double FinalResidual { get; private set; }

        public bool Solve(SparseMatrixCSR A, DenseVector b, DenseVector x, double[] diagInv = null)
        {
            int n = A.N;
            int m = KrylovDim;
            if (diagInv == null) diagInv = new double[n];
            for (int i = 0; i < n; i++) if (diagInv[i] == 0.0) diagInv[i] = 1.0;

            DenseVector r = new DenseVector(n);
            DenseVector w = new DenseVector(n);

            // r = M^{-1}(b - Ax)
            A.MatVec(x, w);
            for (int i = 0; i < n; i++) r[i] = diagInv[i] * (b[i] - w[i]);

            double beta = r.Norm2();
            double bnorm = b.Norm2();
            if (bnorm == 0.0) bnorm = 1.0;

            if (beta / bnorm < Tolerance) { FinalResidual = beta / bnorm; FinalIterations = 0; return true; }

            DenseVector[] V = new DenseVector[m + 1];
            for (int i = 0; i <= m; i++) V[i] = new DenseVector(n);
            double[,] H = new double[m + 1, m];
            double[] cs = new double[m];
            double[] sn = new double[m];
            double[] g = new double[m + 1];

            for (int restart = 0; restart < MaxRestart; restart++)
            {
                for (int i = 0; i < n; i++) V[0][i] = r[i] / beta;
                Array.Clear(g, 0, g.Length);
                g[0] = beta;

                int j;
                for (j = 0; j < m; j++)
                {
                    // w = M^{-1} A V[j]
                    A.MatVec(V[j], w);
                    for (int i = 0; i < n; i++) w[i] *= diagInv[i];

                    // Arnoldi 修正 Gram-Schmidt
                    for (int i = 0; i <= j; i++)
                    {
                        H[i, j] = w.Dot(V[i]);
                        w.AddScaled(V[i], -H[i, j]);
                    }
                    H[j + 1, j] = w.Norm2();

                    if (Math.Abs(H[j + 1, j]) < 1e-14)
                    {
                        // Happy breakdown
                        for (int i = 0; i <= j; i++) V[j + 1][i] = 0.0;
                    }
                    else
                    {
                        for (int i = 0; i < n; i++) V[j + 1][i] = w[i] / H[j + 1, j];
                    }

                    // Apply previous Givens rotations
                    for (int i = 0; i < j; i++)
                    {
                        double h1 = H[i, j];
                        double h2 = H[i + 1, j];
                        H[i, j] = cs[i] * h1 + sn[i] * h2;
                        H[i + 1, j] = -sn[i] * h1 + cs[i] * h2;
                    }

                    // Compute new Givens rotation
                    double h_jj = H[j, j];
                    double h_j1j = H[j + 1, j];
                    double denom = Math.Sqrt(h_jj * h_jj + h_j1j * h_j1j);
                    if (denom < 1e-30) { cs[j] = 1.0; sn[j] = 0.0; }
                    else { cs[j] = h_jj / denom; sn[j] = h_j1j / denom; }
                    H[j, j] = cs[j] * h_jj + sn[j] * h_j1j;
                    H[j + 1, j] = 0.0;

                    double g_j = g[j];
                    g[j] = cs[j] * g_j;
                    g[j + 1] = -sn[j] * g_j;

                    double res = Math.Abs(g[j + 1]);
                    FinalIterations = restart * m + j + 1;
                    FinalResidual = res / bnorm;
                    if (PrintLevel > 0 && (FinalIterations % 10 == 0))
                    {
                        Console.WriteLine($"GMRES iter {FinalIterations}, relres = {FinalResidual:E4}");
                    }
                    if (FinalResidual < Tolerance)
                    {
                        UpdateSolution(x, H, V, g, j, n);
                        return true;
                    }
                }
                j = m;
                UpdateSolution(x, H, V, g, j - 1, n);

                A.MatVec(x, w);
                for (int i = 0; i < n; i++) r[i] = diagInv[i] * (b[i] - w[i]);
                beta = r.Norm2();
                FinalResidual = beta / bnorm;
                if (FinalResidual < Tolerance) return true;
            }
            return false;
        }

        private static void UpdateSolution(DenseVector x, double[,] H, DenseVector[] V, double[] g, int j, int n)
        {
            // 回代
            double[] y = new double[j + 1];
            for (int i = j; i >= 0; i--)
            {
                double s = g[i];
                for (int k = i + 1; k <= j; k++) s -= H[i, k] * y[k];
                y[i] = s / H[i, i];
            }
            for (int i = 0; i <= j; i++)
            {
                x.AddScaled(V[i], y[i]);
            }
        }
    }

    /// <summary>
    /// 共轭梯度法，用于对称正定系统（结构刚度）。
    /// </summary>
    public sealed class CGSolver
    {
        public int MaxIterations { get; set; } = 5000;
        public double Tolerance { get; set; } = 1e-8;
        public int FinalIterations { get; private set; }
        public double FinalResidual { get; private set; }

        public bool Solve(SparseMatrixCSR A, DenseVector b, DenseVector x, double[] diagInv = null)
        {
            int n = A.N;
            DenseVector r = new DenseVector(n);
            DenseVector p = new DenseVector(n);
            DenseVector Ap = new DenseVector(n);

            A.MatVec(x, Ap);
            for (int i = 0; i < n; i++)
            {
                r[i] = b[i] - Ap[i];
                if (diagInv != null) r[i] *= diagInv[i];
                p[i] = r[i];
            }
            double rsold = r.Dot(r);
            double bnorm = b.Norm2();
            if (bnorm == 0.0) bnorm = 1.0;

            for (int it = 0; it < MaxIterations; it++)
            {
                A.MatVec(p, Ap);
                if (diagInv != null) for (int i = 0; i < n; i++) Ap[i] *= diagInv[i];
                double pAp = p.Dot(Ap);
                if (Math.Abs(pAp) < 1e-30) return false;
                double alpha = rsold / pAp;
                x.AddScaled(p, alpha);
                r.AddScaled(Ap, -alpha);
                double rsnew = r.Dot(r);
                FinalResidual = Math.Sqrt(rsnew) / bnorm;
                FinalIterations = it + 1;
                if (FinalResidual < Tolerance) return true;
                double beta = rsnew / rsold;
                for (int i = 0; i < n; i++) p[i] = r[i] + beta * p[i];
                rsold = rsnew;
            }
            return FinalResidual < Tolerance;
        }
    }
}
