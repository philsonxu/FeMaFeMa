using System;

namespace NonlinearFEM2D.Core
{
    public sealed class CGSolver : ILinearSolver
    {
        public double Tolerance { get; set; } = 1.0e-8;
        public int MaxIterations { get; set; } = 5000;
        public int FinalIterations { get; private set; }
        public double FinalResidual { get; private set; }
        public Action<int, double> OnIteration { get; set; }

        public int Solve(SparseMatrixCSR A, DenseVector b, DenseVector x, IPreconditioner M)
        {
            int n = A.Rows;
            DenseVector r = new DenseVector(n);
            DenseVector z = new DenseVector(n);
            DenseVector p = new DenseVector(n);
            DenseVector Ap = new DenseVector(n);
            A.MatVec(x, r);
            for (int i = 0; i < n; i++) r[i] = b[i] - r[i];
            double bnorm = b.Norm2();
            if (bnorm < 1.0e-30) bnorm = 1.0;
            M.Apply(r, z);
            p.CopyFrom(z);
            double rz = r.Dot(z);
            int it = 0;
            for (; it < this.MaxIterations; it++)
            {
                A.MatVec(p, Ap);
                double pAp = p.Dot(Ap);
                if (Math.Abs(pAp) < 1.0e-30) break;
                double alpha = rz / pAp;
                x.Axpy(alpha, p);
                r.Axpy(-alpha, Ap);
                double rn = r.Norm2() / bnorm;
                if (this.OnIteration != null) this.OnIteration(it + 1, rn);
                if (rn < this.Tolerance) { it++; break; }
                M.Apply(r, z);
                double rznew = r.Dot(z);
                double beta = rznew / rz;
                for (int i = 0; i < n; i++) p[i] = z[i] + beta * p[i];
                rz = rznew;
            }
            this.FinalIterations = it;
            this.FinalResidual = r.Norm2() / bnorm;
            return it;
        }
    }

    public sealed class GMRESSolver : ILinearSolver
    {
        public int Restart { get; set; } = 80;
        public double Tolerance { get; set; } = 1.0e-8;
        public int MaxIterations { get; set; } = 5000;
        public int FinalIterations { get; private set; }
        public double FinalResidual { get; private set; }
        public Action<int, double> OnIteration { get; set; }

        public int Solve(SparseMatrixCSR A, DenseVector b, DenseVector x, IPreconditioner M)
        {
            int n = A.Rows;
            int m = this.Restart;
            double bnorm = b.Norm2();
            if (bnorm < 1.0e-30) bnorm = 1.0;
            DenseVector r = new DenseVector(n);
            DenseVector w = new DenseVector(n);
            int totalIt = 0;
            double[,] H = new double[m + 1, m];
            DenseVector[] V = new DenseVector[m + 1];
            for (int k = 0; k <= m; k++) V[k] = new DenseVector(n);
            double[] cs = new double[m];
            double[] sn = new double[m];
            double[] g = new double[m + 1];
            int outer;
            for (outer = 0; outer < 20; outer++)
            {
                A.MatVec(x, r);
                for (int i = 0; i < n; i++) r[i] = b[i] - r[i];
                M.Apply(r, V[0]);
                double beta = V[0].Norm2();
                for (int i = 0; i <= m; i++) g[i] = 0.0;
                g[0] = beta;
                V[0].Scale(1.0 / beta);
                for (int k = 0; k < m; k++)
                    for (int l = 0; l <= m; l++) H[l, k] = 0.0;
                int j;
                for (j = 0; j < m; j++)
                {
                    A.MatVec(V[j], w);
                    M.Apply(w, V[j + 1]);
                    for (int i = 0; i <= j; i++)
                    {
                        H[i, j] = V[j + 1].Dot(V[i]);
                        V[j + 1].Axpy(-H[i, j], V[i]);
                    }
                    H[j + 1, j] = V[j + 1].Norm2();
                    V[j + 1].Scale(1.0 / (H[j + 1, j] + 1.0e-30));
                    for (int i = 0; i < j; i++)
                    {
                        double h1 = H[i, j];
                        double h2 = H[i + 1, j];
                        double c = cs[i], s = sn[i];
                        H[i, j] = c * h1 + s * h2;
                        H[i + 1, j] = -s * h1 + c * h2;
                    }
                    double den = Math.Sqrt(H[j, j] * H[j, j] + H[j + 1, j] * H[j + 1, j]) + 1.0e-30;
                    cs[j] = H[j, j] / den;
                    sn[j] = H[j + 1, j] / den;
                    H[j, j] = cs[j] * H[j, j] + sn[j] * H[j + 1, j];
                    H[j + 1, j] = 0.0;
                    double g0 = g[j];
                    g[j] = cs[j] * g0;
                    g[j + 1] = -sn[j] * g0;
                    totalIt++;
                    double res = Math.Abs(g[j + 1]) / bnorm;
                    if (this.OnIteration != null) this.OnIteration(totalIt, res);
                    if (res < this.Tolerance)
                    {
                        double[] y = new double[j + 1];
                        for (int i = j; i >= 0; i--)
                        {
                            y[i] = g[i];
                            for (int l = i + 1; l <= j; l++) y[i] -= H[i, l] * y[l];
                            y[i] /= H[i, i];
                        }
                        for (int i = 0; i <= j; i++) x.Axpy(y[i], V[i]);
                        this.FinalIterations = totalIt;
                        this.FinalResidual = res;
                        return totalIt;
                    }
                }
                double[] yy = new double[m];
                for (int i = m - 1; i >= 0; i--)
                {
                    yy[i] = g[i];
                    for (int l = i + 1; l < m; l++) yy[i] -= H[i, l] * yy[l];
                    yy[i] /= H[i, i];
                }
                for (int i = 0; i < m; i++) x.Axpy(yy[i], V[i]);
                if (totalIt >= this.MaxIterations) break;
            }
            A.MatVec(x, r);
            for (int i = 0; i < n; i++) r[i] = b[i] - r[i];
            this.FinalIterations = totalIt;
            this.FinalResidual = r.Norm2() / bnorm;
            return totalIt;
        }
    }

    public sealed class BiCGSTABSolver : ILinearSolver
    {
        public double Tolerance { get; set; } = 1.0e-8;
        public int MaxIterations { get; set; } = 5000;
        public int FinalIterations { get; private set; }
        public double FinalResidual { get; private set; }
        public Action<int, double> OnIteration { get; set; }

        public int Solve(SparseMatrixCSR A, DenseVector b, DenseVector x, IPreconditioner M)
        {
            int n = A.Rows;
            DenseVector r = new DenseVector(n);
            DenseVector r0 = new DenseVector(n);
            DenseVector p = new DenseVector(n);
            DenseVector v = new DenseVector(n);
            DenseVector s = new DenseVector(n);
            DenseVector t = new DenseVector(n);
            DenseVector phat = new DenseVector(n);
            DenseVector shat = new DenseVector(n);
            A.MatVec(x, r);
            for (int i = 0; i < n; i++) r[i] = b[i] - r[i];
            r0.CopyFrom(r);
            p.CopyFrom(r);
            double bnorm = b.Norm2();
            if (bnorm < 1.0e-30) bnorm = 1.0;
            double rho = 1.0, alpha = 1.0, omega = 1.0;
            int it = 0;
            for (; it < this.MaxIterations; it++)
            {
                double rhonew = r0.Dot(r);
                if (Math.Abs(rhonew) < 1.0e-30) break;
                double beta = (rhonew / rho) * (alpha / omega);
                for (int i = 0; i < n; i++) p[i] = r[i] + beta * (p[i] - omega * v[i]);
                M.Apply(p, phat);
                A.MatVec(phat, v);
                double tv = r0.Dot(v);
                if (Math.Abs(tv) < 1.0e-30) break;
                alpha = rhonew / tv;
                for (int i = 0; i < n; i++) s[i] = r[i] - alpha * v[i];
                M.Apply(s, shat);
                A.MatVec(shat, t);
                double tt = t.Dot(t);
                omega = (tt < 1.0e-30) ? 0.0 : t.Dot(s) / tt;
                for (int i = 0; i < n; i++) x[i] = x[i] + alpha * phat[i] + omega * shat[i];
                for (int i = 0; i < n; i++) r[i] = s[i] - omega * t[i];
                rho = rhonew;
                double res = r.Norm2() / bnorm;
                if (this.OnIteration != null) this.OnIteration(it + 1, res);
                if (res < this.Tolerance) { it++; break; }
                if (Math.Abs(omega) < 1.0e-30) break;
            }
            this.FinalIterations = it;
            this.FinalResidual = r.Norm2() / bnorm;
            return it;
        }
    }
}
