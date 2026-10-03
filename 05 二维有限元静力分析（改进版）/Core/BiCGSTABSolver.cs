// BiCGSTABSolver.cs - 稳定双共轭梯度法（适用于非对称矩阵，可用于NS动量方程）
using System;

namespace FEM2D.Core
{
    public sealed class BiCGSTABSolver
    {
        public int MaxIterations = 3000;
        public double Tolerance = 1e-9;
        public IPreconditioner Preconditioner;
        public int LastIterations;
        public double LastResidual;
        public event Action<int, double> OnProgress;

        public bool Solve(SparseMatrixCSR A, DenseVector b, DenseVector x)
        {
            int n = A.N;
            IPreconditioner M = Preconditioner ?? new JacobiPreconditioner(A);
            double bnorm = b.Norm2();
            if (bnorm < 1e-30) { x.SetZero(); LastIterations = 0; LastResidual = 0.0; return true; }

            DenseVector r = new DenseVector(n);
            DenseVector r0hat = new DenseVector(n);
            DenseVector v = new DenseVector(n);
            DenseVector p = new DenseVector(n);
            DenseVector s = new DenseVector(n);
            DenseVector t = new DenseVector(n);
            DenseVector p_hat = new DenseVector(n);
            DenseVector s_hat = new DenseVector(n);

            A.Multiply(x, r);
            for (int i = 0; i < n; i++) { r.Values[i] = b.Values[i] - r.Values[i]; r0hat.Values[i] = r.Values[i]; }
            double rho = 1.0, alpha = 1.0, omega = 1.0;
            double res = r.Norm2() / bnorm;

            for (int k = 0; k < MaxIterations; k++)
            {
                double rho_new = r0hat.Dot(r);
                if (Math.Abs(rho_new) < 1e-30) break;
                double beta = (rho_new / rho) * (alpha / omega);
                rho = rho_new;
                for (int i = 0; i < n; i++)
                    p.Values[i] = r.Values[i] + beta * (p.Values[i] - omega * v.Values[i]);
                M.Apply(p, p_hat);
                A.Multiply(p_hat, v);
                double vdotr = r0hat.Dot(v);
                if (Math.Abs(vdotr) < 1e-30) break;
                alpha = rho / vdotr;
                // s = r - alpha*v
                for (int i = 0; i < n; i++) s.Values[i] = r.Values[i] - alpha * v.Values[i];
                res = s.Norm2() / bnorm;
                OnProgress?.Invoke(k, res);
                if (res < Tolerance)
                {
                    x.Axpy(alpha, p_hat);
                    LastIterations = k + 1; LastResidual = res; return true;
                }
                M.Apply(s, s_hat);
                A.Multiply(s_hat, t);
                double tt = t.Dot(t);
                if (tt < 1e-30) { omega = 0.0; }
                else { omega = t.Dot(s) / tt; }
                x.Axpy(alpha, p_hat);
                x.Axpy(omega, s_hat);
                // r = s - omega*t
                for (int i = 0; i < n; i++) r.Values[i] = s.Values[i] - omega * t.Values[i];
                res = r.Norm2() / bnorm;
                OnProgress?.Invoke(k, res);
                LastIterations = k + 1; LastResidual = res;
                if (res < Tolerance) return true;
                if (Math.Abs(omega) < 1e-30) break;
            }
            return res < Tolerance;
        }
    }
}
