// CGSolver.cs - 预处理共轭梯度法（对称正定系统）
using System;

namespace FEM2D.Core
{
    public sealed class CGSolver
    {
        public int MaxIterations = 5000;
        public double Tolerance = 1e-10;
        public IPreconditioner Preconditioner;
        public int LastIterations;
        public double LastResidual;
        public event Action<int, double> OnProgress;
        public Action<string> Log;

        public bool Solve(SparseMatrixCSR A, DenseVector b, DenseVector x)
        {
            int n = A.N;
            IPreconditioner M = Preconditioner ?? new JacobiPreconditioner(A);
            DenseVector r = new DenseVector(n);
            DenseVector p = new DenseVector(n);
            DenseVector Ap = new DenseVector(n);
            DenseVector z = new DenseVector(n);
            A.Multiply(x, Ap);
            for (int i = 0; i < n; i++) r.Values[i] = b.Values[i] - Ap.Values[i];
            double bnorm = b.Norm2();
            if (bnorm < 1e-30) { x.SetZero(); LastIterations = 0; LastResidual = 0.0; return true; }
            M.Apply(r, z);
            double rz = r.Dot(z), rzold = 0.0;
            p.CopyFrom(z);
            double res = Math.Sqrt(r.Dot(r)) / bnorm;

            for (int it = 0; it < MaxIterations; it++)
            {
                A.Multiply(p, Ap);
                double pAp = p.Dot(Ap);
                if (Math.Abs(pAp) < 1e-30) break;
                double alpha = rz / pAp;
                x.Axpy(alpha, p);
                r.Axpy(-alpha, Ap);
                M.Apply(r, z);
                rzold = rz;
                rz = r.Dot(z);
                double rnorm = Math.Sqrt(Math.Abs(r.Dot(r)));
                res = rnorm / bnorm;
                OnProgress?.Invoke(it, res);
                if (Log != null && it % 50 == 0)
                    Log(string.Format("CG iter={0}, relRes={1:E3}", it, res));
                if (res < Tolerance) { LastIterations = it + 1; LastResidual = res; return true; }
                double beta = rz / rzold;
                for (int i = 0; i < n; i++) p.Values[i] = z.Values[i] + beta * p.Values[i];
            }
            LastIterations = MaxIterations; LastResidual = res;
            return res < Tolerance;
        }
    }
}
