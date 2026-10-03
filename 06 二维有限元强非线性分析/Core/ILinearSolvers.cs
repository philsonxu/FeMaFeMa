using System;

namespace NonlinearFEM2D.Core
{
    public interface ILinearSolver
    {
        int Solve(SparseMatrixCSR A, DenseVector b, DenseVector x, IPreconditioner M);
        double Tolerance { get; set; }
        int MaxIterations { get; set; }
        int FinalIterations { get; }
        double FinalResidual { get; }
        Action<int, double> OnIteration { get; set; }
    }

    public interface IPreconditioner
    {
        void Build(SparseMatrixCSR A);
        void Apply(DenseVector r, DenseVector z);
    }

    public sealed class IdentityPreconditioner : IPreconditioner
    {
        public void Build(SparseMatrixCSR A) { }
        public void Apply(DenseVector r, DenseVector z) { z.CopyFrom(r); }
    }

    public sealed class JacobiPreconditioner : IPreconditioner
    {
        private double[] _diag;
        private int _n;
        public void Build(SparseMatrixCSR A)
        {
            this._n = A.Rows;
            this._diag = A.GetDiagonal();
            for (int i = 0; i < this._n; i++)
            {
                if (Math.Abs(this._diag[i]) < 1.0e-30) this._diag[i] = 1.0;
            }
        }
        public void Apply(DenseVector r, DenseVector z)
        {
            for (int i = 0; i < this._n; i++) z[i] = r[i] / this._diag[i];
        }
    }
}
