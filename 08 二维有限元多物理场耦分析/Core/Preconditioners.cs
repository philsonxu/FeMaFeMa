using System;

namespace MultiPhysicsFEM2D.Core
{
    public interface IPreconditioner
    {
        void Build(SparseMatrixCSR A);
        void Apply(DenseVector r, DenseVector z);
    }

    public class IdentityPreconditioner : IPreconditioner
    {
        public void Build(SparseMatrixCSR A) { }
        public void Apply(DenseVector r, DenseVector z) { z.CopyFrom(r); }
    }

    public class JacobiPreconditioner : IPreconditioner
    {
        private double[] _diagInv;
        private int _n;

        public void Build(SparseMatrixCSR A)
        {
            _n = A.Nrows;
            _diagInv = new double[_n];
            for (int i = 0; i < _n; i++)
            {
                double d = 1.0;
                for (int k = A.RowPtr[i]; k < A.RowPtr[i + 1]; k++)
                    if (A.ColIdx[k] == i) { d = A.Values[k]; break; }
                _diagInv[i] = Math.Abs(d) < 1e-30 ? 1.0 : 1.0 / d;
            }
        }

        public void Apply(DenseVector r, DenseVector z)
        {
            for (int i = 0; i < _n; i++) z.Values[i] = _diagInv[i] * r.Values[i];
        }
    }
}
