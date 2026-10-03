// IPreconditioner.cs - 预条件子统一接口 + Jacobi 预条件
using System;

namespace FEM2D.Core
{
    public interface IPreconditioner
    {
        void Apply(DenseVector r, DenseVector z);
    }

    public sealed class IdentityPreconditioner : IPreconditioner
    {
        public void Apply(DenseVector r, DenseVector z) { z.CopyFrom(r); }
    }

    public sealed class JacobiPreconditioner : IPreconditioner
    {
        private readonly double[] _diagInv;
        private readonly int _n;

        public JacobiPreconditioner(SparseMatrixCSR A)
        {
            _n = A.N;
            _diagInv = new double[_n];
            Array.Copy(A.DiagInv, _diagInv, _n);
        }

        public JacobiPreconditioner(double[] diagInv)
        {
            _n = diagInv.Length;
            _diagInv = new double[_n];
            Array.Copy(diagInv, _diagInv, _n);
        }

        public void Apply(DenseVector r, DenseVector z)
        {
            int n = _n;
            for (int i = 0; i < n; i++) z.Values[i] = r.Values[i] * _diagInv[i];
        }
    }
}
