namespace ModalFEM2D.Core
{
    /// <summary>Preconditioner interface: applies z = M^{-1} r.</summary>
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
        private readonly SparseMatrixCSR _mat;
        public JacobiPreconditioner(SparseMatrixCSR mat) { _mat = mat; }
        public void Apply(DenseVector r, DenseVector z)
        {
            int n = _mat.NumRows;
            double[] di = _mat.DiagInv;
            for (int i = 0; i < n; i++) z.Values[i] = r.Values[i] * di[i];
        }
    }
}
