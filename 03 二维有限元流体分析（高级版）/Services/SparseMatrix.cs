using System;
using System.Collections.Generic;

namespace Fem2DFluid.Services
{
    /// <summary>
    /// Coordinate (triplet) sparse matrix used during assembly. Call Add(i,j,v)
    /// repeatedly; duplicates are summed when converting to CSR.
    /// </summary>
    public class CooMatrix
    {
        private List<int> _rows;
        private List<int> _cols;
        private List<double> _vals;
        private int _n;

        public CooMatrix(int n)
        {
            _n = n;
            _rows = new List<int>();
            _cols = new List<int>();
            _vals = new List<double>();
        }

        public int N
        {
            get { return _n; }
        }

        public List<int> Rows
        {
            get { return _rows; }
        }

        public List<int> Cols
        {
            get { return _cols; }
        }

        public List<double> Vals
        {
            get { return _vals; }
        }

        public void Add(int i, int j, double v)
        {
            if (i < 0 || j < 0 || i >= _n || j >= _n)
            {
                return;
            }
            _rows.Add(i);
            _cols.Add(j);
            _vals.Add(v);
        }

        public void AddBlock3x3(int[] idx, double[,] ke)
        {
            for (int a = 0; a < 3; a++)
            {
                for (int b = 0; b < 3; b++)
                {
                    Add(idx[a], idx[b], ke[a, b]);
                }
            }
        }

        public void AddBlock6x6(int[] idx, double[,] ke)
        {
            for (int a = 0; a < 6; a++)
            {
                for (int b = 0; b < 6; b++)
                {
                    Add(idx[a], idx[b], ke[a, b]);
                }
            }
        }
    }

    /// <summary>
    /// Compressed Sparse Row matrix with ILU0 preconditioner and restarted GMRES(m)
    /// solver suitable for symmetric positive definite (potential) and mildly
    /// non-symmetric (convection-diffusion / NS momentum) systems.
    /// </summary>
    public class SparseMatrix
    {
        private int _n;
        private int[] _rowPtr;
        private int[] _cols;
        private double[] _vals;

        // ILU0 factors (same sparsity pattern, stored as separate arrays parallel to vals)
        private double[] _luVals;

        public SparseMatrix(int n, int[] rowPtr, int[] cols, double[] vals)
        {
            _n = n;
            _rowPtr = rowPtr;
            _cols = cols;
            _vals = vals;
            _luVals = null;
        }

        public int N
        {
            get { return _n; }
        }

        public int[] RowPtr
        {
            get { return _rowPtr; }
        }

        public int[] Cols
        {
            get { return _cols; }
        }

        public double[] Vals
        {
            get { return _vals; }
        }

        /// <summary>
        /// Build a CSR matrix from COO triplets. Duplicate (i,j) entries are summed.
        /// Rows and columns must be within [0,n-1].
        /// </summary>
        public static SparseMatrix FromCoo(CooMatrix coo)
        {
            int n = coo.N;
            List<int> cr = coo.Rows;
            List<int> cc = coo.Cols;
            List<double> cv = coo.Vals;
            int nnz = cr.Count;

            // Count entries per row
            int[] rowCount = new int[n + 1];
            for (int k = 0; k < nnz; k++)
            {
                rowCount[cr[k] + 1]++;
            }
            // Prefix sum for rowPtr (provisional)
            int[] rowPtr = new int[n + 1];
            for (int i = 0; i < n; i++)
            {
                rowPtr[i + 1] = rowPtr[i] + rowCount[i + 1];
            }
            // Bucket sort by row
            int[] cols = new int[nnz];
            double[] vals = new double[nnz];
            int[] cursor = new int[n];
            Array.Copy(rowPtr, cursor, n);
            for (int k = 0; k < nnz; k++)
            {
                int r = cr[k];
                int pos = cursor[r];
                cols[pos] = cc[k];
                vals[pos] = cv[k];
                cursor[r]++;
            }
            // Within each row: sort by column, sum duplicates
            List<int> fcols = new List<int>(nnz);
            List<double> fvals = new List<double>(nnz);
            int[] frowPtr = new int[n + 1];
            for (int i = 0; i < n; i++)
            {
                frowPtr[i] = fcols.Count;
                int start = rowPtr[i];
                int end = rowPtr[i + 1];
                if (end <= start) { continue; }
                // Extract and sort by column
                int cnt = end - start;
                int[] lcols = new int[cnt];
                double[] lvals = new double[cnt];
                Array.Copy(cols, start, lcols, 0, cnt);
                Array.Copy(vals, start, lvals, 0, cnt);
                Array.Sort(lcols, lvals);
                // Merge duplicates
                int curCol = lcols[0];
                double curVal = lvals[0];
                for (int k = 1; k < cnt; k++)
                {
                    if (lcols[k] == curCol)
                    {
                        curVal += lvals[k];
                    }
                    else
                    {
                        fcols.Add(curCol);
                        fvals.Add(curVal);
                        curCol = lcols[k];
                        curVal = lvals[k];
                    }
                }
                fcols.Add(curCol);
                fvals.Add(curVal);
            }
            frowPtr[n] = fcols.Count;
            return new SparseMatrix(n, frowPtr, fcols.ToArray(), fvals.ToArray());
        }

        /// <summary>
        /// y = A * x (CSR matrix-vector product).
        /// </summary>
        public void Mult(double[] x, double[] y)
        {
            for (int i = 0; i < _n; i++)
            {
                double sum = 0.0;
                int end = _rowPtr[i + 1];
                for (int k = _rowPtr[i]; k < end; k++)
                {
                    sum += _vals[k] * x[_cols[k]];
                }
                y[i] = sum;
            }
        }

        /// <summary>
        /// In-place ILU(0) factorization: L and U share the same CSR sparsity pattern
        /// as A (zero fill-in). After calling this method, <see cref="ApplyILU0"/> can
        /// be used as a preconditioner.
        /// </summary>
        public void BuildILU0()
        {
            _luVals = new double[_vals.Length];
            Array.Copy(_vals, _luVals, _vals.Length);
            // Helper: find column index within a row
            for (int i = 0; i < _n; i++)
            {
                int p;
                int rowIStart = _rowPtr[i];
                int rowIEnd = _rowPtr[i + 1];
                for (p = rowIStart; p < rowIEnd; p++)
                {
                    int jp = _cols[p];
                    if (jp >= i) { break; }
                    // Find diagonal position for row jp
                    int rowJStart = _rowPtr[jp];
                    int rowJEnd = _rowPtr[jp + 1];
                    int kk = -1;
                    for (int q = rowJStart; q < rowJEnd; q++)
                    {
                        if (_cols[q] == jp) { kk = q; break; }
                    }
                    if (kk < 0) { continue; }
                    double diagJ = _luVals[kk];
                    if (Math.Abs(diagJ) < 1.0e-30) { continue; }
                    _luVals[p] /= diagJ;
                    double aik = _luVals[p];
                    // Multiply row jp entries beyond diagonal, subtract from row i
                    for (int q = kk + 1; q < rowJEnd; q++)
                    {
                        int jw = _cols[q];
                        // Find position jw in row i (binary search)
                        int pos = FindCol(i, jw);
                        if (pos >= 0)
                        {
                            _luVals[pos] -= aik * _luVals[q];
                        }
                    }
                }
            }
        }

        private int FindCol(int row, int col)
        {
            int lo = _rowPtr[row];
            int hi = _rowPtr[row + 1] - 1;
            while (lo <= hi)
            {
                int mid = (lo + hi) >> 1;
                int c = _cols[mid];
                if (c == col) { return mid; }
                if (c < col) { lo = mid + 1; } else { hi = mid - 1; }
            }
            return -1;
        }

        /// <summary>
        /// Solve (LU) z = r using the stored ILU0 factors (forward + backward).
        /// </summary>
        public void ApplyILU0(double[] r, double[] z)
        {
            if (_luVals == null)
            {
                Array.Copy(r, z, _n);
                return;
            }
            // Forward: L y = r (L unit lower triangular with diagonal=1)
            double[] y = new double[_n];
            for (int i = 0; i < _n; i++)
            {
                double sum = r[i];
                int rowStart = _rowPtr[i];
                int rowEnd = _rowPtr[i + 1];
                int diag = -1;
                for (int k = rowStart; k < rowEnd; k++)
                {
                    int c = _cols[k];
                    if (c < i)
                    {
                        sum -= _luVals[k] * y[c];
                    }
                    else if (c == i)
                    {
                        diag = k;
                        break;
                    }
                }
                y[i] = sum;
                if (diag >= 0)
                {
                    // Nothing: L has unit diagonal
                }
            }
            // Backward: U z = y
            for (int i = _n - 1; i >= 0; i--)
            {
                double sum = y[i];
                double diag = 1.0;
                int rowStart = _rowPtr[i];
                int rowEnd = _rowPtr[i + 1];
                int diagPos = -1;
                for (int k = rowStart; k < rowEnd; k++)
                {
                    if (_cols[k] == i) { diagPos = k; diag = _luVals[k]; break; }
                }
                if (diagPos >= 0)
                {
                    for (int k = diagPos + 1; k < rowEnd; k++)
                    {
                        sum -= _luVals[k] * z[_cols[k]];
                    }
                }
                z[i] = (Math.Abs(diag) > 1.0e-30) ? sum / diag : sum;
            }
        }

        /// <summary>
        /// Restarted GMRES(m), m=30. Solves A x = b with initial guess x.
        /// Returns number of iterations performed; convergence achieved if
        /// return < maxIter and final residual < tol * ||b||.
        /// GMRES（Generalized Minimal Residual，广义最小残差法）是一种用于求解大型稀疏非对称线性方程组的迭代算法‌。
        /// 该方法由 Yousef Saad 和 Martin H. Schultz 于 1986 年提出。‌‌
        /// </summary>
        public int GMRES(double[] b, double[] x, double tol, int maxIter, bool useILU0)
        {
            int m = 30;
            int n = _n;
            double[] r = new double[n];
            double[] w = new double[n];
            double[,] V = new double[m + 1, n];
            double[,] H = new double[m + 1, m];
            double[] cs = new double[m];
            double[] sn = new double[m];
            double[] g = new double[m + 1];
            double[] z = new double[n];
            Mult(x, w);
            for (int i = 0; i < n; i++) { r[i] = b[i] - w[i]; }
            if (useILU0 && _luVals != null)
            {
                ApplyILU0(r, z);
            }
            else
            {
                Array.Copy(r, z, n);
            }
            double bnorm = VecNorm(b);
            if (bnorm < 1.0e-30) { bnorm = 1.0; }
            double beta = VecNorm(z);
            if (beta < tol * bnorm) { return 0; }
            int totalIter = 0;
            for (int outer = 0; outer < maxIter; outer += m)
            {
                for (int i = 0; i < n; i++) { V[0, i] = z[i] / beta; }
                for (int jj = 0; jj < m + 1; jj++) { g[jj] = 0.0; }
                g[0] = beta;
                int j;
                for (j = 0; j < m; j++)
                {
                    // w = A * V[j]
                    double[] vj = new double[n];
                    for (int i = 0; i < n; i++) { vj[i] = V[j, i]; }
                    Mult(vj, w);
                    if (useILU0 && _luVals != null)
                    {
                        ApplyILU0(w, z);
                    }
                    else { Array.Copy(w, z, n); }
                    // Arnoldi
                    for (int i = 0; i <= j; i++)
                    {
                        double[] vi = new double[n];
                        for (int k = 0; k < n; k++) { vi[k] = V[i, k]; }
                        H[i, j] = VecDot(vi, z);
                        for (int k = 0; k < n; k++) { z[k] -= H[i, j] * vi[k]; }
                    }
                    H[j + 1, j] = VecNorm(z);
                    if (H[j + 1, j] < 1.0e-30) { break; }
                    for (int k = 0; k < n; k++) { V[j + 1, k] = z[k] / H[j + 1, j]; }
                    // Apply previous Givens rotations
                    for (int i = 0; i < j; i++)
                    {
                        double tmp = cs[i] * H[i, j] + sn[i] * H[i + 1, j];
                        H[i + 1, j] = -sn[i] * H[i, j] + cs[i] * H[i + 1, j];
                        H[i, j] = tmp;
                    }
                    // Compute new Givens rotation
                    double denom = Math.Sqrt(H[j, j] * H[j, j] + H[j + 1, j] * H[j + 1, j]);
                    if (denom < 1.0e-30) { cs[j] = 1.0; sn[j] = 0.0; }
                    else { cs[j] = H[j, j] / denom; sn[j] = H[j + 1, j] / denom; }
                    H[j, j] = cs[j] * H[j, j] + sn[j] * H[j + 1, j];
                    H[j + 1, j] = 0.0;
                    double gtmp = cs[j] * g[j];
                    g[j + 1] = -sn[j] * g[j];
                    g[j] = gtmp;
                    totalIter++;
                    double res = Math.Abs(g[j + 1]);
                    if (res < tol * bnorm) { j++; break; }
                }
                // Back-solve upper triangular H*y = g (up to j-1)
                int kk = j - 1;
                double[] y = new double[j];
                for (int i = kk; i >= 0; i--)
                {
                    double s = g[i];
                    for (int p = i + 1; p <= kk; p++) { s -= H[i, p] * y[p]; }
                    y[i] = s / H[i, i];
                }
                // Update x
                for (int i = 0; i <= kk; i++)
                {
                    for (int k = 0; k < n; k++) { x[k] += y[i] * V[i, k]; }
                }
                // Recompute residual
                Mult(x, w);
                for (int i = 0; i < n; i++) { r[i] = b[i] - w[i]; }
                if (useILU0 && _luVals != null) { ApplyILU0(r, z); }
                else { Array.Copy(r, z, n); }
                beta = VecNorm(z);
                if (beta < tol * bnorm) { return totalIter; }
            }
            return totalIter;
        }

        /// <summary>
        /// Apply Dirichlet BC: set row i to e_i, rhs[i] = value; symmetrically for column.
        /// </summary>
        public void ApplyDirichlet(int i, double value, double[] rhs)
        {
            int rs = _rowPtr[i];
            int re = _rowPtr[i + 1];
            for (int k = rs; k < re; k++)
            {
                int c = _cols[k];
                if (c == i) { _vals[k] = 1.0; }
                else
                {
                    rhs[c] -= _vals[k] * value;
                    _vals[k] = 0.0;
                }
            }
            rhs[i] = value;
            // Now zero column i in other rows (diagonal is kept, already handled)
            for (int r = 0; r < _n; r++)
            {
                if (r == i) { continue; }
                int pos = FindCol(r, i);
                if (pos >= 0) { _vals[pos] = 0.0; }
            }
        }

        public static double VecDot(double[] a, double[] b)
        {
            double s = 0.0;
            for (int i = 0; i < a.Length; i++) { s += a[i] * b[i]; }
            return s;
        }

        public static double VecNorm(double[] a)
        {
            return Math.Sqrt(VecDot(a, a));
        }
    }
}
