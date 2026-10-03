using System;
using System.Collections.Generic;
using System.Threading.Tasks;

namespace NonlinearFEM2D.Core
{
    public sealed class AMGPreconditioner : IPreconditioner
    {
        private const int MAX_LEVELS = 8;
        private const int COARSE_MAX = 64;
        private const double THETA = 0.08;
        private const int NU1 = 2;
        private const int NU2 = 2;
        private const double OMEGA = 0.6;

        private sealed class Level
        {
            public SparseMatrixCSR A;
            public SparseMatrixCSR P, R;
            public double[] Diag;
            public double[] DiagInv;
            public double[,] Ainv;
            public int N;
        }

        private readonly List<Level> _levels = new List<Level>();

        public void Build(SparseMatrixCSR A)
        {
            this._levels.Clear();
            SparseMatrixCSR cur = A;
            for (int lvl = 0; lvl < MAX_LEVELS; lvl++)
            {
                Level L = new Level();
                L.A = cur;
                L.N = cur.Rows;
                L.Diag = cur.GetDiagonal();
                L.DiagInv = new double[L.N];
                for (int i = 0; i < L.N; i++)
                    L.DiagInv[i] = (Math.Abs(L.Diag[i]) < 1.0e-30) ? 1.0 : 1.0 / L.Diag[i];
                if (cur.Rows <= COARSE_MAX || lvl == MAX_LEVELS - 1)
                {
                    int n = cur.Rows;
                    double[,] M = new double[n, n];
                    for (int i = 0; i < n; i++)
                        for (int k = cur.RowPtr[i]; k < cur.RowPtr[i + 1]; k++)
                            M[i, cur.ColIdx[k]] = cur.Values[k];
                    L.Ainv = GaussJordanInverse(M, n);
                    this._levels.Add(L);
                    break;
                }
                int[] aggr = Aggregate(cur);
                int nc = 0;
                for (int i = 0; i < cur.Rows; i++) if (aggr[i] >= nc) nc = aggr[i] + 1;
                SparseMatrixCSR Pt = BuildPTent(cur, aggr, nc);
                SparseMatrixCSR Psm = SmoothP(cur, Pt, L.DiagInv);
                SparseMatrixCSR R = Transpose(Psm);
                SparseMatrixCSR Ac = RAP(R, cur, Psm);
                L.P = Psm; L.R = R;
                this._levels.Add(L);
                cur = Ac;
            }
        }

        private static int[] Aggregate(SparseMatrixCSR A)
        {
            int n = A.Rows;
            int[] aggr = new int[n];
            for (int i = 0; i < n; i++) aggr[i] = -1;
            double[] d = A.GetDiagonal();
            int[] order = new int[n];
            for (int i = 0; i < n; i++) order[i] = i;
            Array.Sort(order, (a, b) => Math.Abs(d[b]).CompareTo(Math.Abs(d[a])));
            int cur = 0;
            for (int ii = 0; ii < n; ii++)
            {
                int i = order[ii];
                if (aggr[i] >= 0) continue;
                Queue<int> q = new Queue<int>();
                q.Enqueue(i); aggr[i] = cur;
                while (q.Count > 0)
                {
                    int u = q.Dequeue();
                    double th = THETA * Math.Sqrt(Math.Abs(d[u]));
                    for (int k = A.RowPtr[u]; k < A.RowPtr[u + 1]; k++)
                    {
                        int v = A.ColIdx[k];
                        if (aggr[v] >= 0) continue;
                        if (Math.Abs(A.Values[k]) >= th)
                        {
                            aggr[v] = cur;
                            q.Enqueue(v);
                        }
                    }
                }
                cur++;
            }
            return aggr;
        }

        private static SparseMatrixCSR BuildPTent(SparseMatrixCSR A, int[] aggr, int nc)
        {
            int n = A.Rows;
            CooBuilder cb = new CooBuilder(n, nc);
            for (int i = 0; i < n; i++) cb.Add(i, aggr[i], 1.0);
            return cb.Build();
        }

        private static SparseMatrixCSR SmoothP(SparseMatrixCSR A, SparseMatrixCSR P, double[] diagInv)
        {
            int n = A.Rows;
            int nc = P.Cols;
            CooBuilder cb = new CooBuilder(n, nc);
            double[] tmp = new double[nc];
            for (int i = 0; i < n; i++)
            {
                DenseVector ei = new DenseVector(n); ei[i] = 1.0;
                DenseVector ap = new DenseVector(n);
                A.MatVec(ei, ap);
                double w = OMEGA * diagInv[i];
                for (int k = P.RowPtr[i]; k < P.RowPtr[i + 1]; k++)
                {
                    int col = P.ColIdx[k];
                    tmp[col] = P.Values[k] - w * ap[i];
                }
                double[] rowsum = new double[nc];
                for (int k = A.RowPtr[i]; k < A.RowPtr[i + 1]; k++)
                {
                    int j = A.ColIdx[k];
                    double aij = A.Values[k];
                    for (int kk = P.RowPtr[j]; kk < P.RowPtr[j + 1]; kk++)
                        rowsum[P.ColIdx[kk]] += w * aij * P.Values[kk];
                }
                for (int c = 0; c < nc; c++)
                {
                    double v = tmp[c] + rowsum[c];
                    if (Math.Abs(v) > 1.0e-20) cb.Add(i, c, v);
                }
                for (int k = P.RowPtr[i]; k < P.RowPtr[i + 1]; k++) tmp[P.ColIdx[k]] = 0.0;
            }
            return cb.Build();
        }

        private static SparseMatrixCSR Transpose(SparseMatrixCSR A)
        {
            CooBuilder cb = new CooBuilder(A.Cols, A.Rows);
            for (int i = 0; i < A.Rows; i++)
                for (int k = A.RowPtr[i]; k < A.RowPtr[i + 1]; k++)
                    cb.Add(A.ColIdx[k], i, A.Values[k]);
            return cb.Build();
        }

        private static SparseMatrixCSR RAP(SparseMatrixCSR R, SparseMatrixCSR A, SparseMatrixCSR P)
        {
            int rows = R.Rows;
            int cols = P.Cols;
            CooBuilder cb = new CooBuilder(rows, cols);
            int cores = Math.Max(1, Environment.ProcessorCount);
            object lk = new object();
            Parallel.For(0, rows, i =>
            {
                Dictionary<int, double> local = new Dictionary<int, double>();
                for (int kr = R.RowPtr[i]; kr < R.RowPtr[i + 1]; kr++)
                {
                    int k = R.ColIdx[kr];
                    double rik = R.Values[kr];
                    for (int ka = A.RowPtr[k]; ka < A.RowPtr[k + 1]; ka++)
                    {
                        int j = A.ColIdx[ka];
                        double aik = rik * A.Values[ka];
                        for (int kp = P.RowPtr[j]; kp < P.RowPtr[j + 1]; kp++)
                        {
                            int c = P.ColIdx[kp];
                            double v = aik * P.Values[kp];
                            double old;
                            if (local.TryGetValue(c, out old)) local[c] = old + v;
                            else local[c] = v;
                        }
                    }
                }
                lock (lk)
                {
                    foreach (KeyValuePair<int, double> kv in local) cb.Add(i, kv.Key, kv.Value);
                }
            });
            return cb.Build();
        }

        private static double[,] GaussJordanInverse(double[,] M, int n)
        {
            double[,] A = (double[,])M.Clone();
            double[,] inv = new double[n, n];
            for (int i = 0; i < n; i++) inv[i, i] = 1.0;
            for (int i = 0; i < n; i++)
            {
                int piv = i;
                double pv = Math.Abs(A[i, i]);
                for (int r = i + 1; r < n; r++)
                    if (Math.Abs(A[r, i]) > pv) { pv = Math.Abs(A[r, i]); piv = r; }
                if (piv != i)
                {
                    for (int c = 0; c < n; c++)
                    {
                        double t = A[i, c]; A[i, c] = A[piv, c]; A[piv, c] = t;
                        t = inv[i, c]; inv[i, c] = inv[piv, c]; inv[piv, c] = t;
                    }
                }
                double d = A[i, i];
                if (Math.Abs(d) < 1.0e-30) d = 1.0e-30;
                for (int c = 0; c < n; c++) { A[i, c] /= d; inv[i, c] /= d; }
                for (int r = 0; r < n; r++)
                {
                    if (r == i) continue;
                    double f = A[r, i];
                    if (Math.Abs(f) < 1.0e-30) continue;
                    for (int c = 0; c < n; c++)
                    {
                        A[r, c] -= f * A[i, c];
                        inv[r, c] -= f * inv[i, c];
                    }
                }
            }
            return inv;
        }

        private void Smooth(Level L, DenseVector x, DenseVector b, int steps)
        {
            int n = L.N;
            for (int s = 0; s < steps; s++)
            {
                for (int i = 0; i < n; i++)
                {
                    double sum = b[i];
                    for (int k = L.A.RowPtr[i]; k < L.A.RowPtr[i + 1]; k++)
                    {
                        int j = L.A.ColIdx[k];
                        if (j != i) sum -= L.A.Values[k] * x[j];
                    }
                    x[i] += OMEGA * L.DiagInv[i] * (sum - L.Diag[i] * x[i]);
                }
                for (int i = n - 1; i >= 0; i--)
                {
                    double sum = b[i];
                    for (int k = L.A.RowPtr[i]; k < L.A.RowPtr[i + 1]; k++)
                    {
                        int j = L.A.ColIdx[k];
                        if (j != i) sum -= L.A.Values[k] * x[j];
                    }
                    x[i] += OMEGA * L.DiagInv[i] * (sum - L.Diag[i] * x[i]);
                }
            }
        }

        private void Vcycle(int l, DenseVector b, DenseVector x)
        {
            Level L = this._levels[l];
            if (L.Ainv != null)
            {
                int n = L.N;
                for (int i = 0; i < n; i++)
                {
                    double s = 0.0;
                    for (int j = 0; j < n; j++) s += L.Ainv[i, j] * b[j];
                    x[i] = s;
                }
                return;
            }
            Smooth(L, x, b, NU1);
            DenseVector r = new DenseVector(L.N);
            L.A.MatVec(x, r);
            for (int i = 0; i < L.N; i++) r[i] = b[i] - r[i];
            DenseVector rc = new DenseVector(L.R.Rows);
            L.R.MatVec(r, rc);
            DenseVector ec = new DenseVector(L.R.Rows);
            Vcycle(l + 1, rc, ec);
            DenseVector e = new DenseVector(L.N);
            L.P.MatVec(ec, e);
            x.Axpy(1.0, e);
            Smooth(L, x, b, NU2);
        }

        public void Apply(DenseVector r, DenseVector z)
        {
            z.SetZero();
            Vcycle(0, r, z);
        }
    }
}
