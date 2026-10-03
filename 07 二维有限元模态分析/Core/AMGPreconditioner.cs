using System;
using System.Collections.Generic;
using System.Threading;

namespace ModalFEM2D.Core
{

    /// <summary>
    /// Algebraic multigrid preconditioner (Smoothed Aggregation variant).
    /// Coarsens via strong-connection aggregation; uses weighted-Jacobi smoothing on each level;
    /// exact Gauss-Jordan solve on the coarsest level.
    /// </summary>
    public sealed class AMGPreconditioner : IPreconditioner
    {
        private readonly List<AMGLevel> _levels = new List<AMGLevel>();
        private const int CoarseMax = 80;
        private const int SmoothSweeps = 2;
        private const double Omega = 0.6;
        private const double Theta = 0.08;

        private sealed class AMGLevel
        {
            public SparseMatrixCSR A;
            public SparseMatrixCSR P;       // prolongation  (fine -> coarse)
            public SparseMatrixCSR R;       // restriction   (coarse -> fine = P^T)
            public double[,] CoarseInv;     // dense inverse (only set for the coarsest)
        }

        public AMGPreconditioner(SparseMatrixCSR A)
        {
            SparseMatrixCSR cur = A;
            while (cur.NumRows > CoarseMax)
            {
                Tuple<SparseMatrixCSR, SparseMatrixCSR> pr = BuildProlongation(cur);
                SparseMatrixCSR P = pr.Item1;
                SparseMatrixCSR R = pr.Item2;
                SparseMatrixCSR Ac = GalerkinRAP(cur, P);
                if (Ac.NumRows >= cur.NumRows) break; // no further coarsening
                _levels.Add(new AMGLevel { A = cur, P = P, R = R });
                cur = Ac;
            }
            _levels.Add(new AMGLevel { A = cur, CoarseInv = DenseInverse(cur) });
        }

        private static Tuple<SparseMatrixCSR, SparseMatrixCSR> BuildProlongation(SparseMatrixCSR A)
        {
            int n = A.NumRows;
            double[] aii = new double[n];
            for (int i = 0; i < n; i++)
                for (int k = A.RowPtr[i]; k < A.RowPtr[i + 1]; k++)
                    if (A.ColIdx[k] == i) { aii[i] = Math.Abs(A.Values[k]); break; }

            int[] agg = new int[n];
            for (int i = 0; i < n; i++) agg[i] = -1;
            int nAgg = 0;
            int[] order = new int[n];
            for (int i = 0; i < n; i++) order[i] = i;
            Array.Sort(aii, order); Array.Reverse(order);
            foreach (int seed in order)
            {
                if (agg[seed] >= 0) continue;
                Queue<int> q = new Queue<int>();
                q.Enqueue(seed); agg[seed] = nAgg;
                while (q.Count > 0)
                {
                    int i = q.Dequeue();
                    for (int k = A.RowPtr[i]; k < A.RowPtr[i + 1]; k++)
                    {
                        int j = A.ColIdx[k]; if (j == i) continue;
                        double aij = Math.Abs(A.Values[k]);
                        double thr = Theta * Math.Sqrt(aii[i] * aii[j] + 1e-30);
                        if (aij >= thr && agg[j] < 0) { agg[j] = nAgg; q.Enqueue(j); }
                    }
                }
                nAgg++;
            }
            for (int i = 0; i < n; i++) if (agg[i] < 0) agg[i] = nAgg++;

            // Tentative P0 (piecewise constant) + one damped-Jacobi smoothing pass.
            int[] pRP = new int[n + 1];
            List<int> pCI = new List<int>();
            List<double> pVA = new List<double>();
            for (int i = 0; i < n; i++)
            {
                pRP[i] = pVA.Count;
                Dictionary<int, double> contrib = new Dictionary<int, double>();
                contrib[agg[i]] = 1.0;
                double di = Math.Abs(aii[i]) > 1e-30 ? aii[i] : 1.0;
                double omegaDi = Omega / di;
                // P_smooth = (I - omega D^{-1}A) P0
                for (int k = A.RowPtr[i]; k < A.RowPtr[i + 1]; k++)
                {
                    int j = A.ColIdx[k];
                    AddTo(contrib, agg[j], -omegaDi * A.Values[k]);
                    if (j == i) AddTo(contrib, agg[j], omegaDi * A.Values[k]); // compensate diagonal part of I
                }
                foreach (KeyValuePair<int, double> kv in contrib)
                {
                    if (Math.Abs(kv.Value) < 1e-14) continue;
                    pCI.Add(kv.Key); pVA.Add(kv.Value);
                }
            }
            pRP[n] = pVA.Count;
            SparseMatrixCSR P = new SparseMatrixCSR(n, nAgg, pRP, pCI.ToArray(), pVA.ToArray());
            SparseMatrixCSR R = Transpose(P);
            return Tuple.Create(P, R);
        }

        private static void AddTo(Dictionary<int, double> d, int k, double v)
        {
            double cur; d[k] = d.TryGetValue(k, out cur) ? cur + v : v;
        }

        private static SparseMatrixCSR Transpose(SparseMatrixCSR A)
        {
            int m = A.NumRows, n = A.NumCols, nnz = A.Values.Length;
            int[] rp = new int[n + 1];
            int[] cnt = new int[n];
            for (int k = 0; k < nnz; k++) cnt[A.ColIdx[k]]++;
            int p = 0; for (int j = 0; j < n; j++) { rp[j] = p; p += cnt[j]; }
            rp[n] = p;
            Array.Copy(rp, cnt, n);
            int[] ci = new int[nnz]; double[] va = new double[nnz];
            for (int i = 0; i < m; i++)
                for (int k = A.RowPtr[i]; k < A.RowPtr[i + 1]; k++)
                {
                    int j = A.ColIdx[k]; int d = cnt[j]++; ci[d] = i; va[d] = A.Values[k];
                }
            return new SparseMatrixCSR(n, m, rp, ci, va);
        }

        /// <summary>Compute Ac = R * A * P = P^T * A * P (symmetric Galerkin).</summary>
        private static SparseMatrixCSR GalerkinRAP(SparseMatrixCSR A, SparseMatrixCSR P)
        {
            SparseMatrixCSR R = Transpose(P);
            int m = R.NumRows; int n = P.NumCols;
            Dictionary<int, double>[] rows = new Dictionary<int, double>[m];
            for (int i = 0; i < m; i++) rows[i] = new Dictionary<int, double>();
            for (int i = 0; i < m; i++)
            {
                Dictionary<int, double> tmp = new Dictionary<int, double>();
                for (int k = R.RowPtr[i]; k < R.RowPtr[i + 1]; k++)
                {
                    int row = R.ColIdx[k]; double rv = R.Values[k];
                    for (int a = A.RowPtr[row]; a < A.RowPtr[row + 1]; a++)
                    {
                        int col = A.ColIdx[a]; AddTo(tmp, col, rv * A.Values[a]);
                    }
                }
                foreach (KeyValuePair<int, double> kv in tmp)
                {
                    int k = kv.Key; double ak = kv.Value;
                    for (int a = P.RowPtr[k]; a < P.RowPtr[k + 1]; a++)
                    {
                        int c = P.ColIdx[a]; AddTo(rows[i], c, ak * P.Values[a]);
                    }
                }
            }
            List<Tuple<int, int, double>> triples = new List<Tuple<int, int, double>>();
            for (int i = 0; i < m; i++)
                foreach (KeyValuePair<int, double> kv in rows[i])
                    triples.Add(Tuple.Create(i, kv.Key, kv.Value));
            return SparseMatrixCSR.FromCoo(m, n, triples);
        }

        private static double[,] DenseInverse(SparseMatrixCSR A)
        {
            int n = A.NumRows;
            double[,] M = new double[n, n];
            for (int i = 0; i < n; i++)
                for (int k = A.RowPtr[i]; k < A.RowPtr[i + 1]; k++)
                    M[i, A.ColIdx[k]] = A.Values[k];
            double[,] I = new double[n, n];
            for (int i = 0; i < n; i++) I[i, i] = 1.0;
            for (int i = 0; i < n; i++)
            {
                int pv = i; double mval = Math.Abs(M[i, i]);
                for (int r = i + 1; r < n; r++) if (Math.Abs(M[r, i]) > mval) { pv = r; mval = Math.Abs(M[r, i]); }
                if (pv != i)
                {
                    for (int c = 0; c < n; c++)
                    {
                        double t = M[i, c]; M[i, c] = M[pv, c]; M[pv, c] = t;
                        t = I[i, c]; I[i, c] = I[pv, c]; I[pv, c] = t;
                    }
                }
                double d = M[i, i]; if (Math.Abs(d) < 1e-14) d = 1e-14;
                for (int c = 0; c < n; c++) { M[i, c] /= d; I[i, c] /= d; }
                for (int r = 0; r < n; r++)
                {
                    if (r == i) continue;
                    double f = M[r, i]; if (Math.Abs(f) < 1e-30) continue;
                    for (int c = 0; c < n; c++) { M[r, c] -= f * M[i, c]; I[r, c] -= f * I[i, c]; }
                }
            }
            return I;
        }

        public void Apply(DenseVector r, DenseVector z)
        {
            z.Clear();
            VCycle(0, r, z);
        }

        private void VCycle(int lvl, DenseVector rb, DenseVector xb)
        {
            AMGLevel L = _levels[lvl];
            if (L.CoarseInv != null)
            {
                int n = L.CoarseInv.GetLength(0);
                for (int i = 0; i < n; i++)
                {
                    double s = 0.0;
                    for (int j = 0; j < n; j++) s += L.CoarseInv[i, j] * rb.Values[j];
                    xb.Values[i] = s;
                }
                return;
            }
            for (int s = 0; s < SmoothSweeps; s++) JacobiSmooth(L.A, xb, rb);
            DenseVector Ax = new DenseVector(L.A.NumRows);
            L.A.MatVec(xb, Ax);
            DenseVector res = new DenseVector(L.A.NumRows);
            for (int i = 0; i < L.A.NumRows; i++) res.Values[i] = rb.Values[i] - Ax.Values[i];
            int nc = L.P.NumCols;
            DenseVector rc = new DenseVector(nc);
            L.R.MatVec(res, rc);
            DenseVector ec = new DenseVector(nc);
            VCycle(lvl + 1, rc, ec);
            DenseVector ef = new DenseVector(L.A.NumRows);
            L.P.MatVec(ec, ef);
            xb.Axpy(1.0, ef);
            for (int s = 0; s < SmoothSweeps; s++) JacobiSmooth(L.A, xb, rb);
        }

        private static void JacobiSmooth(SparseMatrixCSR A, DenseVector x, DenseVector b)
        {
            int n = A.NumRows;
            DenseVector Ax = new DenseVector(n);
            A.MatVec(x, Ax);
            for (int i = 0; i < n; i++)
            {
                double di = 0.0;
                for (int k = A.RowPtr[i]; k < A.RowPtr[i + 1]; k++)
                    if (A.ColIdx[k] == i) { di = A.Values[k]; break; }
                if (Math.Abs(di) < 1e-30) di = 1.0;
                x.Values[i] += Omega * (b.Values[i] - Ax.Values[i]) / di;
            }
        }
    }
}
