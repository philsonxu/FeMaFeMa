using System;
using System.Collections.Generic;
using System.Threading;

namespace MultiPhysicsFEM2D.Core
{
    public class AMGPreconditioner : IPreconditioner
    {
        private class Level
        {
            public SparseMatrixCSR A;
            public SparseMatrixCSR P;
            public SparseMatrixCSR R;
            public double[] DiagInv;
            public double[,] AcDenseInv;
            public int CoarseN;
        }

        private List<Level> _levels;
        private const int MAX_LEVELS = 8;
        private const int COARSE_DIRECT = 64;
        private const double THETA = 0.08;
        private const double OMEGA_JAC = 0.7;
        private const int NU1 = 2, NU2 = 2;

        public void Build(SparseMatrixCSR A)
        {
            _levels = new List<Level>();
            SparseMatrixCSR cur = A;
            while (true)
            {
                Level lv = new Level();
                lv.A = cur;
                int n = cur.Nrows;
                lv.DiagInv = new double[n];
                for (int i = 0; i < n; i++)
                {
                    double d = 1.0;
                    for (int k = cur.RowPtr[i]; k < cur.RowPtr[i + 1]; k++)
                        if (cur.ColIdx[k] == i) { d = cur.Values[k]; break; }
                    lv.DiagInv[i] = Math.Abs(d) < 1e-30 ? 1.0 : 1.0 / d;
                }
                if (n <= COARSE_DIRECT || _levels.Count >= MAX_LEVELS - 1)
                {
                    lv.AcDenseInv = DenseInverse(cur);
                    lv.CoarseN = n;
                    _levels.Add(lv);
                    break;
                }
                int[] agg = Aggregate(cur, THETA);
                int nc = 0;
                for (int i = 0; i < n; i++) if (agg[i] >= nc) nc = agg[i] + 1;
                lv.CoarseN = nc;
                SparseMatrixCSR P0 = BuildTentativeP(n, nc, agg);
                SparseMatrixCSR P = SmoothProlongation(cur, P0, lv.DiagInv);
                SparseMatrixCSR R = Transpose(P);
                SparseMatrixCSR Ac = RAP(R, cur, P);
                lv.P = P; lv.R = R;
                _levels.Add(lv);
                cur = Ac;
            }
        }

        public void Apply(DenseVector r, DenseVector z)
        {
            z.SetZero();
            Vcycle(0, r, z);
        }

        private void Vcycle(int lv, DenseVector b, DenseVector x)
        {
            Level L = _levels[lv];
            int n = L.A.Nrows;
            for (int s = 0; s < NU1; s++) SmoothJ(L.A, L.DiagInv, b, x);
            if (L.AcDenseInv != null)
            {
                DenseVector res = new DenseVector(n);
                L.A.Multiply(x, res);
                for (int i = 0; i < n; i++) res.Values[i] = b.Values[i] - res.Values[i];
                double[] corr = new double[n];
                for (int i = 0; i < n; i++)
                    for (int j = 0; j < n; j++) corr[i] += L.AcDenseInv[i, j] * res.Values[j];
                for (int i = 0; i < n; i++) x.Values[i] += corr[i];
            }
            else
            {
                DenseVector res = new DenseVector(n);
                L.A.Multiply(x, res);
                for (int i = 0; i < n; i++) res.Values[i] = b.Values[i] - res.Values[i];
                int nc = L.CoarseN;
                DenseVector rc = new DenseVector(nc);
                L.R.MultiplyAdd(1.0, res, 0.0, rc);
                DenseVector xc = new DenseVector(nc);
                Vcycle(lv + 1, rc, xc);
                L.P.MultiplyAdd(1.0, xc, 1.0, x);
            }
            for (int s = 0; s < NU2; s++) SmoothJ(L.A, L.DiagInv, b, x);
        }

        private static void SmoothJ(SparseMatrixCSR A, double[] di, DenseVector b, DenseVector x)
        {
            int n = A.Nrows;
            DenseVector Ax = new DenseVector(n);
            A.Multiply(x, Ax);
            for (int i = 0; i < n; i++)
            {
                double r = b.Values[i] - Ax.Values[i];
                x.Values[i] += OMEGA_JAC * di[i] * r;
            }
        }

        private static double[,] DenseInverse(SparseMatrixCSR A)
        {
            int n = A.Nrows;
            double[,] M = new double[n, n];
            for (int i = 0; i < n; i++)
                for (int k = A.RowPtr[i]; k < A.RowPtr[i + 1]; k++)
                    M[i, A.ColIdx[k]] = A.Values[k];
            double[,] aug = new double[n, 2 * n];
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++) aug[i, j] = M[i, j];
                aug[i, n + i] = 1.0;
            }
            for (int col = 0; col < n; col++)
            {
                int piv = col;
                double mx = Math.Abs(aug[col, col]);
                for (int r = col + 1; r < n; r++)
                    if (Math.Abs(aug[r, col]) > mx) { mx = Math.Abs(aug[r, col]); piv = r; }
                if (piv != col)
                    for (int c = 0; c < 2 * n; c++) { double t = aug[col, c]; aug[col, c] = aug[piv, c]; aug[piv, c] = t; }
                double d = aug[col, col];
                if (Math.Abs(d) < 1e-30) d = 1e-30;
                for (int c = 0; c < 2 * n; c++) aug[col, c] /= d;
                for (int r = 0; r < n; r++)
                {
                    if (r == col) continue;
                    double f = aug[r, col];
                    if (Math.Abs(f) < 1e-30) continue;
                    for (int c = 0; c < 2 * n; c++) aug[r, c] -= f * aug[col, c];
                }
            }
            double[,] inv = new double[n, n];
            for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) inv[i, j] = aug[i, n + j];
            return inv;
        }

        private static int[] Aggregate(SparseMatrixCSR A, double theta)
        {
            int n = A.Nrows;
            double[] diag = new double[n];
            for (int i = 0; i < n; i++)
                for (int k = A.RowPtr[i]; k < A.RowPtr[i + 1]; k++)
                    if (A.ColIdx[k] == i) { diag[i] = Math.Abs(A.Values[k]); break; }
            int[] order = new int[n];
            for (int i = 0; i < n; i++) order[i] = i;
            Array.Sort<int>(order, delegate (int a, int b) { return diag[b].CompareTo(diag[a]); });
            int[] agg = new int[n];
            for (int i = 0; i < n; i++) agg[i] = -1;
            int cur = 0;
            Queue<int> q = new Queue<int>();
            bool[] inQ = new bool[n];
            for (int ii = 0; ii < n; ii++)
            {
                int seed = order[ii];
                if (agg[seed] != -1) continue;
                q.Clear(); Array.Clear(inQ, 0, n);
                q.Enqueue(seed); agg[seed] = cur; inQ[seed] = true;
                while (q.Count > 0)
                {
                    int i = q.Dequeue();
                    for (int k = A.RowPtr[i]; k < A.RowPtr[i + 1]; k++)
                    {
                        int j = A.ColIdx[k];
                        if (j == i || agg[j] != -1) continue;
                        double th = theta * Math.Sqrt(Math.Abs(diag[i] * diag[j]));
                        if (Math.Abs(A.Values[k]) > th && !inQ[j])
                        {
                            agg[j] = cur; inQ[j] = true; q.Enqueue(j);
                        }
                    }
                }
                cur++;
            }
            return agg;
        }

        private static SparseMatrixCSR BuildTentativeP(int n, int nc, int[] agg)
        {
            int[] rp = new int[n + 1];
            List<int> cl = new List<int>();
            List<double> vl = new List<double>();
            int p = 0;
            for (int i = 0; i < n; i++)
            {
                rp[i] = p;
                cl.Add(agg[i]); vl.Add(1.0); p++;
            }
            rp[n] = p;
            return new SparseMatrixCSR(n, nc, rp, cl.ToArray(), vl.ToArray());
        }

        private static SparseMatrixCSR SmoothProlongation(SparseMatrixCSR A, SparseMatrixCSR P0, double[] diagInv)
        {
            int n = A.Nrows, nc = P0.Ncols;
            int threads = Environment.ProcessorCount;
            WaitHandle[] whs = new WaitHandle[threads];
            int block = (n + threads - 1) / threads;
            List<Tuple<int, int, double>>[] per = new List<Tuple<int, int, double>>[threads];
            for (int t = 0; t < threads; t++) per[t] = new List<Tuple<int, int, double>>();
            for (int t = 0; t < threads; t++)
            {
                int start = t * block, end = Math.Min(start + block, n);
                int tt = t;
                ManualResetEvent mre = new ManualResetEvent(false); whs[t] = mre;
                ThreadPool.QueueUserWorkItem(delegate
                {
                    for (int i = start; i < end; i++)
                    {
                        Dictionary<int, double> row = new Dictionary<int, double>();
                        for (int k = P0.RowPtr[i]; k < P0.RowPtr[i + 1]; k++) row[P0.ColIdx[k]] = P0.Values[k];
                        double coef = -OMEGA_JAC * diagInv[i];
                        for (int a = A.RowPtr[i]; a < A.RowPtr[i + 1]; a++)
                        {
                            int ja = A.ColIdx[a]; double av = A.Values[a];
                            for (int pk = P0.RowPtr[ja]; pk < P0.RowPtr[ja + 1]; pk++)
                            {
                                int c = P0.ColIdx[pk];
                                double v = coef * av * P0.Values[pk];
                                double cur;
                                if (row.TryGetValue(c, out cur)) row[c] = cur + v;
                                else row[c] = v;
                            }
                        }
                        foreach (KeyValuePair<int, double> kv in row)
                            if (Math.Abs(kv.Value) > 1e-15) per[tt].Add(Tuple.Create(i, kv.Key, kv.Value));
                    }
                    mre.Set();
                });
            }
            WaitHandle.WaitAll(whs);
            List<Tuple<int, int, double>> all = new List<Tuple<int, int, double>>();
            for (int t = 0; t < threads; t++) all.AddRange(per[t]);
            all.Sort(delegate (Tuple<int, int, double> a, Tuple<int, int, double> b)
            {
                if (a.Item1 != b.Item1) return a.Item1.CompareTo(b.Item1);
                return a.Item2.CompareTo(b.Item2);
            });
            int nnz = all.Count;
            int[] rp = new int[n + 1], ci = new int[nnz];
            double[] vv = new double[nnz];
            int curR = -1, p = 0;
            foreach (Tuple<int, int, double> e in all)
            {
                while (curR < e.Item1) { curR++; rp[curR] = p; }
                ci[p] = e.Item2; vv[p] = e.Item3; p++;
            }
            while (curR < n) { curR++; rp[curR] = nnz; }
            return new SparseMatrixCSR(n, nc, rp, ci, vv);
        }

        private static SparseMatrixCSR Transpose(SparseMatrixCSR A)
        {
            int m = A.Nrows, n = A.Ncols, nnz = A.Nnz;
            int[] cnt = new int[n + 1];
            for (int k = 0; k < nnz; k++) cnt[A.ColIdx[k] + 1]++;
            for (int i = 1; i <= n; i++) cnt[i] += cnt[i - 1];
            int[] ptr = (int[])cnt.Clone();
            int[] rp = new int[n + 1];
            Array.Copy(cnt, rp, n + 1);
            int[] ci = new int[nnz]; double[] vv = new double[nnz];
            for (int i = 0; i < m; i++)
                for (int k = A.RowPtr[i]; k < A.RowPtr[i + 1]; k++)
                {
                    int col = A.ColIdx[k];
                    int pos = ptr[col]++;
                    ci[pos] = i; vv[pos] = A.Values[k];
                }
            return new SparseMatrixCSR(n, m, rp, ci, vv);
        }

        private static SparseMatrixCSR RAP(SparseMatrixCSR R, SparseMatrixCSR A, SparseMatrixCSR P)
        {
            int nF = A.Nrows, nC = P.Ncols;
            CooAssembler asm = new CooAssembler(nC, nC);
            int threads = Environment.ProcessorCount;
            WaitHandle[] whs = new WaitHandle[threads];
            int block = (nC + threads - 1) / threads;
            for (int t = 0; t < threads; t++)
            {
                int start = t * block, end = Math.Min(start + block, nC);
                ManualResetEvent mre = new ManualResetEvent(false); whs[t] = mre;
                ThreadPool.QueueUserWorkItem(delegate
                {
                    // 源代码错误
                    //List<Tuple<int, int, double>> local = new List<Tuple<int, int, double>>();
                    List<Tuple<int[], double[,]>> local2 = new List<Tuple<int[], double[,]>>();
                    double[] Pcol = new double[nF];
                    double[] Ap = new double[nF];
                    for (int j = start; j < end; j++)
                    {
                        Array.Clear(Pcol, 0, nF);
                        for (int ip = 0; ip < nF; ip++)
                            for (int kp = P.RowPtr[ip]; kp < P.RowPtr[ip + 1]; kp++)
                                if (P.ColIdx[kp] == j) Pcol[ip] = P.Values[kp];
                        Array.Clear(Ap, 0, nF);
                        for (int i = 0; i < nF; i++)
                        {
                            double s = 0;
                            for (int k = A.RowPtr[i]; k < A.RowPtr[i + 1]; k++)
                                s += A.Values[k] * Pcol[A.ColIdx[k]];
                            Ap[i] = s;
                        }
                        for (int i = 0; i < nC; i++)
                        {
                            double s = 0;
                            for (int rk = R.RowPtr[i]; rk < R.RowPtr[i + 1]; rk++)
                                s += R.Values[rk] * Ap[R.ColIdx[rk]];
                            if (Math.Abs(s) > 1e-15)
                            {
                                //local.Add(Tuple.Create(i, j, s));
                                local2.Add(Tuple.Create(new int[] { i, j }, new double[1, 1] { { s } }));
                            }
                        }
                    }
                    //asm.AddLocalBulk(local);
                    asm.AddLocalBulk(local2);
                    mre.Set();
                });
            }
            WaitHandle.WaitAll(whs);
            return asm.BuildCSR();
        }
    }
}
