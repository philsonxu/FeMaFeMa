using System;
using System.Collections.Generic;
using System.Numerics;
using System.Threading;

namespace MultiPhysicsFEM2D.Core
{
    public class SparseMatrixCSR
    {
        public int Nrows;
        public int Ncols;
        public int Nnz;
        public int[] RowPtr;
        public int[] ColIdx;
        public double[] Values;
        public bool IsSquare { get { return Nrows == Ncols; } }

        public SparseMatrixCSR() { }
        public SparseMatrixCSR(int nr, int nc, int[] rp, int[] ci, double[] v)
        {
            Nrows = nr; Ncols = nc; RowPtr = rp; ColIdx = ci; Values = v; Nnz = v.Length;
        }

        public void MultiplyAdd(double alpha, DenseVector x, double beta, DenseVector y)
        {
            int n = Nrows;
            int vlen = Vector<double>.Count;
            if (n < 2000)
            {
                for (int i = 0; i < n; i++)
                {
                    int r0 = RowPtr[i], r1 = RowPtr[i + 1];
                    double s = 0.0;
                    int k;
                    for (k = r0; k + vlen <= r1; k += vlen)
                    {
                        Span<double> xseg = stackalloc double[vlen];
                        for (int j = 0; j < vlen; j++) xseg[j] = x.Values[ColIdx[k + j]];
                        Vector<double> xv = new Vector<double>(xseg);
                        Vector<double> vv = new Vector<double>(Values, k);
                        s += Vector.Dot(vv, xv);
                    }
                    for (; k < r1; k++) s += Values[k] * x.Values[ColIdx[k]];
                    y.Values[i] = beta * y.Values[i] + alpha * s;
                }
                return;
            }
            int threads = Environment.ProcessorCount;
            WaitHandle[] whs = new WaitHandle[threads];
            int block = (n + threads - 1) / threads;
            for (int t = 0; t < threads; t++)
            {
                int start = t * block, end = Math.Min(start + block, n);
                ManualResetEvent mre = new ManualResetEvent(false); whs[t] = mre;
                ThreadPool.QueueUserWorkItem(delegate
                {
                    for (int i = start; i < end; i++)
                    {
                        int r0 = RowPtr[i], r1 = RowPtr[i + 1];
                        double s = 0.0;
                        int k;
                        for (k = r0; k + vlen <= r1; k += vlen)
                        {
                            Span<double> xseg = stackalloc double[vlen];
                            for (int j = 0; j < vlen; j++) xseg[j] = x.Values[ColIdx[k + j]];
                            Vector<double> xv = new Vector<double>(xseg);
                            Vector<double> vv = new Vector<double>(Values, k);
                            s += Vector.Dot(vv, xv);
                        }
                        for (; k < r1; k++) s += Values[k] * x.Values[ColIdx[k]];
                        y.Values[i] = beta * y.Values[i] + alpha * s;
                    }
                    mre.Set();
                });
            }
            WaitHandle.WaitAll(whs);
        }

        public void Multiply(DenseVector x, DenseVector y)
        {
            y.SetZero();
            MultiplyAdd(1.0, x, 0.0, y);
        }

        public DenseVector Multiply(DenseVector x)
        {
            DenseVector y = new DenseVector(Nrows);
            Multiply(x, y);
            return y;
        }

        public DenseVector ExtractDiagonal()
        {
            DenseVector d = new DenseVector(Nrows);
            for (int i = 0; i < Nrows; i++)
                for (int k = RowPtr[i]; k < RowPtr[i + 1]; k++)
                    if (ColIdx[k] == i) { d.Values[i] = Values[k]; break; }
            return d;
        }

        public void ApplyDirichlet(int dof, double penalty, double value, DenseVector b)
        {
            int r0 = RowPtr[dof], r1 = RowPtr[dof + 1];
            for (int k = r0; k < r1; k++)
            {
                if (ColIdx[k] == dof) Values[k] = penalty;
                else Values[k] = 0.0;
            }
            b.Values[dof] = penalty * value;
        }
    }

    public class CooAssembler
    {
        private int _nr, _nc;
        private List<Tuple<int, int, double>> _entries = new List<Tuple<int, int, double>>();
        private object _lock = new object();

        public CooAssembler(int nr, int nc) { _nr = nr; _nc = nc; }

        public void AddEntry(int r, int c, double v)
        {
            lock (_lock) _entries.Add(Tuple.Create(r, c, v));
        }

        public void AddLocal(int[] dofMap, double[,] localM)
        {
            int ldof = dofMap.Length;
            lock (_lock)
            {
                for (int i = 0; i < ldof; i++)
                {
                    int gi = dofMap[i]; if (gi < 0 || gi >= _nr) continue;
                    for (int j = 0; j < ldof; j++)
                    {
                        int gj = dofMap[j]; if (gj < 0 || gj >= _nc) continue;
                        _entries.Add(Tuple.Create(gi, gj, localM[i, j]));
                    }
                }
            }
        }

        public void AddLocalBulk(List<Tuple<int[], double[,]>> batch)
        {
            lock (_lock)
            {
                foreach (Tuple<int[], double[,]> it in batch)
                {
                    int[] map = it.Item1; double[,] lm = it.Item2;
                    int ldof = map.Length;
                    for (int i = 0; i < ldof; i++)
                    {
                        int gi = map[i]; if (gi < 0 || gi >= _nr) continue;
                        for (int j = 0; j < ldof; j++)
                        {
                            int gj = map[j]; if (gj < 0 || gj >= _nc) continue;
                            _entries.Add(Tuple.Create(gi, gj, lm[i, j]));
                        }
                    }
                }
            }
        }

        public SparseMatrixCSR BuildCSR()
        {
            Dictionary<long, double> dict = new Dictionary<long, double>();
            foreach (Tuple<int, int, double> e in _entries)
            {
                long key = ((long)e.Item1 << 32) | (uint)e.Item2;
                double cur;
                if (dict.TryGetValue(key, out cur)) dict[key] = cur + e.Item3;
                else dict[key] = e.Item3;
            }
            List<KeyValuePair<long, double>> lst = new List<KeyValuePair<long, double>>(dict);
            lst.Sort(delegate (KeyValuePair<long, double> a, KeyValuePair<long, double> b)
            {
                int ra = (int)(a.Key >> 32), rb = (int)(b.Key >> 32);
                if (ra != rb) return ra.CompareTo(rb);
                int ca = (int)(a.Key & 0xFFFFFFFF), cb = (int)(b.Key & 0xFFFFFFFF);
                return ca.CompareTo(cb);
            });
            int nnz = lst.Count;
            int[] rp = new int[_nr + 1];
            int[] ci = new int[nnz];
            double[] vv = new double[nnz];
            int curRow = -1, p = 0;
            foreach (KeyValuePair<long, double> e in lst)
            {
                int r = (int)(e.Key >> 32), c = (int)(e.Key & 0xFFFFFFFF);
                while (curRow < r) { curRow++; rp[curRow] = p; }
                ci[p] = c; vv[p] = e.Value; p++;
            }
            while (curRow < _nr) { curRow++; rp[curRow] = nnz; }
            return new SparseMatrixCSR(_nr, _nc, rp, ci, vv);
        }
    }
}
