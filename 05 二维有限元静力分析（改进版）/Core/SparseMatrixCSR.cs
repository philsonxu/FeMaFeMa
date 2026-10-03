using System;
using System.Collections.Generic;
using System.Numerics;
using System.Threading;

namespace FEM2D.Core
{
    public sealed class SparseMatrixCSR
    {
        public int N;
        public int[] RowPtr;
        public int[] ColIdx;
        public double[] Values;
        public double[] DiagInv;

        public SparseMatrixCSR(int n, int[] rowptr, int[] colidx, double[] vals)
        {
            N = n; RowPtr = rowptr; ColIdx = colidx; Values = vals;
            DiagInv = new double[n];
            for (int i = 0; i < n; i++)
            {
                double d = 0.0;
                for (int p = RowPtr[i]; p < RowPtr[i + 1]; p++)
                    if (ColIdx[p] == i) { d = Values[p]; break; }
                DiagInv[i] = (Math.Abs(d) > 1e-30) ? 1.0 / d : 1.0;
            }
        }

        public int Rows { get { return N; } }

        /// <summary>
        /// 多核 + SIMD 向量化 SpMV：y = A * x。
        /// 策略：按行分块，每块交给线程池；每行使用 Vector&lt;double&gt; 4/8 路同时累加非零项，
        /// 尾部不足向量宽度的用标量收尾。
        /// </summary>
        public void Multiply(DenseVector x, DenseVector y)
        {
            int n = N;
            Array.Clear(y.Values, 0, n);
            int vecLen = Vector<double>.Count;
            int procCount = Environment.ProcessorCount;
            int threads = Math.Min(procCount <= 0 ? 4 : procCount, 16);
            if (threads < 1) threads = 1;
            if (n < 256) threads = 1; // 小规模串行避免调度开销
            int chunk = (n + threads - 1) / threads;
            int[] rp = RowPtr, ci = ColIdx;
            double[] va = Values, xv = x.Values, yv = y.Values;

            if (threads == 1)
            {
                SpmvRange(0, n, rp, ci, va, xv, yv, vecLen);
                return;
            }

            ManualResetEvent[] done = new ManualResetEvent[threads];
            int launched = 0;
            for (int t = 0; t < threads; t++)
            {
                int from = t * chunk;
                int to = Math.Min(from + chunk, n);
                if (from >= to) break;
                done[t] = new ManualResetEvent(false);
                launched++;
                int tid = t;
                ThreadPool.QueueUserWorkItem(delegate(object s)
                {
                    int idx = (int)s;
                    int a2 = idx * chunk, b2 = Math.Min(a2 + chunk, n);
                    SpmvRange(a2, b2, rp, ci, va, xv, yv, vecLen);
                    done[idx].Set();
                }, tid);
            }
            if (launched > 0)
            {
                WaitHandle[] wait = new WaitHandle[launched];
                for (int t = 0; t < launched; t++) wait[t] = done[t];
                WaitHandle.WaitAll(wait);
                for (int t = 0; t < launched; t++) done[t].Close();
            }
        }

        private static void SpmvRange(int a, int b, int[] rp, int[] ci, double[] va, double[] xv, double[] yv, int vecLen)
        {
            // 使用栈上缓冲避免每行分配；vecLen 最大 8（AVX512 下为8，SSE为2，AVX2为4）
            Span<double> tmp = stackalloc double[8];
            for (int i = a; i < b; i++)
            {
                int pStart = rp[i];
                int pEnd = rp[i + 1];
                int p = pStart;
                int pSimdEnd = pStart + ((pEnd - pStart) / vecLen) * vecLen;
                Vector<double> acc = Vector<double>.Zero;
                for (; p < pSimdEnd; p += vecLen)
                {
                    Vector<double> valsVec = new Vector<double>(va, p);
                    Span<double> slice = tmp.Slice(0, vecLen);
                    for (int k = 0; k < vecLen; k++) slice[k] = xv[ci[p + k]];
                    Vector<double> xvec = new Vector<double>(slice);
                    acc += valsVec * xvec;
                }
                double sum = 0.0;
                for (int k = 0; k < vecLen; k++) sum += acc[k];
                for (; p < pEnd; p++) sum += va[p] * xv[ci[p]];
                yv[i] = sum;
            }
        }

        public void MultiplyTranspose(DenseVector x, DenseVector y)
        {
            int n = N;
            Array.Clear(y.Values, 0, n);
            for (int i = 0; i < n; i++)
            {
                double xi = x.Values[i];
                int end = RowPtr[i + 1];
                for (int p = RowPtr[i]; p < end; p++)
                    y.Values[ColIdx[p]] += Values[p] * xi;
            }
        }

        public double[] GetDiagonals()
        {
            double[] d = new double[N];
            for (int i = 0; i < N; i++)
            {
                for (int p = RowPtr[i]; p < RowPtr[i + 1]; p++)
                    if (ColIdx[p] == i) { d[i] = Values[p]; break; }
            }
            return d;
        }

        public void ApplyDirichlet(int dof, double value, DenseVector rhs, double penalty = 1e16)
        {
            int pdiag = -1;
            for (int p = RowPtr[dof]; p < RowPtr[dof + 1]; p++)
                if (ColIdx[p] == dof) { pdiag = p; break; }
            if (pdiag < 0) return;
            Values[pdiag] = penalty;
            rhs.Values[dof] = penalty * value;
            DiagInv[dof] = 1.0 / penalty;
            for (int p = RowPtr[dof]; p < RowPtr[dof + 1]; p++)
                if (ColIdx[p] != dof) Values[p] = 0.0;
            for (int i = 0; i < N; i++)
            {
                if (i == dof) continue;
                for (int p = RowPtr[i]; p < RowPtr[i + 1]; p++)
                    if (ColIdx[p] == dof) Values[p] = 0.0;
            }
        }

        public void ApplyJacobi(DenseVector r, DenseVector z)
        {
            int n = N;
            int vecLen = Vector<double>.Count;
            int i = 0;
            for (; i + vecLen <= n; i += vecLen)
            {
                Vector<double> rv = new Vector<double>(r.Values, i);
                Vector<double> dv = new Vector<double>(DiagInv, i);
                (rv * dv).CopyTo(z.Values, i);
            }
            for (; i < n; i++) z.Values[i] = r.Values[i] * DiagInv[i];
        }

        public static SparseMatrixCSR FromCOOMaps(int n, Dictionary<long, double>[] maps)
        {
            Dictionary<long, double> g = new Dictionary<long, double>(n * 12);
            for (int t = 0; t < maps.Length; t++)
            {
                if (maps[t] == null) continue;
                foreach (KeyValuePair<long, double> kv in maps[t])
                {
                    double old;
                    if (g.TryGetValue(kv.Key, out old)) g[kv.Key] = old + kv.Value;
                    else g[kv.Key] = kv.Value;
                }
            }
            int nnz = g.Count;
            int[] rowptr = new int[n + 1];
            int[] colidx = new int[nnz];
            double[] vals = new double[nnz];
            foreach (KeyValuePair<long, double> kv in g) { int i = (int)(kv.Key >> 32); rowptr[i + 1]++; }
            for (int i = 0; i < n; i++) rowptr[i + 1] += rowptr[i];
            int[] cur = new int[n]; Array.Copy(rowptr, cur, n);
            foreach (KeyValuePair<long, double> kv in g)
            {
                int i = (int)(kv.Key >> 32);
                int j = (int)(kv.Key & 0xFFFFFFFFL);
                int pos = cur[i]++;
                colidx[pos] = j; vals[pos] = kv.Value;
            }
            for (int i = 0; i < n; i++)
            {
                int a = rowptr[i], b = rowptr[i + 1];
                for (int p = a + 1; p < b; p++)
                {
                    int c = colidx[p]; double vv = vals[p]; int q = p - 1;
                    while (q >= a && colidx[q] > c)
                    { colidx[q + 1] = colidx[q]; vals[q + 1] = vals[q]; q--; }
                    colidx[q + 1] = c; vals[q + 1] = vv;
                }
            }
            return new SparseMatrixCSR(n, rowptr, colidx, vals);
        }
    }
}
