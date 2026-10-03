using System;
using System.Collections.Generic;
using System.Numerics;
using System.Threading;
using System.Threading.Tasks;

namespace ModalFEM2D.Core
{

    /// <summary>
    /// CSR sparse matrix for linear-algebra operations. MatVec is multi-threaded + SIMD.
    /// </summary>
    public sealed class SparseMatrixCSR
    {
        public readonly int NumRows;
        public readonly int NumCols;
        public readonly int[] RowPtr;   // length N+1
        public readonly int[] ColIdx;   // length NNZ
        public readonly double[] Values;// length NNZ
        public double[] DiagInv;        // cached Jacobi preconditioner

        public SparseMatrixCSR(int nRows, int nCols, int[] rowPtr, int[] colIdx, double[] values)
        {
            NumRows = nRows;
            NumCols = nCols;
            RowPtr = rowPtr;
            ColIdx = colIdx;
            Values = values;
        }

        public static SparseMatrixCSR FromCoo(int nRows, int nCols, List<Tuple<int, int, double>> triples)
        {
            int nnz = triples.Count;
            int[] rp = new int[nRows + 1];
            int[] ci = new int[nnz];
            double[] va = new double[nnz];
            int[] cnt = new int[nRows];
            for (int k = 0; k < nnz; k++) cnt[triples[k].Item1]++;
            int pos = 0;
            for (int i = 0; i < nRows; i++) { rp[i] = pos; pos += cnt[i]; }
            rp[nRows] = pos;
            Array.Copy(rp, cnt, nRows);
            for (int k = 0; k < nnz; k++)
            {
                int r = triples[k].Item1;
                int dst = cnt[r]++;
                ci[dst] = triples[k].Item2;
                va[dst] = triples[k].Item3;
            }
            // sort each row by column
            for (int i = 0; i < nRows; i++)
            {
                int s = rp[i], e = rp[i + 1];
                for (int a = s + 1; a < e; a++)
                {
                    int j = ci[a]; double v = va[a]; int b = a - 1;
                    while (b >= s && ci[b] > j) { ci[b + 1] = ci[b]; va[b + 1] = va[b]; b--; }
                    ci[b + 1] = j; va[b + 1] = v;
                }
            }
            // merge duplicates
            int w = 0;
            for (int i = 0; i < nRows; i++)
            {
                int s = rp[i]; int e = rp[i + 1]; int newStart = w;
                for (int k = s; k < e; k++)
                {
                    if (k > s && ci[k] == ci[w - 1]) { va[w - 1] += va[k]; }
                    else { ci[w] = ci[k]; va[w] = va[k]; w++; }
                }
                rp[i] = newStart;
            }
            rp[nRows] = w;
            Array.Resize(ref ci, w);
            Array.Resize(ref va, w);
            SparseMatrixCSR mat = new SparseMatrixCSR(nRows, nCols, rp, ci, va);
            mat.BuildDiagInv();
            return mat;
        }

        public void BuildDiagInv()
        {
            DiagInv = new double[NumRows];
            for (int i = 0; i < NumRows; i++)
            {
                double d = 0.0;
                for (int k = RowPtr[i]; k < RowPtr[i + 1]; k++)
                    if (ColIdx[k] == i) { d = Values[k]; break; }
                DiagInv[i] = Math.Abs(d) > 1e-30 ? 1.0 / d : 1.0;
            }
        }

        /// <summary>y = this * x (clears y first). Multi-threaded + SIMD.</summary>
        public void MatVec(DenseVector x, DenseVector y)
        {
            if (x.Length != NumCols || y.Length != NumRows) throw new ArgumentException("dim");
            Array.Clear(y.Values, 0, NumRows);
            int vlen = Vector<double>.Count;
            if (NumRows < 2000 || Environment.ProcessorCount == 1)
            {
                MatVecSegment(0, NumRows, x, y, vlen);
            }
            else
            {
                int threads = Environment.ProcessorCount;
                int chunk = (NumRows + threads - 1) / threads;
                WaitHandle[] handles = new WaitHandle[threads];
                for (int t = 0; t < threads; t++)
                {
                    int s = t * chunk; int e = Math.Min(s + chunk, NumRows);
                    if (s >= e) break;
                    ManualResetEvent mre = new ManualResetEvent(false);
                    handles[t] = mre;
                    ThreadPool.QueueUserWorkItem(_ =>
                    {
                        MatVecSegment(s, e, x, y, vlen);
                        mre.Set();
                    });
                }
                WaitHandle.WaitAll(handles);
                for (int t = 0; t < handles.Length; t++) if (handles[t] != null) handles[t].Dispose();
            }
        }

        private void MatVecSegment(int s, int e, DenseVector x, DenseVector y, int vlen)
        {
            double[] xv = x.Values;
            Span<double> xs = stackalloc double[vlen];
            for (int i = s; i < e; i++)
            {
                int rs = RowPtr[i], re = RowPtr[i + 1]; double sum = 0.0; int k;
                for (k = rs; k + vlen <= re; k += vlen)
                {
                    for (int j = 0; j < vlen; j++) xs[j] = xv[ColIdx[k + j]];
                    Vector<double> a = new Vector<double>(Values, k);
                    Vector<double> b = new Vector<double>(xs);
                    sum += Vector.Dot(a, b);
                }
                for (; k < re; k++) sum += Values[k] * xv[ColIdx[k]];
                y.Values[i] = sum;
            }
        }

        /// <summary>Apply Dirichlet: set row/col i to zero, diagonal = big.</summary>
        public void ApplyDirichlet(int i, double penalty)
        {
            for (int k = RowPtr[i]; k < RowPtr[i + 1]; k++)
            {
                Values[k] = (ColIdx[k] == i) ? penalty : 0.0;
            }
            if (DiagInv != null) DiagInv[i] = 1.0 / penalty;
        }

        /// <summary>Extract diagonal.</summary>
        public DenseVector GetDiagonal()
        {
            DenseVector d = new DenseVector(NumRows);
            for (int i = 0; i < NumRows; i++)
                for (int k = RowPtr[i]; k < RowPtr[i + 1]; k++)
                    if (ColIdx[k] == i) { d.Values[i] = Values[k]; break; }
            return d;
        }
    }

    /// <summary>
    /// Thread-safe COO triple assembler. Supports per-thread local builders merged at end.
    /// </summary>
    public sealed class CooAssembler
    {
        public readonly int NumRows;
        public readonly int NumCols;
        private readonly object _lock = new object();
        private readonly List<Tuple<int, int, double>> _triples = new List<Tuple<int, int, double>>();

        public CooAssembler(int nRows, int nCols) { NumRows = nRows; NumCols = nCols; }

        public void Add(int r, int c, double v)
        {
            if (r < 0 || c < 0 || r >= NumRows || c >= NumCols) return;
            lock (_lock) { _triples.Add(Tuple.Create(r, c, v)); }
        }

        public void AddLocal(int[] rmap, int[] cmap, double[,] ke)
        {
            int nr = rmap.Length, nc = cmap.Length;
            lock (_lock)
            {
                for (int i = 0; i < nr; i++)
                {
                    int rr = rmap[i];
                    if (rr < 0 || rr >= NumRows) continue;
                    for (int j = 0; j < nc; j++)
                    {
                        int cc = cmap[j];
                        if (cc < 0 || cc >= NumCols) continue;
                        _triples.Add(Tuple.Create(rr, cc, ke[i, j]));
                    }
                }
            }
        }

        public void AddSymmLocal(int[] dofMap, double[,] ke)
        {
            int n = dofMap.Length;
            lock (_lock)
            {
                for (int i = 0; i < n; i++)
                {
                    int rr = dofMap[i];
                    if (rr < 0 || rr >= NumRows) continue;
                    for (int j = 0; j < n; j++)
                    {
                        int cc = dofMap[j];
                        if (cc < 0 || cc >= NumCols) continue;
                        _triples.Add(Tuple.Create(rr, cc, ke[i, j]));
                    }
                }
            }
        }

        public void AddMassLocal(int[] dofMap, double[,] me) { AddSymmLocal(dofMap, me); }

        public SparseMatrixCSR Build()
        {
            return SparseMatrixCSR.FromCoo(NumRows, NumCols, _triples);
        }

        public void Merge(CooAssembler other)
        {
            lock (_lock) { _triples.AddRange(other._triples); }
        }
    }
}
