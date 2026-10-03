using System;
using System.Collections.Generic;
using System.Numerics;
using System.Threading.Tasks;

namespace NonlinearFEM2D.Core
{
    public class SparseMatrixCSR
    {
        public int Rows, Cols;
        public int[] RowPtr;
        public int[] ColIdx;
        public double[] Values;

        public SparseMatrixCSR(int rows, int cols, int[] rowptr, int[] colidx, double[] vals)
        {
            this.Rows = rows;
            this.Cols = cols;
            this.RowPtr = rowptr;
            this.ColIdx = colidx;
            this.Values = vals;
        }

        public void MatVec(DenseVector x, DenseVector y)
        {
            int n = this.Rows;
            double[] xv = x.Values;
            double[] yv = new double[n];
            int vs = Vector<double>.Count;
            if (n < 2000)
            {
                for (int i = 0; i < n; i++)
                {
                    int s = this.RowPtr[i];
                    int e = this.RowPtr[i + 1];
                    double sum = 0.0;
                    int k = s;
                    for (; k + vs <= e; k += vs)
                    {
                        Span<double> tv = stackalloc double[vs];
                        for (int t = 0; t < vs; t++) tv[t] = xv[this.ColIdx[k + t]];
                        Vector<double> av = new Vector<double>(this.Values, k);
                        Vector<double> xv2 = new Vector<double>(tv);
                        sum += Vector.Dot(av, xv2);
                    }
                    for (; k < e; k++) sum += this.Values[k] * xv[this.ColIdx[k]];
                    yv[i] = sum;
                }
            }
            else
            {
                int cores = Math.Max(1, Environment.ProcessorCount);
                int chunk = (n + cores - 1) / cores;
                Parallel.For(0, cores, tid =>
                {
                    int begin = tid * chunk;
                    int end = Math.Min(begin + chunk, n);
                    for (int i = begin; i < end; i++)
                    {
                        int s = this.RowPtr[i];
                        int e = this.RowPtr[i + 1];
                        double sum = 0.0;
                        int k = s;
                        for (; k + vs <= e; k += vs)
                        {
                            Span<double> tv = stackalloc double[vs];
                            for (int t = 0; t < vs; t++) tv[t] = xv[this.ColIdx[k + t]];
                            Vector<double> av = new Vector<double>(this.Values, k);
                            Vector<double> xv2 = new Vector<double>(tv);
                            sum += Vector.Dot(av, xv2);
                        }
                        for (; k < e; k++) sum += this.Values[k] * xv[this.ColIdx[k]];
                        yv[i] = sum;
                    }
                });
            }
            Array.Copy(yv, y.Values, n);
        }

        public double[] GetDiagonal()
        {
            double[] d = new double[this.Rows];
            for (int i = 0; i < this.Rows; i++)
            {
                for (int k = this.RowPtr[i]; k < this.RowPtr[i + 1]; k++)
                {
                    if (this.ColIdx[k] == i)
                    {
                        d[i] = this.Values[k];
                        break;
                    }
                }
            }
            return d;
        }

        public void ApplyDirichletPenalty(int dof, double penalty)
        {
            for (int k = this.RowPtr[dof]; k < this.RowPtr[dof + 1]; k++)
            {
                if (this.ColIdx[k] == dof) { this.Values[k] = penalty; break; }
            }
        }
    }

    public class CooEntry
    {
        public int I, J;
        public double V;
        public CooEntry(int i, int j, double v) { this.I = i; this.J = j; this.V = v; }
    }

    public class CooBuilder
    {
        private readonly Dictionary<long, double> _dict;
        private readonly object _lk;
        public int Rows, Cols;

        public CooBuilder(int rows, int cols)
        {
            this.Rows = rows;
            this.Cols = cols;
            this._dict = new Dictionary<long, double>();
            this._lk = new object();
        }

        public void Add(int i, int j, double v)
        {
            if (Math.Abs(v) < 1.0e-30) return;
            long key = ((long)i << 32) | ((uint)j);
            lock (this._lk)
            {
                double old;
                if (this._dict.TryGetValue(key, out old))
                    this._dict[key] = old + v;
                else
                    this._dict[key] = v;
            }
        }

        public void AddLocal(int[] dofs, double[,] Ke)
        {
            int ne = dofs.Length;
            for (int a = 0; a < ne; a++)
            {
                int ia = dofs[a];
                for (int b = 0; b < ne; b++)
                {
                    double v = Ke[a, b];
                    if (Math.Abs(v) < 1.0e-30) continue;
                    long key = ((long)ia << 32) | ((uint)dofs[b]);
                    lock (this._lk)
                    {
                        double old;
                        if (this._dict.TryGetValue(key, out old)) this._dict[key] = old + v;
                        else this._dict[key] = v;
                    }
                }
            }
        }

        public SparseMatrixCSR Build()
        {
            int nnz = this._dict.Count;
            CooEntry[] entries = new CooEntry[nnz];
            int idx = 0;
            foreach (KeyValuePair<long, double> kv in this._dict)
            {
                long key = kv.Key;
                int i = (int)(key >> 32);
                int j = (int)(key & 0xFFFFFFFFL);
                entries[idx++] = new CooEntry(i, j, kv.Value);
            }
            Array.Sort(entries, (a, b) =>
            {
                int c = a.I.CompareTo(b.I);
                if (c != 0) return c;
                return a.J.CompareTo(b.J);
            });
            int[] rp = new int[this.Rows + 1];
            int[] ci = new int[nnz];
            double[] vs = new double[nnz];
            int currow = 0;
            rp[0] = 0;
            for (int k = 0; k < nnz; k++)
            {
                while (currow < entries[k].I) { currow++; rp[currow] = k; }
                ci[k] = entries[k].J;
                vs[k] = entries[k].V;
            }
            for (int r = currow + 1; r <= this.Rows; r++) rp[r] = nnz;
            return new SparseMatrixCSR(this.Rows, this.Cols, rp, ci, vs);
        }
    }
}
