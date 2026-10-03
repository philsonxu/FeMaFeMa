using System;
using System.Collections.Generic;
using System.Threading.Tasks;

namespace FEM2D.Core
{
    /// <summary>
    /// CSR 格式稀疏矩阵。支持 COO 装配、并行化 Mat-Vec。
    /// </summary>
    public sealed class SparseMatrixCSR
    {
        public int N { get; private set; }
        public int[] RowPtr { get; private set; }
        public int[] ColIdx { get; private set; }
        public double[] Values { get; private set; }

        private SparseMatrixCSR() { }

        /// <summary>
        /// COO 三元组构造器（线程安全累加）
        /// </summary>
        public sealed class CooBuilder
        {
            private readonly object _lock = new object();
            private readonly List<(int r, int c, double v)> _entries = new List<(int, int, double)>();
            private readonly int _n;

            public CooBuilder(int n) { _n = n; }

            public void Add(int r, int c, double v)
            {
                if (v == 0.0) return;
                lock (_lock)
                {
                    _entries.Add((r, c, v));
                }
            }

            public void AddLocal(int[] dofMap, double[,] ke)
            {
                int m = dofMap.Length;
                lock (_lock)
                {
                    for (int i = 0; i < m; i++)
                    {
                        int gi = dofMap[i];
                        for (int j = 0; j < m; j++)
                        {
                            int gj = dofMap[j];
                            double v = ke[i, j];
                            if (v != 0.0) _entries.Add((gi, gj, v));
                        }
                    }
                }
            }

            public SparseMatrixCSR Build(bool parallelSort = true)
            {
                int n = _n;
                int nnz = _entries.Count;

                int[] rowPtr = new int[n + 1];
                int[] colIdx = new int[nnz];
                double[] vals = new double[nnz];

                // 统计每行非零元个数
                foreach ((int r, int c, double v) in _entries)
                {
                    rowPtr[r + 1]++;
                }
                for (int i = 0; i < n; i++)
                {
                    rowPtr[i + 1] += rowPtr[i];
                }

                int[] cursor = new int[n];
                Array.Copy(rowPtr, cursor, n);

                // 列索引数组填装
                foreach ((int r, int c, double v) in _entries)
                {
                    int pos = cursor[r]++;
                    colIdx[pos] = c;
                    vals[pos] = v;
                }

                // 每行内按列号排序并累加重复项
                int newNnz = 0;
                if (parallelSort && n > 500)
                {
                    int procs = Environment.ProcessorCount;
                    int chunk = (n + procs - 1) / procs;
                    int[] segStart = new int[procs + 1];
                    for (int p = 0; p < procs; p++)
                    {
                        int r0 = Math.Min(p * chunk, n);
                        segStart[p] = rowPtr[r0];
                    }
                    segStart[procs] = nnz;

                    Parallel.For(0, procs, p =>
                    {
                        int r0 = Math.Min(p * chunk, n);
                        int r1 = Math.Min((p + 1) * chunk, n);
                        for (int r = r0; r < r1; r++)
                        {
                            SortAndMergeRow(rowPtr, colIdx, vals, r);
                        }
                    });

                    // 压缩
                    int[] newCol = new int[nnz];
                    double[] newVal = new double[nnz];
                    int[] newRowPtr = new int[n + 1];
                    int k = 0;
                    for (int r = 0; r < n; r++)
                    {
                        newRowPtr[r] = k;
                        int start = rowPtr[r];
                        int end = rowPtr[r + 1];
                        int ci = colIdx[start];
                        double vs = vals[start];
                        for (int i = start + 1; i < end; i++)
                        {
                            if (colIdx[i] == ci) { vs += vals[i]; }
                            else
                            {
                                if (vs != 0.0) { newCol[k] = ci; newVal[k] = vs; k++; }
                                ci = colIdx[i]; vs = vals[i];
                            }
                        }
                        if (end > start && vs != 0.0) { newCol[k] = ci; newVal[k] = vs; k++; }
                    }
                    newRowPtr[n] = k;
                    int[] finalCol = new int[k];
                    double[] finalVal = new double[k];
                    Array.Copy(newCol, 0, finalCol, 0, k);
                    Array.Copy(newVal, 0, finalVal, 0, k);

                    return new SparseMatrixCSR { N = n, RowPtr = newRowPtr, ColIdx = finalCol, Values = finalVal };
                }
                else
                {
                    // 串行版本
                    int[] newCol = new int[nnz];
                    double[] newVal = new double[nnz];
                    int[] newRowPtr = new int[n + 1];
                    int k = 0;
                    for (int r = 0; r < n; r++)
                    {
                        newRowPtr[r] = k;
                        int start = rowPtr[r];
                        int end = rowPtr[r + 1];
                        if (end <= start) continue;
                        // 简易排序
                        for (int i = start; i < end - 1; i++)
                        {
                            for (int j = i + 1; j < end; j++)
                            {
                                if (colIdx[j] < colIdx[i])
                                {
                                    int tc = colIdx[i]; colIdx[i] = colIdx[j]; colIdx[j] = tc;
                                    double tv = vals[i]; vals[i] = vals[j]; vals[j] = tv;
                                }
                            }
                        }
                        int ci = colIdx[start]; double vs = vals[start];
                        for (int i = start + 1; i < end; i++)
                        {
                            if (colIdx[i] == ci) { vs += vals[i]; }
                            else
                            {
                                if (vs != 0.0) { newCol[k] = ci; newVal[k] = vs; k++; }
                                ci = colIdx[i]; vs = vals[i];
                            }
                        }
                        if (vs != 0.0) { newCol[k] = ci; newVal[k] = vs; k++; }
                    }
                    newRowPtr[n] = k;
                    int[] finalCol = new int[k];
                    double[] finalVal = new double[k];
                    Array.Copy(newCol, 0, finalCol, 0, k);
                    Array.Copy(newVal, 0, finalVal, 0, k);
                    return new SparseMatrixCSR { N = n, RowPtr = newRowPtr, ColIdx = finalCol, Values = finalVal };
                }
            }

            private static void SortAndMergeRow(int[] rowPtr, int[] colIdx, double[] vals, int r)
            {
                int start = rowPtr[r];
                int end = rowPtr[r + 1];
                int len = end - start;
                if (len <= 1) return;
                Array.Sort(colIdx, vals, start, len);
            }
        }

        public void MatVec(DenseVector x, DenseVector y)
        {
            int n = N;
            int[] rp = RowPtr;
            int[] ci = ColIdx;
            double[] va = Values;
            double[] xv = x.Values;
            double[] yv = y.Values;

            int procs = Environment.ProcessorCount;
            if (n < 2000 || procs <= 1)
            {
                for (int i = 0; i < n; i++)
                {
                    double s = 0.0;
                    int e = rp[i + 1];
                    for (int k = rp[i]; k < e; k++)
                    {
                        s += va[k] * xv[ci[k]];
                    }
                    yv[i] = s;
                }
                return;
            }
            Parallel.For(0, n, i =>
            {
                double s = 0.0;
                int e = rp[i + 1];
                for (int k = rp[i]; k < e; k++)
                {
                    s += va[k] * xv[ci[k]];
                }
                yv[i] = s;
            });
        }

        /// <summary>
        /// 雅可比对角预处理矩阵逆（M^{-1}）
        /// </summary>
        public double[] DiagonalInverse()
        {
            double[] diag = new double[N];
            for (int i = 0; i < N; i++)
            {
                double d = 0.0;
                for (int k = RowPtr[i]; k < RowPtr[i + 1]; k++)
                {
                    if (ColIdx[k] == i) { d = Values[k]; break; }
                }
                diag[i] = (Math.Abs(d) > 1e-30) ? 1.0 / d : 1.0;
            }
            return diag;
        }
    }
}
