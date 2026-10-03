// CooAssembler.cs - 线程安全的 COO 三元组装配器（多单元并行累加单元矩阵）
using System;
using System.Collections.Generic;
using System.Threading;

namespace FEM2D.Core
{
    public sealed class CooAssembler
    {
        private readonly Dictionary<long, double>[] _maps;

        public CooAssembler()
        {
            int threads = Math.Min(Environment.ProcessorCount, 16);
            if (threads < 1) threads = 1;
            _maps = new Dictionary<long, double>[threads];
            for (int t = 0; t < threads; t++) _maps[t] = new Dictionary<long, double>();
        }

        public void Add(int i, int j, double v)
        {
            int tid = Thread.CurrentThread.ManagedThreadId % _maps.Length;
            long key = ((long)i << 32) | (uint)j;
            double old;
            if (_maps[tid].TryGetValue(key, out old)) _maps[tid][key] = old + v;
            else _maps[tid][key] = v;
        }

        public void AddBlock(int[] dofMap, double[,] ke)
        {
            int n = dofMap.Length;
            int tid = Thread.CurrentThread.ManagedThreadId % _maps.Length;
            Dictionary<long, double> M = _maps[tid];
            for (int a = 0; a < n; a++)
            {
                int ia = dofMap[a];
                for (int b = 0; b < n; b++)
                {
                    double v = ke[a, b];
                    if (v == 0.0) continue;
                    long key = ((long)ia << 32) | (uint)dofMap[b];
                    double old;
                    if (M.TryGetValue(key, out old)) M[key] = old + v;
                    else M[key] = v;
                }
            }
        }

        public SparseMatrixCSR ToCSR(int n)
        {
            return SparseMatrixCSR.FromCOOMaps(n, _maps);
        }
    }
}
