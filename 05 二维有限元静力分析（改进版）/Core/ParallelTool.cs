// ParallelTool.cs - 多线程并行工具封装
using System;
using System.Collections.Generic;
using System.Threading.Tasks;

namespace FEM2D.Core
{
    public static class ParallelTool
    {
        public static int MaxDegreeOfParallelism = Environment.ProcessorCount;

        public static void ForEach<T>(IList<T> list, Action<T> body)
        {
            ParallelOptions opts = new ParallelOptions();
            opts.MaxDegreeOfParallelism = MaxDegreeOfParallelism;
            Parallel.ForEach(list, opts, body);
        }

        public static void For(int fromInclusive, int toExclusive, Action<int> body)
        {
            ParallelOptions opts = new ParallelOptions();
            opts.MaxDegreeOfParallelism = MaxDegreeOfParallelism;
            Parallel.For(fromInclusive, toExclusive, opts, body);
        }
    }
}
