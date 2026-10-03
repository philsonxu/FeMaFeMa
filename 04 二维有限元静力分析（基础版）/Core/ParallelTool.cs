using System;
using System.Threading.Tasks;

namespace FEM2D.Core
{
    public static class ParallelTool
    {
        public static int GetOptimalChunk(int total)
        {
            int procs = Environment.ProcessorCount;
            return Math.Max(1, (total + procs - 1) / procs);
        }

        public static void For(int start, int end, Action<int> body)
        {
            int count = end - start;
            if (count < 64)
            {
                for (int i = start; i < end; i++) body(i);
                return;
            }
            Parallel.For(start, end, body);
        }
    }
}
