using System;
using System.Threading;

namespace MultiPhysicsFEM2D.Core
{
    public class DenseVector
    {
        public double[] Values;
        public int Length { get { return Values != null ? Values.Length : 0; } }

        public DenseVector(int n)
        {
            Values = new double[n];
        }

        public DenseVector(double[] src)
        {
            Values = new double[src.Length];
            Array.Copy(src, Values, src.Length);
        }

        public DenseVector Clone()
        {
            DenseVector r = new DenseVector(Length);
            Array.Copy(Values, r.Values, Length);
            return r;
        }

        public void SetZero() { Array.Clear(Values, 0, Length); }
        public void SetValue(double v) { for (int i = 0; i < Length; i++) Values[i] = v; }
        public void CopyFrom(DenseVector o) { Array.Copy(o.Values, Values, Math.Min(Length, o.Length)); }
        public void Axpy(double a, DenseVector x) { for (int i = 0; i < Length; i++) Values[i] += a * x.Values[i]; }
        public void Scale(double s) { for (int i = 0; i < Length; i++) Values[i] *= s; }
        public void Add(DenseVector x) { for (int i = 0; i < Length; i++) Values[i] += x.Values[i]; }
        public void Subtract(DenseVector x) { for (int i = 0; i < Length; i++) Values[i] -= x.Values[i]; }

        public double Dot(DenseVector other)
        {
            int n = Length;
            if (n < 2000)
            {
                double s = 0.0;
                for (int i = 0; i < n; i++) s += Values[i] * other.Values[i];
                return s;
            }
            int threads = Environment.ProcessorCount;
            double[] parts = new double[threads];
            int block = (n + threads - 1) / threads;
            WaitHandle[] whs = new WaitHandle[threads];
            for (int t = 0; t < threads; t++)
            {
                int start = t * block, end = Math.Min(start + block, n);
                int tt = t;
                ManualResetEvent mre = new ManualResetEvent(false);
                whs[t] = mre;
                ThreadPool.QueueUserWorkItem(delegate
                {
                    double s = 0.0;
                    for (int i = start; i < end; i++) s += Values[i] * other.Values[i];
                    parts[tt] = s;
                    mre.Set();
                });
            }
            WaitHandle.WaitAll(whs);
            double total = 0; for (int t = 0; t < threads; t++) total += parts[t];
            return total;
        }

        public double Norm2() { return Math.Sqrt(Dot(this)); }
        public double NormInf()
        {
            double m = 0;
            for (int i = 0; i < Length; i++) { double a = Math.Abs(Values[i]); if (a > m) m = a; }
            return m;
        }
    }
}
