using System;
using System.Collections.Concurrent;
using System.Threading.Tasks;

namespace NonlinearFEM2D.Core
{
    public class DenseVector
    {
        public double[] Values;
        public int Length;

        public DenseVector(int n)
        {
            this.Length = n;
            this.Values = new double[n];
        }

        public DenseVector(double[] v)
        {
            this.Length = v.Length;
            this.Values = new double[v.Length];
            Array.Copy(v, this.Values, v.Length);
        }

        public double this[int i]
        {
            get { return this.Values[i]; }
            set { this.Values[i] = value; }
        }

        public void SetZero() { Array.Clear(this.Values, 0, this.Length); }

        public void SetValue(double v)
        {
            for (int i = 0; i < this.Length; i++) this.Values[i] = v;
        }

        public void CopyFrom(DenseVector other)
        {
            Array.Copy(other.Values, this.Values, this.Length);
        }

        public DenseVector Clone()
        {
            DenseVector r = new DenseVector(this.Length);
            Array.Copy(this.Values, r.Values, this.Length);
            return r;
        }

        public double Norm2()
        {
            double s = 0.0;
            for (int i = 0; i < this.Length; i++) s += this.Values[i] * this.Values[i];
            return Math.Sqrt(s);
        }

        public double Dot(DenseVector other)
        {
            double s = 0.0;
            int n = this.Length;
            if (n < 2000)
            {
                for (int i = 0; i < n; i++) s += this.Values[i] * other.Values[i];
                return s;
            }
            object lk = new object();
            int cores = Math.Max(1, Environment.ProcessorCount);
            int chunk = (n + cores - 1) / cores;
            Parallel.For(0, cores, tid =>
            {
                int begin = tid * chunk;
                int end = Math.Min(begin + chunk, n);
                double local = 0.0;
                for (int i = begin; i < end; i++) local += this.Values[i] * other.Values[i];
                lock (lk) { s += local; }
            });
            return s;
        }

        public void Axpy(double a, DenseVector x)
        {
            int n = this.Length;
            for (int i = 0; i < n; i++) this.Values[i] += a * x.Values[i];
        }

        public void Scale(double a)
        {
            for (int i = 0; i < this.Length; i++) this.Values[i] *= a;
        }

        public double NormInf()
        {
            double m = 0.0;
            for (int i = 0; i < this.Length; i++)
            {
                double a = Math.Abs(this.Values[i]);
                if (a > m) m = a;
            }
            return m;
        }
    }
}
