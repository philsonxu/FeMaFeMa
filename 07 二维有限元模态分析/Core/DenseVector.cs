using System;
using System.Threading.Tasks;
using System.Collections.Concurrent;

namespace ModalFEM2D.Core
{

    /// <summary>
    /// Dense double vector with BLAS-1 operations (parallel Dot/Axpy).
    /// </summary>
    public sealed class DenseVector
    {
        public readonly double[] Values;
        public readonly int Length;

        public DenseVector(int n)
        {
            Length = n;
            Values = new double[n];
        }

        public DenseVector(double[] v)
        {
            Values = v;
            Length = v.Length;
        }

        public double this[int i]
        {
            get { return Values[i]; }
            set { Values[i] = value; }
        }

        public void Clear()
        {
            Array.Clear(Values, 0, Length);
        }

        public void SetAll(double v)
        {
            for (int i = 0; i < Length; i++) Values[i] = v;
        }

        public void CopyFrom(DenseVector other)
        {
            Array.Copy(other.Values, Values, Length);
        }

        public void Scale(double a)
        {
            for (int i = 0; i < Length; i++) Values[i] *= a;
        }

        /// <summary>this = this + a*v</summary>
        public void Axpy(double a, DenseVector v)
        {
            if (Length < 2000)
            {
                for (int i = 0; i < Length; i++) Values[i] += a * v.Values[i];
            }
            else
            {
                Parallel.For(0, Length, i => { Values[i] += a * v.Values[i]; });
            }
        }

        /// <summary>this = a*this + b*v</summary>
        public void ScaleAxpy(double a, double b, DenseVector v)
        {
            for (int i = 0; i < Length; i++) Values[i] = a * Values[i] + b * v.Values[i];
        }

        public double Dot(DenseVector other)
        {
            double s = 0.0;
            if (Length < 4000)
            {
                for (int i = 0; i < Length; i++) s += Values[i] * other.Values[i];
            }
            else
            {
                object locker = new object();
                Parallel.ForEach(Partitioner.Create(0, Length), range =>
                {
                    double ls = 0.0;
                    for (int i = range.Item1; i < range.Item2; i++) ls += Values[i] * other.Values[i];
                    lock (locker) s += ls;
                });
            }
            return s;
        }

        public double Norm2()
        {
            return Math.Sqrt(Dot(this));
        }

        public double NormInf()
        {
            double m = 0.0;
            for (int i = 0; i < Length; i++)
            {
                double a = Math.Abs(Values[i]);
                if (a > m) m = a;
            }
            return m;
        }

        public DenseVector Clone()
        {
            DenseVector c = new DenseVector(Length);
            Array.Copy(Values, c.Values, Length);
            return c;
        }
    }
}
