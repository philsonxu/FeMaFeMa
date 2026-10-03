using System;

namespace FEM2D.Core
{
    public sealed class DenseVector
    {
        public double[] Values;
        public int Length;

        public DenseVector(int n) { Length = n; Values = new double[n]; }
        public DenseVector(double[] data) { Length = data.Length; Values = (double[])data.Clone(); }

        public double this[int i]
        {
            get { return Values[i]; }
            set { Values[i] = value; }
        }

        public void SetZero() { Array.Clear(Values, 0, Length); }
        public void CopyFrom(DenseVector other) { Array.Copy(other.Values, Values, Length); }
        public DenseVector Clone() { return new DenseVector(Values); }

        public double Dot(DenseVector other)
        {
            double s = 0.0;
            for (int i = 0; i < Length; i++) s += Values[i] * other.Values[i];
            return s;
        }
        public double Norm2() { return Math.Sqrt(Dot(this)); }

        public void Axpy(double alpha, DenseVector x)
        { for (int i = 0; i < Length; i++) Values[i] += alpha * x.Values[i]; }

        public void Scale(double alpha)
        { for (int i = 0; i < Length; i++) Values[i] *= alpha; }

        public static DenseVector operator +(DenseVector a, DenseVector b)
        {
            DenseVector r = new DenseVector(a.Length);
            for (int i = 0; i < a.Length; i++) r.Values[i] = a.Values[i] + b.Values[i];
            return r;
        }
        public static DenseVector operator -(DenseVector a, DenseVector b)
        {
            DenseVector r = new DenseVector(a.Length);
            for (int i = 0; i < a.Length; i++) r.Values[i] = a.Values[i] - b.Values[i];
            return r;
        }
        public static DenseVector operator *(double a, DenseVector x)
        {
            DenseVector r = new DenseVector(x.Length);
            for (int i = 0; i < x.Length; i++) r.Values[i] = a * x.Values[i];
            return r;
        }
    }
}
