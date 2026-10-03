using System;
using System.Threading.Tasks;

namespace FEM2D.Core
{
    /// <summary>
    /// 稠密向量，封装 BLAS-1 级基本运算
    /// </summary>
    public sealed class DenseVector
    {
        private readonly double[] _values;

        public DenseVector(int n)
        {
            _values = new double[n];
        }

        public DenseVector(double[] values)
        {
            _values = (double[])values.Clone();
        }

        public int Length => _values.Length;

        public double this[int i]
        {
            get => _values[i];
            set => _values[i] = value;
        }

        public double[] Values => _values;

        public DenseVector Clone()
        {
            return new DenseVector((double[])_values.Clone());
        }

        public void Zero()
        {
            Array.Clear(_values, 0, _values.Length);
        }

        public double Dot(DenseVector other)
        {
            double sum = 0.0;
            int n = _values.Length;
            for (int i = 0; i < n; i++)
            {
                sum += _values[i] * other._values[i];
            }
            return sum;
        }

        public double Norm2()
        {
            return Math.Sqrt(Dot(this));
        }

        public void Add(DenseVector a, DenseVector b)
        {
            int n = _values.Length;
            for (int i = 0; i < n; i++)
            {
                _values[i] = a._values[i] + b._values[i];
            }
        }

        public void AddScaled(DenseVector v, double alpha)
        {
            int n = _values.Length;
            for (int i = 0; i < n; i++)
            {
                _values[i] += alpha * v._values[i];
            }
        }

        public void Scale(double alpha)
        {
            int n = _values.Length;
            for (int i = 0; i < n; i++)
            {
                _values[i] *= alpha;
            }
        }

        public void CopyFrom(DenseVector src)
        {
            Array.Copy(src._values, _values, _values.Length);
        }
    }
}
