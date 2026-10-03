using System;
using System.Globalization;

namespace Fem2DFluid.Models
{
    [Serializable]
    public class Element
    {
        private int _id;
        private int _a;
        private int _b;
        private int _c;
        private int _n4;
        private int _n5;
        private int _n6;
        private int _matId;
        private int _order;
        private double _area;
        private double[] _barycenter;
        private double _vx;
        private double _vy;
        private double _vmag;
        private double _p;

        public Element()
        {
            _id = -1;
            _a = -1;
            _b = -1;
            _c = -1;
            _n4 = -1;
            _n5 = -1;
            _n6 = -1;
            _matId = 0;
            _order = 1;
            _area = 0.0;
            _barycenter = new double[2];
            _vx = 0.0;
            _vy = 0.0;
            _vmag = 0.0;
            _p = 0.0;
        }

        public Element(int id, int a, int b, int c)
            : this()
        {
            _id = id;
            _a = a;
            _b = b;
            _c = c;
        }

        public int Id
        {
            get { return _id; }
            set { _id = value; }
        }

        public int A
        {
            get { return _a; }
            set { _a = value; }
        }

        public int B
        {
            get { return _b; }
            set { _b = value; }
        }

        public int C
        {
            get { return _c; }
            set { _c = value; }
        }

        /// <summary>Mid-edge node opposite A (between B and C), -1 if linear</summary>
        public int N4
        {
            get { return _n4; }
            set { _n4 = value; }
        }

        /// <summary>Mid-edge node opposite B (between C and A), -1 if linear</summary>
        public int N5
        {
            get { return _n5; }
            set { _n5 = value; }
        }

        /// <summary>Mid-edge node opposite C (between A and B), -1 if linear</summary>
        public int N6
        {
            get { return _n6; }
            set { _n6 = value; }
        }

        public int MatId
        {
            get { return _matId; }
            set { _matId = value; }
        }

        /// <summary>1=linear (CST, 3 nodes), 2=quadratic (LT6, 6 nodes)</summary>
        public int Order
        {
            get { return _order; }
            set { _order = value; }
        }

        public double Area
        {
            get { return _area; }
            set { _area = value; }
        }

        public double[] Barycenter
        {
            get { return _barycenter; }
            set { _barycenter = value; }
        }

        public double Vx
        {
            get { return _vx; }
            set { _vx = value; }
        }

        public double Vy
        {
            get { return _vy; }
            set { _vy = value; }
        }

        public double Vmag
        {
            get { return _vmag; }
            set { _vmag = value; }
        }

        public double P
        {
            get { return _p; }
            set { _p = value; }
        }

        public int[] NodeIds
        {
            get
            {
                if (_order == 2)
                {
                    return new int[] { _a, _b, _c, _n4, _n5, _n6 };
                }
                return new int[] { _a, _b, _c };
            }
        }

        public override string ToString()
        {
            return string.Format(CultureInfo.InvariantCulture,
                "Element({0},{1},{2},{3},area={4:F4})", _id, _a, _b, _c, _area);
        }
    }
}
