using System;

namespace Fem2DFluid.Models
{
    public enum NodeBCType
    {
        Neumann = 0,
        Dirichlet = 1
    }

    public enum NodeType
    {
        Corner = 0,
        MidEdge = 1
    }

    [Serializable]
    public class Node
    {
        private int _id;
        private double _x;
        private double _y;
        private int _bcType;
        private double _bcValue;
        private int _tag;
        private double _vx;
        private double _vy;
        private int _nodeType;
        private double _phi;
        private double _p;
        private double _u;
        private double _v;

        public Node()
        {
            _id = -1;
            _x = 0.0;
            _y = 0.0;
            _bcType = 0;
            _bcValue = 0.0;
            _tag = 0;
            _vx = 0.0;
            _vy = 0.0;
            _nodeType = 0;
            _phi = 0.0;
            _p = 0.0;
            _u = 0.0;
            _v = 0.0;
        }

        public Node(int id, double x, double y)
            : this()
        {
            _id = id;
            _x = x;
            _y = y;
        }

        public int Id
        {
            get { return _id; }
            set { _id = value; }
        }

        public double X
        {
            get { return _x; }
            set { _x = value; }
        }

        public double Y
        {
            get { return _y; }
            set { _y = value; }
        }

        /// <summary>0=Neumann, 1=Dirichlet</summary>
        public int BCType
        {
            get { return _bcType; }
            set { _bcType = value; }
        }

        public double BCValue
        {
            get { return _bcValue; }
            set { _bcValue = value; }
        }

        /// <summary>Physical group tag (0=interior, 1=inlet, 2=outlet, 3=wall, 4=cylinder)</summary>
        public int Tag
        {
            get { return _tag; }
            set { _tag = value; }
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

        /// <summary>0=corner, 1=mid-edge</summary>
        public int NodeType
        {
            get { return _nodeType; }
            set { _nodeType = value; }
        }

        /// <summary>Velocity potential (potential flow) or transported scalar</summary>
        public double Phi
        {
            get { return _phi; }
            set { _phi = value; }
        }

        /// <summary>Pressure</summary>
        public double P
        {
            get { return _p; }
            set { _p = value; }
        }

        /// <summary>Velocity x component (NS)</summary>
        public double U
        {
            get { return _u; }
            set { _u = value; }
        }

        /// <summary>Velocity y component (NS)</summary>
        public double V
        {
            get { return _v; }
            set { _v = value; }
        }

        public override string ToString()
        {
            return string.Format("Node({0},{1:F4},{2:F4})", _id, _x, _y);
        }
    }
}
