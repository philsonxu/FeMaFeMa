using System;

namespace Fem2DFluid.Models
{
    /// <summary>
    /// Explicit boundary edge segment (ordered endpoints). Unlike the unordered
    /// <see cref="Edge"/> used as a hash key for topology, BoundaryEdge preserves
    /// orientation to support oriented Neumann flux integration in the future.
    /// </summary>
    [Serializable]
    public class BoundaryEdge
    {
        private int _p1;
        private int _p2;
        private int _tag;
        private double _value;

        public BoundaryEdge()
        {
            _p1 = -1;
            _p2 = -1;
            _tag = 0;
            _value = 0.0;
        }

        public BoundaryEdge(int p1, int p2)
            : this()
        {
            _p1 = p1;
            _p2 = p2;
        }

        public BoundaryEdge(int p1, int p2, int tag)
            : this(p1, p2)
        {
            _tag = tag;
        }

        public int P1
        {
            get { return _p1; }
            set { _p1 = value; }
        }

        public int P2
        {
            get { return _p2; }
            set { _p2 = value; }
        }

        /// <summary>Physical group tag (1=inlet, 2=outlet, 3=wall, 4=cylinder, etc.)</summary>
        public int Tag
        {
            get { return _tag; }
            set { _tag = value; }
        }

        /// <summary>Boundary condition value (e.g. inlet velocity head, pressure head)</summary>
        public double Value
        {
            get { return _value; }
            set { _value = value; }
        }

        /// <summary>Convert to canonical unordered Edge for topology lookup.</summary>
        public Edge ToEdge()
        {
            return new Edge(_p1, _p2);
        }

        public override string ToString()
        {
            return string.Format("BoundaryEdge({0}-{1},tag={2})", _p1, _p2, _tag);
        }
    }
}
