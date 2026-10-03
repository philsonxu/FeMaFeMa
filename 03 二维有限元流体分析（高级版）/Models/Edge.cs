using System;

namespace Fem2DFluid.Models
{
    /// <summary>
    /// Undirected edge between two nodes. Equals and GetHashCode are symmetric
    /// so that Edge(a,b) == Edge(b,a) for use in HashSet/ Dictionary keys.
    /// </summary>
    [Serializable]
    public struct Edge : IEquatable<Edge>
    {
        private int _p1;
        private int _p2;

        public Edge(int a, int b)
        {
            // Canonical ordering: smaller id first
            if (a <= b)
            {
                _p1 = a;
                _p2 = b;
            }
            else
            {
                _p1 = b;
                _p2 = a;
            }
        }

        public int P1
        {
            get { return _p1; }
        }

        public int P2
        {
            get { return _p2; }
        }

        public bool Equals(Edge other)
        {
            return _p1 == other._p1 && _p2 == other._p2;
        }

        public override bool Equals(object obj)
        {
            if (!(obj is Edge))
            {
                return false;
            }
            Edge other = (Edge)obj;
            return Equals(other);
        }

        public override int GetHashCode()
        {
            // Symmetric hash: addition commutes, so Edge(a,b) and Edge(b,a) collide intentionally
            int h1 = _p1.GetHashCode();
            int h2 = _p2.GetHashCode();
            return (int)(h1 ^ h2 ^ ((h1 + h2) * 2654435761));
        }

        public static bool operator ==(Edge left, Edge right)
        {
            return left.Equals(right);
        }

        public static bool operator !=(Edge left, Edge right)
        {
            return !left.Equals(right);
        }

        public override string ToString()
        {
            return string.Format("Edge({0}-{1})", _p1, _p2);
        }
    }
}
