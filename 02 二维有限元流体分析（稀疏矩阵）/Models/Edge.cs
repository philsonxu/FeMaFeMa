using System;
using System.Collections.Generic;

namespace Fem2DFluid.Models
{
    /// <summary>
    /// 边（用于边界条件识别、网格剖分边界）
    /// </summary>
    public class Edge
    {
        public int P1;
        public int P2;
        public int BCType;
        public double BCValue;

        public Edge(int p1, int p2)
        {
            // 保留调用方传入的方向，不强制排序；
            // Equals/GetHashCode 使用无序比较，保证哈希表中 (a,b) 与 (b,a) 视为同一条边。
            P1 = p1;
            P2 = p2;
            BCType = 0;
            BCValue = 0.0;
        }

        public Edge(int p1, int p2, int bcType, double bcValue) : this(p1, p2)
        {
            BCType = bcType;
            BCValue = bcValue;
        }

        public override bool Equals(object obj)
        {
            Edge other = obj as Edge;
            if (other == null) return false;
            return P1 == other.P1 && P2 == other.P2;
        }

        public override int GetHashCode()
        {
            return P1 ^ (P2 * 179424673);
        }
    }
}
