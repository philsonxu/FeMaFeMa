using System;

namespace Fem2DFluid.Models
{
    /// <summary>
    /// 二维节点
    /// </summary>
    public class Node
    {
        public int Id;
        public double X;
        public double Y;
        /// <summary>边界条件类型：0无，1第一类Dirichlet(给定值)，2第二类Neumann(给定法向导数/通量)</summary>
        public int BCType;
        /// <summary>第一类边界：压力/速度势值；第二类：法向通量</summary>
        public double BCValue;
        /// <summary>节点类型：0=普通/角点, 1=单元边中点（二次单元用）</summary>
        public int NodeType;
        /// <summary>计算结果：势函数值</summary>
        public double Phi;
        /// <summary>节点速度分量（二次单元在节点上计算）</summary>
        public double Vx;
        public double Vy;

        public Node()
        {
        }

        public Node(int id, double x, double y)
        {
            Id = id;
            X = x;
            Y = y;
            BCType = 0;
            BCValue = 0.0;
            Phi = 0.0;
            NodeType = 0;
            Vx = 0.0; Vy = 0.0;
        }

        public override string ToString()
        {
            return string.Format("Node{0}({1:F4},{2:F4})", Id, X, Y);
        }
    }
}
