using System;
using System.Collections.Generic;

namespace Fem2DFluid.Models
{
    /// <summary>
    /// 三角单元
    /// </summary>
    /// <summary>三角单元：线性CST(3节点)或二次LT6(6节点，N4边1-2中点、N5边2-3中点、N6边3-1中点)</summary>
    public class Element
    {
        public int Id;
        public int N1;
        public int N2;
        public int N3;
        /// <summary>二次边中点：边1-2（0表示线性单元无中点）</summary>
        public int N4;
        /// <summary>二次边中点：边2-3</summary>
        public int N5;
        /// <summary>二次边中点：边3-1</summary>
        public int N6;
        /// <summary>单元类型：0=线性CST, 1=二次LT6</summary>
        public int Order;
        /// <summary>单元材料编号</summary>
        public int MatId;
        /// <summary>单元中心X</summary>
        public double Cx;
        /// <summary>单元中心Y</summary>
        public double Cy;
        /// <summary>单元面积</summary>
        public double Area;

        /// <summary>单元内速度分量（由单元常梯度计算——CST版本；LT6在节点上计算速度，此处作为单元平均）</summary>
        public double Vx;
        public double Vy;
        /// <summary>速度大小</summary>
        public double Vmag;

        public Element()
        {
        }

        public Element(int id, int n1, int n2, int n3, int matId)
        {
            Id = id;
            N1 = n1;
            N2 = n2;
            N3 = n3;
            N4 = 0; N5 = 0; N6 = 0; Order = 0;
            MatId = matId;
        }

        public List<int> NodeList()
        {
            List<int> list = new List<int>();
            list.Add(N1); list.Add(N2); list.Add(N3);
            if (Order >= 1)
            {
                if (N4 > 0) list.Add(N4);
                if (N5 > 0) list.Add(N5);
                if (N6 > 0) list.Add(N6);
            }
            return list;
        }

        /// <summary>返回单元的6个节点索引（LT6），线性单元用0占位</summary>
        public int[] NodeIds6()
        {
            return new int[] { N1, N2, N3, N4, N5, N6 };
        }
    }
}
