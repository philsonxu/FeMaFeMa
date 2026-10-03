using System;
using System.Collections.Generic;

namespace Fem2DFluid.Models
{
    /// <summary>
    /// 三角单元
    /// </summary>
    public class Element
    {
        public int Id;
        public int N1;
        public int N2;
        public int N3;
        /// <summary>单元材料编号</summary>
        public int MatId;
        /// <summary>单元中心X</summary>
        public double Cx;
        /// <summary>单元中心Y</summary>
        public double Cy;
        /// <summary>单元面积</summary>
        public double Area;

        /// <summary>单元内速度分量（由单元常梯度计算）</summary>
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
            MatId = matId;
        }

        public List<int> NodeList()
        {
            List<int> list = new List<int>();
            list.Add(N1);
            list.Add(N2);
            list.Add(N3);
            return list;
        }
    }
}
