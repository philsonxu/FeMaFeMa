using System;
using System.Collections.Generic;

namespace Fem2DFluid.Models
{
    /// <summary>
    /// 有限元模型：节点+单元+材料+边界条件
    /// </summary>
    public class FemModel
    {
        public List<Node> Nodes;
        public List<Element> Elements;
        public List<Material> Materials;
        /// <summary>模型标题</summary>
        public string Title;
        /// <summary>求解类型：0=势流φ（速度势），1=渗流水头H，2=稳态导热T</summary>
        public int ProblemType;

        public FemModel()
        {
            Nodes = new List<Node>();
            Elements = new List<Element>();
            Materials = new List<Material>();
            Title = "Untitled";
            ProblemType = 0;
        }

        public Node GetNode(int id)
        {
            for (int i = 0; i < Nodes.Count; i++)
            {
                if (Nodes[i].Id == id) return Nodes[i];
            }
            return null;
        }

        public Element GetElement(int id)
        {
            for (int i = 0; i < Elements.Count; i++)
            {
                if (Elements[i].Id == id) return Elements[i];
            }
            return null;
        }

        public Material GetMaterial(int id)
        {
            for (int i = 0; i < Materials.Count; i++)
            {
                if (Materials[i].Id == id) return Materials[i];
            }
            return null;
        }

        public void ComputeBounds(out double xmin, out double ymin, out double xmax, out double ymax)
        {
            xmin = double.MaxValue; ymin = double.MaxValue;
            xmax = double.MinValue; ymax = double.MinValue;
            for (int i = 0; i < Nodes.Count; i++)
            {
                Node n = Nodes[i];
                if (n.X < xmin) xmin = n.X;
                if (n.X > xmax) xmax = n.X;
                if (n.Y < ymin) ymin = n.Y;
                if (n.Y > ymax) ymax = n.Y;
            }
        }

        public void ResetResults()
        {
            for (int i = 0; i < Nodes.Count; i++)
            {
                Nodes[i].Phi = 0.0;
            }
            for (int i = 0; i < Elements.Count; i++)
            {
                Elements[i].Vx = 0.0;
                Elements[i].Vy = 0.0;
                Elements[i].Vmag = 0.0;
            }
        }
    }
}
