using System;
using System.Collections.Generic;

namespace FEM2D.Mesh
{
    public enum ElementType
    {
        CST3,
        LT6,
        Q4
    }

    public sealed class FEMesh
    {
        public List<double[]> Nodes { get; } = new List<double[]>();
        public List<int[]> Elements { get; } = new List<int[]>();
        public List<ElementType> ElementTypes { get; } = new List<ElementType>();
        public List<int> MaterialIDs { get; } = new List<int>();

        // Dirichlet 约束：(nodeLocalIndex, dofLocalIndex, value)  dof: 0=ux(or u),1=uy(or v)
        public List<(int node, int dof, double val)> DirichletBCs { get; } = new List<(int, int, double)>();
        // Neumann：边荷载 / 集中力 (nodeLocalIndex, dofLocalIndex, value)
        public List<(int node, int dof, double val)> ConcentratedForces { get; } = new List<(int, int, double)>();
        // 边压力（由 reader 转成等效节点力）：保留元信息
        public List<(int e1, int e2, double pressure)> EdgePressures { get; } = new List<(int, int, double)>();

        public Dictionary<string, object> Metadata { get; } = new Dictionary<string, object>();

        public int NumNodes => Nodes.Count;
        public int NumElements => Elements.Count;

        public int NDOF
        {
            get
            {
                int dofPerNode = 2;
                if (Metadata.ContainsKey("ProblemType") && (string)Metadata["ProblemType"] == "NS")
                    dofPerNode = 3; // u,v,p
                return Nodes.Count * dofPerNode;
            }
        }

        public double X(int n) => Nodes[n][0];
        public double Y(int n) => Nodes[n][1];
    }
}
