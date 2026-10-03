// FEMesh.cs - 二维有限元网格容器
using System;
using System.Collections.Generic;
using FEM2D.Elements;

namespace FEM2D.Mesh
{
    public sealed class FEMesh
    {
        // ---- 基本拓扑 ----
        public int NumNodes;
        public int NumElements;
        public ElementType ElementType = ElementType.CST3;
        public double[] X;
        public double[] Y;
        public List<FiniteElement> Elements = new List<FiniteElement>();

        // ---- 边界条件（节点编号、dof(0=u/1=v)、数值） ----
        public List<Tuple<int, int, double>> DirichletBCs = new List<Tuple<int, int, double>>();
        public List<Tuple<int, int, double>> NodalForces = new List<Tuple<int, int, double>>();
        public List<Tuple<int, int, double>> VelocityBCs = new List<Tuple<int, int, double>>();

        // ---- 边界包围盒 ----
        public double XMin, XMax, YMin, YMax;

        // ---- 派生快访视图（由 BuildConnectivity 填充） ----
        public double[,] Coords;        // [nN,2]
        public int[,] Connectivity;     // [nE, nen]
        public int[] NodeBC;            // [nN] 位掩码：bit0=u固定，bit1=v固定
        public double[] DirichletU;     // [nN] 速度/位移 u 的 Dirichlet 值
        public double[] DirichletV;     // [nN] v 的 Dirichlet 值

        // ---- 结果数组（由求解器填充） ----
        public double[] DisplacementU;
        public double[] DisplacementV;
        public double[] VelocityU;
        public double[] VelocityV;
        public double[] Pressure;
        public double[] StressXX;   // 单元级
        public double[] StressYY;
        public double[] StressXY;
        public double[] VonMises;

        public void UpdateBounds()
        {
            XMin = double.MaxValue; XMax = double.MinValue;
            YMin = double.MaxValue; YMax = double.MinValue;
            for (int i = 0; i < NumNodes; i++)
            {
                if (X[i] < XMin) XMin = X[i];
                if (X[i] > XMax) XMax = X[i];
                if (Y[i] < YMin) YMin = Y[i];
                if (Y[i] > YMax) YMax = Y[i];
            }
        }

        public void AssignNodeCoordinatesToElements()
        {
            foreach (FiniteElement e in Elements)
            {
                int nen = e.NodesPerElement;
                e.X = new double[nen];
                e.Y = new double[nen];
                for (int k = 0; k < nen; k++)
                {
                    int gid = e.NodeIds[k];
                    e.X[k] = X[gid];
                    e.Y[k] = Y[gid];
                }
            }
        }

        /// <summary>
        /// 构造 Coords / Connectivity / NodeBC / DirichletUV 快访数组，
        /// 必须在节点/单元/BC 全部添加完成后调用。
        /// </summary>
        public void BuildConnectivity()
        {
            Coords = new double[NumNodes, 2];
            for (int i = 0; i < NumNodes; i++) { Coords[i, 0] = X[i]; Coords[i, 1] = Y[i]; }

            if (NumElements > 0) ElementType = Elements[0].Type;
            int nen = NodesPerElement(ElementType);
            Connectivity = new int[NumElements, nen];
            for (int e = 0; e < NumElements; e++)
            {
                FiniteElement el = Elements[e];
                for (int k = 0; k < nen; k++) Connectivity[e, k] = el.NodeIds[k];
            }

            NodeBC = new int[NumNodes];
            DirichletU = new double[NumNodes];
            DirichletV = new double[NumNodes];
            foreach (Tuple<int, int, double> bc in DirichletBCs)
            {
                if (bc.Item2 == 0) { NodeBC[bc.Item1] |= 0b01; DirichletU[bc.Item1] = bc.Item3; }
                else { NodeBC[bc.Item1] |= 0b10; DirichletV[bc.Item1] = bc.Item3; }
            }
            foreach (Tuple<int, int, double> bc in VelocityBCs)
            {
                if (bc.Item2 == 0) { NodeBC[bc.Item1] |= 0b01; DirichletU[bc.Item1] = bc.Item3; }
                else { NodeBC[bc.Item1] |= 0b10; DirichletV[bc.Item1] = bc.Item3; }
            }

            // 初始化结果数组
            DisplacementU = new double[NumNodes];
            DisplacementV = new double[NumNodes];
            VelocityU = new double[NumNodes];
            VelocityV = new double[NumNodes];
            Pressure = new double[NumNodes];
            StressXX = new double[NumElements];
            StressYY = new double[NumElements];
            StressXY = new double[NumElements];
            VonMises = new double[NumElements];
        }

        public double[,] ElementCoords(int e)
        {
            int nen = NodesPerElement(ElementType);
            double[,] c = new double[nen, 2];
            for (int k = 0; k < nen; k++) { c[k, 0] = X[Connectivity[e, k]]; c[k, 1] = Y[Connectivity[e, k]]; }
            return c;
        }

        public int PressureDofsPerElement()
        {
            return (ElementType == ElementType.Q4) ? 4 : 3;
        }

        public static int NodesPerElement(ElementType t)
        {
            switch (t)
            {
                case ElementType.CST3: return 3;
                case ElementType.LT6: return 6;
                case ElementType.Q4: return 4;
                default: return 3;
            }
        }
    }
}
