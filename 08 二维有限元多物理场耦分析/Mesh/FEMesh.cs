using System;
using System.Collections.Generic;
using MultiPhysicsFEM2D.Elements;

namespace MultiPhysicsFEM2D.Mesh
{
    public enum BCType
    {
        FixX, FixY, FixBoth,            // 结构 Dirichlet
        TractionX, TractionY,           // 结构 Neumann (集中力)
        VelocityX, VelocityY,           // 流体 Dirichlet
        Pressure,                       // 流体压力 Dirichlet
        Temperature,                    // 热 Dirichlet
        HeatFlux                        // 热 Neumann
    }

    public class BoundaryCondition
    {
        public BCType Type;
        public int NodeId;
        public double Value;
    }

    public class FEMesh
    {
        public List<FiniteElement> Elements = new List<FiniteElement>();
        public List<BoundaryCondition> BCs = new List<BoundaryCondition>();

        // 节点
        public int NumNodes;
        public double[] X;
        public double[] Y;

        // 结果（节点）
        public double[] DisplacementU;
        public double[] DisplacementV;
        public double[] VelocityU;
        public double[] VelocityV;
        public double[] Pressure;
        public double[] Temperature;
        public double[] StressXX;
        public double[] StressYY;
        public double[] StressXY;
        public double[] VonMises;
        public double[] PlasticStrain;

        // 材料默认参数（生成网格时拷贝到单元）
        public double YoungModulus = 2.1e11;
        public double PoissonRatio = 0.3;
        public double Density = 7800.0;
        public double Thickness = 1.0;
        public double YieldStress = 2.5e8;
        public double HardeningModulus = 2.1e9;
        public double Conductivity = 50.0;
        public double HeatCapacity = 500.0;
        public double ThermalExpansion = 1.2e-5;
        public double ReferenceTemperature = 293.0;
        public double Viscosity = 0.001;

        public void AllocateNodes(int n)
        {
            NumNodes = n;
            X = new double[n]; Y = new double[n];
        }

        public void AllocateResults()
        {
            DisplacementU = new double[NumNodes];
            DisplacementV = new double[NumNodes];
            VelocityU = new double[NumNodes];
            VelocityV = new double[NumNodes];
            Pressure = new double[NumNodes];
            Temperature = new double[NumNodes];
            StressXX = new double[NumNodes];
            StressYY = new double[NumNodes];
            StressXY = new double[NumNodes];
            VonMises = new double[NumNodes];
            PlasticStrain = new double[NumNodes];
        }

        /// <summary>
        /// 每个单元更新自身 X/Y 副本（拷贝自全局 X[NodeIds[i]]）
        /// </summary>
        public void UpdateElementCoords()
        {
            foreach (FiniteElement e in Elements) e.UpdateCoords(X, Y);
        }

        public void ApplyDefaultMaterial()
        {
            foreach (FiniteElement e in Elements)
            {
                e.YoungModulus = YoungModulus;
                e.PoissonRatio = PoissonRatio;
                e.Density = Density;
                e.Thickness = Thickness;
                e.YieldStress = YieldStress;
                e.HardeningModulus = HardeningModulus;
                e.Conductivity = Conductivity;
                e.HeatCapacity = HeatCapacity;
                e.ThermalExpansion = ThermalExpansion;
                e.ReferenceTemperature = ReferenceTemperature;
                e.Viscosity = Viscosity;
            }
        }

        /// <summary>
        /// 把每个单元的形心应力投影到节点（按面积加权平均）。
        /// </summary>
        public void ProjectStressFromElements()
        {
            double[] wsum = new double[NumNodes];
            Array.Clear(StressXX, 0, NumNodes);
            Array.Clear(StressYY, 0, NumNodes);
            Array.Clear(StressXY, 0, NumNodes);
            Array.Clear(VonMises, 0, NumNodes);
            Array.Clear(PlasticStrain, 0, NumNodes);
            foreach (FiniteElement e in Elements)
            {
                int nn = e.NodesPerElement;
                double a = e.ComputeArea();
                double w = a * e.Thickness / nn;
                for (int i = 0; i < nn; i++)
                {
                    int nid = e.NodeIds[i];
                    StressXX[nid] += w * e.StressXX;
                    StressYY[nid] += w * e.StressYY;
                    StressXY[nid] += w * e.StressXY;
                    VonMises[nid] += w * e.VonMises;
                    PlasticStrain[nid] += w * e.PlasticStrain;
                    wsum[nid] += w;
                }
            }
            for (int i = 0; i < NumNodes; i++)
            {
                if (wsum[i] > 1e-30)
                {
                    double inv = 1.0 / wsum[i];
                    StressXX[i] *= inv;
                    StressYY[i] *= inv;
                    StressXY[i] *= inv;
                    VonMises[i] *= inv;
                    PlasticStrain[i] *= inv;
                }
            }
        }
    }
}
