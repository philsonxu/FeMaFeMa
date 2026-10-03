using System;

namespace MultiPhysicsFEM2D.Elements
{
    public enum ElementType
    {
        CST3,   // 3节点常应变三角
        LT6,    // 6节点二次三角
        Q4      // 4节点双线性四边形
    }

    /// <summary>
    /// 二维有限单元基类。所有派生单元必须严格实现下面的抽象方法。
    /// </summary>
    public abstract class FiniteElement
    {
        public int Id;
        public ElementType Type;
        public int[] NodeIds;
        public double[] X;
        public double[] Y;

        // 结构材料
        public double YoungModulus = 2.1e11;
        public double PoissonRatio = 0.3;
        public double Density = 7800.0;
        public double Thickness = 1.0;
        public double YieldStress = 2.5e8;
        public double HardeningModulus = 2.1e9;

        // 热材料
        public double Conductivity = 50.0;      // k [W/(m·K)]
        public double HeatCapacity = 500.0;     // c
        public double ThermalExpansion = 1.2e-5; // α [1/K]
        public double ReferenceTemperature = 293.0;

        // 流体
        public double Viscosity = 0.001;        // μ [Pa·s]

        // 单元应力（形心/平均）
        public double StressXX;
        public double StressYY;
        public double StressXY;
        public double VonMises;
        public double PlasticStrain;

        public void AllocateNodes(int n)
        {
            NodeIds = new int[n];
            X = new double[n];
            Y = new double[n];
        }

        public void UpdateCoords(double[] gX, double[] gY)
        {
            int n = NodeIds.Length;
            for (int i = 0; i < n; i++) { X[i] = gX[NodeIds[i]]; Y[i] = gY[NodeIds[i]]; }
        }

        // —— 基本信息 ——
        public abstract int NodesPerElement { get; }
        public abstract int PressureNodesPerElement { get; }
        public abstract double ComputeArea();
        public abstract int[] GetPressureNodeMap();

        // —— 结构矩阵 ——
        public abstract void BuildStiffnessMatrix(double[,] K);
        public abstract void BuildMassMatrix(double[,] M, bool lumped);
        public abstract void ComputeStress(double[] u, double[] v,
            out double sxx, out double syy, out double sxy, out double vonMises);
        public abstract void BuildThermalLoadVector(double[] TNodal, double TRef, double[] f);

        // —— 热传导矩阵 ——
        public abstract void BuildConductivityMatrix(double[,] Kt);
        public abstract void BuildCapacityMatrix(double[,] C, bool lumped);

        // —— NS 矩阵 ——
        public abstract void BuildViscousMatrix(double[,] Visc);
        public abstract void BuildConvectionMatrix(double[] uNodal, double[] vNodal, double[,] Conv);
        public abstract void BuildGradientMatrix(double[,] Gx, double[,] Gy);
        public abstract void BuildPressureMassMatrix(double[,] Mp);
        public abstract void BuildEnergyConvection(double[] uNodal, double[] vNodal, double[,] Ec);

        // 形函数工具
        public abstract double[] ShapeAt(double xi, double eta);
        public abstract double[,] ShapeDerivLocal(double xi, double eta);  // dN/dxi, dN/deta
    }
}
