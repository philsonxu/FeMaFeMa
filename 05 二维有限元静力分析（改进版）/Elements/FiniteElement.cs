// FiniteElement.cs - 二维有限元单元抽象基类（含所有单元共享的实例字段与辅助方法）
using System;

namespace FEM2D.Elements
{
    public enum ElementType { CST3, LT6, Q4 }

    public abstract class FiniteElement
    {
        // 单元标识与拓扑
        public int Id;
        public int[] NodeIds;
        // 节点坐标（由 FEMesh.AssignNodeCoordinatesToElements 填充）
        public double[] X;
        public double[] Y;
        // 材料参数
        public double YoungModulus;
        public double PoissonRatio;
        public double Density;
        public double Viscosity;
        public double Thickness;
        // 单元级结果（单元中心/平均应力，由求解器填充）
        public double StressXX;
        public double StressYY;
        public double StressXY;
        public double VonMises;

        public abstract ElementType Type { get; }
        public abstract int NodesPerElement { get; }
        public abstract int DofPerNode { get; }
        public abstract int NumGaussPoints { get; }

        public abstract void ComputeStiffness(double[,] coords, double E, double nu, double thickness, double[,] ke);
        public abstract void ComputeMass(double[,] coords, double rho, double thickness, double[,] me);
        public abstract void ComputeConvection(double[,] coords, double[] uNodal, double thickness, double[,] ce);
        public abstract void ComputeViscous(double[,] coords, double mu, double thickness, double[,] ve);
        public abstract void ComputePressureGradient(double[,] coords, double thickness, double[,] gx, double[,] gy);
        public abstract double ComputeArea(double[,] coords);
        public abstract void ComputeStress(double[,] coords, double[] ue, double E, double nu,
            out double sxx, out double syy, out double sxy, out double von);
        public abstract void ShapeFunctions(double[,] coords, double xi, double eta,
            double[] N, double[] dNdx, double[] dNdy, out double detJ);

        // ========== 便捷辅助方法（基于抽象方法构造） ==========
        public double[,] GetCoordsMatrix()
        {
            int nen = NodesPerElement;
            double[,] c = new double[nen, 2];
            for (int i = 0; i < nen; i++) { c[i, 0] = X[i]; c[i, 1] = Y[i]; }
            return c;
        }

        public double[,] BuildStiffnessMatrix()
        {
            int nd = NodesPerElement * DofPerNode;
            double[,] ke = new double[nd, nd];
            double[,] c = GetCoordsMatrix();
            ComputeStiffness(c, YoungModulus, PoissonRatio, Thickness, ke);
            return ke;
        }

        public double[,] BuildMassMatrix()
        {
            int nd = NodesPerElement * DofPerNode;
            double[,] me = new double[nd, nd];
            double[,] c = GetCoordsMatrix();
            ComputeMass(c, Density, Thickness, me);
            return me;
        }

        public int PressureDofsCount
        {
            get
            {
                // CST3/LT6 压力 P1（角点 3 个），Q4 压力 P1（等阶 4 个）
                return (Type == ElementType.Q4) ? 4 : 3;
            }
        }
    }
}
