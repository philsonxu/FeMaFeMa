using System;
using NonlinearFEM2D.Materials;

namespace NonlinearFEM2D.Elements
{
    public enum ElementType { CST3, LT6, Q4 }

    // Gauss 点历史（用于非线性）
    public sealed class GPHistory
    {
        public double[] Stress = new double[4];      // [xx,yy,zz,xy]
        public double[] PlasticStrain = new double[4];
        public double EqPlasticStrain;
        public double VonMises;
        public bool Yielded;
        public double J; // 上次 Jacobian 行列式（后处理用）
        public double Wt;
        public void CloneFrom(GPHistory o)
        {
            Array.Copy(o.Stress, this.Stress, 4);
            Array.Copy(o.PlasticStrain, this.PlasticStrain, 4);
            this.EqPlasticStrain = o.EqPlasticStrain;
            this.VonMises = o.VonMises;
            this.Yielded = o.Yielded;
            this.J = o.J;
            this.Wt = o.Wt;
        }
    }

    public abstract class FiniteElement
    {
        public int Id;
        public int[] NodeIds;
        public double[] X, Y;
        public double YoungModulus = 2.1e11;
        public double PoissonRatio = 0.3;
        public double Density = 7850.0;
        public double Viscosity = 1.0e-3;
        public double Thickness = 0.01;
        public J2PlasticityMaterial Material;
        public ElementType Type;
        public int NodesPerElement { get; protected set; }
        public int NumGaussPoints { get; protected set; }
        public GPHistory[] Histories;

        // 位移场：u[i], v[i] 为 NodeIds[i] 节点的 u,v；返回一致切线 Ke(npe*2, npe*2) 和内力向量 fint(2*npe)
        public abstract void BuildTangentAndInternal(double[] u, double[] v, double[,] Ke, double[] fint, bool geometricNonlinearity);
        // 线弹性刚度矩阵（用于线性静力/模态/L 矩阵）
        public abstract void BuildLinearStiffness(double[,] Ke);
        // 一致质量矩阵
        public abstract void BuildMassMatrix(double[,] Me);
        // 对流矩阵（NS）
        public abstract void BuildConvection(double[] uc, double[] vc, double[,] C);
        // 粘性矩阵（NS）
        public abstract void BuildViscousMatrix(double[,] Kv);
        // 压力梯度矩阵（NS），npePressure=3/6/4
        public abstract void BuildPressureGradient(double[,] G);
        // 单元面积
        public abstract double ComputeArea();
        // 按给定位移计算中心处应力（线性后处理）
        public abstract void ComputeStress(double[] u, double[] v, out double sxx, out double syy, out double sxy);
        // 内部初始化（如 Gauss 点历史数组）
        public virtual void Initialize()
        {
            this.Material = new J2PlasticityMaterial(this.YoungModulus, this.PoissonRatio, 2.5e8, this.YoungModulus / 100.0);
            this.Histories = new GPHistory[this.NumGaussPoints];
            for (int g = 0; g < this.NumGaussPoints; g++)
            {
                this.Histories[g] = new GPHistory();
            }
        }
        // 克隆当前历史到目标（载荷步/迭代回退用）
        public GPHistory[] CloneHistories()
        {
            GPHistory[] c = new GPHistory[this.NumGaussPoints];
            for (int g = 0; g < this.NumGaussPoints; g++)
            {
                c[g] = new GPHistory();
                c[g].CloneFrom(this.Histories[g]);
            }
            return c;
        }
        public void RestoreHistories(GPHistory[] snap)
        {
            for (int g = 0; g < this.NumGaussPoints; g++) this.Histories[g].CloneFrom(snap[g]);
        }
    }
}
