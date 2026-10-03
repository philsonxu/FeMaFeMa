using System;
using FEM2D.Mesh;

namespace FEM2D.Elements
{
    /// <summary>
    /// 单元公共接口。支持刚度矩阵、质量矩阵、对流矩阵、粘性矩阵、压力梯度耦合矩阵、面积、应力计算、速度插值。
    /// 坐标数组均按节点顺序为 [x0,y0, x1,y1, ...]。
    /// </summary>
    public abstract class FiniteElement
    {
        public abstract ElementType Type { get; }
        public abstract int NodesPerElement { get; }
        public abstract int SpatialDim { get; }

        /// <summary>
        /// 单元刚度矩阵（线弹性平面应变），按 2*NPE 顺序排列自由度 (u0,v0,u1,v1,...)
        /// </summary>
        public abstract void ComputeStiffness(FEMesh mesh, int elemId, double E, double nu, double thickness, double[,] ke);

        /// <summary>
        /// 一致质量矩阵（密度 rho）
        /// </summary>
        public abstract void ComputeMass(FEMesh mesh, int elemId, double rho, double thickness, double[,] me);

        /// <summary>
        /// 粘性（扩散）矩阵（标量粘度 mu），输出速度自由度大小 NPE*2 的方阵
        /// </summary>
        public abstract void ComputeViscous(FEMesh mesh, int elemId, double mu, double thickness, double[,] ke);

        /// <summary>
        /// 对流矩阵 C(u)：C_ij = ∫ N_i (u·∇) N_j dΩ ，速度由当前 ue,ve（逐节点）传入；输出 2*NPE 方阵
        /// </summary>
        public abstract void ComputeConvection(FEMesh mesh, int elemId, double[] ue, double[] ve, double thickness, double[,] ce);

        /// <summary>
        /// 压力梯度散度耦合：
        ///   B1[i,j] = -∫ ∂N_i/∂x M_j dΩ     (u-p)
        ///   B2[i,j] = -∫ ∂N_i/∂y M_j dΩ     (v-p)
        ///   D[j,i]  = -∫ M_j ∂N_i/∂x dΩ 等  (p-u), 对连续压力场 M_j=N_j 时 D = -B^T
        /// 这里等阶单元用 N_i 作为压力插值，因此输出 B 为 (NPE*2) x NPE 矩阵。
        /// </summary>
        public abstract void ComputePressureCoupling(FEMesh mesh, int elemId, double thickness, double[,] B);

        public abstract double ComputeArea(FEMesh mesh, int elemId);

        /// <summary>
        /// 给定单元位移 ue=[u0,v0,...]，返回单元中心处应力 σxx,σyy,σxy,vonMises
        /// </summary>
        public abstract double[] ComputeStress(FEMesh mesh, int elemId, double E, double nu, double[] ue);

        /// <summary>
        /// 节点坐标数组
        /// </summary>
        protected double[] GetCoords(FEMesh mesh, int elemId)
        {
            int[] conn = mesh.Elements[elemId];
            int npe = conn.Length;
            double[] xy = new double[npe * 2];
            for (int k = 0; k < npe; k++)
            {
                xy[2 * k] = mesh.X(conn[k]);
                xy[2 * k + 1] = mesh.Y(conn[k]);
            }
            return xy;
        }
    }
}
