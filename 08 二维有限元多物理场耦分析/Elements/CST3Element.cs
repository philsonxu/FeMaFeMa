using System;

namespace MultiPhysicsFEM2D.Elements
{
    /// <summary>
    /// 3 节点常应变三角形单元 (CST, CPS3)。
    /// 面积坐标形函数，1 点 Gauss 积分（形心）。
    /// u/v: 3 DOF × 2 = 6 DOF; pressure: 线性 P1（3 节点共节点，带压力质量稳定化）。
    /// </summary>
    public class CST3Element : FiniteElement
    {
        public CST3Element()
        {
            Type = ElementType.CST3;
            AllocateNodes(3);
        }

        public override int NodesPerElement { get { return 3; } }
        public override int PressureNodesPerElement { get { return 3; } }

        public override int[] GetPressureNodeMap()
        {
            int[] m = new int[3];
            for (int i = 0; i < 3; i++) m[i] = i;
            return m;
        }

        // 常应变：B 矩阵与坐标无关
        private void BMatrix(out double bx0, out double by0, out double gx1,
                              out double bx1, out double by1, out double gx2,
                              out double bx2, out double by2, out double gx3,
                              out double area)
        {
            double x0 = X[0], y0 = Y[0];
            double x1 = X[1], y1 = Y[1];
            double x2 = X[2], y2 = Y[2];
            double b0 = y1 - y2;
            double b1 = y2 - y0;
            double b2 = y0 - y1;
            double c0 = x2 - x1;
            double c1 = x0 - x2;
            double c2 = x1 - x0;
            area = 0.5 * Math.Abs(b0 * c1 - b1 * c0);
            double inv2A = 1.0 / (2.0 * area);
            bx0 = b0 * inv2A; by0 = c0 * inv2A;
            bx1 = b1 * inv2A; by1 = c1 * inv2A;
            bx2 = b2 * inv2A; by2 = c2 * inv2A;
            gx1 = inv2A; gx2 = inv2A; gx3 = inv2A;
            // gx1,gx2,gx3 placeholder; γ_xy uses B3 row
        }

        public override double ComputeArea()
        {
            double x0 = X[0], y0 = Y[0];
            double x1 = X[1], y1 = Y[1];
            double x2 = X[2], y2 = Y[2];
            return 0.5 * Math.Abs((x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0));
        }

        // 平面应力 D 矩阵（结构用平面应力，NS 用单元粘性矩阵）
        private void DPlaneStress(out double d00, out double d01, out double d11, out double d22)
        {
            double e = YoungModulus, nu = PoissonRatio;
            double factor = e / (1.0 - nu * nu);
            d00 = factor;
            d01 = factor * nu;
            d11 = factor;
            d22 = factor * (1.0 - nu) * 0.5;
        }

        public override void BuildStiffnessMatrix(double[,] K)
        {
            double area;
            double bx0, by0, bx1, by1, bx2, by2, g;
            BMatrix(out bx0, out by0, out g, out bx1, out by1, out g, out bx2, out by2, out g, out area);
            double d00, d01, d11, d22;
            DPlaneStress(out d00, out d01, out d11, out d22);
            double t = Thickness;
            // B 矩阵 3×6
            // 行0(εx): [bx0  0   bx1  0   bx2  0]
            // 行1(εy): [0   by0  0   by1  0   by2]
            // 行2(γxy):[by0 bx0 by1 bx1 by2 bx2]
            double[,] B = new double[3, 6];
            B[0, 0] = bx0; B[0, 2] = bx1; B[0, 4] = bx2;
            B[1, 1] = by0; B[1, 3] = by1; B[1, 5] = by2;
            B[2, 0] = by0; B[2, 1] = bx0; B[2, 2] = by1; B[2, 3] = bx1; B[2, 4] = by2; B[2, 5] = bx2;
            double[,] D = new double[3, 3];
            D[0, 0] = d00; D[0, 1] = d01;
            D[1, 0] = d01; D[1, 1] = d11;
            D[2, 2] = d22;
            // K = t*A * B^T D B
            double coef = t * area;
            // DB
            double[,] DB = new double[3, 6];
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 6; j++)
                {
                    double s = 0.0;
                    for (int k = 0; k < 3; k++) s += D[i, k] * B[k, j];
                    DB[i, j] = s;
                }
            for (int i = 0; i < 6; i++)
                for (int j = 0; j < 6; j++)
                {
                    double s = 0.0;
                    for (int k = 0; k < 3; k++) s += B[k, i] * DB[k, j];
                    K[i, j] = coef * s;
                }
        }

        public override void BuildMassMatrix(double[,] M, bool lumped)
        {
            double area = ComputeArea();
            double rho = Density, t = Thickness;
            double m = rho * t * area;
            if (lumped)
            {
                double ml = m / 3.0;
                for (int i = 0; i < 6; i++) M[i, i] = ml;
            }
            else
            {
                // 一致质量：∫ N_i N_j = A/12 (1+δij)
                for (int i = 0; i < 3; i++)
                    for (int j = 0; j < 3; j++)
                    {
                        double v = m * ((i == j) ? 1.0 / 6.0 : 1.0 / 12.0);
                        M[2 * i, 2 * j] += v;
                        M[2 * i + 1, 2 * j + 1] += v;
                    }
            }
        }

        public override void ComputeStress(double[] u, double[] v,
            out double sxx, out double syy, out double sxy, out double vonMises)
        {
            double area;
            double bx0, by0, bx1, by1, bx2, by2, g;
            BMatrix(out bx0, out by0, out g, out bx1, out by1, out g, out bx2, out by2, out g, out area);
            double ex = bx0 * u[0] + bx1 * u[1] + bx2 * u[2];
            double ey = by0 * v[0] + by1 * v[1] + by2 * v[2];
            double gxy = by0 * u[0] + bx0 * v[0] + by1 * u[1] + bx1 * v[1] + by2 * u[2] + bx2 * v[2];
            double d00, d01, d11, d22;
            DPlaneStress(out d00, out d01, out d11, out d22);
            sxx = d00 * ex + d01 * ey;
            syy = d01 * ex + d11 * ey;
            sxy = d22 * gxy;
            vonMises = Math.Sqrt(sxx * sxx - sxx * syy + syy * syy + 3.0 * sxy * sxy);
        }

        public override void BuildThermalLoadVector(double[] TNodal, double TRef, double[] f)
        {
            double area = ComputeArea();
            double dT0 = TNodal[0] - TRef, dT1 = TNodal[1] - TRef, dT2 = TNodal[2] - TRef;
            // 形心 ΔT（CST 取平均）
            double dT = (dT0 + dT1 + dT2) / 3.0;
            double e = YoungModulus, nu = PoissonRatio, alpha = ThermalExpansion, t = Thickness;
            double s = e * alpha * dT * t * area / (1.0 - nu);
            // 热应变 {αΔT, αΔT, 0} -> F_th = B^T D ε_th
            // 形函数梯度系数
            double x0 = X[0], y0 = Y[0], x1 = X[1], y1 = Y[1], x2 = X[2], y2 = Y[2];
            double b0 = y1 - y2, b1 = y2 - y0, b2 = y0 - y1;
            double c0 = x2 - x1, c1 = x0 - x2, c2 = x1 - x0;
            // B^T (D * ε_th) = B^T * [EαΔT/(1-ν), EαΔT/(1-ν), 0]^T
            // B 第一列对应 bx0=... 其实为 b_i/(2A)，于是 Fx_i = s * b_i / (2A) * 2A ???
            // 精确：B^T * [s0,s1,0] 中 B 行0 系数 = b_i/(2A); 乘以 s = 2A*EαΔT*t/(1-ν) → Fx_i = 0.5*t*b_i*EαΔT/(1-ν)
            // 用 area = 0.5 |b0c1-b1c0|，取符号一致版本：Fx_i = 0.5*(b_i)*EαΔT*t/(1-ν)
            double coef = 0.5 * e * alpha * dT * t / (1.0 - nu);
            f[0] = coef * b0; f[1] = coef * c0;
            f[2] = coef * b1; f[3] = coef * c1;
            f[4] = coef * b2; f[5] = coef * c2;
        }

        public override void BuildConductivityMatrix(double[,] Kt)
        {
            double area;
            double bx0, by0, bx1, by1, bx2, by2, g;
            BMatrix(out bx0, out by0, out g, out bx1, out by1, out g, out bx2, out by2, out g, out area);
            double k = Conductivity, t = Thickness;
            // Kt = k*t*A * (Bth^T Bth), Bth 行 = [b_i/(2A), c_i/(2A)] 即 [∂N/∂x, ∂N/∂y]
            double[,] G = new double[2, 3];
            G[0, 0] = bx0; G[0, 1] = bx1; G[0, 2] = bx2;
            G[1, 0] = by0; G[1, 1] = by1; G[1, 2] = by2;
            double coef = k * t * area;
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                    Kt[i, j] = coef * (G[0, i] * G[0, j] + G[1, i] * G[1, j]);
        }

        public override void BuildCapacityMatrix(double[,] C, bool lumped)
        {
            double area = ComputeArea();
            double rho = Density, cp = HeatCapacity, t = Thickness;
            double m = rho * cp * t * area;
            if (lumped)
            {
                double ml = m / 3.0;
                for (int i = 0; i < 3; i++) C[i, i] = ml;
            }
            else
            {
                for (int i = 0; i < 3; i++)
                    for (int j = 0; j < 3; j++)
                        C[i, j] = m * ((i == j) ? 1.0 / 6.0 : 1.0 / 12.0);
            }
        }

        public override void BuildViscousMatrix(double[,] Visc)
        {
            // μ*(∇u:∇v) —— 对 u 和 v 两个分量分别组装，2×3 DOF
            double area;
            double bx0, by0, bx1, by1, bx2, by2, g;
            BMatrix(out bx0, out by0, out g, out bx1, out by1, out g, out bx2, out by2, out g, out area);
            double mu = Viscosity * Thickness;
            double[,] G = new double[2, 3];
            G[0, 0] = bx0; G[0, 1] = bx1; G[0, 2] = bx2;
            G[1, 0] = by0; G[1, 1] = by1; G[1, 2] = by2;
            // 单分量 Laplacian 部分（含 ∂/∂x, ∂/∂y）：Kuu = μ*A*(Gx^T Gx + Gy^T Gy)
            // 再加 ∂u_i/∂x ∂v_i/∂x...  Stokes 粘性 = μ*(∇u + ∇u^T):∇v /2 *2 = μ∇u:∇v
            double coef = mu * area;
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                {
                    double s = coef * (G[0, i] * G[0, j] + G[1, i] * G[1, j]);
                    Visc[2 * i, 2 * j] += s;
                    Visc[2 * i + 1, 2 * j + 1] += s;
                }
        }

        public override void BuildConvectionMatrix(double[] uNodal, double[] vNodal, double[,] Conv)
        {
            // u·∇u : ∫ N_i (u·∇) N_j dΩ
            double area = ComputeArea();
            double rho = Density, t = Thickness;
            // u, v 取单元平均
            double uc = (uNodal[0] + uNodal[1] + uNodal[2]) / 3.0;
            double vc = (vNodal[0] + vNodal[1] + vNodal[2]) / 3.0;
            double x0 = X[0], y0 = Y[0], x1 = X[1], y1 = Y[1], x2 = X[2], y2 = Y[2];
            double b0 = y1 - y2, b1 = y2 - y0, b2 = y0 - y1;
            double c0 = x2 - x1, c1 = x0 - x2, c2 = x1 - x0;
            double inv2A = 1.0 / (2.0 * area);
            double[] dNdx = new double[] { b0 * inv2A, b1 * inv2A, b2 * inv2A };
            double[] dNdy = new double[] { c0 * inv2A, c1 * inv2A, c2 * inv2A };
            // ∫ N_i dN_j/dx = A/3 * dN_j/dx   (N 积分和)
            double coef = rho * t * area / 3.0;
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                {
                    double v = coef * (uc * dNdx[j] + vc * dNdy[j]);
                    Conv[2 * i, 2 * j] += v;
                    Conv[2 * i + 1, 2 * j + 1] += v;
                }
        }

        public override void BuildGradientMatrix(double[,] Gx, double[,] Gy)
        {
            // (∇·v, p): Gx[i,j] = ∫ N^p_i ∂N^u_j/∂x dΩ
            double area = ComputeArea();
            double t = Thickness;
            double x0 = X[0], y0 = Y[0], x1 = X[1], y1 = Y[1], x2 = X[2], y2 = Y[2];
            double b0 = y1 - y2, b1 = y2 - y0, b2 = y0 - y1;
            double c0 = x2 - x1, c1 = x0 - x2, c2 = x1 - x0;
            double inv2A = 1.0 / (2.0 * area);
            double ddx0 = b0 * inv2A, ddx1 = b1 * inv2A, ddx2 = b2 * inv2A;
            double ddy0 = c0 * inv2A, ddy1 = c1 * inv2A, ddy2 = c2 * inv2A;
            double a = t * area / 3.0; // ∫ N_i dΩ = A/3
            // u 分量: j=0,2,4 (节点 0,1,2 u)
            int[] uMap = new int[] { 0, 2, 4 };
            int[] vMap = new int[] { 1, 3, 5 };
            for (int i = 0; i < 3; i++)
            {
                Gx[i, uMap[0]] = a * ddx0; Gx[i, uMap[1]] = a * ddx1; Gx[i, uMap[2]] = a * ddx2;
                Gy[i, vMap[0]] = a * ddy0; Gy[i, vMap[1]] = a * ddy1; Gy[i, vMap[2]] = a * ddy2;
            }
        }

        public override void BuildPressureMassMatrix(double[,] Mp)
        {
            double area = ComputeArea();
            double t = Thickness;
            double m = t * area;
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                    Mp[i, j] = m * ((i == j) ? 1.0 / 6.0 : 1.0 / 12.0);
        }

        public override void BuildEnergyConvection(double[] uNodal, double[] vNodal, double[,] Ec)
        {
            // ρ c_p ∫ N_i u·∇N_j T dΩ （单分量温度）
            double area = ComputeArea();
            double rho = Density, cp = HeatCapacity, t = Thickness;
            double uc = (uNodal[0] + uNodal[1] + uNodal[2]) / 3.0;
            double vc = (vNodal[0] + vNodal[1] + vNodal[2]) / 3.0;
            double x0 = X[0], y0 = Y[0], x1 = X[1], y1 = Y[1], x2 = X[2], y2 = Y[2];
            double b0 = y1 - y2, b1 = y2 - y0, b2 = y0 - y1;
            double c0 = x2 - x1, c1 = x0 - x2, c2 = x1 - x0;
            double inv2A = 1.0 / (2.0 * area);
            double[] dNdx = new double[] { b0 * inv2A, b1 * inv2A, b2 * inv2A };
            double[] dNdy = new double[] { c0 * inv2A, c1 * inv2A, c2 * inv2A };
            double coef = rho * cp * t * area / 3.0;
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                    Ec[i, j] = coef * (uc * dNdx[j] + vc * dNdy[j]);
        }

        public override double[] ShapeAt(double xi, double eta)
        {
            double zeta = 1.0 - xi - eta;
            return new double[] { zeta, xi, eta };
        }

        public override double[,] ShapeDerivLocal(double xi, double eta)
        {
            // dN/dL1, dN/dL2 在面积坐标下; 这里返回 dN/dxi, dN/deta
            double[,] d = new double[3, 2];
            d[0, 0] = -1; d[0, 1] = -1;
            d[1, 0] = 1;  d[1, 1] = 0;
            d[2, 0] = 0;  d[2, 1] = 1;
            return d;
        }
    }
}
