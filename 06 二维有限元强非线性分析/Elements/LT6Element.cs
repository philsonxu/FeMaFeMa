using System;

namespace NonlinearFEM2D.Elements
{
    public sealed class LT6Element : FiniteElement
    {
        // 6节点三角形（角点 0,1,2 + 边中点 3,4,5），3 点 Hammer 积分
        public LT6Element()
        {
            this.NodesPerElement = 6;
            this.NumGaussPoints = 3;
            this.Type = ElementType.LT6;
        }

        // 面积坐标上的形函数在点(L1,L2,L3)
        private static void Shape(double L1, double L2, double L3, double[] N)
        {
            N[0] = L1 * (2.0 * L1 - 1.0);
            N[1] = L2 * (2.0 * L2 - 1.0);
            N[2] = L3 * (2.0 * L3 - 1.0);
            N[3] = 4.0 * L1 * L2;
            N[4] = 4.0 * L2 * L3;
            N[5] = 4.0 * L3 * L1;
        }
        private static void ShapeDeriv(double L1, double L2, double L3, double[,] dNdL)
        {
            dNdL[0, 0] = 4.0 * L1 - 1.0; dNdL[0, 1] = 0.0;          dNdL[0, 2] = 0.0;
            dNdL[1, 0] = 0.0;          dNdL[1, 1] = 4.0 * L2 - 1.0; dNdL[1, 2] = 0.0;
            dNdL[2, 0] = 0.0;          dNdL[2, 1] = 0.0;          dNdL[2, 2] = 4.0 * L3 - 1.0;
            dNdL[3, 0] = 4.0 * L2;     dNdL[3, 1] = 4.0 * L1;     dNdL[3, 2] = 0.0;
            dNdL[4, 0] = 0.0;          dNdL[4, 1] = 4.0 * L3;     dNdL[4, 2] = 4.0 * L2;
            dNdL[5, 0] = 4.0 * L3;     dNdL[5, 1] = 0.0;          dNdL[5, 2] = 4.0 * L1;
        }

        private bool Jacobian(double L1, double L2, double L3, out double J00, out double J01, out double J10, out double J11, out double detJ, out double[,] dNdx, out double[,] dNdy)
        {
            double[,] dNdL = new double[6, 3]; ShapeDeriv(L1, L2, L3, dNdL);
            J00 = 0; J01 = 0; J10 = 0; J11 = 0;
            for (int i = 0; i < 6; i++)
            {
                double dL1 = dNdL[i, 0], dL2 = dNdL[i, 1], dL3 = dNdL[i, 2];
                double dx = dL1 * this.X[0] + dL2 * this.X[1] + dL3 * this.X[2];
                double dy = dL1 * this.Y[0] + dL2 * this.Y[1] + dL3 * this.Y[2];
                J00 += dx; J01 += dy; J10 += dx; J11 += dy;
            }
            detJ = J00 * J11 - J01 * J10;
            dNdx = new double[6, 1]; dNdy = new double[6, 1];
            if (Math.Abs(detJ) < 1.0e-20) return false;
            double invDet = 1.0 / detJ;
            dNdx = new double[6, 1]; dNdy = new double[6, 1];
            for (int i = 0; i < 6; i++)
            {
                double dL1 = dNdL[i, 0], dL2 = dNdL[i, 1], dL3 = dNdL[i, 2];
                double dx = dL1 * this.X[0] + dL2 * this.X[1] + dL3 * this.X[2];
                double dy = dL1 * this.Y[0] + dL2 * this.Y[1] + dL3 * this.Y[2];
                // dN/dx = dN/dL * dL/dx ; 由 J = [[xL1,yL1],[xL2,yL2]] (L3=1-L1-L2)
                double xL1 = this.X[0] - this.X[2], yL1 = this.Y[0] - this.Y[2];
                double xL2 = this.X[1] - this.X[2], yL2 = this.Y[1] - this.Y[2];
                double detT = xL1 * yL2 - xL2 * yL1;
                double dL1dx = yL2 / detT, dL1dy = -xL2 / detT;
                double dL2dx = -yL1 / detT, dL2dy = xL1 / detT;
                double dL3dx = -dL1dx - dL2dx, dL3dy = -dL1dy - dL2dy;
                double dNxI = dL1 * dL1dx + dL2 * dL2dx + dL3 * dL3dx;
                double dNyI = dL1 * dL1dy + dL2 * dL2dy + dL3 * dL3dy;
                dNdx[i, 0] = dNxI;
                dNdy[i, 0] = dNyI;
            }
            return true;
        }

        public override double ComputeArea()
        {
            double x0 = this.X[0], y0 = this.Y[0], x1 = this.X[1], y1 = this.Y[1], x2 = this.X[2], y2 = this.Y[2];
            return 0.5 * Math.Abs(x0 * (y1 - y2) + x1 * (y2 - y0) + x2 * (y0 - y1));
        }

        private double[,] Dmatrix()
        {
            double E = this.YoungModulus, nu = this.PoissonRatio;
            double c = E / (1.0 - nu * nu);
            double[,] D = new double[3, 3];
            D[0, 0] = c; D[0, 1] = c * nu;
            D[1, 0] = c * nu; D[1, 1] = c;
            D[2, 2] = c * (1.0 - nu) / 2.0;
            return D;
        }

        // 3 点 Hammer 积分
        private static readonly double[] GP_L1 = { 1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0 };
        private static readonly double[] GP_L2 = { 1.0 / 6.0, 1.0 / 6.0, 2.0 / 3.0 };
        private const double GP_W = 1.0 / 3.0;

        public override void BuildLinearStiffness(double[,] Ke)
        {
            Array.Clear(Ke, 0, Ke.Length);
            double[,] D = Dmatrix();
            double t = this.Thickness, area = this.ComputeArea();
            for (int g = 0; g < 3; g++)
            {
                double L1 = GP_L1[g], L2 = GP_L2[g], L3 = 1.0 - L1 - L2;
                double J00, J01, J10, J11, detJ;
                double[,] dNdx, dNdy;
                if (!Jacobian(L1, L2, L3, out J00, out J01, out J10, out J11, out detJ, out dNdx, out dNdy)) continue;
                double w = GP_W * Math.Abs(detJ) * t;
                double[,] B = new double[3, 12];
                for (int i = 0; i < 6; i++)
                {
                    B[0, 2 * i] = dNdx[i, 0];
                    B[1, 2 * i + 1] = dNdy[i, 0];
                    B[2, 2 * i] = dNdy[i, 0]; B[2, 2 * i + 1] = dNdx[i, 0];
                }
                for (int i = 0; i < 12; i++) for (int j = 0; j < 12; j++)
                    {
                        double s = 0.0;
                        for (int m = 0; m < 3; m++) for (int n = 0; n < 3; n++) s += B[m, i] * D[m, n] * B[n, j];
                        Ke[i, j] += s * w;
                    }
            }
        }

        // LT6 非线性：对非线性求解退化为线性切线 + 内力由线性应力计算（几何非线性可独立开启）
        public override void BuildTangentAndInternal(double[] u, double[] v, double[,] Ke, double[] fint, bool geomNonlin)
        {
            BuildLinearStiffness(Ke);
            Array.Clear(fint, 0, fint.Length);
            // 内力 ∫ B^T σ dΩ 由线性弹性本构从位移求得
            double[,] D = Dmatrix();
            double t = this.Thickness;
            for (int g = 0; g < 3; g++)
            {
                double L1 = GP_L1[g], L2 = GP_L2[g], L3 = 1.0 - L1 - L2;
                double J00, J01, J10, J11, detJ;
                double[,] dNdx, dNdy;
                if (!Jacobian(L1, L2, L3, out J00, out J01, out J10, out J11, out detJ, out dNdx, out dNdy)) continue;
                double w = GP_W * Math.Abs(detJ) * t;
                double exx = 0, eyy = 0, gxy = 0;
                for (int i = 0; i < 6; i++) { exx += dNdx[i, 0] * u[i]; eyy += dNdy[i, 0] * v[i]; gxy += dNdy[i, 0] * u[i] + dNdx[i, 0] * v[i]; }
                double sxx = D[0, 0] * exx + D[0, 1] * eyy;
                double syy = D[1, 0] * exx + D[1, 1] * eyy;
                double sxy = D[2, 2] * gxy;
                this.Histories[g].Stress[0] = sxx; this.Histories[g].Stress[1] = syy; this.Histories[g].Stress[3] = sxy;
                this.Histories[g].VonMises = Math.Sqrt(Math.Max(0.0, sxx * sxx - sxx * syy + syy * syy + 3.0 * sxy * sxy));
                for (int i = 0; i < 12; i++)
                {
                    double f = 0.0;
                    f += (i % 2 == 0 ? dNdx[i / 2, 0] : 0.0) * sxx;
                    f += (i % 2 == 1 ? dNdy[i / 2, 0] : 0.0) * syy;
                    f += (i % 2 == 0 ? dNdy[i / 2, 0] : dNdx[i / 2, 0]) * sxy;
                    fint[i] += f * w;
                }
            }
        }

        public override void BuildMassMatrix(double[,] Me)
        {
            Array.Clear(Me, 0, Me.Length);
            double t = this.Thickness, rho = this.Density;
            // Row-sum lumped + 行一致混合（简化：一致质量通过 Hammer 积分）
            double[] N = new double[6];
            for (int g = 0; g < 3; g++)
            {
                double L1 = GP_L1[g], L2 = GP_L2[g], L3 = 1.0 - L1 - L2;
                double J00, J01, J10, J11, detJ;
                double[,] dNdx, dNdy;
                if (!Jacobian(L1, L2, L3, out J00, out J01, out J10, out J11, out detJ, out dNdx, out dNdy)) continue;
                Shape(L1, L2, L3, N);
                double w = GP_W * Math.Abs(detJ) * t * rho;
                for (int a = 0; a < 6; a++) for (int b = 0; b < 6; b++)
                    { Me[2 * a, 2 * b] += N[a] * N[b] * w; Me[2 * a + 1, 2 * b + 1] += N[a] * N[b] * w; }
            }
        }

        public override void BuildConvection(double[] uc, double[] vc, double[,] C)
        {
            Array.Clear(C, 0, C.Length);
            double t = this.Thickness;
            double[] N = new double[6];
            for (int g = 0; g < 3; g++)
            {
                double L1 = GP_L1[g], L2 = GP_L2[g], L3 = 1.0 - L1 - L2;
                double J00, J01, J10, J11, detJ;
                double[,] dNdx, dNdy;
                if (!Jacobian(L1, L2, L3, out J00, out J01, out J10, out J11, out detJ, out dNdx, out dNdy)) continue;
                Shape(L1, L2, L3, N);
                double uavg = 0, vavg = 0;
                for (int i = 0; i < 6; i++) { uavg += N[i] * uc[i]; vavg += N[i] * vc[i]; }
                double w = GP_W * Math.Abs(detJ) * t;
                for (int a = 0; a < 6; a++) for (int b = 0; b < 6; b++)
                    {
                        double cv = N[a] * (uavg * dNdx[b, 0] + vavg * dNdy[b, 0]) * w;
                        C[2 * a, 2 * b] += cv; C[2 * a + 1, 2 * b + 1] += cv;
                    }
            }
        }

        public override void BuildViscousMatrix(double[,] Kv)
        {
            Array.Clear(Kv, 0, Kv.Length);
            double mu = this.Viscosity, t = this.Thickness;
            for (int g = 0; g < 3; g++)
            {
                double L1 = GP_L1[g], L2 = GP_L2[g], L3 = 1.0 - L1 - L2;
                double J00, J01, J10, J11, detJ;
                double[,] dNdx, dNdy;
                if (!Jacobian(L1, L2, L3, out J00, out J01, out J10, out J11, out detJ, out dNdx, out dNdy)) continue;
                double w = GP_W * Math.Abs(detJ) * t * mu;
                for (int a = 0; a < 6; a++) for (int b = 0; b < 6; b++)
                    {
                        double v = (dNdx[a, 0] * dNdx[b, 0] + dNdy[a, 0] * dNdy[b, 0]) * w;
                        Kv[2 * a, 2 * b] += v; Kv[2 * a + 1, 2 * b + 1] += v;
                    }
            }
        }

        public override void BuildPressureGradient(double[,] G)
        {
            Array.Clear(G, 0, G.Length);
            double t = this.Thickness;
            double[] N = new double[6];
            for (int g = 0; g < 3; g++)
            {
                double L1 = GP_L1[g], L2 = GP_L2[g], L3 = 1.0 - L1 - L2;
                double J00, J01, J10, J11, detJ;
                double[,] dNdx, dNdy;
                if (!Jacobian(L1, L2, L3, out J00, out J01, out J10, out J11, out detJ, out dNdx, out dNdy)) continue;
                Shape(L1, L2, L3, N);
                double w = GP_W * Math.Abs(detJ) * t;
                for (int a = 0; a < 6; a++) for (int bp = 0; bp < 3; bp++) // pressure P2-P1? 这里简单使用角点3个压力自由度
                    {
                        double Np = (bp == 0) ? L1 : (bp == 1 ? L2 : L3);
                        G[2 * a, bp] -= N[a] * dNdx[bp, 0] * w; // 简化
                        G[2 * a + 1, bp] -= N[a] * dNdy[bp, 0] * w;
                    }
            }
        }

        public override void ComputeStress(double[] u, double[] v, out double sxx, out double syy, out double sxy)
        {
            sxx = 0; syy = 0; sxy = 0;
            double total = 0;
            for (int g = 0; g < 3; g++)
            {
                sxx += this.Histories[g].Stress[0]; syy += this.Histories[g].Stress[1]; sxy += this.Histories[g].Stress[3]; total += 1.0;
            }
            sxx /= total; syy /= total; sxy /= total;
        }
    }
}
