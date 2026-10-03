using System;
using FEM2D.Mesh;

namespace FEM2D.Elements
{
    /// <summary>
    /// 六节点二次三角形单元（LT6 / CPS6）
    /// 面积坐标 (L1,L2,L3), L1+L2+L3=1；自然坐标取 (L1,L2)∈三角形。
    /// 顶点节点 (1,0,0),(0,1,0),(0,0,1)，边中点 (0.5,0.5,0),(0,0.5,0.5),(0.5,0,0.5)。
    /// 采用 3 点 Hammer 积分（对二次多项式精确）。
    /// </summary>
    public sealed class LT6Element : FiniteElement
    {
        public override ElementType Type => ElementType.LT6;
        public override int NodesPerElement => 6;
        public override int SpatialDim => 2;

        // 3 点面积坐标高斯：权 1/3 位于 (2/3,1/6,1/6) 轮换
        private static readonly double[] Wt = new double[] { 1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0 };
        private static readonly double[,] Gp = new double[,]
        {
            { 2.0/3.0, 1.0/6.0, 1.0/6.0 },
            { 1.0/6.0, 2.0/3.0, 1.0/6.0 },
            { 1.0/6.0, 1.0/6.0, 2.0/3.0 }
        };

        // 节点顺序：顶点0,1,2，边中点3(0-1),4(1-2),5(2-0)
        private static void ShapeAt(double L1, double L2, double L3, double[] N, double[,] dNdL)
        {
            N[0] = L1 * (2.0 * L1 - 1.0);
            N[1] = L2 * (2.0 * L2 - 1.0);
            N[2] = L3 * (2.0 * L3 - 1.0);
            N[3] = 4.0 * L1 * L2;
            N[4] = 4.0 * L2 * L3;
            N[5] = 4.0 * L3 * L1;
            // dN/dL1, dN/dL2, dN/dL3 —— 实际用到 d/dL1, d/dL2
            dNdL[0, 0] = 4.0 * L1 - 1.0; dNdL[0, 1] = 0;              dNdL[0, 2] = 0;
            dNdL[1, 0] = 0;              dNdL[1, 1] = 4.0 * L2 - 1.0; dNdL[1, 2] = 0;
            dNdL[2, 0] = -1 * (4.0 * L3 - 1.0); dNdL[2, 1] = -1 * (4.0 * L3 - 1.0); dNdL[2, 2] = 4.0 * L3 - 1.0;
            dNdL[3, 0] = 4.0 * L2;       dNdL[3, 1] = 4.0 * L1;       dNdL[3, 2] = 0;
            dNdL[4, 0] = 0;              dNdL[4, 1] = 4.0 * L3;       dNdL[4, 2] = 4.0 * L2;
            dNdL[5, 0] = 4.0 * L3;       dNdL[5, 1] = 0;              dNdL[5, 2] = 4.0 * L1;
        }

        private delegate void Kernel(double[] N, double[,] dNdx, double detJ, double w, double thickness,
            double[] xy, double[] ue, double[] ve, double[,] outMat, int rows, int cols);

        private void Integrate(FEMesh mesh, int elemId, double thickness, double[] xy,
            double[] ue, double[] ve, double[,] outMat, int rows, int cols, Kernel fn)
        {
            double[] N = new double[6];
            double[,] dNdL = new double[6, 3];
            double[,] dNdx = new double[6, 2];

            double x0 = xy[0], y0 = xy[2], x1 = xy[2], y1 = xy[3], x2 = xy[4], y2 = xy[5];
            // 参考三角形到物理三角形的线性变换：
            //   x = x0 + (x1-x0) L1 + (x2-x0) L2？ 不对，LT6 用三个顶点坐标线性映射顶点，边中点是线性中点。
            //   因为二次单元含边中点，真实映射是等参：x = Σ N_i x_i 对所有 6 节点
            for (int ip = 0; ip < 3; ip++)
            {
                double L1 = Gp[ip, 0], L2 = Gp[ip, 1], L3 = Gp[ip, 2];
                ShapeAt(L1, L2, L3, N, dNdL);
                // Jacobian J = Σ dN_i/d(L1,L2) * (x_i,y_i)
                double dxL1 = 0, dxL2 = 0, dyL1 = 0, dyL2 = 0;
                for (int k = 0; k < 6; k++)
                {
                    double xk = xy[2 * k], yk = xy[2 * k + 1];
                    double d1 = dNdL[k, 0] - dNdL[k, 2];
                    double d2 = dNdL[k, 1] - dNdL[k, 2];
                    dxL1 += d1 * xk; dyL1 += d1 * yk;
                    dxL2 += d2 * xk; dyL2 += d2 * yk;
                }
                double detJ = dxL1 * dyL2 - dxL2 * dyL1;
                double inv00 = dyL2 / detJ;
                double inv01 = -dxL2 / detJ;
                double inv10 = -dyL1 / detJ;
                double inv11 = dxL1 / detJ;
                for (int k = 0; k < 6; k++)
                {
                    double d1 = dNdL[k, 0] - dNdL[k, 2];
                    double d2 = dNdL[k, 1] - dNdL[k, 2];
                    dNdx[k, 0] = inv00 * d1 + inv01 * d2;
                    dNdx[k, 1] = inv10 * d1 + inv11 * d2;
                }
                // 参考三角形积分面积为 0.5（在 L1+L2<=1, L1,L2>=0 上）
                double w = Wt[ip] * Math.Abs(detJ) * thickness * 0.5;
                fn(N, dNdx, Math.Abs(detJ), w, thickness, xy, ue, ve, outMat, rows, cols);
            }
        }

        private static double[,] DMatrixPS(double E, double nu)
        {
            double[,] D = new double[3, 3];
            double c = E / ((1.0 + nu) * (1.0 - 2.0 * nu));
            D[0, 0] = c * (1 - nu); D[0, 1] = c * nu;
            D[1, 0] = c * nu; D[1, 1] = c * (1 - nu);
            D[2, 2] = c * (1 - 2 * nu) * 0.5;
            return D;
        }

        public override void ComputeStiffness(FEMesh mesh, int elemId, double E, double nu, double thickness, double[,] ke)
        {
            double[] xy = GetCoords(mesh, elemId);
            double[,] D = DMatrixPS(E, nu);
            Array.Clear(ke, 0, ke.Length);
            int ndof = 12;
            Integrate(mesh, elemId, thickness, xy, null, null, ke, ndof, ndof,
                delegate (double[] N, double[,] dNdx, double detJ, double w, double t, double[] _xy, double[] _ue, double[] _ve, double[,] outM, int or_, int oc_)
                {
                    double[,] B = new double[3, 12];
                    for (int k = 0; k < 6; k++)
                    {
                        B[0, 2 * k] = dNdx[k, 0];
                        B[1, 2 * k + 1] = dNdx[k, 1];
                        B[2, 2 * k] = dNdx[k, 1];
                        B[2, 2 * k + 1] = dNdx[k, 0];
                    }
                    for (int i = 0; i < ndof; i++)
                        for (int j = 0; j < ndof; j++)
                        {
                            double s = 0;
                            for (int m = 0; m < 3; m++) for (int n = 0; n < 3; n++)
                                    s += B[n, i] * D[n, m] * B[m, j];
                            outM[i, j] += s * w;
                        }
                });
        }

        public override void ComputeMass(FEMesh mesh, int elemId, double rho, double thickness, double[,] me)
        {
            double[] xy = GetCoords(mesh, elemId);
            Array.Clear(me, 0, me.Length);
            Integrate(mesh, elemId, thickness, xy, null, null, me, 12, 12,
                delegate (double[] N, double[,] dNdx, double detJ, double w, double t, double[] _xy, double[] _ue, double[] _ve, double[,] outM, int or_, int oc_)
                {
                    double wi = w * rho;
                    for (int i = 0; i < 6; i++)
                        for (int j = 0; j < 6; j++)
                        {
                            double m = N[i] * N[j] * wi;
                            outM[2 * i, 2 * j] += m;
                            outM[2 * i + 1, 2 * j + 1] += m;
                        }
                });
        }

        public override void ComputeViscous(FEMesh mesh, int elemId, double mu, double thickness, double[,] ke)
        {
            double[] xy = GetCoords(mesh, elemId);
            Array.Clear(ke, 0, ke.Length);
            Integrate(mesh, elemId, thickness, xy, null, null, ke, 12, 12,
                delegate (double[] N, double[,] dNdx, double detJ, double w, double t, double[] _xy, double[] _ue, double[] _ve, double[,] outM, int or_, int oc_)
                {
                    double wi = w * mu;
                    for (int i = 0; i < 6; i++)
                        for (int j = 0; j < 6; j++)
                        {
                            double g = dNdx[i, 0] * dNdx[j, 0] + dNdx[i, 1] * dNdx[j, 1];
                            outM[2 * i, 2 * j] += g * wi;
                            outM[2 * i + 1, 2 * j + 1] += g * wi;
                        }
                });
        }

        public override void ComputeConvection(FEMesh mesh, int elemId, double[] ue, double[] ve, double thickness, double[,] ce)
        {
            double[] xy = GetCoords(mesh, elemId);
            Array.Clear(ce, 0, ce.Length);
            Integrate(mesh, elemId, thickness, xy, ue, ve, ce, 12, 12,
                delegate (double[] N, double[,] dNdx, double detJ, double w, double t, double[] _xy, double[] uE, double[] vE, double[,] outM, int or_, int oc_)
                {
                    double u = 0, v = 0;
                    for (int k = 0; k < 6; k++) { u += N[k] * uE[k]; v += N[k] * vE[k]; }
                    for (int i = 0; i < 6; i++)
                        for (int j = 0; j < 6; j++)
                        {
                            double adv = u * dNdx[j, 0] + v * dNdx[j, 1];
                            outM[2 * i, 2 * j] += N[i] * adv * w;
                            outM[2 * i + 1, 2 * j + 1] += N[i] * adv * w;
                        }
                });
        }

        public override void ComputePressureCoupling(FEMesh mesh, int elemId, double thickness, double[,] Bmat)
        {
            double[] xy = GetCoords(mesh, elemId);
            Array.Clear(Bmat, 0, Bmat.Length);
            Integrate(mesh, elemId, thickness, xy, null, null, Bmat, 12, 6,
                delegate (double[] N, double[,] dNdx, double detJ, double w, double t, double[] _xy, double[] _ue, double[] _ve, double[,] outM, int or_, int oc_)
                {
                    for (int i = 0; i < 6; i++)
                        for (int j = 0; j < 6; j++)
                        {
                            outM[2 * i, j] -= dNdx[i, 0] * N[j] * w;
                            outM[2 * i + 1, j] -= dNdx[i, 1] * N[j] * w;
                        }
                });
        }

        public override double ComputeArea(FEMesh mesh, int elemId)
        {
            double[] xy = GetCoords(mesh, elemId);
            double area = 0;
            Integrate(mesh, elemId, 1.0, xy, null, null, new double[1, 1], 1, 1,
                delegate (double[] N, double[,] dNdx, double detJ, double w, double t, double[] _xy, double[] _ue, double[] _ve, double[,] outM, int or_, int oc_)
                { area += w / t; });
            return area;
        }

        public override double[] ComputeStress(FEMesh mesh, int elemId, double E, double nu, double[] ue)
        {
            // 取单元形心（面积重心 L1=L2=L3=1/3）
            double[] xy = GetCoords(mesh, elemId);
            double[] N = new double[6];
            double[,] dNdL = new double[6, 3];
            ShapeAt(1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0, N, dNdL);
            double dxL1 = 0, dxL2 = 0, dyL1 = 0, dyL2 = 0;
            for (int k = 0; k < 6; k++)
            {
                double xk = xy[2 * k], yk = xy[2 * k + 1];
                double d1 = dNdL[k, 0] - dNdL[k, 2];
                double d2 = dNdL[k, 1] - dNdL[k, 2];
                dxL1 += d1 * xk; dyL1 += d1 * yk;
                dxL2 += d2 * xk; dyL2 += d2 * yk;
            }
            double detJ = dxL1 * dyL2 - dxL2 * dyL1;
            double inv00 = dyL2 / detJ, inv01 = -dxL2 / detJ;
            double inv10 = -dyL1 / detJ, inv11 = dxL1 / detJ;
            double[,] dNdx = new double[6, 2];
            for (int k = 0; k < 6; k++)
            {
                double d1 = dNdL[k, 0] - dNdL[k, 2];
                double d2 = dNdL[k, 1] - dNdL[k, 2];
                dNdx[k, 0] = inv00 * d1 + inv01 * d2;
                dNdx[k, 1] = inv10 * d1 + inv11 * d2;
            }
            double ex = 0, ey = 0, gxy = 0;
            for (int k = 0; k < 6; k++)
            {
                ex += dNdx[k, 0] * ue[2 * k];
                ey += dNdx[k, 1] * ue[2 * k + 1];
                gxy += dNdx[k, 1] * ue[2 * k] + dNdx[k, 0] * ue[2 * k + 1];
            }
            double[,] D = DMatrixPS(E, nu);
            double sx = D[0, 0] * ex + D[0, 1] * ey;
            double sy = D[1, 0] * ex + D[1, 1] * ey;
            double sxy = D[2, 2] * gxy;
            double vm = Math.Sqrt(sx * sx - sx * sy + sy * sy + 3.0 * sxy * sxy);
            return new double[] { sx, sy, sxy, vm };
        }
    }
}
