using System;
using FEM2D.Mesh;

namespace FEM2D.Elements
{
    /// <summary>
    /// 4 节点双线性等参四边形单元（Q4）
    /// 自然坐标 ξ,η ∈ [-1,1]^2；形函数 N_i = 0.25(1+ξξ_i)(1+ηη_i)
    /// 采用 2x2 高斯积分（4 点）。
    /// </summary>
    public sealed class Q4Element : FiniteElement
    {
        public override ElementType Type => ElementType.Q4;
        public override int NodesPerElement => 4;
        public override int SpatialDim => 2;

        private static readonly double[] gp = new double[] { -1.0 / Math.Sqrt(3.0), 1.0 / Math.Sqrt(3.0) };
        private static readonly double gw = 1.0;

        // 等参形函数 N_i、导数（对 ξ/η）
        private static void ShapeAt(double xi, double eta, double[] N, double[,] dNdxi)
        {
            double[] s = new double[] { -1, 1, 1, -1 };
            double[] t = new double[] { -1, -1, 1, 1 };
            for (int i = 0; i < 4; i++)
            {
                N[i] = 0.25 * (1 + s[i] * xi) * (1 + t[i] * eta);
                dNdxi[i, 0] = 0.25 * s[i] * (1 + t[i] * eta);
                dNdxi[i, 1] = 0.25 * t[i] * (1 + s[i] * xi);
            }
        }

        private delegate void Kernel(double xi, double eta, double[] N, double[,] dNdx, double detJ, double weight, double thickness,
            double[] xy, double[] ue, double[] ve, double[,] outMat, int outRows, int outCols);

        private void Integrate(FEMesh mesh, int elemId, double thickness, double[] xy,
            double[] ue, double[] ve, double[,] outMat, int outRows, int outCols, Kernel fn)
        {
            double[] N = new double[4];
            double[,] dNdxi = new double[4, 2];
            double[,] dNdx = new double[4, 2];
            double[,] J = new double[2, 2];
            for (int a = 0; a < 2; a++)
            {
                for (int b = 0; b < 2; b++)
                {
                    double xi = gp[a], eta = gp[b];
                    ShapeAt(xi, eta, N, dNdxi);
                    J[0, 0] = 0; J[0, 1] = 0; J[1, 0] = 0; J[1, 1] = 0;
                    for (int k = 0; k < 4; k++)
                    {
                        double xk = xy[2 * k], yk = xy[2 * k + 1];
                        J[0, 0] += dNdxi[k, 0] * xk;
                        J[0, 1] += dNdxi[k, 0] * yk;
                        J[1, 0] += dNdxi[k, 1] * xk;
                        J[1, 1] += dNdxi[k, 1] * yk;
                    }
                    double detJ = J[0, 0] * J[1, 1] - J[0, 1] * J[1, 0];
                    if (detJ < 0) { /* 反向编号取绝对值 */ }
                    double inv00 = J[1, 1] / detJ;
                    double inv01 = -J[0, 1] / detJ;
                    double inv10 = -J[1, 0] / detJ;
                    double inv11 = J[0, 0] / detJ;
                    for (int k = 0; k < 4; k++)
                    {
                        dNdx[k, 0] = inv00 * dNdxi[k, 0] + inv01 * dNdxi[k, 1];
                        dNdx[k, 1] = inv10 * dNdxi[k, 0] + inv11 * dNdxi[k, 1];
                    }
                    fn(xi, eta, N, dNdx, Math.Abs(detJ), gw * gw, thickness, xy, ue, ve, outMat, outRows, outCols);
                }
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
            Integrate(mesh, elemId, thickness, xy, null, null, ke, 8, 8,
                delegate (double xi, double eta, double[] N, double[,] dNdx, double detJ, double wt, double t, double[] _xy, double[] _ue, double[] _ve, double[,] outM, int or_, int oc_)
                {
                    double w = wt * detJ * t;
                    // B (3x8)
                    double[,] B = new double[3, 8];
                    for (int k = 0; k < 4; k++)
                    {
                        B[0, 2 * k] = dNdx[k, 0];
                        B[1, 2 * k + 1] = dNdx[k, 1];
                        B[2, 2 * k] = dNdx[k, 1];
                        B[2, 2 * k + 1] = dNdx[k, 0];
                    }
                    for (int i = 0; i < 8; i++)
                        for (int j = 0; j < 8; j++)
                        {
                            double s = 0;
                            for (int m = 0; m < 3; m++)
                                for (int n = 0; n < 3; n++)
                                    s += B[n, i] * D[n, m] * B[m, j];
                            outM[i, j] += s * w;
                        }
                });
        }

        public override void ComputeMass(FEMesh mesh, int elemId, double rho, double thickness, double[,] me)
        {
            double[] xy = GetCoords(mesh, elemId);
            Array.Clear(me, 0, me.Length);
            Integrate(mesh, elemId, thickness, xy, null, null, me, 8, 8,
                delegate (double xi, double eta, double[] N, double[,] dNdx, double detJ, double wt, double t, double[] _xy, double[] _ue, double[] _ve, double[,] outM, int or_, int oc_)
                {
                    double w = wt * detJ * t * rho;
                    for (int i = 0; i < 4; i++)
                        for (int j = 0; j < 4; j++)
                        {
                            outM[2 * i, 2 * j] += N[i] * N[j] * w;
                            outM[2 * i + 1, 2 * j + 1] += N[i] * N[j] * w;
                        }
                });
        }

        public override void ComputeViscous(FEMesh mesh, int elemId, double mu, double thickness, double[,] ke)
        {
            double[] xy = GetCoords(mesh, elemId);
            Array.Clear(ke, 0, ke.Length);
            Integrate(mesh, elemId, thickness, xy, null, null, ke, 8, 8,
                delegate (double xi, double eta, double[] N, double[,] dNdx, double detJ, double wt, double t, double[] _xy, double[] _ue, double[] _ve, double[,] outM, int or_, int oc_)
                {
                    double w = wt * detJ * t * mu;
                    for (int i = 0; i < 4; i++)
                        for (int j = 0; j < 4; j++)
                        {
                            double grad = dNdx[i, 0] * dNdx[j, 0] + dNdx[i, 1] * dNdx[j, 1];
                            outM[2 * i, 2 * j] += grad * w;
                            outM[2 * i + 1, 2 * j + 1] += grad * w;
                        }
                });
        }

        public override void ComputeConvection(FEMesh mesh, int elemId, double[] ue, double[] ve, double thickness, double[,] ce)
        {
            double[] xy = GetCoords(mesh, elemId);
            Array.Clear(ce, 0, ce.Length);
            Integrate(mesh, elemId, thickness, xy, ue, ve, ce, 8, 8,
                delegate (double xi, double eta, double[] N, double[,] dNdx, double detJ, double wt, double t, double[] _xy, double[] uE, double[] vE, double[,] outM, int or_, int oc_)
                {
                    double u = 0, v = 0;
                    for (int k = 0; k < 4; k++) { u += N[k] * uE[k]; v += N[k] * vE[k]; }
                    double w = wt * detJ * t;
                    for (int i = 0; i < 4; i++)
                        for (int j = 0; j < 4; j++)
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
            Integrate(mesh, elemId, thickness, xy, null, null, Bmat, 8, 4,
                delegate (double xi, double eta, double[] N, double[,] dNdx, double detJ, double wt, double t, double[] _xy, double[] _ue, double[] _ve, double[,] outM, int or_, int oc_)
                {
                    double w = wt * detJ * t;
                    for (int i = 0; i < 4; i++)
                        for (int j = 0; j < 4; j++)
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
                delegate (double xi, double eta, double[] N, double[,] dNdx, double detJ, double wt, double t, double[] _xy, double[] _ue, double[] _ve, double[,] outM, int or_, int oc_)
                { area += wt * Math.Abs(detJ) * t; });
            return area;
        }

        public override double[] ComputeStress(FEMesh mesh, int elemId, double E, double nu, double[] ue)
        {
            // 在单元中心 (0,0) 处计算
            double[] xy = GetCoords(mesh, elemId);
            double[] N = new double[4];
            double[,] dNdxi = new double[4, 2];
            ShapeAt(0, 0, N, dNdxi);
            double[,] J = new double[2, 2];
            for (int k = 0; k < 4; k++)
            {
                double xk = xy[2 * k], yk = xy[2 * k + 1];
                J[0, 0] += dNdxi[k, 0] * xk; J[0, 1] += dNdxi[k, 0] * yk;
                J[1, 0] += dNdxi[k, 1] * xk; J[1, 1] += dNdxi[k, 1] * yk;
            }
            double detJ = J[0, 0] * J[1, 1] - J[0, 1] * J[1, 0];
            double inv00 = J[1, 1] / detJ, inv01 = -J[0, 1] / detJ;
            double inv10 = -J[1, 0] / detJ, inv11 = J[0, 0] / detJ;
            double[,] dNdx = new double[4, 2];
            for (int k = 0; k < 4; k++)
            {
                dNdx[k, 0] = inv00 * dNdxi[k, 0] + inv01 * dNdxi[k, 1];
                dNdx[k, 1] = inv10 * dNdxi[k, 0] + inv11 * dNdxi[k, 1];
            }
            double ex = 0, ey = 0, gxy = 0;
            for (int k = 0; k < 4; k++)
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
