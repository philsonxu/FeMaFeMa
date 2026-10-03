using System;
using FEM2D.Mesh;

namespace FEM2D.Elements
{
    /// <summary>
    /// 三节点常应变三角单元（CST / CPS3）
    /// 形函数：N_i = (a_i + b_i x + c_i y) / (2A)
    /// 应变矩阵 B 为常数矩阵；使用 1 点高斯精确积分。
    /// </summary>
    public sealed class CST3Element : FiniteElement
    {
        public override ElementType Type => ElementType.CST3;
        public override int NodesPerElement => 3;
        public override int SpatialDim => 2;

        private struct CSTGeom
        {
            public double Area;
            public double b0, b1, b2;
            public double c0, c1, c2;
            public double x0, y0, x1, y1, x2, y2;
        }

        private CSTGeom ComputeGeom(FEMesh mesh, int elemId)
        {
            int[] c = mesh.Elements[elemId];
            double x0 = mesh.X(c[0]), y0 = mesh.Y(c[0]);
            double x1 = mesh.X(c[1]), y1 = mesh.Y(c[1]);
            double x2 = mesh.X(c[2]), y2 = mesh.Y(c[2]);
            double b0 = y1 - y2, b1 = y2 - y0, b2 = y0 - y1;
            double c0 = x2 - x1, c1 = x0 - x2, c2 = x1 - x0;
            double A2 = b0 * c1 - b1 * c0; // = 2A * det(sign)
            double A = 0.5 * Math.Abs(b0 * (x1 - x2) + b1 * (x2 - x0) + b2 * (x0 - x1));
            // 正确的 2A（含符号）：
            double twoA = (x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0);
            if (Math.Abs(twoA) < 1e-30) twoA = 1e-30;
            return new CSTGeom { Area = 0.5 * Math.Abs(twoA), b0 = b0, b1 = b1, b2 = b2, c0 = c0, c1 = c1, c2 = c2, x0 = x0, y0 = y0, x1 = x1, y1 = y1, x2 = x2, y2 = y2 };
        }

        private static double[,] DMatrix(double E, double nu, bool planeStress)
        {
            double[,] D = new double[3, 3];
            if (planeStress)
            {
                double c = E / (1.0 - nu * nu);
                D[0, 0] = c; D[0, 1] = c * nu; D[0, 2] = 0;
                D[1, 0] = c * nu; D[1, 1] = c; D[1, 2] = 0;
                D[2, 0] = 0; D[2, 1] = 0; D[2, 2] = c * (1 - nu) * 0.5;
            }
            else
            {
                double c = E / ((1.0 + nu) * (1.0 - 2.0 * nu));
                D[0, 0] = c * (1 - nu); D[0, 1] = c * nu; D[0, 2] = 0;
                D[1, 0] = c * nu; D[1, 1] = c * (1 - nu); D[1, 2] = 0;
                D[2, 0] = 0; D[2, 1] = 0; D[2, 2] = c * (1 - 2 * nu) * 0.5;
            }
            return D;
        }

        public override void ComputeStiffness(FEMesh mesh, int elemId, double E, double nu, double thickness, double[,] ke)
        {
            CSTGeom g = ComputeGeom(mesh, elemId);
            double A = g.Area;
            double inv2A = 1.0 / (2.0 * A);
            double b0 = g.b0 * inv2A, b1 = g.b1 * inv2A, b2 = g.b2 * inv2A;
            double c0 = g.c0 * inv2A, c1 = g.c1 * inv2A, c2 = g.c2 * inv2A;

            double[,] B = new double[3, 6]
            {
                { b0, 0,  b1, 0,  b2, 0 },
                { 0,  c0, 0,  c1, 0,  c2 },
                { c0, b0, c1, b1, c2, b2 }
            };
            double[,] D = DMatrix(E, nu, false);
            double[,] BTDB = new double[6, 6];
            for (int i = 0; i < 6; i++)
                for (int j = 0; j < 6; j++)
                {
                    double s = 0;
                    for (int m = 0; m < 3; m++) for (int k = 0; k < 3; k++)
                            s += B[k, i] * D[k, m] * B[m, j];
                    BTDB[i, j] = s;
                }
            double wt = A * thickness;
            for (int i = 0; i < 6; i++)
                for (int j = 0; j < 6; j++)
                    ke[i, j] = BTDB[i, j] * wt;
        }

        public override void ComputeMass(FEMesh mesh, int elemId, double rho, double thickness, double[,] me)
        {
            double A = ComputeArea(mesh, elemId);
            double c = rho * thickness * A / 12.0;
            int ndof = 6;
            for (int i = 0; i < ndof; i++)
                for (int j = 0; j < ndof; j++)
                    me[i, j] = c;
            for (int k = 0; k < 3; k++)
            {
                me[2 * k, 2 * k] = 2 * c;
                me[2 * k + 1, 2 * k + 1] = 2 * c;
            }
        }

        public override void ComputeViscous(FEMesh mesh, int elemId, double mu, double thickness, double[,] ke)
        {
            // Stokes 粘性: mu ∫ ∇u : ∇v   （包含耦合项：简化为 2*mu*(∇_s u : ∇_s v)）
            // 这里使用矢量拉普拉斯形式 mu ∫ ∇u·∇v，为 NS 标准 Galerkin 离散。
            CSTGeom g = ComputeGeom(mesh, elemId);
            double A = g.Area;
            double inv2A = 1.0 / (2.0 * A);
            double b0 = g.b0 * inv2A, b1 = g.b1 * inv2A, b2 = g.b2 * inv2A;
            double c0 = g.c0 * inv2A, c1 = g.c1 * inv2A, c2 = g.c2 * inv2A;

            double[] dNdx = new double[] { b0, b1, b2 };
            double[] dNdy = new double[] { c0, c1, c2 };

            double coef = mu * thickness * A;
            Array.Clear(ke, 0, ke.Length);
            for (int i = 0; i < 3; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    double val = coef * (dNdx[i] * dNdx[j] + dNdy[i] * dNdy[j]);
                    ke[2 * i, 2 * j] = val;
                    ke[2 * i + 1, 2 * j + 1] = val;
                }
            }
        }

        public override void ComputeConvection(FEMesh mesh, int elemId, double[] ue, double[] ve, double thickness, double[,] ce)
        {
            // u·∇N_j 在常应变单元上为常数；形函数导数常数，N 在形心处=1/3
            CSTGeom g = ComputeGeom(mesh, elemId);
            double A = g.Area;
            double inv2A = 1.0 / (2.0 * A);
            double b0 = g.b0 * inv2A, b1 = g.b1 * inv2A, b2 = g.b2 * inv2A;
            double c0 = g.c0 * inv2A, c1 = g.c1 * inv2A, c2 = g.c2 * inv2A;

            double uc = (ue[0] + ue[1] + ue[2]) / 3.0;
            double vc = (ve[0] + ve[1] + ve[2]) / 3.0;

            double[] dNdx = new double[] { b0, b1, b2 };
            double[] dNdy = new double[] { c0, c1, c2 };

            // C_ij = N_i(centroid) * (u dN_j/dx + v dN_j/dy) * A * t;
            Array.Clear(ce, 0, ce.Length);
            double wt = thickness * A;
            for (int i = 0; i < 3; i++)
            {
                double Ni = 1.0 / 3.0;
                for (int j = 0; j < 3; j++)
                {
                    double adv = uc * dNdx[j] + vc * dNdy[j];
                    ce[2 * i, 2 * j] = Ni * adv * wt;
                    ce[2 * i + 1, 2 * j + 1] = Ni * adv * wt;
                }
            }
        }

        public override void ComputePressureCoupling(FEMesh mesh, int elemId, double thickness, double[,] B)
        {
            CSTGeom g = ComputeGeom(mesh, elemId);
            double A = g.Area;
            double inv2A = 1.0 / (2.0 * A);
            double b0 = g.b0 * inv2A, b1 = g.b1 * inv2A, b2 = g.b2 * inv2A;
            double c0 = g.c0 * inv2A, c1 = g.c1 * inv2A, c2 = g.c2 * inv2A;
            double[] dNdx = new double[] { b0, b1, b2 };
            double[] dNdy = new double[] { c0, c1, c2 };
            double wt = thickness * A;
            // 1 点积分：-∫ ∂N_i/∂x * N_j dΩ = -∂N_i/∂x * N_j(c) * wt
            Array.Clear(B, 0, B.Length);
            for (int i = 0; i < 3; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    B[2 * i, j] = -dNdx[i] * (1.0 / 3.0) * wt;
                    B[2 * i + 1, j] = -dNdy[i] * (1.0 / 3.0) * wt;
                }
            }
        }

        public override double ComputeArea(FEMesh mesh, int elemId)
        {
            CSTGeom g = ComputeGeom(mesh, elemId);
            return g.Area;
        }

        public override double[] ComputeStress(FEMesh mesh, int elemId, double E, double nu, double[] ue)
        {
            CSTGeom g = ComputeGeom(mesh, elemId);
            double A = g.Area;
            double inv2A = 1.0 / (2.0 * A);
            double b0 = g.b0 * inv2A, b1 = g.b1 * inv2A, b2 = g.b2 * inv2A;
            double c0 = g.c0 * inv2A, c1 = g.c1 * inv2A, c2 = g.c2 * inv2A;

            double ex = b0 * ue[0] + b1 * ue[2] + b2 * ue[4];
            double ey = c0 * ue[1] + c1 * ue[3] + c2 * ue[5];
            double gxy = c0 * ue[0] + b0 * ue[1] + c1 * ue[2] + b1 * ue[3] + c2 * ue[4] + b2 * ue[5];

            double[,] D = DMatrix(E, nu, false);
            double sx = D[0, 0] * ex + D[0, 1] * ey;
            double sy = D[1, 0] * ex + D[1, 1] * ey;
            double sxy = D[2, 2] * gxy;
            double vm = Math.Sqrt(sx * sx - sx * sy + sy * sy + 3.0 * sxy * sxy);
            return new double[] { sx, sy, sxy, vm };
        }
    }
}
