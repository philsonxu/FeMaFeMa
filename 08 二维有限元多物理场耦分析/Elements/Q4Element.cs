using System;

namespace MultiPhysicsFEM2D.Elements
{
    /// <summary>
    /// 4节点双线性四边形等参单元 (Q4)，2×2 Gauss 积分。
    /// 局部坐标 ξ,η ∈ [-1,1]。
    /// 形状函数: N_i = (1+ξξ_i)(1+ηη_i)/4
    /// </summary>
    public class Q4Element : FiniteElement
    {
        private const int NN = 4;
        // 节点顺序: (0)-1,-1 (1)+1,-1 (2)+1,+1 (3)-1,+1
        private static readonly double[] _xi = new double[] { -1, 1, 1, -1 };
        private static readonly double[] _eta = new double[] { -1, -1, 1, 1 };
        // 2×2 Gauss
        private static readonly double[] _g = new double[] { -1.0 / Math.Sqrt(3.0), 1.0 / Math.Sqrt(3.0) };
        private static readonly double _gw = 1.0;

        public Q4Element()
        {
            Type = ElementType.Q4;
            AllocateNodes(NN);
        }

        public override int NodesPerElement { get { return NN; } }
        public override int PressureNodesPerElement { get { return NN; } }
        public override int[] GetPressureNodeMap()
        {
            int[] m = new int[NN];
            for (int i = 0; i < NN; i++) m[i] = i;
            return m;
        }

        private void ShapeAt(double xi, double eta, double[] N, double[] dNdxi, double[] dNdeta)
        {
            for (int i = 0; i < NN; i++)
            {
                double xi_i = _xi[i], et_i = _eta[i];
                N[i] = 0.25 * (1.0 + xi * xi_i) * (1.0 + eta * et_i);
                dNdxi[i] = 0.25 * xi_i * (1.0 + eta * et_i);
                dNdeta[i] = 0.25 * et_i * (1.0 + xi * xi_i);
            }
        }

        private double Jacobian(double xi, double eta, double[] dNdxi, double[] dNdeta,
            double[] dNdx, double[] dNdy)
        {
            double dxdxi = 0, dydxi = 0, dxdeta = 0, dydeta = 0;
            for (int i = 0; i < NN; i++)
            {
                dxdxi += dNdxi[i] * X[i];
                dydxi += dNdxi[i] * Y[i];
                dxdeta += dNdeta[i] * X[i];
                dydeta += dNdeta[i] * Y[i];
            }
            double J = dxdxi * dydeta - dxdeta * dydxi;
            double invJ = 1.0 / J;
            for (int i = 0; i < NN; i++)
            {
                dNdx[i] = (dNdxi[i] * dydeta - dNdeta[i] * dydxi) * invJ;
                dNdy[i] = (dNdeta[i] * dxdxi - dNdxi[i] * dxdeta) * invJ;
            }
            return J;
        }

        private delegate void GaussAction(double[] N, double[] dNdx, double[] dNdy, double wdet);

        private void ForEachGaussPoint(GaussAction act)
        {
            double[] N = new double[NN];
            double[] dNdxi = new double[NN], dNdeta = new double[NN];
            double[] dNdx = new double[NN], dNdy = new double[NN];
            for (int gi = 0; gi < 2; gi++)
                for (int gj = 0; gj < 2; gj++)
                {
                    double xi = _g[gi], eta = _g[gj];
                    ShapeAt(xi, eta, N, dNdxi, dNdeta);
                    double J = Jacobian(xi, eta, dNdxi, dNdeta, dNdx, dNdy);
                    double wdet = _gw * _gw * Math.Abs(J) * Thickness;
                    double[] Nc = (double[])N.Clone();
                    double[] dxc = (double[])dNdx.Clone();
                    double[] dyc = (double[])dNdy.Clone();
                    act(Nc, dxc, dyc, wdet);
                }
        }

        public override double ComputeArea()
        {
            double a = 0.0;
            ForEachGaussPoint(delegate (double[] N, double[] dx, double[] dy, double w)
            {
                a += w / Thickness;
            });
            return a;
        }

        private void DPlaneStress(out double d00, out double d01, out double d11, out double d22)
        {
            double e = YoungModulus, nu = PoissonRatio;
            double f = e / (1.0 - nu * nu);
            d00 = f; d01 = f * nu; d11 = f; d22 = f * (1.0 - nu) * 0.5;
        }

        public override void BuildStiffnessMatrix(double[,] K)
        {
            double d00, d01, d11, d22;
            DPlaneStress(out d00, out d01, out d11, out d22);
            double[,] Dm = new double[3, 3] {
                { d00, d01, 0 }, { d01, d11, 0 }, { 0, 0, d22 }
            };
            ForEachGaussPoint(delegate (double[] N, double[] dNdx, double[] dNdy, double w)
            {
                double[,] B = new double[3, 2 * NN];
                for (int a = 0; a < NN; a++)
                {
                    B[0, 2 * a] = dNdx[a];
                    B[1, 2 * a + 1] = dNdy[a];
                    B[2, 2 * a] = dNdy[a];
                    B[2, 2 * a + 1] = dNdx[a];
                }
                for (int i = 0; i < 2 * NN; i++)
                    for (int j = 0; j < 2 * NN; j++)
                    {
                        double s = 0.0;
                        for (int p = 0; p < 3; p++)
                            for (int q = 0; q < 3; q++)
                                s += B[p, i] * Dm[p, q] * B[q, j];
                        K[i, j] += w * s;
                    }
            });
        }

        public override void BuildMassMatrix(double[,] M, bool lumped)
        {
            if (lumped)
            {
                double[] rs = new double[NN];
                ForEachGaussPoint(delegate (double[] N, double[] dx, double[] dy, double w)
                {
                    for (int a = 0; a < NN; a++)
                        for (int b = 0; b < NN; b++) rs[a] += Density * w * N[a] * N[b];
                });
                for (int a = 0; a < NN; a++) { M[2 * a, 2 * a] = rs[a]; M[2 * a + 1, 2 * a + 1] = rs[a]; }
            }
            else
            {
                ForEachGaussPoint(delegate (double[] N, double[] dx, double[] dy, double w)
                {
                    for (int a = 0; a < NN; a++)
                        for (int b = 0; b < NN; b++)
                        {
                            double v = Density * w * N[a] * N[b];
                            M[2 * a, 2 * b] += v; M[2 * a + 1, 2 * b + 1] += v;
                        }
                });
            }
        }

        public override void ComputeStress(double[] u, double[] v,
            out double txx, out double tyy, out double txy, out double vonMises)
        {
            double d00, d01, d11, d22;
            DPlaneStress(out d00, out d01, out d11, out d22);
            double sxx = 0;
            double syy = 0;
            double sxy = 0;
            vonMises = 0;
            double tw = 0;
            ForEachGaussPoint(delegate (double[] N, double[] dNdx, double[] dNdy, double w)
            {
                double ex = 0, ey = 0, gxy = 0;
                for (int a = 0; a < NN; a++)
                {
                    ex += dNdx[a] * u[a]; ey += dNdy[a] * v[a];
                    gxy += dNdy[a] * u[a] + dNdx[a] * v[a];
                }
                sxx += w * (d00 * ex + d01 * ey);
                syy += w * (d01 * ex + d11 * ey);
                sxy += w * d22 * gxy;
                tw += w;
            });
            sxx /= tw; syy /= tw; sxy /= tw;
            txx = sxx;
            tyy = syy;
            txy = sxy;
            vonMises = Math.Sqrt(sxx * sxx - sxx * syy + syy * syy + 3.0 * sxy * sxy);
        }

        public override void BuildThermalLoadVector(double[] TNodal, double TRef, double[] f)
        {
            double d00, d01, d11, d22;
            DPlaneStress(out d00, out d01, out d11, out d22);
            double alpha = ThermalExpansion;
            ForEachGaussPoint(delegate (double[] N, double[] dNdx, double[] dNdy, double w)
            {
                double dT = 0;
                for (int a = 0; a < NN; a++) dT += N[a] * (TNodal[a] - TRef);
                double et = alpha * dT;
                double sx = (d00 + d01) * et;
                double sy = (d01 + d11) * et;
                for (int a = 0; a < NN; a++)
                {
                    f[2 * a] += w * dNdx[a] * sx;
                    f[2 * a + 1] += w * dNdy[a] * sy;
                }
            });
        }

        public override void BuildConductivityMatrix(double[,] Kt)
        {
            double k = Conductivity;
            ForEachGaussPoint(delegate (double[] N, double[] dNdx, double[] dNdy, double w)
            {
                double wk = k * w;
                for (int a = 0; a < NN; a++)
                    for (int b = 0; b < NN; b++)
                        Kt[a, b] += wk * (dNdx[a] * dNdx[b] + dNdy[a] * dNdy[b]);
            });
        }

        public override void BuildCapacityMatrix(double[,] C, bool lumped)
        {
            double m0 = Density * HeatCapacity;
            if (lumped)
            {
                double[] rs = new double[NN];
                ForEachGaussPoint(delegate (double[] N, double[] dx, double[] dy, double w)
                {
                    for (int a = 0; a < NN; a++)
                        for (int b = 0; b < NN; b++) rs[a] += m0 * w * N[a] * N[b];
                });
                for (int a = 0; a < NN; a++) C[a, a] = rs[a];
            }
            else
            {
                ForEachGaussPoint(delegate (double[] N, double[] dx, double[] dy, double w)
                {
                    for (int a = 0; a < NN; a++)
                        for (int b = 0; b < NN; b++)
                            C[a, b] += m0 * w * N[a] * N[b];
                });
            }
        }

        public override void BuildViscousMatrix(double[,] Visc)
        {
            double mu = Viscosity;
            ForEachGaussPoint(delegate (double[] N, double[] dNdx, double[] dNdy, double w)
            {
                double wm = mu * w;
                for (int a = 0; a < NN; a++)
                    for (int b = 0; b < NN; b++)
                    {
                        double lap = dNdx[a] * dNdx[b] + dNdy[a] * dNdy[b];
                        Visc[2 * a, 2 * b] += wm * lap;
                        Visc[2 * a + 1, 2 * b + 1] += wm * lap;
                    }
            });
        }

        public override void BuildConvectionMatrix(double[] uNodal, double[] vNodal, double[,] Conv)
        {
            double rho = Density;
            double[] us = new double[NN], vs = new double[NN];
            Array.Copy(uNodal, us, NN); Array.Copy(vNodal, vs, NN);
            ForEachGaussPoint(delegate (double[] N, double[] dNdx, double[] dNdy, double w)
            {
                double uc = 0, vc = 0;
                for (int a = 0; a < NN; a++) { uc += N[a] * us[a]; vc += N[a] * vs[a]; }
                double wr = rho * w;
                for (int a = 0; a < NN; a++)
                    for (int b = 0; b < NN; b++)
                    {
                        double val = wr * N[a] * (uc * dNdx[b] + vc * dNdy[b]);
                        Conv[2 * a, 2 * b] += val;
                        Conv[2 * a + 1, 2 * b + 1] += val;
                    }
            });
        }

        public override void BuildGradientMatrix(double[,] Gx, double[,] Gy)
        {
            ForEachGaussPoint(delegate (double[] N, double[] dNdx, double[] dNdy, double w)
            {
                for (int i = 0; i < NN; i++)
                {
                    double ni = N[i] * w;
                    for (int j = 0; j < NN; j++)
                    {
                        Gx[i, 2 * j] += ni * dNdx[j];
                        Gy[i, 2 * j + 1] += ni * dNdy[j];
                    }
                }
            });
        }

        public override void BuildPressureMassMatrix(double[,] Mp)
        {
            ForEachGaussPoint(delegate (double[] N, double[] dx, double[] dy, double w)
            {
                for (int a = 0; a < NN; a++)
                    for (int b = 0; b < NN; b++)
                        Mp[a, b] += w * N[a] * N[b];
            });
        }

        public override void BuildEnergyConvection(double[] uNodal, double[] vNodal, double[,] Ec)
        {
            double rhoCp = Density * HeatCapacity;
            double[] us = new double[NN], vs = new double[NN];
            Array.Copy(uNodal, us, NN); Array.Copy(vNodal, vs, NN);
            ForEachGaussPoint(delegate (double[] N, double[] dNdx, double[] dNdy, double w)
            {
                double uc = 0, vc = 0;
                for (int a = 0; a < NN; a++) { uc += N[a] * us[a]; vc += N[a] * vs[a]; }
                double wc = rhoCp * w;
                for (int a = 0; a < NN; a++)
                    for (int b = 0; b < NN; b++)
                        Ec[a, b] += wc * N[a] * (uc * dNdx[b] + vc * dNdy[b]);
            });
        }

        public override double[] ShapeAt(double xi, double eta)
        {
            double[] N = new double[NN];
            double[] dx = new double[NN], dy = new double[NN];
            ShapeAt(xi, eta, N, dx, dy);
            return N;
        }

        public override double[,] ShapeDerivLocal(double xi, double eta)
        {
            double[] N = new double[NN];
            double[] dxi = new double[NN], deta = new double[NN];
            ShapeAt(xi, eta, N, dxi, deta);
            double[,] dd = new double[NN, 2];
            for (int i = 0; i < NN; i++) { dd[i, 0] = dxi[i]; dd[i, 1] = deta[i]; }
            return dd;
        }
    }
}
