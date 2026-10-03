using System;

namespace MultiPhysicsFEM2D.Elements
{
    /// <summary>
    /// 6节点二次三角形单元 (LT6/P2)。
    /// 3点Hammer面积坐标Gauss积分。u/v/P同阶P2（12个位移DOF，6个压力DOF）。
    /// 面积坐标形函数：
    ///   N1 = λ1(2λ1-1), N2 = λ2(2λ2-1), N3 = λ3(2λ3-1)
    ///   N4 = 4λ1λ2, N5 = 4λ2λ3, N6 = 4λ3λ1   (边中节点 1-2, 2-3, 3-1)
    /// </summary>
    public class LT6Element : FiniteElement
    {
        private const int NN = 6;
        // Hammer 3点积分
        private static readonly double[] _gp1 = new double[] { 0.5, 0.5, 0.0, 1.0 / 3.0 };
        private static readonly double[] _gp2 = new double[] { 0.0, 0.5, 0.5, 1.0 / 3.0 };
        private static readonly double[] _gp3 = new double[] { 0.5, 0.0, 0.5, 1.0 / 3.0 };
        private static readonly double[][] _gps = new double[][] { _gp1, _gp2, _gp3 };

        public LT6Element()
        {
            Type = ElementType.LT6;
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

        private static void ShapeAndDeriv(double L1, double L2, double L3,
            double[] N, double[] dNdL1, double[] dNdL2)
        {
            N[0] = L1 * (2.0 * L1 - 1.0);
            N[1] = L2 * (2.0 * L2 - 1.0);
            N[2] = L3 * (2.0 * L3 - 1.0);
            N[3] = 4.0 * L1 * L2;
            N[4] = 4.0 * L2 * L3;
            N[5] = 4.0 * L3 * L1;
            dNdL1[0] = 4.0 * L1 - 1.0; dNdL2[0] = 0.0;
            dNdL1[1] = 0.0; dNdL2[1] = 4.0 * L2 - 1.0;
            dNdL1[2] = -(4.0 * L3 - 1.0); dNdL2[2] = -(4.0 * L3 - 1.0);
            dNdL1[3] = 4.0 * L2; dNdL2[3] = 4.0 * L1;
            dNdL1[4] = -4.0 * L2; dNdL2[4] = 4.0 * (L3 - L2);
            dNdL1[5] = 4.0 * (L3 - L1); dNdL2[5] = -4.0 * L1;
        }

        private double Jacobian(double L1, double L2, double L3,
            double[] dNdL1, double[] dNdL2, double[] dNdx, double[] dNdy)
        {
            double dxdl1 = 0, dydl1 = 0, dxdl2 = 0, dydl2 = 0;
            for (int i = 0; i < NN; i++)
            {
                dxdl1 += dNdL1[i] * X[i];
                dydl1 += dNdL1[i] * Y[i];
                dxdl2 += dNdL2[i] * X[i];
                dydl2 += dNdL2[i] * Y[i];
            }
            double J = dxdl1 * dydl2 - dxdl2 * dydl1;
            double invJ = 1.0 / J;
            for (int i = 0; i < NN; i++)
            {
                dNdx[i] = (dNdL1[i] * dydl2 - dNdL2[i] * dydl1) * invJ;
                dNdy[i] = (dNdL2[i] * dxdl1 - dNdL1[i] * dxdl2) * invJ;
            }
            return J;
        }

        public override double ComputeArea()
        {
            // 参考3节点三角形（角点 0,1,2）
            double x0 = X[0], y0 = Y[0], x1 = X[1], y1 = Y[1], x2 = X[2], y2 = Y[2];
            return 0.5 * Math.Abs((x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0));
        }

        private void DPlaneStress(out double d00, out double d01, out double d11, out double d22)
        {
            double e = YoungModulus, nu = PoissonRatio;
            double f = e / (1.0 - nu * nu);
            d00 = f; d01 = f * nu; d11 = f; d22 = f * (1.0 - nu) * 0.5;
        }

        private delegate void GaussAction(double[] N, double[] dNdx, double[] dNdy, double wdet);

        private void ForEachGaussPoint(GaussAction act)
        {
            double[] N = new double[NN];
            double[] dNdL1 = new double[NN], dNdL2 = new double[NN];
            double[] dNdx = new double[NN], dNdy = new double[NN];
            for (int g = 0; g < 3; g++)
            {
                double L1 = _gps[g][0], L2 = _gps[g][1], L3 = _gps[g][2], w = _gps[g][3];
                ShapeAndDeriv(L1, L2, L3, N, dNdL1, dNdL2);
                double J = Jacobian(L1, L2, L3, dNdL1, dNdL2, dNdx, dNdy);
                double wdet = w * Math.Abs(J) * Thickness;
                double[] Nc = (double[])N.Clone();
                double[] dxc = (double[])dNdx.Clone();
                double[] dyc = (double[])dNdy.Clone();
                act(Nc, dxc, dyc, wdet);
            }
        }

        public override void BuildStiffnessMatrix(double[,] K)
        {
            double d00, d01, d11, d22;
            DPlaneStress(out d00, out d01, out d11, out d22);
            ForEachGaussPoint(delegate (double[] N, double[] dNdx, double[] dNdy, double w)
            {
                // B (3×12)
                for (int a = 0; a < NN; a++)
                {
                    int ua = 2 * a, va = 2 * a + 1;
                    double bx = dNdx[a], by = dNdy[a];
                    // εx = bx*u_a, εy = by*v_a, γxy = by*u_a + bx*v_a
                    // D·B:
                    double s0 = d00 * bx; double s1 = d01 * by;
                    double t0 = d01 * bx; double t1 = d11 * by;
                    double g0 = d22 * by; double g1 = d22 * bx;
                    for (int b = 0; b < NN; b++)
                    {
                        int ub = 2 * b, vb = 2 * b + 1;
                        double cx = dNdx[b], cy = dNdy[b];
                        K[ua, ub] += w * (s0 * cx + g0 * cy);
                        K[ua, vb] += w * (s1 * cy + g1 * cx); // ?
                    }
                }
                // 直接按 B^T D B 全量计算更清晰：
                // 先构造 B(3×2N), 再计算 K += w*B^T D B
            });
            // 上面的简化式容易出错，重写：
            // 清空再重新计算
            for (int i = 0; i < 2 * NN; i++)
                for (int j = 0; j < 2 * NN; j++) K[i, j] = 0.0;
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
                double totArea = ComputeArea() * Thickness * Density;
                // 角点 1/12，中点 0? 简化：行和法
                double[] rowSum = new double[NN];
                ForEachGaussPoint(delegate (double[] N, double[] dx, double[] dy, double w)
                {
                    for (int a = 0; a < NN; a++)
                        for (int b = 0; b < NN; b++)
                            rowSum[a] += Density * w * N[a] * N[b];
                });
                double total = 0;
                for (int i = 0; i < NN; i++) total += rowSum[i];
                for (int a = 0; a < NN; a++)
                {
                    double m = rowSum[a] * (totArea / total);
                    M[2 * a, 2 * a] = m;
                    M[2 * a + 1, 2 * a + 1] = m;
                }
            }
            else
            {
                ForEachGaussPoint(delegate (double[] N, double[] dx, double[] dy, double w)
                {
                    for (int a = 0; a < NN; a++)
                        for (int b = 0; b < NN; b++)
                        {
                            double v = Density * w * N[a] * N[b];
                            M[2 * a, 2 * b] += v;
                            M[2 * a + 1, 2 * b + 1] += v;
                        }
                });
            }
        }

        public override void ComputeStress(double[] u, double[] v,
            out double sxx, out double syy, out double sxy, out double vonMises)
        {
            double d00, d01, d11, d22;
            DPlaneStress(out d00, out d01, out d11, out d22);
            sxx = 0; syy = 0; sxy = 0; vonMises = 0;
            double totalW = 0;
            double[] N = new double[NN];
            double[] dNdL1 = new double[NN], dNdL2 = new double[NN];
            double[] dNdx = new double[NN], dNdy = new double[NN];
            for (int g = 0; g < 3; g++)
            {
                double L1 = _gps[g][0], L2 = _gps[g][1], L3 = _gps[g][2], wg = _gps[g][3];
                ShapeAndDeriv(L1, L2, L3, N, dNdL1, dNdL2);
                double J = Jacobian(L1, L2, L3, dNdL1, dNdL2, dNdx, dNdy);
                double wdet = wg * Math.Abs(J) * Thickness;
                double ex = 0, ey = 0, gxy = 0;
                for (int a = 0; a < NN; a++)
                {
                    ex += dNdx[a] * u[a];
                    ey += dNdy[a] * v[a];
                    gxy += dNdy[a] * u[a] + dNdx[a] * v[a];
                }
                sxx += wdet * (d00 * ex + d01 * ey);
                syy += wdet * (d01 * ex + d11 * ey);
                sxy += wdet * d22 * gxy;
                totalW += wdet;
            }
            sxx /= totalW; syy /= totalW; sxy /= totalW;
            vonMises = Math.Sqrt(sxx * sxx - sxx * syy + syy * syy + 3.0 * sxy * sxy);
        }

        public override void BuildThermalLoadVector(double[] TNodal, double TRef, double[] f)
        {
            double d00, d01, d11, d22;
            DPlaneStress(out d00, out d01, out d11, out d22);
            double alpha = ThermalExpansion;
            ForEachGaussPoint(delegate (double[] N, double[] dNdx, double[] dNdy, double w)
            {
                double dT = 0.0;
                for (int a = 0; a < NN; a++) dT += N[a] * (TNodal[a] - TRef);
                double et0 = alpha * dT;
                // σ0 = D·[εt,εt,0]^T
                double sx = d00 * et0 + d01 * et0;
                double sy = d01 * et0 + d11 * et0;
                for (int a = 0; a < NN; a++)
                {
                    f[2 * a] += w * (dNdx[a] * sx + dNdy[a] * 0 /*sxy*/);
                    f[2 * a + 1] += w * (dNdy[a] * sy + dNdx[a] * 0);
                }
            });
        }

        public override void BuildConductivityMatrix(double[,] Kt)
        {
            double k = Conductivity, t = Thickness;
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
                        for (int b = 0; b < NN; b++)
                            rs[a] += m0 * w * N[a] * N[b];
                });
                double total = 0, rowS = 0;
                for (int i = 0; i < NN; i++) { total += rs[i]; rowS += 1.0; }
                double area = ComputeArea() * Thickness;
                double factor = m0 * area / total;
                for (int a = 0; a < NN; a++) C[a, a] = rs[a] * factor * NN;
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
            double L1 = xi, L2 = eta, L3 = 1.0 - xi - eta;
            double[] N = new double[NN];
            N[0] = L1 * (2 * L1 - 1);
            N[1] = L2 * (2 * L2 - 1);
            N[2] = L3 * (2 * L3 - 1);
            N[3] = 4 * L1 * L2;
            N[4] = 4 * L2 * L3;
            N[5] = 4 * L3 * L1;
            return N;
        }

        public override double[,] ShapeDerivLocal(double xi, double eta)
        {
            double L1 = xi, L2 = eta, L3 = 1.0 - xi - eta;
            double[] N = new double[NN];
            double[] d1 = new double[NN], d2 = new double[NN];
            ShapeAndDeriv(L1, L2, L3, N, d1, d2);
            double[,] dd = new double[NN, 2];
            for (int i = 0; i < NN; i++) { dd[i, 0] = d1[i]; dd[i, 1] = d2[i]; }
            return dd;
        }
    }
}
