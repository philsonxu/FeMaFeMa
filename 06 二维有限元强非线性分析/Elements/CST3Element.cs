using System;
using NonlinearFEM2D.Materials;

namespace NonlinearFEM2D.Elements
{
    public sealed class CST3Element : FiniteElement
    {
        public CST3Element()
        {
            this.NodesPerElement = 3;
            this.NumGaussPoints = 1;
            this.Type = ElementType.CST3;
        }

        // 形状函数对(x,y)偏导（CST为常数）
        private void GradN(out double[] b, out double[] c, out double A)
        {
            double x0 = this.X[0], y0 = this.Y[0];
            double x1 = this.X[1], y1 = this.Y[1];
            double x2 = this.X[2], y2 = this.Y[2];
            b = new double[3]; c = new double[3];
            b[0] = y1 - y2; b[1] = y2 - y0; b[2] = y0 - y1;
            c[0] = x2 - x1; c[1] = x0 - x2; c[2] = x1 - x0;
            A = 0.5 * (b[0] * c[1] - b[1] * c[0]);
        }

        public override double ComputeArea()
        {
            double[] b, c; double A; GradN(out b, out c, out A); return Math.Abs(A);
        }

        public override void BuildLinearStiffness(double[,] Ke)
        {
            int npe = 3;
            Array.Clear(Ke, 0, Ke.Length);
            double[] b, c; double A; GradN(out b, out c, out A);
            double E = this.YoungModulus, nu = this.PoissonRatio, t = this.Thickness;
            double cf = E * t / (4.0 * A * (1.0 - nu * nu));
            // 平面应变 D 矩阵 (3x3) [xx,yy,xy]
            double[,] D = new double[3, 3];
            D[0, 0] = 1.0; D[0, 1] = nu; D[1, 0] = nu; D[1, 1] = 1.0; D[2, 2] = (1.0 - nu) / 2.0;
            for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) D[i, j] *= (E / (1.0 - nu * nu));
            for (int i = 0; i < npe; i++)
            {
                double[] Bi = new double[] { b[i], 0, c[i] / 2.0 };
                double[] Bi2 = new double[] { 0, c[i], b[i] / 2.0 };
                for (int j = 0; j < npe; j++)
                {
                    double[] Bj = new double[] { b[j], 0, c[j] };
                    double[] Bj2 = new double[] { 0, c[j], b[j] };
                    double k1 = 0.0, k2 = 0.0;
                    for (int m = 0; m < 3; m++) { k1 += Bi[m] * (D[m, 0] * Bj[0] + D[m, 1] * Bj[1] + D[m, 2] * Bj[2]); k2 += Bi2[m] * (D[m, 0] * Bj[0] + D[m, 1] * Bj[1] + D[m, 2] * Bj[2]); }
                    Ke[2 * i, 2 * j] += k1 * t * Math.Abs(A);
                    Ke[2 * i, 2 * j + 1] += k2 * t * Math.Abs(A);
                }
                for (int j = 0; j < npe; j++)
                {
                    double[] Bj = new double[] { b[j], 0, c[j] };
                    double[] Bj2 = new double[] { 0, c[j], b[j] };
                    double k1 = 0.0, k2 = 0.0;
                    for (int m = 0; m < 3; m++) { k1 += Bi[m] * (D[m, 0] * Bj2[0] + D[m, 1] * Bj2[1] + D[m, 2] * Bj2[2]); k2 += Bi2[m] * (D[m, 0] * Bj2[0] + D[m, 1] * Bj2[1] + D[m, 2] * Bj2[2]); }
                    Ke[2 * i + 1, 2 * j] += k1 * t * Math.Abs(A);
                    Ke[2 * i + 1, 2 * j + 1] += k2 * t * Math.Abs(A);
                }
            }
        }

        // 非线性一致切线 + 内力（J2 弹塑性，可选几何非线性）
        public override void BuildTangentAndInternal(double[] u, double[] v, double[,] Ke, double[] fint, bool geomNonlin)
        {
            int npe = 3;
            Array.Clear(Ke, 0, Ke.Length);
            Array.Clear(fint, 0, fint.Length);
            double[] b, c; double A; GradN(out b, out c, out A);
            double t = this.Thickness;
            double absA = Math.Abs(A);
            double w = absA * t;
            // 位移梯度：H = sum_i grad(N_i) (u_i, v_i)
            double ux = 0.0, uy = 0.0, vx = 0.0, vy = 0.0;
            double invA = 1.0 / (2.0 * A);
            double[] dNdx = new double[3], dNdy = new double[3];
            for (int i = 0; i < 3; i++) { dNdx[i] = b[i] * invA; dNdy[i] = c[i] * invA; }
            for (int i = 0; i < npe; i++)
            {
                ux += dNdx[i] * u[i]; uy += dNdy[i] * u[i];
                vx += dNdx[i] * v[i]; vy += dNdy[i] * v[i];
            }
            // 小应变
            double exx = ux, eyy = vy, gxy = uy + vx;
            double ezz = -this.PoissonRatio / (1.0 - this.PoissonRatio) * (exx + eyy); // 平面应变弹性估计（用于材料增量应变输入）
            double[] deps = new double[4] { exx, eyy, ezz, gxy / 2.0 };
            // 更新应力与切线 Dalg[4,4]
            double[,] Dalg = new double[4, 4];
            GPHistory h = this.Histories[0];
            this.Material.UpdateStressAndTangent(h.Stress, deps, ref h.EqPlasticStrain, Dalg);
            h.VonMises = J2PlasticityMaterial.VonMises(h.Stress);
            h.Yielded = (h.EqPlasticStrain > 1.0e-12);
            double sxx = h.Stress[0], syy = h.Stress[1], szz = h.Stress[2], sxy = h.Stress[3];
            // 将 Dalg[4,4] 收缩到 2D 平面应变 D3[3,3]，忽略 zz 对 Ke 的影响（平面应变已在材料中处理）
            // 实际使用 xx/yy/xy 三行三列：
            double[,] D3 = new double[3, 3];
            D3[0, 0] = Dalg[0, 0]; D3[0, 1] = Dalg[0, 1]; D3[0, 2] = Dalg[0, 3];
            D3[1, 0] = Dalg[1, 0]; D3[1, 1] = Dalg[1, 1]; D3[1, 2] = Dalg[1, 3];
            D3[2, 0] = Dalg[3, 0]; D3[2, 1] = Dalg[3, 1]; D3[2, 2] = Dalg[3, 3];
            // B 矩阵 (3 x 2npe): B = [[dNdx,0],[0,dNdy],[dNdy,dNdx]]
            double[,] B = new double[3, 6];
            for (int i = 0; i < npe; i++)
            {
                B[0, 2 * i] = dNdx[i];
                B[1, 2 * i + 1] = dNdy[i];
                B[2, 2 * i] = dNdy[i]; B[2, 2 * i + 1] = dNdx[i];
            }
            // 内力 fint = ∫ B^T σ dΩ
            double[] sv = new double[3] { sxx, syy, sxy };
            for (int i = 0; i < 2 * npe; i++)
            {
                double s = 0.0;
                for (int m = 0; m < 3; m++) s += B[m, i] * sv[m];
                fint[i] = s * w;
            }
            // 切线 Kt = ∫ B^T D3 B dΩ + Kgeo（可选）
            double[,] Kt = new double[6, 6];
            for (int i = 0; i < 6; i++)
                for (int j = 0; j < 6; j++)
                {
                    double s = 0.0;
                    for (int m = 0; m < 3; m++) for (int n = 0; n < 3; n++) s += B[m, i] * D3[m, n] * B[n, j];
                    Kt[i, j] = s * w;
                }
            if (geomNonlin)
            {
                // 几何刚度：Kσ = ∫ G^T [σ] G dΩ, G(2x6) = [[dNdx],[dNdy]] 对每节点, [σ] = [[sxx,sxy],[sxy,syy]]
                for (int ia = 0; ia < npe; ia++)
                {
                    for (int ib = 0; ib < npe; ib++)
                    {
                        double gx = dNdx[ia] * dNdx[ib] * sxx + dNdx[ia] * dNdy[ib] * sxy
                                  + dNdy[ia] * dNdx[ib] * sxy + dNdy[ia] * dNdy[ib] * syy;
                        Kt[2 * ia, 2 * ib] += gx * w;
                        Kt[2 * ia + 1, 2 * ib + 1] += gx * w;
                    }
                }
            }
            for (int i = 0; i < 6; i++) for (int j = 0; j < 6; j++) Ke[i, j] = Kt[i, j];
        }

        public override void BuildMassMatrix(double[,] Me)
        {
            Array.Clear(Me, 0, Me.Length);
            double A = this.ComputeArea();
            double t = this.Thickness, rho = this.Density;
            double cm = rho * t * A / 12.0;
            for (int i = 0; i < 3; i++)
            {
                Me[2 * i, 2 * i] = 2.0 * cm;
                Me[2 * i + 1, 2 * i + 1] = 2.0 * cm;
                for (int j = 0; j < 3; j++)
                {
                    if (i == j) continue;
                    Me[2 * i, 2 * j] = cm;
                    Me[2 * i + 1, 2 * j + 1] = cm;
                }
            }
        }

        public override void BuildConvection(double[] uc, double[] vc, double[,] C)
        {
            // 集中近似：C[a,b] = ∫ N_a (u·∇N_b) dΩ，u 单元平均
            Array.Clear(C, 0, C.Length);
            double[] b, c; double A; GradN(out b, out c, out A);
            double absA = Math.Abs(A), inv2A = 1.0 / (2.0 * A);
            double uavg = (uc[0] + uc[1] + uc[2]) / 3.0;
            double vavg = (vc[0] + vc[1] + vc[2]) / 3.0;
            double t = this.Thickness;
            double[] dNx = new double[3], dNy = new double[3];
            for (int i = 0; i < 3; i++) { dNx[i] = b[i] * inv2A; dNy[i] = c[i] * inv2A; }
            double[] Ni = new double[3] { 1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0 }; // 单点高斯
            for (int a = 0; a < 3; a++)
            {
                for (int bb = 0; bb < 3; bb++)
                {
                    double cval = Ni[a] * (uavg * dNx[bb] + vavg * dNy[bb]) * t * absA;
                    C[2 * a, 2 * bb] += cval;
                    C[2 * a + 1, 2 * bb + 1] += cval;
                }
            }
        }

        public override void BuildViscousMatrix(double[,] Kv)
        {
            Array.Clear(Kv, 0, Kv.Length);
            double[] b, c; double A; GradN(out b, out c, out A);
            double absA = Math.Abs(A), inv2A = 1.0 / (2.0 * A);
            double mu = this.Viscosity, t = this.Thickness;
            double[] dNx = new double[3], dNy = new double[3];
            for (int i = 0; i < 3; i++) { dNx[i] = b[i] * inv2A; dNy[i] = c[i] * inv2A; }
            for (int a = 0; a < 3; a++) for (int bb = 0; bb < 3; bb++)
            {
                double v = mu * t * absA * (dNx[a] * dNx[bb] + dNy[a] * dNy[bb]);
                Kv[2 * a, 2 * bb] += v;
                Kv[2 * a + 1, 2 * bb + 1] += v;
            }
        }

        public override void BuildPressureGradient(double[,] G)
        {
            Array.Clear(G, 0, G.Length);
            double[] b, c; double A; GradN(out b, out c, out A);
            double absA = Math.Abs(A), inv2A = 1.0 / (2.0 * A);
            double t = this.Thickness;
            double[] dNx = new double[3], dNy = new double[3];
            for (int i = 0; i < 3; i++) { dNx[i] = b[i] * inv2A; dNy[i] = c[i] * inv2A; }
            for (int a = 0; a < 3; a++) for (int bp = 0; bp < 3; bp++)
            {
                double Ni = 1.0 / 3.0;
                G[2 * a, bp] = -Ni * dNx[bp] * t * absA;
                G[2 * a + 1, bp] = -Ni * dNy[bp] * t * absA;
            }
        }

        public override void ComputeStress(double[] u, double[] v, out double sxx, out double syy, out double sxy)
        {
            sxx = this.Histories[0].Stress[0];
            syy = this.Histories[0].Stress[1];
            sxy = this.Histories[0].Stress[3];
        }
    }
}
