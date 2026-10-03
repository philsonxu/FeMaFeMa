using System;

namespace NonlinearFEM2D.Elements
{
    public sealed class Q4Element : FiniteElement
    {
        // 4-node bilinear quad, 2x2 Gauss
        private static readonly double[] GP = { -1.0 / Math.Sqrt(3.0), 1.0 / Math.Sqrt(3.0) };

        public Q4Element()
        {
            this.NodesPerElement = 4;
            this.NumGaussPoints = 4;
            this.Type = ElementType.Q4;
        }

        private static double N1(double xi, double eta) { return 0.25 * (1.0 - xi) * (1.0 - eta); }
        private static double N2(double xi, double eta) { return 0.25 * (1.0 + xi) * (1.0 - eta); }
        private static double N3(double xi, double eta) { return 0.25 * (1.0 + xi) * (1.0 + eta); }
        private static double N4(double xi, double eta) { return 0.25 * (1.0 - xi) * (1.0 + eta); }
        private static void EvalShape(double xi, double eta, double[] N, double[] dNdxi, double[] dNdeta)
        {
            N[0] = 0.25 * (1.0 - xi) * (1.0 - eta);
            N[1] = 0.25 * (1.0 + xi) * (1.0 - eta);
            N[2] = 0.25 * (1.0 + xi) * (1.0 + eta);
            N[3] = 0.25 * (1.0 - xi) * (1.0 + eta);
            dNdxi[0] = -0.25 * (1.0 - eta); dNdxi[1] = 0.25 * (1.0 - eta); dNdxi[2] = 0.25 * (1.0 + eta); dNdxi[3] = -0.25 * (1.0 + eta);
            dNdeta[0] = -0.25 * (1.0 - xi); dNdeta[1] = -0.25 * (1.0 + xi); dNdeta[2] = 0.25 * (1.0 + xi); dNdeta[3] = 0.25 * (1.0 - xi);
        }

        private bool Jacobian(double xi, double eta, out double J00, out double J01, out double J10, out double J11, out double detJ, out double[] dNdx, out double[] dNdy)
        {
            double[] N = new double[4], dNdxi = new double[4], dNdeta = new double[4];
            EvalShape(xi, eta, N, dNdxi, dNdeta);
            J00 = 0; J01 = 0; J10 = 0; J11 = 0;
            for (int i = 0; i < 4; i++)
            {
                J00 += dNdxi[i] * this.X[i]; J01 += dNdxi[i] * this.Y[i];
                J10 += dNdeta[i] * this.X[i]; J11 += dNdeta[i] * this.Y[i];
            }
            detJ = J00 * J11 - J01 * J10;
            dNdx = new double[4]; dNdy = new double[4];
            if (Math.Abs(detJ) < 1.0e-20) return false;
            double invDet = 1.0 / detJ;
            for (int i = 0; i < 4; i++)
            {
                dNdx[i] = (J11 * dNdxi[i] - J01 * dNdeta[i]) * invDet;
                dNdy[i] = (-J10 * dNdxi[i] + J00 * dNdeta[i]) * invDet;
            }
            return true;
        }

        public override double ComputeArea()
        {
            double a = 0.0;
            for (int ig = 0; ig < 2; ig++) for (int jg = 0; jg < 2; jg++)
                {
                    double xi = GP[ig], eta = GP[jg];
                    double J00, J01, J10, J11, detJ; double[] dx, dy;
                    if (Jacobian(xi, eta, out J00, out J01, out J10, out J11, out detJ, out dx, out dy)) a += Math.Abs(detJ);
                }
            return a;
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

        public override void BuildLinearStiffness(double[,] Ke)
        {
            Array.Clear(Ke, 0, Ke.Length);
            double[,] D = Dmatrix();
            double t = this.Thickness;
            for (int ig = 0; ig < 2; ig++) for (int jg = 0; jg < 2; jg++)
                {
                    double xi = GP[ig], eta = GP[jg];
                    double J00, J01, J10, J11, detJ; double[] dNdx, dNdy;
                    if (!Jacobian(xi, eta, out J00, out J01, out J10, out J11, out detJ, out dNdx, out dNdy)) continue;
                    double w = Math.Abs(detJ) * t;
                    double[,] B = new double[3, 8];
                    for (int i = 0; i < 4; i++)
                    {
                        B[0, 2 * i] = dNdx[i];
                        B[1, 2 * i + 1] = dNdy[i];
                        B[2, 2 * i] = dNdy[i]; B[2, 2 * i + 1] = dNdx[i];
                    }
                    for (int i = 0; i < 8; i++) for (int j = 0; j < 8; j++)
                        {
                            double s = 0.0;
                            for (int m = 0; m < 3; m++) for (int n = 0; n < 3; n++) s += B[m, i] * D[m, n] * B[n, j];
                            Ke[i, j] += s * w;
                        }
                }
        }

        public override void BuildTangentAndInternal(double[] u, double[] v, double[,] Ke, double[] fint, bool geomNonlin)
        {
            BuildLinearStiffness(Ke);
            Array.Clear(fint, 0, fint.Length);
            double[,] D = Dmatrix();
            double t = this.Thickness;
            int gidx = 0;
            for (int ig = 0; ig < 2; ig++) for (int jg = 0; jg < 2; jg++)
                {
                    double xi = GP[ig], eta = GP[jg];
                    double J00, J01, J10, J11, detJ; double[] dNdx, dNdy;
                    if (!Jacobian(xi, eta, out J00, out J01, out J10, out J11, out detJ, out dNdx, out dNdy)) continue;
                    double w = Math.Abs(detJ) * t;
                    double exx = 0, eyy = 0, gxy = 0;
                    for (int i = 0; i < 4; i++) { exx += dNdx[i] * u[i]; eyy += dNdy[i] * v[i]; gxy += dNdy[i] * u[i] + dNdx[i] * v[i]; }
                    double sxx = D[0, 0] * exx + D[0, 1] * eyy;
                    double syy = D[1, 0] * exx + D[1, 1] * eyy;
                    double sxy = D[2, 2] * gxy;
                    this.Histories[gidx].Stress[0] = sxx; this.Histories[gidx].Stress[1] = syy; this.Histories[gidx].Stress[3] = sxy;
                    this.Histories[gidx].VonMises = Math.Sqrt(Math.Max(0.0, sxx * sxx - sxx * syy + syy * syy + 3.0 * sxy * sxy));
                    for (int i = 0; i < 8; i++)
                    {
                        double f = 0.0;
                        int ii = i / 2;
                        f += (i % 2 == 0 ? dNdx[ii] : 0.0) * sxx;
                        f += (i % 2 == 1 ? dNdy[ii] : 0.0) * syy;
                        f += (i % 2 == 0 ? dNdy[ii] : dNdx[ii]) * sxy;
                        fint[i] += f * w;
                    }
                    gidx++;
                }
        }

        public override void BuildMassMatrix(double[,] Me)
        {
            Array.Clear(Me, 0, Me.Length);
            double t = this.Thickness, rho = this.Density;
            double[] N = new double[4], dx = new double[4], dy = new double[4];
            for (int ig = 0; ig < 2; ig++) for (int jg = 0; jg < 2; jg++)
                {
                    double xi = GP[ig], eta = GP[jg];
                    double J00, J01, J10, J11, detJ; double[] dNdx, dNdy;
                    if (!Jacobian(xi, eta, out J00, out J01, out J10, out J11, out detJ, out dNdx, out dNdy)) continue;
                    double w = Math.Abs(detJ) * t * rho;
                    // 重新计算N
                    N[0] = 0.25 * (1.0 - xi) * (1.0 - eta);
                    N[1] = 0.25 * (1.0 + xi) * (1.0 - eta);
                    N[2] = 0.25 * (1.0 + xi) * (1.0 + eta);
                    N[3] = 0.25 * (1.0 - xi) * (1.0 + eta);
                    for (int a = 0; a < 4; a++) for (int b = 0; b < 4; b++)
                        { Me[2 * a, 2 * b] += N[a] * N[b] * w; Me[2 * a + 1, 2 * b + 1] += N[a] * N[b] * w; }
                }
        }

        public override void BuildConvection(double[] uc, double[] vc, double[,] C)
        {
            Array.Clear(C, 0, C.Length);
            double t = this.Thickness;
            double[] N = new double[4];
            for (int ig = 0; ig < 2; ig++) for (int jg = 0; jg < 2; jg++)
                {
                    double xi = GP[ig], eta = GP[jg];
                    double J00, J01, J10, J11, detJ; double[] dNdx, dNdy;
                    if (!Jacobian(xi, eta, out J00, out J01, out J10, out J11, out detJ, out dNdx, out dNdy)) continue;
                    N[0] = 0.25 * (1.0 - xi) * (1.0 - eta);
                    N[1] = 0.25 * (1.0 + xi) * (1.0 - eta);
                    N[2] = 0.25 * (1.0 + xi) * (1.0 + eta);
                    N[3] = 0.25 * (1.0 - xi) * (1.0 + eta);
                    double uavg = 0, vavg = 0;
                    for (int i = 0; i < 4; i++) { uavg += N[i] * uc[i]; vavg += N[i] * vc[i]; }
                    double w = Math.Abs(detJ) * t;
                    for (int a = 0; a < 4; a++) for (int b = 0; b < 4; b++)
                        {
                            double cv = N[a] * (uavg * dNdx[b] + vavg * dNdy[b]) * w;
                            C[2 * a, 2 * b] += cv; C[2 * a + 1, 2 * b + 1] += cv;
                        }
                }
        }

        public override void BuildViscousMatrix(double[,] Kv)
        {
            Array.Clear(Kv, 0, Kv.Length);
            double mu = this.Viscosity, t = this.Thickness;
            for (int ig = 0; ig < 2; ig++) for (int jg = 0; jg < 2; jg++)
                {
                    double xi = GP[ig], eta = GP[jg];
                    double J00, J01, J10, J11, detJ; double[] dNdx, dNdy;
                    if (!Jacobian(xi, eta, out J00, out J01, out J10, out J11, out detJ, out dNdx, out dNdy)) continue;
                    double w = Math.Abs(detJ) * mu * t;
                    for (int a = 0; a < 4; a++) for (int b = 0; b < 4; b++)
                        {
                            double v = (dNdx[a] * dNdx[b] + dNdy[a] * dNdy[b]) * w;
                            Kv[2 * a, 2 * b] += v; Kv[2 * a + 1, 2 * b + 1] += v;
                        }
                }
        }

        public override void BuildPressureGradient(double[,] G)
        {
            Array.Clear(G, 0, G.Length);
            double t = this.Thickness;
            double[] N = new double[4];
            for (int ig = 0; ig < 2; ig++) for (int jg = 0; jg < 2; jg++)
                {
                    double xi = GP[ig], eta = GP[jg];
                    double J00, J01, J10, J11, detJ; double[] dNdx, dNdy;
                    if (!Jacobian(xi, eta, out J00, out J01, out J10, out J11, out detJ, out dNdx, out dNdy)) continue;
                    N[0] = 0.25 * (1.0 - xi) * (1.0 - eta);
                    N[1] = 0.25 * (1.0 + xi) * (1.0 - eta);
                    N[2] = 0.25 * (1.0 + xi) * (1.0 + eta);
                    N[3] = 0.25 * (1.0 - xi) * (1.0 + eta);
                    double w = Math.Abs(detJ) * t;
                    for (int a = 0; a < 4; a++) for (int bp = 0; bp < 4; bp++)
                        {
                            G[2 * a, bp] -= N[a] * dNdx[bp] * w;
                            G[2 * a + 1, bp] -= N[a] * dNdy[bp] * w;
                        }
                }
        }

        public override void ComputeStress(double[] u, double[] v, out double sxx, out double syy, out double sxy)
        {
            sxx = 0; syy = 0; sxy = 0;
            for (int g = 0; g < 4; g++)
            {
                sxx += this.Histories[g].Stress[0]; syy += this.Histories[g].Stress[1]; sxy += this.Histories[g].Stress[3];
            }
            sxx /= 4; syy /= 4; sxy /= 4;
        }
    }
}
