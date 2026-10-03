using System;

namespace ModalFEM2D.Elements
{
    /// <summary>
    /// Element type tag.
    /// </summary>
    public enum ElementType
    {
        CST3,
        LT6,
        Q4
    }

    /// <summary>
    /// Abstract base class for 2D finite elements. Coordinates and material are stored per instance.
    /// All assembly uses dense local matrices returned by Build*Matrix(), so elements can plug into
    /// structural, modal, or NS solvers through a uniform API.
    /// </summary>
    public abstract class FiniteElement
    {
        public int Id;
        public int[] NodeIds;
        public double[,] Coords;    // [npe,2]
        public double YoungModulus;
        public double PoissonRatio;
        public double Density;
        public double Viscosity;
        public double Thickness;
        public double StressXX;
        public double StressYY;
        public double StressXY;
        public double VonMises;

        public abstract ElementType Type { get; }
        public abstract int NodesPerElement { get; }
        public abstract int PressureNodes { get; }   // number of pressure DOFs per element (0 for structural)

        protected int Npe { get { return NodesPerElement; } }

        protected double[,] GetCoords()
        {
            if (Coords != null && Coords.GetLength(0) == Npe) return Coords;
            double[,] c = new double[Npe, 2];
            for (int i = 0; i < Npe; i++) { c[i, 0] = Coords[i, 0]; c[i, 1] = Coords[i, 1]; }
            return c;
        }

        // ===== Shape functions / geometry =====
        // points: [ngp,2] in natural coords; weights: [ngp]; dNdxi: [ngp,npe,2]
        protected abstract void GetIntegrationPoints(out double[,] points, out double[] weights);
        protected abstract void ShapeAndDerivs(double xi, double eta, out double[] N, out double[,] dN);

        protected static double Det2(double a, double b, double c, double d) { return a * d - b * c; }

        // Jacobian J = dN^T * coords; returns detJ
        protected static double Jacobian(double[,] dNdx, double[,] coords, out double[,] J, out double[,] invJ)
        {
            int npe = dNdx.GetLength(0);
            J = new double[2, 2];
            for (int i = 0; i < npe; i++)
            {
                J[0, 0] += dNdx[i, 0] * coords[i, 0]; J[0, 1] += dNdx[i, 0] * coords[i, 1];
                J[1, 0] += dNdx[i, 1] * coords[i, 0]; J[1, 1] += dNdx[i, 1] * coords[i, 1];
            }
            double det = J[0, 0] * J[1, 1] - J[0, 1] * J[1, 0];
            invJ = new double[2, 2];
            double s = 1.0 / det;
            invJ[0, 0] = J[1, 1] * s; invJ[0, 1] = -J[0, 1] * s;
            invJ[1, 0] = -J[1, 0] * s; invJ[1, 1] = J[0, 0] * s;
            return det;
        }

        // global B for plane strain elasticity at a Gauss point: B[3,2*npe]; dN is dN/dxi [npe,2]
        protected void StrainDisplacement(double[,] dN, double[,] invJ, out double[,] B)
        {
            int npe = dN.GetLength(0);
            double[,] dNdxy = new double[npe, 2];
            for (int i = 0; i < npe; i++)
            {
                dNdxy[i, 0] = invJ[0, 0] * dN[i, 0] + invJ[0, 1] * dN[i, 1];
                dNdxy[i, 1] = invJ[1, 0] * dN[i, 0] + invJ[1, 1] * dN[i, 1];
            }
            B = new double[3, 2 * npe];
            for (int i = 0; i < npe; i++)
            {
                B[0, 2 * i] = dNdxy[i, 0];
                B[1, 2 * i + 1] = dNdxy[i, 1];
                B[2, 2 * i] = dNdxy[i, 1];
                B[2, 2 * i + 1] = dNdxy[i, 0];
            }
        }

        // Plane-strain constitutive matrix D[3,3]
        protected void PlaneStrainD(out double[,] D)
        {
            double E = YoungModulus; double nu = PoissonRatio;
            double c = E / ((1.0 + nu) * (1.0 - 2.0 * nu));
            D = new double[3, 3];
            D[0, 0] = c * (1.0 - nu); D[0, 1] = c * nu;
            D[1, 0] = c * nu; D[1, 1] = c * (1.0 - nu);
            D[2, 2] = c * (0.5 - nu);
        }

        // ===== Public API =====
        public double ComputeArea()
        {
            double[,] c = GetCoords();
            double[,] pts; double[] w;
            GetIntegrationPoints(out pts, out w);
            double area = 0.0;
            for (int g = 0; g < w.Length; g++)
            {
                double[] N; double[,] dN;
                ShapeAndDerivs(pts[g, 0], pts[g, 1], out N, out dN);
                double[,] J, invJ;
                double detJ = Jacobian(dN, c, out J, out invJ);
                area += Math.Abs(detJ) * w[g] * Thickness;
            }
            return area;
        }

        /// <summary>Stiffness matrix Ke[2*npe, 2*npe].</summary>
        public double[,] BuildStiffnessMatrix()
        {
            double[,] c = GetCoords();
            int nd = 2 * Npe;
            double[,] Ke = new double[nd, nd];
            double[,] D; PlaneStrainD(out D);
            double[,] pts; double[] w; GetIntegrationPoints(out pts, out w);
            for (int g = 0; g < w.Length; g++)
            {
                double[] N; double[,] dN;
                ShapeAndDerivs(pts[g, 0], pts[g, 1], out N, out dN);
                double[,] J, invJ; double detJ = Jacobian(dN, c, out J, out invJ);
                double[,] B; StrainDisplacement(dN, invJ, out B);
                // B^T D B * w*detJ*t
                double wt = w[g] * detJ * Thickness;
                // B^T D: [2npe,3]
                double[,] BtD = new double[nd, 3];
                for (int i = 0; i < nd; i++)
                    for (int k = 0; k < 3; k++)
                        for (int j = 0; j < 3; j++)
                            BtD[i, k] += B[j, i] * D[j, k];
                for (int i = 0; i < nd; i++)
                    for (int j = 0; j < nd; j++)
                    {
                        double s = 0.0;
                        for (int k = 0; k < 3; k++) s += BtD[i, k] * B[k, j];
                        Ke[i, j] += s * wt;
                    }
            }
            return Ke;
        }

        /// <summary>Consistent mass matrix Me[2*npe,2*npe].</summary>
        public double[,] BuildMassMatrix()
        {
            double[,] c = GetCoords();
            int nd = 2 * Npe;
            double[,] Me = new double[nd, nd];
            double[,] pts; double[] w; GetIntegrationPoints(out pts, out w);
            double rho = Density, t = Thickness;
            for (int g = 0; g < w.Length; g++)
            {
                double[] N; double[,] dN;
                ShapeAndDerivs(pts[g, 0], pts[g, 1], out N, out dN);
                double[,] J, invJ; double detJ = Jacobian(dN, c, out J, out invJ);
                double wt = w[g] * detJ * t * rho;
                for (int i = 0; i < Npe; i++)
                    for (int j = 0; j < Npe; j++)
                    {
                        double v = N[i] * N[j] * wt;
                        Me[2 * i, 2 * j] += v;
                        Me[2 * i + 1, 2 * j + 1] += v;
                    }
            }
            return Me;
        }

        /// <summary>Lumped (diagonal) mass vector length 2*npe.</summary>
        public double[] BuildLumpedMass()
        {
            double[,] Me = BuildMassMatrix();
            int nd = 2 * Npe;
            double[] m = new double[nd];
            for (int i = 0; i < nd; i++)
            {
                double s = 0.0;
                for (int j = 0; j < nd; j++) s += Me[i, j];
                m[i] = s;
            }
            return m;
        }

        /// <summary>Viscous (diffusion) matrix for Stokes/NS: Ce[nd,nd] (velocity-velocity).</summary>
        public double[,] BuildViscousMatrix()
        {
            double[,] c = GetCoords();
            int nd = 2 * Npe;
            double[,] Ce = new double[nd, nd];
            double mu = Viscosity, t = Thickness;
            double[,] pts; double[] w; GetIntegrationPoints(out pts, out w);
            for (int g = 0; g < w.Length; g++)
            {
                double[] N; double[,] dN;
                ShapeAndDerivs(pts[g, 0], pts[g, 1], out N, out dN);
                double[,] J, invJ; double detJ = Jacobian(dN, c, out J, out invJ);
                int npe = Npe;
                double[,] dNdxy = new double[npe, 2];
                for (int i = 0; i < npe; i++)
                {
                    dNdxy[i, 0] = invJ[0, 0] * dN[i, 0] + invJ[0, 1] * dN[i, 1];
                    dNdxy[i, 1] = invJ[1, 0] * dN[i, 0] + invJ[1, 1] * dN[i, 1];
                }
                double wt = mu * w[g] * detJ * t;
                // grad u . grad u  +  grad v . grad v (simplified Newtonian 2D 2x2 block diagonal)
                for (int i = 0; i < npe; i++)
                    for (int j = 0; j < npe; j++)
                    {
                        double d = (dNdxy[i, 0] * dNdxy[j, 0] + dNdxy[i, 1] * dNdxy[j, 1]) * wt;
                        Ce[2 * i, 2 * j] += d;
                        Ce[2 * i + 1, 2 * j + 1] += d;
                    }
            }
            return Ce;
        }

        /// <summary>
        /// Convection matrix for velocity (u·∇)u using advecting velocity at nodes (Picard linearization):
        /// KeConv[nd,nd] with K_ij = sum_g ( u_h · ∇N_i ) N_j |J| w t
        /// velu, velv are per-node velocities of length Npe (element-local nodes).
        /// </summary>
        public double[,] BuildConvectionMatrix(double[] velu, double[] velv)
        {
            double[,] c = GetCoords();
            int nd = 2 * Npe;
            int npe = Npe;
            double[,] Kc = new double[nd, nd];
            double rho = Density, t = Thickness;
            double[,] pts; double[] w; GetIntegrationPoints(out pts, out w);
            for (int g = 0; g < w.Length; g++)
            {
                double[] N; double[,] dN;
                ShapeAndDerivs(pts[g, 0], pts[g, 1], out N, out dN);
                double[,] J, invJ; double detJ = Jacobian(dN, c, out J, out invJ);
                double[,] dNdxy = new double[npe, 2];
                for (int i = 0; i < npe; i++)
                {
                    dNdxy[i, 0] = invJ[0, 0] * dN[i, 0] + invJ[0, 1] * dN[i, 1];
                    dNdxy[i, 1] = invJ[1, 0] * dN[i, 0] + invJ[1, 1] * dN[i, 1];
                }
                double uh = 0.0, vh = 0.0;
                for (int k = 0; k < npe; k++) { uh += N[k] * velu[k]; vh += N[k] * velv[k]; }
                double wt = rho * w[g] * detJ * t;
                double[] gN = new double[npe];
                for (int i = 0; i < npe; i++) gN[i] = uh * dNdxy[i, 0] + vh * dNdxy[i, 1];
                for (int i = 0; i < npe; i++)
                    for (int j = 0; j < npe; j++)
                    {
                        double v = gN[i] * N[j] * wt;
                        Kc[2 * i, 2 * j] += v;
                        Kc[2 * i + 1, 2 * j + 1] += v;
                    }
            }
            return Kc;
        }

        /// <summary>Pressure-gradient coupling G = ∫ B_p^T N_u dΩ for velocity block in SIMPLE (D [nd, pnpe]).
        /// For collocated elements, pressure DOFs coincide with velocity nodes (pnpe = Npe), but we stabilise separately.</summary>
        public virtual double[,] BuildPressureGradient()
        {
            int nd = 2 * Npe; int pn = Npe;
            double[,] G = new double[nd, pn];
            double[,] c = GetCoords(); double t = Thickness;
            double[,] pts; double[] w; GetIntegrationPoints(out pts, out w);
            for (int g = 0; g < w.Length; g++)
            {
                double[] N; double[,] dN;
                ShapeAndDerivs(pts[g, 0], pts[g, 1], out N, out dN);
                double[,] J, invJ; double detJ = Jacobian(dN, c, out J, out invJ);
                int npe = Npe;
                double[,] dNdxy = new double[npe, 2];
                for (int i = 0; i < npe; i++)
                {
                    dNdxy[i, 0] = invJ[0, 0] * dN[i, 0] + invJ[0, 1] * dN[i, 1];
                    dNdxy[i, 1] = invJ[1, 0] * dN[i, 0] + invJ[1, 1] * dN[i, 1];
                }
                double wt = w[g] * detJ * t;
                for (int i = 0; i < npe; i++)
                    for (int j = 0; j < npe; j++)
                    {
                        G[2 * i, j] += dNdxy[i, 0] * N[j] * wt;
                        G[2 * i + 1, j] += dNdxy[i, 1] * N[j] * wt;
                    }
            }
            return G;
        }

        /// <summary>Divergence operator D = -G^T for pressure equation (G was ∇N_p · N_u? For equal-order we use -∫ ∇·u N_p).</summary>
        public virtual double[,] BuildDivergence()
        {
            double[,] G = BuildPressureGradient();
            int nd = G.GetLength(0); int pn = G.GetLength(1);
            double[,] D = new double[pn, nd];
            for (int i = 0; i < pn; i++)
                for (int j = 0; j < nd; j++)
                    D[i, j] = -G[j, i];
            return D;
        }

        /// <summary>Pressure Laplacian stabilisation Lp[pn,pn] = ∫ ∇N_p · ∇N_p dΩ; scaled by small coefficient in SIMPLE.</summary>
        public virtual double[,] BuildPressureLaplacian()
        {
            int pn = Npe;
            double[,] L = new double[pn, pn];
            double[,] c = GetCoords(); double t = Thickness;
            double[,] pts; double[] w; GetIntegrationPoints(out pts, out w);
            for (int g = 0; g < w.Length; g++)
            {
                double[] N; double[,] dN;
                ShapeAndDerivs(pts[g, 0], pts[g, 1], out N, out dN);
                double[,] J, invJ; double detJ = Jacobian(dN, c, out J, out invJ);
                int npe = Npe;
                double[,] dNdxy = new double[npe, 2];
                for (int i = 0; i < npe; i++)
                {
                    dNdxy[i, 0] = invJ[0, 0] * dN[i, 0] + invJ[0, 1] * dN[i, 1];
                    dNdxy[i, 1] = invJ[1, 0] * dN[i, 0] + invJ[1, 1] * dN[i, 1];
                }
                double wt = w[g] * detJ * t;
                for (int i = 0; i < npe; i++)
                    for (int j = 0; j < npe; j++)
                        L[i, j] += (dNdxy[i, 0] * dNdxy[j, 0] + dNdxy[i, 1] * dNdxy[j, 1]) * wt;
            }
            return L;
        }

        /// <summary>Compute element-center stresses from displacement vector u[2*npe] (element local).</summary>
        public void ComputeStress(double[] ue)
        {
            double[,] c = GetCoords();
            double[,] D; PlaneStrainD(out D);
            // Evaluate at element center (natural coords (0,0) or centroid depending on element).
            double xic = 0.0, etac = 0.0;
            if (Type == ElementType.CST3) { xic = 1.0 / 3.0; etac = 1.0 / 3.0; }
            double[] N; double[,] dN;
            ShapeAndDerivs(xic, etac, out N, out dN);
            double[,] J, invJ; Jacobian(dN, c, out J, out invJ);
            double[,] B; StrainDisplacement(dN, invJ, out B);
            double[] eps = new double[3];
            int nd = 2 * Npe;
            for (int k = 0; k < 3; k++)
                for (int j = 0; j < nd; j++)
                    eps[k] += B[k, j] * ue[j];
            double[] sig = new double[3];
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                    sig[i] += D[i, j] * eps[j];
            StressXX = sig[0]; StressYY = sig[1]; StressXY = sig[2];
            VonMises = Math.Sqrt(sig[0] * sig[0] - sig[0] * sig[1] + sig[1] * sig[1] + 3.0 * sig[2] * sig[2]);
        }

        /// <summary>Interpolate field to a natural coordinate given nodal values field[npe].</summary>
        public double Interpolate(double xi, double eta, double[] f)
        {
            double[] N; double[,] dN;
            ShapeAndDerivs(xi, eta, out N, out dN);
            double s = 0.0;
            for (int i = 0; i < Npe; i++) s += N[i] * f[i];
            return s;
        }
    }
}
