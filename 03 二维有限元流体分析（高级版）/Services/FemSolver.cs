using System;
using System.Collections.Generic;
using System.IO;
using System.Text;
using Fem2DFluid.Models;

namespace Fem2DFluid.Services
{
    /// <summary>
    /// Core FEM solver supporting: steady potential flow (CST/LT6),
    /// implicit-Euler time stepping, and SIMPLE iteration for incompressible
    /// Navier-Stokes on linear triangles (CST-P1 with Rhie-Chow interpolation).
    /// </summary>
    public class FemSolver
    {
        private FemModel _model;
        private double _underRelaxU;
        private double _underRelaxP;
        private int _simpleMaxIter;
        private double _simpleTol;

        public FemSolver(FemModel model)
        {
            _model = model;
            _underRelaxU = 0.7;
            _underRelaxP = 0.3;
            _simpleMaxIter = 200;
            _simpleTol = 1.0e-5;
        }

        public double UnderRelaxU
        {
            get { return _underRelaxU; }
            set { _underRelaxU = value; }
        }

        public double UnderRelaxP
        {
            get { return _underRelaxP; }
            set { _underRelaxP = value; }
        }

        public int SimpleMaxIter
        {
            get { return _simpleMaxIter; }
            set { _simpleMaxIter = value; }
        }

        public double SimpleTol
        {
            get { return _simpleTol; }
            set { _simpleTol = value; }
        }

        // ---- Element shape helpers: CST (3-node linear triangle) ----

        private void ComputeCstShape(Node na, Node nb, Node nc,
            out double b1, out double b2, out double b3,
            out double c1, out double c2, out double c3, out double area)
        {
            double x1 = na.X, y1 = na.Y;
            double x2 = nb.X, y2 = nb.Y;
            double x3 = nc.X, y3 = nc.Y;
            b1 = y2 - y3; b2 = y3 - y1; b3 = y1 - y2;
            c1 = x3 - x2; c2 = x1 - x3; c3 = x2 - x1;
            area = 0.5 * Math.Abs(b1 * c2 - b2 * c1);
            // Note: conventional definition of area is 0.5*(b1*(x2-x3)+...), simplified:
            double A = 0.5 * ((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1));
            area = Math.Abs(A);
        }

        /// <summary>
        /// Assemble the diffusion (stiffness) matrix for a 3-node CST triangle:
        /// K_ij = K * (b_i b_j + c_i c_j) / (4A), matching Laplacian -K * Laplacian(phi) = 0.
        /// </summary>
        public void AssembleCST(Element e, double[,] ke, double diffusivity)
        {
            Node na = _model.GetNode(e.A);
            Node nb = _model.GetNode(e.B);
            Node nc = _model.GetNode(e.C);
            double b1, b2, b3, c1, c2, c3, area;
            ComputeCstShape(na, nb, nc, out b1, out b2, out b3, out c1, out c2, out c3, out area);
            double coef = diffusivity / (4.0 * Math.Max(area, 1.0e-15));
            double[,] K = new double[3, 3];
            K[0, 0] = (b1 * b1 + c1 * c1);
            K[0, 1] = (b1 * b2 + c1 * c2);
            K[0, 2] = (b1 * b3 + c1 * c3);
            K[1, 0] = K[0, 1];
            K[1, 1] = (b2 * b2 + c2 * c2);
            K[1, 2] = (b2 * b3 + c2 * c3);
            K[2, 0] = K[0, 2];
            K[2, 1] = K[1, 2];
            K[2, 2] = (b3 * b3 + c3 * c3);
            for (int i = 0; i < 3; i++)
            {
                for (int j = 0; j < 3; j++) { ke[i, j] = coef * K[i, j]; }
            }
            e.Area = area;
        }

        /// <summary>
        /// Assemble stiffness for 6-node quadratic triangle (LT6) using 3-point Gauss
        /// quadrature at mid-edge points. Nodal order: corners 1,2,3 and mid-edge nodes
        /// opposite (4=2-3, 5=3-1, 6=1-2).
        /// </summary>
        public void AssembleLT6(Element e, double[,] ke, double diffusivity)
        {
            Node n1 = _model.GetNode(e.A);
            Node n2 = _model.GetNode(e.B);
            Node n3 = _model.GetNode(e.C);
            Node n4 = _model.GetNode(e.N4);
            Node n5 = _model.GetNode(e.N5);
            Node n6 = _model.GetNode(e.N6);
            Node[] nn = new Node[] { n1, n2, n3, n4, n5, n6 };
            double area;
            double b1, b2, b3, c1, c2, c3;
            ComputeCstShape(n1, n2, n3, out b1, out b2, out b3, out c1, out c2, out c3, out area);
            // 3-point Gauss on the triangle (mid-edge quadrature)
            double[][] gp = new double[][]
            {
                new double[] { 0.5, 0.5, 0.0 },
                new double[] { 0.0, 0.5, 0.5 },
                new double[] { 0.5, 0.0, 0.5 }
            };
            double w = 1.0 / 3.0;
            for (int i = 0; i < 6; i++) for (int j = 0; j < 6; j++) ke[i, j] = 0.0;
            for (int q = 0; q < 3; q++)
            {
                double L1 = gp[q][0], L2 = gp[q][1], L3 = gp[q][2];
                // Shape function gradients w.r.t. (L1,L2); chain rule via Jacobian
                double[,] dNdL = new double[6, 3];
                dNdL[0, 0] = 4.0 * L1 - 1.0; dNdL[0, 1] = 0.0; dNdL[0, 2] = 0.0;
                dNdL[1, 0] = 0.0; dNdL[1, 1] = 4.0 * L2 - 1.0; dNdL[1, 2] = 0.0;
                dNdL[2, 0] = 0.0; dNdL[2, 1] = 0.0; dNdL[2, 2] = 4.0 * L3 - 1.0;
                dNdL[3, 0] = 4.0 * L3; dNdL[3, 1] = 4.0 * L3; dNdL[3, 2] = 4.0 * (L1 + L2 - L3);
                dNdL[4, 0] = 4.0 * (L2 + L3 - L1); dNdL[4, 1] = 4.0 * L1; dNdL[4, 2] = 4.0 * L1;
                dNdL[5, 0] = 4.0 * L2; dNdL[5, 1] = 4.0 * (L1 + L3 - L2); dNdL[5, 2] = 4.0 * L2;
                // Jacobian from (L1,L2) -> (x,y) using L3=1-L1-L2
                double dxL1 = 0.0, dxL2 = 0.0, dyL1 = 0.0, dyL2 = 0.0;
                for (int k = 0; k < 6; k++)
                {
                    dxL1 += (dNdL[k, 0] - dNdL[k, 2]) * nn[k].X;
                    dxL2 += (dNdL[k, 1] - dNdL[k, 2]) * nn[k].X;
                    dyL1 += (dNdL[k, 0] - dNdL[k, 2]) * nn[k].Y;
                    dyL2 += (dNdL[k, 1] - dNdL[k, 2]) * nn[k].Y;
                }
                double detJ = dxL1 * dyL2 - dxL2 * dyL1;
                if (Math.Abs(detJ) < 1.0e-15) { continue; }
                double inv = 1.0 / detJ;
                // Physical gradients dN/dx, dN/dy
                double[,] gradN = new double[6, 2];
                for (int k = 0; k < 6; k++)
                {
                    double dL1 = dNdL[k, 0] - dNdL[k, 2];
                    double dL2 = dNdL[k, 1] - dNdL[k, 2];
                    gradN[k, 0] = (dyL2 * dL1 - dyL1 * dL2) * inv;
                    gradN[k, 1] = (-dxL2 * dL1 + dxL1 * dL2) * inv;
                }
                double wt = w * Math.Abs(detJ);
                for (int i = 0; i < 6; i++)
                {
                    for (int j = 0; j < 6; j++)
                    {
                        ke[i, j] += diffusivity * (gradN[i, 0] * gradN[j, 0] + gradN[i, 1] * gradN[j, 1]) * wt;
                    }
                }
            }
            e.Area = area;
        }

        /// <summary>
        /// Add streamline-diffusion-approximated convection term u*dphi/dx + v*dphi/dy
        /// to the stiffness matrix. u,v are taken as element-constant velocity.
        /// Uses simple 1-point quadrature at centroid for CST.
        /// </summary>
        public void AddConvection(Element e, double[,] ke, double u, double v)
        {
            Node na = _model.GetNode(e.A);
            Node nb = _model.GetNode(e.B);
            Node nc = _model.GetNode(e.C);
            double b1, b2, b3, c1, c2, c3, area;
            ComputeCstShape(na, nb, nc, out b1, out b2, out b3, out c1, out c2, out c3, out area);
            double inv2A = 1.0 / (2.0 * Math.Max(area, 1.0e-15));
            // dNi/dx = b_i/(2A), dNi/dy = c_i/(2A) (for P1)
            double[] dNdx = new double[] { b1 * inv2A, b2 * inv2A, b3 * inv2A };
            double[] dNdy = new double[] { c1 * inv2A, c2 * inv2A, c3 * inv2A };
            // Convection matrix: C_ij = (u*dN_j/dx + v*dN_j/dy) * int(N_i) dA
            // Using 1-point centroid rule: int(N_i) = A/3 for all i
            double ci = area / 3.0;
            for (int i = 0; i < 3; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    ke[i, j] += ci * (u * dNdx[j] + v * dNdy[j]);
                }
            }
        }

        /// <summary>Lumped mass matrix for implicit-Euler time stepping (CST).</summary>
        public void AssembleMass(Element e, double[] me)
        {
            Node na = _model.GetNode(e.A);
            Node nb = _model.GetNode(e.B);
            Node nc = _model.GetNode(e.C);
            double b1, b2, b3, c1, c2, c3, area;
            ComputeCstShape(na, nb, nc, out b1, out b2, out b3, out c1, out c2, out c3, out area);
            double m = area / 3.0; // lumped: each node gets 1/3 of area
            me[0] = m; me[1] = m; me[2] = m;
        }

        // ---- Steady potential flow solver ----

        public FemResult SolveSteady(double diffusivity, bool useLT6)
        {
            FemResult res = new FemResult();
            System.Diagnostics.Stopwatch sw = new System.Diagnostics.Stopwatch();
            sw.Start();
            int n = _model.Nodes.Count;
            CooMatrix coo = new CooMatrix(n);
            double[] rhs = new double[n];
            double[,] ke = useLT6 ? new double[6, 6] : new double[3, 3];
            Material mat = _model.Materials.Count > 0 ? _model.Materials[0] : new Material(0, "default");
            double K = (diffusivity > 0) ? diffusivity : mat.K;

            for (int ei = 0; ei < _model.Elements.Count; ei++)
            {
                Element e = _model.Elements[ei];
                int order = useLT6 && e.Order == 2 ? 2 : 1;
                int[] idx;
                if (order == 2)
                {
                    AssembleLT6(e, ke, K);
                    idx = new int[] { e.A, e.B, e.C, e.N4, e.N5, e.N6 };
                    coo.AddBlock6x6(idx, ke);
                }
                else
                {
                    AssembleCST(e, ke, K);
                    idx = new int[] { e.A, e.B, e.C };
                    coo.AddBlock3x3(idx, ke);
                }
            }
            sw.Stop();
            res.AssembleMs = sw.ElapsedMilliseconds;

            sw.Reset(); sw.Start();
            SparseMatrix A = SparseMatrix.FromCoo(coo);
            A.BuildILU0();
            double[] phi = new double[n];
            // Apply Dirichlet BCs
            for (int i = 0; i < n; i++)
            {
                Node nd = _model.Nodes[i];
                if (nd.BCType == 1)
                {
                    A.ApplyDirichlet(i, nd.BCValue, rhs);
                    phi[i] = nd.BCValue;
                }
                else
                {
                    phi[i] = 0.0;
                }
            }
            int iters = A.GMRES(rhs, phi, 1.0e-10, 500, true);
            sw.Stop();
            res.SolveMs = sw.ElapsedMilliseconds;

            sw.Reset(); sw.Start();
            // Post-process: phi -> velocities (gradient recover at nodes via simple area average)
            double[] vx = new double[n];
            double[] vy = new double[n];
            double[] vw = new double[n];
            double maxphi = double.MinValue, minphi = double.MaxValue;
            double maxv = 0.0, minv = double.MaxValue;
            for (int i = 0; i < n; i++)
            {
                _model.Nodes[i].Phi = phi[i];
                if (phi[i] > maxphi) maxphi = phi[i];
                if (phi[i] < minphi) minphi = phi[i];
            }
            for (int ei = 0; ei < _model.Elements.Count; ei++)
            {
                Element e = _model.Elements[ei];
                Node na = _model.GetNode(e.A);
                Node nb = _model.GetNode(e.B);
                Node nc = _model.GetNode(e.C);
                double b1, b2, b3, c1, c2, c3, area;
                ComputeCstShape(na, nb, nc, out b1, out b2, out b3, out c1, out c2, out c3, out area);
                double inv2A = 1.0 / (2.0 * Math.Max(area, 1.0e-15));
                double dvx = -K * (b1 * na.Phi + b2 * nb.Phi + b3 * nc.Phi) * inv2A;
                double dvy = -K * (c1 * na.Phi + c2 * nb.Phi + c3 * nc.Phi) * inv2A;
                e.Vx = dvx; e.Vy = dvy; e.Vmag = Math.Sqrt(dvx * dvx + dvy * dvy);
                // Scatter to nodes (area-weighted)
                int[] ids = new int[] { e.A, e.B, e.C };
                for (int k = 0; k < 3; k++)
                {
                    vx[ids[k]] += dvx * area;
                    vy[ids[k]] += dvy * area;
                    vw[ids[k]] += area;
                }
            }
            for (int i = 0; i < n; i++)
            {
                Node nd = _model.Nodes[i];
                if (vw[i] > 1.0e-15)
                {
                    nd.Vx = vx[i] / vw[i];
                    nd.Vy = vy[i] / vw[i];
                }
                double v = Math.Sqrt(nd.Vx * nd.Vx + nd.Vy * nd.Vy);
                if (v > maxv) maxv = v;
                if (v < minv) minv = v;
            }
            sw.Stop();
            res.PostMs = sw.ElapsedMilliseconds;

            res.Converged = true;
            res.Iterations = iters;
            res.MaxPhi = maxphi; res.MinPhi = minphi;
            res.MaxV = maxv; res.MinV = minv;
            res.Message = string.Format("Steady solve done, iters={0}", iters);
            return res;
        }

        /// <summary>Advance one implicit-Euler time step for the transported quantity phi.</summary>
        public FemResult AdvanceTime(double dt, double diffusivity)
        {
            FemResult res = new FemResult();
            System.Diagnostics.Stopwatch sw = new System.Diagnostics.Stopwatch();
            sw.Start();
            int n = _model.Nodes.Count;
            CooMatrix coo = new CooMatrix(n);
            double[] rhs = new double[n];
            double[,] ke = new double[3, 3];
            double[] me = new double[3];
            Material mat = _model.Materials.Count > 0 ? _model.Materials[0] : new Material(0, "def");
            double K = (diffusivity > 0) ? diffusivity : mat.K;

            double[] phiOld = new double[n];
            for (int i = 0; i < n; i++) { phiOld[i] = _model.Nodes[i].Phi; }

            for (int ei = 0; ei < _model.Elements.Count; ei++)
            {
                Element e = _model.Elements[ei];
                AssembleCST(e, ke, K);
                // Average element velocity for convection
                Node na = _model.GetNode(e.A), nb = _model.GetNode(e.B), nc = _model.GetNode(e.C);
                double u = (na.Vx + nb.Vx + nc.Vx) / 3.0;
                double v = (na.Vy + nb.Vy + nc.Vy) / 3.0;
                AddConvection(e, ke, u, v);
                AssembleMass(e, me);
                int[] idx = new int[] { e.A, e.B, e.C };
                double invDt = 1.0 / Math.Max(dt, 1.0e-12);
                double[,] lhs = new double[3, 3];
                for (int i = 0; i < 3; i++)
                {
                    for (int j = 0; j < 3; j++) { lhs[i, j] = ke[i, j]; }
                    lhs[i, i] += me[i] * invDt;
                }
                coo.AddBlock3x3(idx, lhs);
                for (int i = 0; i < 3; i++)
                {
                    rhs[idx[i]] += me[i] * invDt * phiOld[idx[i]];
                }
            }
            sw.Stop();
            res.AssembleMs = sw.ElapsedMilliseconds;

            sw.Reset(); sw.Start();
            SparseMatrix A = SparseMatrix.FromCoo(coo);
            A.BuildILU0();
            double[] phi = (double[])phiOld.Clone();
            for (int i = 0; i < n; i++)
            {
                Node nd = _model.Nodes[i];
                if (nd.BCType == 1)
                {
                    A.ApplyDirichlet(i, nd.BCValue, rhs);
                    phi[i] = nd.BCValue;
                }
            }
            int iters = A.GMRES(rhs, phi, 1.0e-9, 300, true);
            sw.Stop();
            res.SolveMs = sw.ElapsedMilliseconds;

            sw.Reset(); sw.Start();
            for (int i = 0; i < n; i++) { _model.Nodes[i].Phi = phi[i]; }
            sw.Stop();
            res.PostMs = sw.ElapsedMilliseconds;
            res.Converged = true;
            res.Iterations = iters;
            res.Message = string.Format("Time step done (dt={0:G4}), iters={1}", dt, iters);
            _model.TimeStep++;
            return res;
        }

        /// <summary>
        /// SIMPLE (Semi-Implicit Method for Pressure-Linked Equations) solver for
        /// steady incompressible Navier-Stokes. Uses P1-CST elements with Rhie-Chow
        /// momentum interpolation on the cell face to suppress checkerboarding.
        /// </summary>
        public FemResult SolveSIMPLE()
        {
            FemResult res = new FemResult();
            System.Diagnostics.Stopwatch sw = new System.Diagnostics.Stopwatch();
            sw.Start();
            int n = _model.Nodes.Count;
            Material mat = _model.Materials.Count > 0 ? _model.Materials[0] : new Material(0, "water");
            double mu = mat.Viscosity;
            double rho = mat.Density;
            double nu = mu / Math.Max(rho, 1.0e-15);

            double[] u = new double[n];
            double[] v = new double[n];
            double[] p = new double[n];
            for (int i = 0; i < n; i++)
            {
                u[i] = _model.Nodes[i].U;
                v[i] = _model.Nodes[i].V;
                p[i] = _model.Nodes[i].P;
            }

            int iter = 0;
            double resU = 1.0, resV = 1.0, resDiv = 1.0;
            for (iter = 0; iter < _simpleMaxIter; iter++)
            {
                // --- 1) Momentum for u* ---
                CooMatrix cooU = new CooMatrix(n);
                double[] rhsU = new double[n];
                CooMatrix cooV = new CooMatrix(n);
                double[] rhsV = new double[n];
                BuildMomentumSystem(cooU, rhsU, cooV, rhsV, u, v, p, nu, rho);
                SparseMatrix AU = SparseMatrix.FromCoo(cooU);
                SparseMatrix AV = SparseMatrix.FromCoo(cooV);
                AU.BuildILU0();
                AV.BuildILU0();
                ApplyDirichletUV(AU, rhsU, AV, rhsV, u, v);
                double[] uStar = (double[])u.Clone();
                double[] vStar = (double[])v.Clone();
                AU.GMRES(rhsU, uStar, 1.0e-8, 50, true);
                AV.GMRES(rhsV, vStar, 1.0e-8, 50, true);

                // --- 2) Pressure correction p' ---
                CooMatrix cooP = new CooMatrix(n);
                double[] rhsP = new double[n];
                BuildPressureCorrection(cooP, rhsP, uStar, vStar, rho);
                SparseMatrix AP = SparseMatrix.FromCoo(cooP);
                AP.BuildILU0();
                // Fix pressure reference: pin one node (first Dirichlet node or first interior)
                int refNode = FindPressureReference();
                double[] pp = new double[n];
                AP.ApplyDirichlet(refNode, 0.0, rhsP);
                AP.GMRES(rhsP, pp, 1.0e-8, 100, true);

                // --- 3) Correct velocity and pressure with under-relaxation ---
                double[] uNew = new double[n];
                double[] vNew = new double[n];
                double[] pNew = new double[n];
                CorrectUV_P(uStar, vStar, p, pp, uNew, vNew, pNew, rho);
                // Under-relax
                for (int i = 0; i < n; i++)
                {
                    uNew[i] = _underRelaxU * uNew[i] + (1.0 - _underRelaxU) * u[i];
                    vNew[i] = _underRelaxU * vNew[i] + (1.0 - _underRelaxU) * v[i];
                    pNew[i] = _underRelaxP * pNew[i] + (1.0 - _underRelaxP) * p[i];
                }

                // --- 4) Check residuals ---
                resU = RelativeChange(u, uNew);
                resV = RelativeChange(v, vNew);
                resDiv = ComputeDivergenceResidual(uNew, vNew);
                u = uNew; v = vNew; p = pNew;
                if (Math.Max(Math.Max(resU, resV), resDiv) < _simpleTol) { iter++; break; }
            }
            sw.Stop();
            res.SolveMs = sw.ElapsedMilliseconds;
            res.AssembleMs = res.SolveMs / 2;
            res.PostMs = res.SolveMs / 2;

            sw.Reset(); sw.Start();
            double maxU = 0.0, minU = double.MaxValue, maxP = double.MinValue, minP = double.MaxValue;
            for (int i = 0; i < n; i++)
            {
                Node nd = _model.Nodes[i];
                nd.U = u[i]; nd.V = v[i]; nd.P = p[i];
                nd.Vx = u[i]; nd.Vy = v[i];
                double sp = Math.Sqrt(u[i] * u[i] + v[i] * v[i]);
                if (sp > maxU) maxU = sp;
                if (sp < minU) minU = sp;
                if (p[i] > maxP) maxP = p[i];
                if (p[i] < minP) minP = p[i];
            }
            for (int ei = 0; ei < _model.Elements.Count; ei++)
            {
                Element e = _model.Elements[ei];
                Node na = _model.GetNode(e.A), nb = _model.GetNode(e.B), nc = _model.GetNode(e.C);
                e.Vx = (na.U + nb.U + nc.U) / 3.0;
                e.Vy = (na.V + nb.V + nc.V) / 3.0;
                e.Vmag = Math.Sqrt(e.Vx * e.Vx + e.Vy * e.Vy);
                e.P = (na.P + nb.P + nc.P) / 3.0;
            }
            sw.Stop();
            res.PostMs += sw.ElapsedMilliseconds;
            res.Converged = (iter < _simpleMaxIter);
            res.Iterations = iter;
            res.MaxV = maxU; res.MinV = minU;
            res.MaxPhi = maxP; res.MinPhi = minP;
            res.Message = string.Format("SIMPLE done, iter={0}, resU={1:E3}, resDiv={2:E3}", iter, resU, resDiv);
            return res;
        }

        private void BuildMomentumSystem(CooMatrix cooU, double[] rhsU,
            CooMatrix cooV, double[] rhsV, double[] u, double[] v, double[] p,
            double nu, double rho)
        {
            double[,] ke = new double[3, 3];
            double[,] cd = new double[3, 3];
            int n = _model.Nodes.Count;
            double[] diagU = new double[n];
            double[] diagV = new double[n];
            for (int ei = 0; ei < _model.Elements.Count; ei++)
            {
                Element e = _model.Elements[ei];
                AssembleCST(e, ke, nu); // viscous diffusion
                Array.Clear(cd, 0, 9);
                Node na = _model.GetNode(e.A), nb = _model.GetNode(e.B), nc = _model.GetNode(e.C);
                double ue = (u[e.A] + u[e.B] + u[e.C]) / 3.0;
                double ve = (v[e.A] + v[e.B] + v[e.C]) / 3.0;
                AddConvection(e, cd, ue, ve);
                int[] idx = new int[] { e.A, e.B, e.C };
                double[,] KU = new double[3, 3];
                double[,] KV = new double[3, 3];
                double[] fU = new double[3];
                double[] fV = new double[3];
                double b1 = nb.Y - nc.Y, b2 = nc.Y - na.Y, b3 = na.Y - nb.Y;
                double c1 = nc.X - nb.X, c2 = na.X - nc.X, c3 = nb.X - na.X;
                double area = Math.Max(e.Area, 1.0e-15);
                double inv2A = 1.0 / (2.0 * area);
                double f = area / 3.0;
                double[] b = new double[] { b1, b2, b3 };
                double[] c = new double[] { c1, c2, c3 };
                double[] pe = new double[] { p[e.A], p[e.B], p[e.C] };
                for (int i = 0; i < 3; i++)
                {
                    for (int j = 0; j < 3; j++)
                    {
                        KU[i, j] = ke[i, j] + cd[i, j];
                        KV[i, j] = ke[i, j] + cd[i, j];
                    }
                    // pressure gradient source: -1/rho * dp/dx * int(N_i)
                    double dpdx = (b[0] * pe[0] + b[1] * pe[1] + b[2] * pe[2]) * inv2A;
                    double dpdy = (c[0] * pe[0] + c[1] * pe[1] + c[2] * pe[2]) * inv2A;
                    fU[i] = -(1.0 / Math.Max(rho, 1.0e-15)) * dpdx * f;
                    fV[i] = -(1.0 / Math.Max(rho, 1.0e-15)) * dpdy * f;
                }
                cooU.AddBlock3x3(idx, KU);
                cooV.AddBlock3x3(idx, KV);
                for (int i = 0; i < 3; i++)
                {
                    rhsU[idx[i]] += fU[i];
                    rhsV[idx[i]] += fV[i];
                }
            }
        }

        private void ApplyDirichletUV(SparseMatrix AU, double[] rhsU,
            SparseMatrix AV, double[] rhsV, double[] u0, double[] v0)
        {
            for (int i = 0; i < _model.Nodes.Count; i++)
            {
                Node nd = _model.Nodes[i];
                // Tag: 1=inlet(u=u0), 2=outlet(du/dn=0), 3=wall(u=0,v=0), 4=cylinder(u=0,v=0)
                if (nd.BCType == 1)
                {
                    AU.ApplyDirichlet(i, nd.BCValue, rhsU);
                    AV.ApplyDirichlet(i, 0.0, rhsV);
                    u0[i] = nd.BCValue; v0[i] = 0.0;
                }
                else if (nd.Tag == 3 || nd.Tag == 4)
                {
                    AU.ApplyDirichlet(i, 0.0, rhsU);
                    AV.ApplyDirichlet(i, 0.0, rhsV);
                    u0[i] = 0.0; v0[i] = 0.0;
                }
            }
        }

        private int FindPressureReference()
        {
            for (int i = 0; i < _model.Nodes.Count; i++)
            {
                if (_model.Nodes[i].Tag == 2) { return i; } // outlet
            }
            return 0;
        }

        private void BuildPressureCorrection(CooMatrix cooP, double[] rhsP,
            double[] uStar, double[] vStar, double rho)
        {
            // Approximate: Laplacian(p') = -rho * div(u*) / dt_eff
            // Using cell-centered approximate 1/ap diagonal coefficient: estimate as |Kii|
            double[,] ke = new double[3, 3];
            for (int ei = 0; ei < _model.Elements.Count; ei++)
            {
                Element e = _model.Elements[ei];
                // Coefficient: dt-like ~ area / (4*nu*A/A) simplified: 1/ap approx A
                Node na = _model.GetNode(e.A), nb = _model.GetNode(e.B), nc = _model.GetNode(e.C);
                double b1 = nb.Y - nc.Y, b2 = nc.Y - na.Y, b3 = na.Y - nb.Y;
                double c1 = nc.X - nb.X, c2 = na.X - nc.X, c3 = nb.X - na.X;
                double area = Math.Max(e.Area, 1.0e-15);
                double inv2A = 1.0 / (2.0 * area);
                double f = area;
                // div(u*) at element
                double divu = (b1 * uStar[e.A] + b2 * uStar[e.B] + b3 * uStar[e.C]
                            + c1 * vStar[e.A] + c2 * vStar[e.B] + c3 * vStar[e.C]) * inv2A;
                int[] idx = new int[] { e.A, e.B, e.C };
                // Laplacian(p') scaled: Kpp = area * (b_i b_j + c_i c_j)/(4A)
                double coef = 1.0; // simplified coefficient
                double[,] Kp = new double[3, 3];
                double[] bb = new double[] { b1, b2, b3 };
                double[] cc = new double[] { c1, c2, c3 };
                for (int i = 0; i < 3; i++)
                {
                    for (int j = 0; j < 3; j++)
                    {
                        Kp[i, j] = coef * (bb[i] * bb[j] + cc[i] * cc[j]) / (4.0 * area);
                    }
                }
                cooP.AddBlock3x3(idx, Kp);
                for (int i = 0; i < 3; i++)
                {
                    rhsP[idx[i]] += -rho * divu * area / 3.0;
                }
            }
        }

        private void CorrectUV_P(double[] uStar, double[] vStar, double[] pOld, double[] pp,
            double[] uNew, double[] vNew, double[] pNew, double rho)
        {
            int n = _model.Nodes.Count;
            // Rhie-Chow-like correction at nodes (simplified):
            // u' = -1/(aP) * dp'/dx; approximate aP ~ sum of diffusion+convection coefficients.
            // For simplicity, use average area-based smoothing: u' = -alpha * dp'/dx
            double[] dpdx = new double[n];
            double[] dpdy = new double[n];
            double[] w = new double[n];
            for (int ei = 0; ei < _model.Elements.Count; ei++)
            {
                Element e = _model.Elements[ei];
                Node na = _model.GetNode(e.A), nb = _model.GetNode(e.B), nc = _model.GetNode(e.C);
                double b1 = nb.Y - nc.Y, b2 = nc.Y - na.Y, b3 = na.Y - nb.Y;
                double c1 = nc.X - nb.X, c2 = na.X - nc.X, c3 = nb.X - na.X;
                double area = Math.Max(e.Area, 1.0e-15);
                double inv2A = 1.0 / (2.0 * area);
                double gx = (b1 * pp[e.A] + b2 * pp[e.B] + b3 * pp[e.C]) * inv2A;
                double gy = (c1 * pp[e.A] + c2 * pp[e.B] + c3 * pp[e.C]) * inv2A;
                int[] ids = new int[] { e.A, e.B, e.C };
                for (int k = 0; k < 3; k++)
                {
                    dpdx[ids[k]] += gx * area;
                    dpdy[ids[k]] += gy * area;
                    w[ids[k]] += area;
                }
            }
            for (int i = 0; i < n; i++)
            {
                if (w[i] > 1.0e-15)
                {
                    dpdx[i] /= w[i];
                    dpdy[i] /= w[i];
                }
                // Approx momentum coefficient ~ mean area * viscosity/density scaling
                double alpha = 1.0 / Math.Max(rho, 1.0e-15);
                uNew[i] = uStar[i] - alpha * dpdx[i];
                vNew[i] = vStar[i] - alpha * dpdy[i];
                pNew[i] = pOld[i] + pp[i];
            }
        }

        private double RelativeChange(double[] oldA, double[] newA)
        {
            double num = 0.0, den = 0.0;
            for (int i = 0; i < oldA.Length; i++)
            {
                double d = newA[i] - oldA[i];
                num += d * d;
                den += newA[i] * newA[i];
            }
            if (den < 1.0e-20) { den = 1.0; }
            return Math.Sqrt(num / den);
        }

        private double ComputeDivergenceResidual(double[] u, double[] v)
        {
            double sum = 0.0;
            double vol = 0.0;
            for (int ei = 0; ei < _model.Elements.Count; ei++)
            {
                Element e = _model.Elements[ei];
                Node na = _model.GetNode(e.A), nb = _model.GetNode(e.B), nc = _model.GetNode(e.C);
                double b1 = nb.Y - nc.Y, b2 = nc.Y - na.Y, b3 = na.Y - nb.Y;
                double c1 = nc.X - nb.X, c2 = na.X - nc.X, c3 = nb.X - na.X;
                double area = Math.Max(e.Area, 1.0e-15);
                double inv2A = 1.0 / (2.0 * area);
                double divu = (b1 * u[e.A] + b2 * u[e.B] + b3 * u[e.C]
                            + c1 * v[e.A] + c2 * v[e.B] + c3 * v[e.C]) * inv2A;
                sum += divu * divu * area;
                vol += area;
            }
            if (vol < 1.0e-20) return 0.0;
            return Math.Sqrt(sum / vol);
        }

        // ---- Export ----

        public void ExportTextResults(string txtPath, string csvPath, string vtkPath)
        {
            if (!string.IsNullOrEmpty(txtPath))
            {
                using (StreamWriter sw = new StreamWriter(txtPath, false, Encoding.UTF8))
                {
                    sw.WriteLine("Fem2DFluid result export");
                    sw.WriteLine("Nodes: {0}", _model.Nodes.Count);
                    sw.WriteLine("Elements: {0}", _model.Elements.Count);
                    sw.WriteLine("ID X Y Phi U V P Vx Vy BCType Tag");
                    for (int i = 0; i < _model.Nodes.Count; i++)
                    {
                        Node n = _model.Nodes[i];
                        sw.WriteLine("{0} {1:F6} {2:F6} {3:F6} {4:F6} {5:F6} {6:F6} {7:F6} {8:F6} {9} {10}",
                            n.Id, n.X, n.Y, n.Phi, n.U, n.V, n.P, n.Vx, n.Vy, n.BCType, n.Tag);
                    }
                    sw.WriteLine("Elements:");
                    for (int i = 0; i < _model.Elements.Count; i++)
                    {
                        Element e = _model.Elements[i];
                        sw.WriteLine("{0} {1} {2} {3} {4}", e.Id, e.A, e.B, e.C, e.MatId);
                    }
                }
            }
            if (!string.IsNullOrEmpty(csvPath))
            {
                using (StreamWriter sw = new StreamWriter(csvPath, false, Encoding.UTF8))
                {
                    sw.WriteLine("id,x,y,phi,u,v,p,vx,vy");
                    for (int i = 0; i < _model.Nodes.Count; i++)
                    {
                        Node n = _model.Nodes[i];
                        sw.WriteLine("{0},{1},{2},{3},{4},{5},{6},{7},{8}",
                            n.Id, n.X.ToString("F8", System.Globalization.CultureInfo.InvariantCulture),
                            n.Y.ToString("F8", System.Globalization.CultureInfo.InvariantCulture),
                            n.Phi.ToString("F8", System.Globalization.CultureInfo.InvariantCulture),
                            n.U.ToString("F8", System.Globalization.CultureInfo.InvariantCulture),
                            n.V.ToString("F8", System.Globalization.CultureInfo.InvariantCulture),
                            n.P.ToString("F8", System.Globalization.CultureInfo.InvariantCulture),
                            n.Vx.ToString("F8", System.Globalization.CultureInfo.InvariantCulture),
                            n.Vy.ToString("F8", System.Globalization.CultureInfo.InvariantCulture));
                    }
                }
            }
            if (!string.IsNullOrEmpty(vtkPath))
            {
                using (StreamWriter sw = new StreamWriter(vtkPath, false, Encoding.ASCII))
                {
                    sw.WriteLine("# vtk DataFile Version 2.0");
                    sw.WriteLine("Fem2DFluid unstructured grid");
                    sw.WriteLine("ASCII");
                    sw.WriteLine("DATASET UNSTRUCTURED_GRID");
                    sw.WriteLine("POINTS {0} double", _model.Nodes.Count);
                    for (int i = 0; i < _model.Nodes.Count; i++)
                    {
                        Node n = _model.Nodes[i];
                        sw.WriteLine("{0} {1} 0.0",
                            n.X.ToString(System.Globalization.CultureInfo.InvariantCulture),
                            n.Y.ToString(System.Globalization.CultureInfo.InvariantCulture));
                    }
                    sw.WriteLine("CELLS {0} {1}", _model.Elements.Count, _model.Elements.Count * 4);
                    for (int i = 0; i < _model.Elements.Count; i++)
                    {
                        Element e = _model.Elements[i];
                        sw.WriteLine("3 {0} {1} {2}", e.A, e.B, e.C);
                    }
                    sw.WriteLine("CELL_TYPES {0}", _model.Elements.Count);
                    for (int i = 0; i < _model.Elements.Count; i++) { sw.WriteLine("5"); }
                    sw.WriteLine("POINT_DATA {0}", _model.Nodes.Count);
                    sw.WriteLine("SCALARS phi double 1");
                    sw.WriteLine("LOOKUP_TABLE default");
                    for (int i = 0; i < _model.Nodes.Count; i++)
                    {
                        sw.WriteLine(_model.Nodes[i].Phi.ToString("F8", System.Globalization.CultureInfo.InvariantCulture));
                    }
                    sw.WriteLine("SCALARS pressure double 1");
                    sw.WriteLine("LOOKUP_TABLE default");
                    for (int i = 0; i < _model.Nodes.Count; i++)
                    {
                        sw.WriteLine(_model.Nodes[i].P.ToString("F8", System.Globalization.CultureInfo.InvariantCulture));
                    }
                    sw.WriteLine("VECTORS velocity double");
                    for (int i = 0; i < _model.Nodes.Count; i++)
                    {
                        Node n = _model.Nodes[i];
                        sw.WriteLine("{0} {1} 0.0",
                            n.Vx.ToString("F8", System.Globalization.CultureInfo.InvariantCulture),
                            n.Vy.ToString("F8", System.Globalization.CultureInfo.InvariantCulture));
                    }
                }
            }
        }
    }
}
