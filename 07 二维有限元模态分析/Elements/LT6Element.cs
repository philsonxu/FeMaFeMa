namespace ModalFEM2D.Elements
{
    /// <summary>
    /// 6-node quadratic triangle (LT6 / P2). Three-point Hammer quadrature (exact for quadratic).
    /// Node order in natural coords (xi, eta) with L1=1-xi-eta, L2=xi, L3=eta:
    ///   0: (0,0) = corner L1, 1: (1,0) = corner L2, 2: (0,1) = corner L3,
    ///   3: (0.5,0) = midside 0-1, 4: (0.5,0.5) = midside 1-2, 5: (0,0.5) = midside 2-0.
    /// </summary>
    public sealed class LT6Element : FiniteElement
    {
        public override ElementType Type { get { return ElementType.LT6; } }
        public override int NodesPerElement { get { return 6; } }
        public override int PressureNodes { get { return 6; } }

        protected override void GetIntegrationPoints(out double[,] points, out double[] weights)
        {
            points = new double[3, 2];
            double a = 1.0 / 6.0;
            points[0, 0] = a; points[0, 1] = a;
            points[1, 0] = 2.0 / 3.0; points[1, 1] = a;
            points[2, 0] = a; points[2, 1] = 2.0 / 3.0;
            weights = new double[3] { 1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0 };
        }

        protected override void ShapeAndDerivs(double xi, double eta, out double[] N, out double[,] dN)
        {
            double L1 = 1.0 - xi - eta;
            double L2 = xi;
            double L3 = eta;
            N = new double[6];
            N[0] = L1 * (2.0 * L1 - 1.0);
            N[1] = L2 * (2.0 * L2 - 1.0);
            N[2] = L3 * (2.0 * L3 - 1.0);
            N[3] = 4.0 * L1 * L2;
            N[4] = 4.0 * L2 * L3;
            N[5] = 4.0 * L3 * L1;
            // derivatives w.r.t (xi, eta): dL1/dxi=-1, dL1/deta=-1; dL2/dxi=1,dL2/deta=0; dL3/dxi=0,dL3/deta=1
            dN = new double[6, 2];
            dN[0, 0] = -(4.0 * L1 - 1.0); dN[0, 1] = -(4.0 * L1 - 1.0);
            dN[1, 0] =  (4.0 * L2 - 1.0); dN[1, 1] = 0.0;
            dN[2, 0] = 0.0;               dN[2, 1] =  (4.0 * L3 - 1.0);
            dN[3, 0] = 4.0 * (-L2 + L1);  dN[3, 1] = -4.0 * L2;
            dN[4, 0] = 4.0 * L3;          dN[4, 1] = 4.0 * L2;
            dN[5, 0] = -4.0 * L3;         dN[5, 1] = 4.0 * (L1 - L3);
        }
    }
}
