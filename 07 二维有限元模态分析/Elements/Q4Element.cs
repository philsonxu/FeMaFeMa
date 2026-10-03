namespace ModalFEM2D.Elements
{
    /// <summary>
    /// 4-node bilinear isoparametric quadrilateral (Q4).
    /// Natural coordinates xi, eta in [-1,1]^2. Node order:
    ///   0 (-1,-1), 1 (1,-1), 2 (1,1), 3 (-1,1) (counter-clockwise).
    /// Integration: 2x2 Gauss (exact for bilinear polynomials).
    /// </summary>
    public sealed class Q4Element : FiniteElement
    {
        private static readonly double GP = 1.0 / System.Math.Sqrt(3.0);

        public override ElementType Type { get { return ElementType.Q4; } }
        public override int NodesPerElement { get { return 4; } }
        public override int PressureNodes { get { return 4; } }

        protected override void GetIntegrationPoints(out double[,] points, out double[] weights)
        {
            points = new double[4, 2];
            weights = new double[4];
            int k = 0;
            for (int i = 0; i < 2; i++) for (int j = 0; j < 2; j++)
            {
                points[k, 0] = (i == 0 ? -GP : GP);
                points[k, 1] = (j == 0 ? -GP : GP);
                weights[k] = 1.0;
                k++;
            }
        }

        protected override void ShapeAndDerivs(double xi, double eta, out double[] N, out double[,] dN)
        {
            double[] sx = new double[4] { -1.0, 1.0, 1.0, -1.0 };
            double[] sy = new double[4] { -1.0, -1.0, 1.0, 1.0 };
            N = new double[4];
            dN = new double[4, 2];
            for (int i = 0; i < 4; i++)
            {
                N[i] = 0.25 * (1.0 + sx[i] * xi) * (1.0 + sy[i] * eta);
                dN[i, 0] = 0.25 * sx[i] * (1.0 + sy[i] * eta);
                dN[i, 1] = 0.25 * sy[i] * (1.0 + sx[i] * xi);
            }
        }
    }
}
