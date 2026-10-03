namespace ModalFEM2D.Elements
{
    /// <summary>
    /// 3-node constant-strain triangle (linear triangle, CST). Area-coordinate integration (1 point, exact for linear).
    /// Natural coordinates: L1=1-xi-eta, L2=xi, L3=eta with xi>=0, eta>=0, xi+eta<=1.
    /// </summary>
    public sealed class CST3Element : FiniteElement
    {
        public override ElementType Type { get { return ElementType.CST3; } }
        public override int NodesPerElement { get { return 3; } }
        public override int PressureNodes { get { return 3; } }

        protected override void GetIntegrationPoints(out double[,] points, out double[] weights)
        {
            points = new double[1, 2]; points[0, 0] = 1.0 / 3.0; points[0, 1] = 1.0 / 3.0;
            weights = new double[1] { 0.5 };
        }

        protected override void ShapeAndDerivs(double xi, double eta, out double[] N, out double[,] dN)
        {
            N = new double[3] { 1.0 - xi - eta, xi, eta };
            dN = new double[3, 2] { { -1.0, -1.0 }, { 1.0, 0.0 }, { 0.0, 1.0 } };
        }
    }
}
