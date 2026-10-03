namespace ModalFEM2D.Mesh
{
    using System.Collections.Generic;
    using ModalFEM2D.Elements;

    public enum BCKind { FixedX, FixedY, Fixed, ForceX, ForceY, VelInX, VelInY, Pressure }

    public struct BCEntry
    {
        public int Node;
        public BCKind Kind;
        public double Value;
        public BCEntry(int n, BCKind k, double v) { Node = n; Kind = k; Value = v; }
    }

    public sealed class FEMesh
    {
        public ElementType ElementType;
        public double YoungModulus;
        public double PoissonRatio;
        public double Density;
        public double Viscosity;
        public double Thickness;

        public List<double[]> Nodes = new List<double[]>(); // [2] x,y
        public List<int[]> Elements = new List<int[]>();    // [npe]
        public List<BCEntry> BCs = new List<BCEntry>();

        // Fast-access arrays (built by BuildConnectivity)
        public int NumNodes;
        public int NumElements;
        public int NodesPerElement;
        public double[] X;
        public double[] Y;
        public int[,] Connectivity;
        public int[] DirichletU;  // -1=free, else prescribed value id via DirichletUValue
        public int[] DirichletV;
        public double[] DirichletUValue;
        public double[] DirichletVValue;
        public double[] ForceU;
        public double[] ForceV;

        // Results
        public double[] DisplacementU;
        public double[] DisplacementV;
        public double[] VelocityU;
        public double[] VelocityV;
        public double[] Pressure;
        public double[] StressXX;
        public double[] StressYY;
        public double[] StressXY;
        public double[] VonMises;

        // Modal
        public double[] ModalFreqHz;
        public double[,] ModalShapes; // [mode, 2*n]

        public void AddNode(double x, double y) { Nodes.Add(new double[2] { x, y }); }
        public void AddElement(int[] nodes) { Elements.Add(nodes); }
        public void AddBC(int node, BCKind k, double v) { BCs.Add(new BCEntry(node, k, v)); }

        public void BuildConnectivity()
        {
            NumNodes = Nodes.Count;
            NumElements = Elements.Count;
            if (NumElements == 0) { NodesPerElement = 0; return; }
            NodesPerElement = Elements[0].Length;
            X = new double[NumNodes]; Y = new double[NumNodes];
            for (int i = 0; i < NumNodes; i++) { X[i] = Nodes[i][0]; Y[i] = Nodes[i][1]; }
            Connectivity = new int[NumElements, NodesPerElement];
            for (int e = 0; e < NumElements; e++)
                for (int k = 0; k < NodesPerElement; k++)
                    Connectivity[e, k] = Elements[e][k];
            DirichletU = new int[NumNodes]; DirichletV = new int[NumNodes];
            DirichletUValue = new double[NumNodes]; DirichletVValue = new double[NumNodes];
            ForceU = new double[NumNodes]; ForceV = new double[NumNodes];
            for (int i = 0; i < NumNodes; i++) { DirichletU[i] = -1; DirichletV[i] = -1; }
            foreach (BCEntry bc in BCs)
            {
                int n = bc.Node; if (n < 0 || n >= NumNodes) continue;
                switch (bc.Kind)
                {
                    case BCKind.FixedX: DirichletU[n] = 1; DirichletUValue[n] = bc.Value; break;
                    case BCKind.FixedY: DirichletV[n] = 1; DirichletVValue[n] = bc.Value; break;
                    case BCKind.Fixed:  DirichletU[n] = 1; DirichletV[n] = 1;
                                        DirichletUValue[n] = bc.Value; DirichletVValue[n] = bc.Value; break;
                    case BCKind.ForceX: ForceU[n] += bc.Value; break;
                    case BCKind.ForceY: ForceV[n] += bc.Value; break;
                    case BCKind.VelInX: DirichletU[n] = 1; DirichletUValue[n] = bc.Value; break;
                    case BCKind.VelInY: DirichletV[n] = 1; DirichletVValue[n] = bc.Value; break;
                    case BCKind.Pressure: break; // applied in solver
                }
            }
            DisplacementU = new double[NumNodes]; DisplacementV = new double[NumNodes];
            VelocityU = new double[NumNodes]; VelocityV = new double[NumNodes];
            Pressure = new double[NumNodes];
            StressXX = new double[NumNodes]; StressYY = new double[NumNodes];
            StressXY = new double[NumNodes]; VonMises = new double[NumNodes];
        }

        public FiniteElement CreateElement(int eId)
        {
            FiniteElement el = ElementFactory.Create(ElementType);
            el.Id = eId;
            int npe = NodesPerElement;
            el.NodeIds = new int[npe];
            el.Coords = new double[npe, 2];
            for (int k = 0; k < npe; k++)
            {
                int n = Connectivity[eId, k];
                el.NodeIds[k] = n;
                el.Coords[k, 0] = X[n]; el.Coords[k, 1] = Y[n];
            }
            el.YoungModulus = YoungModulus;
            el.PoissonRatio = PoissonRatio;
            el.Density = Density;
            el.Viscosity = Viscosity;
            el.Thickness = Thickness;
            return el;
        }

        public void ResetResults()
        {
            for (int i = 0; i < NumNodes; i++)
            {
                DisplacementU[i] = 0.0; DisplacementV[i] = 0.0;
                VelocityU[i] = 0.0; VelocityV[i] = 0.0;
                Pressure[i] = 0.0;
                StressXX[i] = 0.0; StressYY[i] = 0.0; StressXY[i] = 0.0; VonMises[i] = 0.0;
            }
            ModalFreqHz = null; ModalShapes = null;
        }
    }
}
