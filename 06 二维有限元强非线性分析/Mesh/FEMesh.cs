using System;
using System.Collections.Generic;
using NonlinearFEM2D.Elements;
using NonlinearFEM2D.Core;

namespace NonlinearFEM2D.Mesh
{
    public class BoundaryCondition
    {
        public int NodeId;
        public int Dof;       // 0=u 1=v
        public double Value;
    }

    public sealed class FEMesh
    {
        public int NumNodes;
        public double[] X, Y;
        public List<FiniteElement> Elements = new List<FiniteElement>();
        public List<BoundaryCondition> Dirichlet = new List<BoundaryCondition>();
        // 节点集中力
        public double[] NodeForceX, NodeForceY;
        // 计算结果
        public double[] DisplacementU, DisplacementV;
        public double[] VelocityU, VelocityV;
        public double[] Pressure;
        public double[] StressXX, StressYY, StressXY;
        public double[] VonMises;
        public double[] EqPlasticStrain;
        public ElementType ElementType;
        // 模态
        public List<double> ModalFrequencies = new List<double>();
        public List<DenseVector> ModalShapes = new List<DenseVector>();

        public void Allocate(int numNodes)
        {
            this.NumNodes = numNodes;
            this.X = new double[numNodes]; this.Y = new double[numNodes];
            this.NodeForceX = new double[numNodes]; this.NodeForceY = new double[numNodes];
            this.DisplacementU = new double[numNodes]; this.DisplacementV = new double[numNodes];
            this.VelocityU = new double[numNodes]; this.VelocityV = new double[numNodes];
            this.Pressure = new double[numNodes];
            this.StressXX = new double[numNodes]; this.StressYY = new double[numNodes]; this.StressXY = new double[numNodes];
            this.VonMises = new double[numNodes];
            this.EqPlasticStrain = new double[numNodes];
        }

        public void ResetFields()
        {
            Array.Clear(this.DisplacementU, 0, this.NumNodes);
            Array.Clear(this.DisplacementV, 0, this.NumNodes);
            Array.Clear(this.VelocityU, 0, this.NumNodes);
            Array.Clear(this.VelocityV, 0, this.NumNodes);
            Array.Clear(this.Pressure, 0, this.NumNodes);
            Array.Clear(this.StressXX, 0, this.NumNodes);
            Array.Clear(this.StressYY, 0, this.NumNodes);
            Array.Clear(this.StressXY, 0, this.NumNodes);
            Array.Clear(this.VonMises, 0, this.NumNodes);
            Array.Clear(this.EqPlasticStrain, 0, this.NumNodes);
            this.ModalFrequencies.Clear(); this.ModalShapes.Clear();
        }
    }
}
