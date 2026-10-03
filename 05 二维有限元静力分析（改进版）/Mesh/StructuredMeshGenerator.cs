// StructuredMeshGenerator.cs - 矩形区域结构化网格自动剖分（悬臂梁/顶盖驱动方腔/纯矩形）
using System;
using System.Collections.Generic;
using FEM2D.Elements;

namespace FEM2D.Mesh
{
    public static class StructuredMeshGenerator
    {
        private const double DefaultE = 2.1e11;
        private const double DefaultNu = 0.3;
        private const double DefaultRho = 7850.0;
        private const double DefaultMu = 1.0e-3;
        private const double DefaultThick = 1.0;

        /// <summary>悬臂梁：左端固支 u=v=0，右端中点加竖向力 loadY；材料为钢。</summary>
        public static FEMesh GenerateCantilever(ElementType type, double L, double H, int nx, int ny, double loadY,
            double E = DefaultE, double nu = DefaultNu, double rho = DefaultRho, double mu = DefaultMu, double thick = DefaultThick)
        {
            FEMesh mesh = BuildRectMesh(type, nx, ny, L, H, E, nu, rho, mu, thick);
            int NNx = (type == ElementType.LT6) ? (2 * nx + 1) : (nx + 1);
            int NNy = (type == ElementType.LT6) ? (2 * ny + 1) : (ny + 1);
            // 左端固支 (x=0)
            for (int j = 0; j < NNy; j++)
            {
                int n = j * NNx + 0;
                mesh.DirichletBCs.Add(Tuple.Create(n, 0, 0.0));
                mesh.DirichletBCs.Add(Tuple.Create(n, 1, 0.0));
            }
            // 右端中点竖向荷载 (x=L)
            int tipMid = (ny / ((type == ElementType.LT6) ? 1 : 1)) * NNx + (NNx - 1);
            if (type == ElementType.LT6) tipMid = ny * NNx + (NNx - 1);
            else tipMid = (ny / 2) * NNx + (NNx - 1);
            mesh.NodalForces.Add(Tuple.Create(tipMid, 1, loadY));
            mesh.BuildConnectivity();
            return mesh;
        }

        /// <summary>顶盖驱动方腔：四周墙无滑移、顶盖 x 方向速度 U0。为 NS 准备。</summary>
        public static FEMesh GenerateLidDrivenCavity(ElementType type, double L, int nx, double U0,
            double rho = 1.0, double mu = 0.01, double thick = 1.0)
        {
            int ny = nx;
            FEMesh mesh = BuildRectMesh(type, nx, ny, L, L, 0.0, 0.0, rho, mu, thick);
            int NNx = (type == ElementType.LT6) ? (2 * nx + 1) : (nx + 1);
            int NNy = NNx;
            // 底部无滑移
            for (int i = 0; i < NNx; i++)
            {
                int n = 0 * NNx + i;
                mesh.VelocityBCs.Add(Tuple.Create(n, 0, 0.0));
                mesh.VelocityBCs.Add(Tuple.Create(n, 1, 0.0));
            }
            // 顶盖 u=U0
            for (int i = 0; i < NNx; i++)
            {
                int n = (NNy - 1) * NNx + i;
                mesh.VelocityBCs.Add(Tuple.Create(n, 0, U0));
                mesh.VelocityBCs.Add(Tuple.Create(n, 1, 0.0));
            }
            // 左右墙无滑移
            for (int j = 0; j < NNy; j++)
            {
                int nL = j * NNx + 0;
                int nR = j * NNx + (NNx - 1);
                mesh.VelocityBCs.Add(Tuple.Create(nL, 0, 0.0));
                mesh.VelocityBCs.Add(Tuple.Create(nL, 1, 0.0));
                mesh.VelocityBCs.Add(Tuple.Create(nR, 0, 0.0));
                mesh.VelocityBCs.Add(Tuple.Create(nR, 1, 0.0));
            }
            mesh.BuildConnectivity();
            return mesh;
        }

        private static FEMesh BuildRectMesh(ElementType type, int nx, int ny, double L, double H,
            double E, double nu, double rho, double mu, double thick)
        {
            FEMesh mesh = new FEMesh();
            List<FiniteElement> elems = new List<FiniteElement>();
            double[] X, Y;
            int numNodes;

            if (type == ElementType.Q4)
            {
                int NNx = nx + 1, NNy = ny + 1;
                numNodes = NNx * NNy;
                X = new double[numNodes]; Y = new double[numNodes];
                for (int j = 0; j < NNy; j++)
                    for (int i = 0; i < NNx; i++)
                    { int n = j * NNx + i; X[n] = i * L / nx; Y[n] = j * H / ny; }
                for (int j = 0; j < ny; j++)
                    for (int i = 0; i < nx; i++)
                    {
                        int n0 = j * NNx + i, n1 = n0 + 1, n2 = n0 + NNx + 1, n3 = n0 + NNx;
                        Q4Element e = new Q4Element();
                        e.Id = elems.Count; e.NodeIds = new int[] { n0, n1, n2, n3 };
                        e.YoungModulus = E; e.PoissonRatio = nu; e.Density = rho; e.Viscosity = mu; e.Thickness = thick;
                        elems.Add(e);
                    }
            }
            else if (type == ElementType.CST3)
            {
                int NNx = nx + 1, NNy = ny + 1;
                numNodes = NNx * NNy;
                X = new double[numNodes]; Y = new double[numNodes];
                for (int j = 0; j < NNy; j++)
                    for (int i = 0; i < NNx; i++)
                    { int n = j * NNx + i; X[n] = i * L / nx; Y[n] = j * H / ny; }
                for (int j = 0; j < ny; j++)
                    for (int i = 0; i < nx; i++)
                    {
                        int n0 = j * NNx + i, n1 = n0 + 1, n2 = n0 + NNx + 1, n3 = n0 + NNx;
                        CST3Element e1 = new CST3Element();
                        e1.Id = elems.Count; e1.NodeIds = new int[] { n0, n1, n2 };
                        e1.YoungModulus = E; e1.PoissonRatio = nu; e1.Density = rho; e1.Viscosity = mu; e1.Thickness = thick;
                        elems.Add(e1);
                        CST3Element e2 = new CST3Element();
                        e2.Id = elems.Count; e2.NodeIds = new int[] { n0, n2, n3 };
                        e2.YoungModulus = E; e2.PoissonRatio = nu; e2.Density = rho; e2.Viscosity = mu; e2.Thickness = thick;
                        elems.Add(e2);
                    }
            }
            else // LT6
            {
                int NNx = 2 * nx + 1, NNy = 2 * ny + 1;
                numNodes = NNx * NNy;
                X = new double[numNodes]; Y = new double[numNodes];
                for (int j = 0; j < NNy; j++)
                    for (int i = 0; i < NNx; i++)
                    { int n = j * NNx + i; X[n] = i * L / (2.0 * nx); Y[n] = j * H / (2.0 * ny); }
                for (int j = 0; j < ny; j++)
                    for (int i = 0; i < nx; i++)
                    {
                        int i0 = 2 * i, j0 = 2 * j;
                        int n00 = j0 * NNx + i0;
                        int n20 = j0 * NNx + i0 + 2;
                        int n02 = (j0 + 2) * NNx + i0;
                        int n10 = j0 * NNx + i0 + 1;
                        int n11 = (j0 + 1) * NNx + i0 + 1;
                        int n01 = (j0 + 1) * NNx + i0;
                        LT6Element t1 = new LT6Element();
                        t1.Id = elems.Count; t1.NodeIds = new int[] { n00, n20, n02, n10, n11, n01 };
                        t1.YoungModulus = E; t1.PoissonRatio = nu; t1.Density = rho; t1.Viscosity = mu; t1.Thickness = thick;
                        elems.Add(t1);
                        int n22 = (j0 + 2) * NNx + i0 + 2;
                        int n21 = (j0 + 1) * NNx + i0 + 2;
                        int n12 = (j0 + 2) * NNx + i0 + 1;
                        LT6Element t2 = new LT6Element();
                        t2.Id = elems.Count; t2.NodeIds = new int[] { n20, n22, n02, n21, n12, n11 };
                        t2.YoungModulus = E; t2.PoissonRatio = nu; t2.Density = rho; t2.Viscosity = mu; t2.Thickness = thick;
                        elems.Add(t2);
                    }
            }

            mesh.X = X; mesh.Y = Y; mesh.NumNodes = numNodes;
            mesh.Elements = elems; mesh.NumElements = elems.Count;
            mesh.ElementType = type;
            mesh.AssignNodeCoordinatesToElements();
            mesh.UpdateBounds();
            return mesh;
        }
    }
}
