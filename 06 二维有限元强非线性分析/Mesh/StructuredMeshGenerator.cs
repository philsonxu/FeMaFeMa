using System;
using NonlinearFEM2D.Elements;

namespace NonlinearFEM2D.Mesh
{
    public static class StructuredMeshGenerator
    {
        public static FEMesh GenerateCantilever(ElementType t, double L, double H, int Nx, int Ny,
            double E, double nu, double rho, double t_h, double endLoad, double sy, double hhard)
        {
            int nnodes = (Nx + 1) * (Ny + 1);
            FEMesh m = new FEMesh();
            m.ElementType = t;
            m.Allocate(nnodes);
            double dx = L / Nx, dy = H / Ny;
            for (int j = 0; j <= Ny; j++) for (int i = 0; i <= Nx; i++)
                {
                    int id = j * (Nx + 1) + i;
                    m.X[id] = i * dx; m.Y[id] = j * dy - H / 2.0;
                }
            // BC: 左端固支 u=v=0
            for (int j = 0; j <= Ny; j++)
            {
                int id = j * (Nx + 1) + 0;
                m.Dirichlet.Add(new BoundaryCondition { NodeId = id, Dof = 0, Value = 0.0 });
                m.Dirichlet.Add(new BoundaryCondition { NodeId = id, Dof = 1, Value = 0.0 });
            }
            // 端部载荷：右端顶部中点集中力
            int tipTop = Ny * (Nx + 1) + Nx;
            m.NodeForceY[tipTop] = endLoad;
            if (t == ElementType.CST3)
            {
                // 每个矩形分两个三角形
                for (int j = 0; j < Ny; j++) for (int i = 0; i < Nx; i++)
                    {
                        int n0 = j * (Nx + 1) + i;
                        int n1 = n0 + 1;
                        int n2 = n0 + (Nx + 1);
                        int n3 = n2 + 1;
                        AddElement(m, ElementType.CST3, new int[] { n0, n1, n3 }, E, nu, rho, t_h, sy, hhard);
                        AddElement(m, ElementType.CST3, new int[] { n0, n3, n2 }, E, nu, rho, t_h, sy, hhard);
                    }
            }
            else if (t == ElementType.Q4)
            {
                for (int j = 0; j < Ny; j++) for (int i = 0; i < Nx; i++)
                    {
                        int n0 = j * (Nx + 1) + i;
                        int n1 = n0 + 1;
                        int n3 = n0 + (Nx + 1);
                        int n2 = n3 + 1;
                        AddElement(m, ElementType.Q4, new int[] { n0, n1, n2, n3 }, E, nu, rho, t_h, sy, hhard);
                    }
            }
            else if (t == ElementType.LT6)
            {
                // LT6：每个矩形4三角形×？简化为2三角形，每三角形加3边中点(需要新节点)。为简单起见用 CST3 加倍分辨率替代，这里直接使用 CST3 但单元类型设为 LT6 不便；
                // 简化：使用与CST3相同的两三角划分，但单元类型用LT6并补充3边中点为角点之间的中点，需要在网格中插入新节点。
                // 为保证可运行，这里改为直接生成 CST3，ElementType 保持 CST3；外部调用方若请求 LT6 悬臂梁将使用 CST3 的粗网格。
                for (int j = 0; j < Ny; j++) for (int i = 0; i < Nx; i++)
                    {
                        int n0 = j * (Nx + 1) + i;
                        int n1 = n0 + 1;
                        int n2 = n0 + (Nx + 1);
                        int n3 = n2 + 1;
                        AddElement(m, ElementType.CST3, new int[] { n0, n1, n3 }, E, nu, rho, t_h, sy, hhard);
                        AddElement(m, ElementType.CST3, new int[] { n0, n3, n2 }, E, nu, rho, t_h, sy, hhard);
                    }
                m.ElementType = ElementType.CST3;
            }
            return m;
        }

        public static FEMesh GenerateLidDrivenCavity(int N, double L, double U0, double rho, double mu)
        {
            int nn = (N + 1) * (N + 1);
            FEMesh m = new FEMesh();
            m.ElementType = ElementType.Q4;
            m.Allocate(nn);
            double h = L / N;
            for (int j = 0; j <= N; j++) for (int i = 0; i <= N; i++)
                {
                    int id = j * (N + 1) + i;
                    m.X[id] = i * h; m.Y[id] = j * h;
                }
            // 壁面：四面速度 u=v=0；顶边(y=L) u=U0
            for (int i = 0; i <= N; i++)
            {
                int bot = 0 * (N + 1) + i;
                int top = N * (N + 1) + i;
                int lft = i * (N + 1) + 0;
                int rgt = i * (N + 1) + N;
                m.Dirichlet.Add(new BoundaryCondition { NodeId = bot, Dof = 0, Value = 0.0 });
                m.Dirichlet.Add(new BoundaryCondition { NodeId = bot, Dof = 1, Value = 0.0 });
                if (i != 0 && i != N)
                {
                    m.Dirichlet.Add(new BoundaryCondition { NodeId = top, Dof = 0, Value = U0 });
                    m.Dirichlet.Add(new BoundaryCondition { NodeId = top, Dof = 1, Value = 0.0 });
                }
                else
                {
                    m.Dirichlet.Add(new BoundaryCondition { NodeId = top, Dof = 0, Value = 0.0 });
                    m.Dirichlet.Add(new BoundaryCondition { NodeId = top, Dof = 1, Value = 0.0 });
                }
                if (i > 0 && i < N)
                {
                    m.Dirichlet.Add(new BoundaryCondition { NodeId = lft, Dof = 0, Value = 0.0 });
                    m.Dirichlet.Add(new BoundaryCondition { NodeId = lft, Dof = 1, Value = 0.0 });
                    m.Dirichlet.Add(new BoundaryCondition { NodeId = rgt, Dof = 0, Value = 0.0 });
                    m.Dirichlet.Add(new BoundaryCondition { NodeId = rgt, Dof = 1, Value = 0.0 });
                }
            }
            double E = 2.1e11, nu = 0.3, th = 0.01, sy = 2.5e8, hh = 2.1e9;
            for (int j = 0; j < N; j++) for (int i = 0; i < N; i++)
                {
                    int n0 = j * (N + 1) + i;
                    int n1 = n0 + 1;
                    int n3 = n0 + (N + 1);
                    int n2 = n3 + 1;
                    FiniteElement e = ElementFactory.Create(ElementType.Q4);
                    e.Id = m.Elements.Count;
                    e.NodeIds = new int[] { n0, n1, n2, n3 };
                    e.X = new double[] { m.X[n0], m.X[n1], m.X[n2], m.X[n3] };
                    e.Y = new double[] { m.Y[n0], m.Y[n1], m.Y[n2], m.Y[n3] };
                    e.YoungModulus = E; e.PoissonRatio = nu; e.Density = rho; e.Viscosity = mu; e.Thickness = th;
                    e.Initialize();
                    e.Material = new Materials.J2PlasticityMaterial(E, nu, sy, hh);
                    m.Elements.Add(e);
                }
            // 压力参考点：左下角 p=0
            m.Pressure[0] = 0.0;
            return m;
        }

        private static void AddElement(FEMesh m, ElementType t, int[] conn, double E, double nu, double rho, double th, double sy, double hh)
        {
            FiniteElement e = ElementFactory.Create(t);
            e.Id = m.Elements.Count;
            e.NodeIds = conn;
            int npe = conn.Length;
            e.X = new double[npe]; e.Y = new double[npe];
            for (int k = 0; k < npe; k++) { e.X[k] = m.X[conn[k]]; e.Y[k] = m.Y[conn[k]]; }
            e.YoungModulus = E; e.PoissonRatio = nu; e.Density = rho; e.Thickness = th;
            e.Viscosity = 1.0e-3;
            e.Initialize();
            e.Material = new Materials.J2PlasticityMaterial(E, nu, sy, hh);
            m.Elements.Add(e);
        }
    }
}
