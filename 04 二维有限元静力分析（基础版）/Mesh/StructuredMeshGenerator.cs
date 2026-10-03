using System;
using System.Collections.Generic;

namespace FEM2D.Mesh
{
    /// <summary>
    /// 矩形区域结构化网格自动剖分器。
    /// 几何：[0, Lx] x [0, Ly]，可生成 CST3、LT6、Q4 三种网格。
    /// 自动施加悬臂梁/固定-压力等常见边界供算例使用。
    /// </summary>
    public static class StructuredMeshGenerator
    {
        public static FEMesh GenerateRectangle(double Lx, double Ly, int nx, int ny,
            ElementType type, bool fixLeft, double topPressure, double endLoad = 0.0)
        {
            FEMesh mesh = new FEMesh();
            if (type == ElementType.Q4) return GenerateQ4(mesh, Lx, Ly, nx, ny, fixLeft, topPressure, endLoad);
            if (type == ElementType.CST3) return GenerateCST(mesh, Lx, Ly, nx, ny, fixLeft, topPressure, endLoad, false);
            if (type == ElementType.LT6) return GenerateLT6(mesh, Lx, Ly, nx, ny, fixLeft, topPressure, endLoad);
            return mesh;
        }

        private static FEMesh GenerateQ4(FEMesh mesh, double Lx, double Ly, int nx, int ny,
            bool fixLeft, double topPressure, double endLoad)
        {
            double dx = Lx / nx;
            double dy = Ly / ny;
            int[,] nodeId = new int[nx + 1, ny + 1];
            for (int j = 0; j <= ny; j++)
            {
                for (int i = 0; i <= nx; i++)
                {
                    nodeId[i, j] = mesh.Nodes.Count;
                    mesh.Nodes.Add(new double[] { i * dx, j * dy });
                }
            }
            for (int j = 0; j < ny; j++)
            {
                for (int i = 0; i < nx; i++)
                {
                    int[] conn = new int[]
                    {
                        nodeId[i, j], nodeId[i+1, j], nodeId[i+1, j+1], nodeId[i, j+1]
                    };
                    mesh.Elements.Add(conn);
                    mesh.ElementTypes.Add(ElementType.Q4);
                    mesh.MaterialIDs.Add(0);
                }
            }
            ApplyBoundary(mesh, nx, ny, nodeId, fixLeft, topPressure, endLoad, dx);
            return mesh;
        }

        private static FEMesh GenerateCST(FEMesh mesh, double Lx, double Ly, int nx, int ny,
            bool fixLeft, double topPressure, double endLoad, bool l6)
        {
            double dx = Lx / nx;
            double dy = Ly / ny;
            int[,] nodeId = new int[nx + 1, ny + 1];
            for (int j = 0; j <= ny; j++)
            {
                for (int i = 0; i <= nx; i++)
                {
                    nodeId[i, j] = mesh.Nodes.Count;
                    mesh.Nodes.Add(new double[] { i * dx, j * dy });
                }
            }
            for (int j = 0; j < ny; j++)
            {
                for (int i = 0; i < nx; i++)
                {
                    int n0 = nodeId[i, j];
                    int n1 = nodeId[i + 1, j];
                    int n2 = nodeId[i + 1, j + 1];
                    int n3 = nodeId[i, j + 1];
                    mesh.Elements.Add(new int[] { n0, n1, n2 });
                    mesh.ElementTypes.Add(ElementType.CST3);
                    mesh.MaterialIDs.Add(0);
                    mesh.Elements.Add(new int[] { n0, n2, n3 });
                    mesh.ElementTypes.Add(ElementType.CST3);
                    mesh.MaterialIDs.Add(0);
                }
            }
            ApplyBoundary(mesh, nx, ny, nodeId, fixLeft, topPressure, endLoad, dx);
            return mesh;
        }

        private static FEMesh GenerateLT6(FEMesh mesh, double Lx, double Ly, int nx, int ny,
            bool fixLeft, double topPressure, double endLoad)
        {
            // 二次三角单元：顶点+边中点；先建立顶点网格，再在每个矩形单元细分两个三角单元并插入中点
            double dx = Lx / nx;
            double dy = Ly / ny;
            int[,] cornerId = new int[nx + 1, ny + 1];
            for (int j = 0; j <= ny; j++)
            {
                for (int i = 0; i <= nx; i++)
                {
                    cornerId[i, j] = mesh.Nodes.Count;
                    mesh.Nodes.Add(new double[] { i * dx, j * dy });
                }
            }
            // 边中点缓存：水平边 (i,j)-(i+1,j): key=(i,j,0)；
            // 竖直边 (i,j)-(i,j+1): key=(i,j,1)；
            // 对角线 (i,j)-(i+1,j+1): key=(i,j,2)
            Dictionary<(int, int, int), int> mid = new Dictionary<(int, int, int), int>();

            int MidNode(int a, int b)
            {
                double mx = 0.5 * (mesh.X(a) + mesh.X(b));
                double my = 0.5 * (mesh.Y(a) + mesh.Y(b));
                int id = mesh.Nodes.Count;
                mesh.Nodes.Add(new double[] { mx, my });
                return id;
            }
            int HMid(int i, int j)
            {
                ValueTuple<int, int, int> k = (i, j, 0);
                if (!mid.ContainsKey(k)) mid[k] = MidNode(cornerId[i, j], cornerId[i + 1, j]);
                return mid[k];
            }
            int VMid(int i, int j)
            {
                ValueTuple<int, int, int> k = (i, j, 1);
                if (!mid.ContainsKey(k)) mid[k] = MidNode(cornerId[i, j], cornerId[i, j + 1]);
                return mid[k];
            }
            int DMid(int i, int j)
            {
                ValueTuple<int, int, int> k = (i, j, 2);
                if (!mid.ContainsKey(k)) mid[k] = MidNode(cornerId[i, j], cornerId[i + 1, j + 1]);
                return mid[k];
            }

            for (int j = 0; j < ny; j++)
            {
                for (int i = 0; i < nx; i++)
                {
                    int n0 = cornerId[i, j];
                    int n1 = cornerId[i + 1, j];
                    int n2 = cornerId[i + 1, j + 1];
                    int n3 = cornerId[i, j + 1];
                    int m01 = HMid(i, j);
                    int m12 = VMid(i + 1, j);
                    int m02 = DMid(i, j);
                    int m23 = HMid(i, j + 1);
                    int m03 = VMid(i, j);
                    // 三角 1: n0,n1,n2 中点：m01,m12,m02
                    mesh.Elements.Add(new int[] { n0, n1, n2, m01, m12, m02 });
                    mesh.ElementTypes.Add(ElementType.LT6);
                    mesh.MaterialIDs.Add(0);
                    // 三角 2: n0,n2,n3 中点：m02,m23,m03
                    mesh.Elements.Add(new int[] { n0, n2, n3, m02, m23, m03 });
                    mesh.ElementTypes.Add(ElementType.LT6);
                    mesh.MaterialIDs.Add(0);
                }
            }
            ApplyBoundary(mesh, nx, ny, cornerId, fixLeft, topPressure, endLoad, dx);
            // 对新增中点节点：若两端顶点皆被约束则中点同样约束
            ApplyEdgeMidBoundary(mesh, cornerId, nx, ny, mid);
            return mesh;
        }

        private static void ApplyBoundary(FEMesh mesh, int nx, int ny, int[,] cornerId,
            bool fixLeft, double topPressure, double endLoad, double dx)
        {
            // Dirichlet：左边界固支 u=v=0
            if (fixLeft)
            {
                for (int j = 0; j <= ny; j++)
                {
                    int n = cornerId[0, j];
                    mesh.DirichletBCs.Add((n, 0, 0.0));
                    mesh.DirichletBCs.Add((n, 1, 0.0));
                }
            }
            // Neumann：顶面压力（向下为正）
            if (Math.Abs(topPressure) > 1e-15)
            {
                double faceLen = dx;
                double fy = -topPressure * faceLen; // 面力：y 向；顶面外法线向上，若 topPressure 为压（正）则 fy 负
                for (int i = 0; i < nx; i++)
                {
                    int na = cornerId[i, ny];
                    int nb = cornerId[i + 1, ny];
                    mesh.ConcentratedForces.Add((na, 1, fy * 0.5));
                    mesh.ConcentratedForces.Add((nb, 1, fy * 0.5));
                }
            }
            // 右端点集中力（悬臂梁端部荷载）
            if (Math.Abs(endLoad) > 1e-15)
            {
                int tipTop = cornerId[nx, ny];
                int tipBot = cornerId[nx, 0];
                mesh.ConcentratedForces.Add((tipTop, 1, endLoad));
                mesh.ConcentratedForces.Add((tipBot, 1, endLoad));
            }
        }

        private static void ApplyEdgeMidBoundary(FEMesh mesh, int[,] cornerId, int nx, int ny,
            Dictionary<(int, int, int), int> mid)
        {
            // 左边界中点固支
            for (int j = 0; j < ny; j++)
            {
                ValueTuple<int, int, int> k = (0, j, 1);
                if (mid.ContainsKey(k))
                {
                    int m = mid[k];
                    mesh.DirichletBCs.Add((m, 0, 0.0));
                    mesh.DirichletBCs.Add((m, 1, 0.0));
                }
            }
        }

        /// <summary>
        /// 生成顶盖驱动方腔 NS 算例网格：四周壁面 u=v=0，顶盖 v=0,u=U0；压力参考点。
        /// </summary>
        public static FEMesh GenerateLidDrivenCavity(double L, int n, ElementType type, double U0)
        {
            FEMesh mesh = GenerateRectangle(L, L, n, n, type, false, 0.0, 0.0);
            mesh.Metadata["ProblemType"] = "NS";
            // 清除原有边界
            mesh.DirichletBCs.Clear();
            mesh.ConcentratedForces.Clear();
            // 四周壁面 u=v=0（底面和两侧），顶面 y=L：u=U0, v=0
            int nn = n;
            if (type == ElementType.LT6)
            {
                // 顶点网格还是 (n+1)*(n+1)，但中间还有中点——通过坐标判边界
                for (int k = 0; k < mesh.Nodes.Count; k++)
                {
                    double x = mesh.X(k); double y = mesh.Y(k);
                    bool onBottom = Math.Abs(y) < 1e-9;
                    bool onTop = Math.Abs(y - L) < 1e-9;
                    bool onLeft = Math.Abs(x) < 1e-9;
                    bool onRight = Math.Abs(x - L) < 1e-9;
                    if (onBottom || onLeft || onRight)
                    {
                        mesh.DirichletBCs.Add((k, 0, 0.0));
                        mesh.DirichletBCs.Add((k, 1, 0.0));
                    }
                    if (onTop)
                    {
                        mesh.DirichletBCs.Add((k, 0, U0));
                        mesh.DirichletBCs.Add((k, 1, 0.0));
                    }
                }
            }
            else
            {
                int[,] nid = new int[n + 1, n + 1];
                // 重建顶点索引（Q4/CST 仅顶点）
                if (type == ElementType.Q4)
                {
                    int idx = 0;
                    for (int j = 0; j <= n; j++)
                        for (int i = 0; i <= n; i++)
                            nid[i, j] = idx++;
                }
                else // CST
                {
                    int idx = 0;
                    for (int j = 0; j <= n; j++)
                        for (int i = 0; i <= n; i++)
                            nid[i, j] = idx++;
                }
                for (int i = 0; i <= n; i++)
                {
                    // bottom
                    int nb = nid[i, 0];
                    mesh.DirichletBCs.Add((nb, 0, 0.0));
                    mesh.DirichletBCs.Add((nb, 1, 0.0));
                    // top
                    int nt = nid[i, n];
                    mesh.DirichletBCs.Add((nt, 0, U0));
                    mesh.DirichletBCs.Add((nt, 1, 0.0));
                }
                for (int j = 1; j < n; j++)
                {
                    int nl = nid[0, j];
                    int nr = nid[n, j];
                    mesh.DirichletBCs.Add((nl, 0, 0.0));
                    mesh.DirichletBCs.Add((nl, 1, 0.0));
                    mesh.DirichletBCs.Add((nr, 0, 0.0));
                    mesh.DirichletBCs.Add((nr, 1, 0.0));
                }
            }
            // 压力参考点：左下角 p=0
            mesh.DirichletBCs.Add((0, 2, 0.0));
            return mesh;
        }
    }
}
