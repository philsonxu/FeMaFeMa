using System;
using System.Collections.Generic;
using System.IO;
using System.Globalization;
using Fem2DFluid.Models;

namespace Fem2DFluid.Services
{
    /// <summary>
    /// 内置示例生成：
    /// 1. 方腔绕流/势流：1x1方腔，左入口φ=0，右出口φ=1，上下壁面齐次Neumann(固壁不可穿透)
    /// 2. 带障碍管道：2x1矩形中间一个圆形障碍物，入口φ=1出口φ=0，上下壁面+圆柱固壁φ=0? 不对：固壁用Neumann q=0
    /// 3. 达西渗流坝体：梯形坝，上游水头H1=10，下游H2=2
    /// </summary>
    public class SampleBuilder
    {
        /// <summary>
        /// 示例1：二维平面势流——方腔均匀流
        /// 1×1正方形区域，左侧φ=0，右侧φ=1，上下壁面∂φ/∂n=0
        /// 解析解：φ=x，速度u=(1,0)
        /// </summary>
        public static FemModel CavityFlow(int nx, int ny)
        {
            FemModel m = new FemModel();
            m.Title = "方腔均匀势流(解析φ=x)";
            m.ProblemType = 0;
            m.Materials.Add(new Material(1, "Fluid", 1.0));
            int id = 1;
            double Lx = 1.0, Ly = 1.0;
            double dx = Lx / (nx - 1), dy = Ly / (ny - 1);
            Node[,] grid = new Node[nx, ny];
            for (int j = 0; j < ny; j++)
            {
                for (int i = 0; i < nx; i++)
                {
                    double x = i * dx, y = j * dy;
                    Node n = new Node(id, x, y);
                    // 左边界φ=0，右边界φ=1（Dirichlet）
                    if (i == 0) { n.BCType = 1; n.BCValue = 0.0; }
                    else if (i == nx - 1) { n.BCType = 1; n.BCValue = 1.0; }
                    // 上下壁面Neumann q=0（自然边界条件，默认即为0，无需额外设置）
                    m.Nodes.Add(n);
                    grid[i, j] = n;
                    id++;
                }
            }
            int eid = 1;
            for (int j = 0; j < ny - 1; j++)
            {
                for (int i = 0; i < nx - 1; i++)
                {
                    Node a = grid[i, j];
                    Node b = grid[i + 1, j];
                    Node c = grid[i + 1, j + 1];
                    Node d = grid[i, j + 1];
                    // 两个三角形：a-b-c 与 a-c-d（保持逆时针为正）
                    m.Elements.Add(new Element(eid++, a.Id, b.Id, c.Id, 1));
                    m.Elements.Add(new Element(eid++, a.Id, c.Id, d.Id, 1));
                }
            }
            // 计算单元几何
            for (int i = 0; i < m.Elements.Count; i++) Triangulator.ComputeElementGeom(m.Elements[i], m);
            return m;
        }

        /// <summary>
        /// 示例2：绕圆柱势流（简单矩形区域内挖圆孔，带障碍物的管道流）
        /// 区域 [0,4]x[0,1]，圆柱圆心(2,0.5)半径0.2
        /// 入口(左)φ=0，出口(右)φ=4（产生u≈1来流），上下壁面+圆柱面Neumann固壁
        /// 使用DXF文件导入后自动剖分；也可直接程序生成粗糙网格作为演示
        /// </summary>
        public static void WriteFlowAroundCylinderDxf(string path)
        {
            // 直接写出DXF文件：矩形外框+圆（使用netDxf）
            netDxf.DxfDocument doc = new netDxf.DxfDocument();
            // 外框（逆时针）：(0,0)->(4,0)->(4,1)->(0,1)->(0,0)
            doc.Entities.Add(new netDxf.Entities.Line(new netDxf.Vector3(0, 0, 0), new netDxf.Vector3(4, 0, 0)));
            doc.Entities.Add(new netDxf.Entities.Line(new netDxf.Vector3(4, 0, 0), new netDxf.Vector3(4, 1, 0)));
            doc.Entities.Add(new netDxf.Entities.Line(new netDxf.Vector3(4, 1, 0), new netDxf.Vector3(0, 1, 0)));
            doc.Entities.Add(new netDxf.Entities.Line(new netDxf.Vector3(0, 1, 0), new netDxf.Vector3(0, 0, 0)));
            // 圆柱障碍（顺时针作为孔）
            doc.Entities.Add(new netDxf.Entities.Circle(new netDxf.Vector3(2, 0.5, 0), 0.2));
            doc.Save(path);
        }

        /// <summary>
        /// 给圆柱绕流几何设置边界条件：左入口φ=0，右出口φ=4，上下壁面/圆柱面Neumann=0（不指定BC即自然）
        /// 注意：在后续剖分添加内部点时，只对落在入口/出口直线上的边界节点重新标记为Dirichlet，
        /// 其他边界节点（上下壁、圆柱）保持BCType=0 作为Neumann自然边界。
        /// 为区分"入口/出口"和"壁面"，我们用 BCType=1 且 BCTag 用 BCValue 的符号来标记：
        ///   BCValue=0.0 且 X≈xmin  => 左入口（Dirichlet φ=0）
        ///   BCValue=4.0 且 X≈xmax  => 右出口（Dirichlet φ=4）
        ///   其他边界节点BCType=0  => Neumann壁面
        /// 为避免DoMesh中逻辑错误，在节点上打Tag：把BCValue的符号标记直接通过BCType区分：
        ///   BCType=3  => 入口/出口（剖分后再转为BCType=1，保留BCValue）
        ///   BCType=4  => 壁面/圆柱（剖分后转为BCType=0 Neumann）
        /// </summary>
        public static void ApplyCylinderBC(FemModel model)
        {
            double xmin, ymin, xmax, ymax;
            model.ComputeBounds(out xmin, out ymin, out xmax, out ymax);
            double cx = (xmin + xmax) / 2.0;
            double cy = (ymin + ymax) / 2.0;
            // 用BCType=100/101/102打Tag（剖分后再重置为正式类型）
            const int TAG_INLET = 100;
            const int TAG_OUTLET = 101;
            const int TAG_WALL = 102;
            for (int i = 0; i < model.Nodes.Count; i++)
            {
                Node n = model.Nodes[i];
                if (Math.Abs(n.X - xmin) < 1e-6)
                {
                    n.BCType = TAG_INLET; n.BCValue = 0.0;
                }
                else if (Math.Abs(n.X - xmax) < 1e-6)
                {
                    n.BCType = TAG_OUTLET; n.BCValue = xmax - xmin;
                }
                else if (Math.Abs(n.Y - ymin) < 1e-6 || Math.Abs(n.Y - ymax) < 1e-6)
                {
                    n.BCType = TAG_WALL; n.BCValue = 0.0;
                }
                else
                {
                    double dx = n.X - cx;
                    double dy = n.Y - cy;
                    double dist = Math.Sqrt(dx * dx + dy * dy);
                    // 圆柱壁面识别：距离圆柱半径0.2的容差取h/3量级，h≈0.2时容差≈0.07
                    // 为稳妥取0.05（h=0.2时边界节点步长0.2，内部点最少离边界0.1远，故0.05不会误判）
                    if (Math.Abs(dist - 0.2) < 0.05) { n.BCType = TAG_WALL; n.BCValue = 0.0; }
                    else { n.BCType = 0; n.BCValue = 0.0; }
                }
            }
        }

        /// <summary>
        /// 示例3：达西渗流——梯形坝体（文本格式网格）
        /// 生成结构化网格的文本文件
        /// </summary>
        public static void WriteDamSeepageInput(string path)
        {
            using (StreamWriter sw = new StreamWriter(path, false, System.Text.Encoding.ASCII))
            {
                sw.WriteLine("# 梯形坝体达西渗流示例：上游φ=10，下游φ=2，坝底不透水Neumann=0");
                sw.WriteLine("TITLE DamSeepage");
                sw.WriteLine("PROBLEM 1");
                sw.WriteLine("MATERIALS 1");
                sw.WriteLine("1 Soil 0.001 0");
                // 构造梯形：上游面(0,0)-(0,10)，下游面(20,0)-(20,4)，中间坡面
                int nx = 21, ny = 11;
                // 节点
                List<Node> nodes = new List<Node>();
                int id = 1;
                for (int j = 0; j < ny; j++)
                {
                    double y = 10.0 * j / (ny - 1);
                    // 上游水平宽度20，下游顶宽因边坡而变：这里简化为矩形20×10
                    // （真正梯形剖分需要更复杂，矩形已足够演示渗流）
                    for (int i = 0; i < nx; i++)
                    {
                        double x = 20.0 * i / (nx - 1);
                        Node n = new Node(id, x, y);
                        // 上游水位(左) φ=10，下游水位(右) φ=2，底部y=0 Neumann=0
                        if (i == 0 && y <= 10) { n.BCType = 1; n.BCValue = 10; }
                        else if (i == nx - 1 && y <= 4) { n.BCType = 1; n.BCValue = 2; }
                        // 顶面y=10为渗流自由面，简化为Dirichlet φ=y(即10)
                        else if (j == ny - 1) { n.BCType = 1; n.BCValue = y; }
                        nodes.Add(n);
                        id++;
                    }
                }
                sw.WriteLine("NODES " + nodes.Count);
                for (int i = 0; i < nodes.Count; i++)
                {
                    Node n = nodes[i];
                    sw.WriteLine(string.Format(CultureInfo.InvariantCulture,
                        "{0} {1:F4} {2:F4} {3} {4:F4}", n.Id, n.X, n.Y, n.BCType, n.BCValue));
                }
                List<Element> elems = new List<Element>();
                int eid = 1;
                for (int j = 0; j < ny - 1; j++)
                {
                    for (int i = 0; i < nx - 1; i++)
                    {
                        int a = j * nx + i + 1;
                        int b = a + 1;
                        int c = a + nx + 1;
                        int d = a + nx;
                        elems.Add(new Element(eid++, a, b, c, 1));
                        elems.Add(new Element(eid++, a, c, d, 1));
                    }
                }
                sw.WriteLine("ELEMENTS " + elems.Count);
                for (int i = 0; i < elems.Count; i++)
                {
                    Element e = elems[i];
                    sw.WriteLine(string.Format("{0} {1} {2} {3} {4}", e.Id, e.N1, e.N2, e.N3, e.MatId));
                }
            }
        }
    }
}
