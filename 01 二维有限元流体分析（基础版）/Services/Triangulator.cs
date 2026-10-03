using System;
using System.Collections.Generic;
using Fem2DFluid.Models;

namespace Fem2DFluid.Services
{
    /// <summary>
    /// 简单Delaunay三角剖分器（Bowyer-Watson算法），支持域内插入点
    /// 可对给定点集进行三角剖分，然后通过边界裁剪得到域内单元
    /// </summary>
    public class Triangulator
    {
        private class Triangle2D
        {
            public int A, B, C;
            public double Cx, Cy, R2; // 外接圆圆心与半径平方
            public bool Bad;
            public int MatId;

            public Triangle2D(int a, int b, int c, IList<Node> nodes)
            {
                A = a; B = b; C = c; Bad = false; MatId = 1;
                CalcCircumcircle(nodes);
            }

            public bool CircumContains(Node p, IList<Node> nodes)
            {
                double dx = p.X - Cx;
                double dy = p.Y - Cy;
                return dx * dx + dy * dy <= R2 * (1 + 1e-12);
            }

            public bool HasVertex(int v)
            {
                return A == v || B == v || C == v;
            }

            private void CalcCircumcircle(IList<Node> nodes)
            {
                // 注意：Triangle2D存的是节点Id，nodes列表是FemModel.Nodes(可能含超级三角形节点)，
                // 节点Id不一定与列表下标-1一致——通过扫描查找对应节点。
                Node pa = FindNodeById(nodes, A);
                Node pb = FindNodeById(nodes, B);
                Node pc = FindNodeById(nodes, C);
                if (pa == null || pb == null || pc == null) { Cx = 0; Cy = 0; R2 = 1e30; return; }
                double ax = pa.X, ay = pa.Y;
                double bx = pb.X, by = pb.Y;
                double cx = pc.X, cy = pc.Y;
                double d = 2 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by));
                if (Math.Abs(d) < 1e-12) { Cx = (ax + bx + cx) / 3; Cy = (ay + by + cy) / 3; R2 = 1e30; return; }
                double ax2 = ax * ax + ay * ay;
                double bx2 = bx * bx + by * by;
                double cx2 = cx * cx + cy * cy;
                Cx = (ax2 * (by - cy) + bx2 * (cy - ay) + cx2 * (ay - by)) / d;
                Cy = (ax2 * (cx - bx) + bx2 * (ax - cx) + cx2 * (bx - ax)) / d;
                double dx2 = ax - Cx, dy2 = ay - Cy;
                R2 = dx2 * dx2 + dy2 * dy2;
            }

            public static Node FindNodeById(IList<Node> nodes, int id)
            {
                for (int i = 0; i < nodes.Count; i++) if (nodes[i].Id == id) return nodes[i];
                return null;
            }
        }

        /// <summary>
        /// 对点集执行Bowyer-Watson Delaunay三角剖分（含超级三角形）
        /// </summary>
        public static void Delaunay(FemModel model)
        {
            double xmin, ymin, xmax, ymax;
            model.ComputeBounds(out xmin, out ymin, out xmax, out ymax);
            double dx = xmax - xmin, dy = ymax - ymin;
            double dmax = Math.Max(dx, dy) * 20.0;
            double cx = (xmax + xmin) / 2.0;
            double cy = (ymax + ymin) / 2.0;

            // 添加超级三角形三个顶点（用非常大的负数/正数ID，标记为临时节点）
            int superId1 = -1, superId2 = -2, superId3 = -3;
            Node s1 = new Node(superId1, cx - dmax, cy - dmax);
            Node s2 = new Node(superId2, cx + dmax, cy - dmax);
            Node s3 = new Node(superId3, cx, cy + dmax);
            // 临时插入到nodes末尾（记住起始索引以便移除）
            int startSuper = model.Nodes.Count;
            model.Nodes.Add(s1); model.Nodes.Add(s2); model.Nodes.Add(s3);

            List<Triangle2D> triangles = new List<Triangle2D>();
            triangles.Add(new Triangle2D(superId1, superId2, superId3, model.Nodes));

            List<Node> pts = new List<Node>();
            for (int i = 0; i < model.Nodes.Count; i++)
            {
                if (model.Nodes[i].Id > 0) pts.Add(model.Nodes[i]);
            }
            for (int ip = 0; ip < pts.Count; ip++)
            {
                Node p = pts[ip];
                List<Triangle2D> bad = new List<Triangle2D>();
                for (int t = 0; t < triangles.Count; t++)
                {
                    if (triangles[t].CircumContains(p, model.Nodes))
                    {
                        triangles[t].Bad = true;
                        bad.Add(triangles[t]);
                    }
                }
                List<Edge> polygon = new List<Edge>();
                for (int b = 0; b < bad.Count; b++)
                {
                    Triangle2D t = bad[b];
                    Edge[] triEdges = new Edge[] { new Edge(t.A, t.B), new Edge(t.B, t.C), new Edge(t.C, t.A) };
                    for (int e = 0; e < 3; e++)
                    {
                        bool shared = false;
                        for (int o = 0; o < bad.Count; o++)
                        {
                            if (o == b) continue;
                            Triangle2D oth = bad[o];
                            Edge[] oEdges = new Edge[] { new Edge(oth.A, oth.B), new Edge(oth.B, oth.C), new Edge(oth.C, oth.A) };
                            for (int oe = 0; oe < 3; oe++)
                            {
                                if (EdgeEqualUnordered(triEdges[e], oEdges[oe])) { shared = true; break; }
                            }
                            if (shared) break;
                        }
                        if (!shared) polygon.Add(triEdges[e]);
                    }
                }
                for (int t = triangles.Count - 1; t >= 0; t--)
                {
                    if (triangles[t].Bad) triangles.RemoveAt(t);
                }
                for (int e = 0; e < polygon.Count; e++)
                {
                    triangles.Add(new Triangle2D(polygon[e].P1, polygon[e].P2, p.Id, model.Nodes));
                }
            }

            // 移除包含超级三角形顶点的单元
            for (int t = triangles.Count - 1; t >= 0; t--)
            {
                if (triangles[t].HasVertex(superId1) || triangles[t].HasVertex(superId2) || triangles[t].HasVertex(superId3))
                    triangles.RemoveAt(t);
            }
            // 移除超级三角形节点
            model.Nodes.RemoveRange(startSuper, 3);

            model.Elements.Clear();
            int eid = 1;
            for (int t = 0; t < triangles.Count; t++)
            {
                Triangle2D tr = triangles[t];
                Node na = Triangle2D.FindNodeById(model.Nodes, tr.A);
                Node nb = Triangle2D.FindNodeById(model.Nodes, tr.B);
                Node nc = Triangle2D.FindNodeById(model.Nodes, tr.C);
                if (na == null || nb == null || nc == null) continue;
                // 保证节点顺序为 CCW（逆时针），面积为正
                double cross = (nb.X - na.X) * (nc.Y - na.Y) - (nc.X - na.X) * (nb.Y - na.Y);
                int n1 = tr.A, n2 = tr.B, n3 = tr.C;
                if (cross < 0) { int tmp = n2; n2 = n3; n3 = tmp; }
                Element elem = new Element(eid, n1, n2, n3, 1);
                ComputeElementGeom(elem, model);
                if (elem.Area < 1e-14) continue; // 零面积单元丢弃
                model.Elements.Add(elem);
                eid++;
            }
        }

        /// <summary>
        /// 在边界包围盒内生成规则网格点，对DXF导入的几何进行自动剖分
        /// h为目标网格尺寸
        /// </summary>
        public static void GenerateStructuredPoints(FemModel model, double h)
        {
            // 第一步：先对过长的边界边按 h 做等距细分，确保边界节点密度与内部一致
            SubdivideBoundaryEdges(model, h);

            double xmin, ymin, xmax, ymax;
            model.ComputeBounds(out xmin, out ymin, out xmax, out ymax);
            int nx = Math.Max(3, (int)Math.Ceiling((xmax - xmin) / h) + 1);
            int ny = Math.Max(3, (int)Math.Ceiling((ymax - ymin) / h) + 1);
            double dx = (xmax - xmin) / (nx - 1);
            double dy = (ymax - ymin) / (ny - 1);

            // 找到已有最大ID
            int maxId = 0;
            for (int i = 0; i < model.Nodes.Count; i++)
            {
                if (model.Nodes[i].Id > maxId) maxId = model.Nodes[i].Id;
            }
            // 用边界边集合进行"简单内外判定"：只保留在边界多边形内部的点
            List<Edge> boundary = BoundaryEdgeCache.Get(model);
            if (boundary == null || boundary.Count == 0)
            {
                // 没边界，就用所有网格点
                for (int i = 1; i < nx - 1; i++)
                {
                    for (int j = 1; j < ny - 1; j++)
                    {
                        double x = xmin + i * dx;
                        double y = ymin + j * dy;
                        if (PointTooCloseToExisting(model, x, y, h * 0.3)) continue;
                        maxId++;
                        model.Nodes.Add(new Node(maxId, x, y));
                    }
                }
                return;
            }
            // 把边界从点构建多边形列表（可能多连通：外+孔）
            List<List<int>> polygons = BuildPolygons(boundary, model);

            // 添加内部网格点，用射线法判定是否在域内（奇数次在内部）
            for (int i = 0; i <= nx - 1; i++)
            {
                for (int j = 0; j <= ny - 1; j++)
                {
                    double x = xmin + i * dx;
                    double y = ymin + j * dy;
                    if (PointTooCloseToExisting(model, x, y, h * 0.3)) continue;
                    if (PointInsidePolygons(x, y, polygons, model))
                    {
                        maxId++;
                        model.Nodes.Add(new Node(maxId, x, y));
                    }
                }
            }
        }

        private static bool PointTooCloseToExisting(FemModel model, double x, double y, double tol)
        {
            double tol2 = tol * tol;
            for (int i = 0; i < model.Nodes.Count; i++)
            {
                Node n = model.Nodes[i];
                double dx = n.X - x;
                double dy = n.Y - y;
                if (dx * dx + dy * dy < tol2) return true;
            }
            return false;
        }

        private static List<List<int>> BuildPolygons(List<Edge> edges, FemModel model)
        {
            List<List<int>> polys = new List<List<int>>();
            // adjacency: key=nodeId, value=list of (neighborId, edgeDirection-from-node-to-neighbor)
            // 我们保留方向：每条边 (p1,p2) 表示从p1到p2的有向边，但在多边形提取阶段我们需要按有向方式行走
            // 由于我们传入的边界边均为 "沿多边形CCW(外环)/CW(孔)" 方向组织（DXF导入按Circle/Polyline顺序建边），
            // 但细分/其他来源未必保持方向一致，因此用 "从当前点出发的所有候选下一顶点"，
            // 每走一条有向边就把该有向边标记为已用；若某点有多个候选，按转角最小（最右转）原则选择，保证沿环闭合行走。
            //
            // 实现：用 HashSet 存"有向边字符串"(cur->next)，每个无向边允许两个方向各走一次；
            // 但边界简单地每个无向边只应属于一个多边形，因此我们同时跟踪无向边是否已走。

            Dictionary<int, List<int>> adj = new Dictionary<int, List<int>>();
            for (int i = 0; i < edges.Count; i++)
            {
                Edge e = edges[i];
                if (!adj.ContainsKey(e.P1)) adj[e.P1] = new List<int>();
                if (!adj.ContainsKey(e.P2)) adj[e.P2] = new List<int>();
                adj[e.P1].Add(e.P2);
                adj[e.P2].Add(e.P1);
            }

            // 记录"无向边已使用"：用p<q序
            HashSet<long> usedUndir = new HashSet<long>();
            // 记录"有向边已使用"：cur -> next
            HashSet<long> usedDir = new HashSet<long>();

            for (int i = 0; i < edges.Count; i++)
            {
                Edge startEdge = edges[i];
                if (IsUndirUsed(usedUndir, startEdge.P1, startEdge.P2)) continue;

                // 选起始方向：尝试两个方向，选取第一个未用的
                int startNode = startEdge.P1;
                int secondNode = startEdge.P2;
                if (IsDirUsed(usedDir, startNode, secondNode))
                {
                    startNode = startEdge.P2;
                    secondNode = startEdge.P1;
                }
                if (IsDirUsed(usedDir, startNode, secondNode)) continue;

                List<int> poly = new List<int>();
                int cur = startNode;
                int prev = -1;
                int next = secondNode;
                poly.Add(cur);

                int safety = 0;
                while (safety < edges.Count * 4 + 10)
                {
                    safety++;
                    MarkDirUsed(usedDir, cur, next);
                    MarkUndirUsed(usedUndir, cur, next);
                    poly.Add(next);

                    if (next == startNode) break; // 回到起点，闭合

                    // 找 next->? 的候选：所有邻居中除cur之外的，且有向边未被使用
                    Node pCur = model.GetNode(next);
                    List<int> nbs = adj.ContainsKey(next) ? adj[next] : new List<int>();
                    double bestCross = -1e30;
                    int bestNext = -1;
                    double bestDot = 0.0;
                    double inDx = pCur.X - model.GetNode(cur).X;
                    double inDy = pCur.Y - model.GetNode(cur).Y;
                    double inLen = Math.Sqrt(inDx * inDx + inDy * inDy) + 1e-14;
                    inDx /= inLen; inDy /= inLen;
                    for (int k = 0; k < nbs.Count; k++)
                    {
                        int cand = nbs[k];
                        if (cand == cur) continue;
                        if (IsDirUsed(usedDir, next, cand)) continue;
                        Node pCand = model.GetNode(cand);
                        double outDx = pCand.X - pCur.X;
                        double outDy = pCand.Y - pCur.Y;
                        double outLen = Math.Sqrt(outDx * outDx + outDy * outDy) + 1e-14;
                        outDx /= outLen; outDy /= outLen;
                        // cross>0 表示左转，cross<0 表示右转；我们选**左转最小（最紧左转，等价于右转最大）**
                        // 对于CCW外环，始终选"最左"转可以沿外环行走闭合
                        double cross = inDx * outDy - inDy * outDx;
                        double dot = inDx * outDx + inDy * outDy;
                        // 优先最左转（cross最大）；若共线则选dot最大（正向）
                        if (bestNext == -1 || cross > bestCross + 1e-12 || (Math.Abs(cross - bestCross) < 1e-12 && dot > bestDot))
                        {
                            bestCross = cross;
                            bestDot = dot;
                            bestNext = cand;
                        }
                    }
                    if (bestNext == -1) break; // 死路
                    prev = cur;
                    cur = next;
                    next = bestNext;
                }
                // 若尾端不等于startNode，说明走了死链，丢弃
                if (poly.Count >= 3 && poly[poly.Count - 1] == startNode)
                {
                    polys.Add(poly);
                }
            }
            return polys;
        }

        private static long EdgeKey(int a, int b)
        {
            // 编码成 64-bit：高32位较小id，低32位较大id，保证无序
            int lo = a < b ? a : b;
            int hi = a < b ? b : a;
            return ((long)lo << 32) | (uint)hi;
        }
        private static long DirKey(int a, int b)
        {
            // 有向边：高32位from，低32位to
            return ((long)a << 32) | (uint)b;
        }
        private static bool IsUndirUsed(HashSet<long> set, int a, int b) { return set.Contains(EdgeKey(a, b)); }
        private static bool IsDirUsed(HashSet<long> set, int a, int b) { return set.Contains(DirKey(a, b)); }
        private static void MarkUndirUsed(HashSet<long> set, int a, int b) { set.Add(EdgeKey(a, b)); }
        private static void MarkDirUsed(HashSet<long> set, int a, int b) { set.Add(DirKey(a, b)); }

        private static bool PointInsidePolygons(double x, double y, List<List<int>> polys, FemModel model)
        {
            // 找到包含该点的最大/最小多边形判定：简单实现——取面积最大的多边形为外边界，其余为孔
            int outer = -1; double outerArea = -1;
            List<double> areas = new List<double>();
            for (int i = 0; i < polys.Count; i++)
            {
                double a = PolygonSignedArea(polys[i], model);
                areas.Add(a);
                if (Math.Abs(a) > outerArea) { outerArea = Math.Abs(a); outer = i; }
            }
            if (outer < 0) return false;
            bool insideOuter = PointInPolygon(x, y, polys[outer], model);
            if (!insideOuter) return false;
            // 是否在孔内
            for (int i = 0; i < polys.Count; i++)
            {
                if (i == outer) continue;
                if (PointInPolygon(x, y, polys[i], model)) return false;
            }
            return true;
        }

        private static double PolygonSignedArea(List<int> poly, FemModel model)
        {
            double a = 0;
            for (int i = 0; i < poly.Count; i++)
            {
                Node p1 = model.GetNode(poly[i]);
                Node p2 = model.GetNode(poly[(i + 1) % poly.Count]);
                a += p1.X * p2.Y - p2.X * p1.Y;
            }
            return a * 0.5;
        }

        private static bool PointInPolygon(double x, double y, List<int> poly, FemModel model)
        {
            bool inside = false;
            int n = poly.Count;
            for (int i = 0, j = n - 1; i < n; j = i++)
            {
                Node pi = model.GetNode(poly[i]);
                Node pj = model.GetNode(poly[j]);
                bool intersect = ((pi.Y > y) != (pj.Y > y)) &&
                                 (x < (pj.X - pi.X) * (y - pi.Y) / (pj.Y - pi.Y + 1e-12) + pi.X);
                if (intersect) inside = !inside;
            }
            return inside;
        }

        /// <summary>
        /// 两条边是否相同（不计方向）
        /// </summary>
        private static bool EdgeEqualUnordered(Edge a, Edge b)
        {
            return (a.P1 == b.P1 && a.P2 == b.P2) || (a.P1 == b.P2 && a.P2 == b.P1);
        }

        /// <summary>
        /// 对缓存中的边界边做等距细分：长度 > 1.1*h 的边按 ceil(L/h) 段插入中间节点，
        /// 并同步替换 BoundaryEdgeCache 中的边；中间节点线性插值两端点的BC。
        /// 这一步保证边界节点密度与内部一致，避免出现"拉面条"异形三角形。
        /// </summary>
        public static void SubdivideBoundaryEdges(FemModel model, double h)
        {
            List<Edge> boundary = BoundaryEdgeCache.Get(model);
            if (boundary == null || boundary.Count == 0) return;

            List<Edge> newEdges = new List<Edge>();
            int maxId = 0;
            for (int i = 0; i < model.Nodes.Count; i++)
            {
                if (model.Nodes[i].Id > maxId) maxId = model.Nodes[i].Id;
            }

            // 仅允许在"当前边自身"的端点附近做节点合并，禁止跨边合并到其他边界上的节点，
            // 否则会在两条不相邻边之间连出边，破坏边界环拓扑。
            for (int ei = 0; ei < boundary.Count; ei++)
            {
                Edge e = boundary[ei];
                Node n1 = model.GetNode(e.P1);
                Node n2 = model.GetNode(e.P2);
                if (n1 == null || n2 == null) { newEdges.Add(e); continue; }
                double dx = n2.X - n1.X;
                double dy = n2.Y - n1.Y;
                double L = Math.Sqrt(dx * dx + dy * dy);
                if (L < 1.1 * h) { newEdges.Add(e); continue; }
                int m = Math.Max(2, (int)Math.Ceiling(L / h));
                int prevId = e.P1;
                for (int k = 1; k < m; k++)
                {
                    double t = (double)k / m;
                    double x = n1.X + dx * t;
                    double y = n1.Y + dy * t;
                    // 仅检测距当前插入点极近（<0.3h）的已有点，若是本边端点或其他边中间点就合并；
                    // 由于细分是顺序执行的，其他边尚未细分，不会意外跨边合并。
                    int nearId = -1; double nearD2 = (0.3 * h) * (0.3 * h);
                    for (int ni = 0; ni < model.Nodes.Count; ni++)
                    {
                        Node nn = model.Nodes[ni];
                        double ndx = nn.X - x, ndy = nn.Y - y;
                        double nd2 = ndx * ndx + ndy * ndy;
                        if (nd2 < nearD2) { nearD2 = nd2; nearId = nn.Id; }
                    }
                    if (nearId > 0 && nearId != prevId)
                    {
                        newEdges.Add(new Edge(prevId, nearId));
                        prevId = nearId;
                        continue;
                    }
                    maxId++;
                    Node nm = new Node(maxId, x, y);
                    if (n1.BCType != 0 && n2.BCType != 0 && n1.BCType == n2.BCType)
                    {
                        nm.BCType = n1.BCType;
                        nm.BCValue = n1.BCValue + (n2.BCValue - n1.BCValue) * t;
                    }
                    else if (n1.BCType != 0 && n2.BCType != 0)
                    {
                        nm.BCType = 1;
                        nm.BCValue = n1.BCValue + (n2.BCValue - n1.BCValue) * t;
                    }
                    else
                    {
                        nm.BCType = 0;
                        nm.BCValue = 0.0;
                    }
                    model.Nodes.Add(nm);
                    newEdges.Add(new Edge(prevId, maxId));
                    prevId = maxId;
                }
                newEdges.Add(new Edge(prevId, e.P2));
            }
            BoundaryEdgeCache.Put(model, newEdges);
        }

        /// <summary>
        /// 计算单元面积与形心
        /// </summary>
        public static void ComputeElementGeom(Element e, FemModel model)
        {
            Node n1 = model.GetNode(e.N1);
            Node n2 = model.GetNode(e.N2);
            Node n3 = model.GetNode(e.N3);
            e.Area = 0.5 * Math.Abs((n2.X - n1.X) * (n3.Y - n1.Y) - (n3.X - n1.X) * (n2.Y - n1.Y));
            e.Cx = (n1.X + n2.X + n3.X) / 3.0;
            e.Cy = (n1.Y + n2.Y + n3.Y) / 3.0;
        }

        /// <summary>
        /// 移除落在域外的单元（质心不在任意多边形内部）
        /// </summary>
        public static void RemoveOutsideElements(FemModel model)
        {
            List<Edge> boundary = BoundaryEdgeCache.Get(model);
            if (boundary == null || boundary.Count == 0) return;
            List<List<int>> polys = BuildPolygons(boundary, model);
            for (int i = model.Elements.Count - 1; i >= 0; i--)
            {
                Element e = model.Elements[i];
                if (!PointInsidePolygons(e.Cx, e.Cy, polys, model))
                    model.Elements.RemoveAt(i);
            }
            // 重新编号ID
            for (int i = 0; i < model.Elements.Count; i++)
            {
                model.Elements[i].Id = i + 1;
            }
        }
    }
}
