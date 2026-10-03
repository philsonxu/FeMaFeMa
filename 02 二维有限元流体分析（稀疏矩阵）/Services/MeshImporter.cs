using System;
using System.Collections;
using System.Collections.Generic;
using System.IO;
using System.Globalization;
using System.Reflection;
using netDxf;
using netDxf.Entities;
using Fem2DFluid.Models;

namespace Fem2DFluid.Services
{
    /// <summary>
    /// 网格/几何导入器：支持自定义节点-单元文本格式、DXF边界线导入
    /// </summary>
    public class MeshImporter
    {
        /// <summary>
        /// 从自定义文本文件导入节点和单元
        /// 文件格式（#开头为注释，空行忽略）：
        /// TITLE 模型标题
        /// MATERIALS 材料数
        ///   id name k source
        /// NODES 节点数
        ///   id x y [bcType bcValue]
        /// ELEMENTS 单元数
        ///   id n1 n2 n3 matId
        /// EDGES 边界边数(可选)
        ///   p1 p2 bcType bcValue
        /// </summary>
        public static FemModel FromTextFile(string filePath)
        {
            FemModel model = new FemModel();
            string[] lines = File.ReadAllLines(filePath);
            string section = "";
            int i = 0;
            while (i < lines.Length)
            {
                string raw = lines[i].Trim();
                i++;
                if (raw.Length == 0 || raw.StartsWith("#")) continue;
                string upper = raw.ToUpper();
                if (upper.StartsWith("TITLE"))
                {
                    model.Title = raw.Substring(5).Trim();
                    continue;
                }
                if (upper.StartsWith("PROBLEM"))
                {
                    string[] t = raw.Split(new char[] { ' ', '\t' }, StringSplitOptions.RemoveEmptyEntries);
                    if (t.Length >= 2)
                    {
                        int pt;
                        if (int.TryParse(t[1], out pt)) model.ProblemType = pt;
                    }
                    continue;
                }
                if (upper.StartsWith("MATERIALS"))
                {
                    section = "MATERIALS";
                    string[] cnt = raw.Split(new char[] { ' ', '\t' }, StringSplitOptions.RemoveEmptyEntries);
                    int n = int.Parse(cnt[1]);
                    for (int j = 0; j < n && i < lines.Length; j++, i++)
                    {
                        string line = lines[i].Trim();
                        while ((line.Length == 0 || line.StartsWith("#")) && i < lines.Length)
                        {
                            i++; if (i < lines.Length) line = lines[i].Trim();
                        }
                        string[] parts = line.Split(new char[] { ' ', '\t', ',' }, StringSplitOptions.RemoveEmptyEntries);
                        Material m = new Material();
                        m.Id = int.Parse(parts[0]);
                        m.Name = parts[1];
                        m.K = double.Parse(parts[2], CultureInfo.InvariantCulture);
                        m.Source = parts.Length >= 4 ? double.Parse(parts[3], CultureInfo.InvariantCulture) : 0.0;
                        model.Materials.Add(m);
                    }
                    continue;
                }
                if (upper.StartsWith("NODES"))
                {
                    section = "NODES";
                    string[] cnt = raw.Split(new char[] { ' ', '\t' }, StringSplitOptions.RemoveEmptyEntries);
                    int n = int.Parse(cnt[1]);
                    for (int j = 0; j < n && i < lines.Length; j++, i++)
                    {
                        string line = lines[i].Trim();
                        while ((line.Length == 0 || line.StartsWith("#")) && i < lines.Length)
                        {
                            i++; if (i < lines.Length) line = lines[i].Trim();
                        }
                        string[] parts = line.Split(new char[] { ' ', '\t', ',' }, StringSplitOptions.RemoveEmptyEntries);
                        Node node = new Node();
                        node.Id = int.Parse(parts[0]);
                        node.X = double.Parse(parts[1], CultureInfo.InvariantCulture);
                        node.Y = double.Parse(parts[2], CultureInfo.InvariantCulture);
                        if (parts.Length >= 4)
                            node.BCType = int.Parse(parts[3]);
                        if (parts.Length >= 5)
                            node.BCValue = double.Parse(parts[4], CultureInfo.InvariantCulture);
                        model.Nodes.Add(node);
                    }
                    continue;
                }
                if (upper.StartsWith("ELEMENTS"))
                {
                    section = "ELEMENTS";
                    string[] cnt = raw.Split(new char[] { ' ', '\t' }, StringSplitOptions.RemoveEmptyEntries);
                    int n = int.Parse(cnt[1]);
                    for (int j = 0; j < n && i < lines.Length; j++, i++)
                    {
                        string line = lines[i].Trim();
                        while ((line.Length == 0 || line.StartsWith("#")) && i < lines.Length)
                        {
                            i++; if (i < lines.Length) line = lines[i].Trim();
                        }
                        string[] parts = line.Split(new char[] { ' ', '\t', ',' }, StringSplitOptions.RemoveEmptyEntries);
                        Element e = new Element();
                        e.Id = int.Parse(parts[0]);
                        e.N1 = int.Parse(parts[1]);
                        e.N2 = int.Parse(parts[2]);
                        e.N3 = int.Parse(parts[3]);
                        e.MatId = int.Parse(parts[4]);
                        model.Elements.Add(e);
                    }
                    continue;
                }
                if (upper.StartsWith("EDGES"))
                {
                    section = "EDGES";
                    string[] cnt = raw.Split(new char[] { ' ', '\t' }, StringSplitOptions.RemoveEmptyEntries);
                    int n = int.Parse(cnt[1]);
                    for (int j = 0; j < n && i < lines.Length; j++, i++)
                    {
                        string line = lines[i].Trim();
                        while ((line.Length == 0 || line.StartsWith("#")) && i < lines.Length)
                        {
                            i++; if (i < lines.Length) line = lines[i].Trim();
                        }
                        string[] parts = line.Split(new char[] { ' ', '\t', ',' }, StringSplitOptions.RemoveEmptyEntries);
                        int p1 = int.Parse(parts[0]);
                        int p2 = int.Parse(parts[1]);
                        int bct = int.Parse(parts[2]);
                        double bcv = double.Parse(parts[3], CultureInfo.InvariantCulture);
                        Node n1 = model.GetNode(p1);
                        Node n2 = model.GetNode(p2);
                        if (n1 != null) { n1.BCType = bct; n1.BCValue = bcv; }
                        if (n2 != null) { n2.BCType = bct; n2.BCValue = bcv; }
                    }
                    continue;
                }
            }
            // 若没给材料，默认加k=1材料
            if (model.Materials.Count == 0)
            {
                model.Materials.Add(new Material(1, "Fluid", 1.0));
            }
            return model;
        }

        /// <summary>
        /// 从DXF文件导入边界线段，返回仅含边界节点（无单元），需后续自动剖分
        /// </summary>
        public static FemModel FromDxf(string filePath)
        {
            FemModel model = new FemModel();
            model.Title = Path.GetFileNameWithoutExtension(filePath);
            model.Materials.Add(new Material(1, "Fluid", 1.0));

            DxfDocument dxf = DxfDocument.Load(filePath);
            Dictionary<string, int> vertexMap = new Dictionary<string, int>();
            List<Edge> edges = new List<Edge>();
            int nodeId = 1;

            foreach (EntityObject ent in dxf.Entities.All)
            {
                if (ent is Line)
                {
                    Line ln = (Line)ent;
                    int p1 = AddVertex(model, vertexMap, ref nodeId, ln.StartPoint.X, ln.StartPoint.Y, 0, 0);
                    int p2 = AddVertex(model, vertexMap, ref nodeId, ln.EndPoint.X, ln.EndPoint.Y, 0, 0);
                    edges.Add(new Edge(p1, p2));
                }
                else if (ent.GetType().Name == "LwPolyline")
                {
                    // 使用反射兼容 netDxf 各版本（避免 LwPolylineVertex 类型依赖）
                    List<int> ids = CollectVertexIds(ent, model, vertexMap, ref nodeId);
                    bool closed = (bool)GetPropValue(ent, "IsClosed");
                    ConnectChain(edges, ids, closed);
                }
                else if (ent.GetType().Name == "Polyline2D")
                {
                    List<int> ids = CollectVertexIds(ent, model, vertexMap, ref nodeId);
                    bool closed = (bool)GetPropValue(ent, "IsClosed");
                    ConnectChain(edges, ids, closed);
                }
                else if (ent is Circle)
                {
                    Circle c = (Circle)ent;
                    int np = 36;
                    List<int> circle = new List<int>();
                    for (int k = 0; k < np; k++)
                    {
                        double theta = 2 * Math.PI * k / np;
                        double x = c.Center.X + c.Radius * Math.Cos(theta);
                        double y = c.Center.Y + c.Radius * Math.Sin(theta);
                        int p = AddVertex(model, vertexMap, ref nodeId, x, y, 0, 0);
                        circle.Add(p);
                    }
                    for (int k = 0; k < np; k++)
                    {
                        edges.Add(new Edge(circle[k], circle[(k + 1) % np]));
                    }
                }
                else if (ent is Arc)
                {
                    Arc a = (Arc)ent;
                    int np = Math.Max(8, (int)Math.Ceiling(Math.Abs(a.EndAngle - a.StartAngle) / 10.0));
                    List<int> arc = new List<int>();
                    for (int k = 0; k <= np; k++)
                    {
                        double ang = a.StartAngle * Math.PI / 180.0 + (a.EndAngle - a.StartAngle) * Math.PI / 180.0 * k / np;
                        double x = a.Center.X + a.Radius * Math.Cos(ang);
                        double y = a.Center.Y + a.Radius * Math.Sin(ang);
                        int p = AddVertex(model, vertexMap, ref nodeId, x, y, 0, 0);
                        arc.Add(p);
                    }
                    for (int k = 0; k < arc.Count - 1; k++)
                    {
                        edges.Add(new Edge(arc[k], arc[k + 1]));
                    }
                }
            }
            // 把边界边存入模型Tag（用Source临时？不，我们给模型加Tags字典太麻烦，这里只是说明需要剖分）
            // 这里直接通过返回的模型，Elements为空表，Edges通过Tag传递。
            // 简化：把边信息写入到一个特殊变量：模型.Materials[0].Source=edges的数量作为标记，并另外通过Tag属性传递
            // 为简单起见，我们把边界作为临时"线框"存在一个附加属性：使用自定义的材料源项标记
            // 因为模型没有Edges列表，我们临时在节点BCType=1标记为边界节点，供剖分器使用
            // 同时我们临时记录边：通过给模型Title添加标记，并在外部通过节点位置重建（更稳）
            // 实际：给边界节点标记BCType=1（默认壁面 Dirichlet=0），然后剖分后再由用户修改边界条件
            HashSet<int> boundaryNodes = new HashSet<int>();
            for (int k = 0; k < edges.Count; k++)
            {
                boundaryNodes.Add(edges[k].P1);
                boundaryNodes.Add(edges[k].P2);
            }
            foreach (int nid in boundaryNodes)
            {
                Node n = model.GetNode(nid);
                // 先不指定具体BC类型（保持0），由调用方在几何加载后调用SetBC函数设置真正边界条件
                if (n != null) { n.BCType = 0; n.BCValue = 0.0; }
            }
            // edges数据：存入一个静态缓存（避免修改模型类）——使用MeshTriangulator可调用的全局方法
            BoundaryEdgeCache.Put(model, edges);
            return model;
        }

        /// <summary>
        /// 反射读取对象属性值——用于跨netDxf版本兼容，避免直接引用顶点类型名
        /// </summary>
        private static object GetPropValue(object obj, string name)
        {
            PropertyInfo pi = obj.GetType().GetProperty(name, BindingFlags.Public | BindingFlags.Instance);
            if (pi == null) return null;
            return pi.GetValue(obj, null);
        }

        /// <summary>
        /// 通用方法：从任意多段线类型（LwPolyline/Polyline/Polyline2D/Polyline3D）反射取顶点坐标
        /// 兼容顶点上的Position.Vector2、Position.X/Y、以及直接X/Y两种结构
        /// </summary>
        private static List<int> CollectVertexIds(object polyline, FemModel model, Dictionary<string, int> vertexMap, ref int nodeId)
        {
            List<int> ids = new List<int>();
            object vxs = GetPropValue(polyline, "Vertexes");
            if (vxs == null) return ids;
            IEnumerable en = vxs as IEnumerable;
            if (en == null) return ids;
            foreach (object vx in en)
            {
                double x = 0.0, y = 0.0;
                object pos = GetPropValue(vx, "Position");
                if (pos != null)
                {
                    object ox = GetPropValue(pos, "X");
                    object oy = GetPropValue(pos, "Y");
                    if (ox != null) x = (double)Convert.ChangeType(ox, typeof(double));
                    if (oy != null) y = (double)Convert.ChangeType(oy, typeof(double));
                }
                else
                {
                    // 某些版本顶点直接是Vector2/X/Y成员
                    object ox = GetPropValue(vx, "X");
                    object oy = GetPropValue(vx, "Y");
                    if (ox != null) x = (double)Convert.ChangeType(ox, typeof(double));
                    if (oy != null) y = (double)Convert.ChangeType(oy, typeof(double));
                }
                int pid = AddVertex(model, vertexMap, ref nodeId, x, y, 0, 0);
                ids.Add(pid);
            }
            return ids;
        }

        /// <summary>
        /// 将顶点id序列连为顺序边；closed为true时首尾相连形成闭合环
        /// </summary>
        private static void ConnectChain(List<Edge> edges, List<int> ids, bool closed)
        {
            for (int k = 0; k < ids.Count - 1; k++)
            {
                edges.Add(new Edge(ids[k], ids[k + 1]));
            }
            if (closed && ids.Count > 2)
            {
                edges.Add(new Edge(ids[ids.Count - 1], ids[0]));
            }
        }

        private static int AddVertex(FemModel model, Dictionary<string, int> vertexMap, ref int nodeId, double x, double y, int bct, double bcv)
        {
            string key = string.Format(CultureInfo.InvariantCulture, "{0:F6},{1:F6}", x, y);
            int id;
            if (vertexMap.TryGetValue(key, out id)) return id;
            id = nodeId++;
            Node n = new Node(id, x, y);
            n.BCType = bct;
            n.BCValue = bcv;
            model.Nodes.Add(n);
            vertexMap.Add(key, id);
            return id;
        }
    }

    /// <summary>
    /// 简单缓存：因为FemModel没有Edges字段，导入DXF时临时保存边界边，供三角剖分使用
    /// </summary>
    public static class BoundaryEdgeCache
    {
        private static Dictionary<FemModel, List<Edge>> _cache = new Dictionary<FemModel, List<Edge>>();

        public static void Put(FemModel model, List<Edge> edges)
        {
            _cache[model] = edges;
        }

        public static List<Edge> Get(FemModel model)
        {
            List<Edge> e;
            if (_cache.TryGetValue(model, out e)) return e;
            return null;
        }

        public static void Clear(FemModel model)
        {
            _cache.Remove(model);
        }
    }
}
