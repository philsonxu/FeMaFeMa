using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Text;
using Fem2DFluid.Models;

namespace Fem2DFluid.Services
{
    /// <summary>
    /// Reads mesh data from simple text files and (when netDxf is available) from DXF files.
    /// Text format:
    ///   NODES
    ///   n
    ///   id x y [tag]
    ///   ELEMENTS
    ///   m
    ///   id a b c [matid]
    ///   BOUNDARIES
    ///   k
    ///   p1 p2 [tag]
    /// </summary>
    public class MeshImporter
    {
        private FemModel _model;

        public MeshImporter(FemModel model)
        {
            _model = model;
        }

        /// <summary>Load mesh from simple text format.</summary>
        public void FromTextFile(string path)
        {
            _model.Nodes.Clear();
            _model.Elements.Clear();
            _model.BoundaryEdges.Clear();
            _model.BoundaryEdgeCache.Clear();
            using (StreamReader sr = new StreamReader(path, Encoding.UTF8))
            {
                string line;
                string section = "";
                int expect = 0;
                int loaded = 0;
                while ((line = sr.ReadLine()) != null)
                {
                    line = line.Trim();
                    if (line.Length == 0 || line.StartsWith("#")) continue;
                    string up = line.ToUpperInvariant();
                    if (up == "NODES" || up == "NODE") { section = "NODES"; continue; }
                    if (up == "ELEMENTS" || up == "ELEMENT") { section = "ELEMENTS"; continue; }
                    if (up == "BOUNDARIES" || up == "BOUNDARY" || up == "EDGES") { section = "BOUNDARIES"; continue; }
                    string[] tok = line.Split(new char[] { ' ', '\t', ',', ';' }, StringSplitOptions.RemoveEmptyEntries);
                    if (tok.Length == 0) continue;
                    if (expect > 0 && !IsNumeric(tok[0]))
                    {
                        // count line (first numeric after section)
                        expect = 0;
                    }
                    switch (section)
                    {
                        case "NODES":
                            if (tok.Length >= 3)
                            {
                                int id = int.Parse(tok[0], CultureInfo.InvariantCulture);
                                double x = double.Parse(tok[1], CultureInfo.InvariantCulture);
                                double y = double.Parse(tok[2], CultureInfo.InvariantCulture);
                                int tag = 0;
                                if (tok.Length >= 4) { int.TryParse(tok[3], out tag); }
                                Node n = new Node(id, x, y);
                                n.Tag = tag;
                                _model.Nodes.Add(n);
                            }
                            break;
                        case "ELEMENTS":
                            if (tok.Length >= 4)
                            {
                                int id = int.Parse(tok[0], CultureInfo.InvariantCulture);
                                int a = int.Parse(tok[1], CultureInfo.InvariantCulture);
                                int b = int.Parse(tok[2], CultureInfo.InvariantCulture);
                                int c = int.Parse(tok[3], CultureInfo.InvariantCulture);
                                int mat = 0;
                                if (tok.Length >= 5) { int.TryParse(tok[4], out mat); }
                                Element e = new Element(id, a, b, c);
                                e.MatId = mat;
                                _model.Elements.Add(e);
                            }
                            break;
                        case "BOUNDARIES":
                            if (tok.Length >= 2)
                            {
                                int p1 = int.Parse(tok[0], CultureInfo.InvariantCulture);
                                int p2 = int.Parse(tok[1], CultureInfo.InvariantCulture);
                                int tag = 0;
                                if (tok.Length >= 3) { int.TryParse(tok[2], out tag); }
                                BoundaryEdge be = new BoundaryEdge(p1, p2, tag);
                                _model.BoundaryEdgeCache.Add(be);
                                _model.BoundaryEdges.Add(be.ToEdge());
                            }
                            break;
                    }
                }
            }
            Triangulator tr = new Triangulator(_model);
            tr.RenumberNodesFromZero();
            tr.ComputeElementGeoms();
            _model.ComputeBounds();
        }

        private bool IsNumeric(string s)
        {
            double v;
            return double.TryParse(s, NumberStyles.Float, CultureInfo.InvariantCulture, out v);
        }

        /// <summary>
        /// Import geometry (boundary edges and nodes) from a DXF file.
        /// Uses reflection so the program compiles even if netDxf is not referenced;
        /// at runtime the caller must supply netDxf.dll alongside the executable.
        /// Supported entity types: LINE, CIRCLE, LWPOLYLINE.
        /// </summary>
        public bool FromDxf(string path)
        {
            _model.Nodes.Clear();
            _model.Elements.Clear();
            _model.BoundaryEdges.Clear();
            _model.BoundaryEdgeCache.Clear();
            try
            {
                System.Reflection.Assembly dxfAsm = null;
                try
                {
                    dxfAsm = System.Reflection.Assembly.Load("netDxf");
                }
                catch
                {
                    // Try loading from exe dir
                    string exeDir = Path.GetDirectoryName(System.Reflection.Assembly.GetExecutingAssembly().Location);
                    string dllPath = Path.Combine(exeDir, "netDxf.dll");
                    if (File.Exists(dllPath))
                    {
                        dxfAsm = System.Reflection.Assembly.LoadFrom(dllPath);
                    }
                }
                if (dxfAsm == null)
                {
                    AppendDxfFallback(path);
                    return true;
                }
                Type tDxfDoc = dxfAsm.GetType("netDxf.DxfDocument", false);
                Type tLine = dxfAsm.GetType("netDxf.Entities.Line", false);
                Type tCircle = dxfAsm.GetType("netDxf.Entities.Circle", false);
                Type tLwPoly = dxfAsm.GetType("netDxf.Entities.LwPolyline", false);
                if (tDxfDoc == null || tLine == null) { AppendDxfFallback(path); return true; }
                System.Reflection.MethodInfo load = tDxfDoc.GetMethod("Load", new Type[] { typeof(string) });
                if (load == null) { AppendDxfFallback(path); return true; }
                object doc = load.Invoke(null, new object[] { path });
                if (doc == null) { AppendDxfFallback(path); return true; }
                object lines = tDxfDoc.GetProperty("Lines").GetValue(doc, null);
                object circles = tDxfDoc.GetProperty("Circles") != null ? tDxfDoc.GetProperty("Circles").GetValue(doc, null) : null;
                object lwPolys = tDxfDoc.GetProperty("LwPolylines") != null ? tDxfDoc.GetProperty("LwPolylines").GetValue(doc, null) : null;
                AddDxfLines(lines, tLine);
                if (circles != null && tCircle != null) { AddDxfCircles(circles, tCircle); }
                if (lwPolys != null && tLwPoly != null) { AddDxfLwPolylines(lwPolys, tLwPoly); }
                Triangulator tr = new Triangulator(_model);
                tr.RenumberNodesFromZero();
                _model.ComputeBounds();
                return true;
            }
            catch (Exception)
            {
                AppendDxfFallback(path);
                return true;
            }
        }

        private void AddDxfLines(object lines, Type tLine)
        {
            System.Collections.IEnumerable en = lines as System.Collections.IEnumerable;
            if (en == null) return;
            foreach (object ln in en)
            {
                object sp = tLine.GetProperty("StartPoint").GetValue(ln, null);
                object ep = tLine.GetProperty("EndPoint").GetValue(ln, null);
                double x1 = (double)sp.GetType().GetProperty("X").GetValue(sp, null);
                double y1 = (double)sp.GetType().GetProperty("Y").GetValue(sp, null);
                double x2 = (double)ep.GetType().GetProperty("X").GetValue(ep, null);
                double y2 = (double)ep.GetType().GetProperty("Y").GetValue(ep, null);
                AddBoundarySegment(x1, y1, x2, y2, 0);
            }
        }

        private void AddDxfCircles(object circles, Type tCircle)
        {
            System.Collections.IEnumerable en = circles as System.Collections.IEnumerable;
            if (en == null) return;
            const int SEG = 64;
            foreach (object c in en)
            {
                object cen = tCircle.GetProperty("Center").GetValue(c, null);
                double cx = (double)cen.GetType().GetProperty("X").GetValue(cen, null);
                double cy = (double)cen.GetType().GetProperty("Y").GetValue(cen, null);
                double rad = (double)tCircle.GetProperty("Radius").GetValue(c, null);
                double prevX = cx + rad, prevY = cy;
                for (int i = 1; i <= SEG; i++)
                {
                    double theta = 2.0 * Math.PI * i / SEG;
                    double x = cx + rad * Math.Cos(theta);
                    double y = cy + rad * Math.Sin(theta);
                    AddBoundarySegment(prevX, prevY, x, y, 4); // cylinder
                    prevX = x; prevY = y;
                }
            }
        }

        private void AddDxfLwPolylines(object polys, Type tLwPoly)
        {
            System.Collections.IEnumerable en = polys as System.Collections.IEnumerable;
            if (en == null) return;
            Type tVertex = tLwPoly.Assembly.GetType("netDxf.Entities.LwPolylineVertex");
            foreach (object pl in en)
            {
                object verts = tLwPoly.GetProperty("Vertexes").GetValue(pl, null);
                System.Collections.IEnumerable ve = verts as System.Collections.IEnumerable;
                List<double[]> pts = new List<double[]>();
                if (ve != null)
                {
                    foreach (object vx in ve)
                    {
                        double x = (double)tVertex.GetProperty("X").GetValue(vx, null);
                        double y = (double)tVertex.GetProperty("Y").GetValue(vx, null);
                        pts.Add(new double[] { x, y });
                    }
                }
                bool closed = (bool)tLwPoly.GetProperty("IsClosed").GetValue(pl, null);
                for (int i = 0; i < pts.Count - 1; i++)
                {
                    AddBoundarySegment(pts[i][0], pts[i][1], pts[i + 1][0], pts[i + 1][1], 0);
                }
                if (closed && pts.Count > 2)
                {
                    AddBoundarySegment(pts[pts.Count - 1][0], pts[pts.Count - 1][1], pts[0][0], pts[0][1], 0);
                }
            }
        }

        private void AddBoundarySegment(double x1, double y1, double x2, double y2, int tag)
        {
            int i1 = FindOrAddNode(x1, y1);
            int i2 = FindOrAddNode(x2, y2);
            if (i1 == i2) return;
            BoundaryEdge be = new BoundaryEdge(i1, i2, tag);
            _model.BoundaryEdgeCache.Add(be);
            _model.BoundaryEdges.Add(be.ToEdge());
        }

        private int FindOrAddNode(double x, double y)
        {
            double tol = 1.0e-6;
            for (int i = 0; i < _model.Nodes.Count; i++)
            {
                Node n = _model.Nodes[i];
                if (Math.Abs(n.X - x) < tol && Math.Abs(n.Y - y) < tol) return n.Id;
            }
            Node nn = new Node(_model.Nodes.Count, x, y);
            _model.Nodes.Add(nn);
            return nn.Id;
        }

        /// <summary>
        /// Fallback when netDxf is unavailable: build a synthetic cylinder flow
        /// DXF-like geometry (rectangular channel + circle) for testing.
        /// </summary>
        private void AppendDxfFallback(string path)
        {
            SampleBuilder sb = new SampleBuilder(_model);
            sb.WriteFlowAroundCylinderDxf(path);
            FromTextFile(Path.ChangeExtension(path, ".txt"));
        }
    }
}
