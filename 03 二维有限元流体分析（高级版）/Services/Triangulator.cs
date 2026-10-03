using System;
using System.Collections.Generic;
using Fem2DFluid.Models;

namespace Fem2DFluid.Services
{
    /// <summary>
    /// Delaunay triangulation (Bowyer-Watson) with boundary refinement,
    /// interior structured-point generation, polygon-cropping, and LT6
    /// mid-edge node insertion/removal. All indices are kept 0-based
    /// and contiguous via <see cref="RenumberNodesFromZero"/>.
    /// </summary>
    public class Triangulator
    {
        private FemModel _model;
        private Random _rand;

        public Triangulator(FemModel model)
        {
            _model = model;
            _rand = new Random(12345);
        }

        // ========== Bowyer-Watson Delaunay ==========

        private class SuperTriangle
        {
            public int A, B, C;
        }

        private bool InCircumcircle(Node p1, Node p2, Node p3, Node p)
        {
            double ax = p1.X - p.X, ay = p1.Y - p.Y;
            double bx = p2.X - p.X, by = p2.Y - p.Y;
            double cx = p3.X - p.X, cy = p3.Y - p.Y;
            double det = (ax * ax + ay * ay) * (bx * cy - cx * by)
                       - (bx * bx + by * by) * (ax * cy - cx * ay)
                       + (cx * cx + cy * cy) * (ax * by - bx * ay);
            double ccw = (p2.X - p1.X) * (p3.Y - p1.Y) - (p3.X - p1.X) * (p2.Y - p1.Y);
            return (det > 0 && ccw > 0) || (det < 0 && ccw < 0);
        }

        private bool PointInTriangle(Node p, Node a, Node b, Node c)
        {
            double d1 = (p.X - b.X) * (a.Y - b.Y) - (a.X - b.X) * (p.Y - b.Y);
            double d2 = (p.X - c.X) * (b.Y - c.Y) - (b.X - c.X) * (p.Y - c.Y);
            double d3 = (p.X - a.X) * (c.Y - a.Y) - (c.X - a.X) * (p.Y - a.Y);
            bool neg = (d1 < 0) || (d2 < 0) || (d3 < 0);
            bool pos = (d1 > 0) || (d2 > 0) || (d3 > 0);
            return !(neg && pos);
        }

        /// <summary>
        /// Run Bowyer-Watson on the current Nodes list. Existing Elements are cleared.
        /// </summary>
        public void Delaunay()
        {
            List<Node> nodes = _model.Nodes;
            if (nodes.Count < 3) { return; }
            _model.ComputeBounds();
            double dx = _model.Xmax - _model.Xmin;
            double dy = _model.Ymax - _model.Ymin;
            double dmax = Math.Max(dx, dy) * 2.0;
            double cx = (_model.Xmin + _model.Xmax) * 0.5;
            double cy = (_model.Ymin + _model.Ymax) * 0.5;

            Node sA = new Node(-10, cx - dmax * 2, cy - dmax);
            Node sB = new Node(-11, cx + dmax * 2, cy - dmax);
            Node sC = new Node(-12, cx, cy + dmax * 2);
            List<Node> pn = new List<Node>(nodes);
            pn.Add(sA); pn.Add(sB); pn.Add(sC);
            int iA = pn.Count - 3, iB = pn.Count - 2, iC = pn.Count - 1;

            List<int[]> tris = new List<int[]>();
            tris.Add(new int[] { iA, iB, iC });

            for (int ip = 0; ip < nodes.Count; ip++)
            {
                Node p = nodes[ip];
                List<int[]> bad = new List<int[]>();
                List<Edge> poly = new List<Edge>();
                for (int t = tris.Count - 1; t >= 0; t--)
                {
                    int[] tri = tris[t];
                    if (InCircumcircle(pn[tri[0]], pn[tri[1]], pn[tri[2]], p))
                    {
                        bad.Add(tri);
                        tris.RemoveAt(t);
                    }
                }
                // Find polygon boundary of bad triangles
                HashSet<Edge> edgeCount = new HashSet<Edge>();
                Dictionary<Edge, int> edgeDups = new Dictionary<Edge, int>();
                for (int b = 0; b < bad.Count; b++)
                {
                    int[] tri = bad[b];
                    Edge e1 = new Edge(tri[0], tri[1]);
                    Edge e2 = new Edge(tri[1], tri[2]);
                    Edge e3 = new Edge(tri[2], tri[0]);
                    IncEdge(edgeDups, e1);
                    IncEdge(edgeDups, e2);
                    IncEdge(edgeDups, e3);
                }
                foreach (KeyValuePair<Edge, int> kv in edgeDups)
                {
                    if (kv.Value == 1) { poly.Add(kv.Key); }
                }
                for (int k = 0; k < poly.Count; k++)
                {
                    Edge e = poly[k];
                    tris.Add(new int[] { e.P1, e.P2, ip });
                }
            }

            // Remove triangles touching super nodes
            List<Element> elements = new List<Element>();
            int eid = 0;
            for (int t = 0; t < tris.Count; t++)
            {
                int[] tri = tris[t];
                if (tri[0] >= nodes.Count || tri[1] >= nodes.Count || tri[2] >= nodes.Count) { continue; }
                Element el = new Element(eid++, tri[0], tri[1], tri[2]);
                el.Order = 1;
                elements.Add(el);
            }
            _model.Elements = elements;
        }

        private void IncEdge(Dictionary<Edge, int> d, Edge e)
        {
            int c;
            if (d.TryGetValue(e, out c)) { d[e] = c + 1; }
            else { d[e] = 1; }
        }

        // ========== Boundary subdivision ==========

        /// <summary>
        /// Subdivide all boundary segments in BoundaryEdgeCache so that no segment
        /// is longer than h. New boundary nodes are appended to Nodes and Tag is
        /// inherited from the parent segment.
        /// </summary>
        public int SubdivideBoundaryEdges(double h)
        {
            if (h <= 0) { return 0; }
            List<BoundaryEdge> newBe = new List<BoundaryEdge>();
            int added = 0;
            for (int i = 0; i < _model.BoundaryEdgeCache.Count; i++)
            {
                BoundaryEdge be = _model.BoundaryEdgeCache[i];
                Node p1 = FindNodeById(be.P1);
                Node p2 = FindNodeById(be.P2);
                if (p1 == null || p2 == null) { newBe.Add(be); continue; }
                double dx = p2.X - p1.X, dy = p2.Y - p1.Y;
                double len = Math.Sqrt(dx * dx + dy * dy);
                int nSeg = Math.Max(1, (int)Math.Ceiling(len / h));
                if (nSeg == 1) { newBe.Add(be); continue; }
                int prevId = be.P1;
                for (int s = 1; s < nSeg; s++)
                {
                    double t = (double)s / nSeg;
                    double x = p1.X + t * dx;
                    double y = p1.Y + t * dy;
                    Node nm = new Node(_model.Nodes.Count, x, y);
                    nm.Tag = be.Tag;
                    nm.NodeType = 1; // mid-edge
                    // If tag indicates Dirichlet, propagate BC
                    if (be.Tag == 1) { nm.BCType = 1; nm.BCValue = 0.0; }
                    _model.Nodes.Add(nm);
                    BoundaryEdge seg = new BoundaryEdge(prevId, nm.Id, be.Tag);
                    seg.Value = be.Value;
                    newBe.Add(seg);
                    prevId = nm.Id;
                    added++;
                }
                BoundaryEdge last = new BoundaryEdge(prevId, be.P2, be.Tag);
                last.Value = be.Value;
                newBe.Add(last);
            }
            _model.BoundaryEdgeCache = newBe;
            return added;
        }

        /// <summary>
        /// Generate a lattice of interior points with spacing h, lying strictly
        /// inside the domain polygon (not on boundary edges).
        /// </summary>
        public int GenerateStructuredPoints(double h)
        {
            if (h <= 0) { return 0; }
            List<double[]> polygon = BuildOuterPolygon();
            if (polygon == null || polygon.Count < 3) { return 0; }
            _model.ComputeBounds();
            int added = 0;
            double x0 = _model.Xmin + h * 0.5;
            double y0 = _model.Ymin + h * 0.5;
            int nx = (int)Math.Ceiling((_model.Xmax - _model.Xmin) / h);
            int ny = (int)Math.Ceiling((_model.Ymax - _model.Ymin) / h);
            for (int j = 0; j < ny; j++)
            {
                double y = y0 + j * h;
                if (y > _model.Ymax - h * 0.1) { continue; }
                for (int i = 0; i < nx; i++)
                {
                    double x = x0 + i * h;
                    if (x > _model.Xmax - h * 0.1) { continue; }
                    if (!PointInPolygon(x, y, polygon)) { continue; }
                    if (NearBoundary(x, y, h * 0.45)) { continue; }
                    Node n = new Node(_model.Nodes.Count, x, y);
                    n.Tag = 0;
                    _model.Nodes.Add(n);
                    added++;
                }
            }
            return added;
        }

        private bool NearBoundary(double x, double y, double tol)
        {
            for (int i = 0; i < _model.BoundaryEdgeCache.Count; i++)
            {
                BoundaryEdge be = _model.BoundaryEdgeCache[i];
                Node p1 = FindNodeById(be.P1);
                Node p2 = FindNodeById(be.P2);
                if (p1 == null || p2 == null) continue;
                double dx = p2.X - p1.X, dy = p2.Y - p1.Y;
                double len2 = dx * dx + dy * dy;
                if (len2 < 1.0e-20) continue;
                double t = ((x - p1.X) * dx + (y - p1.Y) * dy) / len2;
                if (t < 0) t = 0; if (t > 1) t = 1;
                double px = p1.X + t * dx, py = p1.Y + t * dy;
                double d = Math.Sqrt((x - px) * (x - px) + (y - py) * (y - py));
                if (d < tol) return true;
            }
            return false;
        }

        /// <summary>
        /// Remove elements whose centroid lies outside the outer polygon or
        /// inside any hole polygon.
        /// </summary>
        public int RemoveOutsideElements()
        {
            List<double[]> outer = BuildOuterPolygon();
            List<List<double[]>> holes = BuildHolePolygons();
            if (outer == null || outer.Count < 3) { return 0; }
            ComputeElementGeoms();
            List<Element> keep = new List<Element>();
            int removed = 0;
            for (int i = 0; i < _model.Elements.Count; i++)
            {
                Element e = _model.Elements[i];
                double cx = e.Barycenter[0], cy = e.Barycenter[1];
                if (!PointInPolygon(cx, cy, outer)) { removed++; continue; }
                bool insideHole = false;
                for (int h = 0; h < holes.Count; h++)
                {
                    if (PointInPolygon(cx, cy, holes[h])) { insideHole = true; break; }
                }
                if (insideHole) { removed++; continue; }
                keep.Add(e);
            }
            _model.Elements = keep;
            return removed;
        }

        /// <summary>
        /// Insert mid-edge nodes on every element edge that is not on a Dirichlet
        /// boundary where the user specified an essential BC (those remain linear
        /// to preserve the boundary value profile). Used to convert CST mesh to LT6.
        /// </summary>
        public int AddMidEdgeNodes()
        {
            _model.BoundaryEdges.Clear();
            for (int i = 0; i < _model.BoundaryEdgeCache.Count; i++)
            {
                _model.BoundaryEdges.Add(_model.BoundaryEdgeCache[i].ToEdge());
            }
            Dictionary<Edge, int> midNodes = new Dictionary<Edge, int>();
            int added = 0;
            // First pass: create mid-edge nodes
            for (int ei = 0; ei < _model.Elements.Count; ei++)
            {
                Element e = _model.Elements[ei];
                Edge e1 = new Edge(e.B, e.C);
                Edge e2 = new Edge(e.C, e.A);
                Edge e3 = new Edge(e.A, e.B);
                EnsureMidNode(e1, midNodes, ref added);
                EnsureMidNode(e2, midNodes, ref added);
                EnsureMidNode(e3, midNodes, ref added);
                e.N4 = midNodes[e1];
                e.N5 = midNodes[e2];
                e.N6 = midNodes[e3];
                e.Order = 2;
            }
            // Set mid-edge node BC based on whether edge is a BC edge
            for (int i = 0; i < _model.BoundaryEdgeCache.Count; i++)
            {
                BoundaryEdge be = _model.BoundaryEdgeCache[i];
                Edge ek = be.ToEdge();
                int mid;
                if (midNodes.TryGetValue(ek, out mid))
                {
                    Node n1 = FindNodeById(be.P1), n2 = FindNodeById(be.P2), nm = FindNodeById(mid);
                    if (n1 != null && n2 != null && nm != null)
                    {
                        nm.Tag = be.Tag;
                        nm.BCType = n1.BCType == 1 && n2.BCType == 1 ? 1 : 0;
                        nm.BCValue = 0.5 * (n1.BCValue + n2.BCValue);
                    }
                }
            }
            return added;
        }

        private void EnsureMidNode(Edge e, Dictionary<Edge, int> midNodes, ref int added)
        {
            if (midNodes.ContainsKey(e)) return;
            Node p1 = FindNodeById(e.P1);
            Node p2 = FindNodeById(e.P2);
            if (p1 == null || p2 == null) return;
            Node nm = new Node(_model.Nodes.Count, 0.5 * (p1.X + p2.X), 0.5 * (p1.Y + p2.Y));
            nm.NodeType = 1;
            // Interpolate BC conservatively (Tag inherited later if needed)
            _model.Nodes.Add(nm);
            midNodes[e] = nm.Id;
            added++;
        }

        /// <summary>Remove all mid-edge nodes and revert to linear elements.</summary>
        public int RemoveMidEdgeNodes()
        {
            List<Node> newNodes = new List<Node>();
            Dictionary<int, int> remap = new Dictionary<int, int>();
            for (int i = 0; i < _model.Nodes.Count; i++)
            {
                Node n = _model.Nodes[i];
                if (n.NodeType == 1) { continue; }
                remap[n.Id] = newNodes.Count;
                Node nn = new Node(newNodes.Count, n.X, n.Y);
                nn.BCType = n.BCType;
                nn.BCValue = n.BCValue;
                nn.Tag = n.Tag;
                newNodes.Add(nn);
            }
            List<Element> newElems = new List<Element>();
            for (int i = 0; i < _model.Elements.Count; i++)
            {
                Element e = _model.Elements[i];
                Element ne = new Element(newElems.Count, remap[e.A], remap[e.B], remap[e.C]);
                ne.MatId = e.MatId;
                ne.Order = 1;
                ne.Area = e.Area;
                newElems.Add(ne);
            }
            List<BoundaryEdge> nbe = new List<BoundaryEdge>();
            for (int i = 0; i < _model.BoundaryEdgeCache.Count; i++)
            {
                BoundaryEdge be = _model.BoundaryEdgeCache[i];
                if (!remap.ContainsKey(be.P1) || !remap.ContainsKey(be.P2)) { continue; }
                BoundaryEdge nb = new BoundaryEdge(remap[be.P1], remap[be.P2], be.Tag);
                nb.Value = be.Value;
                nbe.Add(nb);
            }
            _model.Nodes = newNodes;
            _model.Elements = newElems;
            _model.BoundaryEdgeCache = nbe;
            _model.BoundaryEdges.Clear();
            for (int i = 0; i < nbe.Count; i++) { _model.BoundaryEdges.Add(nbe[i].ToEdge()); }
            RenumberNodesFromZero();
            return newNodes.Count;
        }

        /// <summary>Renumber nodes so their Ids are 0..n-1 in list order; rebuild element indices.</summary>
        public void RenumberNodesFromZero()
        {
            Dictionary<int, int> remap = new Dictionary<int, int>();
            for (int i = 0; i < _model.Nodes.Count; i++)
            {
                remap[_model.Nodes[i].Id] = i;
                _model.Nodes[i].Id = i;
            }
            for (int i = 0; i < _model.Elements.Count; i++)
            {
                Element e = _model.Elements[i];
                e.Id = i;
                e.A = remap[e.A]; e.B = remap[e.B]; e.C = remap[e.C];
                if (e.Order == 2 && e.N4 >= 0 && e.N5 >= 0 && e.N6 >= 0)
                {
                    e.N4 = remap[e.N4]; e.N5 = remap[e.N5]; e.N6 = remap[e.N6];
                }
            }
            for (int i = 0; i < _model.BoundaryEdgeCache.Count; i++)
            {
                BoundaryEdge be = _model.BoundaryEdgeCache[i];
                be.P1 = remap[be.P1];
                be.P2 = remap[be.P2];
            }
        }

        public void ComputeElementGeom(Element e)
        {
            Node a = _model.GetNode(e.A);
            Node b = _model.GetNode(e.B);
            Node c = _model.GetNode(e.C);
            double ar = Math.Abs(0.5 * ((b.X - a.X) * (c.Y - a.Y) - (c.X - a.X) * (b.Y - a.Y)));
            e.Area = ar;
            e.Barycenter[0] = (a.X + b.X + c.X) / 3.0;
            e.Barycenter[1] = (a.Y + b.Y + c.Y) / 3.0;
        }

        public void ComputeElementGeoms()
        {
            for (int i = 0; i < _model.Elements.Count; i++)
            {
                ComputeElementGeom(_model.Elements[i]);
            }
        }

        public Node FindNodeById(int id)
        {
            for (int i = 0; i < _model.Nodes.Count; i++)
            {
                if (_model.Nodes[i].Id == id) return _model.Nodes[i];
            }
            return null;
        }

        // ========== Polygon construction ==========

        /// <summary>
        /// Build closed polygon rings from BoundaryEdgeCache (simple chain;
        /// assumes boundary is well-formed and nodes connect end-to-end).
        /// Returns list of rings; first is outer, others are holes (sorted by area desc).
        /// </summary>
        public List<List<double[]>> BuildPolygons()
        {
            List<BoundaryEdge> edges = new List<BoundaryEdge>(_model.BoundaryEdgeCache);
            List<List<double[]>> rings = new List<List<double[]>>();
            while (edges.Count > 0)
            {
                BoundaryEdge first = edges[0];
                edges.RemoveAt(0);
                List<int> chain = new List<int>();
                chain.Add(first.P1);
                chain.Add(first.P2);
                bool grew = true;
                int safety = edges.Count + 2;
                while (grew && safety-- > 0)
                {
                    grew = false;
                    for (int i = edges.Count - 1; i >= 0; i--)
                    {
                        BoundaryEdge e = edges[i];
                        if (e.P1 == chain[chain.Count - 1])
                        {
                            chain.Add(e.P2); edges.RemoveAt(i); grew = true;
                        }
                        else if (e.P2 == chain[chain.Count - 1])
                        {
                            chain.Add(e.P1); edges.RemoveAt(i); grew = true;
                        }
                        else if (e.P1 == chain[0])
                        {
                            chain.Insert(0, e.P2); edges.RemoveAt(i); grew = true;
                        }
                        else if (e.P2 == chain[0])
                        {
                            chain.Insert(0, e.P1); edges.RemoveAt(i); grew = true;
                        }
                    }
                }
                List<double[]> ring = new List<double[]>();
                for (int k = 0; k < chain.Count; k++)
                {
                    Node n = FindNodeById(chain[k]);
                    if (n != null) { ring.Add(new double[] { n.X, n.Y }); }
                }
                if (ring.Count > 2) { rings.Add(ring); }
            }
            // Sort: largest area first (outer), rest holes
            rings.Sort(delegate (List<double[]> a, List<double[]> b)
            {
                return Math.Abs(PolygonArea(b)).CompareTo(Math.Abs(PolygonArea(a)));
            });
            return rings;
        }

        private List<double[]> BuildOuterPolygon()
        {
            List<List<double[]>> rings = BuildPolygons();
            if (rings.Count == 0) return null;
            return rings[0];
        }

        private List<List<double[]>> BuildHolePolygons()
        {
            List<List<double[]>> rings = BuildPolygons();
            List<List<double[]>> holes = new List<List<double[]>>();
            for (int i = 1; i < rings.Count; i++) { holes.Add(rings[i]); }
            return holes;
        }

        private double PolygonArea(List<double[]> pts)
        {
            double a = 0;
            for (int i = 0; i < pts.Count; i++)
            {
                int j = (i + 1) % pts.Count;
                a += pts[i][0] * pts[j][1] - pts[j][0] * pts[i][1];
            }
            return 0.5 * a;
        }

        /// <summary>Ray-casting point-in-polygon test.</summary>
        public bool PointInPolygon(double x, double y, List<double[]> pts)
        {
            if (pts == null || pts.Count < 3) return false;
            bool inside = false;
            for (int i = 0, j = pts.Count - 1; i < pts.Count; j = i++)
            {
                double xi = pts[i][0], yi = pts[i][1];
                double xj = pts[j][0], yj = pts[j][1];
                bool intersect = ((yi > y) != (yj > y))
                    && (x < (xj - xi) * (y - yi) / (yj - yi + 1.0e-20) + xi);
                if (intersect) inside = !inside;
            }
            return inside;
        }
    }
}
