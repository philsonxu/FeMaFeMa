using System;
using System.Collections.Generic;
using System.Drawing;
using System.Drawing.Drawing2D;
using System.Drawing.Imaging;
using System.Globalization;
using System.IO;
using System.Text;
using Fem2DFluid.Models;

namespace Fem2DFluid.Services
{
    /// <summary>
    /// Export results to CSV, Tecplot ASCII, plain text, and color PNG cloud-map
    /// with Jet colormap, mesh lines, velocity vectors, legend/colorbar and title.
    /// </summary>
    public class ResultExporter
    {
        private FemModel _model;

        public ResultExporter(FemModel model)
        {
            _model = model;
        }

        public void ExportCsv(string path)
        {
            using (StreamWriter sw = new StreamWriter(path, false, Encoding.UTF8))
            {
                sw.WriteLine("id,x,y,phi,u,v,p,vx,vy,vmag,tag,bctype");
                for (int i = 0; i < _model.Nodes.Count; i++)
                {
                    Node n = _model.Nodes[i];
                    double vm = Math.Sqrt(n.Vx * n.Vx + n.Vy * n.Vy);
                    sw.WriteLine("{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11}",
                        n.Id,
                        n.X.ToString(CultureInfo.InvariantCulture),
                        n.Y.ToString(CultureInfo.InvariantCulture),
                        n.Phi.ToString(CultureInfo.InvariantCulture),
                        n.U.ToString(CultureInfo.InvariantCulture),
                        n.V.ToString(CultureInfo.InvariantCulture),
                        n.P.ToString(CultureInfo.InvariantCulture),
                        n.Vx.ToString(CultureInfo.InvariantCulture),
                        n.Vy.ToString(CultureInfo.InvariantCulture),
                        vm.ToString(CultureInfo.InvariantCulture),
                        n.Tag, n.BCType);
                }
            }
        }

        public void ExportTecplot(string path)
        {
            using (StreamWriter sw = new StreamWriter(path, false, Encoding.ASCII))
            {
                sw.WriteLine("TITLE = \"Fem2DFluid result\"");
                sw.WriteLine("VARIABLES = \"X\", \"Y\", \"PHI\", \"U\", \"V\", \"P\"");
                sw.WriteLine("ZONE N={0}, E={1}, F=FEPOINT, ET=TRIANGLE",
                    _model.Nodes.Count, _model.Elements.Count);
                for (int i = 0; i < _model.Nodes.Count; i++)
                {
                    Node n = _model.Nodes[i];
                    sw.WriteLine("{0} {1} {2} {3} {4} {5}",
                        n.X.ToString(CultureInfo.InvariantCulture),
                        n.Y.ToString(CultureInfo.InvariantCulture),
                        n.Phi.ToString("E6", CultureInfo.InvariantCulture),
                        n.U.ToString("E6", CultureInfo.InvariantCulture),
                        n.V.ToString("E6", CultureInfo.InvariantCulture),
                        n.P.ToString("E6", CultureInfo.InvariantCulture));
                }
                for (int i = 0; i < _model.Elements.Count; i++)
                {
                    Element e = _model.Elements[i];
                    sw.WriteLine("{0} {1} {2}", e.A + 1, e.B + 1, e.C + 1);
                }
            }
        }

        public void ExportText(string path)
        {
            using (StreamWriter sw = new StreamWriter(path, false, Encoding.UTF8))
            {
                sw.WriteLine("=== Fem2DFluid Result ===");
                sw.WriteLine("Model: {0}", _model.Name ?? "");
                sw.WriteLine("Nodes: {0}  Elements: {1}  Boundary edges: {2}",
                    _model.Nodes.Count, _model.Elements.Count, _model.BoundaryEdgeCache.Count);
                sw.WriteLine("Bounds: X[{0:F4},{1:F4}] Y[{2:F4},{3:F4}]",
                    _model.Xmin, _model.Xmax, _model.Ymin, _model.Ymax);
                sw.WriteLine();
                sw.WriteLine("--- Node results ---");
                sw.WriteLine("ID\tX\tY\tPhi\tU\tV\tP\t|V|\tTag");
                for (int i = 0; i < _model.Nodes.Count; i++)
                {
                    Node n = _model.Nodes[i];
                    double vm = Math.Sqrt(n.Vx * n.Vx + n.Vy * n.Vy);
                    sw.WriteLine("{0}\t{1:F4}\t{2:F4}\t{3:F4}\t{4:F4}\t{5:F4}\t{6:F4}\t{7:F4}\t{8}",
                        n.Id, n.X, n.Y, n.Phi, n.U, n.V, n.P, vm, n.Tag);
                }
                sw.WriteLine();
                sw.WriteLine("--- Element results ---");
                sw.WriteLine("ID\tA\tB\tC\tArea\tVx\tVy\tVmag\tP");
                for (int i = 0; i < _model.Elements.Count; i++)
                {
                    Element e = _model.Elements[i];
                    sw.WriteLine("{0}\t{1}\t{2}\t{3}\t{4:F4}\t{5:F4}\t{6:F4}\t{7:F4}\t{8:F4}",
                        e.Id, e.A, e.B, e.C, e.Area, e.Vx, e.Vy, e.Vmag, e.P);
                }
            }
        }

        /// <summary>
        /// Render a color cloud-map PNG with Jet colormap, optional mesh lines,
        /// velocity vectors, right-side colorbar, and title. For LT6 elements the
        /// quadratic triangle is split into four sub-triangles for smoother rendering.
        /// </summary>
        public void SavePng(string path, string field, int width, int height,
            bool showMesh, bool showVectors)
        {
            if (width < 200) width = 800;
            if (height < 200) height = 600;
            int cbWidth = 60;
            int margin = 40;
            int plotW = width - cbWidth - margin * 2;
            int plotH = height - margin * 2 - 30;
            using (Bitmap bmp = new Bitmap(width, height))
            using (Graphics g = Graphics.FromImage(bmp))
            {
                g.SmoothingMode = SmoothingMode.AntiAlias;
                g.Clear(Color.White);
                _model.ComputeBounds();
                double xmin = _model.Xmin, xmax = _model.Xmax;
                double ymin = _model.Ymin, ymax = _model.Ymax;
                double sx = plotW / Math.Max(xmax - xmin, 1.0e-10);
                double sy = plotH / Math.Max(ymax - ymin, 1.0e-10);
                double s = Math.Min(sx, sy);
                double ox = margin + (plotW - (xmax - xmin) * s) * 0.5;
                double oy = margin + (plotH - (ymax - ymin) * s) * 0.5;
                // Compute field range
                double vmin = double.MaxValue, vmax = double.MinValue;
                for (int i = 0; i < _model.Nodes.Count; i++)
                {
                    double v = GetField(_model.Nodes[i], field);
                    if (v < vmin) vmin = v;
                    if (v > vmax) vmax = v;
                }
                if (vmax - vmin < 1.0e-12) { vmax = vmin + 1.0; }
                // Draw elements as filled triangles (sub-triangulate LT6)
                for (int i = 0; i < _model.Elements.Count; i++)
                {
                    Element e = _model.Elements[i];
                    DrawElement(e, field, vmin, vmax, g, s, ox, oy, xmin, ymax);
                }
                if (showMesh)
                {
                    using (Pen p = new Pen(Color.FromArgb(80, Color.Black), 0.5f))
                    {
                        for (int i = 0; i < _model.Elements.Count; i++)
                        {
                            Element e = _model.Elements[i];
                            PointF[] pts = TriPts(e, s, ox, oy, xmin, ymax);
                            g.DrawPolygon(p, pts);
                        }
                    }
                }
                if (showVectors)
                {
                    DrawVelocityVectors(g, s, ox, oy, xmin, ymax, vmax - vmin);
                }
                // Colorbar
                DrawColorBar(g, margin + plotW + 10, margin, cbWidth - 10, plotH, vmin, vmax, field);
                // Title
                using (Font f = new Font("Arial", 12, FontStyle.Bold))
                {
                    string title = string.Format("{0} - {1}", _model.Name ?? "FEM", field);
                    g.DrawString(title, f, Brushes.Black, margin, 8);
                }
                // Frame
                g.DrawRectangle(Pens.Black, (float)ox, (float)oy, (float)plotW, (float)plotH);
                bmp.Save(path, ImageFormat.Png);
            }
        }

        private PointF[] TriPts(Element e, double s, double ox, double oy, double xmin, double ymax)
        {
            Node na = _model.GetNode(e.A);
            Node nb = _model.GetNode(e.B);
            Node nc = _model.GetNode(e.C);
            PointF[] p = new PointF[3];
            p[0] = new PointF((float)(ox + (na.X - xmin) * s), (float)(oy + (ymax - na.Y) * s));
            p[1] = new PointF((float)(ox + (nb.X - xmin) * s), (float)(oy + (ymax - nb.Y) * s));
            p[2] = new PointF((float)(ox + (nc.X - xmin) * s), (float)(oy + (ymax - nc.Y) * s));
            return p;
        }

        private void DrawElement(Element e, string field, double vmin, double vmax,
            Graphics g, double s, double ox, double oy, double xmin, double ymax)
        {
            // Split LT6 into 4 sub-triangles for smoother quadratic visualisation
            if (e.Order == 2 && e.N4 >= 0 && e.N5 >= 0 && e.N6 >= 0)
            {
                SubTriangle[] subs = SubTri(e);
                for (int k = 0; k < subs.Length; k++)
                {
                    FillTriField(subs[k].na, subs[k].nb, subs[k].nc, field, vmin, vmax, g, s, ox, oy, xmin, ymax);
                }
            }
            else
            {
                Node na = _model.GetNode(e.A);
                Node nb = _model.GetNode(e.B);
                Node nc = _model.GetNode(e.C);
                FillTriField(na, nb, nc, field, vmin, vmax, g, s, ox, oy, xmin, ymax);
            }
        }

        private struct SubTriangle
        {
            public Node na, nb, nc;
            public SubTriangle(Node a, Node b, Node c) { na = a; nb = b; nc = c; }
        }

        private SubTriangle[] SubTri(Element e)
        {
            Node a = _model.GetNode(e.A), b = _model.GetNode(e.B), c = _model.GetNode(e.C);
            Node m4 = _model.GetNode(e.N4), m5 = _model.GetNode(e.N5), m6 = _model.GetNode(e.N6);
            return new SubTriangle[]
            {
                new SubTriangle(a, m6, m5),
                new SubTriangle(b, m4, m6),
                new SubTriangle(c, m5, m4),
                new SubTriangle(m4, m5, m6)
            };
        }

        private void FillTriField(Node na, Node nb, Node nc, string field,
            double vmin, double vmax, Graphics g, double s, double ox, double oy,
            double xmin, double ymax)
        {
            PointF[] pts = new PointF[3];
            pts[0] = new PointF((float)(ox + (na.X - xmin) * s), (float)(oy + (ymax - na.Y) * s));
            pts[1] = new PointF((float)(ox + (nb.X - xmin) * s), (float)(oy + (ymax - nb.Y) * s));
            pts[2] = new PointF((float)(ox + (nc.X - xmin) * s), (float)(oy + (ymax - nc.Y) * s));
            double va = GetField(na, field);
            double vb = GetField(nb, field);
            double vc = GetField(nc, field);
            double avg = (va + vb + vc) / 3.0;
            double t = (avg - vmin) / (vmax - vmin);
            if (t < 0) t = 0; if (t > 1) t = 1;
            Color col = JetColor(t);
            using (Brush br = new SolidBrush(col))
            {
                g.FillPolygon(br, pts);
            }
        }

        private double GetField(Node n, string field)
        {
            if (n == null) return 0.0;
            string f = (field ?? "phi").ToLowerInvariant();
            switch (f)
            {
                case "u": return n.U;
                case "v": return n.V;
                case "p": return n.P;
                case "vx": return n.Vx;
                case "vy": return n.Vy;
                case "vmag":
                    //case "v":
                    return Math.Sqrt(n.Vx * n.Vx + n.Vy * n.Vy);
                case "phi":
                default:
                    return n.Phi;
            }
        }

        /// <summary>Jet colormap: blue -> cyan -> green -> yellow -> red.</summary>
        public static Color JetColor(double t)
        {
            if (t < 0) t = 0; if (t > 1) t = 1;
            double r, g, b;
            if (t < 0.25)
            {
                r = 0; g = 4.0 * t; b = 1;
            }
            else if (t < 0.5)
            {
                r = 0; g = 1; b = 1 - 4.0 * (t - 0.25);
            }
            else if (t < 0.75)
            {
                r = 4.0 * (t - 0.5); g = 1; b = 0;
            }
            else
            {
                r = 1; g = 1 - 4.0 * (t - 0.75); b = 0;
            }
            int R = (int)(r * 255); if (R < 0) R = 0; if (R > 255) R = 255;
            int G = (int)(g * 255); if (G < 0) G = 0; if (G > 255) G = 255;
            int B = (int)(b * 255); if (B < 0) B = 0; if (B > 255) B = 255;
            return Color.FromArgb(255, R, G, B);
        }

        private void DrawVelocityVectors(Graphics g, double s, double ox, double oy,
            double xmin, double ymax, double scale0)
        {
            double maxv = 0;
            for (int i = 0; i < _model.Elements.Count; i++)
            {
                Element e = _model.Elements[i];
                if (e.Vmag > maxv) maxv = e.Vmag;
            }
            if (maxv < 1.0e-10) return;
            double refLen = (_model.Xmax - _model.Xmin) * 0.03;
            using (Pen p = new Pen(Color.Black, 0.8f))
            {
                p.EndCap = LineCap.ArrowAnchor;
                int step = Math.Max(1, _model.Elements.Count / 600);
                for (int i = 0; i < _model.Elements.Count; i += step)
                {
                    Element e = _model.Elements[i];
                    double cx = e.Barycenter[0], cy = e.Barycenter[1];
                    double x1 = ox + (cx - xmin) * s;
                    double y1 = oy + (ymax - cy) * s;
                    double len = e.Vmag / maxv * refLen * s;
                    double ang = Math.Atan2(e.Vy, e.Vx);
                    double x2 = x1 + Math.Cos(ang) * len;
                    double y2 = y1 - Math.Sin(ang) * len;
                    g.DrawLine(p, (float)x1, (float)y1, (float)x2, (float)y2);
                }
            }
        }

        private void DrawColorBar(Graphics g, int x, int y, int w, int h,
            double vmin, double vmax, string field)
        {
            int nBins = 64;
            for (int i = 0; i < nBins; i++)
            {
                double t = 1.0 - (double)i / (nBins - 1);
                Color c = JetColor(t);
                using (Brush br = new SolidBrush(c))
                {
                    int y0 = y + (int)((double)i * h / nBins);
                    int y1 = y + (int)((double)(i + 1) * h / nBins);
                    g.FillRectangle(br, x, y0, w, y1 - y0 + 1);
                }
            }
            g.DrawRectangle(Pens.Black, x, y, w, h);
            using (Font f = new Font("Arial", 8))
            {
                for (int k = 0; k <= 4; k++)
                {
                    double t = (double)k / 4;
                    double v = vmax - t * (vmax - vmin);
                    int yp = y + (int)(t * h);
                    g.DrawLine(Pens.Black, x + w, yp, x + w + 5, yp);
                    g.DrawString(v.ToString("G3"), f, Brushes.Black, x + w + 7, yp - 6);
                }
                using (Font ff = new Font("Arial", 9, FontStyle.Bold))
                {
                    g.DrawString(field ?? "", ff, Brushes.Black, x - 5, y - 18);
                }
            }
        }
    }
}
