// ContourRenderer.cs - GDI+ Jet 色阶云图绘制（网格线框 / 应力 / 位移 / 速度 / 压力 / 变形图）
using System;
using System.Drawing;
using System.Drawing.Drawing2D;
using System.Drawing.Imaging;
using FEM2D.Elements;
using FEM2D.Mesh;

namespace FEM2D.Rendering
{
    public static class ContourRenderer
    {
        private enum FieldMode
        { Mesh, Displacement, Deformed, StressXX, StressYY, StressXY, VonMises, Velocity, Pressure }

        public static Bitmap Render(FEMesh mesh, int width, int height, string field, bool showDeformed)
        {
            FieldMode fm = ParseField(field);
            Bitmap bmp = new Bitmap(Math.Max(width, 100), Math.Max(height, 100));
            using (Graphics g = Graphics.FromImage(bmp))
            {
                g.SmoothingMode = SmoothingMode.AntiAlias;
                g.Clear(Color.White);
                double xmin = 1e30, xmax = -1e30, ymin = 1e30, ymax = -1e30;
                for (int i = 0; i < mesh.NumNodes; i++)
                {
                    double x = mesh.Coords[i, 0], y = mesh.Coords[i, 1];
                    if (x < xmin) xmin = x; if (x > xmax) xmax = x;
                    if (y < ymin) ymin = y; if (y > ymax) ymax = y;
                }
                double defScale = 0.0;
                if (showDeformed || fm == FieldMode.Deformed)
                {
                    double du = 0.0, dv = 0.0;
                    for (int i = 0; i < mesh.NumNodes; i++)
                    {
                        if (Math.Abs(mesh.DisplacementU[i]) > du) du = Math.Abs(mesh.DisplacementU[i]);
                        if (Math.Abs(mesh.DisplacementV[i]) > dv) dv = Math.Abs(mesh.DisplacementV[i]);
                    }
                    double maxDisp = Math.Max(du, dv);
                    if (maxDisp > 1e-30) defScale = 0.1 * Math.Max(xmax - xmin, ymax - ymin) / maxDisp;
                }

                double dx = xmax - xmin, dy = ymax - ymin;
                if (dx < 1e-12) dx = 1.0;
                double asp = (double)bmp.Width / bmp.Height;
                int marginLeft = 60, marginRight = 80, marginTop = 40, marginBottom = 40;
                int pw = bmp.Width - marginLeft - marginRight;
                int ph = bmp.Height - marginTop - marginBottom;
                double sc = Math.Min(pw / dx, ph / dy);
                double W = dx * sc, H = dy * sc;
                double ox = marginLeft + (pw - W) * 0.5 - xmin * sc;
                double oy = marginTop + (ph - H) * 0.5 + ymax * sc;

                double vmin = 0.0, vmax = 0.0;
                GetFieldRange(mesh, fm, out vmin, out vmax);
                if (vmax - vmin < 1e-20) { vmin -= 1.0; vmax += 1.0; }

                int nen = FEMesh.NodesPerElement(mesh.ElementType);
                double[,] xDef = new double[mesh.NumNodes, 2];
                for (int i = 0; i < mesh.NumNodes; i++)
                {
                    xDef[i, 0] = mesh.Coords[i, 0] + defScale * mesh.DisplacementU[i];
                    xDef[i, 1] = mesh.Coords[i, 1] + defScale * mesh.DisplacementV[i];
                }
                double[,] nodeField = null;
                if (fm != FieldMode.Mesh && IsNodalField(fm))
                    nodeField = BuildNodalField(mesh, fm);

                for (int e = 0; e < mesh.NumElements; e++)
                {
                    PointF[] pts = new PointF[nen];
                    double[] fn = new double[nen];
                    for (int k = 0; k < nen; k++)
                    {
                        int n = mesh.Connectivity[e, k];
                        double xx, yy;
                        if (showDeformed || fm == FieldMode.Deformed) { xx = xDef[n, 0]; yy = xDef[n, 1]; }
                        else { xx = mesh.Coords[n, 0]; yy = mesh.Coords[n, 1]; }
                        pts[k] = new PointF((float)(ox + xx * sc), (float)(oy - yy * sc));
                        if (nodeField != null) fn[k] = nodeField[n, 0];
                    }
                    if (fm != FieldMode.Mesh)
                    {
                        double v;
                        if (nodeField != null) { v = 0.0; for (int k = 0; k < nen; k++) v += fn[k]; v /= nen; }
                        else v = CellFieldValue(mesh, e, fm);
                        double t = (v - vmin) / (vmax - vmin);
                        Color c = JetColor(t);
                        using (SolidBrush brush = new SolidBrush(c))
                        {
                            FillPoly(g, brush, pts, nen);
                        }
                    }
                    using (Pen pen = new Pen((fm == FieldMode.Mesh) ? Color.Black : Color.FromArgb(80, 40, 40, 40), 0.5f))
                    {
                        if (fm == FieldMode.Deformed)
                        {
                            PointF[] pts0 = new PointF[nen];
                            for (int k = 0; k < nen; k++)
                            {
                                int n = mesh.Connectivity[e, k];
                                pts0[k] = new PointF((float)(ox + mesh.Coords[n, 0] * sc), (float)(oy - mesh.Coords[n, 1] * sc));
                            }
                            using (Pen dashed = new Pen(Color.Gray, 0.8f) { DashStyle = DashStyle.Dash })
                                DrawPoly(g, dashed, pts0, nen);
                        }
                        DrawPoly(g, pen, pts, nen);
                    }
                }

                using (Font font = new Font("Consolas", 9.0f))
                {
                    DrawColorBar(g, bmp.Width - marginRight + 15, marginTop, 30, ph, vmin, vmax, font);
                    g.DrawString("FEM2D  " + fm.ToString() + "  (" + mesh.ElementType + ", N=" + mesh.NumNodes + ", E=" + mesh.NumElements + ")",
                        new Font("Consolas", 11, FontStyle.Bold), Brushes.Black, marginLeft, 8);
                }
            }
            return bmp;
        }

        public static void RenderToFile(FEMesh mesh, string path, int w, int h, string field, bool showDef)
        {
            using (Bitmap bmp = Render(mesh, w, h, field, showDef))
                bmp.Save(path, ImageFormat.Png);
        }

        private static FieldMode ParseField(string f)
        {
            string s = (f ?? "").ToLowerInvariant();
            if (s.Contains("mesh") || s.Contains("网格")) return FieldMode.Mesh;
            if (s.Contains("deform") || s.Contains("变形")) return FieldMode.Deformed;
            if (s.Contains("disp") || s.Contains("位移")) return FieldMode.Displacement;
            if (s.Contains("sxx") || s.Contains("σxx")) return FieldMode.StressXX;
            if (s.Contains("syy") || s.Contains("σyy")) return FieldMode.StressYY;
            if (s.Contains("sxy") || s.Contains("txy") || s.Contains("τxy")) return FieldMode.StressXY;
            if (s.Contains("von") || s.Contains("mises")) return FieldMode.VonMises;
            if (s.Contains("vel") || s.Contains("速度")) return FieldMode.Velocity;
            if (s.Contains("pres") || s.Contains("压力")) return FieldMode.Pressure;
            return FieldMode.Mesh;
        }

        private static bool IsNodalField(FieldMode fm)
        {
            return fm == FieldMode.Displacement || fm == FieldMode.Deformed || fm == FieldMode.Velocity || fm == FieldMode.Pressure;
        }

        private static double[,] BuildNodalField(FEMesh mesh, FieldMode fm)
        {
            int nN = mesh.NumNodes;
            double[,] f = new double[nN, 1];
            int[] cnt = new int[nN];
            for (int i = 0; i < nN; i++) f[i, 0] = 0.0;
            if (fm == FieldMode.Pressure) { for (int i = 0; i < nN; i++) f[i, 0] = mesh.Pressure[i]; return f; }
            int nen = FEMesh.NodesPerElement(mesh.ElementType);
            for (int e = 0; e < mesh.NumElements; e++)
            {
                double v = CellFieldValue(mesh, e, fm);
                for (int k = 0; k < nen; k++)
                {
                    int n = mesh.Connectivity[e, k];
                    f[n, 0] += v; cnt[n]++;
                }
            }
            for (int i = 0; i < nN; i++) if (cnt[i] > 0) f[i, 0] /= cnt[i];
            return f;
        }

        private static void GetFieldRange(FEMesh mesh, FieldMode fm, out double vmin, out double vmax)
        {
            vmin = double.MaxValue; vmax = double.MinValue;
            if (fm == FieldMode.Mesh || fm == FieldMode.Deformed) { vmin = 0.0; vmax = 1.0; return; }
            if (fm == FieldMode.Pressure)
            {
                for (int i = 0; i < mesh.NumNodes; i++) { double v = mesh.Pressure[i]; if (v < vmin) vmin = v; if (v > vmax) vmax = v; }
                return;
            }
            if (fm == FieldMode.Displacement)
            {
                for (int i = 0; i < mesh.NumNodes; i++)
                {
                    double mag = Math.Sqrt(mesh.DisplacementU[i] * mesh.DisplacementU[i] + mesh.DisplacementV[i] * mesh.DisplacementV[i]);
                    if (mag < vmin) vmin = mag; if (mag > vmax) vmax = mag;
                }
                return;
            }
            if (fm == FieldMode.Velocity)
            {
                for (int i = 0; i < mesh.NumNodes; i++)
                {
                    double mag = Math.Sqrt(mesh.VelocityU[i] * mesh.VelocityU[i] + mesh.VelocityV[i] * mesh.VelocityV[i]);
                    if (mag < vmin) vmin = mag; if (mag > vmax) vmax = mag;
                }
                return;
            }
            for (int e = 0; e < mesh.NumElements; e++)
            {
                double v = CellFieldValue(mesh, e, fm);
                if (v < vmin) vmin = v; if (v > vmax) vmax = v;
            }
        }

        private static double CellFieldValue(FEMesh mesh, int e, FieldMode fm)
        {
            int nen = FEMesh.NodesPerElement(mesh.ElementType);
            switch (fm)
            {
                case FieldMode.StressXX: return mesh.StressXX[e];
                case FieldMode.StressYY: return mesh.StressYY[e];
                case FieldMode.StressXY: return mesh.StressXY[e];
                case FieldMode.VonMises: return mesh.VonMises[e];
                case FieldMode.Deformed:
                case FieldMode.Displacement:
                    { double s = 0.0; for (int k = 0; k < nen; k++) { int n = mesh.Connectivity[e, k]; s += Math.Sqrt(mesh.DisplacementU[n] * mesh.DisplacementU[n] + mesh.DisplacementV[n] * mesh.DisplacementV[n]); } return s / nen; }
                case FieldMode.Velocity:
                    { double s = 0.0; for (int k = 0; k < nen; k++) { int n = mesh.Connectivity[e, k]; s += Math.Sqrt(mesh.VelocityU[n] * mesh.VelocityU[n] + mesh.VelocityV[n] * mesh.VelocityV[n]); } return s / nen; }
                case FieldMode.Pressure:
                    { double s = 0.0; for (int k = 0; k < nen; k++) s += mesh.Pressure[mesh.Connectivity[e, k]]; return s / nen; }
                default: return 0.0;
            }
        }

        private static void FillPoly(Graphics g, Brush b, PointF[] pts, int n)
        {
            if (n == 3) g.FillPolygon(b, new PointF[] { pts[0], pts[1], pts[2] });
            else if (n == 4) g.FillPolygon(b, new PointF[] { pts[0], pts[1], pts[2], pts[3] });
            else if (n == 6) g.FillPolygon(b, new PointF[] { pts[0], pts[3], pts[1], pts[4], pts[2], pts[5] });
        }

        private static void DrawPoly(Graphics g, Pen p, PointF[] pts, int n)
        {
            if (n == 3) g.DrawPolygon(p, new PointF[] { pts[0], pts[1], pts[2] });
            else if (n == 4) g.DrawPolygon(p, new PointF[] { pts[0], pts[1], pts[2], pts[3] });
            else if (n == 6)
            {
                g.DrawLines(p, new PointF[] { pts[0], pts[3], pts[1], pts[4], pts[2], pts[5], pts[0] });
            }
        }

        private static void DrawColorBar(Graphics g, int x, int y, int w, int h, double vmin, double vmax, Font font)
        {
            int steps = 64;
            for (int i = 0; i < steps; i++)
            {
                double t = (double)i / (steps - 1);
                Color c = JetColor(t);
                using (SolidBrush b = new SolidBrush(c))
                    g.FillRectangle(b, x, y + h - (int)(t * h), w, (int)(h / (double)steps) + 1);
            }
            g.DrawRectangle(Pens.Black, x, y, w, h);
            StringFormat fmt = new StringFormat();
            fmt.Alignment = StringAlignment.Near;
            for (int k = 0; k <= 5; k++)
            {
                double t = k / 5.0;
                double v = vmin + t * (vmax - vmin);
                int yy = y + h - (int)(t * h);
                g.DrawLine(Pens.Black, x + w, yy, x + w + 5, yy);
                g.DrawString(v.ToString("G3", System.Globalization.CultureInfo.InvariantCulture),
                    font, Brushes.Black, x + w + 8, yy - 7, fmt);
            }
        }

        private static Color JetColor(double t)
        {
            if (t < 0.0) t = 0.0; if (t > 1.0) t = 1.0;
            double r, g, b;
            if (t < 0.125) { r = 0; g = 0; b = 0.5 + 4 * t; }
            else if (t < 0.375) { r = 0; g = 4 * (t - 0.125); b = 1; }
            else if (t < 0.625) { r = 4 * (t - 0.375); g = 1; b = 1 - 4 * (t - 0.375); }
            else if (t < 0.875) { r = 1; g = 1 - 4 * (t - 0.625); b = 0; }
            else { r = 1 - 4 * (t - 0.875); g = 0; b = 0; }
            return Color.FromArgb((int)(r * 255), (int)(g * 255), (int)(b * 255));
        }
    }
}
