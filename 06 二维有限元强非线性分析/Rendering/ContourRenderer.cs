using System;
using System.Drawing;
using System.Drawing.Drawing2D;
using NonlinearFEM2D.Mesh;
using NonlinearFEM2D.Elements;

namespace NonlinearFEM2D.Rendering
{
    public sealed class ContourRenderer
    {
        public Bitmap Render(FEMesh mesh, string field, int width, int height, int modeIdx)
        {
            double xmin = 1e99, xmax = -1e99, ymin = 1e99, ymax = -1e99;
            for (int i = 0; i < mesh.NumNodes; i++)
            {
                if (mesh.X[i] < xmin) xmin = mesh.X[i];
                if (mesh.X[i] > xmax) xmax = mesh.X[i];
                if (mesh.Y[i] < ymin) ymin = mesh.Y[i];
                if (mesh.Y[i] > ymax) ymax = mesh.Y[i];
            }
            double pad = Math.Max((xmax - xmin), (ymax - ymin)) * 0.08;
            double W = xmax - xmin, H = ymax - ymin;
            double scale = Math.Min((width - 180) / W, (height - 40) / H);
            double ox = 90 - (xmin + xmax) / 2 * scale + W / 2 * scale;
            double oy = 20 - (ymin + ymax) / 2 * scale + H / 2 * scale;
            Bitmap bmp = new Bitmap(width, height);
            using (Graphics g = Graphics.FromImage(bmp))
            {
                g.SmoothingMode = SmoothingMode.AntiAlias;
                g.Clear(Color.White);
                // 提取场值
                double[] v = new double[mesh.NumNodes];
                double vmin = 1e99, vmax = -1e99;
                bool meshOnly = false, deform = false;
                double scaleDef = 1.0;
                string label = field;
                string fn = field;
                if (fn.Contains("网格") || fn == "Mesh") meshOnly = true;
                else if (fn.Contains("变形") || fn.Contains("Deform")) deform = true;
                if (fn.Contains("位移") || fn.Contains("Displacement"))
                {
                    for (int i = 0; i < mesh.NumNodes; i++)
                        v[i] = Math.Sqrt(mesh.DisplacementU[i] * mesh.DisplacementU[i] + mesh.DisplacementV[i] * mesh.DisplacementV[i]);
                    label = "位移大小";
                }
                else if (fn.Contains("σxx") || fn.Contains("sxx") || fn.Contains("Sxx")) { for (int i = 0; i < mesh.NumNodes; i++) v[i] = mesh.StressXX[i]; label = "σxx"; }
                else if (fn.Contains("σyy") || fn.Contains("syy") || fn.Contains("Syy")) { for (int i = 0; i < mesh.NumNodes; i++) v[i] = mesh.StressYY[i]; label = "σyy"; }
                else if (fn.Contains("τxy") || fn.Contains("τ") || fn.Contains("sxy")) { for (int i = 0; i < mesh.NumNodes; i++) v[i] = mesh.StressXY[i]; label = "τxy"; }
                else if (fn.Contains("Von") || fn.Contains("Mises")) { for (int i = 0; i < mesh.NumNodes; i++) v[i] = mesh.VonMises[i]; label = "Von Mises"; }
                else if (fn.Contains("εp") || fn.Contains("Plas")) { for (int i = 0; i < mesh.NumNodes; i++) v[i] = mesh.EqPlasticStrain[i]; label = "等效塑性应变"; }
                else if (fn.Contains("速度") || fn.Contains("Velocity") || fn.Contains("|v|"))
                {
                    for (int i = 0; i < mesh.NumNodes; i++)
                        v[i] = Math.Sqrt(mesh.VelocityU[i] * mesh.VelocityU[i] + mesh.VelocityV[i] * mesh.VelocityV[i]);
                    label = "速度大小";
                }
                else if (fn.Contains("压力") || fn.Contains("Pressure")) { for (int i = 0; i < mesh.NumNodes; i++) v[i] = mesh.Pressure[i]; label = "压力"; }
                else
                {
                    meshOnly = true;
                }
                for (int i = 0; i < mesh.NumNodes; i++)
                {
                    if (double.IsNaN(v[i])) v[i] = 0;
                    if (v[i] < vmin) vmin = v[i];
                    if (v[i] > vmax) vmax = v[i];
                }
                if (Math.Abs(vmax - vmin) < 1.0e-20) { vmax = vmin + 1.0; }
                if (deform)
                {
                    double maxD = 0;
                    for (int i = 0; i < mesh.NumNodes; i++)
                    {
                        double d = Math.Sqrt(mesh.DisplacementU[i] * mesh.DisplacementU[i] + mesh.DisplacementV[i] * mesh.DisplacementV[i]);
                        if (d > maxD) maxD = d;
                    }
                    if (maxD > 1.0e-12) scaleDef = 0.1 * Math.Max(W, H) / maxD;
                }
                // 填色
                foreach (FiniteElement e in mesh.Elements)
                {
                    int npe = e.NodeIds.Length;
                    PointF[] pts = new PointF[npe];
                    PointF[] pts0 = new PointF[npe];
                    Color[] cols = new Color[npe];
                    for (int k = 0; k < npe; k++)
                    {
                        int nid = e.NodeIds[k];
                        double xv = mesh.X[nid] + (deform ? mesh.DisplacementU[nid] * scaleDef : 0);
                        double yv = mesh.Y[nid] + (deform ? mesh.DisplacementV[nid] * scaleDef : 0);
                        pts[k] = new PointF((float)(xv * scale + ox), (float)(height - (yv * scale + oy)));
                        pts0[k] = new PointF((float)(mesh.X[nid] * scale + ox), (float)(height - (mesh.Y[nid] * scale + oy)));
                        double t = (v[nid] - vmin) / (vmax - vmin);
                        if (t < 0) t = 0; if (t > 1) t = 1;
                        cols[k] = JetColor(t);
                    }
                    if (!meshOnly)
                    {
                        // 简化：用中心单色填充（CST3），对 LT6/Q4 用三角剖分
                        FillElement(g, e, pts, cols);
                    }
                    using (Pen pen = new Pen(meshOnly ? Color.Black : Color.DarkGray, 1.0f))
                    {
                        g.DrawPolygon(pen, pts);
                    }
                    if (deform)
                    {
                        using (Pen p0 = new Pen(Color.Gray, 1.0f) { DashStyle = DashStyle.Dash })
                        {
                            g.DrawPolygon(p0, pts0);
                        }
                    }
                }
                // 标题
                using (Font f = new Font("微软雅黑", 12, FontStyle.Bold))
                {
                    g.DrawString(label + "  范围: [" + vmin.ToString("E3") + ", " + vmax.ToString("E3") + "]",
                        f, Brushes.Black, 20, 5);
                }
                // 色条
                DrawColorBar(g, width - 150, 40, 30, height - 80, vmin, vmax);
            }
            return bmp;
        }

        private void FillElement(Graphics g, FiniteElement e, PointF[] pts, Color[] cols)
        {
            int npe = pts.Length;
            if (npe == 3)
            {
                using (GraphicsPath gp = new GraphicsPath())
                {
                    gp.AddPolygon(pts);
                    PathGradientBrush pgb = new PathGradientBrush(gp);
                    pgb.CenterColor = Blend3(cols[0], cols[1], cols[2]);
                    pgb.SurroundColors = cols;
                    g.FillPath(pgb, gp);
                    pgb.Dispose();
                }
            }
            else if (npe == 4)
            {
                FillTriangle(g, pts[0], pts[1], pts[2], cols[0], cols[1], cols[2]);
                FillTriangle(g, pts[0], pts[2], pts[3], cols[0], cols[2], cols[3]);
            }
            else if (npe == 6)
            {
                FillTriangle(g, pts[0], pts[3], pts[5], cols[0], cols[3], cols[5]);
                FillTriangle(g, pts[3], pts[1], pts[4], cols[3], cols[1], cols[4]);
                FillTriangle(g, pts[5], pts[4], pts[2], cols[5], cols[4], cols[2]);
                FillTriangle(g, pts[3], pts[4], pts[5], cols[3], cols[4], cols[5]);
            }
        }
        private void FillTriangle(Graphics g, PointF a, PointF b, PointF c, Color ca, Color cb, Color cc)
        {
            PointF[] tri = { a, b, c };
            Color[] cols = { ca, cb, cc };
            using (GraphicsPath gp = new GraphicsPath())
            {
                gp.AddPolygon(tri);
                PathGradientBrush pgb = new PathGradientBrush(gp);
                pgb.CenterColor = Blend3(ca, cb, cc);
                pgb.SurroundColors = cols;
                g.FillPath(pgb, gp);
                pgb.Dispose();
            }
        }
        private Color Blend3(Color a, Color b, Color c)
        {
            return Color.FromArgb((a.R + b.R + c.R) / 3, (a.G + b.G + c.G) / 3, (a.B + c.B + c.B) / 3);
        }
        private Color JetColor(double t)
        {
            double r, g, b;
            if (t < 0.25) { r = 0; g = 4 * t; b = 1; }
            else if (t < 0.5) { r = 0; g = 1; b = 1 - 4 * (t - 0.25); }
            else if (t < 0.75) { r = 4 * (t - 0.5); g = 1; b = 0; }
            else { r = 1; g = 1 - 4 * (t - 0.75); b = 0; }
            int R = (int)Math.Max(0, Math.Min(255, r * 255));
            int G = (int)Math.Max(0, Math.Min(255, g * 255));
            int B = (int)Math.Max(0, Math.Min(255, b * 255));
            return Color.FromArgb(R, G, B);
        }
        private void DrawColorBar(Graphics g, int x, int y, int w, int h, double vmin, double vmax)
        {
            for (int i = 0; i < h; i++)
            {
                double t = 1.0 - (double)i / h;
                using (Pen p = new Pen(JetColor(t)))
                {
                    g.DrawLine(p, x, y + i, x + w, y + i);
                }
            }
            using (Pen p = new Pen(Color.Black))
            {
                g.DrawRectangle(p, x, y, w, h);
            }
            using (Font f = new Font("Consolas", 9))
            {
                for (int i = 0; i <= 6; i++)
                {
                    double t = (double)i / 6;
                    double val = vmin + t * (vmax - vmin);
                    int yy = y + h - (int)(t * h);
                    g.DrawString(val.ToString("E2"), f, Brushes.Black, x + w + 4, yy - 7);
                    g.DrawLine(Pens.Black, x + w, yy, x + w + 3, yy);
                }
            }
        }
    }
}
