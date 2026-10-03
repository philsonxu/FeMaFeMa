using System;
using System.Drawing;
using System.Drawing.Drawing2D;
using System.Drawing.Imaging;
using MultiPhysicsFEM2D.Elements;
using MultiPhysicsFEM2D.Mesh;

namespace MultiPhysicsFEM2D.Rendering
{
    public static class ContourRenderer
    {
        private static Color JetColor(double v)
        {
            // v ∈ [0,1] -> Jet
            v = Math.Max(0.0, Math.Min(1.0, v));
            double r, g, b;
            if (v < 0.125) { r = 0; g = 0; b = 0.5 + 4 * v; }
            else if (v < 0.375) { v -= 0.125; r = 0; g = 4 * v; b = 1.0; }
            else if (v < 0.625) { v -= 0.375; r = 4 * v; g = 1.0; b = 1.0 - 4 * v; }
            else if (v < 0.875) { v -= 0.625; r = 1.0; g = 1.0 - 4 * v; b = 0; }
            else { v -= 0.875; r = 1.0 - 4 * v; g = 0; b = 0; }
            return Color.FromArgb(255,
                (int)Math.Max(0, Math.Min(255, r * 255)),
                (int)Math.Max(0, Math.Min(255, g * 255)),
                (int)Math.Max(0, Math.Min(255, b * 255)));
        }

        public static void RenderToPng(FEMesh mesh, string path, int width, int height,
            string field, double deformScale)
        {
            using (Bitmap bmp = new Bitmap(width, height))
            using (Graphics g = Graphics.FromImage(bmp))
            {
                g.SmoothingMode = SmoothingMode.AntiAlias;
                g.Clear(Color.White);
                // 计算显示坐标与范围
                double xmin = double.MaxValue, xmax = double.MinValue;
                double ymin = double.MaxValue, ymax = double.MinValue;
                for (int i = 0; i < mesh.NumNodes; i++)
                {
                    if (mesh.X[i] < xmin) xmin = mesh.X[i];
                    if (mesh.X[i] > xmax) xmax = mesh.X[i];
                    if (mesh.Y[i] < ymin) ymin = mesh.Y[i];
                    if (mesh.Y[i] > ymax) ymax = mesh.Y[i];
                }
                // 变形范围
                double dxmax = 0, dymax = 0;
                for (int i = 0; i < mesh.NumNodes; i++)
                {
                    if (Math.Abs(mesh.DisplacementU[i]) > dxmax) dxmax = Math.Abs(mesh.DisplacementU[i]);
                    if (Math.Abs(mesh.DisplacementV[i]) > dymax) dymax = Math.Abs(mesh.DisplacementV[i]);
                }
                double scale = deformScale;
                double dpmax = Math.Max(dxmax, dymax);
                double L = Math.Max(xmax - xmin, ymax - ymin);
                if (field.Contains("变形") && dpmax > 1e-30) scale = deformScale * 0.1 * L / dpmax;
                int marginL = 80, marginR = 120, marginT = 40, marginB = 80;
                int drawW = width - marginL - marginR;
                int drawH = height - marginT - marginB;
                double sx = drawW / (xmax - xmin);
                double sy = drawH / (ymax - ymin);
                double s = Math.Min(sx, sy);
                double ox = marginL + (drawW - (xmax - xmin) * s) * 0.5;
                double oy = marginT + (drawH + (ymax - ymin) * s) * 0.5; // Y 轴翻转

                // 获取场
                double[] values; bool isDeformed = false;
                double vmin = 0, vmax = 1;
                ParseField(mesh, field, out values, out vmin, out vmax, out isDeformed);
                if (Math.Abs(vmax - vmin) < 1e-30) vmax = vmin + 1;

                // 逐单元填充彩色三角/四边形
                foreach (FiniteElement e in mesh.Elements)
                {
                    int nn = e.NodesPerElement;
                    PointF[] pts = new PointF[nn];
                    double[] vals = new double[nn];
                    for (int k = 0; k < nn; k++)
                    {
                        int nid = e.NodeIds[k];
                        double px = mesh.X[nid] + scale * mesh.DisplacementU[nid];
                        double py = mesh.Y[nid] + scale * mesh.DisplacementV[nid];
                        pts[k] = new PointF((float)(ox + (px - xmin) * s), (float)(oy - (py - ymin) * s));
                        vals[k] = values[nid];
                    }
                    // 简单按单元中心颜色填充（CST 正好是常数；Q4/LT6 切分三角）
                    if (nn == 3)
                    {
                        FillTriangle(g, pts, vals, vmin, vmax);
                    }
                    else if (nn == 4)
                    {
                        PointF[] t1 = new PointF[] { pts[0], pts[1], pts[2] };
                        PointF[] t2 = new PointF[] { pts[0], pts[2], pts[3] };
                        double[] v1 = new double[] { vals[0], vals[1], vals[2] };
                        double[] v2 = new double[] { vals[0], vals[2], vals[3] };
                        FillTriangle(g, t1, v1, vmin, vmax);
                        FillTriangle(g, t2, v2, vmin, vmax);
                    }
                    else if (nn == 6)
                    {
                        // LT6 拆成 4 个小三角（节点: 0-3-1 / 1-4-2 / 2-5-0 / 3-4-5）
                        int[][] tri = new int[][] {
                            new int[] {0,3,1}, new int[] {1,4,2}, new int[] {2,5,0}, new int[] {3,4,5}
                        };
                        foreach (int[] t in tri)
                        {
                            PointF[] tp = new PointF[] { pts[t[0]], pts[t[1]], pts[t[2]] };
                            double[] tv = new double[] { vals[t[0]], vals[t[1]], vals[t[2]] };
                            FillTriangle(g, tp, tv, vmin, vmax);
                        }
                    }
                }

                // 网格线（细线）
                using (Pen pGrid = new Pen(Color.FromArgb(48, Color.Black), 0.7f))
                {
                    foreach (FiniteElement e in mesh.Elements)
                    {
                        int nn = e.NodesPerElement;
                        PointF[] pts = new PointF[nn];
                        for (int k = 0; k < nn; k++)
                        {
                            int nid = e.NodeIds[k];
                            double px = mesh.X[nid] + scale * mesh.DisplacementU[nid];
                            double py = mesh.Y[nid] + scale * mesh.DisplacementV[nid];
                            pts[k] = new PointF((float)(ox + (px - xmin) * s), (float)(oy - (py - ymin) * s));
                        }
                        if (nn == 3) g.DrawPolygon(pGrid, pts);
                        else if (nn == 4) g.DrawPolygon(pGrid, pts);
                        else if (nn == 6)
                        {
                            g.DrawLine(pGrid, pts[0], pts[3]); g.DrawLine(pGrid, pts[3], pts[1]);
                            g.DrawLine(pGrid, pts[1], pts[4]); g.DrawLine(pGrid, pts[4], pts[2]);
                            g.DrawLine(pGrid, pts[2], pts[5]); g.DrawLine(pGrid, pts[5], pts[0]);
                            g.DrawLine(pGrid, pts[3], pts[4]); g.DrawLine(pGrid, pts[4], pts[5]); g.DrawLine(pGrid, pts[5], pts[3]);
                        }
                    }
                }

                // 变形前轮廓（灰色虚线）
                if (isDeformed)
                {
                    using (Pen pDef = new Pen(Color.Gray, 1.0f) { DashStyle = DashStyle.Dash })
                    {
                        foreach (FiniteElement e in mesh.Elements)
                        {
                            int nn = e.NodesPerElement;
                            PointF[] pts = new PointF[nn];
                            for (int k = 0; k < nn; k++)
                            {
                                int nid = e.NodeIds[k];
                                pts[k] = new PointF(
                                    (float)(ox + (mesh.X[nid] - xmin) * s),
                                    (float)(oy - (mesh.Y[nid] - ymin) * s));
                            }
                            if (nn == 3 || nn == 4) g.DrawPolygon(pDef, pts);
                        }
                    }
                }

                // 色条
                int barW = 24, barH = drawH - 40;
                int barX = width - marginR + 20;
                int barY = marginT + 20;
                using (Font ft = new Font("Consolas", 9f))
                using (SolidBrush fb = new SolidBrush(Color.Black))
                using (StringFormat sf = new StringFormat() { Alignment = StringAlignment.Near })
                {
                    int steps = 200;
                    for (int i = 0; i < steps; i++)
                    {
                        double v = 1.0 - (double)i / (steps - 1);
                        using (Pen pc = new Pen(JetColor(v), barW))
                            g.DrawLine(pc, barX, barY + barH * i / steps, barX + barW, barY + barH * i / steps);
                    }
                    g.DrawRectangle(Pens.Black, barX, barY, barW, barH);
                    // 刻度
                    for (int t = 0; t <= 5; t++)
                    {
                        double v = vmin + (vmax - vmin) * t / 5.0;
                        int yy = barY + barH * (5 - t) / 5;
                        g.DrawLine(Pens.Black, barX + barW, yy, barX + barW + 6, yy);
                        g.DrawString(string.Format("{0:G3}", v), ft, fb, barX + barW + 8, yy - 7, sf);
                    }
                    // 标题
                    using (Font fth = new Font("Arial", 11f, FontStyle.Bold))
                        g.DrawString(field, fth, fb, new PointF(marginL, 10));
                }

                bmp.Save(path, ImageFormat.Png);
            }
        }

        private static void FillTriangle(Graphics g, PointF[] pts, double[] vals, double vmin, double vmax)
        {
            double avg = (vals[0] + vals[1] + vals[2]) / 3.0;
            double t = (avg - vmin) / (vmax - vmin);
            using (SolidBrush br = new SolidBrush(JetColor(t)))
                g.FillPolygon(br, pts);
        }

        private static void ParseField(FEMesh mesh, string field,
            out double[] values, out double vmin, out double vmax, out bool isDeformed)
        {
            values = new double[mesh.NumNodes];
            vmin = double.MaxValue; vmax = double.MinValue;
            isDeformed = false;
            string f = field.ToLowerInvariant();
            if (f.Contains("网格") || f.Contains("mesh"))
            {
                for (int i = 0; i < mesh.NumNodes; i++) values[i] = 0;
                vmin = 0; vmax = 1;
                return;
            }
            double[] src = null;
            if (f.Contains("位移") || f.Contains("ux") || f.Contains("u")) { src = mesh.DisplacementU; isDeformed = true; }
            else if (f.Contains("uy") || f.Contains("vy")) { src = mesh.DisplacementV; isDeformed = true; }
            else if (f.Contains("变形") || f.Contains("deform"))
            {
                isDeformed = true;
                for (int i = 0; i < mesh.NumNodes; i++)
                    values[i] = Math.Sqrt(mesh.DisplacementU[i] * mesh.DisplacementU[i] + mesh.DisplacementV[i] * mesh.DisplacementV[i]);
            }
            else if (f.Contains("sxx") || f.Contains("σx") || f.Contains("应力 σxx")) { src = mesh.StressXX; }
            else if (f.Contains("syy") || f.Contains("σy")) { src = mesh.StressYY; }
            else if (f.Contains("sxy") || f.Contains("τxy") || f.Contains("τ")) { src = mesh.StressXY; }
            else if (f.Contains("von") || f.Contains("vm") || f.Contains("等效应力")) { src = mesh.VonMises; }
            else if (f.Contains("ep") || f.Contains("塑性")) { src = mesh.PlasticStrain; }
            else if (f.Contains("速度") || f.Contains("vel"))
            {
                for (int i = 0; i < mesh.NumNodes; i++)
                    values[i] = Math.Sqrt(mesh.VelocityU[i] * mesh.VelocityU[i] + mesh.VelocityV[i] * mesh.VelocityV[i]);
            }
            else if (f.Contains("vx")) src = mesh.VelocityU;
            else if (f.Contains("vy")) src = mesh.VelocityV;
            else if (f.Contains("压力") || f.Contains("press")) { src = mesh.Pressure; }
            else if (f.Contains("温度") || f.Contains("temp") || f.Contains("t")) { src = mesh.Temperature; }
            else src = mesh.VonMises;
            if (src != null) Array.Copy(src, values, mesh.NumNodes);
            if (f.Contains("变形") || f.Contains("速度")) { }
            else if (src == null) { }
            for (int i = 0; i < mesh.NumNodes; i++)
            {
                if (values[i] < vmin) vmin = values[i];
                if (values[i] > vmax) vmax = values[i];
            }
            if (vmin == vmax) { vmin -= 1e-12; vmax += 1e-12; }
        }
    }
}
