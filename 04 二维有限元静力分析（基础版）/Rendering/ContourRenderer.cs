using System;
using System.Drawing;
using System.Drawing.Drawing2D;
using FEM2D.Mesh;
using FEM2D.Solvers;

namespace FEM2D.Rendering
{
    public enum FieldType
    {
        Mesh,
        Displacement,
        StressXX,
        StressYY,
        StressXY,
        VonMises,
        Velocity,
        Pressure,
        Deformed
    }

    public static class ContourRenderer
    {
        public static Bitmap Render(FEMesh mesh, object result, FieldType field, int width, int height,
            double deformScale = 0.0, bool showDeformedOutline = false)
        {
            Bitmap bmp = new Bitmap(width, height);
            using (Graphics g = Graphics.FromImage(bmp))
            {
                g.SmoothingMode = SmoothingMode.AntiAlias;
                g.Clear(Color.White);

                // 计算坐标范围
                double xmin = double.PositiveInfinity, xmax = double.NegativeInfinity;
                double ymin = double.PositiveInfinity, ymax = double.NegativeInfinity;
                for (int i = 0; i < mesh.NumNodes; i++)
                {
                    double x = mesh.X(i); double y = mesh.Y(i);
                    if (x < xmin) xmin = x; if (x > xmax) xmax = x;
                    if (y < ymin) ymin = y; if (y > ymax) ymax = y;
                }
                double margin = Math.Max((xmax - xmin), (ymax - ymin)) * 0.08;
                xmin -= margin; xmax += margin; ymin -= margin; ymax += margin;
                double w = xmax - xmin, h = ymax - ymin;
                double scale = Math.Min((width - 60) / w, (height - 60) / h);
                double ox = (width - scale * w) * 0.5 - scale * xmin;
                double oy = (height - scale * h) * 0.5 + scale * ymax; // y 翻转
                Func<double, double, PointF> toPx = (x, y) => new PointF((float)(ox + scale * x), (float)(oy - scale * y));

                // 准备节点值（用于节点型字段）
                double[] nodeValues = null;
                double[] nodeDx = null, nodeDy = null;
                double vmin = 0, vmax = 0;
                if (result is StructuralResult sr)
                {
                    nodeDx = new double[mesh.NumNodes];
                    nodeDy = new double[mesh.NumNodes];
                    nodeValues = new double[mesh.NumNodes];
                    for (int i = 0; i < mesh.NumNodes; i++)
                    {
                        nodeDx[i] = sr.Displacement[2 * i];
                        nodeDy[i] = sr.Displacement[2 * i + 1];
                        switch (field)
                        {
                            case FieldType.Displacement:
                                nodeValues[i] = Math.Sqrt(nodeDx[i] * nodeDx[i] + nodeDy[i] * nodeDy[i]); break;
                            case FieldType.Deformed:
                                nodeValues[i] = Math.Sqrt(nodeDx[i] * nodeDx[i] + nodeDy[i] * nodeDy[i]); break;
                        }
                    }
                    if (field == FieldType.Deformed && deformScale == 0.0)
                    {
                        double maxd = sr.MaxDisp;
                        double L = Math.Max(w, h);
                        deformScale = 0.1 * L / Math.Max(maxd, 1e-12);
                    }
                }
                else if (result is NSResult nr)
                {
                    nodeValues = new double[mesh.NumNodes];
                    for (int i = 0; i < mesh.NumNodes; i++)
                    {
                        switch (field)
                        {
                            case FieldType.Velocity: nodeValues[i] = nr.Speed[i]; break;
                            case FieldType.Pressure: nodeValues[i] = nr.P[i]; break;
                        }
                    }
                }
                // 单元值（应力型）
                double[] cellValues = null;
                if (result is StructuralResult sr2 && (field == FieldType.StressXX || field == FieldType.StressYY || field == FieldType.StressXY || field == FieldType.VonMises))
                {
                    cellValues = new double[mesh.NumElements];
                    for (int e = 0; e < mesh.NumElements; e++)
                    {
                        switch (field)
                        {
                            case FieldType.StressXX: cellValues[e] = sr2.StressXX[e]; break;
                            case FieldType.StressYY: cellValues[e] = sr2.StressYY[e]; break;
                            case FieldType.StressXY: cellValues[e] = sr2.StressXY[e]; break;
                            case FieldType.VonMises: cellValues[e] = sr2.VonMises[e]; break;
                        }
                    }
                    // 将单元值扩散到节点（平均到节点）
                    nodeValues = new double[mesh.NumNodes];
                    int[] cnt = new int[mesh.NumNodes];
                    for (int e = 0; e < mesh.NumElements; e++)
                    {
                        foreach (int n in mesh.Elements[e])
                        {
                            nodeValues[n] += cellValues[e];
                            cnt[n]++;
                        }
                    }
                    for (int i = 0; i < mesh.NumNodes; i++) if (cnt[i] > 0) nodeValues[i] /= cnt[i];
                }

                if (nodeValues != null)
                {
                    vmin = double.PositiveInfinity; vmax = double.NegativeInfinity;
                    for (int i = 0; i < nodeValues.Length; i++)
                    {
                        if (nodeValues[i] < vmin) vmin = nodeValues[i];
                        if (nodeValues[i] > vmax) vmax = nodeValues[i];
                    }
                    if (Math.Abs(vmax - vmin) < 1e-30) { vmax = vmin + 1e-6; }
                }

                // 绘制单元填充
                for (int e = 0; e < mesh.NumElements; e++)
                {
                    int[] c = mesh.Elements[e];
                    int npe = c.Length;
                    // 简化：对于二次单元 T6，仅使用角点绘制三角面（中点用于插值颜色）
                    int nCorners;
                    int[] corners;
                    if (npe == 6) { nCorners = 3; corners = new int[] { c[0], c[1], c[2] }; }
                    else if (npe == 4) { nCorners = 4; corners = new int[] { c[0], c[1], c[2], c[3] }; }
                    else { nCorners = 3; corners = new int[] { c[0], c[1], c[2] }; }

                    PointF[] pts = new PointF[nCorners];
                    Color[] colors = new Color[nCorners];
                    for (int k = 0; k < nCorners; k++)
                    {
                        int nd = corners[k];
                        double x = mesh.X(nd), y = mesh.Y(nd);
                        if (nodeDx != null && deformScale > 0 && (field == FieldType.Deformed || field == FieldType.Displacement || showDeformedOutline))
                        {
                            x += nodeDx[nd] * deformScale;
                            y += nodeDy[nd] * deformScale;
                        }
                        pts[k] = toPx(x, y);
                        if (nodeValues != null)
                            colors[k] = JetColor((nodeValues[nd] - vmin) / (vmax - vmin));
                        else
                            colors[k] = Color.White;
                    }

                    if (field == FieldType.Mesh)
                    {
                        using (Pen p = new Pen(Color.DarkSlateGray, 0.7f)) g.DrawPolygon(p, pts);
                    }
                    else
                    {
                        // 三角形/四边形线性插值填充：分两三角做渐变
                        using (PathGradientBrush pgb = new PathGradientBrush(pts))
                        {
                            pgb.CenterColor = Color.FromArgb(
                                (colors[0].R + colors[1].R + colors[2].R) / 3,
                                (colors[0].G + colors[1].G + colors[2].G) / 3,
                                (colors[0].B + colors[1].B + colors[2].B) / 3);
                            if (nCorners == 3)
                            {
                                pgb.SurroundColors = new Color[] { colors[0], colors[1], colors[2] };
                                g.FillPolygon(pgb, pts);
                            }
                            else
                            {
                                // Q4：分两个三角
                                PointF[] tri1 = new PointF[] { pts[0], pts[1], pts[2] };
                                Color[] c1 = new Color[] { colors[0], colors[1], colors[2] };
                                using (PathGradientBrush b1 = new PathGradientBrush(tri1))
                                {
                                    b1.CenterColor = Color.FromArgb((c1[0].R + c1[1].R + c1[2].R) / 3, (c1[0].G + c1[1].G + c1[2].G) / 3, (c1[0].B + c1[1].B + c1[2].B) / 3);
                                    b1.SurroundColors = c1;
                                    g.FillPolygon(b1, tri1);
                                }
                                PointF[] tri2 = new PointF[] { pts[0], pts[2], pts[3] };
                                Color[] c2 = new Color[] { colors[0], colors[2], colors[3] };
                                using (PathGradientBrush b2 = new PathGradientBrush(tri2))
                                {
                                    b2.CenterColor = Color.FromArgb((c2[0].R + c2[1].R + c2[2].R) / 3, (c2[0].G + c2[1].G + c2[2].G) / 3, (c2[0].B + c2[1].B + c2[2].B) / 3);
                                    b2.SurroundColors = c2;
                                    g.FillPolygon(b2, tri2);
                                }
                            }
                        }
                    }
                }
                // 网格线
                if (field != FieldType.Mesh)
                {
                    using (Pen pen = new Pen(Color.FromArgb(120, Color.Black), 0.4f))
                    {
                        for (int e = 0; e < mesh.NumElements; e++)
                        {
                            int[] c = mesh.Elements[e];
                            int npe = c.Length;
                            int[] corners;
                            if (npe == 6) corners = new int[] { c[0], c[1], c[2] };
                            else if (npe == 4) corners = c;
                            else corners = c;
                            PointF[] pts = new PointF[corners.Length];
                            for (int k = 0; k < corners.Length; k++)
                            {
                                int nd = corners[k];
                                double x = mesh.X(nd), y = mesh.Y(nd);
                                if (nodeDx != null && deformScale > 0 && (field == FieldType.Deformed || showDeformedOutline))
                                { x += nodeDx[nd] * deformScale; y += nodeDy[nd] * deformScale; }
                                pts[k] = toPx(x, y);
                            }
                            g.DrawPolygon(pen, pts);
                        }
                    }
                }

                // 色条图例
                if (field != FieldType.Mesh && nodeValues != null)
                {
                    int barW = 20, barH = height - 100;
                    int bx0 = width - 50, by0 = 50;
                    for (int i = 0; i < barH; i++)
                    {
                        double t = 1.0 - (double)i / barH;
                        Color cl = JetColor(t);
                        using (Pen p = new Pen(cl)) g.DrawLine(p, bx0, by0 + i, bx0 + barW, by0 + i);
                    }
                    g.DrawRectangle(Pens.Black, bx0, by0, barW, barH);
                    using (Font f = new Font("Consolas", 9))
                    {
                        StringFormat sf = new StringFormat();
                        g.DrawString(vmax.ToString("G3"), f, Brushes.Black, bx0 + barW + 4, by0 - 6);
                        g.DrawString(vmin.ToString("G3"), f, Brushes.Black, bx0 + barW + 4, by0 + barH - 10);
                        g.DrawString(((vmax + vmin) * 0.5).ToString("G3"), f, Brushes.Black, bx0 + barW + 4, by0 + barH / 2 - 6);
                        string name = field.ToString();
                        g.DrawString(name, f, Brushes.Black, bx0 - 10, by0 - 25);
                    }
                }
                // 标题
                using (Font titleF = new Font("Arial", 12, FontStyle.Bold))
                {
                    string title = "FEM2D - " + field.ToString();
                    g.DrawString(title, titleF, Brushes.DarkBlue, 10, 10);
                }
            }
            return bmp;
        }

        private static Color JetColor(double t)
        {
            if (t < 0) t = 0; if (t > 1) t = 1;
            double r, g, b;
            if (t < 0.25) { r = 0; g = 4 * t; b = 1; }
            else if (t < 0.5) { r = 0; g = 1; b = 1 - 4 * (t - 0.25); }
            else if (t < 0.75) { r = 4 * (t - 0.5); g = 1; b = 0; }
            else { r = 1; g = 1 - 4 * (t - 0.75); b = 0; }
            return Color.FromArgb((int)(255 * Clamp(r)), (int)(255 * Clamp(g)), (int)(255 * Clamp(b)));
        }

        private static double Clamp(double v) { return v < 0 ? 0 : v > 1 ? 1 : v; }
    }
}
