using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Drawing;
using System.Drawing.Imaging;
using Fem2DFluid.Models;

namespace Fem2DFluid.Services
{
    /// <summary>
    /// 结果输出：数据文件(CSV/Tecplot ASCII DAT)+彩色云图PNG
    /// </summary>
    public class ResultExporter
    {
        /// <summary>
        /// 导出节点解 CSV
        /// </summary>
        public static void ExportNodeCsv(FemModel model, string path)
        {
            using (StreamWriter sw = new StreamWriter(path, false, System.Text.Encoding.UTF8))
            {
                sw.WriteLine("Id,X,Y,Phi,BCType,BCValue");
                for (int i = 0; i < model.Nodes.Count; i++)
                {
                    Node n = model.Nodes[i];
                    sw.WriteLine(string.Format(CultureInfo.InvariantCulture,
                        "{0},{1:F6},{2:F6},{3:F6},{4},{5:F6}",
                        n.Id, n.X, n.Y, n.Phi, n.BCType, n.BCValue));
                }
            }
        }

        /// <summary>
        /// 导出单元解 CSV
        /// </summary>
        public static void ExportElementCsv(FemModel model, string path)
        {
            using (StreamWriter sw = new StreamWriter(path, false, System.Text.Encoding.UTF8))
            {
                sw.WriteLine("Id,N1,N2,N3,MatId,Cx,Cy,Area,Vx,Vy,Vmag");
                for (int i = 0; i < model.Elements.Count; i++)
                {
                    Element e = model.Elements[i];
                    sw.WriteLine(string.Format(CultureInfo.InvariantCulture,
                        "{0},{1},{2},{3},{4},{5:F6},{6:F6},{7:F6},{8:F6},{9:F6},{10:F6}",
                        e.Id, e.N1, e.N2, e.N3, e.MatId, e.Cx, e.Cy, e.Area, e.Vx, e.Vy, e.Vmag));
                }
            }
        }

        /// <summary>
        /// 导出Tecplot格式DAT（可被Tecplot/Paraview读取）
        /// </summary>
        public static void ExportTecplot(FemModel model, string path)
        {
            using (StreamWriter sw = new StreamWriter(path, false, System.Text.Encoding.ASCII))
            {
                sw.WriteLine("TITLE = \"" + model.Title + "\"");
                sw.WriteLine("VARIABLES = \"X\", \"Y\", \"Phi\"");
                sw.WriteLine(string.Format("ZONE T=\"Result\", N={0}, E={1}, DATAPACKING=POINT, ZONETYPE=FETRIANGLE",
                    model.Nodes.Count, model.Elements.Count));
                for (int i = 0; i < model.Nodes.Count; i++)
                {
                    Node n = model.Nodes[i];
                    sw.WriteLine(string.Format(CultureInfo.InvariantCulture, "{0:F6} {1:F6} {2:F6}", n.X, n.Y, n.Phi));
                }
                for (int i = 0; i < model.Elements.Count; i++)
                {
                    Element e = model.Elements[i];
                    sw.WriteLine(string.Format("{0} {1} {2}", e.N1, e.N2, e.N3));
                }
            }
        }

        /// <summary>
        /// 绘制彩色云图（势函数填充等值色+网格线+速度矢量+图例色条）
        /// scalarField: "Phi" 或 "Vmag"
        /// </summary>
        public static Bitmap RenderContour(FemModel model, string scalarField, int width, int height,
            bool showMesh, bool showVectors, out double vmin, out double vmax, out List<LegendItem> legend)
        {
            legend = new List<LegendItem>();
            vmin = double.MaxValue; vmax = double.MinValue;
            double[] values = new double[model.Elements.Count];
            for (int i = 0; i < model.Elements.Count; i++)
            {
                Element e = model.Elements[i];
                double v = 0;
                if (scalarField == "Vmag") v = e.Vmag;
                else v = (model.GetNode(e.N1).Phi + model.GetNode(e.N2).Phi + model.GetNode(e.N3).Phi) / 3.0;
                values[i] = v;
                if (v < vmin) vmin = v;
                if (v > vmax) vmax = v;
            }
            if (Math.Abs(vmax - vmin) < 1e-12) vmax = vmin + 1e-9;

            Bitmap bmp = new Bitmap(width, height);
            using (Graphics g = Graphics.FromImage(bmp))
            {
                g.SmoothingMode = System.Drawing.Drawing2D.SmoothingMode.AntiAlias;
                g.Clear(Color.White);
                double xmin, ymin, xmax, ymax;
                model.ComputeBounds(out xmin, out ymin, out xmax, out ymax);
                double margin = 40;
                double drawW = width - margin * 2 - 100;
                double drawH = height - margin * 2;
                double sx = drawW / Math.Max(1e-9, (xmax - xmin));
                double sy = drawH / Math.Max(1e-9, (ymax - ymin));
                double scale = Math.Min(sx, sy);
                double ox = margin + (drawW - (xmax - xmin) * scale) / 2.0;
                double oy = margin + (drawH - (ymax - ymin) * scale) / 2.0;
                // Y轴翻转（屏幕坐标Y向下）
                System.Func<double, double, PointF> map = (x, y) => new PointF(
                    (float)(ox + (x - xmin) * scale),
                    (float)(oy + (ymax - y) * scale)
                );

                // 1. 填充三角形
                int levels = 24;
                Color[] cmap = BuildJetColormap(levels);
                for (int i = 0; i < model.Elements.Count; i++)
                {
                    Element e = model.Elements[i];
                    Node n1 = model.GetNode(e.N1);
                    Node n2 = model.GetNode(e.N2);
                    Node n3 = model.GetNode(e.N3);
                    PointF p1 = map(n1.X, n1.Y);
                    PointF p2 = map(n2.X, n2.Y);
                    PointF p3 = map(n3.X, n3.Y);
                    double t = (values[i] - vmin) / (vmax - vmin);
                    int li = (int)Math.Min(levels - 1, Math.Max(0, t * levels));
                    using (Brush br = new SolidBrush(cmap[li]))
                    {
                        g.FillPolygon(br, new PointF[] { p1, p2, p3 });
                    }
                }

                // 2. 网格线
                if (showMesh)
                {
                    using (Pen p = new Pen(Color.FromArgb(80, Color.Black), 0.5f))
                    {
                        for (int i = 0; i < model.Elements.Count; i++)
                        {
                            Element e = model.Elements[i];
                            Node n1 = model.GetNode(e.N1);
                            Node n2 = model.GetNode(e.N2);
                            Node n3 = model.GetNode(e.N3);
                            PointF p1 = map(n1.X, n1.Y);
                            PointF p2 = map(n2.X, n2.Y);
                            PointF p3 = map(n3.X, n3.Y);
                            g.DrawLine(p, p1, p2);
                            g.DrawLine(p, p2, p3);
                            g.DrawLine(p, p3, p1);
                        }
                    }
                }

                // 3. 边界条件节点标记
                for (int i = 0; i < model.Nodes.Count; i++)
                {
                    Node n = model.Nodes[i];
                    if (n.BCType == 1)
                    {
                        PointF pt = map(n.X, n.Y);
                        g.FillEllipse(Brushes.Red, pt.X - 2, pt.Y - 2, 4, 4);
                    }
                    else if (n.BCType == 2)
                    {
                        PointF pt = map(n.X, n.Y);
                        g.FillEllipse(Brushes.Blue, pt.X - 2, pt.Y - 2, 4, 4);
                    }
                }

                // 4. 速度矢量
                if (showVectors)
                {
                    using (Pen vp = new Pen(Color.Black, 1f))
                    {
                        vp.EndCap = System.Drawing.Drawing2D.LineCap.ArrowAnchor;
                        double refLen = 0.04 * Math.Max(xmax - xmin, ymax - ymin);
                        double vmaxGlobal = vmax;
                        if (vmaxGlobal < 1e-9) vmaxGlobal = 1.0;
                        for (int i = 0; i < model.Elements.Count; i++)
                        {
                            Element e = model.Elements[i];
                            PointF c = map(e.Cx, e.Cy);
                            double vmag = e.Vmag;
                            double vx = e.Vx, vy = e.Vy;
                            if (vmag < 1e-9) continue;
                            double len = refLen * vmag / vmaxGlobal;
                            double dx = vx / vmag * len;
                            double dy = -vy / vmag * len; // Y翻转
                            PointF end = new PointF(c.X + (float)(dx * scale), c.Y + (float)(dy * scale));
                            g.DrawLine(vp, c, end);
                        }
                    }
                }

                // 5. 坐标边框
                using (Pen frame = new Pen(Color.Black, 1f))
                {
                    g.DrawRectangle(frame, (float)ox, (float)oy, (float)((xmax - xmin) * scale), (float)((ymax - ymin) * scale));
                    Font sf = new Font("Arial", 7f);
                    g.DrawString(string.Format(CultureInfo.InvariantCulture, "Xmin={0:F2}", xmin), sf, Brushes.Black, (float)ox, height - (float)margin + 5);
                    g.DrawString(string.Format(CultureInfo.InvariantCulture, "Xmax={0:F2}", xmax), sf, Brushes.Black, (float)(ox + (xmax - xmin) * scale - 50), height - (float)margin + 5);
                    g.DrawString(string.Format(CultureInfo.InvariantCulture, "Ymin={0:F2}", ymin), sf, Brushes.Black, 2, (float)(oy + (ymax - ymin) * scale - 8));
                    g.DrawString(string.Format(CultureInfo.InvariantCulture, "Ymax={0:F2}", ymax), sf, Brushes.Black, 2, (float)oy);
                }

                // 6. 图例色条（右侧）
                int barX = (int)(ox + (xmax - xmin) * scale + 20);
                int barY = (int)oy;
                int barW = 20;
                int barH = (int)((ymax - ymin) * scale);
                for (int k = 0; k < barH; k++)
                {
                    double t = 1.0 - (double)k / barH;
                    int li = (int)Math.Min(levels - 1, Math.Max(0, t * levels));
                    using (Pen pc = new Pen(cmap[li]))
                    {
                        g.DrawLine(pc, barX, barY + k, barX + barW, barY + k);
                    }
                }
                using (Pen bp = new Pen(Color.Black)) { g.DrawRectangle(bp, barX, barY, barW, barH); }
                Font lf = new Font("Arial", 8f);
                StringFormat rfmt = new StringFormat();
                rfmt.Alignment = StringAlignment.Near;
                string fieldName = (scalarField == "Vmag") ? "|V|" : "φ";
                g.DrawString(fieldName + " max", lf, Brushes.Black, barX + barW + 4, barY - 10);
                g.DrawString(string.Format(CultureInfo.InvariantCulture, "{0:F4}", vmax), lf, Brushes.Black, barX + barW + 4, barY);
                g.DrawString(string.Format(CultureInfo.InvariantCulture, "{0:F4}", (vmax + vmin) / 2), lf, Brushes.Black, barX + barW + 4, barY + barH / 2 - 6);
                g.DrawString(string.Format(CultureInfo.InvariantCulture, "{0:F4}", vmin), lf, Brushes.Black, barX + barW + 4, barY + barH - 12);
                g.DrawString(fieldName + " min", lf, Brushes.Black, barX + barW + 4, barY + barH);

                // 7. 标题
                Font tf = new Font("Arial", 12f, FontStyle.Bold);
                string title = model.Title + "  [" + fieldName + "]  节点=" + model.Nodes.Count + " 单元=" + model.Elements.Count;
                SizeF ts = g.MeasureString(title, tf);
                g.DrawString(title, tf, Brushes.Black, (width - ts.Width) / 2, 8);

                // 图例数据返回
                for (int lv = 0; lv < 6; lv++)
                {
                    double t = 1.0 - (double)lv / 5.0;
                    double val = vmin + t * (vmax - vmin);
                    int li = (int)Math.Min(levels - 1, Math.Max(0, t * levels));
                    LegendItem it = new LegendItem();
                    it.Value = val;
                    it.Color = cmap[li];
                    legend.Add(it);
                }
            }
            return bmp;
        }

        private static Color[] BuildJetColormap(int levels)
        {
            Color[] cmap = new Color[levels];
            for (int i = 0; i < levels; i++)
            {
                double v = (double)i / (levels - 1); // 0~1
                double r, g, b;
                // Jet colormap
                if (v < 0.25) { r = 0; g = 4 * v; b = 1; }
                else if (v < 0.5) { r = 0; g = 1; b = 1 - 4 * (v - 0.25); }
                else if (v < 0.75) { r = 4 * (v - 0.5); g = 1; b = 0; }
                else { r = 1; g = 1 - 4 * (v - 0.75); b = 0; }
                int ri = (int)Math.Max(0, Math.Min(255, r * 255));
                int gi = (int)Math.Max(0, Math.Min(255, g * 255));
                int bi = (int)Math.Max(0, Math.Min(255, b * 255));
                cmap[i] = Color.FromArgb(ri, gi, bi);
            }
            return cmap;
        }
    }

    public class LegendItem
    {
        public double Value;
        public Color Color;
    }
}
