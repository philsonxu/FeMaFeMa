namespace ModalFEM2D.Rendering
{
    using System;
    using System.Collections.Generic;
    using System.Drawing;
    using System.Drawing.Imaging;
    using System.Globalization;
    using ModalFEM2D.Elements;
    using ModalFEM2D.Mesh;

    public enum DisplayField
    {
        Mesh,
        Displacement,
        Deformed,
        StressXX,
        StressYY,
        StressXY,
        VonMises,
        Velocity,
        Pressure,
        ModeShape
    }

    public static class ContourRenderer
    {
        public static Bitmap Render(FEMesh mesh, DisplayField field, int modeIndex, double deformScale,
                                     int width, int height)
        {
            Bitmap bmp = new Bitmap(width, height);
            using (Graphics g = Graphics.FromImage(bmp))
            {
                g.Clear(Color.Black);
                if (mesh == null || mesh.NumNodes == 0) { return bmp; }
                // bounds
                double xmin = mesh.X[0], xmax = mesh.X[0], ymin = mesh.Y[0], ymax = mesh.Y[0];
                for (int i = 0; i < mesh.NumNodes; i++)
                {
                    if (mesh.X[i] < xmin) xmin = mesh.X[i]; if (mesh.X[i] > xmax) xmax = mesh.X[i];
                    if (mesh.Y[i] < ymin) ymin = mesh.Y[i]; if (mesh.Y[i] > ymax) ymax = mesh.Y[i];
                }
                double dx = xmax - xmin; double dy = ymax - ymin;
                if (dx < 1e-12) dx = 1.0; if (dy < 1e-12) dy = 1.0;
                int margin = 60;
                int plotW = width - margin - 10; int plotH = height - margin - 60;
                if (plotW < 20) plotW = 20; if (plotH < 20) plotH = 20;
                double sx = plotW / dx; double sy = plotH / dy;
                double s = Math.Min(sx, sy);
                int offX = margin + (int)((plotW - dx * s) * 0.5);
                int offY = (int)((plotH - dy * s) * 0.5) + 5;
                Func<double, double, PointF> map = (x, y) => new PointF((float)(offX + (x - xmin) * s), (float)(offY + (ymax - y) * s));

                // gather field values
                double[] fld = new double[mesh.NumNodes];
                if (field == DisplayField.Mesh) { for (int i = 0; i < mesh.NumNodes; i++) fld[i] = 0.0; }
                else if (field == DisplayField.Deformed || field == DisplayField.Displacement)
                { for (int i = 0; i < mesh.NumNodes; i++) fld[i] = Math.Sqrt(mesh.DisplacementU[i] * mesh.DisplacementU[i] + mesh.DisplacementV[i] * mesh.DisplacementV[i]); }
                else if (field == DisplayField.StressXX) fld = (double[])mesh.StressXX.Clone();
                else if (field == DisplayField.StressYY) fld = (double[])mesh.StressYY.Clone();
                else if (field == DisplayField.StressXY) fld = (double[])mesh.StressXY.Clone();
                else if (field == DisplayField.VonMises) fld = (double[])mesh.VonMises.Clone();
                else if (field == DisplayField.Velocity)
                { for (int i = 0; i < mesh.NumNodes; i++) fld[i] = Math.Sqrt(mesh.VelocityU[i] * mesh.VelocityU[i] + mesh.VelocityV[i] * mesh.VelocityV[i]); }
                else if (field == DisplayField.Pressure) fld = (double[])mesh.Pressure.Clone();
                else if (field == DisplayField.ModeShape)
                {
                    if (mesh.ModalShapes != null && modeIndex >= 0 && modeIndex < mesh.ModalFreqHz.Length)
                        for (int i = 0; i < mesh.NumNodes; i++)
                            fld[i] = Math.Sqrt(mesh.ModalShapes[modeIndex, 2 * i] * mesh.ModalShapes[modeIndex, 2 * i] +
                                               mesh.ModalShapes[modeIndex, 2 * i + 1] * mesh.ModalShapes[modeIndex, 2 * i + 1]);
                }
                double vmin = fld[0], vmax = fld[0];
                for (int i = 1; i < mesh.NumNodes; i++) { if (fld[i] < vmin) vmin = fld[i]; if (fld[i] > vmax) vmax = fld[i]; }
                if (vmax - vmin < 1e-14) vmax = vmin + 1.0;

                // get disp/shape for deformed overlay
                double[] ux = new double[mesh.NumNodes]; double[] uy = new double[mesh.NumNodes];
                bool deform = (field == DisplayField.Deformed || field == DisplayField.ModeShape);
                if (deform)
                {
                    double maxDisp = 0.0;
                    double[] srcU = mesh.DisplacementU, srcV = mesh.DisplacementV;
                    if (field == DisplayField.ModeShape && mesh.ModalShapes != null && modeIndex >= 0)
                    {
                        for (int i = 0; i < mesh.NumNodes; i++)
                        {
                            ux[i] = mesh.ModalShapes[modeIndex, 2 * i];
                            uy[i] = mesh.ModalShapes[modeIndex, 2 * i + 1];
                            double a = Math.Max(Math.Abs(ux[i]), Math.Abs(uy[i]));
                            if (a > maxDisp) maxDisp = a;
                        }
                    }
                    else
                    {
                        for (int i = 0; i < mesh.NumNodes; i++)
                        { ux[i] = mesh.DisplacementU[i]; uy[i] = mesh.DisplacementV[i]; double a = Math.Max(Math.Abs(ux[i]), Math.Abs(uy[i])); if (a > maxDisp) maxDisp = a; }
                    }
                    // scale so max displacement ~ deformScale * smaller of dx/dy
                    double target = deformScale * 0.15 * Math.Min(dx, dy);
                    double fact = maxDisp > 1e-14 ? target / maxDisp : 0.0;
                    for (int i = 0; i < mesh.NumNodes; i++) { ux[i] *= fact; uy[i] *= fact; }
                }

                // Draw filled cells with per-nodal color via Gouraud-like approach using GDI Polygon + average cell color (simple and fast)
                for (int e = 0; e < mesh.NumElements; e++)
                {
                    int npe = mesh.NodesPerElement;
                    PointF[] pts = new PointF[npe];
                    double avgF = 0.0;
                    for (int k = 0; k < npe; k++)
                    {
                        int gid = mesh.Connectivity[e, k];
                        double x = mesh.X[gid], y = mesh.Y[gid];
                        if (deform) { x += ux[gid]; y += uy[gid]; }
                        pts[k] = map(x, y);
                        avgF += fld[gid];
                    }
                    avgF /= npe;
                    double t = (avgF - vmin) / (vmax - vmin); if (t < 0) t = 0; if (t > 1) t = 1;
                    Color c = JetColor(t);
                    using (Brush br = new SolidBrush(c))
                    {
                        if (npe == 3) g.FillPolygon(br, new PointF[] { pts[0], pts[1], pts[2] });
                        else if (npe == 4) g.FillPolygon(br, new PointF[] { pts[0], pts[1], pts[2], pts[3] });
                        else g.FillPolygon(br, pts); // 6-node LT6: 6-point polygon
                    }
                }
                // wireframe
                using (Pen pen = new Pen(Color.FromArgb(70, 255, 255, 255), 0.6f))
                {
                    for (int e = 0; e < mesh.NumElements; e++)
                    {
                        int npe = mesh.NodesPerElement;
                        PointF[] pts = new PointF[npe];
                        for (int k = 0; k < npe; k++)
                        {
                            int gid = mesh.Connectivity[e, k];
                            double x = mesh.X[gid], y = mesh.Y[gid];
                            if (deform) { x += ux[gid]; y += uy[gid]; }
                            pts[k] = map(x, y);
                        }
                        if (npe == 3) { g.DrawLine(pen, pts[0], pts[1]); g.DrawLine(pen, pts[1], pts[2]); g.DrawLine(pen, pts[2], pts[0]); }
                        else if (npe == 4) { g.DrawLine(pen, pts[0], pts[1]); g.DrawLine(pen, pts[1], pts[2]); g.DrawLine(pen, pts[2], pts[3]); g.DrawLine(pen, pts[3], pts[0]); }
                        else { for (int k = 0; k < npe; k++) g.DrawLine(pen, pts[k], pts[(k + 1) % npe]); }
                    }
                }
                // undeformed wireframe (dashed gray) when deformed
                if (deform)
                {
                    using (Pen pd = new Pen(Color.FromArgb(120, 200, 200, 200)) { DashStyle = System.Drawing.Drawing2D.DashStyle.Dash })
                    {
                        for (int e = 0; e < mesh.NumElements; e++)
                        {
                            int npe = mesh.NodesPerElement;
                            PointF[] pts = new PointF[npe];
                            for (int k = 0; k < npe; k++)
                            {
                                int gid = mesh.Connectivity[e, k];
                                pts[k] = map(mesh.X[gid], mesh.Y[gid]);
                            }
                            if (npe == 3) { g.DrawLine(pd, pts[0], pts[1]); g.DrawLine(pd, pts[1], pts[2]); g.DrawLine(pd, pts[2], pts[0]); }
                            else if (npe == 4) { g.DrawLine(pd, pts[0], pts[1]); g.DrawLine(pd, pts[1], pts[2]); g.DrawLine(pd, pts[2], pts[3]); g.DrawLine(pd, pts[3], pts[0]); }
                        }
                    }
                }

                // color bar
                int barX = width - 24; int barY = margin; int barW = 14; int barH = plotH;
                if (barH < 20) barH = 20;
                using (Pen bp = new Pen(Color.White))
                {
                    for (int y = 0; y < barH; y++)
                    {
                        double t = 1.0 - (double)y / (barH - 1);
                        using (Brush br = new SolidBrush(JetColor(t))) g.FillRectangle(br, barX, barY + y, barW, 1);
                    }
                    g.DrawRectangle(bp, barX, barY, barW, barH);
                    using (Font ft = new Font("Consolas", 8))
                    {
                        StringFormat fmt = new StringFormat(); fmt.Alignment = StringAlignment.Far;
                        for (int k = 0; k <= 5; k++)
                        {
                            double t = 1.0 - k / 5.0; double v = vmin + t * (vmax - vmin);
                            int y = barY + (int)(k / 5.0 * (barH - 1));
                            g.DrawLine(bp, barX - 3, y, barX, y);
                            g.DrawString(v.ToString("G3", CultureInfo.InvariantCulture), ft, Brushes.White, barX - 5, y - 6, fmt);
                        }
                    }
                }
                // title
                using (Font ft = new Font("Arial", 11, FontStyle.Bold))
                {
                    string title = field.ToString();
                    if (field == DisplayField.ModeShape && mesh.ModalFreqHz != null)
                        title = string.Format("Mode {0}  f = {1:F2} Hz", modeIndex + 1, mesh.ModalFreqHz[modeIndex]);
                    g.DrawString(title, ft, Brushes.White, margin, 8);
                    g.DrawString(string.Format("Nodes={0}  Elements={1}  Type={2}", mesh.NumNodes, mesh.NumElements, mesh.ElementType),
                                 new Font("Arial", 8), Brushes.LightGray, margin, 28);
                }
            }
            return bmp;
        }

        private static Color JetColor(double t)
        {
            if (t < 0) t = 0; if (t > 1) t = 1;
            double r, g, b;
            if (t < 0.25) { r = 0.0; g = t / 0.25 * 0.5; b = 1.0; }
            else if (t < 0.5) { r = 0.0; g = 0.5 + (t - 0.25) / 0.25 * 0.5; b = 1.0 - (t - 0.25) / 0.25; }
            else if (t < 0.75) { r = (t - 0.5) / 0.25; g = 1.0; b = 0.0; }
            else { r = 1.0; g = 1.0 - (t - 0.75) / 0.25; b = 0.0; }
            return Color.FromArgb(255, (int)(255 * r), (int)(255 * g), (int)(255 * b));
        }
    }
}
