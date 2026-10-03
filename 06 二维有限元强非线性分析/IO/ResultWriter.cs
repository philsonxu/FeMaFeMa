using System;
using System.Globalization;
using System.IO;
using System.Drawing;
using System.Drawing.Imaging;
using NonlinearFEM2D.Mesh;
using NonlinearFEM2D.Elements;
using NonlinearFEM2D.Rendering;

namespace NonlinearFEM2D.IO
{
    public static class ResultWriter
    {
        public static void WriteCSV(FEMesh mesh, string path)
        {
            using (StreamWriter sw = new StreamWriter(Path.Combine(path, "nodes.csv")))
            {
                sw.WriteLine("id,x,y,u,v,ux,uy,p,sxx,syy,sxy,vm,eqp");
                for (int i = 0; i < mesh.NumNodes; i++)
                {
                    sw.WriteLine(string.Format(CultureInfo.InvariantCulture,
                        "{0},{1:F6},{2:F6},{3:F8},{4:F8},{5:F8},{6:F8},{7:F8},{8:F6},{9:F6},{10:F6},{11:F6},{12:E4}",
                        i, mesh.X[i], mesh.Y[i],
                        mesh.DisplacementU[i], mesh.DisplacementV[i],
                        mesh.VelocityU[i], mesh.VelocityV[i],
                        mesh.Pressure[i],
                        mesh.StressXX[i], mesh.StressYY[i], mesh.StressXY[i],
                        mesh.VonMises[i], mesh.EqPlasticStrain[i]));
                }
            }
            using (StreamWriter sw = new StreamWriter(Path.Combine(path, "elements.csv")))
            {
                sw.WriteLine("id,type,nodes");
                foreach (FiniteElement e in mesh.Elements)
                {
                    sw.Write("{0},{1},", e.Id, e.Type);
                    for (int i = 0; i < e.NodeIds.Length; i++)
                    {
                        sw.Write(e.NodeIds[i]);
                        if (i + 1 < e.NodeIds.Length) sw.Write(" ");
                    }
                    sw.WriteLine();
                }
            }
        }

        public static void WriteVTK(FEMesh mesh, string path)
        {
            using (StreamWriter sw = new StreamWriter(Path.Combine(path, "result.vtk")))
            {
                sw.WriteLine("# vtk DataFile Version 3.0");
                sw.WriteLine("NonlinearFEM2D result");
                sw.WriteLine("ASCII");
                sw.WriteLine("DATASET UNSTRUCTURED_GRID");
                sw.WriteLine("POINTS {0} double", mesh.NumNodes);
                for (int i = 0; i < mesh.NumNodes; i++)
                    sw.WriteLine(string.Format(CultureInfo.InvariantCulture, "{0:F6} {1:F6} 0.0", mesh.X[i], mesh.Y[i]));
                int totConn = 0;
                foreach (FiniteElement e in mesh.Elements) totConn += e.NodeIds.Length + 1;
                sw.WriteLine("CELLS {0} {1}", mesh.Elements.Count, totConn);
                foreach (FiniteElement e in mesh.Elements)
                {
                    sw.Write(e.NodeIds.Length);
                    for (int i = 0; i < e.NodeIds.Length; i++) sw.Write(" " + e.NodeIds[i]);
                    sw.WriteLine();
                }
                sw.WriteLine("CELL_TYPES {0}", mesh.Elements.Count);
                foreach (FiniteElement e in mesh.Elements)
                {
                    if (e.Type == ElementType.CST3) sw.WriteLine("5");
                    else if (e.Type == ElementType.LT6) sw.WriteLine("22");
                    else sw.WriteLine("9");
                }
                sw.WriteLine("POINT_DATA {0}", mesh.NumNodes);
                sw.WriteLine("VECTORS displacement double");
                for (int i = 0; i < mesh.NumNodes; i++)
                    sw.WriteLine(string.Format(CultureInfo.InvariantCulture, "{0:F8} {1:F8} 0.0", mesh.DisplacementU[i], mesh.DisplacementV[i]));
                sw.WriteLine("VECTORS velocity double");
                for (int i = 0; i < mesh.NumNodes; i++)
                    sw.WriteLine(string.Format(CultureInfo.InvariantCulture, "{0:F8} {1:F8} 0.0", mesh.VelocityU[i], mesh.VelocityV[i]));
                sw.WriteLine("SCALARS pressure double 1"); sw.WriteLine("LOOKUP_TABLE default");
                for (int i = 0; i < mesh.NumNodes; i++) sw.WriteLine(string.Format(CultureInfo.InvariantCulture, "{0:F8}", mesh.Pressure[i]));
                sw.WriteLine("SCALARS vonMises double 1"); sw.WriteLine("LOOKUP_TABLE default");
                for (int i = 0; i < mesh.NumNodes; i++) sw.WriteLine(string.Format(CultureInfo.InvariantCulture, "{0:F6}", mesh.VonMises[i]));
                sw.WriteLine("SCALARS eqPlasticStrain double 1"); sw.WriteLine("LOOKUP_TABLE default");
                for (int i = 0; i < mesh.NumNodes; i++) sw.WriteLine(string.Format(CultureInfo.InvariantCulture, "{0:E4}", mesh.EqPlasticStrain[i]));
                sw.WriteLine("SCALARS sxx double 1"); sw.WriteLine("LOOKUP_TABLE default");
                for (int i = 0; i < mesh.NumNodes; i++) sw.WriteLine(string.Format(CultureInfo.InvariantCulture, "{0:F6}", mesh.StressXX[i]));
            }
        }

        public static void WritePNG(FEMesh mesh, string field, string path, int modeIdx)
        {
            ContourRenderer cr = new ContourRenderer();
            Bitmap bmp = cr.Render(mesh, field, 1600, 1200, modeIdx);
            bmp.Save(Path.Combine(path, "contour_" + field + ".png"), ImageFormat.Png);
            bmp.Dispose();
        }
    }
}
