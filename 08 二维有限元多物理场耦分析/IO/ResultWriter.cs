using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using MultiPhysicsFEM2D.Elements;
using MultiPhysicsFEM2D.Mesh;
using MultiPhysicsFEM2D.Rendering;

namespace MultiPhysicsFEM2D.IO
{
    public static class ResultWriter
    {
        public static void WriteAll(FEMesh mesh, string directory, string prefix, string field, double deformScale)
        {
            Directory.CreateDirectory(directory);
            WriteCSV(mesh, Path.Combine(directory, prefix + "_nodes.csv"), Path.Combine(directory, prefix + "_elements.csv"));
            WriteVTK(mesh, Path.Combine(directory, prefix + ".vtk"));
            string pngPath = Path.Combine(directory, prefix + "_" + SafeField(field) + ".png");
            ContourRenderer.RenderToPng(mesh, pngPath, 1600, 1200, field, deformScale);
        }

        private static string SafeField(string f)
        {
            string s = f.Replace("σ", "s").Replace("τ", "t").Replace(" ", "").Replace("(", "").Replace(")", "");
            return s;
        }

        public static void WriteCSV(FEMesh mesh, string nodePath, string elemPath)
        {
            using (StreamWriter w = new StreamWriter(nodePath))
            {
                w.WriteLine("id,x,y,ux,uy,vx,vy,p,T,sxx,syy,sxy,vm,ep");
                for (int i = 0; i < mesh.NumNodes; i++)
                {
                    w.WriteLine(string.Format(CultureInfo.InvariantCulture,
                        "{0},{1:F6},{2:F6},{3:E6},{4:E6},{5:E6},{6:E6},{7:E6},{8:F4},{9:E4},{10:E4},{11:E4},{12:E4},{13:E4}",
                        i, mesh.X[i], mesh.Y[i],
                        mesh.DisplacementU[i], mesh.DisplacementV[i],
                        mesh.VelocityU[i], mesh.VelocityV[i],
                        mesh.Pressure[i], mesh.Temperature[i],
                        mesh.StressXX[i], mesh.StressYY[i], mesh.StressXY[i],
                        mesh.VonMises[i], mesh.PlasticStrain[i]));
                }
            }
            using (StreamWriter w = new StreamWriter(elemPath))
            {
                w.WriteLine("id,type,nodes...,sxx,syy,sxy,vm");
                foreach (FiniteElement e in mesh.Elements)
                {
                    string nodes = "";
                    for (int k = 0; k < e.NodeIds.Length; k++) nodes += "," + e.NodeIds[k];
                    w.WriteLine(string.Format(CultureInfo.InvariantCulture,
                        "{0},{1}{2},{3:E4},{4:E4},{5:E4},{6:E4}",
                        e.Id, e.Type, nodes, e.StressXX, e.StressYY, e.StressXY, e.VonMises));
                }
            }
        }

        public static void WriteVTK(FEMesh mesh, string path)
        {
            using (StreamWriter w = new StreamWriter(path))
            {
                w.WriteLine("# vtk DataFile Version 3.0");
                w.WriteLine("FEM2D output");
                w.WriteLine("ASCII");
                w.WriteLine("DATASET UNSTRUCTURED_GRID");
                w.WriteLine("POINTS " + mesh.NumNodes + " double");
                for (int i = 0; i < mesh.NumNodes; i++)
                    w.WriteLine(string.Format(CultureInfo.InvariantCulture, "{0:F6} {1:F6} 0.0", mesh.X[i], mesh.Y[i]));
                int totalConn = 0;
                foreach (FiniteElement e in mesh.Elements) totalConn += e.NodeIds.Length + 1;
                w.WriteLine("CELLS " + mesh.Elements.Count + " " + totalConn);
                foreach (FiniteElement e in mesh.Elements)
                {
                    w.Write(e.NodeIds.Length);
                    for (int k = 0; k < e.NodeIds.Length; k++) w.Write(" " + e.NodeIds[k]);
                    w.WriteLine();
                }
                w.WriteLine("CELL_TYPES " + mesh.Elements.Count);
                foreach (FiniteElement e in mesh.Elements)
                {
                    int t = 5;
                    if (e.Type == ElementType.CST3) t = 5;        // VTK_TRIANGLE
                    else if (e.Type == ElementType.Q4) t = 9;   // VTK_QUAD
                    else if (e.Type == ElementType.LT6) t = 22; // VTK_QUADRATIC_TRIANGLE
                    w.WriteLine(t);
                }
                w.WriteLine("POINT_DATA " + mesh.NumNodes);
                w.WriteLine("VECTORS displacement double");
                for (int i = 0; i < mesh.NumNodes; i++)
                    w.WriteLine(string.Format(CultureInfo.InvariantCulture, "{0:E6} {1:E6} 0", mesh.DisplacementU[i], mesh.DisplacementV[i]));
                w.WriteLine("VECTORS velocity double");
                for (int i = 0; i < mesh.NumNodes; i++)
                    w.WriteLine(string.Format(CultureInfo.InvariantCulture, "{0:E6} {1:E6} 0", mesh.VelocityU[i], mesh.VelocityV[i]));
                WriteScalar(w, "pressure", mesh.Pressure);
                WriteScalar(w, "temperature", mesh.Temperature);
                WriteScalar(w, "stress_xx", mesh.StressXX);
                WriteScalar(w, "stress_yy", mesh.StressYY);
                WriteScalar(w, "stress_xy", mesh.StressXY);
                WriteScalar(w, "von_mises", mesh.VonMises);
                WriteScalar(w, "plastic_strain", mesh.PlasticStrain);
            }
        }

        private static void WriteScalar(StreamWriter w, string name, double[] a)
        {
            w.WriteLine("SCALARS " + name + " double 1");
            w.WriteLine("LOOKUP_TABLE default");
            for (int i = 0; i < a.Length; i++)
                w.WriteLine(string.Format(CultureInfo.InvariantCulture, "{0:E6}", a[i]));
        }
    }
}
