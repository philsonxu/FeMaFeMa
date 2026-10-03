namespace ModalFEM2D.IO
{
    using System;
    using System.Globalization;
    using System.IO;
    using System.Text;
    using ModalFEM2D.Mesh;

    public static class ResultWriter
    {
        public static void WriteCSV(FEMesh mesh, string outDir)
        {
            Directory.CreateDirectory(outDir);
            using (StreamWriter sw = new StreamWriter(Path.Combine(outDir, "nodes.csv"), false, Encoding.UTF8))
            {
                sw.WriteLine("id,x,y,ux,uy,vx,vy,p,sxx,syy,sxy,vm");
                for (int i = 0; i < mesh.NumNodes; i++)
                {
                    sw.WriteLine(string.Format(CultureInfo.InvariantCulture,
                        "{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11}",
                        i, mesh.X[i], mesh.Y[i],
                        mesh.DisplacementU[i], mesh.DisplacementV[i],
                        mesh.VelocityU[i], mesh.VelocityV[i], mesh.Pressure[i],
                        mesh.StressXX[i], mesh.StressYY[i], mesh.StressXY[i], mesh.VonMises[i]));
                }
            }
            using (StreamWriter sw = new StreamWriter(Path.Combine(outDir, "elements.csv"), false, Encoding.UTF8))
            {
                sw.Write("id,type");
                for (int k = 0; k < mesh.NodesPerElement; k++) sw.Write(",n" + k);
                sw.WriteLine();
                string tname = mesh.ElementType.ToString();
                for (int e = 0; e < mesh.NumElements; e++)
                {
                    sw.Write(e); sw.Write(","); sw.Write(tname);
                    for (int k = 0; k < mesh.NodesPerElement; k++) { sw.Write(","); sw.Write(mesh.Connectivity[e, k]); }
                    sw.WriteLine();
                }
            }
            if (mesh.ModalFreqHz != null)
            {
                using (StreamWriter sw = new StreamWriter(Path.Combine(outDir, "modes.txt"), false, Encoding.UTF8))
                {
                    sw.WriteLine("# Mode, f [Hz], omega [rad/s]");
                    for (int m = 0; m < mesh.ModalFreqHz.Length; m++)
                    {
                        double f = mesh.ModalFreqHz[m];
                        sw.WriteLine(string.Format(CultureInfo.InvariantCulture, "{0},{1},{2}", m + 1, f, 2 * Math.PI * f));
                    }
                }
            }
        }

        public static void WriteVTK(FEMesh mesh, string outDir)
        {
            Directory.CreateDirectory(outDir);
            string path = Path.Combine(outDir, "result.vtk");
            using (StreamWriter sw = new StreamWriter(path, false, Encoding.ASCII))
            {
                sw.WriteLine("# vtk DataFile Version 3.0");
                sw.WriteLine("FEM2D result");
                sw.WriteLine("ASCII");
                sw.WriteLine("DATASET UNSTRUCTURED_GRID");
                sw.WriteLine("POINTS " + mesh.NumNodes + " double");
                for (int i = 0; i < mesh.NumNodes; i++)
                    sw.WriteLine(string.Format(CultureInfo.InvariantCulture, "{0} {1} 0", mesh.X[i], mesh.Y[i]));
                int cellType = 5; int vtkNpe = mesh.NodesPerElement;
                if (mesh.ElementType == Elements.ElementType.LT6) { cellType = 22; vtkNpe = 6; }
                else if (mesh.ElementType == Elements.ElementType.Q4) { cellType = 9; vtkNpe = 4; }
                else { cellType = 5; vtkNpe = 3; }
                sw.WriteLine("CELLS " + mesh.NumElements + " " + (mesh.NumElements * (vtkNpe + 1)));
                for (int e = 0; e < mesh.NumElements; e++)
                {
                    sw.Write(vtkNpe);
                    for (int k = 0; k < vtkNpe; k++) sw.Write(" " + mesh.Connectivity[e, k]);
                    sw.WriteLine();
                }
                sw.WriteLine("CELL_TYPES " + mesh.NumElements);
                for (int e = 0; e < mesh.NumElements; e++) sw.WriteLine(cellType);
                sw.WriteLine("POINT_DATA " + mesh.NumNodes);
                WriteScalar(sw, "displacement", Mag(mesh.DisplacementU, mesh.DisplacementV));
                WriteVector(sw, "displacement_vec", mesh.DisplacementU, mesh.DisplacementV);
                WriteScalar(sw, "velocity", Mag(mesh.VelocityU, mesh.VelocityV));
                WriteScalar(sw, "pressure", mesh.Pressure);
                WriteScalar(sw, "stress_xx", mesh.StressXX);
                WriteScalar(sw, "stress_yy", mesh.StressYY);
                WriteScalar(sw, "stress_xy", mesh.StressXY);
                WriteScalar(sw, "von_mises", mesh.VonMises);
                if (mesh.ModalFreqHz != null)
                {
                    for (int m = 0; m < mesh.ModalFreqHz.Length; m++)
                    {
                        double[] ux = new double[mesh.NumNodes]; double[] uy = new double[mesh.NumNodes];
                        for (int i = 0; i < mesh.NumNodes; i++) { ux[i] = mesh.ModalShapes[m, 2 * i]; uy[i] = mesh.ModalShapes[m, 2 * i + 1]; }
                        WriteScalar(sw, "mode" + (m + 1) + "_f=" + mesh.ModalFreqHz[m].ToString("F3", CultureInfo.InvariantCulture) + "Hz", Mag(ux, uy));
                    }
                }
            }
        }

        private static double[] Mag(double[] ux, double[] uy)
        {
            double[] m = new double[ux.Length];
            for (int i = 0; i < ux.Length; i++) m[i] = Math.Sqrt(ux[i] * ux[i] + uy[i] * uy[i]);
            return m;
        }

        private static void WriteScalar(StreamWriter sw, string name, double[] arr)
        {
            sw.WriteLine("SCALARS " + name + " double 1");
            sw.WriteLine("LOOKUP_TABLE default");
            for (int i = 0; i < arr.Length; i++)
                sw.WriteLine(arr[i].ToString("G6", CultureInfo.InvariantCulture));
        }

        private static void WriteVector(StreamWriter sw, string name, double[] ux, double[] uy)
        {
            sw.WriteLine("VECTORS " + name + " double");
            for (int i = 0; i < ux.Length; i++)
                sw.WriteLine(string.Format(CultureInfo.InvariantCulture, "{0} {1} 0", ux[i], uy[i]));
        }
    }
}
