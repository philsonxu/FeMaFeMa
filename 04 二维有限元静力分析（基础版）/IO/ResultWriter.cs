using System;
using System.IO;
using System.Globalization;
using System.Text;
using FEM2D.Mesh;
using FEM2D.Solvers;

namespace FEM2D.IO
{
    public static class ResultWriter
    {
        public static void WriteStructuralCsv(FEMesh mesh, StructuralResult result, string outDir)
        {
            Directory.CreateDirectory(outDir);
            string nodeFile = Path.Combine(outDir, "nodes_displacement.csv");
            StringBuilder sb = new StringBuilder();
            sb.AppendLine("NodeId,X,Y,Ux,Uy,DispMagnitude");
            for (int i = 0; i < mesh.NumNodes; i++)
            {
                double x = mesh.X(i);
                double y = mesh.Y(i);
                double ux = result.Displacement[2 * i];
                double uy = result.Displacement[2 * i + 1];
                double m = Math.Sqrt(ux * ux + uy * uy);
                sb.AppendLine($"{i},{x.ToString("G6", CultureInfo.InvariantCulture)},{y.ToString("G6", CultureInfo.InvariantCulture)},{ux:G6},{uy:G6},{m:G6}");
            }
            File.WriteAllText(nodeFile, sb.ToString());

            string elemFile = Path.Combine(outDir, "elements_stress.csv");
            sb = new StringBuilder();
            sb.AppendLine("ElementId,Type,Sxx,Syy,Sxy,VonMises");
            for (int e = 0; e < mesh.NumElements; e++)
            {
                sb.AppendLine($"{e},{mesh.ElementTypes[e]},{result.StressXX[e]:G6},{result.StressYY[e]:G6},{result.StressXY[e]:G6},{result.VonMises[e]:G6}");
            }
            File.WriteAllText(elemFile, sb.ToString());
        }

        public static void WriteNSCsv(FEMesh mesh, NSResult result, string outDir)
        {
            Directory.CreateDirectory(outDir);
            string nodeFile = Path.Combine(outDir, "nodes_ns.csv");
            StringBuilder sb = new StringBuilder();
            sb.AppendLine("NodeId,X,Y,U,V,P,Speed");
            for (int i = 0; i < mesh.NumNodes; i++)
            {
                double x = mesh.X(i);
                double y = mesh.Y(i);
                sb.AppendLine($"{i},{x:G6},{y:G6},{result.U[i]:G6},{result.V[i]:G6},{result.P[i]:G6},{result.Speed[i]:G6}");
            }
            File.WriteAllText(nodeFile, sb.ToString());
        }

        /// <summary>
        /// 写入 VTK Legacy 非结构化网格（ParaView 可读）
        /// </summary>
        public static void WriteVtk(FEMesh mesh, object result, string path, string fieldName = "Structural")
        {
            StringBuilder sb = new StringBuilder();
            sb.AppendLine("# vtk DataFile Version 3.0");
            sb.AppendLine("FEM2D output");
            sb.AppendLine("ASCII");
            sb.AppendLine("DATASET UNSTRUCTURED_GRID");
            sb.AppendLine($"POINTS {mesh.NumNodes} double");
            for (int i = 0; i < mesh.NumNodes; i++)
            {
                sb.AppendLine($"{mesh.X(i).ToString("G6", CultureInfo.InvariantCulture)} {mesh.Y(i).ToString("G6", CultureInfo.InvariantCulture)} 0.0");
            }
            // Cells
            int totalEntries = 0;
            for (int e = 0; e < mesh.NumElements; e++) totalEntries += 1 + mesh.Elements[e].Length;
            sb.AppendLine($"CELLS {mesh.NumElements} {totalEntries}");
            for (int e = 0; e < mesh.NumElements; e++)
            {
                int[] c = mesh.Elements[e];
                sb.Append(c.Length);
                for (int k = 0; k < c.Length; k++) sb.Append(" " + c[k]);
                sb.AppendLine();
            }
            sb.AppendLine($"CELL_TYPES {mesh.NumElements}");
            for (int e = 0; e < mesh.NumElements; e++)
            {
                int[] c = mesh.Elements[e];
                int vtkType;
                switch (c.Length)
                {
                    case 3: vtkType = 5; break;   // VTK_TRIANGLE
                    case 4: vtkType = 9; break;   // VTK_QUAD
                    case 6: vtkType = 22; break;  // VTK_QUADRATIC_TRIANGLE
                    default: vtkType = 7; break;
                }
                sb.AppendLine(vtkType.ToString());
            }

            sb.AppendLine($"POINT_DATA {mesh.NumNodes}");
            if (result is StructuralResult s)
            {
                sb.AppendLine("VECTORS displacement double");
                for (int i = 0; i < mesh.NumNodes; i++)
                {
                    sb.AppendLine($"{s.Displacement[2 * i]:G6} {s.Displacement[2 * i + 1]:G6} 0.0");
                }
                sb.AppendLine("SCALARS disp_mag double 1");
                sb.AppendLine("LOOKUP_TABLE default");
                for (int i = 0; i < mesh.NumNodes; i++)
                {
                    double ux = s.Displacement[2 * i]; double uy = s.Displacement[2 * i + 1];
                    sb.AppendLine($"{Math.Sqrt(ux * ux + uy * uy):G6}");
                }
            }
            else if (result is NSResult ns)
            {
                sb.AppendLine("VECTORS velocity double");
                for (int i = 0; i < mesh.NumNodes; i++)
                {
                    sb.AppendLine($"{ns.U[i]:G6} {ns.V[i]:G6} 0.0");
                }
                sb.AppendLine("SCALARS pressure double 1");
                sb.AppendLine("LOOKUP_TABLE default");
                for (int i = 0; i < mesh.NumNodes; i++) sb.AppendLine($"{ns.P[i]:G6}");
                sb.AppendLine("SCALARS speed double 1");
                sb.AppendLine("LOOKUP_TABLE default");
                for (int i = 0; i < mesh.NumNodes; i++) sb.AppendLine($"{ns.Speed[i]:G6}");
            }

            // Cell data: stress
            if (result is StructuralResult se)
            {
                sb.AppendLine($"CELL_DATA {mesh.NumElements}");
                sb.AppendLine("SCALARS vonMises double 1");
                sb.AppendLine("LOOKUP_TABLE default");
                for (int e = 0; e < mesh.NumElements; e++) sb.AppendLine($"{se.VonMises[e]:G6}");
                sb.AppendLine("SCALARS Sxx double 1");
                sb.AppendLine("LOOKUP_TABLE default");
                for (int e = 0; e < mesh.NumElements; e++) sb.AppendLine($"{se.StressXX[e]:G6}");
                sb.AppendLine("SCALARS Syy double 1");
                sb.AppendLine("LOOKUP_TABLE default");
                for (int e = 0; e < mesh.NumElements; e++) sb.AppendLine($"{se.StressYY[e]:G6}");
                sb.AppendLine("SCALARS Sxy double 1");
                sb.AppendLine("LOOKUP_TABLE default");
                for (int e = 0; e < mesh.NumElements; e++) sb.AppendLine($"{se.StressXY[e]:G6}");
            }
            File.WriteAllText(path, sb.ToString());
        }
    }
}
