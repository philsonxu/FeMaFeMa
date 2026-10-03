// ResultWriter.cs - 结果输出（CSV / VTK Legacy / PNG 云图）
using System;
using System.Globalization;
using System.IO;
using System.Text;
using FEM2D.Elements;
using FEM2D.Mesh;
using FEM2D.Rendering;

namespace FEM2D.IO
{
    public static class ResultWriter
    {
        public static void WriteAll(FEMesh mesh, string outputDir, string caseName, string fieldName)
        {
            if (!Directory.Exists(outputDir)) Directory.CreateDirectory(outputDir);
            string nodesCsv = Path.Combine(outputDir, caseName + "_nodes.csv");
            string elemsCsv = Path.Combine(outputDir, caseName + "_elements.csv");
            string vtk = Path.Combine(outputDir, caseName + ".vtk");
            string png = Path.Combine(outputDir, caseName + "_" + SafeField(fieldName) + ".png");
            WriteNodeCSV(mesh, nodesCsv);
            WriteElementCSV(mesh, elemsCsv);
            WriteVTK(mesh, vtk);
            ContourRenderer.RenderToFile(mesh, png, 1600, 1200, fieldName, true);
        }

        private static string SafeField(string f)
        {
            string s = (f ?? "mesh").ToLowerInvariant();
            foreach (char c in Path.GetInvalidFileNameChars()) s = s.Replace(c, '_');
            return s.Replace(' ', '_');
        }

        public static void WriteNodeCSV(FEMesh mesh, string path)
        {
            StringBuilder sb = new StringBuilder();
            sb.AppendLine("id,x,y,disp_u,disp_v,vel_u,vel_v,pressure");
            for (int i = 0; i < mesh.NumNodes; i++)
                sb.Append(i.ToString(CultureInfo.InvariantCulture)).Append(',')
                    .Append(mesh.X[i].ToString("G6", CultureInfo.InvariantCulture)).Append(',')
                    .Append(mesh.Y[i].ToString("G6", CultureInfo.InvariantCulture)).Append(',')
                    .Append(mesh.DisplacementU[i].ToString("G6", CultureInfo.InvariantCulture)).Append(',')
                    .Append(mesh.DisplacementV[i].ToString("G6", CultureInfo.InvariantCulture)).Append(',')
                    .Append(mesh.VelocityU[i].ToString("G6", CultureInfo.InvariantCulture)).Append(',')
                    .Append(mesh.VelocityV[i].ToString("G6", CultureInfo.InvariantCulture)).Append(',')
                    .Append(mesh.Pressure[i].ToString("G6", CultureInfo.InvariantCulture)).AppendLine();
            File.WriteAllText(path, sb.ToString());
        }

        public static void WriteElementCSV(FEMesh mesh, string path)
        {
            StringBuilder sb = new StringBuilder();
            sb.Append("id,type,node_ids,sxx,syy,sxy,von_mises");
            for (int k = 0; k < FEMesh.NodesPerElement(mesh.ElementType); k++) sb.Append(",node").Append(k);
            sb.AppendLine();
            for (int e = 0; e < mesh.NumElements; e++)
            {
                sb.Append(e.ToString(CultureInfo.InvariantCulture)).Append(',')
                    .Append(mesh.ElementType.ToString()).Append(',');
                StringBuilder ids = new StringBuilder();
                int nen = FEMesh.NodesPerElement(mesh.ElementType);
                for (int k = 0; k < nen; k++)
                {
                    if (k > 0) ids.Append(';');
                    ids.Append(mesh.Connectivity[e, k].ToString(CultureInfo.InvariantCulture));
                }
                sb.Append(ids.ToString()).Append(',')
                    .Append(mesh.StressXX[e].ToString("G6", CultureInfo.InvariantCulture)).Append(',')
                    .Append(mesh.StressYY[e].ToString("G6", CultureInfo.InvariantCulture)).Append(',')
                    .Append(mesh.StressXY[e].ToString("G6", CultureInfo.InvariantCulture)).Append(',')
                    .Append(mesh.VonMises[e].ToString("G6", CultureInfo.InvariantCulture));
                for (int k = 0; k < nen; k++)
                    sb.Append(',').Append(mesh.Connectivity[e, k].ToString(CultureInfo.InvariantCulture));
                sb.AppendLine();
            }
            File.WriteAllText(path, sb.ToString());
        }

        public static void WriteVTK(FEMesh mesh, string path)
        {
            StringBuilder sb = new StringBuilder();
            sb.AppendLine("# vtk DataFile Version 3.0");
            sb.AppendLine("FEM2D results");
            sb.AppendLine("ASCII");
            sb.AppendLine("DATASET UNSTRUCTURED_GRID");
            sb.Append("POINTS ").Append(mesh.NumNodes).AppendLine(" double");
            for (int i = 0; i < mesh.NumNodes; i++)
                sb.Append(mesh.X[i].ToString("G6", CultureInfo.InvariantCulture)).Append(' ')
                    .Append(mesh.Y[i].ToString("G6", CultureInfo.InvariantCulture)).AppendLine(" 0.0");
            int nen = FEMesh.NodesPerElement(mesh.ElementType);
            int vtkType;
            switch (mesh.ElementType)
            {
                case ElementType.CST3: vtkType = 5; break;   // VTK_TRIANGLE
                case ElementType.Q4: vtkType = 9; break;     // VTK_QUAD
                case ElementType.LT6: vtkType = 22; break;   // VTK_QUADRATIC_TRIANGLE
                default: vtkType = 5; break;
            }
            sb.Append("CELLS ").Append(mesh.NumElements).Append(' ')
                .AppendLine((mesh.NumElements * (nen + 1)).ToString(CultureInfo.InvariantCulture));
            for (int e = 0; e < mesh.NumElements; e++)
            {
                sb.Append(nen.ToString(CultureInfo.InvariantCulture));
                for (int k = 0; k < nen; k++) sb.Append(' ').Append(mesh.Connectivity[e, k]);
                sb.AppendLine();
            }
            sb.Append("CELL_TYPES ").AppendLine(mesh.NumElements.ToString(CultureInfo.InvariantCulture));
            for (int e = 0; e < mesh.NumElements; e++) sb.AppendLine(vtkType.ToString(CultureInfo.InvariantCulture));
            sb.Append("POINT_DATA ").AppendLine(mesh.NumNodes.ToString(CultureInfo.InvariantCulture));
            sb.AppendLine("VECTORS displacement double");
            for (int i = 0; i < mesh.NumNodes; i++)
                sb.Append(mesh.DisplacementU[i].ToString("G6", CultureInfo.InvariantCulture)).Append(' ')
                    .Append(mesh.DisplacementV[i].ToString("G6", CultureInfo.InvariantCulture)).AppendLine(" 0.0");
            sb.AppendLine("VECTORS velocity double");
            for (int i = 0; i < mesh.NumNodes; i++)
                sb.Append(mesh.VelocityU[i].ToString("G6", CultureInfo.InvariantCulture)).Append(' ')
                    .Append(mesh.VelocityV[i].ToString("G6", CultureInfo.InvariantCulture)).AppendLine(" 0.0");
            sb.AppendLine("SCALARS pressure double 1"); sb.AppendLine("LOOKUP_TABLE default");
            for (int i = 0; i < mesh.NumNodes; i++) sb.AppendLine(mesh.Pressure[i].ToString("G6", CultureInfo.InvariantCulture));
            sb.Append("CELL_DATA ").AppendLine(mesh.NumElements.ToString(CultureInfo.InvariantCulture));
            sb.AppendLine("SCALARS vonMises double 1"); sb.AppendLine("LOOKUP_TABLE default");
            for (int e = 0; e < mesh.NumElements; e++) sb.AppendLine(mesh.VonMises[e].ToString("G6", CultureInfo.InvariantCulture));
            sb.AppendLine("SCALARS stressXX double 1"); sb.AppendLine("LOOKUP_TABLE default");
            for (int e = 0; e < mesh.NumElements; e++) sb.AppendLine(mesh.StressXX[e].ToString("G6", CultureInfo.InvariantCulture));
            sb.AppendLine("SCALARS stressYY double 1"); sb.AppendLine("LOOKUP_TABLE default");
            for (int e = 0; e < mesh.NumElements; e++) sb.AppendLine(mesh.StressYY[e].ToString("G6", CultureInfo.InvariantCulture));
            sb.AppendLine("SCALARS stressXY double 1"); sb.AppendLine("LOOKUP_TABLE default");
            for (int e = 0; e < mesh.NumElements; e++) sb.AppendLine(mesh.StressXY[e].ToString("G6", CultureInfo.InvariantCulture));
            File.WriteAllText(path, sb.ToString());
        }
    }
}
