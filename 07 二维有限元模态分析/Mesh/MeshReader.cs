namespace ModalFEM2D.Mesh
{
    using System;
    using System.Collections.Generic;
    using System.Globalization;
    using System.IO;
    using ModalFEM2D.Elements;

    /// <summary>
    /// Reads a simple Abaqus-INP-like text format:
    ///   *HEADING
    ///   *NODE
    ///   id, x, y
    ///   ...
    ///   *ELEMENT, TYPE=CST3 | LT6 | S4 (Q4)
    ///   eid, n1, n2, n3 [, n4, n5, n6]
    ///   ...
    ///   *BOUNDARY           (fixed, ux=uy=0; or ux=val)
    ///   nid, 1, 1, val  (ux=val)
    ///   nid, 2, 2, val  (uy=val)
    ///   nid, 1, 2, val  (ux=uy=val)
    ///   *CLOAD
    ///   nid, dir(1=x,2=y), value
    ///   *MATERIAL
    ///   E, nu, rho, mu, thickness
    /// Node IDs and element IDs in the file are 1-based; they are converted to 0-based on import.
    /// </summary>
    public static class MeshReader
    {
        public static FEMesh Read(string path)
        {
            FEMesh mesh = new FEMesh();
            mesh.Thickness = 1.0; mesh.YoungModulus = 2.1e11; mesh.PoissonRatio = 0.3;
            mesh.Density = 7800.0; mesh.Viscosity = 1.0e-3;
            mesh.ElementType = ElementType.CST3;
            Dictionary<int, int> nodeMap = new Dictionary<int, int>();
            string section = "";
            string[] lines = File.ReadAllLines(path);
            foreach (string rawLine in lines)
            {
                string line = rawLine.Trim();
                if (line.Length == 0 || line.StartsWith("**")) continue;
                if (line.StartsWith("*"))
                {
                    string up = line.ToUpperInvariant();
                    if (up.Contains("NODE")) section = "NODE";
                    else if (up.Contains("ELEMENT"))
                    {
                        section = "ELEMENT";
                        if (up.Contains("CST3") || up.Contains("CPS3") || up.Contains("T2D2")) mesh.ElementType = ElementType.CST3;
                        else if (up.Contains("LT6") || up.Contains("CPS6") || up.Contains("T2D6")) mesh.ElementType = ElementType.LT6;
                        else if (up.Contains("S4") || up.Contains("Q4") || up.Contains("CPS4")) mesh.ElementType = ElementType.Q4;
                    }
                    else if (up.Contains("BOUNDARY")) section = "BOUNDARY";
                    else if (up.Contains("CLOAD") || up.Contains("FORCE")) section = "CLOAD";
                    else if (up.Contains("MATERIAL")) section = "MATERIAL";
                    else if (up.Contains("HEADING") || up.Contains("END")) section = "";
                    continue;
                }
                string[] parts = line.Split(new char[] { ',', ' ', '\t' }, StringSplitOptions.RemoveEmptyEntries);
                switch (section)
                {
                    case "NODE":
                        {
                            int id = int.Parse(parts[0]);
                            double x = double.Parse(parts[1], CultureInfo.InvariantCulture);
                            double y = double.Parse(parts[2], CultureInfo.InvariantCulture);
                            nodeMap[id] = mesh.Nodes.Count;
                            mesh.AddNode(x, y);
                        }
                        break;
                    case "ELEMENT":
                        {
                            int npe = (mesh.ElementType == ElementType.CST3) ? 3 : (mesh.ElementType == ElementType.Q4 ? 4 : 6);
                            int[] nodes = new int[npe];
                            for (int k = 0; k < npe; k++)
                            {
                                int fid = int.Parse(parts[k + 1]);
                                nodes[k] = nodeMap[fid];
                            }
                            mesh.AddElement(nodes);
                        }
                        break;
                    case "BOUNDARY":
                        {
                            int nid = nodeMap[int.Parse(parts[0])];
                            int d1 = int.Parse(parts[1]);
                            int d2 = (parts.Length >= 3) ? int.Parse(parts[2]) : d1;
                            double val = (parts.Length >= 4) ? double.Parse(parts[3], CultureInfo.InvariantCulture) : 0.0;
                            if (d1 <= 1 && d2 >= 1) mesh.AddBC(nid, BCKind.FixedX, val);
                            if (d1 <= 2 && d2 >= 2) mesh.AddBC(nid, BCKind.FixedY, val);
                        }
                        break;
                    case "CLOAD":
                        {
                            int nid = nodeMap[int.Parse(parts[0])];
                            int dir = int.Parse(parts[1]);
                            double v = double.Parse(parts[2], CultureInfo.InvariantCulture);
                            mesh.AddBC(nid, dir == 1 ? BCKind.ForceX : BCKind.ForceY, v);
                        }
                        break;
                    case "MATERIAL":
                        if (parts[0].ToUpperInvariant().StartsWith("E") && parts.Length >= 6)
                        {
                            mesh.YoungModulus = double.Parse(parts[1], CultureInfo.InvariantCulture);
                            mesh.PoissonRatio = double.Parse(parts[2], CultureInfo.InvariantCulture);
                            mesh.Density = double.Parse(parts[3], CultureInfo.InvariantCulture);
                            mesh.Viscosity = double.Parse(parts[4], CultureInfo.InvariantCulture);
                            mesh.Thickness = double.Parse(parts[5], CultureInfo.InvariantCulture);
                        }
                        break;
                }
            }
            mesh.BuildConnectivity();
            return mesh;
        }
    }
}
