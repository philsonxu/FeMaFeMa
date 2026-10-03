using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using NonlinearFEM2D.Elements;

namespace NonlinearFEM2D.Mesh
{
    public static class MeshReader
    {
        public static FEMesh ReadInp(string path)
        {
            string[] lines = File.ReadAllLines(path);
            List<double[]> nodes = new List<double[]>();
            List<int[]> elems = new List<int[]>();
            ElementType eType = ElementType.CST3;
            List<BoundaryCondition> bcs = new List<BoundaryCondition>();
            string section = "";
            bool firstElem = true;
            double E = 2.1e11, nu = 0.3, rho = 7850.0, mu = 1.0e-3, t = 0.01, sy = 2.5e8, h = 2.1e9;
            for (int li = 0; li < lines.Length; li++)
            {
                string raw = lines[li].Trim();
                if (raw.Length == 0 || raw.StartsWith("**") || raw.StartsWith("#")) continue;
                if (raw.StartsWith("*"))
                {
                    section = raw.ToLower();
                    if (section.Contains("cps3") || section.Contains("tri3") || section.Contains("cst")) eType = ElementType.CST3;
                    else if (section.Contains("cps6") || section.Contains("tri6") || section.Contains("lt6")) eType = ElementType.LT6;
                    else if (section.Contains("cpe4") || section.Contains("quad") || section.Contains("q4")) eType = ElementType.Q4;
                    continue;
                }
                string[] parts = raw.Split(new char[] { ',', ' ', '\t' }, StringSplitOptions.RemoveEmptyEntries);
                if (section.Contains("node"))
                {
                    double x = double.Parse(parts[1], CultureInfo.InvariantCulture);
                    double y = double.Parse(parts[2], CultureInfo.InvariantCulture);
                    nodes.Add(new double[] { x, y });
                }
                else if (section.Contains("element"))
                {
                    int npe = (eType == ElementType.CST3) ? 3 : (eType == ElementType.LT6 ? 6 : 4);
                    int[] conn = new int[npe];
                    for (int i = 0; i < npe; i++) conn[i] = int.Parse(parts[i + 1]) - 1;
                    elems.Add(conn);
                }
                else if (section.Contains("boundary") || section.Contains("bc"))
                {
                    // id,dof(1/2),val
                    if (parts.Length >= 3)
                    {
                        int nid = int.Parse(parts[0]) - 1;
                        int dof = int.Parse(parts[1]) - 1;
                        double val = double.Parse(parts[2], CultureInfo.InvariantCulture);
                        bcs.Add(new BoundaryCondition { NodeId = nid, Dof = dof, Value = val });
                    }
                }
                else if (section.Contains("material"))
                {
                    for (int i = 0; i < parts.Length - 1; i += 2)
                    {
                        string k = parts[i].ToLower();
                        double v = double.Parse(parts[i + 1], CultureInfo.InvariantCulture);
                        if (k == "e") E = v; else if (k == "nu") nu = v;
                        else if (k == "rho") rho = v; else if (k == "mu") mu = v;
                        else if (k == "t") t = v; else if (k == "sy") sy = v; else if (k == "h") h = v;
                    }
                }
            }
            FEMesh m = new FEMesh();
            m.ElementType = eType;
            m.Allocate(nodes.Count);
            for (int i = 0; i < nodes.Count; i++) { m.X[i] = nodes[i][0]; m.Y[i] = nodes[i][1]; }
            for (int ie = 0; ie < elems.Count; ie++)
            {
                FiniteElement e = ElementFactory.Create(eType);
                e.Id = ie;
                int[] c = elems[ie];
                int npe = c.Length;
                e.NodeIds = c;
                e.X = new double[npe]; e.Y = new double[npe];
                for (int k = 0; k < npe; k++) { e.X[k] = m.X[c[k]]; e.Y[k] = m.Y[c[k]]; }
                e.YoungModulus = E; e.PoissonRatio = nu; e.Density = rho; e.Viscosity = mu; e.Thickness = t;
                e.Initialize();
                e.Material = new Materials.J2PlasticityMaterial(E, nu, sy, h);
                m.Elements.Add(e);
            }
            foreach (BoundaryCondition bc in bcs) m.Dirichlet.Add(bc);
            return m;
        }
    }
}
