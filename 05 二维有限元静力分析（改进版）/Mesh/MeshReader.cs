// MeshReader.cs - 类 Abaqus INP 文本网格导入器
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using FEM2D.Elements;

namespace FEM2D.Mesh
{
    public static class MeshReader
    {
        public static FEMesh Read(string path)
        {
            FEMesh mesh = new FEMesh();
            List<double> xs = new List<double>();
            List<double> ys = new List<double>();
            Dictionary<int, int> nodeIdMap = new Dictionary<int, int>();

            string[] lines = File.ReadAllLines(path);
            string section = "";
            ElementType curType = ElementType.CST3;
            double E = 2.1e11, nu = 0.3, rho = 7800, mu = 1.0e-3, thick = 1.0;

            foreach (string rawLine in lines)
            {
                string line = rawLine.Trim();
                if (line.Length == 0 || line.StartsWith("**")) continue;
                if (line.StartsWith("*"))
                {
                    string up = line.ToUpper();
                    if (up.StartsWith("*NODE")) section = "NODE";
                    else if (up.StartsWith("*ELEMENT"))
                    {
                        section = "ELEMENT";
                        if (up.Contains("TYPE=CST3") || up.Contains("TYPE=S3")) curType = ElementType.CST3;
                        else if (up.Contains("TYPE=LT6") || up.Contains("TYPE=S6")) curType = ElementType.LT6;
                        else if (up.Contains("TYPE=Q4") || up.Contains("TYPE=S4")) curType = ElementType.Q4;
                        else curType = ElementType.CST3;
                    }
                    else if (up.StartsWith("*BOUNDARY")) section = "BC";
                    else if (up.StartsWith("*CLOAD")) section = "CLOAD";
                    else if (up.StartsWith("*VELOCITY")) section = "VELO";
                    else if (up.StartsWith("*MATERIAL"))
                    {
                        section = "MAT";
                        string[] parts = up.Split(new char[] { ',', ' ' }, StringSplitOptions.RemoveEmptyEntries);
                        foreach (string p in parts)
                        {
                            if (p.StartsWith("E=")) double.TryParse(p.Substring(2), NumberStyles.Float, CultureInfo.InvariantCulture, out E);
                            else if (p.StartsWith("NU=")) double.TryParse(p.Substring(3), NumberStyles.Float, CultureInfo.InvariantCulture, out nu);
                            else if (p.StartsWith("RHO=")) double.TryParse(p.Substring(4), NumberStyles.Float, CultureInfo.InvariantCulture, out rho);
                            else if (p.StartsWith("MU=")) double.TryParse(p.Substring(3), NumberStyles.Float, CultureInfo.InvariantCulture, out mu);
                            else if (p.StartsWith("T=")) double.TryParse(p.Substring(2), NumberStyles.Float, CultureInfo.InvariantCulture, out thick);
                        }
                    }
                    else section = "";
                    continue;
                }

                string[] toks = line.Split(new char[] { ',', ' ', '\t' }, StringSplitOptions.RemoveEmptyEntries);
                if (toks.Length == 0) continue;

                if (section == "NODE")
                {
                    int id = int.Parse(toks[0]);
                    double x = double.Parse(toks[1], CultureInfo.InvariantCulture);
                    double y = double.Parse(toks[2], CultureInfo.InvariantCulture);
                    nodeIdMap[id] = xs.Count;
                    xs.Add(x); ys.Add(y);
                }
                else if (section == "ELEMENT")
                {
                    int eid = int.Parse(toks[0]);
                    FiniteElement e = CreateElement(curType);
                    e.Id = eid;
                    int nen = e.NodesPerElement;
                    e.NodeIds = new int[nen];
                    for (int k = 0; k < nen; k++)
                    {
                        int gid = int.Parse(toks[1 + k]);
                        e.NodeIds[k] = nodeIdMap[gid];
                    }
                    e.YoungModulus = E; e.PoissonRatio = nu; e.Density = rho; e.Viscosity = mu; e.Thickness = thick;
                    mesh.Elements.Add(e);
                }
                else if (section == "BC")
                {
                    int nid = int.Parse(toks[0]);
                    int dof = int.Parse(toks[1]) - 1;
                    double val = double.Parse(toks[2], CultureInfo.InvariantCulture);
                    mesh.DirichletBCs.Add(Tuple.Create(nodeIdMap[nid], dof, val));
                }
                else if (section == "CLOAD")
                {
                    int nid = int.Parse(toks[0]);
                    int dof = int.Parse(toks[1]) - 1;
                    double val = double.Parse(toks[2], CultureInfo.InvariantCulture);
                    mesh.NodalForces.Add(Tuple.Create(nodeIdMap[nid], dof, val));
                }
                else if (section == "VELO")
                {
                    int nid = int.Parse(toks[0]);
                    int dof = int.Parse(toks[1]) - 1;
                    double val = double.Parse(toks[2], CultureInfo.InvariantCulture);
                    mesh.VelocityBCs.Add(Tuple.Create(nodeIdMap[nid], dof, val));
                }
            }

            mesh.NumNodes = xs.Count;
            mesh.NumElements = mesh.Elements.Count;
            mesh.X = xs.ToArray();
            mesh.Y = ys.ToArray();
            mesh.AssignNodeCoordinatesToElements();
            mesh.UpdateBounds();
            mesh.BuildConnectivity();
            return mesh;
        }

        private static FiniteElement CreateElement(ElementType t)
        {
            switch (t)
            {
                case ElementType.CST3: return new CST3Element();
                case ElementType.LT6: return new LT6Element();
                case ElementType.Q4: return new Q4Element();
                default: return new CST3Element();
            }
        }
    }
}
