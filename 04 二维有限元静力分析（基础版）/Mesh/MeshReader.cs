using System;
using System.Globalization;
using System.IO;
using System.Collections.Generic;
using FEM2D.Elements;

namespace FEM2D.Mesh
{
    /// <summary>
    /// 类 Abaqus INP 格式读取器
    /// 支持：
    ///   *NODE
    ///   编号, x, y
    ///   *ELEMENT, TYPE=CPS3|CPS6|CPS4
    ///   编号, n1, n2, ...
    ///   *BOUNDARY
    ///   node, dof(1/2), value
    ///   *CLOAD
    ///   node, dof(1/2), value
    ///   *DSLOAD
    ///   element, face, pressure
    /// 节点编号与单元编号从1开始。
    /// </summary>
    public static class MeshReader
    {
        public static FEMesh ReadInp(string path)
        {
            FEMesh mesh = new FEMesh();
            Dictionary<int, int> nodeMap = new Dictionary<int, int>();
            Dictionary<int, int> elemMap = new Dictionary<int, int>();

            string[] lines = File.ReadAllLines(path);
            string section = "";
            ElementType currentElemType = ElementType.CST3;

            foreach (string rawLine in lines)
            {
                string line = rawLine.Trim();
                if (string.IsNullOrEmpty(line) || line.StartsWith("**")) continue;

                if (line.StartsWith("*"))
                {
                    string up = line.ToUpper();
                    if (up.StartsWith("*NODE")) { section = "NODE"; continue; }
                    if (up.StartsWith("*ELEMENT"))
                    {
                        section = "ELEMENT";
                        if (up.Contains("CPS3") || up.Contains("T3")) currentElemType = ElementType.CST3;
                        else if (up.Contains("CPS6") || up.Contains("T6")) currentElemType = ElementType.LT6;
                        else if (up.Contains("CPS4") || up.Contains("Q4")) currentElemType = ElementType.Q4;
                        else currentElemType = ElementType.CST3;
                        continue;
                    }
                    if (up.StartsWith("*BOUNDARY")) { section = "BOUNDARY"; continue; }
                    if (up.StartsWith("*CLOAD")) { section = "CLOAD"; continue; }
                    if (up.StartsWith("*DSLOAD") || up.StartsWith("*DLOAD")) { section = "DSLOAD"; continue; }
                    if (up.StartsWith("*MATERIAL") || up.StartsWith("*ELASTIC") || up.StartsWith("*STEP"))
                    { section = "OTHER"; continue; }
                    section = "OTHER";
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
                            mesh.Nodes.Add(new double[] { x, y });
                        }
                        break;
                    case "ELEMENT":
                        {
                            int id = int.Parse(parts[0]);
                            int nn;
                            if (currentElemType == ElementType.CST3) nn = 3;
                            else if (currentElemType == ElementType.LT6) nn = 6;
                            else nn = 4;
                            int[] conn = new int[nn];
                            for (int k = 0; k < nn; k++)
                            {
                                int gid = int.Parse(parts[k + 1]);
                                conn[k] = nodeMap[gid];
                            }
                            elemMap[id] = mesh.Elements.Count;
                            mesh.Elements.Add(conn);
                            mesh.ElementTypes.Add(currentElemType);
                            mesh.MaterialIDs.Add(0);
                        }
                        break;
                    case "BOUNDARY":
                        {
                            int gid = int.Parse(parts[0]);
                            int dof = int.Parse(parts[1]) - 1;
                            double val = parts.Length >= 3 ? double.Parse(parts[2], CultureInfo.InvariantCulture) : 0.0;
                            mesh.DirichletBCs.Add((nodeMap[gid], dof, val));
                        }
                        break;
                    case "CLOAD":
                        {
                            int gid = int.Parse(parts[0]);
                            int dof = int.Parse(parts[1]) - 1;
                            double val = double.Parse(parts[2], CultureInfo.InvariantCulture);
                            mesh.ConcentratedForces.Add((nodeMap[gid], dof, val));
                        }
                        break;
                    case "DSLOAD":
                        {
                            // 简化：elemId, faceLabel, pressure -> 按边均分到两端节点
                            int eid = int.Parse(parts[0]);
                            string face = parts[1];
                            double p = double.Parse(parts[2], CultureInfo.InvariantCulture);
                            int eLocal = elemMap[eid];
                            int[] conn = mesh.Elements[eLocal];
                            ElementType et = mesh.ElementTypes[eLocal];
                            int npe = conn.Length;
                            int fIdx = ParseFaceIndex(face);
                            int a = conn[fIdx % npe];
                            int b = conn[(fIdx + 1) % npe];
                            double dx = mesh.X(b) - mesh.X(a);
                            double dy = mesh.Y(b) - mesh.Y(a);
                            double len = Math.Sqrt(dx * dx + dy * dy);
                            double nx = -dy / len;
                            double ny = dx / len;
                            // 合力 p*len，均分到两端
                            double fx = p * nx * len * 0.5;
                            double fy = p * ny * len * 0.5;
                            mesh.ConcentratedForces.Add((a, 0, fx));
                            mesh.ConcentratedForces.Add((a, 1, fy));
                            mesh.ConcentratedForces.Add((b, 0, fx));
                            mesh.ConcentratedForces.Add((b, 1, fy));
                        }
                        break;
                }
            }
            // 自动判断问题类型：若存在 dof=3 的约束（压力），视为 NS
            bool hasPDof = false;
            foreach ((int node, int dof, double val) in mesh.DirichletBCs)
            {
                if (dof >= 2) { hasPDof = true; break; }
            }
            if (hasPDof) mesh.Metadata["ProblemType"] = "NS";
            return mesh;
        }

        private static int ParseFaceIndex(string face)
        {
            if (string.IsNullOrEmpty(face)) return 0;
            string s = face.Trim().ToUpper();
            if (s.StartsWith("S") || s.StartsWith("E") || s.StartsWith("F") || s.StartsWith("P"))
            {
                s = s.Substring(1);
            }
            int v;
            if (int.TryParse(s, out v)) return Math.Max(0, v - 1);
            return 0;
        }
    }
}
