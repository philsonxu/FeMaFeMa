using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using MultiPhysicsFEM2D.Elements;

namespace MultiPhysicsFEM2D.Mesh
{
    public static class MeshReader
    {
        /// <summary>
        /// 读取类 Abaqus INP 文件。支持以下关键字：
        /// *NODE     节点: id, x, y
        /// *ELEMENT, TYPE=CPS3/CPS6/CPS4   单元: id, n1, n2, ... (1-based)
        /// *BOUNDARY 节点集边界: nid, dof(1=x/2=y/11=T), value  （或 *BOUNDARY node,1,1,0.0 固定）
        /// *CLOAD, *DFLUX, *TEMPERATURE 等简单荷载
        ///
        /// 简单约定（易于手写）：每行一个边界条件：
        ///   节点号, BCType, value
        ///   BCType 字符串：FixX FixY FixBoth TractionX TractionY Temp HeatFlux Vx Vy P
        /// </summary>
        public static FEMesh ReadINP(string path)
        {
            FEMesh mesh = new FEMesh();
            List<string[]> nodeLines = new List<string[]>();
            List<Tuple<ElementType, int[]>> elemDefs = new List<Tuple<ElementType, int[]>>();
            List<string[]> bcLines = new List<string[]>();

            string[] lines = File.ReadAllLines(path);
            string section = "";
            for (int i = 0; i < lines.Length; i++)
            {
                string line = lines[i].Trim();
                if (line.Length == 0 || line.StartsWith("**")) continue;
                if (line.StartsWith("*"))
                {
                    string up = line.ToUpperInvariant();
                    if (up.StartsWith("*NODE")) section = "NODE";
                    else if (up.StartsWith("*ELEMENT"))
                    {
                        section = "ELEMENT";
                        if (up.Contains("CPS3") || up.Contains("CST3")) elemDefs.Clear(); // just marker
                    }
                    else if (up.StartsWith("*BOUNDARY") || up.StartsWith("*BC")) section = "BC";
                    else if (up.StartsWith("*CLOAD") || up.StartsWith("*DSLOAD")) section = "BC";
                    else if (up.StartsWith("*TEMPERATURE") || up.StartsWith("*TEMP")) section = "BC";
                    else if (up.StartsWith("*MATERIAL")) section = "MAT";
                    else if (up.StartsWith("*")) section = "";
                    continue;
                }
                string[] parts = line.Split(new char[] { ',', ' ', '\t' }, StringSplitOptions.RemoveEmptyEntries);
                if (parts.Length == 0) continue;
                if (section == "NODE") nodeLines.Add(parts);
                else if (section == "ELEMENT")
                {
                    ElementType et = ElementType.CST3;
                    int start = 1;
                    // 判断关键字是否含 CPS6/CPS4
                    if (lines[i - 1].ToUpperInvariant().Contains("CPS6") || lines[i - 1].ToUpperInvariant().Contains("LT6")) et = ElementType.LT6;
                    else if (lines[i - 1].ToUpperInvariant().Contains("CPS4") || lines[i - 1].ToUpperInvariant().Contains("Q4")) et = ElementType.Q4;
                    else et = ElementType.CST3;
                    int nNode = et == ElementType.CST3 ? 3 : (et == ElementType.LT6 ? 6 : 4);
                    int[] nodes = new int[nNode];
                    for (int k = 0; k < nNode; k++) nodes[k] = int.Parse(parts[start + k]) - 1; // 0-based
                    elemDefs.Add(Tuple.Create(et, nodes));
                }
                else if (section == "BC") bcLines.Add(parts);
            }

            // 建节点
            int maxId = -1;
            foreach (string[] p in nodeLines) { int id = int.Parse(p[0]) - 1; if (id > maxId) maxId = id; }
            mesh.AllocateNodes(maxId + 1);
            Dictionary<int, int> nodeRemap = new Dictionary<int, int>();
            foreach (string[] p in nodeLines)
            {
                int id = int.Parse(p[0]) - 1;
                double x = double.Parse(p[1], CultureInfo.InvariantCulture);
                double y = double.Parse(p[2], CultureInfo.InvariantCulture);
                mesh.X[id] = x; mesh.Y[id] = y;
                nodeRemap[id] = id;
            }

            // 建单元
            int eid = 0;
            foreach (Tuple<ElementType, int[]> ed in elemDefs)
            {
                FiniteElement e = ElementFactory.Create(ed.Item1);
                e.Id = eid++;
                int[] nids = ed.Item2;
                for (int k = 0; k < nids.Length; k++) e.NodeIds[k] = nids[k];
                mesh.Elements.Add(e);
            }

            // 边界条件
            foreach (string[] p in bcLines)
            {
                int nid;
                BCType bt;
                double val = 0.0;
                if (p.Length >= 3)
                {
                    nid = int.Parse(p[0]) - 1;
                    string s = p[1];
                    switch (s.ToUpperInvariant())
                    {
                        case "FIXX": case "UX": case "X": bt = BCType.FixX; break;
                        case "FIXY": case "UY": case "Y": bt = BCType.FixY; break;
                        case "FIX": case "FIXBOTH": case "ENCASTRE": bt = BCType.FixBoth; break;
                        case "FORCEX": case "FX": bt = BCType.TractionX; break;
                        case "FORCEY": case "FY": bt = BCType.TractionY; break;
                        case "TEMP": case "T": bt = BCType.Temperature; break;
                        case "HEATFLUX": case "Q": bt = BCType.HeatFlux; break;
                        case "VX": case "VELOX": bt = BCType.VelocityX; break;
                        case "VY": case "VELOY": bt = BCType.VelocityY; break;
                        case "P": case "PRESSURE": bt = BCType.Pressure; break;
                        default: bt = BCType.FixX; break;
                    }
                    val = double.Parse(p[2], CultureInfo.InvariantCulture);
                }
                else continue;
                mesh.BCs.Add(new BoundaryCondition { Type = bt, NodeId = nid, Value = val });
            }
            mesh.ApplyDefaultMaterial();
            mesh.UpdateElementCoords();
            mesh.AllocateResults();
            return mesh;
        }
    }
}
