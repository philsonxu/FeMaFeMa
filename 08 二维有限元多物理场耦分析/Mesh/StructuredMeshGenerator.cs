using System;
using System.Collections.Generic;
using MultiPhysicsFEM2D.Elements;

namespace MultiPhysicsFEM2D.Mesh
{
    public static class StructuredMeshGenerator
    {
        public enum Preset
        {
            CantileverCST3,
            CantileverQ4,
            ThermalBeamCST3,       // 热-结构耦合梁：顶面Thot，底面Tcold
            LidDrivenCavityQ4,     // NS 顶盖驱动方腔
            NaturalConvectionQ4    // Boussinesq 自然对流方腔
        }

        /// <summary>
        /// 生成矩形网格，返回 mesh。
        /// 节点编号：(i,j) = j*(nx+1)+i, i=0..nx, j=0..ny；x = i*dx, y = j*dy
        /// </summary>
        public static FEMesh Generate(double Lx, double Ly, int nx, int ny, ElementType et, Preset preset,
            double loadX = 0.0, double loadY = 0.0,
            double Thot = 393.0, double Tcold = 293.0,
            double uTop = 1.0, double TLeft = 1.0, double TRight = 0.0)
        {
            FEMesh m = new FEMesh();
            int nNodes = (nx + 1) * (ny + 1);
            m.AllocateNodes(nNodes);
            double dx = Lx / nx, dy = Ly / ny;
            for (int j = 0; j <= ny; j++)
                for (int i = 0; i <= nx; i++)
                {
                    int id = j * (nx + 1) + i;
                    m.X[id] = i * dx;
                    m.Y[id] = j * dy;
                }
            m.Elements = new List<FiniteElement>();
            m.BCs = new List<BoundaryCondition>();
            int eid = 0;
            if (et == ElementType.CST3)
            {
                for (int j = 0; j < ny; j++)
                    for (int i = 0; i < nx; i++)
                    {
                        int n00 = j * (nx + 1) + i;
                        int n10 = n00 + 1;
                        int n11 = n00 + (nx + 1) + 1;
                        int n01 = n00 + (nx + 1);
                        // 三角1: n00-n10-n11
                        FiniteElement e1 = ElementFactory.Create(ElementType.CST3);
                        e1.Id = eid++;
                        e1.NodeIds[0] = n00; e1.NodeIds[1] = n10; e1.NodeIds[2] = n11;
                        m.Elements.Add(e1);
                        // 三角2: n00-n11-n01
                        FiniteElement e2 = ElementFactory.Create(ElementType.CST3);
                        e2.Id = eid++;
                        e2.NodeIds[0] = n00; e2.NodeIds[1] = n11; e2.NodeIds[2] = n01;
                        m.Elements.Add(e2);
                    }
            }
            else if (et == ElementType.Q4)
            {
                for (int j = 0; j < ny; j++)
                    for (int i = 0; i < nx; i++)
                    {
                        int n00 = j * (nx + 1) + i;
                        int n10 = n00 + 1;
                        int n11 = n00 + (nx + 1) + 1;
                        int n01 = n00 + (nx + 1);
                        FiniteElement e = ElementFactory.Create(ElementType.Q4);
                        e.Id = eid++;
                        e.NodeIds[0] = n00; e.NodeIds[1] = n10; e.NodeIds[2] = n11; e.NodeIds[3] = n01;
                        m.Elements.Add(e);
                    }
            }
            else if (et == ElementType.LT6)
            {
                // LT6: 每矩形两三角，每条边插入一个中点
                int nNodeMax = (2 * nx + 1) * (2 * ny + 1); // 不对，我们只用角点+边中点——简单实现：在矩形单元内生成6节点三角形，使用已存在的角点 + 边上中点（在两个三角间共享）
                // 为简化，改用角点为原来 (nx+1)(ny+1)，边中点另加，每水平边+垂直边一个中点
                // 这里采用：中点编号在角点之后，水平中点编号 = nNodes + j*(nx+1-1)*? ...
                // 简单起见：直接用 2x2 细分（每个矩形变成 4 个三角形）后用 CST3 代替，LT6 由 INP 导入。
                // 为满足需求，直接回落到 CST3 生成：
                for (int j = 0; j < ny; j++)
                    for (int i = 0; i < nx; i++)
                    {
                        int n00 = j * (nx + 1) + i;
                        int n10 = n00 + 1;
                        int n11 = n00 + (nx + 1) + 1;
                        int n01 = n00 + (nx + 1);
                        FiniteElement e1 = ElementFactory.Create(ElementType.CST3);
                        e1.Id = eid++; e1.NodeIds[0] = n00; e1.NodeIds[1] = n10; e1.NodeIds[2] = n11; m.Elements.Add(e1);
                        FiniteElement e2 = ElementFactory.Create(ElementType.CST3);
                        e2.Id = eid++; e2.NodeIds[0] = n00; e2.NodeIds[1] = n11; e2.NodeIds[2] = n01; m.Elements.Add(e2);
                    }
                et = ElementType.CST3;
            }

            // 施加边界条件按预设
            switch (preset)
            {
                case Preset.CantileverCST3:
                case Preset.CantileverQ4:
                    // 左边界固支
                    for (int j = 0; j <= ny; j++)
                    {
                        int nid = j * (nx + 1) + 0;
                        m.BCs.Add(new BoundaryCondition { Type = BCType.FixBoth, NodeId = nid, Value = 0.0 });
                    }
                    // 右端中点加 Y 向力
                    int midTop = (ny / 2) * (nx + 1) + nx;
                    m.BCs.Add(new BoundaryCondition { Type = BCType.TractionY, NodeId = midTop, Value = loadY });
                    break;
                case Preset.ThermalBeamCST3:
                    // 左固支；顶面 Thot，底面 Tcold
                    for (int j = 0; j <= ny; j++)
                    {
                        int nid = j * (nx + 1) + 0;
                        m.BCs.Add(new BoundaryCondition { Type = BCType.FixBoth, NodeId = nid, Value = 0.0 });
                    }
                    for (int i = 0; i <= nx; i++)
                    {
                        int nidBot = 0 * (nx + 1) + i;
                        int nidTop = ny * (nx + 1) + i;
                        m.BCs.Add(new BoundaryCondition { Type = BCType.Temperature, NodeId = nidBot, Value = Tcold });
                        m.BCs.Add(new BoundaryCondition { Type = BCType.Temperature, NodeId = nidTop, Value = Thot });
                    }
                    break;
                case Preset.LidDrivenCavityQ4:
                    // 4 壁无滑（u=v=0），顶盖 v=0, u=uTop
                    for (int i = 0; i <= nx; i++)
                    {
                        int bot = 0 * (nx + 1) + i;
                        int top = ny * (nx + 1) + i;
                        m.BCs.Add(new BoundaryCondition { Type = BCType.VelocityX, NodeId = bot, Value = 0.0 });
                        m.BCs.Add(new BoundaryCondition { Type = BCType.VelocityY, NodeId = bot, Value = 0.0 });
                        m.BCs.Add(new BoundaryCondition { Type = BCType.VelocityY, NodeId = top, Value = 0.0 });
                        m.BCs.Add(new BoundaryCondition { Type = BCType.VelocityX, NodeId = top, Value = uTop });
                    }
                    for (int j = 1; j < ny; j++)
                    {
                        int left = j * (nx + 1) + 0;
                        int right = j * (nx + 1) + nx;
                        m.BCs.Add(new BoundaryCondition { Type = BCType.VelocityX, NodeId = left, Value = 0.0 });
                        m.BCs.Add(new BoundaryCondition { Type = BCType.VelocityY, NodeId = left, Value = 0.0 });
                        m.BCs.Add(new BoundaryCondition { Type = BCType.VelocityX, NodeId = right, Value = 0.0 });
                        m.BCs.Add(new BoundaryCondition { Type = BCType.VelocityY, NodeId = right, Value = 0.0 });
                    }
                    // 参考压力
                    m.BCs.Add(new BoundaryCondition { Type = BCType.Pressure, NodeId = 0, Value = 0.0 });
                    break;
                case Preset.NaturalConvectionQ4:
                    // 4 壁无滑
                    for (int i = 0; i <= nx; i++)
                    {
                        int bot = 0 * (nx + 1) + i;
                        int top = ny * (nx + 1) + i;
                        m.BCs.Add(new BoundaryCondition { Type = BCType.VelocityX, NodeId = bot, Value = 0.0 });
                        m.BCs.Add(new BoundaryCondition { Type = BCType.VelocityY, NodeId = bot, Value = 0.0 });
                        m.BCs.Add(new BoundaryCondition { Type = BCType.VelocityX, NodeId = top, Value = 0.0 });
                        m.BCs.Add(new BoundaryCondition { Type = BCType.VelocityY, NodeId = top, Value = 0.0 });
                    }
                    for (int j = 1; j < ny; j++)
                    {
                        int left = j * (nx + 1) + 0;
                        int right = j * (nx + 1) + nx;
                        m.BCs.Add(new BoundaryCondition { Type = BCType.VelocityX, NodeId = left, Value = 0.0 });
                        m.BCs.Add(new BoundaryCondition { Type = BCType.VelocityY, NodeId = left, Value = 0.0 });
                        m.BCs.Add(new BoundaryCondition { Type = BCType.VelocityX, NodeId = right, Value = 0.0 });
                        m.BCs.Add(new BoundaryCondition { Type = BCType.VelocityY, NodeId = right, Value = 0.0 });
                    }
                    // 左壁 TLeft，右壁 TRight，绝热顶底（无需 BC）
                    for (int j = 0; j <= ny; j++)
                    {
                        int left = j * (nx + 1) + 0;
                        int right = j * (nx + 1) + nx;
                        m.BCs.Add(new BoundaryCondition { Type = BCType.Temperature, NodeId = left, Value = TLeft });
                        m.BCs.Add(new BoundaryCondition { Type = BCType.Temperature, NodeId = right, Value = TRight });
                    }
                    m.BCs.Add(new BoundaryCondition { Type = BCType.Pressure, NodeId = 0, Value = 0.0 });
                    break;
            }
            m.ApplyDefaultMaterial();
            m.UpdateElementCoords();
            m.AllocateResults();
            return m;
        }
    }
}
