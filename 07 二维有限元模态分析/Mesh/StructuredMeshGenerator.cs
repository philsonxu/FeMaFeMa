namespace ModalFEM2D.Mesh
{
    using System;
    using System.Collections.Generic;
    using ModalFEM2D.Elements;

    public enum MeshPreset
    {
        CantileverBeam,
        SimplySupportedBeam,
        LidDrivenCavity
    }

    public static class StructuredMeshGenerator
    {
        /// <summary>Generate a cantilever beam: left edge fixed (x=0), tip (x=L) loaded in -y direction.</summary>
        public static FEMesh CantileverBeam(double L, double H, int Nx, int Ny, ElementType type,
                                            double E, double nu, double rho, double thick, double tipLoad)
        {
            FEMesh m = GenerateRect(L, H, Nx, Ny, type);
            m.YoungModulus = E; m.PoissonRatio = nu; m.Density = rho; m.Thickness = thick; m.Viscosity = 1e-3;
            m.ResetResults();
            m.BCs.Clear();
            // left edge fixed
            for (int j = 0; j <= Ny; j++)
            {
                int nid = NodeIndex(0, j, Nx, Ny, type);
                m.AddBC(nid, BCKind.Fixed, 0.0);
            }
            // tip load: concentrated force at middle right node
            int tip = NodeIndex(Nx, Ny / 2, Nx, Ny, type);
            m.AddBC(tip, BCKind.ForceY, -tipLoad);
            // For LT6, tip middle node may be a corner (if Ny even). Just add force to that; otherwise to top/bottom average.
            m.BuildConnectivity();
            return m;
        }

        /// <summary>Simply supported beam: y=0 at x=0 (both dir) and x=L (y only), uniform transverse load. Used for modal demo.</summary>
        public static FEMesh SimplySupportedBeam(double L, double H, int Nx, int Ny, ElementType type,
                                                 double E, double nu, double rho, double thick, double udl)
        {
            FEMesh m = GenerateRect(L, H, Nx, Ny, type);
            m.YoungModulus = E; m.PoissonRatio = nu; m.Density = rho; m.Thickness = thick; m.Viscosity = 1e-3;
            m.ResetResults(); m.BCs.Clear();
            int nLeftBot = NodeIndex(0, 0, Nx, Ny, type);
            int nRightBot = NodeIndex(Nx, 0, Nx, Ny, type);
            m.AddBC(nLeftBot, BCKind.Fixed, 0.0);
            m.AddBC(nRightBot, BCKind.FixedY, 0.0);
            // uniform load on top edge (y=H)
            List<int> top = GetEdgeNodes(Nx, Ny, 2, type); // top edge = eta=1
            foreach (int n in top) m.AddBC(n, BCKind.ForceY, -udl);
            m.BuildConnectivity();
            return m;
        }

        /// <summary>Lid-driven cavity: Q4 grid on [0,L]^2, all walls u=v=0 except top u=LidVel.</summary>
        public static FEMesh LidDrivenCavity(double L, int N, ElementType type,
                                             double rho, double mu, double lidVel)
        {
            FEMesh m = GenerateRect(L, L, N, N, type);
            m.YoungModulus = 1e9; m.PoissonRatio = 0.3; m.Density = rho; m.Viscosity = mu; m.Thickness = 1.0;
            m.ResetResults(); m.BCs.Clear();
            List<int> bot = GetEdgeNodes(N, N, 0, type);
            List<int> lft = GetEdgeNodes(N, N, 3, type);
            List<int> rgt = GetEdgeNodes(N, N, 1, type);
            List<int> top = GetEdgeNodes(N, N, 2, type);
            foreach (int n in bot) { m.AddBC(n, BCKind.VelInX, 0.0); m.AddBC(n, BCKind.VelInY, 0.0); }
            foreach (int n in lft) { m.AddBC(n, BCKind.VelInX, 0.0); m.AddBC(n, BCKind.VelInY, 0.0); }
            foreach (int n in rgt) { m.AddBC(n, BCKind.VelInX, 0.0); m.AddBC(n, BCKind.VelInY, 0.0); }
            foreach (int n in top) { m.AddBC(n, BCKind.VelInX, lidVel); m.AddBC(n, BCKind.VelInY, 0.0); }
            // pressure reference at (0,0)
            m.AddBC(0, BCKind.Pressure, 0.0);
            m.BuildConnectivity();
            return m;
        }

        // ================ rectangular grid generation ================
        private static FEMesh GenerateRect(double Lx, double Ly, int Nx, int Ny, ElementType type)
        {
            FEMesh m = new FEMesh();
            m.ElementType = type;
            double dx = Lx / Nx; double dy = Ly / Ny;
            if (type == ElementType.LT6)
            {
                // P2 mesh: vertices (i*dx, j*dy) + mid-edge nodes. Each cell is split into two P2 triangles.
                // Build a regular P2 grid by adding midpoints to edges of Nx*Ny quads -> 2 P2 triangles per quad.
                int vx = Nx + 1, vy = Ny + 1;
                // Vertex indices: (i,j) -> i*vy + j
                Func<int, int, int> vIndex = (i, j) => i * vy + j;
                for (int i = 0; i <= Nx; i++) for (int j = 0; j <= Ny; j++) m.AddNode(i * dx, j * dy);
                // horizontal mid-edge nodes (i+0.5, j): (vx*vy) + i*vy + j  for i=0..Nx-1, j=0..Ny
                int hMidStart = vx * vy;
                for (int i = 0; i < Nx; i++) for (int j = 0; j <= Ny; j++) m.AddNode((i + 0.5) * dx, j * dy);
                // vertical mid-edge nodes (i, j+0.5): hMidStart + Nx*vy + i*(Ny+1) + j
                int vMidStart = hMidStart + Nx * vy;
                for (int i = 0; i <= Nx; i++) for (int j = 0; j < Ny; j++) m.AddNode(i * dx, (j + 0.5) * dy);
                // Diagonal mid-edge nodes inside each quad (i+0.5, j+0.5) where diagonal goes from (i,j+1)-(i+1,j):
                // we split quad into two LT6: triangle A = (i,j),(i+1,j),(i,j+1); triangle B = (i+1,j+1),(i,j+1),(i+1,j)
                // Diagonal edge of triangle A is (i+1,j)-(i,j+1): midpoint at (i+0.5,j+0.5)
                // Diagonal edge of triangle B is (i,j+1)-(i+1,j): same midpoint; both triangles share that midpoint -> one diagonal node per quad.
                int dMidStart = vMidStart + (Nx + 1) * Ny;
                for (int i = 0; i < Nx; i++) for (int j = 0; j < Ny; j++) m.AddNode((i + 0.5) * dx, (j + 0.5) * dy);
                // Build 2 LT6 elements per quad
                for (int i = 0; i < Nx; i++)
                {
                    for (int j = 0; j < Ny; j++)
                    {
                        int v0 = vIndex(i, j);
                        int v1 = vIndex(i + 1, j);
                        int v2 = vIndex(i, j + 1);
                        int v3 = vIndex(i + 1, j + 1);
                        int h0 = hMidStart + i * vy + j;           // (i+0.5, j)
                        int h1 = hMidStart + i * vy + (j + 1);     // (i+0.5, j+1)
                        int vA = vMidStart + i * Ny + j;           // (i, j+0.5)
                        int vB = vMidStart + (i + 1) * Ny + j;     // (i+1, j+0.5)
                        int dM = dMidStart + i * Ny + j;           // (i+0.5, j+0.5) diagonal midpoint
                        // Triangle A: (v0, v1, v2) corners -> mids (h0, dM, vA)
                        int[] eA = new int[6] { v0, v1, v2, h0, dM, vA };
                        // Triangle B: (v3, v2, v1) corners -> mids (h1, vA, vB) wait ordering:
                        // corners: v3(1,1)=node2, v2(0,1)=node1? Need to keep CCW around each LT6.
                        // Let B corners CCW be v3=(i+1,j+1), v2=(i,j+1), v1=(i+1,j).
                        // mids: v3-v2 = h1 (top edge mid), v2-v1 = diagonal dM, v1-v3 = vB (right edge mid)
                        int[] eB = new int[6] { v3, v2, v1, h1, dM, vB };
                        m.AddElement(eA); m.AddElement(eB);
                    }
                }
                m.ElementType = ElementType.LT6;
                return m;
            }
            if (type == ElementType.Q4)
            {
                for (int i = 0; i <= Nx; i++) for (int j = 0; j <= Ny; j++) m.AddNode(i * dx, j * dy);
                for (int i = 0; i < Nx; i++) for (int j = 0; j < Ny; j++)
                {
                    int n0 = i * (Ny + 1) + j;
                    int n1 = (i + 1) * (Ny + 1) + j;
                    int n2 = (i + 1) * (Ny + 1) + (j + 1);
                    int n3 = i * (Ny + 1) + (j + 1);
                    m.AddElement(new int[4] { n0, n1, n2, n3 });
                }
                return m;
            }
            // CST3: rectangular grid split into 2 triangles per cell
            for (int i = 0; i <= Nx; i++) for (int j = 0; j <= Ny; j++) m.AddNode(i * dx, j * dy);
            for (int i = 0; i < Nx; i++) for (int j = 0; j < Ny; j++)
            {
                int n0 = i * (Ny + 1) + j;
                int n1 = (i + 1) * (Ny + 1) + j;
                int n2 = (i + 1) * (Ny + 1) + (j + 1);
                int n3 = i * (Ny + 1) + (j + 1);
                m.AddElement(new int[3] { n0, n1, n2 });
                m.AddElement(new int[3] { n0, n2, n3 });
            }
            return m;
        }

        private static int NodeIndex(int i, int j, int Nx, int Ny, ElementType type)
        {
            if (type == ElementType.LT6)
            {
                int vx = Nx + 1, vy = Ny + 1;
                if (i <= Nx && j <= Ny) return i * vy + j;
                // should not be called for midnodes in boundary setups
                return 0;
            }
            return i * (Ny + 1) + j;
        }

        // edge: 0=bottom(y=0),1=right(x=L),2=top(y=H),3=left(x=0). Only returns corner nodes on the edge (for simplicity).
        private static List<int> GetEdgeNodes(int Nx, int Ny, int edge, ElementType type)
        {
            List<int> res = new List<int>();
            if (type == ElementType.LT6)
            {
                int vx = Nx + 1, vy = Ny + 1;
                int vStart = 0;
                int hStart = vx * vy;
                int vStartMid = hStart + Nx * vy;
                switch (edge)
                {
                    case 0:
                        for (int i = 0; i <= Nx; i++) res.Add(vStart + i * vy + 0);
                        for (int i = 0; i < Nx; i++) res.Add(hStart + i * vy + 0);
                        break;
                    case 2:
                        for (int i = 0; i <= Nx; i++) res.Add(vStart + i * vy + Ny);
                        for (int i = 0; i < Nx; i++) res.Add(hStart + i * vy + Ny);
                        break;
                    case 3:
                        for (int j = 0; j <= Ny; j++) res.Add(vStart + 0 + j);
                        for (int j = 0; j < Ny; j++) res.Add(vStartMid + 0 * Ny + j);
                        break;
                    case 1:
                        for (int j = 0; j <= Ny; j++) res.Add(vStart + Nx * vy + j);
                        for (int j = 0; j < Ny; j++) res.Add(vStartMid + Nx * Ny + j);
                        break;
                }
                return res;
            }
            // Q4/CST
            switch (edge)
            {
                case 0: for (int i = 0; i <= Nx; i++) res.Add(i * (Ny + 1) + 0); break;
                case 2: for (int i = 0; i <= Nx; i++) res.Add(i * (Ny + 1) + Ny); break;
                case 3: for (int j = 0; j <= Ny; j++) res.Add(j); break;
                case 1: for (int j = 0; j <= Ny; j++) res.Add(Nx * (Ny + 1) + j); break;
            }
            return res;
        }
    }
}
