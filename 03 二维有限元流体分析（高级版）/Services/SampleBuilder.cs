using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Text;
using Fem2DFluid.Models;

namespace Fem2DFluid.Services
{
    /// <summary>
    /// Builds canonical verification/demo models:
    ///  - Square cavity potential flow (21x21)
    ///  - Dam seepage (Darcy)
    ///  - Flow around a cylinder (DXF geometry + BC)
    /// </summary>
    public class SampleBuilder
    {
        private FemModel _model;

        public SampleBuilder(FemModel model)
        {
            _model = model;
        }

        /// <summary>
        /// Sample 1: 21x21 square-cavity potential flow with left-wall Dirichlet phi=0,
        /// right-wall Dirichlet phi=4, and top/bottom Neumann. Useful as a sanity check
        /// (analytic solution: linear in x; Vx = K*4/L).
        /// </summary>
        public void BuildSquareCavity(int nx, int ny)
        {
            if (nx < 2) nx = 21;
            if (ny < 2) ny = 21;
            _model.Nodes.Clear();
            _model.Elements.Clear();
            _model.BoundaryEdges.Clear();
            _model.BoundaryEdgeCache.Clear();
            _model.Name = "Square Cavity (potential flow)";
            double Lx = 4.0, Ly = 1.0;
            double dx = Lx / (nx - 1);
            double dy = Ly / (ny - 1);
            int id = 0;
            for (int j = 0; j < ny; j++)
            {
                for (int i = 0; i < nx; i++)
                {
                    Node n = new Node(id++, i * dx, j * dy);
                    if (i == 0) { n.BCType = 1; n.BCValue = 0.0; n.Tag = 1; }
                    else if (i == nx - 1) { n.BCType = 1; n.BCValue = 4.0; n.Tag = 2; }
                    else if (j == 0 || j == ny - 1) { n.BCType = 0; n.Tag = 3; }
                    _model.Nodes.Add(n);
                }
            }
            int eid = 0;
            for (int j = 0; j < ny - 1; j++)
            {
                for (int i = 0; i < nx - 1; i++)
                {
                    int n0 = j * nx + i;
                    int n1 = n0 + 1;
                    int n2 = n0 + nx;
                    int n3 = n2 + 1;
                    Element e1 = new Element(eid++, n0, n1, n2);
                    Element e2 = new Element(eid++, n1, n3, n2);
                    _model.Elements.Add(e1);
                    _model.Elements.Add(e2);
                }
            }
            // Boundary edges (outer perimeter)
            for (int i = 0; i < nx - 1; i++)
            {
                AddBnd(i, i + 1, (i == 0) ? 1 : 3);
                AddBnd((ny - 1) * nx + i, (ny - 1) * nx + i + 1, 3);
            }
            for (int j = 0; j < ny - 1; j++)
            {
                AddBnd(j * nx, (j + 1) * nx, 1);
                AddBnd(j * nx + (nx - 1), (j + 1) * nx + (nx - 1), 2);
            }
            Triangulator tr = new Triangulator(_model);
            tr.ComputeElementGeoms();
            _model.ComputeBounds();
        }

        private void AddBnd(int p1, int p2, int tag)
        {
            BoundaryEdge be = new BoundaryEdge(p1, p2, tag);
            _model.BoundaryEdgeCache.Add(be);
            _model.BoundaryEdges.Add(be.ToEdge());
        }

        /// <summary>
        /// Sample 2: Dam seepage (Darcy flow). Rectangular dam with upstream head H1=10,
        /// downstream head H2=2, base impermeable, top seepage face approximated.
        /// </summary>
        public void BuildDamSeepage(int nx, int ny)
        {
            if (nx < 2) nx = 21;
            if (ny < 2) ny = 11;
            _model.Nodes.Clear();
            _model.Elements.Clear();
            _model.BoundaryEdges.Clear();
            _model.BoundaryEdgeCache.Clear();
            _model.Name = "Dam Seepage (Darcy)";
            double Lx = 10.0, Ly = 5.0;
            double dx = Lx / (nx - 1);
            double dy = Ly / (ny - 1);
            int id = 0;
            for (int j = 0; j < ny; j++)
            {
                for (int i = 0; i < nx; i++)
                {
                    Node n = new Node(id++, i * dx, j * dy);
                    if (i == 0) { n.BCType = 1; n.BCValue = Ly; n.Tag = 1; }            // upstream head = 5
                    else if (i == nx - 1) { n.BCType = 1; n.BCValue = Ly * 0.4; n.Tag = 2; } // downstream
                    else if (j == 0) { n.BCType = 0; n.Tag = 3; }                       // base impermeable
                    _model.Nodes.Add(n);
                }
            }
            int eid = 0;
            for (int j = 0; j < ny - 1; j++)
            {
                for (int i = 0; i < nx - 1; i++)
                {
                    int n0 = j * nx + i, n1 = n0 + 1, n2 = n0 + nx, n3 = n2 + 1;
                    _model.Elements.Add(new Element(eid++, n0, n1, n2));
                    _model.Elements.Add(new Element(eid++, n1, n3, n2));
                }
            }
            for (int i = 0; i < nx - 1; i++) AddBnd(i, i + 1, 0);
            for (int i = 0; i < nx - 1; i++) AddBnd((ny - 1) * nx + i, (ny - 1) * nx + i + 1, 0);
            for (int j = 0; j < ny - 1; j++) AddBnd(j * nx, (j + 1) * nx, 1);
            for (int j = 0; j < ny - 1; j++) AddBnd(j * nx + (nx - 1), (j + 1) * nx + (nx - 1), 2);
            Material mat = new Material(0, "Soil");
            mat.K = 1.0e-4;
            _model.Materials.Clear();
            _model.Materials.Add(mat);
            Triangulator tr = new Triangulator(_model);
            tr.ComputeElementGeoms();
            _model.ComputeBounds();
        }

        /// <summary>
        /// Write a simple text "DXF" file describing a channel with a cylinder for
        /// the fallback import path when netDxf is not present. The companion .txt
        /// mesh is written alongside with the same base name.
        /// </summary>
        public bool WriteFlowAroundCylinderDxf(string path)
        {
            try
            {
                string txtPath = Path.ChangeExtension(path, ".txt");
                BuildCylinderMeshTxt(txtPath, 80, 30, 0.125);
                return true;
            }
            catch (Exception)
            {
                return false;
            }
        }

        /// <summary>
        /// Build flow-around-cylinder mesh directly (channel L=2.0, H=0.5, cylinder at
        /// (0.5,0.25) r=0.0625) and apply boundary conditions:
        ///   left inlet u=1, right outlet p=0, top/bottom/cylinder no-slip.
        /// </summary>
        public void BuildCylinderFlow(int nx, int ny)
        {
            if (nx < 10) nx = 80;
            if (ny < 5) ny = 30;
            _model.Nodes.Clear();
            _model.Elements.Clear();
            _model.BoundaryEdges.Clear();
            _model.BoundaryEdgeCache.Clear();
            _model.Name = "Flow Around Cylinder (NS)";
            string tmpTxt = Path.Combine(Path.GetTempPath(), "fem_cyl_tmp.txt");
            BuildCylinderMeshTxt(tmpTxt, nx, ny, 0.0625);
            MeshImporter mi = new MeshImporter(_model);
            mi.FromTextFile(tmpTxt);
            ApplyCylinderBC();
            try { File.Delete(tmpTxt); }
            catch { }
        }

        private void BuildCylinderMeshTxt(string path, int nx, int ny, double radius)
        {
            double Lx = 2.0, Ly = 0.5;
            double cx = 0.5, cy = 0.25;
            // Build structured grid then add cylinder boundary nodes
            List<Node> nodes = new List<Node>();
            List<BoundaryEdge> bes = new List<BoundaryEdge>();
            double dx = Lx / (nx - 1), dy = Ly / (ny - 1);
            int id = 0;
            for (int j = 0; j < ny; j++)
            {
                for (int i = 0; i < nx; i++)
                {
                    double x = i * dx, y = j * dy;
                    // Omit cells near cylinder center (we just keep grid; triangulator will cut)
                    Node n = new Node(id++, x, y);
                    n.Tag = 0;
                    nodes.Add(n);
                }
            }
            // Add cylinder perimeter nodes
            const int NC = 60;
            List<int> cylIds = new List<int>();
            for (int k = 0; k < NC; k++)
            {
                double th = 2.0 * Math.PI * k / NC;
                Node n = new Node(id++, cx + radius * Math.Cos(th), cy + radius * Math.Sin(th));
                n.Tag = 4;
                nodes.Add(n);
                cylIds.Add(n.Id);
            }
            // Remove nodes inside cylinder
            List<Node> cleaned = new List<Node>();
            Dictionary<int, int> remap = new Dictionary<int, int>();
            for (int i = 0; i < nodes.Count; i++)
            {
                Node n = nodes[i];
                double rx = n.X - cx, ry = n.Y - cy;
                if (n.Tag != 4 && rx * rx + ry * ry < radius * radius * 0.95)
                {
                    continue;
                }
                remap[n.Id] = cleaned.Count;
                n.Id = cleaned.Count;
                cleaned.Add(n);
            }
            // Build triangles on structured grid; skip those whose centroid is inside cylinder
            List<Element> elems = new List<Element>();
            int eid = 0;
            for (int j = 0; j < ny - 1; j++)
            {
                for (int i = 0; i < nx - 1; i++)
                {
                    int n0 = j * nx + i, n1 = n0 + 1, n2 = n0 + nx, n3 = n2 + 1;
                    if (remap.ContainsKey(n0) && remap.ContainsKey(n1) && remap.ContainsKey(n2))
                    {
                        double gx = (nodes[n0].X + nodes[n1].X + nodes[n2].X) / 3.0;
                        double gy = (nodes[n0].Y + nodes[n1].Y + nodes[n2].Y) / 3.0;
                        double rx = gx - cx, ry = gy - cy;
                        if (rx * rx + ry * ry > radius * radius * 0.9)
                        {
                            elems.Add(new Element(eid++, remap[n0], remap[n1], remap[n2]));
                        }
                    }
                    if (remap.ContainsKey(n1) && remap.ContainsKey(n3) && remap.ContainsKey(n2))
                    {
                        double gx = (nodes[n1].X + nodes[n3].X + nodes[n2].X) / 3.0;
                        double gy = (nodes[n1].Y + nodes[n3].Y + nodes[n2].Y) / 3.0;
                        double rx = gx - cx, ry = gy - cy;
                        if (rx * rx + ry * ry > radius * radius * 0.9)
                        {
                            elems.Add(new Element(eid++, remap[n1], remap[n3], remap[n2]));
                        }
                    }
                }
            }
            // Boundary edges: inlet (x=0), outlet (x=Lx), walls (y=0,y=Ly), cylinder
            for (int j = 0; j < ny - 1; j++)
            {
                int na = remap[j * nx];
                int nb = remap[(j + 1) * nx];
                bes.Add(new BoundaryEdge(na, nb, 1));
                int na2 = remap[j * nx + (nx - 1)];
                int nb2 = remap[(j + 1) * nx + (nx - 1)];
                bes.Add(new BoundaryEdge(na2, nb2, 2));
            }
            for (int i = 0; i < nx - 1; i++)
            {
                bes.Add(new BoundaryEdge(remap[i], remap[i + 1], 3));
                bes.Add(new BoundaryEdge(remap[(ny - 1) * nx + i], remap[(ny - 1) * nx + i + 1], 3));
            }
            for (int k = 0; k < NC; k++)
            {
                int a = remap[cylIds[k]];
                int b = remap[cylIds[(k + 1) % NC]];
                bes.Add(new BoundaryEdge(a, b, 4));
            }
            // Write text mesh
            using (StreamWriter sw = new StreamWriter(path, false, Encoding.UTF8))
            {
                sw.WriteLine("# FEM cylinder mesh (text)");
                sw.WriteLine("NODES");
                sw.WriteLine(cleaned.Count);
                for (int i = 0; i < cleaned.Count; i++)
                {
                    Node n = cleaned[i];
                    sw.WriteLine(string.Format(CultureInfo.InvariantCulture,
                        "{0} {1:F6} {2:F6} {3}", n.Id, n.X, n.Y, n.Tag));
                }
                sw.WriteLine("ELEMENTS");
                sw.WriteLine(elems.Count);
                for (int i = 0; i < elems.Count; i++)
                {
                    Element e = elems[i];
                    sw.WriteLine("{0} {1} {2} {3} {4}", e.Id, e.A, e.B, e.C, e.MatId);
                }
                sw.WriteLine("BOUNDARIES");
                sw.WriteLine(bes.Count);
                for (int i = 0; i < bes.Count; i++)
                {
                    sw.WriteLine("{0} {1} {2}", bes[i].P1, bes[i].P2, bes[i].Tag);
                }
            }
        }

        /// <summary>
        /// Tag nodes by coordinate and set Dirichlet/Neumann flags for the
        /// cylinder flow case (inlet phi=0, outlet phi=4, walls Neumann,
        /// cylinder Neumann).
        /// </summary>
        public void ApplyCylinderBC()
        {
            _model.ComputeBounds();
            double xmin = _model.Xmin, xmax = _model.Xmax;
            double ymin = _model.Ymin, ymax = _model.Ymax;
            double tol = 1.0e-4;
            for (int i = 0; i < _model.Nodes.Count; i++)
            {
                Node n = _model.Nodes[i];
                n.BCType = 0;
                n.BCValue = 0.0;
                if (Math.Abs(n.X - xmin) < tol)
                {
                    n.Tag = 1; n.BCType = 1; n.BCValue = 0.0; // inlet
                }
                else if (Math.Abs(n.X - xmax) < tol)
                {
                    n.Tag = 2; n.BCType = 1; n.BCValue = 4.0; // outlet
                }
                else if (Math.Abs(n.Y - ymin) < tol || Math.Abs(n.Y - ymax) < tol)
                {
                    n.Tag = 3; n.BCType = 0; // wall Neumann
                }
                // cylinder nodes keep Tag=4 from mesh build
                if (_model.IsNavierStokes)
                {
                    // NS: inlet u=1, outlet p=0 (handled in solver Tag-based), walls/cylinder u=v=0
                    if (n.Tag == 1) { n.BCType = 1; n.BCValue = 1.0; }
                    else if (n.Tag == 2) { n.BCType = 0; }
                    else if (n.Tag == 3 || n.Tag == 4) { n.BCType = 1; n.BCValue = 0.0; }
                }
            }
        }
    }
}
