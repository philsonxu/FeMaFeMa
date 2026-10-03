using System;
using System.Drawing;
using System.IO;
using System.Windows.Forms;
using NonlinearFEM2D.Core;
using NonlinearFEM2D.Elements;
using NonlinearFEM2D.IO;
using NonlinearFEM2D.Mesh;
using NonlinearFEM2D.Rendering;
using NonlinearFEM2D.Solvers;

namespace NonlinearFEM2D
{
    public class AnalParams
    {
        public double Length = 2.0;
        public double Height = 0.5;
        public int Nx = 24;
        public int Ny = 6;
        public double E = 2.1e11;
        public double Nu = 0.3;
        public double Rho = 7850.0;
        public double Thick = 0.01;
        public double EndLoad = -2.0e7;
        public double SigmaY0 = 2.5e8;
        public double Harden = 2.1e9;
        public int LoadSteps = 30;
        public double Viscosity = 1.0e-3;
        public double LidVelocity = 1.0;
        public int CavityN = 20;
        public double Dt = 0.01;
        public int TimeSteps = 50;
        public double AlphaU = 0.7;
        public double AlphaP = 0.3;
        public double DefScale = 20.0;
        public int ModalModes = 4;
    }

    public class MainForm : Form
    {
        private FEMesh _mesh;
        private readonly AnalParams _params = new AnalParams();
        private PropertyGrid _pg;
        private PictureBox _pb;
        private TextBox _log;
        private ComboBox _cbField;
        private ComboBox _cbKrylov;
        private ComboBox _cbPrecond;
        private Label _status;

        public MainForm()
        {
            this.Text = "二维非线性有限元桌面程序 NonlinearFEM2D (.NET 10)";
            this.Width = 1400; this.Height = 900;
            this.StartPosition = FormStartPosition.CenterScreen;
            BuildMenu();
            BuildUI();
        }

        private void BuildMenu()
        {
            MenuStrip ms = new MenuStrip();
            ToolStripMenuItem miFile = new ToolStripMenuItem("文件(&F)");
            ToolStripMenuItem miImport = new ToolStripMenuItem("导入网格 INP...");
            miImport.Click += delegate { DoImport(); };
            ToolStripMenuItem miExport = new ToolStripMenuItem("导出结果 (CSV+VTK+PNG)...");
            miExport.Click += delegate { DoExport(); };
            ToolStripMenuItem miExit = new ToolStripMenuItem("退出");
            miExit.Click += delegate { this.Close(); };
            miFile.DropDownItems.Add(miImport);
            miFile.DropDownItems.Add(miExport);
            miFile.DropDownItems.Add(new ToolStripSeparator());
            miFile.DropDownItems.Add(miExit);

            ToolStripMenuItem miMesh = new ToolStripMenuItem("网格(&M)");
            ToolStripMenuItem miCST = new ToolStripMenuItem("悬臂梁 CST3 (J2 弹塑性)");
            miCST.Click += delegate { BuildCantilever(ElementType.CST3); };
            ToolStripMenuItem miQ4 = new ToolStripMenuItem("悬臂梁 Q4");
            miQ4.Click += delegate { BuildCantilever(ElementType.Q4); };
            ToolStripMenuItem miCav = new ToolStripMenuItem("顶盖驱动方腔 Q4 (NS)");
            miCav.Click += delegate { BuildLidCavity(); };
            miMesh.DropDownItems.Add(miCST);
            miMesh.DropDownItems.Add(miQ4);
            miMesh.DropDownItems.Add(miCav);

            ToolStripMenuItem miSolve = new ToolStripMenuItem("求解(&S)");
            ToolStripMenuItem miLinCG = new ToolStripMenuItem("线性静力 CG+Jacobi");
            miLinCG.Click += delegate { RunLinear(KrylovKind.CG, PrecondKind.Jacobi); };
            ToolStripMenuItem miLinGMRES = new ToolStripMenuItem("线性静力 GMRES+Jacobi");
            miLinGMRES.Click += delegate { RunLinear(KrylovKind.GMRES, PrecondKind.Jacobi); };
            ToolStripMenuItem miLinBCG = new ToolStripMenuItem("线性静力 BiCGSTAB+Jacobi");
            miLinBCG.Click += delegate { RunLinear(KrylovKind.BiCGSTAB, PrecondKind.Jacobi); };
            ToolStripMenuItem miLinAMG = new ToolStripMenuItem("线性静力 CG+AMG ★");
            miLinAMG.Click += delegate { RunLinear(KrylovKind.CG, PrecondKind.AMG); };
            ToolStripMenuItem miNLF = new ToolStripMenuItem("★ 强非线性 J2 塑性（力加载-Newton+AMG）");
            miNLF.Click += delegate { RunNonlinear(LoadControlType.Force); };
            ToolStripMenuItem miNLD = new ToolStripMenuItem("★ 强非线性 J2 塑性（位移加载）");
            miNLD.Click += delegate { RunNonlinear(LoadControlType.Displacement); };
            ToolStripMenuItem miNLGeo = new ToolStripMenuItem("★ 几何+材料双重非线性");
            miNLGeo.Click += delegate { RunNonlinear(LoadControlType.Force, true); };
            ToolStripMenuItem miNS = new ToolStripMenuItem("NS 定常 SIMPLE");
            miNS.Click += delegate { RunNS(true); };
            ToolStripMenuItem miNSU = new ToolStripMenuItem("NS 非稳态");
            miNSU.Click += delegate { RunNS(false); };
            miSolve.DropDownItems.Add(miLinCG); miSolve.DropDownItems.Add(miLinGMRES); miSolve.DropDownItems.Add(miLinBCG);
            miSolve.DropDownItems.Add(miLinAMG);
            miSolve.DropDownItems.Add(new ToolStripSeparator());
            miSolve.DropDownItems.Add(miNLF); miSolve.DropDownItems.Add(miNLD); miSolve.DropDownItems.Add(miNLGeo);
            miSolve.DropDownItems.Add(new ToolStripSeparator());
            miSolve.DropDownItems.Add(miNS); miSolve.DropDownItems.Add(miNSU);

            ToolStripMenuItem miHelp = new ToolStripMenuItem("帮助(&H)");
            ToolStripMenuItem miAbout = new ToolStripMenuItem("关于");
            miAbout.Click += delegate
            {
                MessageBox.Show(
                    "二维强非线性有限元桌面程序\n\n" +
                    "• 单元：CST3（J2 弹塑性, 径向返回）、LT6、Q4\n" +
                    "• 存储：CSR 稀疏 + 线程池并行装配 + SIMD SpMV (Vector<double>)\n" +
                    "• Krylov：CG / GMRES(m) / BiCGSTAB\n" +
                    "• 预条件：Jacobi / AMG (Smoothed Aggregation + V-cycle)\n" +
                    "• 非线性：Newton-Raphson 增量 + 力/位移控制 + 几何非线性可选\n" +
                    "• NS：Picard 对流 + 向后欧拉非稳态 + SIMPLE 压力修正\n" +
                    "• 输出：CSV / VTK (ParaView) / Jet 云图 PNG\n\n" +
                    ".NET 10 · WinForms · C# (无 var)",
                    "关于", MessageBoxButtons.OK, MessageBoxIcon.Information);
            };
            miHelp.DropDownItems.Add(miAbout);

            ms.Items.Add(miFile); ms.Items.Add(miMesh); ms.Items.Add(miSolve); ms.Items.Add(miHelp);
            this.MainMenuStrip = ms;
            this.Controls.Add(ms);
        }

        private void BuildUI()
        {
            SplitContainer split = new SplitContainer();
            split.Dock = DockStyle.Fill;
            split.SplitterDistance = 280;
            split.FixedPanel = FixedPanel.Panel1;

            Panel left = new Panel(); left.Dock = DockStyle.Fill;
            _pg = new PropertyGrid();
            _pg.Dock = DockStyle.Fill;
            _pg.SelectedObject = this._params;
            _pg.PropertySort = PropertySort.Categorized;
            left.Controls.Add(_pg);

            FlowLayoutPanel top = new FlowLayoutPanel();
            top.Dock = DockStyle.Top;
            top.Height = 40;
            Label l1 = new Label() { Text = "显示字段：", AutoSize = true, Margin = new Padding(6, 8, 0, 0) };
            _cbField = new ComboBox() { Width = 160, DropDownStyle = ComboBoxStyle.DropDownList };
            _cbField.Items.AddRange(new object[] { "网格 Mesh", "位移大小 Displacement", "变形 Deform", "σxx", "σyy", "τxy", "Von Mises", "等效塑性应变 εp", "速度 |v|", "压力 Pressure" });
            _cbField.SelectedIndex = 6;
            _cbField.SelectedIndexChanged += delegate { RefreshView(); };
            Label l2 = new Label() { Text = "Krylov：", AutoSize = true, Margin = new Padding(6, 8, 0, 0) };
            _cbKrylov = new ComboBox() { Width = 100, DropDownStyle = ComboBoxStyle.DropDownList };
            _cbKrylov.Items.AddRange(new object[] { "CG", "GMRES", "BiCGSTAB" });
            _cbKrylov.SelectedIndex = 0;
            Label l3 = new Label() { Text = "预条件：", AutoSize = true, Margin = new Padding(6, 8, 0, 0) };
            _cbPrecond = new ComboBox() { Width = 80, DropDownStyle = ComboBoxStyle.DropDownList };
            _cbPrecond.Items.AddRange(new object[] { "None", "Jacobi", "AMG" });
            _cbPrecond.SelectedIndex = 1;
            top.Controls.Add(l1); top.Controls.Add(_cbField); top.Controls.Add(l2); top.Controls.Add(_cbKrylov); top.Controls.Add(l3); top.Controls.Add(_cbPrecond);

            _pb = new PictureBox();
            _pb.Dock = DockStyle.Fill;
            _pb.BackColor = Color.White;
            _pb.SizeMode = PictureBoxSizeMode.Zoom;

            _log = new TextBox();
            _log.Dock = DockStyle.Bottom;
            _log.Height = 220;
            _log.Multiline = true; _log.ScrollBars = ScrollBars.Vertical;
            _log.Font = new Font("Consolas", 9); _log.BackColor = Color.Black; _log.ForeColor = Color.Lime;

            Panel right = new Panel(); right.Dock = DockStyle.Fill;
            right.Controls.Add(_pb);
            right.Controls.Add(_log);
            right.Controls.Add(top);

            _status = new Label(); _status.Dock = DockStyle.Bottom; _status.Height = 22; _status.Text = "就绪";
            this.Controls.Add(_status);

            split.Panel1.Controls.Add(left);
            split.Panel2.Controls.Add(right);
            this.Controls.Add(split);
        }

        private void Log(string msg)
        {
            if (this.InvokeRequired) { this.BeginInvoke(new Action<string>(Log), msg); return; }
            this._log.AppendText(msg + Environment.NewLine);
        }

        private void BuildCantilever(ElementType t)
        {
            this._mesh = StructuredMeshGenerator.GenerateCantilever(
                t, this._params.Length, this._params.Height, this._params.Nx, this._params.Ny,
                this._params.E, this._params.Nu, this._params.Rho, this._params.Thick,
                this._params.EndLoad, this._params.SigmaY0, this._params.Harden);
            RefreshView();
            Log("[网格] 悬臂梁已生成：" + this._mesh.Elements.Count + " 个单元, " + this._mesh.NumNodes + " 节点.");
        }
        private void BuildLidCavity()
        {
            this._mesh = StructuredMeshGenerator.GenerateLidDrivenCavity(this._params.CavityN, this._params.Length,
                this._params.LidVelocity, this._params.Rho, this._params.Viscosity);
            RefreshView();
            Log("[网格] 顶盖方腔 Q4 已生成：" + this._mesh.Elements.Count + " 单元, " + this._mesh.NumNodes + " 节点.");
        }
        private void DoImport()
        {
            OpenFileDialog fd = new OpenFileDialog();
            fd.Filter = "Abaqus INP (*.inp)|*.inp|所有文件 (*.*)|*.*";
            if (fd.ShowDialog() == DialogResult.OK)
            {
                try
                {
                    this._mesh = MeshReader.ReadInp(fd.FileName);
                    RefreshView();
                    Log("[导入] " + fd.FileName + " | 节点=" + this._mesh.NumNodes + " 单元=" + this._mesh.Elements.Count);
                }
                catch (Exception ex) { Log("导入错误: " + ex.Message); }
            }
        }
        private void DoExport()
        {
            if (this._mesh == null) { MessageBox.Show("请先生成网格并求解."); return; }
            using (FolderBrowserDialog fd = new FolderBrowserDialog())
            {
                if (fd.ShowDialog() == DialogResult.OK)
                {
                    try
                    {
                        ResultWriter.WriteCSV(this._mesh, fd.SelectedPath);
                        ResultWriter.WriteVTK(this._mesh, fd.SelectedPath);
                        string fld = GetSelectedField();
                        ResultWriter.WritePNG(this._mesh, fld, fd.SelectedPath, 0);
                        MessageBox.Show("已导出到: " + fd.SelectedPath);
                        Log("[导出] CSV+VTK+PNG -> " + fd.SelectedPath);
                    }
                    catch (Exception ex) { Log("导出错误: " + ex.Message); }
                }
            }
        }

        private string GetSelectedField()
        {
            return (string)this._cbField.SelectedItem;
        }
        private KrylovKind GetKrylov()
        {
            string s = (string)this._cbKrylov.SelectedItem;
            if (s == "GMRES") return KrylovKind.GMRES;
            if (s == "BiCGSTAB") return KrylovKind.BiCGSTAB;
            return KrylovKind.CG;
        }
        private PrecondKind GetPrecond()
        {
            string s = (string)this._cbPrecond.SelectedItem;
            if (s == "AMG") return PrecondKind.AMG;
            if (s == "None") return PrecondKind.None;
            return PrecondKind.Jacobi;
        }

        private void RunLinear(KrylovKind ks, PrecondKind pk)
        {
            if (this._mesh == null) { MessageBox.Show("请先生成网格"); return; }
            this._status.Text = "正在求解(线性)...";
            System.Threading.ThreadPool.QueueUserWorkItem(delegate
            {
                try
                {
                    this._mesh.ResetFields();
                    StructuralSolver s = new StructuralSolver();
                    s.Krylov = ks; s.Precond = pk;
                    s.Tolerance = 1.0e-8; s.MaxIterations = 5000;
                    s.Log = Log;
                    s.OnIter = delegate (int it, double r) { if (it % 50 == 0) Log("    KSP it=" + it + " res=" + r.ToString("E3")); };
                    DateTime t0 = DateTime.Now;
                    s.Solve(this._mesh);
                    Log("[线性求解] 耗时 " + (DateTime.Now - t0).TotalSeconds.ToString("F2") + "s");
                    this.BeginInvoke(new Action(RefreshView));
                }
                catch (Exception ex) { Log("求解错误: " + ex.Message + "\n" + ex.StackTrace); }
                finally { this.BeginInvoke(new Action(delegate { this._status.Text = "就绪"; })); }
            });
        }

        private void RunNonlinear(LoadControlType lc, bool geom = false)
        {
            if (this._mesh == null) { MessageBox.Show("请先生成网格"); return; }
            this._status.Text = "正在求解(非线性)...";
            System.Threading.ThreadPool.QueueUserWorkItem(delegate
            {
                try
                {
                    this._mesh.ResetFields();
                    NonlinearStaticSolver s = new NonlinearStaticSolver();
                    s.Krylov = GetKrylov(); s.Precond = GetPrecond();
                    s.Tolerance = 1.0e-5; s.MaxIterations = 50;
                    s.LoadSteps = this._params.LoadSteps;
                    s.LoadControl = lc;
                    s.UseGeometricNonlinearity = geom;
                    s.DisplacementTarget = this._params.Height * 0.1;
                    int tipId = -1; double tipX = -1e99;
                    for (int i = 0; i < this._mesh.NumNodes; i++) if (this._mesh.X[i] > tipX) { tipX = this._mesh.X[i]; tipId = i; }
                    s.DisplacementNodeId = tipId;
                    s.DisplacementDof = 1;
                    s.Log = Log;
                    s.OnIter = null;
                    DateTime t0 = DateTime.Now;
                    s.Solve(this._mesh);
                    Log("[强非线性] 耗时 " + (DateTime.Now - t0).TotalSeconds.ToString("F2") + "s");
                    this.BeginInvoke(new Action(RefreshView));
                }
                catch (Exception ex) { Log("求解错误: " + ex.Message + "\n" + ex.StackTrace); }
                finally { this.BeginInvoke(new Action(delegate { this._status.Text = "就绪"; })); }
            });
        }

        private void RunNS(bool steady)
        {
            if (this._mesh == null) { MessageBox.Show("请先生成网格"); return; }
            this._status.Text = "正在求解(NS SIMPLE)...";
            System.Threading.ThreadPool.QueueUserWorkItem(delegate
            {
                try
                {
                    this._mesh.ResetFields();
                    NavierStokesSolver s = new NavierStokesSolver();
                    s.Density = this._params.Rho;
                    s.Viscosity = this._params.Viscosity;
                    s.Dt = this._params.Dt;
                    s.TimeSteps = this._params.TimeSteps;
                    s.AlphaU = this._params.AlphaU;
                    s.AlphaP = this._params.AlphaP;
                    s.Steady = steady;
                    s.Tolerance = 1.0e-3;
                    s.Log = Log;
                    DateTime t0 = DateTime.Now;
                    s.Solve(this._mesh);
                    Log("[NS] 耗时 " + (DateTime.Now - t0).TotalSeconds.ToString("F2") + "s");
                    this.BeginInvoke(new Action(RefreshView));
                }
                catch (Exception ex) { Log("求解错误: " + ex.Message + "\n" + ex.StackTrace); }
                finally { this.BeginInvoke(new Action(delegate { this._status.Text = "就绪"; })); }
            });
        }

        private void RefreshView()
        {
            if (this._mesh == null) return;
            string fld = GetSelectedField();
            ContourRenderer cr = new ContourRenderer();
            Bitmap bmp = cr.Render(this._mesh, fld, Math.Max(800, this._pb.Width), Math.Max(600, this._pb.Height), 0);
            Image old = this._pb.Image;
            this._pb.Image = bmp;
            if (old != null) old.Dispose();
        }
    }
}
