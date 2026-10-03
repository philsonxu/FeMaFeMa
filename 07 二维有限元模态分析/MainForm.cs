namespace ModalFEM2D
{
    using System;
    using System.Drawing;
    using System.IO;
    using System.Threading.Tasks;
    using System.Windows.Forms;
    using ModalFEM2D.Elements;
    using ModalFEM2D.IO;
    using ModalFEM2D.Mesh;
    using ModalFEM2D.Rendering;
    using ModalFEM2D.Solvers;

    public enum MeshPreset
    {
        CantileverBeam,
        SimplySupportedBeam,
        LidDrivenCavity
    }

    public class Parameters
    {
        // geometry
        public double Length { get; set; } = 6.0;
        public double Height { get; set; } = 1.0;
        public int Nx { get; set; } = 24;
        public int Ny { get; set; } = 6;
        // material
        public double YoungModulus { get; set; } = 2.1e11;
        public double PoissonRatio { get; set; } = 0.3;
        public double Density { get; set; } = 7800.0;
        public double Thickness { get; set; } = 0.05;
        // structural loads
        public double TipLoad { get; set; } = 2000.0;
        // fluid
        public double Viscosity { get; set; } = 1e-3;
        public double LidVelocity { get; set; } = 1.0;
        // modal
        public int NumModes { get; set; } = 6;
        // SIMPLE
        public double AlphaU { get; set; } = 0.7;
        public double AlphaP { get; set; } = 0.3;
        public double NSTol { get; set; } = 1e-3;
        public int NSMaxOuter { get; set; } = 200;
        public double TimeStep { get; set; } = 0.05;
        public int NumTimeSteps { get; set; } = 20;
        // display
        public double DeformScale { get; set; } = 10.0;
    }

    public class MainForm : Form
    {
        private FEMesh _mesh;
        private Parameters _params = new Parameters();
        private readonly PropertyGrid _pg;
        private readonly PictureBox _pic;
        private readonly TextBox _log;
        private readonly ComboBox _fieldBox;
        private readonly ComboBox _modeBox;
        private readonly ComboBox _krylovBox;
        private readonly ComboBox _precBox;
        private MenuStrip _menu;
        private ToolStripMenuItem _miFile, _miImport, _miExport, _miExit,
            _miMesh, _miCstCant, _miLt6Cant, _miCstSS, _miQ4Cavity,
            _miSolve, _miStructCG, _miStructGMRES, _miStructBiCG, _miStructAMG, _miModal, _miNSSteady, _miNSUnsteady,
            _miHelp, _miAbout;

        public MainForm()
        {
            Text = "ModalFEM2D - 二维有限元计算 (模态分析/静力/NS)";
            Width = 1280; Height = 800; StartPosition = FormStartPosition.CenterScreen;

            _menu = new MenuStrip();
            BuildMenu();
            Controls.Add(_menu);

            SplitContainer split = new SplitContainer { Dock = DockStyle.Fill, SplitterDistance = 300 };

            // Left panel
            Panel left = new Panel { Dock = DockStyle.Fill, Padding = new Padding(6) };
            GroupBox gbP = new GroupBox { Text = "参数", Dock = DockStyle.Top, Height = 430 };
            _pg = new PropertyGrid { Dock = DockStyle.Fill, SelectedObject = _params };
            gbP.Controls.Add(_pg);
            GroupBox gbD = new GroupBox { Text = "显示", Dock = DockStyle.Top, Height = 120, Padding = new Padding(6) };
            _fieldBox = new ComboBox { Dock = DockStyle.Top, DropDownStyle = ComboBoxStyle.DropDownList };
            _fieldBox.Items.AddRange(new object[] { "网格", "位移大小", "变形图", "应力 σxx", "应力 σyy", "应力 τxy", "Von Mises", "速度大小", "压力", "模态振型" });
            _fieldBox.SelectedIndex = 0;
            _fieldBox.SelectedIndexChanged += (s, e) => { RefreshRender(); };
            _modeBox = new ComboBox { Dock = DockStyle.Top, DropDownStyle = ComboBoxStyle.DropDownList, Top = 30 };
            _krylovBox = new ComboBox { Dock = DockStyle.Top, DropDownStyle = ComboBoxStyle.DropDownList, Top = 60 };
            _precBox = new ComboBox { Dock = DockStyle.Top, DropDownStyle = ComboBoxStyle.DropDownList, Top = 90 };
            _krylovBox.Items.AddRange(new object[] { "CG", "GMRES(60)", "BiCGSTAB" }); _krylovBox.SelectedIndex = 0;
            _precBox.Items.AddRange(new object[] { "Jacobi", "AMG" }); _precBox.SelectedIndex = 0;
            gbD.Controls.Add(_precBox); gbD.Controls.Add(_krylovBox); gbD.Controls.Add(_modeBox); gbD.Controls.Add(_fieldBox);
            left.Controls.Add(gbD); left.Controls.Add(gbP);

            // Right panel
            Panel right = new Panel { Dock = DockStyle.Fill };
            _pic = new PictureBox { Dock = DockStyle.Fill, BackColor = Color.Black, SizeMode = PictureBoxSizeMode.Zoom };
            SplitContainer rs = new SplitContainer { Dock = DockStyle.Fill, Orientation = Orientation.Horizontal, SplitterDistance = 560 };
            rs.Panel1.Controls.Add(_pic);
            _log = new TextBox { Dock = DockStyle.Fill, Multiline = true, ScrollBars = ScrollBars.Vertical, BackColor = Color.Black, ForeColor = Color.Lime, Font = new Font("Consolas", 9) };
            rs.Panel2.Controls.Add(_log);
            right.Controls.Add(rs);

            split.Panel1.Controls.Add(left); split.Panel2.Controls.Add(right);
            Controls.Add(split);
            MainMenuStrip = _menu;
            Log("就绪。请从 [网格] 菜单选择算例生成网格，或 [文件]→[导入 INP...]。");
        }

        private void BuildMenu()
        {
            _miFile = new ToolStripMenuItem("文件(&F)");
            _miImport = new ToolStripMenuItem("导入 INP..."); _miImport.Click += (s, e) => DoImport();
            _miExport = new ToolStripMenuItem("导出结果..."); _miExport.Click += (s, e) => DoExport();
            _miExit = new ToolStripMenuItem("退出"); _miExit.Click += (s, e) => Close();
            _miFile.DropDownItems.AddRange(new ToolStripItem[] { _miImport, _miExport, new ToolStripSeparator(), _miExit });

            _miMesh = new ToolStripMenuItem("网格(&M)");
            _miCstCant = new ToolStripMenuItem("悬臂梁 CST3"); _miCstCant.Click += (s, e) => Generate(MeshPreset.CantileverBeam, ElementType.CST3);
            _miLt6Cant = new ToolStripMenuItem("悬臂梁 LT6"); _miLt6Cant.Click += (s, e) => Generate(MeshPreset.CantileverBeam, ElementType.LT6);
            _miCstSS = new ToolStripMenuItem("简支梁 CST3 (模态)"); _miCstSS.Click += (s, e) => Generate(MeshPreset.SimplySupportedBeam, ElementType.CST3);
            _miQ4Cavity = new ToolStripMenuItem("顶盖驱动方腔 Q4"); _miQ4Cavity.Click += (s, e) => Generate(MeshPreset.LidDrivenCavity, ElementType.Q4);
            _miMesh.DropDownItems.AddRange(new ToolStripItem[] { _miCstCant, _miLt6Cant, _miCstSS, new ToolStripSeparator(), _miQ4Cavity });

            _miSolve = new ToolStripMenuItem("求解(&S)");
            _miStructCG = new ToolStripMenuItem("结构静力 (CG+Jacobi)"); _miStructCG.Click += async (s, e) => await RunStruct(KrylovKind.CG, PrecondKind.Jacobi);
            _miStructGMRES = new ToolStripMenuItem("结构静力 (GMRES+Jacobi)"); _miStructGMRES.Click += async (s, e) => await RunStruct(KrylovKind.GMRES, PrecondKind.Jacobi);
            _miStructBiCG = new ToolStripMenuItem("结构静力 (BiCGSTAB+Jacobi)"); _miStructBiCG.Click += async (s, e) => await RunStruct(KrylovKind.BiCGSTAB, PrecondKind.Jacobi);
            _miStructAMG = new ToolStripMenuItem("结构静力 (CG+AMG)  ⭐"); _miStructAMG.Click += async (s, e) => await RunStruct(KrylovKind.CG, PrecondKind.AMG);
            _miModal = new ToolStripMenuItem("模态分析 (子空间迭代 + AMG) ⭐⭐"); _miModal.Click += async (s, e) => await RunModal();
            _miNSSteady = new ToolStripMenuItem("NS 定常 (SIMPLE)"); _miNSSteady.Click += async (s, e) => await RunNS(false);
            _miNSUnsteady = new ToolStripMenuItem("NS 非稳态 (SIMPLE + 向后欧拉)"); _miNSUnsteady.Click += async (s, e) => await RunNS(true);
            _miSolve.DropDownItems.AddRange(new ToolStripItem[] { _miStructCG, _miStructGMRES, _miStructBiCG, _miStructAMG, new ToolStripSeparator(), _miModal, new ToolStripSeparator(), _miNSSteady, _miNSUnsteady });

            _miHelp = new ToolStripMenuItem("帮助(&H)");
            _miAbout = new ToolStripMenuItem("关于"); _miAbout.Click += (s, e) => MessageBox.Show(
                "ModalFEM2D v1.0\n\n" +
                "二维有限元桌面程序\n" +
                "· 结构静力 (CST3/LT6/Q4, 并行 CSR, CG/GMRES/BiCGSTAB+Jacobi/AMG)\n" +
                "· 模态分析 (子空间逆迭代 + Rayleigh-Ritz + CG/AMG)\n" +
                "· Navier-Stokes (SIMPLE + Picard 对流 + 向后欧拉非稳态)\n" +
                "· 多核 SIMD SpMV, VTK/CSV/PNG 输出\n\n" +
                ".NET 10.0 / WinForms / C#", "关于", MessageBoxButtons.OK, MessageBoxIcon.Information);
            _miHelp.DropDownItems.Add(_miAbout);

            _menu.Items.AddRange(new ToolStripItem[] { _miFile, _miMesh, _miSolve, _miHelp });
        }

        private void Log(string msg)
        {
            if (InvokeRequired) { BeginInvoke(new Action<string>(Log), msg); return; }
            _log.AppendText("[" + DateTime.Now.ToString("HH:mm:ss") + "] " + msg + Environment.NewLine);
        }

        private void Generate(MeshPreset preset, ElementType type)
        {
            Parameters p = _params;
            if (preset == MeshPreset.CantileverBeam)
                _mesh = StructuredMeshGenerator.CantileverBeam(p.Length, p.Height, p.Nx, p.Ny, type, p.YoungModulus, p.PoissonRatio, p.Density, p.Thickness, p.TipLoad);
            else if (preset == MeshPreset.SimplySupportedBeam)
                _mesh = StructuredMeshGenerator.SimplySupportedBeam(p.Length, p.Height, p.Nx, p.Ny, type, p.YoungModulus, p.PoissonRatio, p.Density, p.Thickness, 1.0);
            else
                _mesh = StructuredMeshGenerator.LidDrivenCavity(p.Length, p.Nx, type, p.Density, p.Viscosity, p.LidVelocity);
            _mesh.Viscosity = p.Viscosity;
            _fieldBox.SelectedIndex = 0;
            RefreshModes();
            RefreshRender();
            Log(string.Format("网格生成完成: {0} 节点, {1} 单元, 类型={2}", _mesh.NumNodes, _mesh.NumElements, type));
        }

        private void DoImport()
        {
            using (OpenFileDialog ofd = new OpenFileDialog())
            {
                ofd.Filter = "INP 文件 (*.inp)|*.inp|所有文件 (*.*)|*.*";
                string exDir = Path.Combine(AppDomain.CurrentDomain.BaseDirectory, "Examples");
                if (Directory.Exists(exDir)) ofd.InitialDirectory = exDir;
                if (ofd.ShowDialog(this) == DialogResult.OK)
                {
                    try { _mesh = MeshReader.Read(ofd.FileName); _mesh.Viscosity = _params.Viscosity; RefreshModes(); RefreshRender(); Log("已导入: " + ofd.FileName); }
                    catch (Exception ex) { MessageBox.Show("导入失败: " + ex.Message); }
                }
            }
        }

        private void DoExport()
        {
            if (_mesh == null) { MessageBox.Show("请先生成或导入网格并运行求解。"); return; }
            using (FolderBrowserDialog fd = new FolderBrowserDialog())
            {
                fd.Description = "选择输出目录";
                if (fd.ShowDialog(this) == DialogResult.OK)
                {
                    try
                    {
                        string dir = fd.SelectedPath;
                        ResultWriter.WriteCSV(_mesh, dir);
                        ResultWriter.WriteVTK(_mesh, dir);
                        DisplayField df = CurrentField();
                        int mi = _modeBox.SelectedIndex >= 0 ? _modeBox.SelectedIndex : 0;
                        using (Bitmap bmp = ContourRenderer.Render(_mesh, df, mi, _params.DeformScale, 1600, 1200))
                            bmp.Save(Path.Combine(dir, "contour.png"), System.Drawing.Imaging.ImageFormat.Png);
                        Log("结果已导出到 " + dir);
                        MessageBox.Show("已导出到:\n" + dir);
                    }
                    catch (Exception ex) { MessageBox.Show("导出失败: " + ex.Message); }
                }
            }
        }

        private async Task RunStruct(KrylovKind ks, PrecondKind pk)
        {
            if (_mesh == null) { MessageBox.Show("请先建立网格。"); return; }
            _menu.Enabled = false;
            try
            {
                StructuralResult r = null;
                await Task.Run(() => { r = StructuralSolver.Solve(_mesh, ks, pk, Log); });
                Log(string.Format("静力求解完成: {0}, {1} 次迭代", r.SolverUsed, r.Iterations));
                _fieldBox.SelectedIndex = 6; // Von Mises
                RefreshRender();
            }
            catch (Exception ex) { Log("求解失败: " + ex.Message); MessageBox.Show("求解异常: " + ex.Message); }
            finally { _menu.Enabled = true; }
        }

        private async Task RunModal()
        {
            if (_mesh == null) { MessageBox.Show("请先建立网格。"); return; }
            _menu.Enabled = false;
            try
            {
                ModalResult r = null;
                await Task.Run(() => { r = ModalSolver.Solve(_mesh, _params.NumModes, Log); });
                RefreshModes();
                Log(string.Format("模态分析完成: 共提取 {0} 阶, 子空间迭代 {1} 次", r.TargetModes, r.SubspaceIterations));
                for (int i = 0; i < r.FrequenciesHz.Length; i++)
                    Log(string.Format("  第 {0} 阶: f = {1:F3} Hz", i + 1, r.FrequenciesHz[i]));
                _fieldBox.SelectedIndex = 9; // ModeShape
                _modeBox.SelectedIndex = 0;
                RefreshRender();
            }
            catch (Exception ex) { Log("模态求解失败: " + ex.Message); MessageBox.Show("求解异常: " + ex.Message); }
            finally { _menu.Enabled = true; }
        }

        private async Task RunNS(bool unsteady)
        {
            if (_mesh == null) { MessageBox.Show("请先建立网格。"); return; }
            _menu.Enabled = false;
            try
            {
                NSResult r = null;
                Parameters p = _params;
                await Task.Run(() =>
                {
                    if (unsteady)
                        r = NavierStokesSolver.SolveUnsteady(_mesh, p.NSMaxOuter, p.NSTol, p.AlphaU, p.AlphaP, p.Density, p.Viscosity, p.TimeStep, p.NumTimeSteps, Log);
                    else
                        r = NavierStokesSolver.SolveSteady(_mesh, p.NSMaxOuter, p.NSTol, p.AlphaU, p.AlphaP, p.Density, p.Viscosity, Log);
                });
                Log(string.Format("NS 求解完成: outer={0}, |div|={1:E3}", r.OuterIterations, r.FinalResidual));
                _fieldBox.SelectedIndex = 7; // Velocity
                RefreshRender();
            }
            catch (Exception ex) { Log("NS 求解失败: " + ex.Message); MessageBox.Show("求解异常: " + ex.Message); }
            finally { _menu.Enabled = true; }
        }

        private DisplayField CurrentField()
        {
            switch (_fieldBox.SelectedIndex)
            {
                case 0: return DisplayField.Mesh;
                case 1: return DisplayField.Displacement;
                case 2: return DisplayField.Deformed;
                case 3: return DisplayField.StressXX;
                case 4: return DisplayField.StressYY;
                case 5: return DisplayField.StressXY;
                case 6: return DisplayField.VonMises;
                case 7: return DisplayField.Velocity;
                case 8: return DisplayField.Pressure;
                case 9: return DisplayField.ModeShape;
                default: return DisplayField.Mesh;
            }
        }

        private void RefreshModes()
        {
            _modeBox.Items.Clear();
            if (_mesh != null && _mesh.ModalFreqHz != null)
            {
                for (int i = 0; i < _mesh.ModalFreqHz.Length; i++)
                    _modeBox.Items.Add(string.Format("第 {0} 阶  f={1:F2} Hz", i + 1, _mesh.ModalFreqHz[i]));
                if (_modeBox.Items.Count > 0) _modeBox.SelectedIndex = 0;
            }
            _modeBox.SelectedIndexChanged -= ModeChanged;
            _modeBox.SelectedIndexChanged += ModeChanged;
        }

        private void ModeChanged(object s, EventArgs e) { RefreshRender(); }

        private void RefreshRender()
        {
            if (_mesh == null) return;
            try
            {
                int mi = _modeBox.SelectedIndex;
                DisplayField f = CurrentField();
                Bitmap bmp = ContourRenderer.Render(_mesh, f, mi, _params.DeformScale, _pic.Width, _pic.Height);
                Image old = _pic.Image; _pic.Image = bmp; if (old != null) old.Dispose();
            }
            catch (Exception ex) { Log("渲染失败: " + ex.Message); }
        }
    }
}
