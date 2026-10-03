using System;
using System.Drawing;
using System.IO;
using System.Threading;
using System.Threading.Tasks;
using System.Windows.Forms;
using MultiPhysicsFEM2D.Core;
using MultiPhysicsFEM2D.Elements;
using MultiPhysicsFEM2D.IO;
using MultiPhysicsFEM2D.Mesh;
using MultiPhysicsFEM2D.Rendering;
using MultiPhysicsFEM2D.Solvers;

namespace MultiPhysicsFEM2D
{
    public class SimulationParams
    {
        // 几何
        public double Length = 2.0;
        public double Height = 0.5;
        public int Nx = 20;
        public int Ny = 5;
        // 结构材料
        public double E = 2.1e11;
        public double Nu = 0.3;
        public double Density = 7800.0;
        public double Thickness = 0.01;
        public double EndLoad = -2e6;
        // 热
        public double K = 50.0;
        public double Cp = 500.0;
        public double Alpha = 1.2e-5;
        public double T0 = 293.0;
        public double Thot = 393.0;
        public double Tcold = 293.0;
        // 流体
        public double RhoFluid = 1.0;
        public double MuFluid = 0.01;
        public double Beta = 0.00034;
        public double UTop = 1.0;
        public double TLeft = 1.0;
        public double TRight = 0.0;
        public double Gravity = 10.0;
        public double AlphaU = 0.7;
        public double AlphaP = 0.3;
        public int NSOuter = 300;
        public double NSTol = 1e-4;
        public double DT = 0.01;
        public int TimeSteps = 20;
        // 模态/显示
        public int NumModes = 4;
        public double DeformScale = 1.0;
        public KrylovKind Krylov = KrylovKind.CG;
        public PrecondKind Precond = PrecondKind.Jacobi;
    }

    public class MainForm : Form
    {
        private FEMesh _mesh;
        private SimulationParams _params = new SimulationParams();
        private PropertyGrid _pg;
        private PictureBox _pb;
        private TextBox _logBox;
        private ComboBox _fieldCombo;
        private ComboBox _krylovCombo;
        private ComboBox _pcCombo;
        private Label _status;

        public MainForm()
        {
            Text = "二维多物理场耦合有限元 (Multi-Physics FEM 2D)";
            Width = 1400; Height = 880;
            StartPosition = FormStartPosition.CenterScreen;
            BuildUi();
            BuildMenu();
        }

        private void BuildUi()
        {
            SplitContainer sc = new SplitContainer();
            sc.Dock = DockStyle.Fill;
            sc.SplitterDistance = 300;
            sc.Panel1.AutoScroll = true;
            Controls.Add(sc);

            // 左侧参数
            _pg = new PropertyGrid();
            _pg.Dock = DockStyle.Fill;
            _pg.SelectedObject = _params;
            _pg.PropertySort = PropertySort.Categorized;
            sc.Panel1.Controls.Add(_pg);

            // 左侧底部下拉
            Panel bot = new Panel();
            bot.Dock = DockStyle.Bottom;
            bot.Height = 110;
            sc.Panel1.Controls.Add(bot);
            Label l1 = new Label() { Text = "显示字段:", Top = 8, Left = 8, AutoSize = true };
            _fieldCombo = new ComboBox() { Top = 6, Left = 80, Width = 180, DropDownStyle = ComboBoxStyle.DropDownList };
            _fieldCombo.Items.AddRange(new object[] {
                "网格 Mesh", "位移大小 (变形)", "变形图",
                "应力 σxx", "应力 σyy", "剪应力 τxy", "等效应力 Von Mises",
                "速度大小 |V|", "速度 Vx", "速度 Vy", "压力 P", "温度 T"
            });
            _fieldCombo.SelectedIndex = 0;
            _fieldCombo.SelectedIndexChanged += delegate { RefreshPlot(); };
            Label l2 = new Label() { Text = "Krylov:", Top = 38, Left = 8, AutoSize = true };
            _krylovCombo = new ComboBox() { Top = 36, Left = 80, Width = 180, DropDownStyle = ComboBoxStyle.DropDownList };
            _krylovCombo.Items.AddRange(Enum.GetNames(typeof(KrylovKind)));
            _krylovCombo.SelectedItem = _params.Krylov.ToString();
            _krylovCombo.SelectedIndexChanged += delegate
            {
                _params.Krylov = (KrylovKind)Enum.Parse(typeof(KrylovKind), (string)_krylovCombo.SelectedItem);
            };
            Label l3 = new Label() { Text = "预条件:", Top = 68, Left = 8, AutoSize = true };
            _pcCombo = new ComboBox() { Top = 66, Left = 80, Width = 180, DropDownStyle = ComboBoxStyle.DropDownList };
            _pcCombo.Items.AddRange(Enum.GetNames(typeof(PrecondKind)));
            _pcCombo.SelectedItem = _params.Precond.ToString();
            _pcCombo.SelectedIndexChanged += delegate
            {
                _params.Precond = (PrecondKind)Enum.Parse(typeof(PrecondKind), (string)_pcCombo.SelectedItem);
            };
            bot.Controls.Add(l1); bot.Controls.Add(_fieldCombo);
            bot.Controls.Add(l2); bot.Controls.Add(_krylovCombo);
            bot.Controls.Add(l3); bot.Controls.Add(_pcCombo);

            // 右侧
            SplitContainer rsc = new SplitContainer();
            rsc.Dock = DockStyle.Fill;
            rsc.Orientation = Orientation.Horizontal;
            rsc.SplitterDistance = 560;
            sc.Panel2.Controls.Add(rsc);
            _pb = new PictureBox();
            _pb.Dock = DockStyle.Fill;
            _pb.BackColor = Color.White;
            _pb.SizeMode = PictureBoxSizeMode.Zoom;
            rsc.Panel1.Controls.Add(_pb);
            _logBox = new TextBox();
            _logBox.Dock = DockStyle.Fill;
            _logBox.Multiline = true;
            _logBox.ScrollBars = ScrollBars.Vertical;
            _logBox.Font = new Font("Consolas", 9f);
            _logBox.BackColor = Color.Black;
            _logBox.ForeColor = Color.Lime;
            rsc.Panel2.Controls.Add(_logBox);

            _status = new Label() { Dock = DockStyle.Bottom, Height = 22, Text = "就绪" };
            Controls.Add(_status);
        }

        private void BuildMenu()
        {
            MenuStrip ms = new MenuStrip();
            ToolStripMenuItem file = new ToolStripMenuItem("文件(&F)");
            ToolStripMenuItem fiImport = new ToolStripMenuItem("导入网格 INP...");
            fiImport.Click += delegate { ImportMesh(); };
            ToolStripMenuItem fiExport = new ToolStripMenuItem("导出结果...");
            fiExport.Click += delegate { ExportResults(); };
            ToolStripMenuItem fiExit = new ToolStripMenuItem("退出");
            fiExit.Click += delegate { Close(); };
            file.DropDownItems.Add(fiImport);
            file.DropDownItems.Add(fiExport);
            file.DropDownItems.Add(new ToolStripSeparator());
            file.DropDownItems.Add(fiExit);

            ToolStripMenuItem mesh = new ToolStripMenuItem("网格(&M)");
            ToolStripMenuItem miCantCST = new ToolStripMenuItem("悬臂梁 CST3");
            miCantCST.Click += delegate { GenMesh(StructuredMeshGenerator.Preset.CantileverCST3, ElementType.CST3); };
            ToolStripMenuItem miCantQ4 = new ToolStripMenuItem("悬臂梁 Q4");
            miCantQ4.Click += delegate { GenMesh(StructuredMeshGenerator.Preset.CantileverQ4, ElementType.Q4); };
            ToolStripMenuItem miTherm = new ToolStripMenuItem("★ 热-结构耦合梁 CST3");
            miTherm.Click += delegate { GenMesh(StructuredMeshGenerator.Preset.ThermalBeamCST3, ElementType.CST3); };
            ToolStripMenuItem miLid = new ToolStripMenuItem("顶盖驱动方腔 Q4");
            miLid.Click += delegate { GenMesh(StructuredMeshGenerator.Preset.LidDrivenCavityQ4, ElementType.Q4); };
            ToolStripMenuItem miNat = new ToolStripMenuItem("★ 自然对流方腔 Q4 (Boussinesq)");
            miNat.Click += delegate { GenMesh(StructuredMeshGenerator.Preset.NaturalConvectionQ4, ElementType.Q4); };
            mesh.DropDownItems.Add(miCantCST);
            mesh.DropDownItems.Add(miCantQ4);
            mesh.DropDownItems.Add(miTherm);
            mesh.DropDownItems.Add(miLid);
            mesh.DropDownItems.Add(miNat);

            ToolStripMenuItem solve = new ToolStripMenuItem("求解(&S)");
            ToolStripMenuItem sStructCG = new ToolStripMenuItem("结构静力");
            sStructCG.Click += delegate { RunStructural(); };
            ToolStripMenuItem sThermal = new ToolStripMenuItem("热传导(稳态)");
            sThermal.Click += delegate { RunThermal(true); };
            ToolStripMenuItem sThermalT = new ToolStripMenuItem("热传导(瞬态)");
            sThermalT.Click += delegate { RunThermal(false); };
            ToolStripMenuItem sTS = new ToolStripMenuItem("★ 热-结构耦合(单向)");
            sTS.Click += delegate { RunThermalStructural(); };
            ToolStripMenuItem sNS = new ToolStripMenuItem("NS 定常 SIMPLE");
            sNS.Click += delegate { RunNS(false); };
            ToolStripMenuItem sNST = new ToolStripMenuItem("NS 非稳态");
            sNST.Click += delegate { RunNS(true); };
            ToolStripMenuItem sBq = new ToolStripMenuItem("★ 流-热 Boussinesq 自然对流");
            sBq.Click += delegate { RunBoussinesq(); };
            solve.DropDownItems.Add(sStructCG);
            solve.DropDownItems.Add(sThermal);
            solve.DropDownItems.Add(sThermalT);
            solve.DropDownItems.Add(sTS);
            solve.DropDownItems.Add(sNS);
            solve.DropDownItems.Add(sNST);
            solve.DropDownItems.Add(sBq);

            ToolStripMenuItem help = new ToolStripMenuItem("帮助(&H)");
            ToolStripMenuItem hAbout = new ToolStripMenuItem("关于");
            hAbout.Click += delegate
            {
                MessageBox.Show(
                    "二维多物理场耦合有限元\n" +
                    "Multi-Physics FEM 2D\n\n" +
                    "物理场: 结构静力 / 稳态+瞬态热 / NS SIMPLE / Boussinesq 流-热 / 热-结构单向耦合\n" +
                    "单元: CST3 / LT6 / Q4\n" +
                    "求解器: CG / GMRES / BiCGSTAB + Jacobi/AMG 预条件\n" +
                    "存储: CSR + COO 多线程装配 + SIMD SpMV\n" +
                    "规范: .NET 10.0 WinForms, 全项目无 var 关键字",
                    "关于", MessageBoxButtons.OK, MessageBoxIcon.Information);
            };
            help.DropDownItems.Add(hAbout);

            ms.Items.Add(file); ms.Items.Add(mesh); ms.Items.Add(solve); ms.Items.Add(help);
            MainMenuStrip = ms;
            Controls.Add(ms);
        }

        private void Log(string msg)
        {
            if (InvokeRequired) { BeginInvoke(new Action<string>(Log), msg); return; }
            _logBox.AppendText(msg + Environment.NewLine);
        }

        private void ApplyParamsToMesh()
        {
            if (_mesh == null) return;
            _mesh.YoungModulus = _params.E;
            _mesh.PoissonRatio = _params.Nu;
            _mesh.Density = _params.RhoFluid; // 流体密度用于 NS
            _mesh.Thickness = _params.Thickness;
            _mesh.Conductivity = _params.K;
            _mesh.HeatCapacity = _params.Cp;
            _mesh.ThermalExpansion = _params.Alpha;
            _mesh.ReferenceTemperature = _params.T0;
            _mesh.Viscosity = _params.MuFluid;
            _mesh.ApplyDefaultMaterial();
            // 结构密度单独覆盖（结构默认 7800）
            foreach (FiniteElement e in _mesh.Elements) e.Density = _params.RhoFluid;
            _mesh.UpdateElementCoords();
        }

        private void GenMesh(StructuredMeshGenerator.Preset preset, ElementType et)
        {
            _mesh = StructuredMeshGenerator.Generate(
                _params.Length, _params.Height, _params.Nx, _params.Ny, et, preset,
                0, _params.EndLoad, _params.Thot, _params.Tcold,
                _params.UTop, _params.TLeft, _params.TRight);
            ApplyParamsToMesh();
            Log(string.Format("网格已生成: {0} 节点, {1} 单元 ({2})",
                _mesh.NumNodes, _mesh.Elements.Count, et));
            RefreshPlot();
        }

        private void ImportMesh()
        {
            using (OpenFileDialog dlg = new OpenFileDialog())
            {
                dlg.Filter = "INP mesh|*.inp|All|*.*";
                if (dlg.ShowDialog() == DialogResult.OK)
                {
                    try
                    {
                        _mesh = MeshReader.ReadINP(dlg.FileName);
                        ApplyParamsToMesh();
                        Log("导入成功: " + dlg.FileName + "  nodes=" + _mesh.NumNodes + " elems=" + _mesh.Elements.Count);
                        RefreshPlot();
                    }
                    catch (Exception ex) { Log("导入失败: " + ex.Message); }
                }
            }
        }

        private void EnsureMesh()
        {
            if (_mesh == null) GenMesh(StructuredMeshGenerator.Preset.CantileverCST3, ElementType.CST3);
        }

        private void RunStructural()
        {
            EnsureMesh();
            ApplyParamsToMesh();
            _status.Text = "结构分析中...";
            ThreadPool.QueueUserWorkItem(delegate
            {
                try
                {
                    StructuralSolver ss = new StructuralSolver();
                    ss.Krylov = _params.Krylov; ss.Precond = _params.Precond;
                    ss.Log = Log;
                    StructuralResult r = ss.Solve(_mesh);
                    BeginInvoke(new Action(delegate
                    {
                        Log(string.Format("结构求解完成: iter={0}, res={1:E4}, solver={2}", r.Iterations, r.Residual, r.SolverUsed));
                        _fieldCombo.SelectedItem = "等效应力 Von Mises";
                        RefreshPlot(); _status.Text = "完成";
                    }));
                }
                catch (Exception ex) { BeginInvoke(new Action(() => { Log("求解失败: " + ex.Message); _status.Text = "错误"; })); }
            });
        }

        private void RunThermal(bool steady)
        {
            EnsureMesh(); ApplyParamsToMesh();
            _status.Text = "热分析中...";
            ThreadPool.QueueUserWorkItem(delegate
            {
                try
                {
                    ThermalSolver ts = new ThermalSolver();
                    ts.Steady = steady; ts.TimeStep = _params.DT; ts.TimeSteps = _params.TimeSteps;
                    ts.Krylov = _params.Krylov; ts.Precond = _params.Precond; ts.Log = Log;
                    ThermalResult r;
                    if (steady) r = ts.Solve(_mesh); else r = ts.SolveTransient(_mesh);
                    BeginInvoke(new Action(delegate
                    {
                        Log(string.Format("热传导完成: iter={0}, res={1:E4}", r.Iterations, r.Residual));
                        _fieldCombo.SelectedItem = "温度 T";
                        RefreshPlot(); _status.Text = "完成";
                    }));
                }
                catch (Exception ex) { BeginInvoke(new Action(() => { Log("热分析失败: " + ex.Message); _status.Text = "错误"; })); }
            });
        }

        private void RunThermalStructural()
        {
            EnsureMesh(); ApplyParamsToMesh();
            _status.Text = "热-结构耦合分析中...";
            ThreadPool.QueueUserWorkItem(delegate
            {
                try
                {
                    ThermalStructuralSolver tss = new ThermalStructuralSolver();
                    tss.Krylov = _params.Krylov; tss.Precond = _params.Precond;
                    tss.Log = Log;
                    ThermalStructuralResult r = tss.Solve(_mesh);
                    BeginInvoke(new Action(delegate
                    {
                        Log(string.Format("热-结构耦合完成: thermal iter={0}, struct iter={1}",
                            r.Thermal.Iterations, r.Structural.Iterations));
                        _fieldCombo.SelectedItem = "等效应力 Von Mises";
                        RefreshPlot(); _status.Text = "完成";
                    }));
                }
                catch (Exception ex) { BeginInvoke(new Action(() => { Log("耦合求解失败: " + ex.Message); _status.Text = "错误"; })); }
            });
        }

        private void RunNS(bool transient)
        {
            EnsureMesh(); ApplyParamsToMesh();
            _status.Text = "NS 求解中...";
            ThreadPool.QueueUserWorkItem(delegate
            {
                try
                {
                    NavierStokesSolver ns = new NavierStokesSolver();
                    ns.Transient = transient;
                    ns.TimeStep = _params.DT; ns.TimeSteps = _params.TimeSteps;
                    ns.OuterIterations = _params.NSOuter;
                    ns.AlphaU = _params.AlphaU; ns.AlphaP = _params.AlphaP;
                    ns.Tol = _params.NSTol;
                    ns.Gravity = _params.Gravity;
                    ns.BetaExpansion = _params.Beta;
                    ns.IncludeThermal = false;
                    ns.Log = Log;
                    NSResult r = ns.Solve(_mesh);
                    BeginInvoke(new Action(delegate
                    {
                        Log(string.Format("NS 完成: outer={0}, resU={1:E3}", r.OuterIterations, r.FinalResU));
                        _fieldCombo.SelectedItem = "速度大小 |V|";
                        RefreshPlot(); _status.Text = "完成";
                    }));
                }
                catch (Exception ex) { BeginInvoke(new Action(() => { Log("NS 失败: " + ex.Message); _status.Text = "错误"; })); }
            });
        }

        private void RunBoussinesq()
        {
            EnsureMesh(); ApplyParamsToMesh();
            _status.Text = "Boussinesq 流-热耦合中...";
            ThreadPool.QueueUserWorkItem(delegate
            {
                try
                {
                    NavierStokesSolver ns = new NavierStokesSolver();
                    ns.Transient = false;
                    ns.OuterIterations = _params.NSOuter;
                    ns.AlphaU = _params.AlphaU; ns.AlphaP = _params.AlphaP;
                    ns.Tol = _params.NSTol;
                    ns.Gravity = _params.Gravity;
                    ns.BetaExpansion = _params.Beta;
                    ns.TRef = _params.T0;
                    ns.IncludeThermal = true;
                    // 初始温度为参考温度
                    for (int i = 0; i < _mesh.NumNodes; i++) _mesh.Temperature[i] = _params.T0;
                    ns.Log = Log;
                    NSResult r = ns.Solve(_mesh);
                    BeginInvoke(new Action(delegate
                    {
                        Log(string.Format("Boussinesq 完成: outer={0}, resU={1:E3}", r.OuterIterations, r.FinalResU));
                        _fieldCombo.SelectedItem = "温度 T";
                        RefreshPlot(); _status.Text = "完成";
                    }));
                }
                catch (Exception ex) { BeginInvoke(new Action(() => { Log("Boussinesq 失败: " + ex.Message); _status.Text = "错误"; })); }
            });
        }

        private void RefreshPlot()
        {
            if (_mesh == null) return;
            try
            {
                string tmp = Path.Combine(Path.GetTempPath(), "fem2d_preview.png");
                ContourRenderer.RenderToPng(_mesh, tmp, 1200, 900,
                    (string)_fieldCombo.SelectedItem, _params.DeformScale);
                using (FileStream fs = new FileStream(tmp, FileMode.Open, FileAccess.Read))
                {
                    Image old = _pb.Image;
                    _pb.Image = Image.FromStream(fs);
                    if (old != null) old.Dispose();
                }
            }
            catch (Exception ex) { Log("绘图失败: " + ex.Message); }
        }

        private void ExportResults()
        {
            if (_mesh == null) { MessageBox.Show("请先生成或导入网格"); return; }
            using (FolderBrowserDialog dlg = new FolderBrowserDialog())
            {
                dlg.Description = "选择导出目录";
                if (dlg.ShowDialog() == DialogResult.OK)
                {
                    try
                    {
                        string field = (string)_fieldCombo.SelectedItem;
                        ResultWriter.WriteAll(_mesh, dlg.SelectedPath, "fem2d", field, _params.DeformScale);
                        Log("已导出到: " + dlg.SelectedPath);
                    }
                    catch (Exception ex) { Log("导出失败: " + ex.Message); }
                }
            }
        }
    }
}
