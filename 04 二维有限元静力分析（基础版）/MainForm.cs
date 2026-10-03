using System;
using System.Drawing;
using System.IO;
using System.Windows.Forms;
using FEM2D.IO;
using FEM2D.Mesh;
using FEM2D.Rendering;
using FEM2D.Solvers;

namespace FEM2D
{
    public sealed class MainForm : Form
    {
        private FEMesh _mesh;
        private StructuralResult _structResult;
        private NSResult _nsResult;
        private object _currentResult;

        private PictureBox _pictureBox;
        private TextBox _logBox;
        private ComboBox _fieldCombo;
        private PropertyGrid _paramGrid;
        private SplitContainer _split;
        private SplitContainer _splitRight;
        private ToolStripMenuItem _menuExport;
        private Label _statusLabel;

        private FEMParams _params = new FEMParams();

        public MainForm()
        {
            Text = "二维有限元计算平台 FEM2D (C# / .NET 10.0 / WinForms)";
            Width = 1280; Height = 820;
            StartPosition = FormStartPosition.CenterScreen;

            BuildMenu();
            BuildLayout();
        }

        private void BuildMenu()
        {
            MenuStrip ms = new MenuStrip();
            ToolStripMenuItem fileMenu = new ToolStripMenuItem("文件(&F)");
            ToolStripMenuItem openItem = new ToolStripMenuItem("导入网格 INP...", null, OnOpenInp);
            _menuExport = new ToolStripMenuItem("导出结果...", null, OnExport) { Enabled = false };
            ToolStripMenuItem exitItem = new ToolStripMenuItem("退出", null, (s, e) => Close());
            fileMenu.DropDownItems.Add(openItem);
            fileMenu.DropDownItems.Add(_menuExport);
            fileMenu.DropDownItems.Add(new ToolStripSeparator());
            fileMenu.DropDownItems.Add(exitItem);

            ToolStripMenuItem meshMenu = new ToolStripMenuItem("网格(&M)");
            ToolStripMenuItem genCst = new ToolStripMenuItem("自动剖分 - CST3 三角单元", null, (s, e) => GenerateMesh(ElementType.CST3));
            ToolStripMenuItem genLt6 = new ToolStripMenuItem("自动剖分 - LT6 二次三角单元", null, (s, e) => GenerateMesh(ElementType.LT6));
            ToolStripMenuItem genQ4 = new ToolStripMenuItem("自动剖分 - Q4 四边形单元", null, (s, e) => GenerateMesh(ElementType.Q4));
            meshMenu.DropDownItems.Add(genCst);
            meshMenu.DropDownItems.Add(genLt6);
            meshMenu.DropDownItems.Add(genQ4);

            ToolStripMenuItem solveMenu = new ToolStripMenuItem("求解(&S)");
            ToolStripMenuItem structuralItem = new ToolStripMenuItem("结构静力分析 (CG)", null, (s, e) => RunStructural(false));
            ToolStripMenuItem structuralGmres = new ToolStripMenuItem("结构静力分析 (GMRES)", null, (s, e) => RunStructural(true));
            ToolStripMenuItem nsSteadyItem = new ToolStripMenuItem("NS 定常 (SIMPLE)", null, (s, e) => RunNS(false));
            ToolStripMenuItem nsUnsteadyItem = new ToolStripMenuItem("NS 非稳态 (SIMPLE + 时间项)", null, (s, e) => RunNS(true));
            solveMenu.DropDownItems.Add(structuralItem);
            solveMenu.DropDownItems.Add(structuralGmres);
            solveMenu.DropDownItems.Add(new ToolStripSeparator());
            solveMenu.DropDownItems.Add(nsSteadyItem);
            solveMenu.DropDownItems.Add(nsUnsteadyItem);

            ToolStripMenuItem helpMenu = new ToolStripMenuItem("帮助(&H)");
            ToolStripMenuItem aboutItem = new ToolStripMenuItem("关于", null, (s, e) =>
                MessageBox.Show(
                    "二维有限元计算平台 FEM2D\n" +
                    "支持单元：CST3 / LT6 / Q4\n" +
                    "求解器：CG / GMRES；NS 方程 SIMPLE 算法\n" +
                    "稀疏矩阵 CSR 并行组装；.NET 10.0 WinForms\n",
                    "关于 FEM2D", MessageBoxButtons.OK, MessageBoxIcon.Information));
            helpMenu.DropDownItems.Add(aboutItem);

            ms.Items.Add(fileMenu);
            ms.Items.Add(meshMenu);
            ms.Items.Add(solveMenu);
            ms.Items.Add(helpMenu);
            MainMenuStrip = ms;
            Controls.Add(ms);
        }

        private void BuildLayout()
        {
            _split = new SplitContainer { Dock = DockStyle.Fill, Orientation = Orientation.Vertical, SplitterDistance = 30 };
            _splitRight = new SplitContainer { Dock = DockStyle.Fill, Orientation = Orientation.Horizontal, SplitterDistance = 560 };

            // 左侧参数面板
            Panel left = new Panel { Dock = DockStyle.Fill, Padding = new Padding(6) };
            Label l1 = new Label { Text = "参数设置", Dock = DockStyle.Top, Font = new Font("Arial", 10, FontStyle.Bold), Height = 24 };
            _paramGrid = new PropertyGrid { Dock = DockStyle.Fill, SelectedObject = _params, ToolbarVisible = false };
            _fieldCombo = new ComboBox { Dock = DockStyle.Bottom, DropDownStyle = ComboBoxStyle.DropDownList };
            _fieldCombo.Items.AddRange(new object[]
            {
                "网格", "位移大小", "变形图", "σxx", "σyy", "σxy", "Von Mises 应力", "速度大小", "压力"
            });
            _fieldCombo.SelectedIndex = 0;
            _fieldCombo.SelectedIndexChanged += (s, e) => RenderField();
            Label l2 = new Label { Text = "显示字段：", Dock = DockStyle.Bottom, Height = 20 };
            left.Controls.Add(_paramGrid);
            left.Controls.Add(l1);
            left.Controls.Add(l2);
            left.Controls.Add(_fieldCombo);

            // 右上图形区
            _pictureBox = new PictureBox { Dock = DockStyle.Fill, BackColor = Color.White, SizeMode = PictureBoxSizeMode.Zoom };
            Panel top = new Panel { Dock = DockStyle.Fill };
            top.Controls.Add(_pictureBox);

            // 右下日志区
            _logBox = new TextBox { Dock = DockStyle.Fill, Multiline = true, ScrollBars = ScrollBars.Vertical, Font = new Font("Consolas", 9), ReadOnly = true, BackColor = Color.Black, ForeColor = Color.Lime };
            Panel logPanel = new Panel { Dock = DockStyle.Fill, Padding = new Padding(3) };
            Label logLabel = new Label { Text = "计算日志", Dock = DockStyle.Top, Height = 20, Font = new Font("Arial", 9, FontStyle.Bold) };
            logPanel.Controls.Add(_logBox);
            logPanel.Controls.Add(logLabel);

            _splitRight.Panel1.Controls.Add(top);
            _splitRight.Panel2.Controls.Add(logPanel);

            _split.Panel1.Controls.Add(left);
            _split.Panel2.Controls.Add(_splitRight);

            _statusLabel = new Label { Dock = DockStyle.Bottom, Height = 22, Text = "就绪", BorderStyle = BorderStyle.FixedSingle };
            Controls.Add(_statusLabel);
            Controls.Add(_split);
        }

        private void Log(string msg)
        {
            if (InvokeRequired) { Invoke(new Action<string>(Log), msg); return; }
            _logBox.AppendText($"[{DateTime.Now:HH:mm:ss}] {msg}{Environment.NewLine}");
            _statusLabel.Text = msg;
        }

        private void OnOpenInp(object sender, EventArgs e)
        {
            using (OpenFileDialog dlg = new OpenFileDialog())
            {
                dlg.Filter = "Abaqus INP 文件 (*.inp)|*.inp|所有文件 (*.*)|*.*";
                if (dlg.ShowDialog() == DialogResult.OK)
                {
                    try
                    {
                        _mesh = MeshReader.ReadInp(dlg.FileName);
                        _structResult = null; _nsResult = null; _currentResult = null;
                        Log($"已导入网格：{Path.GetFileName(dlg.FileName)}，节点={_mesh.NumNodes}, 单元={_mesh.NumElements}");
                        RenderField();
                    }
                    catch (Exception ex)
                    {
                        MessageBox.Show("导入失败：" + ex.Message, "错误", MessageBoxButtons.OK, MessageBoxIcon.Error);
                    }
                }
            }
        }

        private void GenerateMesh(ElementType type)
        {
            try
            {
                bool ns = _params.ProblemType == ProblemTypeEnum.NS;
                if (ns)
                {
                    _mesh = StructuredMeshGenerator.GenerateLidDrivenCavity(_params.LengthX, _params.Nx, type, _params.LidVelocity);
                    Log($"已生成 NS 顶盖驱动方腔网格：{type}, {_mesh.NumNodes} 节点 / {_mesh.NumElements} 单元");
                }
                else
                {
                    _mesh = StructuredMeshGenerator.GenerateRectangle(
                        _params.LengthX, _params.LengthY, _params.Nx, _params.Ny, type,
                        _params.FixLeft, _params.TopPressure, _params.EndLoad);
                    Log($"已生成结构网格：{type}, {_mesh.NumNodes} 节点 / {_mesh.NumElements} 单元");
                }
                _structResult = null; _nsResult = null; _currentResult = null;
                _fieldCombo.SelectedIndex = 0;
                RenderField();
            }
            catch (Exception ex)
            {
                MessageBox.Show("网格生成失败：" + ex.Message, "错误", MessageBoxButtons.OK, MessageBoxIcon.Error);
                Log("异常：" + ex.Message);
            }
        }

        private void RunStructural(bool useGmres)
        {
            if (_mesh == null) { MessageBox.Show("请先生成或导入网格。", "提示", MessageBoxButtons.OK, MessageBoxIcon.Warning); return; }
            _menuExport.Enabled = false;
            _params.ProblemType = ProblemTypeEnum.Structural;
            _structResult = null; _nsResult = null; _currentResult = null;
            Log("======== 开始结构静力分析 ========");
            System.Threading.ThreadPool.QueueUserWorkItem(_ =>
            {
                try
                {
                    StructuralSolver solver = new StructuralSolver
                    {
                        YoungModulus = _params.E,
                        PoissonRatio = _params.Nu,
                        Thickness = _params.Thickness,
                        UseGMRES = useGmres
                    };
                    StructuralResult r = solver.Solve(_mesh, msg => Log(msg));
                    _structResult = r;
                    _currentResult = r;
                    Log($"分析完成！最大 Von Mises = {r.MaxVM:G4} Pa, 最大位移 = {r.MaxDisp:G4} m, 迭代 {solver.Iterations}");
                    BeginInvoke(new Action(() =>
                    {
                        _menuExport.Enabled = true;
                        _fieldCombo.SelectedIndex = 6; // Von Mises
                        RenderField();
                    }));
                }
                catch (Exception ex)
                {
                    Log("求解失败：" + ex.Message);
                    BeginInvoke(new Action(() => MessageBox.Show("求解失败：" + ex.Message, "错误", MessageBoxButtons.OK, MessageBoxIcon.Error)));
                }
            });
        }

        private void RunNS(bool unsteady)
        {
            if (_mesh == null) { MessageBox.Show("请先生成或导入网格。", "提示", MessageBoxButtons.OK, MessageBoxIcon.Warning); return; }
            _menuExport.Enabled = false;
            _params.ProblemType = ProblemTypeEnum.NS;
            _structResult = null; _nsResult = null; _currentResult = null;
            Log(unsteady ? "======== 开始 NS 非稳态计算 ========" : "======== 开始 NS 定常 SIMPLE 计算 ========");
            System.Threading.ThreadPool.QueueUserWorkItem(_ =>
            {
                try
                {
                    NavierStokesSolver solver = new NavierStokesSolver
                    {
                        Density = _params.Rho,
                        Viscosity = _params.Mu,
                        MaxOuterIter = _params.MaxOuterIter,
                        ConvergenceTol = _params.NSTol,
                        TimeStep = unsteady ? _params.Dt : 0.0,
                        TimeSteps = unsteady ? _params.NTimeSteps : 1,
                        AlphaU = _params.AlphaU,
                        AlphaP = _params.AlphaP
                    };
                    NSResult r = solver.Solve(_mesh, msg => Log(msg));
                    _nsResult = r;
                    _currentResult = r;
                    Log($"NS 完成！外迭代次数 = {r.Iterations}, 最终残差 = {r.FinalResidual:E4}");
                    BeginInvoke(new Action(() =>
                    {
                        _menuExport.Enabled = true;
                        _fieldCombo.SelectedIndex = 7; // Velocity
                        RenderField();
                    }));
                }
                catch (Exception ex)
                {
                    Log("求解失败：" + ex.Message);
                    BeginInvoke(new Action(() => MessageBox.Show("求解失败：" + ex.Message, "错误", MessageBoxButtons.OK, MessageBoxIcon.Error)));
                }
            });
        }

        private void RenderField()
        {
            if (_mesh == null) return;
            FieldType ft;
            switch (_fieldCombo.SelectedIndex)
            {
                case 0: ft = FieldType.Mesh; break;
                case 1: ft = FieldType.Displacement; break;
                case 2: ft = FieldType.Deformed; break;
                case 3: ft = FieldType.StressXX; break;
                case 4: ft = FieldType.StressYY; break;
                case 5: ft = FieldType.StressXY; break;
                case 6: ft = FieldType.VonMises; break;
                case 7: ft = FieldType.Velocity; break;
                case 8: ft = FieldType.Pressure; break;
                default: ft = FieldType.Mesh; break;
            }
            if ((ft == FieldType.StressXX || ft == FieldType.StressYY || ft == FieldType.StressXY || ft == FieldType.VonMises || ft == FieldType.Displacement || ft == FieldType.Deformed) && _structResult == null)
            { ft = FieldType.Mesh; }
            if ((ft == FieldType.Velocity || ft == FieldType.Pressure) && _nsResult == null)
            { ft = FieldType.Mesh; }

            int w = Math.Max(400, _pictureBox.Width);
            int h = Math.Max(400, _pictureBox.Height);
            try
            {
                Bitmap bmp = ContourRenderer.Render(_mesh, _currentResult, ft, w, h);
                Image old = _pictureBox.Image;
                _pictureBox.Image = bmp;
                if (old != null) old.Dispose();
            }
            catch (Exception ex)
            {
                Log("渲染失败：" + ex.Message);
            }
        }

        private void OnExport(object sender, EventArgs e)
        {
            if (_currentResult == null) { MessageBox.Show("没有可导出的结果。", "提示"); return; }
            using (FolderBrowserDialog dlg = new FolderBrowserDialog())
            {
                dlg.Description = "选择结果导出目录";
                if (dlg.ShowDialog() == DialogResult.OK)
                {
                    string dir = dlg.SelectedPath;
                    try
                    {
                        if (_currentResult is StructuralResult sr)
                        {
                            ResultWriter.WriteStructuralCsv(_mesh, sr, dir);
                            ResultWriter.WriteVtk(_mesh, sr, Path.Combine(dir, "structural.vtk"), "Structural");
                            Bitmap bmp = ContourRenderer.Render(_mesh, sr, FieldType.VonMises, 1600, 1200);
                            bmp.Save(Path.Combine(dir, "vonmises.png"), System.Drawing.Imaging.ImageFormat.Png);
                            Bitmap bmp2 = ContourRenderer.Render(_mesh, sr, FieldType.Deformed, 1600, 1200);
                            bmp2.Save(Path.Combine(dir, "deformed.png"), System.Drawing.Imaging.ImageFormat.Png);
                            Log($"结果已导出到：{dir}");
                        }
                        else if (_currentResult is NSResult nr)
                        {
                            ResultWriter.WriteNSCsv(_mesh, nr, dir);
                            ResultWriter.WriteVtk(_mesh, nr, Path.Combine(dir, "ns.vtk"), "NS");
                            Bitmap bmp = ContourRenderer.Render(_mesh, nr, FieldType.Velocity, 1600, 1200);
                            bmp.Save(Path.Combine(dir, "velocity.png"), System.Drawing.Imaging.ImageFormat.Png);
                            Bitmap bmp2 = ContourRenderer.Render(_mesh, nr, FieldType.Pressure, 1600, 1200);
                            bmp2.Save(Path.Combine(dir, "pressure.png"), System.Drawing.Imaging.ImageFormat.Png);
                            Log($"结果已导出到：{dir}");
                        }
                        MessageBox.Show("结果导出成功！\n目录：" + dir, "导出完成", MessageBoxButtons.OK, MessageBoxIcon.Information);
                    }
                    catch (Exception ex)
                    {
                        MessageBox.Show("导出失败：" + ex.Message, "错误", MessageBoxButtons.OK, MessageBoxIcon.Error);
                    }
                }
            }
        }
    }

    public enum ProblemTypeEnum { Structural, NS }

    public sealed class FEMParams
    {
        [System.ComponentModel.Category("几何")]
        [System.ComponentModel.Description("矩形域 X 向长度")]
        public double LengthX { get; set; } = 10.0;
        [System.ComponentModel.Category("几何")]
        public double LengthY { get; set; } = 1.0;
        [System.ComponentModel.Category("网格")]
        public int Nx { get; set; } = 30;
        [System.ComponentModel.Category("网格")]
        public int Ny { get; set; } = 4;
        [System.ComponentModel.Category("网格")]
        public ProblemTypeEnum ProblemType { get; set; } = ProblemTypeEnum.Structural;

        [System.ComponentModel.Category("材料-结构")]
        public double E { get; set; } = 210e9;
        [System.ComponentModel.Category("材料-结构")]
        public double Nu { get; set; } = 0.3;
        [System.ComponentModel.Category("材料-结构")]
        public double Thickness { get; set; } = 0.1;
        [System.ComponentModel.Category("荷载")]
        public bool FixLeft { get; set; } = true;
        [System.ComponentModel.Category("荷载")]
        public double TopPressure { get; set; } = 0.0;
        [System.ComponentModel.Category("荷载")]
        public double EndLoad { get; set; } = -1e5;

        [System.ComponentModel.Category("流体-NS")]
        public double Rho { get; set; } = 1.0;
        [System.ComponentModel.Category("流体-NS")]
        public double Mu { get; set; } = 0.01;
        [System.ComponentModel.Category("流体-NS")]
        public double LidVelocity { get; set; } = 1.0;
        [System.ComponentModel.Category("SIMPLE")]
        public double AlphaU { get; set; } = 0.7;
        [System.ComponentModel.Category("SIMPLE")]
        public double AlphaP { get; set; } = 0.3;
        [System.ComponentModel.Category("SIMPLE")]
        public int MaxOuterIter { get; set; } = 200;
        [System.ComponentModel.Category("SIMPLE")]
        public double NSTol { get; set; } = 1e-4;
        [System.ComponentModel.Category("非稳态")]
        public double Dt { get; set; } = 0.01;
        [System.ComponentModel.Category("非稳态")]
        public int NTimeSteps { get; set; } = 20;
    }
}
