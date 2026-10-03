using System;
using System.Drawing;
using System.Drawing.Drawing2D;
using System.IO;
using System.Text;
using System.Runtime;
using System.ComponentModel;
using System.Windows.Forms;
using Fem2DFluid.Models;
using Fem2DFluid.Services;

namespace Fem2DFluid.Forms
{
    public class MainForm : Form
    {
        private MenuStrip _menu;
        private ToolStripMenuItem _miFile, _miOpenTxt, _miOpenDxf, _miExportCsv,
            _miExportTxt, _miExportVtk, _miExportPng, _miExit;
        private ToolStripMenuItem _miSample, _miSamp1, _miSamp2, _miSamp3;
        private ToolStripMenuItem _miHelp, _miAbout;

        private ToolStrip _tool;
        private ToolStripButton _btnOpen, _btnSave;
        private ToolStripComboBox _cbSample;
        private ToolStripTextBox _tbH;
        private ToolStripButton _btnMesh;
        private ToolStripCheckBox _chkLT6;
        private ToolStripCheckBox _chkNS;
        private ToolStripTextBox _tbDt;
        private ToolStripButton _btnSolve, _btnTimeStep, _btnExport;

        private SplitContainer _split;
        private TreeView _tree;
        private TableLayoutPanel _rightPanel;
        private PictureBox _picture;
        private RichTextBox _logBox;
        private StatusStrip _status;
        private ToolStripStatusLabel _slMain;

        private FemModel _model;
        private FemResult _lastResult;
        private string _currentField;

        public MainForm()
        {
            Text = "Fem2DFluid - 2D Finite Element Fluid";
            Width = 1200;
            Height = 800;
            StartPosition = FormStartPosition.CenterScreen;
            BuildMenu();
            BuildToolstrip();
            BuildBody();
            BuildStatus();
            MainMenuStrip = _menu;
            _model = new FemModel();
            _currentField = "phi";
            Log("Ready. Use File -> Open or sample menu to begin.");
        }

        private void BuildMenu()
        {
            _menu = new MenuStrip();
            _miFile = new ToolStripMenuItem("文件(&F)");
            _miOpenTxt = new ToolStripMenuItem("打开TXT...");
            _miOpenDxf = new ToolStripMenuItem("打开DXF...");
            _miExportCsv = new ToolStripMenuItem("导出CSV...");
            _miExportTxt = new ToolStripMenuItem("导出TXT...");
            _miExportVtk = new ToolStripMenuItem("导出VTK...");
            _miExportPng = new ToolStripMenuItem("导出PNG...");
            _miExit = new ToolStripMenuItem("退出");
            _miFile.DropDownItems.Add(_miOpenTxt);
            _miFile.DropDownItems.Add(_miOpenDxf);
            _miFile.DropDownItems.Add(new ToolStripSeparator());
            _miFile.DropDownItems.Add(_miExportCsv);
            _miFile.DropDownItems.Add(_miExportTxt);
            _miFile.DropDownItems.Add(_miExportVtk);
            _miFile.DropDownItems.Add(_miExportPng);
            _miFile.DropDownItems.Add(new ToolStripSeparator());
            _miFile.DropDownItems.Add(_miExit);

            _miSample = new ToolStripMenuItem("示例(&S)");
            _miSamp1 = new ToolStripMenuItem("示例1 方腔势流");
            _miSamp2 = new ToolStripMenuItem("示例2 圆柱绕流");
            _miSamp3 = new ToolStripMenuItem("示例3 坝体渗流");
            _miSample.DropDownItems.Add(_miSamp1);
            _miSample.DropDownItems.Add(_miSamp2);
            _miSample.DropDownItems.Add(_miSamp3);

            _miHelp = new ToolStripMenuItem("帮助(&H)");
            _miAbout = new ToolStripMenuItem("关于");
            _miHelp.DropDownItems.Add(_miAbout);

            _menu.Items.Add(_miFile);
            _menu.Items.Add(_miSample);
            _menu.Items.Add(_miHelp);

            _miOpenTxt.Click += delegate (object s, EventArgs e) { OpenTxt(); };
            _miOpenDxf.Click += delegate (object s, EventArgs e) { OpenDxf(); };
            _miExportCsv.Click += delegate (object s, EventArgs e) { ExportCsv(); };
            _miExportTxt.Click += delegate (object s, EventArgs e) { ExportTxt(); };
            _miExportVtk.Click += delegate (object s, EventArgs e) { ExportVtk(); };
            _miExportPng.Click += delegate (object s, EventArgs e) { ExportPng(); };
            _miExit.Click += delegate (object s, EventArgs e) { Close(); };
            _miSamp1.Click += delegate (object s, EventArgs e) { LoadSample(1); };
            _miSamp2.Click += delegate (object s, EventArgs e) { LoadSample(2); };
            _miSamp3.Click += delegate (object s, EventArgs e) { LoadSample(3); };
            _miAbout.Click += delegate (object s, EventArgs e)
            {
                MessageBox.Show("Fem2DFluid v1.0\n2D FEM fluid solver (C# 3.0 / WinForm .NET 10)\n"
                    + "Potential flow, Darcy seepage, SIMPLE Navier-Stokes.",
                    "About", MessageBoxButtons.OK, MessageBoxIcon.Information);
            };
        }

        private void BuildToolstrip()
        {
            _tool = new ToolStrip();
            _btnOpen = new ToolStripButton("打开");
            _btnSave = new ToolStripButton("保存");
            _cbSample = new ToolStripComboBox();
            _cbSample.Items.AddRange(new object[] { "示例1 方腔", "示例2 圆柱", "示例3 渗流" });
            _cbSample.SelectedIndex = 0;
            _cbSample.DropDownStyle = ComboBoxStyle.DropDownList;
            ToolStripLabel lblH = new ToolStripLabel("网格h:");
            _tbH = new ToolStripTextBox();
            _tbH.Text = "0.05";
            _tbH.Width = 50;
            _btnMesh = new ToolStripButton("剖分");
            _chkLT6 = new ToolStripCheckBox("LT6二次");
            _chkNS = new ToolStripCheckBox("NS流动");
            ToolStripLabel lblDt = new ToolStripLabel("Δt:");
            _tbDt = new ToolStripTextBox();
            _tbDt.Text = "0.01";
            _tbDt.Width = 50;
            _btnSolve = new ToolStripButton("求解");
            _btnTimeStep = new ToolStripButton("时间步");
            _btnExport = new ToolStripButton("导出");

            _tool.Items.Add(_btnOpen);
            _tool.Items.Add(_btnSave);
            _tool.Items.Add(new ToolStripSeparator());
            _tool.Items.Add(_cbSample);
            _tool.Items.Add(lblH);
            _tool.Items.Add(_tbH);
            _tool.Items.Add(_btnMesh);
            _tool.Items.Add(_chkLT6);
            _tool.Items.Add(_chkNS);
            _tool.Items.Add(lblDt);
            _tool.Items.Add(_tbDt);
            _tool.Items.Add(new ToolStripSeparator());
            _tool.Items.Add(_btnSolve);
            _tool.Items.Add(_btnTimeStep);
            _tool.Items.Add(_btnExport);

            _btnOpen.Click += delegate (object s, EventArgs e) { OpenTxt(); };
            _btnMesh.Click += delegate (object s, EventArgs e) { DoMesh(); };
            _btnSolve.Click += delegate (object s, EventArgs e) { DoSolve(); };
            _btnTimeStep.Click += delegate (object s, EventArgs e) { DoTimeStep(); };
            _btnExport.Click += delegate (object s, EventArgs e) { ExportAll(); };
        }

        private void BuildBody()
        {
            _split = new SplitContainer();
            _split.Dock = DockStyle.Fill;
            _split.Orientation = Orientation.Vertical;
            _split.SplitterDistance = 240;
            _tree = new TreeView();
            _tree.Dock = DockStyle.Fill;
            _split.Panel1.Controls.Add(_tree);

            _rightPanel = new TableLayoutPanel();
            _rightPanel.Dock = DockStyle.Fill;
            _rightPanel.RowCount = 2;
            _rightPanel.ColumnCount = 1;
            _rightPanel.RowStyles.Add(new RowStyle(SizeType.Percent, 70));
            _rightPanel.RowStyles.Add(new RowStyle(SizeType.Percent, 30));
            _picture = new PictureBox();
            _picture.Dock = DockStyle.Fill;
            _picture.BackColor = Color.White;
            _picture.SizeMode = PictureBoxSizeMode.Zoom;
            _logBox = new RichTextBox();
            _logBox.Dock = DockStyle.Fill;
            _logBox.BackColor = Color.Black;
            _logBox.ForeColor = Color.LightGreen;
            _logBox.Font = new Font("Consolas", 9);
            _logBox.ReadOnly = true;
            _rightPanel.Controls.Add(_picture, 0, 0);
            _rightPanel.Controls.Add(_logBox, 0, 1);
            _split.Panel2.Controls.Add(_rightPanel);

            Controls.Add(_split);
            Controls.Add(_tool);
            Controls.Add(_menu);
        }

        private void BuildStatus()
        {
            _status = new StatusStrip();
            _slMain = new ToolStripStatusLabel("Ready");
            _slMain.Spring = true;
            _slMain.TextAlign = ContentAlignment.MiddleLeft;
            _status.Items.Add(_slMain);
            Controls.Add(_status);
        }

        // ========== Actions ==========

        private void Log(string msg)
        {
            if (_logBox == null) return;
            _logBox.AppendText(string.Format("[{0:HH:mm:ss}] {1}\n", DateTime.Now, msg));
            _slMain.Text = msg;
        }

        private void OpenTxt()
        {
            using (OpenFileDialog dlg = new OpenFileDialog())
            {
                dlg.Filter = "Text Mesh (*.txt)|*.txt|All (*.*)|*.*";
                if (dlg.ShowDialog(this) != DialogResult.OK) return;
                _model.ResetResults();
                MeshImporter mi = new MeshImporter(_model);
                mi.FromTextFile(dlg.FileName);
                Log(string.Format("Loaded TXT: {0} nodes, {1} elements", _model.Nodes.Count, _model.Elements.Count));
                RefreshTree();
                DrawResult();
            }
        }

        private void OpenDxf()
        {
            using (OpenFileDialog dlg = new OpenFileDialog())
            {
                dlg.Filter = "DXF (*.dxf)|*.dxf|All (*.*)|*.*";
                if (dlg.ShowDialog(this) != DialogResult.OK) return;
                _model.ResetResults();
                MeshImporter mi = new MeshImporter(_model);
                bool ok = mi.FromDxf(dlg.FileName);
                Log(ok ? "DXF geometry loaded." : "DXF load failed; using fallback mesh.");
                RefreshTree();
                DrawResult();
            }
        }

        private void LoadSample(int idx)
        {
            _model = new FemModel();
            SampleBuilder sb = new SampleBuilder(_model);
            bool doMeshDirect = false;
            switch (idx)
            {
                case 1:
                    sb.BuildSquareCavity(41, 21);
                    _currentField = "phi";
                    doMeshDirect = true;
                    Log("Loaded sample 1: square cavity (structured, no re-mesh needed).");
                    break;
                case 2:
                    sb.BuildCylinderFlow(100, 40);
                    _model.IsNavierStokes = _chkNS.Checked;
                    sb.ApplyCylinderBC();
                    _currentField = _chkNS.Checked ? "vmag" : "phi";
                    doMeshDirect = true;
                    Log("Loaded sample 2: flow around cylinder.");
                    break;
                case 3:
                    sb.BuildDamSeepage(41, 21);
                    _currentField = "phi";
                    doMeshDirect = true;
                    Log("Loaded sample 3: dam seepage.");
                    break;
            }
            _cbSample.SelectedIndex = idx - 1;
            if (doMeshDirect)
            {
                DoMeshDirect();
            }
            else
            {
                DoMesh();
            }
        }

        private void DoMesh()
        {
            double h;
            if (!double.TryParse(_tbH.Text, System.Globalization.NumberStyles.Float,
                System.Globalization.CultureInfo.InvariantCulture, out h) || h <= 0)
            {
                MessageBox.Show("Invalid mesh size h."); return;
            }
            System.Diagnostics.Stopwatch sw = new System.Diagnostics.Stopwatch();
            sw.Start();
            Triangulator tr = new Triangulator(_model);
            tr.SubdivideBoundaryEdges(h);
            tr.GenerateStructuredPoints(h);
            tr.Delaunay();
            int removed = tr.RemoveOutsideElements();
            tr.ComputeElementGeoms();
            tr.RenumberNodesFromZero();
            if (_chkLT6.Checked)
            {
                tr.AddMidEdgeNodes();
            }
            else
            {
                tr.RemoveMidEdgeNodes();
            }
            _model.ComputeBounds();
            sw.Stop();
            Log(string.Format("Meshing done: {0} nodes, {1} elements, removed={2}, {3}ms",
                _model.Nodes.Count, _model.Elements.Count, removed, sw.ElapsedMilliseconds));
            RefreshTree();
            DrawResult();
        }

        private void DoMeshDirect()
        {
            Triangulator tr = new Triangulator(_model);
            tr.ComputeElementGeoms();
            if (_chkLT6.Checked)
            {
                tr.AddMidEdgeNodes();
            }
            _model.ComputeBounds();
            RefreshTree();
            DrawResult();
        }

        private void DoSolve()
        {
            FemSolver solver = new FemSolver(_model);
            _model.IsNavierStokes = _chkNS.Checked;
            if (_chkNS.Checked)
            {
                SampleBuilder sb = new SampleBuilder(_model);
                sb.ApplyCylinderBC();
                _lastResult = solver.SolveSIMPLE();
                _currentField = "vmag";
            }
            else
            {
                _lastResult = solver.SolveSteady(1.0, _chkLT6.Checked);
                _currentField = "phi";
            }
            Log(_lastResult == null ? "Solve failed." : _lastResult.ToString());
            RefreshTree();
            DrawResult();
        }

        private void DoTimeStep()
        {
            double dt;
            if (!double.TryParse(_tbDt.Text, System.Globalization.NumberStyles.Float,
                System.Globalization.CultureInfo.InvariantCulture, out dt) || dt <= 0)
            {
                MessageBox.Show("Invalid dt."); return;
            }
            FemSolver solver = new FemSolver(_model);
            _lastResult = solver.AdvanceTime(dt, 0.01);
            Log(_lastResult == null ? "Time step failed." : _lastResult.ToString());
            DrawResult();
        }

        private void ExportCsv()
        {
            using (SaveFileDialog dlg = new SaveFileDialog())
            {
                dlg.Filter = "CSV (*.csv)|*.csv";
                dlg.FileName = "result.csv";
                if (dlg.ShowDialog(this) != DialogResult.OK) return;
                ResultExporter ex = new ResultExporter(_model);
                ex.ExportCsv(dlg.FileName);
                Log("CSV exported: " + dlg.FileName);
            }
        }

        private void ExportTxt()
        {
            using (SaveFileDialog dlg = new SaveFileDialog())
            {
                dlg.Filter = "Text (*.txt)|*.txt";
                dlg.FileName = "result.txt";
                if (dlg.ShowDialog(this) != DialogResult.OK) return;
                ResultExporter ex = new ResultExporter(_model);
                ex.ExportText(dlg.FileName);
                Log("TXT exported: " + dlg.FileName);
            }
        }

        private void ExportVtk()
        {
            using (SaveFileDialog dlg = new SaveFileDialog())
            {
                dlg.Filter = "VTK (*.vtk)|*.vtk";
                dlg.FileName = "result.vtk";
                if (dlg.ShowDialog(this) != DialogResult.OK) return;
                FemSolver solver = new FemSolver(_model);
                solver.ExportTextResults(null, null, dlg.FileName);
                Log("VTK exported: " + dlg.FileName);
            }
        }

        private void ExportPng()
        {
            using (SaveFileDialog dlg = new SaveFileDialog())
            {
                dlg.Filter = "PNG (*.png)|*.png";
                dlg.FileName = "result.png";
                if (dlg.ShowDialog(this) != DialogResult.OK) return;
                ResultExporter ex = new ResultExporter(_model);
                int w = Math.Max(800, _picture.Width);
                int h = Math.Max(600, _picture.Height);
                ex.SavePng(dlg.FileName, _currentField, w, h, true, true);
                Log("PNG exported: " + dlg.FileName);
            }
        }

        private void ExportAll()
        {
            using (FolderBrowserDialog dlg = new FolderBrowserDialog())
            {
                if (dlg.ShowDialog(this) != DialogResult.OK) return;
                string dir = dlg.SelectedPath;
                ResultExporter ex = new ResultExporter(_model);
                ex.ExportCsv(Path.Combine(dir, "result.csv"));
                ex.ExportText(Path.Combine(dir, "result.txt"));
                FemSolver solver = new FemSolver(_model);
                solver.ExportTextResults(null, null, Path.Combine(dir, "result.vtk"));
                ex.SavePng(Path.Combine(dir, "result.png"), _currentField, 1200, 800, true, true);
                Log("All exports written to " + dir);
            }
        }

        private void RefreshTree()
        {
            _tree.BeginUpdate();
            _tree.Nodes.Clear();
            TreeNode root = new TreeNode(string.Format("Model ({0})", _model.Name ?? ""));
            TreeNode nNodes = new TreeNode(string.Format("Nodes ({0})", _model.Nodes.Count));
            int bnd = 0, dir = 0;
            for (int i = 0; i < _model.Nodes.Count; i++)
            {
                Node n = _model.Nodes[i];
                nNodes.Nodes.Add(new TreeNode(string.Format(
                    "#{0} ({1:F3},{2:F3}) tag={3} bc={4}",
                    n.Id, n.X, n.Y, n.Tag, n.BCType)));
                if (n.Tag > 0) bnd++;
                if (n.BCType == 1) dir++;
            }
            root.Nodes.Add(nNodes);
            TreeNode nElem = new TreeNode(string.Format("Elements ({0})", _model.Elements.Count));
            int show = Math.Min(_model.Elements.Count, 200);
            for (int i = 0; i < show; i++)
            {
                Element e = _model.Elements[i];
                nElem.Nodes.Add(new TreeNode(string.Format(
                    "#{0} ({1},{2},{3}) area={4:F4}", e.Id, e.A, e.B, e.C, e.Area)));
            }
            root.Nodes.Add(nElem);
            TreeNode nMat = new TreeNode(string.Format("Materials ({0})", _model.Materials.Count));
            for (int i = 0; i < _model.Materials.Count; i++)
            {
                Material m = _model.Materials[i];
                nMat.Nodes.Add(new TreeNode(string.Format("{0} K={1}, mu={2}, rho={3}",
                    m.Name, m.K, m.Viscosity, m.Density)));
            }
            root.Nodes.Add(nMat);
            TreeNode nBnd = new TreeNode(string.Format("Boundary (edges={0}, Dirichlet nodes={1})",
                _model.BoundaryEdgeCache.Count, dir));
            root.Nodes.Add(nBnd);
            if (_lastResult != null)
            {
                TreeNode nr = new TreeNode("Last result");
                nr.Nodes.Add(new TreeNode("Converged: " + _lastResult.Converged));
                nr.Nodes.Add(new TreeNode("Iterations: " + _lastResult.Iterations));
                nr.Nodes.Add(new TreeNode("Assemble ms: " + _lastResult.AssembleMs));
                nr.Nodes.Add(new TreeNode("Solve ms: " + _lastResult.SolveMs));
                nr.Nodes.Add(new TreeNode("Post ms: " + _lastResult.PostMs));
                nr.Nodes.Add(new TreeNode(string.Format("Phi range: [{0:F4}, {1:F4}]", _lastResult.MinPhi, _lastResult.MaxPhi)));
                nr.Nodes.Add(new TreeNode(string.Format("|V| range: [{0:F4}, {1:F4}]", _lastResult.MinV, _lastResult.MaxV)));
                root.Nodes.Add(nr);
            }
            _tree.Nodes.Add(root);
            root.Expand();
            _tree.EndUpdate();
        }

        private void DrawResult()
        {
            if (_model == null || _model.Elements.Count == 0)
            {
                _picture.Image = null;
                return;
            }
            int w = Math.Max(400, _picture.ClientSize.Width);
            int h = Math.Max(300, _picture.ClientSize.Height);
            ResultExporter ex = new ResultExporter(_model);
            string tmp = Path.Combine(Path.GetTempPath(), "fem2d_view.png");
            ex.SavePng(tmp, _currentField, w, h, true, !_chkNS.Checked);
            Image old = _picture.Image;
            try
            {
                FileStream fs = new FileStream(tmp, FileMode.Open, FileAccess.Read);
                _picture.Image = Image.FromStream(fs);
                fs.Close();
            }
            catch { }
            if (old != null) { old.Dispose(); }
            try { File.Delete(tmp); } catch { }
        }

        protected override void OnResize(EventArgs e)
        {
            base.OnResize(e);
            if (_picture != null && _model != null && _model.Elements.Count > 0)
            {
                DrawResult();
            }
        }
    }

    /// <summary>Minimal ToolStripCheckBox wrapper (since base ToolStrip doesn't include CheckBox directly).</summary>
    public class ToolStripCheckBox : ToolStripControlHost
    {
        private CheckBox _cb;
        public ToolStripCheckBox(string text) : base(new CheckBox())
        {
            _cb = (CheckBox)Control;
            _cb.Text = text;
            _cb.AutoSize = true;
        }
        [DefaultValue(false)]
        [DesignerSerializationVisibility(DesignerSerializationVisibility.Hidden)]
        public bool Checked
        {
            get { return _cb.Checked; }
            set { _cb.Checked = value; }
        }
        public event EventHandler CheckedChanged
        {
            add { _cb.CheckedChanged += value; }
            remove { _cb.CheckedChanged -= value; }
        }
    }
}
