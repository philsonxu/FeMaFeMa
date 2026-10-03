using System;
using System.Collections.Generic;
using System.Drawing;
using System.IO;
using System.Windows.Forms;
using Fem2DFluid.Models;
using Fem2DFluid.Services;

namespace Fem2DFluid.Forms
{
    public class MainForm : Form
    {
        private FemModel _model;
        private FemResult _lastResult;
        private FemSolver _solver;

        private MenuStrip _menu;
        private ToolStrip _toolbar;
        private StatusStrip _status;
        private ToolStripStatusLabel _statusLabel;
        private SplitContainer _split;
        private TreeView _tree;
        private Panel _canvasPanel;
        private PictureBox _canvas;
        private ComboBox _cbField;
        private CheckBox _chkMesh;
        private CheckBox _chkVectors;
        private CheckBox _chkLT6;
        private NumericUpDown _numMeshSize;
        private NumericUpDown _numDt;
        private Button _btnMesh;
        private Button _btnSolve;
        private Button _btnTimeStep;
        private TextBox _logBox;

        public MainForm()
        {
            Text = "二维有限元流体分析 — Fem2DFluid (C# 3.0 / .NET 10)";
            Width = 1280;
            Height = 820;
            StartPosition = FormStartPosition.CenterScreen;
            BuildUi();
            // 默认载入示例1（方腔流）
            LoadSample1();
        }

        private void BuildUi()
        {
            _menu = new MenuStrip();
            ToolStripMenuItem fileMenu = new ToolStripMenuItem("文件(&F)");
            ToolStripMenuItem miOpenTxt = new ToolStripMenuItem("打开节点单元文件(&T)...");
            miOpenTxt.Click += delegate { OpenTextFile(); };
            ToolStripMenuItem miOpenDxf = new ToolStripMenuItem("导入DXF几何(&D)...");
            miOpenDxf.Click += delegate { OpenDxfFile(); };
            ToolStripMenuItem miExportCsv = new ToolStripMenuItem("导出结果CSV...");
            miExportCsv.Click += delegate { ExportCsv(); };
            ToolStripMenuItem miExportTec = new ToolStripMenuItem("导出Tecplot DAT...");
            miExportTec.Click += delegate { ExportTec(); };
            ToolStripMenuItem miExportTxt = new ToolStripMenuItem("导出文本报告/VTK...");
            miExportTxt.Click += delegate { ExportTextReport(); };
            ToolStripMenuItem miExportPng = new ToolStripMenuItem("导出云图PNG...");
            miExportPng.Click += delegate { ExportPng(); };
            ToolStripMenuItem miExit = new ToolStripMenuItem("退出");
            miExit.Click += delegate { Close(); };
            fileMenu.DropDownItems.Add(miOpenTxt);
            fileMenu.DropDownItems.Add(miOpenDxf);
            fileMenu.DropDownItems.Add(new ToolStripSeparator());
            fileMenu.DropDownItems.Add(miExportCsv);
            fileMenu.DropDownItems.Add(miExportTec);
            fileMenu.DropDownItems.Add(miExportTxt);
            fileMenu.DropDownItems.Add(miExportPng);
            fileMenu.DropDownItems.Add(new ToolStripSeparator());
            fileMenu.DropDownItems.Add(miExit);

            ToolStripMenuItem sampleMenu = new ToolStripMenuItem("示例(&S)");
            ToolStripMenuItem miS1 = new ToolStripMenuItem("示例1：方腔均匀势流");
            miS1.Click += delegate { LoadSample1(); };
            ToolStripMenuItem miS2 = new ToolStripMenuItem("示例2：绕圆柱势流");
            miS2.Click += delegate { LoadSample2(); };
            ToolStripMenuItem miS3 = new ToolStripMenuItem("示例3：达西渗流坝体");
            miS3.Click += delegate { LoadSample3(); };
            sampleMenu.DropDownItems.Add(miS1);
            sampleMenu.DropDownItems.Add(miS2);
            sampleMenu.DropDownItems.Add(miS3);

            ToolStripMenuItem helpMenu = new ToolStripMenuItem("帮助(&H)");
            ToolStripMenuItem miAbout = new ToolStripMenuItem("关于");
            miAbout.Click += delegate
            {
                MessageBox.Show("二维有限元流体分析程序\n" +
                                "C# 3.0 / .NET 10 / WinForm\n" +
                                "线性三角单元(CST) + Bowyer-Watson Delaunay剖分\n" +
                                "支持势流/渗流/导热的二维椭圆方程\n" +
                                "∇·(k∇φ)=f", "关于", MessageBoxButtons.OK, MessageBoxIcon.Information);
            };
            helpMenu.DropDownItems.Add(miAbout);

            _menu.Items.Add(fileMenu);
            _menu.Items.Add(sampleMenu);
            _menu.Items.Add(helpMenu);
            MainMenuStrip = _menu;

            _toolbar = new ToolStrip();
            ToolStripButton btnOpen = new ToolStripButton("📂 打开");
            btnOpen.Click += delegate { OpenTextFile(); };
            ToolStripButton btnS1 = new ToolStripButton("示例1");
            btnS1.Click += delegate { LoadSample1(); };
            ToolStripButton btnS2 = new ToolStripButton("示例2");
            btnS2.Click += delegate { LoadSample2(); };
            ToolStripButton btnS3 = new ToolStripButton("示例3");
            btnS3.Click += delegate { LoadSample3(); };
            ToolStripSeparator sep1 = new ToolStripSeparator();
            ToolStripLabel lblMesh = new ToolStripLabel("网格尺寸:");
            _numMeshSize = new NumericUpDown();
            _numMeshSize.Width = 60;
            _numMeshSize.DecimalPlaces = 3;
            _numMeshSize.Increment = 0.05M;
            _numMeshSize.Minimum = 0.001M;
            _numMeshSize.Maximum = 100M;
            _numMeshSize.Value = 0.10M;
            ToolStripControlHost hostMesh = new ToolStripControlHost(_numMeshSize);
            _btnMesh = new Button();
            _btnMesh.Text = "🔺 剖分";
            _btnMesh.Width = 70;
            _btnMesh.Click += delegate { DoMesh(); };
            ToolStripControlHost hostBtnMesh = new ToolStripControlHost(_btnMesh);
            ToolStripSeparator sep2 = new ToolStripSeparator();
            _btnSolve = new Button();
            _btnSolve.Text = "🧮 求解";
            _btnSolve.Width = 70;
            _btnSolve.Font = new Font(_btnSolve.Font, FontStyle.Bold);
            _btnSolve.BackColor = Color.LightSteelBlue;
            _btnSolve.Click += delegate { DoSolve(); };
            ToolStripControlHost hostBtnSolve = new ToolStripControlHost(_btnSolve);
            ToolStripSeparator sep3 = new ToolStripSeparator();
            ToolStripLabel lblField = new ToolStripLabel("云图:");
            _cbField = new ComboBox();
            _cbField.Width = 90;
            _cbField.DropDownStyle = ComboBoxStyle.DropDownList;
            _cbField.Items.Add("Phi(势函数)");
            _cbField.Items.Add("Vmag(速度)");
            _cbField.SelectedIndex = 0;
            _cbField.SelectedIndexChanged += delegate { RenderCanvas(); };
            ToolStripControlHost hostField = new ToolStripControlHost(_cbField);
            _chkMesh = new CheckBox();
            _chkMesh.Text = "网格";
            _chkMesh.Width = 55;
            _chkMesh.Checked = true;
            _chkMesh.CheckedChanged += delegate { RenderCanvas(); };
            ToolStripControlHost hostChkMesh = new ToolStripControlHost(_chkMesh);
            _chkVectors = new CheckBox();
            _chkVectors.Text = "矢量";
            _chkVectors.Width = 55;
            _chkVectors.Checked = true;
            _chkVectors.CheckedChanged += delegate { RenderCanvas(); };
            ToolStripControlHost hostChkV = new ToolStripControlHost(_chkVectors);
            _chkLT6 = new CheckBox();
            _chkLT6.Text = "LT6";
            _chkLT6.Width = 50;
            _chkLT6.Checked = true;
            //_chkLT6.ToolTipText = "使用六节点二次三角形(LT6)，精度比CST提高一个量级";
            ToolStripControlHost hostLT6 = new ToolStripControlHost(_chkLT6);
            ToolStripLabel lblDt = new ToolStripLabel("Δt:");
            _numDt = new NumericUpDown();
            _numDt.Width = 60;
            _numDt.DecimalPlaces = 3;
            _numDt.Increment = 0.01M;
            _numDt.Minimum = 0.001M;
            _numDt.Maximum = 1000M;
            _numDt.Value = 0.05M;
            ToolStripControlHost hostDt = new ToolStripControlHost(_numDt);
            _btnTimeStep = new Button();
            _btnTimeStep.Text = "▶ 时间步";
            _btnTimeStep.Width = 78;
            _btnTimeStep.Click += delegate { DoTimeStep(); };
            ToolStripControlHost hostBtnStep = new ToolStripControlHost(_btnTimeStep);

            _toolbar.Items.Add(btnOpen);
            _toolbar.Items.Add(btnS1);
            _toolbar.Items.Add(btnS2);
            _toolbar.Items.Add(btnS3);
            _toolbar.Items.Add(sep1);
            _toolbar.Items.Add(lblMesh);
            _toolbar.Items.Add(hostMesh);
            _toolbar.Items.Add(hostBtnMesh);
            _toolbar.Items.Add(hostLT6);
            _toolbar.Items.Add(sep2);
            _toolbar.Items.Add(hostBtnSolve);
            _toolbar.Items.Add(hostBtnStep);
            _toolbar.Items.Add(lblDt);
            _toolbar.Items.Add(hostDt);
            _toolbar.Items.Add(sep3);
            _toolbar.Items.Add(lblField);
            _toolbar.Items.Add(hostField);
            _toolbar.Items.Add(hostChkMesh);
            _toolbar.Items.Add(hostChkV);

            _split = new SplitContainer();
            _split.Dock = DockStyle.Fill;
            _split.SplitterDistance = 30;
            _split.Orientation = Orientation.Vertical;

            _tree = new TreeView();
            _tree.Dock = DockStyle.Fill;
            _tree.Font = new Font("Consolas", 9f);
            _tree.AfterSelect += delegate { RenderCanvas(); };

            _canvasPanel = new Panel();
            _canvasPanel.Dock = DockStyle.Fill;
            _canvasPanel.AutoScroll = true;
            _canvasPanel.BackColor = Color.LightGray;

            SplitContainer rightSplit = new SplitContainer();
            rightSplit.Dock = DockStyle.Fill;
            rightSplit.Orientation = Orientation.Horizontal;
            rightSplit.SplitterDistance = 520;

            _canvas = new PictureBox();
            _canvas.Dock = DockStyle.Fill;
            _canvas.SizeMode = PictureBoxSizeMode.Zoom;
            _canvas.BackColor = Color.White;

            _logBox = new TextBox();
            _logBox.Dock = DockStyle.Fill;
            _logBox.Multiline = true;
            _logBox.ScrollBars = ScrollBars.Vertical;
            _logBox.Font = new Font("Consolas", 9f);
            _logBox.ReadOnly = true;
            _logBox.BackColor = Color.Black;
            _logBox.ForeColor = Color.LightGreen;

            rightSplit.Panel1.Controls.Add(_canvas);
            rightSplit.Panel2.Controls.Add(_logBox);

            _canvasPanel.Controls.Add(rightSplit);
            _split.Panel1.Controls.Add(_tree);
            _split.Panel2.Controls.Add(_canvasPanel);

            _status = new StatusStrip();
            _statusLabel = new ToolStripStatusLabel("就绪");
            _status.Items.Add(_statusLabel);

            Controls.Add(_split);
            Controls.Add(_toolbar);
            Controls.Add(_menu);
            Controls.Add(_status);
        }

        private void Log(string msg)
        {
            _logBox.AppendText("[" + DateTime.Now.ToString("HH:mm:ss") + "] " + msg + Environment.NewLine);
        }

        private void SetStatus(string s)
        {
            _statusLabel.Text = s;
        }

        private void RefreshTree()
        {
            _tree.Nodes.Clear();
            if (_model == null) return;
            TreeNode root = new TreeNode(_model.Title);
            TreeNode info = new TreeNode(string.Format("问题类型: {0}",
                _model.ProblemType == 0 ? "势流(速度势φ)" :
                _model.ProblemType == 1 ? "渗流(水头H)" : "导热(温度T)"));
            info.Nodes.Add("节点: " + _model.Nodes.Count);
            info.Nodes.Add("单元: " + _model.Elements.Count);
            info.Nodes.Add("材料: " + _model.Materials.Count);
            root.Nodes.Add(info);

            TreeNode matNode = new TreeNode("材料");
            for (int i = 0; i < _model.Materials.Count; i++)
            {
                Material m = _model.Materials[i];
                matNode.Nodes.Add(string.Format("#{0} {1}  k={2}  f={3}", m.Id, m.Name, m.K, m.Source));
            }
            root.Nodes.Add(matNode);

            TreeNode bcNode = new TreeNode("边界条件");
            int dCount = 0, nCount = 0;
            for (int i = 0; i < _model.Nodes.Count; i++)
            {
                if (_model.Nodes[i].BCType == 1) dCount++;
                else if (_model.Nodes[i].BCType == 2) nCount++;
            }
            bcNode.Nodes.Add("Dirichlet节点: " + dCount);
            bcNode.Nodes.Add("Neumann节点: " + nCount);
            root.Nodes.Add(bcNode);

            if (_lastResult != null && _lastResult.Success)
            {
                TreeNode res = new TreeNode("计算结果");
                res.Nodes.Add(string.Format("φ 范围: [{0:F4}, {1:F4}]", _lastResult.PhiMin, _lastResult.PhiMax));
                res.Nodes.Add(string.Format("|V| 范围: [{0:F4}, {1:F4}]", _lastResult.VmagMin, _lastResult.VmagMax));
                res.Nodes.Add("求解耗时: " + _lastResult.SolveTime.TotalMilliseconds.ToString("F1") + " ms");
                root.Nodes.Add(res);
            }
            _tree.Nodes.Add(root);
            root.ExpandAll();
        }

        private void RenderCanvas()
        {
            if (_model == null || _model.Elements.Count == 0 && _model.Nodes.Count == 0) return;
            try
            {
                string field = _cbField.SelectedIndex == 1 ? "Vmag" : "Phi";
                int w = Math.Max(400, _canvas.Width);
                int h = Math.Max(400, _canvas.Height);
                double vmin, vmax;
                List<LegendItem> legend;
                Bitmap bmp;
                if (_model.Elements.Count == 0)
                {
                    // 未剖分：只画节点和边界
                    bmp = new Bitmap(w, h);
                    using (Graphics g = Graphics.FromImage(bmp))
                    {
                        g.Clear(Color.White);
                        double xmin, ymin, xmax, ymax;
                        _model.ComputeBounds(out xmin, out ymin, out xmax, out ymax);
                        double margin = 40;
                        double dw = w - margin * 2;
                        double dh = h - margin * 2;
                        double scale = Math.Min(dw / Math.Max(1e-9, xmax - xmin), dh / Math.Max(1e-9, ymax - ymin));
                        double ox = margin + (dw - (xmax - xmin) * scale) / 2;
                        double oy = margin + (dh - (ymax - ymin) * scale) / 2;
                        // 连线（边界缓存）
                        List<Edge> edges = BoundaryEdgeCache.Get(_model);
                        if (edges != null)
                        {
                            using (Pen p = new Pen(Color.DarkBlue, 1.5f))
                            {
                                for (int i = 0; i < edges.Count; i++)
                                {
                                    Node n1 = _model.GetNode(edges[i].P1);
                                    Node n2 = _model.GetNode(edges[i].P2);
                                    if (n1 == null || n2 == null) continue;
                                    float x1 = (float)(ox + (n1.X - xmin) * scale);
                                    float y1 = (float)(oy + (ymax - n1.Y) * scale);
                                    float x2 = (float)(ox + (n2.X - xmin) * scale);
                                    float y2 = (float)(oy + (ymax - n2.Y) * scale);
                                    g.DrawLine(p, x1, y1, x2, y2);
                                }
                            }
                        }
                        for (int i = 0; i < _model.Nodes.Count; i++)
                        {
                            Node n = _model.Nodes[i];
                            float x = (float)(ox + (n.X - xmin) * scale);
                            float y = (float)(oy + (ymax - n.Y) * scale);
                            Brush br = Brushes.Black;
                            if (n.BCType == 1) br = Brushes.Red;
                            else if (n.BCType == 2) br = Brushes.Blue;
                            g.FillEllipse(br, x - 2, y - 2, 4, 4);
                        }
                        Font f = new Font("Arial", 12f, FontStyle.Bold);
                        g.DrawString(_model.Title + "（几何未剖分）", f, Brushes.Black, 20, 10);
                    }
                }
                else
                {
                    bmp = ResultExporter.RenderContour(_model, field, w, h, _chkMesh.Checked, _chkVectors.Checked, out vmin, out vmax, out legend);
                }
                Image old = _canvas.Image;
                _canvas.Image = bmp;
                if (old != null) old.Dispose();
            }
            catch (Exception ex)
            {
                Log("绘制失败: " + ex.Message);
            }
        }

        private void LoadSample1()
        {
            _model = SampleBuilder.CavityFlow(21, 21);
            _lastResult = null;
            _solver = null;
            Log("已加载示例1：方腔均匀势流 21×21 结构化网格");
            Log("左边界φ=0，右边界φ=1，上下壁面Neumann固壁。解析解：φ=x, u=(1,0)");
            RefreshTree();
            RenderCanvas();
            SetStatus("示例1已加载，节点=" + _model.Nodes.Count + " 单元=" + _model.Elements.Count);
        }

        private void LoadSample2()
        {
            string tmp = Path.Combine(Path.GetTempPath(), "fem_cylinder.dxf");
            SampleBuilder.WriteFlowAroundCylinderDxf(tmp);
            _model = MeshImporter.FromDxf(tmp);
            SampleBuilder.ApplyCylinderBC(_model);
            _solver = null;
            _numMeshSize.Value = 0.20M;
            Log("已加载示例2：绕圆柱势流DXF，网格尺寸建议0.2");
            Log("已自动设置：入口φ=0 出口φ=4 壁面/圆柱面Neumann固壁");
            Log("请点击 [🔺剖分] 生成三角网格，再点击 [🧮求解]");
            _lastResult = null;
            RefreshTree();
            RenderCanvas();
            SetStatus("示例2几何已加载，节点=" + _model.Nodes.Count);
        }

        private void LoadSample3()
        {
            string tmp = Path.Combine(Path.GetTempPath(), "fem_dam.txt");
            SampleBuilder.WriteDamSeepageInput(tmp);
            _model = MeshImporter.FromTextFile(tmp);
            _lastResult = null;
            _solver = null;
            Log("已加载示例3：达西渗流坝体 21×11 结构化网格");
            Log("上游φ=10 下游φ=2 k=0.001 m/s");
            RefreshTree();
            RenderCanvas();
            SetStatus("示例3已加载，节点=" + _model.Nodes.Count + " 单元=" + _model.Elements.Count);
        }

        private void DoMesh()
        {
            if (_model == null) { MessageBox.Show("请先加载模型"); return; }
            double h = (double)_numMeshSize.Value;
            try
            {
                Cursor = Cursors.WaitCursor;
                SetStatus("剖分中...");
                Log("开始网格剖分，目标网格尺寸 h=" + h.ToString("F3"));
                // 对只有边界/没有单元的几何，生成内部点再Delaunay
                if (_model.Elements.Count == 0)
                {
                    Triangulator.GenerateStructuredPoints(_model, h);
                    Triangulator.Delaunay(_model);
                    Triangulator.RemoveOutsideElements(_model);
                    // 重新计算单元几何（Delaunay结束时已算过，这里保险再算一次）
                    for (int ei = 0; ei < _model.Elements.Count; ei++)
                        Triangulator.ComputeElementGeom(_model.Elements[ei], _model);
                    // 剖分后所有边界节点（包括细分新增的中间点）按实际坐标重新打tag，
                    // 避免因边界细分遗漏入口/出口Dirichlet节点导致矩阵奇异
                    SampleBuilder.ApplyCylinderBC(_model);
                    const int TAG_INLET = 100;
                    const int TAG_OUTLET = 101;
                    const int TAG_WALL = 102;
                    int cntInlet = 0, cntOutlet = 0, cntWall = 0;
                    for (int i = 0; i < _model.Nodes.Count; i++)
                    {
                        Node n = _model.Nodes[i];
                        if (n.BCType == TAG_INLET) { n.BCType = 1; cntInlet++; }
                        else if (n.BCType == TAG_OUTLET) { n.BCType = 1; cntOutlet++; }
                        else if (n.BCType == TAG_WALL) { n.BCType = 0; cntWall++; }
                        // 内部点保持BCType=0（自然Neumann壁面）
                    }
                    Log("生成结构化内部点 + Delaunay三角剖分完成，单元数=" + _model.Elements.Count +
                        " 入口节点=" + cntInlet + " 出口节点=" + cntOutlet + " 壁面节点=" + cntWall);
                }
                else
                {
                    // 已有网格，仅重新计算几何信息
                    for (int i = 0; i < _model.Elements.Count; i++) Triangulator.ComputeElementGeom(_model.Elements[i], _model);
                    Log("已重新计算单元几何");
                }
                // 根据 LT6 复选框决定是否生成边中点（二次单元）
                bool useLT6 = _chkLT6.Checked;
                if (useLT6)
                {
                    Triangulator.AddMidEdgeNodes(_model);
                    for (int i = 0; i < _model.Elements.Count; i++) Triangulator.ComputeElementGeom(_model.Elements[i], _model);
                    Log("已为所有单元添加边中点，升级为LT6二次三角形（节点数=" + _model.Nodes.Count + "）");
                }
                else
                {
                    // 若之前已生成LT6中点则移除（简单处理：重置Order标记，不实际删除节点——再次剖分会重算）
                    for (int i = 0; i < _model.Elements.Count; i++) { _model.Elements[i].Order = 0; _model.Elements[i].N4 = 0; _model.Elements[i].N5 = 0; _model.Elements[i].N6 = 0; }
                }
                _lastResult = null;
                RefreshTree();
                RenderCanvas();
                SetStatus("剖分完成，节点=" + _model.Nodes.Count + " 单元=" + _model.Elements.Count);
            }
            catch (Exception ex)
            {
                Log("剖分失败: " + ex.Message);
                MessageBox.Show("剖分失败: " + ex.Message);
            }
            finally { Cursor = Cursors.Default; }
        }

        private void DoSolve()
        {
            if (_model == null) { MessageBox.Show("请先加载模型"); return; }
            if (_model.Elements.Count == 0)
            {
                if (MessageBox.Show("还没有网格，是否先自动剖分再求解？", "提示", MessageBoxButtons.YesNo) == DialogResult.Yes)
                {
                    DoMesh();
                }
                else return;
            }
            try
            {
                Cursor = Cursors.WaitCursor;
                SetStatus("求解中...");
                Log("开始有限元求解，节点=" + _model.Nodes.Count + " 单元=" + _model.Elements.Count);
                _model.ResetResults();
                FemResult r = FemSolver.Solve(_model);
                _lastResult = r;
                if (r.Success)
                {
                    Log("✓ 求解成功! 耗时=" + r.SolveTime.TotalMilliseconds.ToString("F1") + "ms");
                    Log(string.Format("  φ: [{0:F4}, {1:F4}]", r.PhiMin, r.PhiMax));
                    Log(string.Format("  |V|: [{0:F4}, {1:F4}]", r.VmagMin, r.VmagMax));
                }
                else
                {
                    Log("✗ 求解失败: " + r.Message);
                    MessageBox.Show(r.Message, "求解失败");
                }
                RefreshTree();
                RenderCanvas();
                SetStatus(r.Success ? "求解完成" : "求解失败");
            }
            catch (Exception ex)
            {
                Log("求解异常: " + ex.Message);
                MessageBox.Show("求解异常: " + ex.Message);
            }
            finally { Cursor = Cursors.Default; }
        }

        private void OpenTextFile()
        {
            OpenFileDialog dlg = new OpenFileDialog();
            dlg.Filter = "节点单元文件(*.txt;*.dat)|*.txt;*.dat|所有文件(*.*)|*.*";
            if (dlg.ShowDialog() == DialogResult.OK)
            {
                try
                {
                    _model = MeshImporter.FromTextFile(dlg.FileName);
                    _lastResult = null;
                    Log("已导入文本文件: " + dlg.FileName);
                    for (int i = 0; i < _model.Elements.Count; i++) Triangulator.ComputeElementGeom(_model.Elements[i], _model);
                    RefreshTree();
                    RenderCanvas();
                    SetStatus("导入完成，节点=" + _model.Nodes.Count + " 单元=" + _model.Elements.Count);
                }
                catch (Exception ex)
                {
                    MessageBox.Show("导入失败: " + ex.Message);
                }
            }
        }

        private void OpenDxfFile()
        {
            OpenFileDialog dlg = new OpenFileDialog();
            dlg.Filter = "DXF文件(*.dxf)|*.dxf|所有文件(*.*)|*.*";
            if (dlg.ShowDialog() == DialogResult.OK)
            {
                try
                {
                    _model = MeshImporter.FromDxf(dlg.FileName);
                    _lastResult = null;
                    Log("已导入DXF几何: " + dlg.FileName + " 节点=" + _model.Nodes.Count);
                    Log("请设置网格尺寸后点 [🔺剖分] 生成三角网格，再设置边界条件并求解");
                    Log("提示：DXF导入后所有边界节点默认标记为壁面(Neumann自然条件)，请在文本中指定入口/出口Dirichlet节点");
                    RefreshTree();
                    RenderCanvas();
                    SetStatus("DXF几何已加载，节点=" + _model.Nodes.Count);
                }
                catch (Exception ex)
                {
                    MessageBox.Show("DXF导入失败: " + ex.Message);
                }
            }
        }

        private void ExportCsv()
        {
            if (_model == null) return;
            using (FolderBrowserDialog dlg = new FolderBrowserDialog())
            {
                dlg.Description = "选择导出目录";
                if (dlg.ShowDialog() == DialogResult.OK)
                {
                    string nodePath = Path.Combine(dlg.SelectedPath, "nodes_result.csv");
                    string elemPath = Path.Combine(dlg.SelectedPath, "elements_result.csv");
                    ResultExporter.ExportNodeCsv(_model, nodePath);
                    ResultExporter.ExportElementCsv(_model, elemPath);
                    Log("已导出CSV:\n  " + nodePath + "\n  " + elemPath);
                    MessageBox.Show("已导出节点/单元结果CSV到:\n" + dlg.SelectedPath);
                }
            }
        }

        private void ExportTec()
        {
            if (_model == null) return;
            SaveFileDialog dlg = new SaveFileDialog();
            dlg.Filter = "Tecplot文件(*.dat)|*.dat";
            dlg.FileName = "result.dat";
            if (dlg.ShowDialog() == DialogResult.OK)
            {
                ResultExporter.ExportTecplot(_model, dlg.FileName);
                Log("已导出Tecplot: " + dlg.FileName);
                MessageBox.Show("导出完成");
            }
        }

        private void ExportPng()
        {
            if (_model == null || _canvas.Image == null) return;
            SaveFileDialog dlg = new SaveFileDialog();
            dlg.Filter = "PNG图片(*.png)|*.png";
            dlg.FileName = "contour.png";
            if (dlg.ShowDialog() == DialogResult.OK)
            {
                _canvas.Image.Save(dlg.FileName, System.Drawing.Imaging.ImageFormat.Png);
                Log("已导出云图PNG: " + dlg.FileName);
                MessageBox.Show("云图已保存");
            }
        }

        private void DoTimeStep()
        {
            if (_model == null) { MessageBox.Show("请先加载模型"); return; }
            if (_model.Elements.Count == 0)
            {
                MessageBox.Show("请先剖分网格再做时间推进");
                return;
            }
            try
            {
                Cursor = Cursors.WaitCursor;
                double dt = (double)_numDt.Value;
                if (_solver == null) _solver = new FemSolver();
                FemResult r = _solver.AdvanceTime(_model, dt);
                _lastResult = r;
                Log("时间步推进 Δt=" + dt.ToString("F3") + "s - " + r.Message);
                RefreshTree();
                RenderCanvas();
                SetStatus(r.Success ? "时间步完成" : "时间步失败");
                if (!r.Success) MessageBox.Show(r.Message, "求解失败");
            }
            catch (Exception ex)
            {
                Log("时间步异常: " + ex.Message);
                MessageBox.Show("时间步异常: " + ex.Message);
            }
            finally { Cursor = Cursors.Default; }
        }

        private void ExportTextReport()
        {
            if (_model == null || _lastResult == null) { MessageBox.Show("请先求解后再导出报告"); return; }
            using (FolderBrowserDialog dlg = new FolderBrowserDialog())
            {
                dlg.Description = "选择文本报告输出目录";
                if (dlg.ShowDialog() == DialogResult.OK)
                {
                    FemSolver.ExportTextResults(_model, _lastResult, dlg.SelectedPath);
                    Log("已导出文本报告+CSV+VTK到: " + dlg.SelectedPath);
                    MessageBox.Show("已导出:" + dlg.SelectedPath);
                }
            }
        }
    }
}
