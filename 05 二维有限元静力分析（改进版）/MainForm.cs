using System;
using System.Drawing;
using System.IO;
using System.Threading;
using System.Windows.Forms;
using FEM2D.Elements;
using FEM2D.IO;
using FEM2D.Mesh;
using FEM2D.Rendering;
using FEM2D.Solvers;

namespace FEM2D
{
    public sealed class MainForm : Form
    {
        private FEMesh _mesh;
        private MenuStrip _menu;
        private ToolStripMenuItem _miFile,_miImport,_miExport,_miExit;
        private ToolStripMenuItem _miMesh,_miGenCST,_miGenLT6,_miGenQ4;
        private ToolStripMenuItem _miSolve,_miStructCG,_miStructGMRES,_miStructBiCGSTAB,_miStructAMG,_miNS,_miNSUnsteady;
        private ToolStripMenuItem _miHelp,_miAbout;
        private SplitContainer _split;
        private PropertyGrid _propGrid;
        private ComboBox _fieldBox;
        private PictureBox _picture;
        private TextBox _logBox;
        private FEMParams _params;

        public MainForm()
        {
            Text="二维有限元计算 (FEM2D) - Visual Studio 2026 / .NET 10.0";
            Width=1400;Height=900;StartPosition=FormStartPosition.CenterScreen;
            BuildMenu();BuildUI();
            _params=new FEMParams();_propGrid.SelectedObject=_params;
        }

        private void BuildMenu()
        {
            _menu=new MenuStrip();
            _miFile=new ToolStripMenuItem("文件(&F)");
            _miImport=new ToolStripMenuItem("导入网格 INP...");
            _miImport.Click+=delegate(object s,EventArgs e){DoImport();};
            _miExport=new ToolStripMenuItem("导出结果...");
            _miExport.Click+=delegate(object s,EventArgs e){DoExport();};
            _miExit=new ToolStripMenuItem("退出");_miExit.Click+=delegate(object s,EventArgs e){Close();};
            _miFile.DropDownItems.Add(_miImport);_miFile.DropDownItems.Add(_miExport);
            _miFile.DropDownItems.Add(new ToolStripSeparator());_miFile.DropDownItems.Add(_miExit);

            _miMesh=new ToolStripMenuItem("网格(&M)");
            _miGenCST=new ToolStripMenuItem("自动生成：悬臂梁 CST3");
            _miGenCST.Click+=delegate(object s,EventArgs e){GenerateMesh(ElementType.CST3);};
            _miGenLT6=new ToolStripMenuItem("自动生成：悬臂梁 LT6");
            _miGenLT6.Click+=delegate(object s,EventArgs e){GenerateMesh(ElementType.LT6);};
            _miGenQ4=new ToolStripMenuItem("自动生成：顶盖驱动方腔 Q4");
            _miGenQ4.Click+=delegate(object s,EventArgs e){GenerateMesh(ElementType.Q4);};
            _miMesh.DropDownItems.Add(_miGenCST);_miMesh.DropDownItems.Add(_miGenLT6);_miMesh.DropDownItems.Add(_miGenQ4);

            _miSolve=new ToolStripMenuItem("求解(&S)");
            _miStructCG=new ToolStripMenuItem("结构静力求解 (CG + Jacobi)");
            _miStructCG.Click+=delegate(object s,EventArgs e){RunStructural(KrylovSolverKind.CG,PreconditionerKind.Jacobi);};
            _miStructGMRES=new ToolStripMenuItem("结构静力求解 (GMRES + Jacobi)");
            _miStructGMRES.Click+=delegate(object s,EventArgs e){RunStructural(KrylovSolverKind.GMRES,PreconditionerKind.Jacobi);};
            _miStructBiCGSTAB=new ToolStripMenuItem("结构静力求解 (BiCGSTAB + Jacobi)");
            _miStructBiCGSTAB.Click+=delegate(object s,EventArgs e){RunStructural(KrylovSolverKind.BiCGSTAB,PreconditionerKind.Jacobi);};
            _miStructAMG=new ToolStripMenuItem("结构静力求解 (CG + AMG 预条件)");
            _miStructAMG.Click+=delegate(object s,EventArgs e){RunStructural(KrylovSolverKind.CG,PreconditionerKind.AMG);};
            _miNS=new ToolStripMenuItem("NS 定常 SIMPLE");
            _miNS.Click+=delegate(object s,EventArgs e){RunNS(false);};
            _miNSUnsteady=new ToolStripMenuItem("NS 非稳态 (向后欧拉+SIMPLE)");
            _miNSUnsteady.Click+=delegate(object s,EventArgs e){RunNS(true);};
            _miSolve.DropDownItems.Add(_miStructCG);_miSolve.DropDownItems.Add(_miStructGMRES);
            _miSolve.DropDownItems.Add(_miStructBiCGSTAB);_miSolve.DropDownItems.Add(_miStructAMG);
            _miSolve.DropDownItems.Add(new ToolStripSeparator());_miSolve.DropDownItems.Add(_miNS);_miSolve.DropDownItems.Add(_miNSUnsteady);

            _miHelp=new ToolStripMenuItem("帮助(&H)");
            _miAbout=new ToolStripMenuItem("关于");
            _miAbout.Click+=delegate(object s,EventArgs e)
            {
                MessageBox.Show("二维有限元计算程序 (FEM2D)\n版本 2.0\n\n支持单元：CST3 / LT6(二次三角) / Q4\n求解器：CG / GMRES / BiCGSTAB\n预条件：Jacobi / AMG (Smoothed Aggregation)\nSpMV：多核 SIMD 向量化\n流体：Navier-Stokes SIMPLE（定常/非稳态，Picard对流、向后欧拉）\n\n输出：CSV / VTK(Legacy) / PNG 云图\n\nVisual Studio 2026 · .NET 10.0 WinForms","关于 FEM2D",MessageBoxButtons.OK,MessageBoxIcon.Information);
            };
            _miHelp.DropDownItems.Add(_miAbout);

            _menu.Items.Add(_miFile);_menu.Items.Add(_miMesh);_menu.Items.Add(_miSolve);_menu.Items.Add(_miHelp);
            MainMenuStrip=_menu;Controls.Add(_menu);
        }

        private void BuildUI()
        {
            _split=new SplitContainer{Dock=DockStyle.Fill,Orientation=Orientation.Vertical,SplitterDistance=310};
            Panel left=new Panel{Dock=DockStyle.Fill,Padding=new Padding(6)};
            _propGrid=new PropertyGrid{Dock=DockStyle.Fill,HelpVisible=false,PropertySort=PropertySort.Categorized};
            Panel top=new Panel{Dock=DockStyle.Top,Height=55,Padding=new Padding(4)};
            Label lbl=new Label{Text="显示字段：",Dock=DockStyle.Top,Height=20};
            _fieldBox=new ComboBox{Dock=DockStyle.Top,DropDownStyle=ComboBoxStyle.DropDownList};
            _fieldBox.Items.AddRange(new object[]{"网格 Mesh","位移 Displacement","变形 Deformed","应力 σxx","应力 σyy","剪应力 τxy","Von Mises 应力","速度 Velocity","压力 Pressure"});
            _fieldBox.SelectedIndex=0;
            _fieldBox.SelectedIndexChanged+=delegate(object s,EventArgs e){RefreshPlot();};
            top.Controls.Add(_fieldBox);top.Controls.Add(lbl);
            left.Controls.Add(_propGrid);left.Controls.Add(top);
            _split.Panel1.Controls.Add(left);

            Panel right=new Panel{Dock=DockStyle.Fill};
            SplitContainer rs=new SplitContainer{Dock=DockStyle.Fill,Orientation=Orientation.Horizontal,SplitterDistance=640};
            _picture=new PictureBox{Dock=DockStyle.Fill,BackColor=Color.White,SizeMode=PictureBoxSizeMode.Zoom};
            _logBox=new TextBox{Dock=DockStyle.Fill,Multiline=true,ScrollBars=ScrollBars.Vertical,ReadOnly=true,Font=new Font("Consolas",9),BackColor=Color.Black,ForeColor=Color.Lime};
            rs.Panel1.Controls.Add(_picture);rs.Panel2.Controls.Add(_logBox);
            right.Controls.Add(rs);
            _split.Panel2.Controls.Add(right);
            Controls.Add(_split);
        }

        private void Log(string msg)
        {
            if(IsDisposed)return;
            if(InvokeRequired){try{Invoke(new Action<string>(Log),msg);}catch{}return;}
            _logBox.AppendText(msg+Environment.NewLine);
        }

        private void DoImport()
        {
            using(OpenFileDialog dlg=new OpenFileDialog())
            {
                dlg.Filter="Abaqus INP (*.inp)|*.inp|所有文件 (*.*)|*.*";
                string ed=Path.Combine(AppDomain.CurrentDomain.BaseDirectory,"Examples");
                if(Directory.Exists(ed))dlg.InitialDirectory=ed;
                if(dlg.ShowDialog(this)==DialogResult.OK)
                {
                    try{_mesh=MeshReader.Read(dlg.FileName);
                        Log("已导入网格："+Path.GetFileName(dlg.FileName));
                        Log(string.Format("  节点 {0}, 单元 {1}, 类型 {2}",_mesh.NumNodes,_mesh.NumElements,_mesh.ElementType));
                        RefreshPlot();}
                    catch(Exception ex){Log("导入失败："+ex.Message);}
                }
            }
        }

        private void GenerateMesh(ElementType type)
        {
            try{
                if(type==ElementType.Q4)
                {_mesh=StructuredMeshGenerator.GenerateLidDrivenCavity(type,_params.Length,_params.Nx,_params.LidVelocity,
                        _params.FluidDensity,_params.Viscosity,_params.Thickness);
                    Log(string.Format("已生成顶盖驱动方腔 Q4：{0}x{1}, L={2}, U0={3}",_params.Nx,_params.Nx,_params.Length,_params.LidVelocity));}
                else
                {_mesh=StructuredMeshGenerator.GenerateCantilever(type,_params.Length,_params.Height,_params.Nx,_params.Ny,_params.LoadY,
                        _params.Young,_params.Poisson,_params.Density,_params.Viscosity,_params.Thickness);
                    Log(string.Format("已生成悬臂梁 ({0})：{1}x{2}, L={3}, H={4}",type,_params.Nx,_params.Ny,_params.Length,_params.Height));}
                Log(string.Format("  节点 {0}, 单元 {1}",_mesh.NumNodes,_mesh.NumElements));
                RefreshPlot();}
            catch(Exception ex){Log("网格生成失败："+ex.Message);}
        }

        private void RunStructural(KrylovSolverKind ks, PreconditionerKind pk)
        {
            if(_mesh==null){MessageBox.Show("请先创建或导入网格。","提示",MessageBoxButtons.OK,MessageBoxIcon.Exclamation);return;}
            SetButtonsEnabled(false);
            string solverName=ks.ToString()+" + "+pk.ToString();
            Log("==== 结构静力分析开始（"+solverName+"） ====");
            ThreadPool.QueueUserWorkItem(delegate(object s)
            {
                try{StructuralSolver solver=new StructuralSolver();
                    solver.SolverKind=ks;
                    solver.PreconditionerKind=pk;
                    solver.Tolerance=1.0e-8;
                    solver.MaxIterations=5000;
                    solver.OnLog+=Log;
                    solver.OnIteration+=delegate(int it,double r){if(it%50==0)Log(string.Format("  迭代 {0}: 残差 {1:E3}",it,r));};
                    DateTime t0=DateTime.Now;
                    solver.Solve(_mesh);
                    TimeSpan dt=DateTime.Now-t0;
                    Log("求解完成，耗时 "+dt.TotalSeconds.ToString("F2")+" 秒");
                    BeginInvoke(new Action(RefreshPlot));}
                catch(Exception ex){Log("求解失败："+ex);}
                finally{BeginInvoke(new Action<bool>(SetButtonsEnabled),true);}
            });
        }

        private void RunNS(bool unsteady)
        {
            if(_mesh==null){MessageBox.Show("请先创建或导入网格。","提示",MessageBoxButtons.OK,MessageBoxIcon.Exclamation);return;}
            SetButtonsEnabled(false);
            Log("==== Navier-Stokes "+(unsteady?"非稳态":"定常")+" SIMPLE 开始 ====");
            ThreadPool.QueueUserWorkItem(delegate(object s)
            {
                try{NavierStokesSolver ns=new NavierStokesSolver();
                    ns.Density=_params.FluidDensity;ns.Viscosity=_params.Viscosity;ns.TimeStep=_params.TimeStep;
                    ns.Unsteady=unsteady;ns.TimeSteps=unsteady?_params.NumTimeSteps:1;
                    ns.MaxSIMPLEIter=_params.MaxSIMPLE;ns.ConvergenceTol=_params.NSTol;
                    ns.AlphaU=_params.AlphaU;ns.AlphaP=_params.AlphaP;
                    ns.OnLog+=Log;
                    DateTime t0=DateTime.Now;
                    ns.Solve(_mesh);
                    TimeSpan dt=DateTime.Now-t0;
                    Log("NS 求解完成，耗时 "+dt.TotalSeconds.ToString("F2")+" 秒");
                    BeginInvoke(new Action(RefreshPlot));}
                catch(Exception ex){Log("NS 求解失败："+ex);}
                finally{BeginInvoke(new Action<bool>(SetButtonsEnabled),true);}
            });
        }

        private void DoExport()
        {
            if(_mesh==null){MessageBox.Show("暂无结果。","提示",MessageBoxButtons.OK,MessageBoxIcon.Exclamation);return;}
            using(FolderBrowserDialog dlg=new FolderBrowserDialog())
            {
                dlg.Description="选择导出目录（CSV/VTK/PNG）";
                if(dlg.ShowDialog(this)==DialogResult.OK)
                {
                    string field=_fieldBox.SelectedItem!=null?_fieldBox.SelectedItem.ToString():"Mesh";
                    string tag=_mesh.ElementType.ToString();
                    try{ResultWriter.WriteAll(_mesh,dlg.SelectedPath,"FEM2D_"+tag,field);
                        Log("结果已导出到："+dlg.SelectedPath);
                        MessageBox.Show("导出成功！\n\n已生成：\n - 节点 CSV\n - 单元 CSV\n - VTK 非结构网格 (ParaView可读)\n - PNG 彩色云图","导出完成",MessageBoxButtons.OK,MessageBoxIcon.Information);}
                    catch(Exception ex){Log("导出失败："+ex.Message);}
                }
            }
        }

        private void RefreshPlot()
        {
            if(_mesh==null)return;
            try{string field=_fieldBox.SelectedItem!=null?_fieldBox.SelectedItem.ToString():"Mesh";
                bool showDeformed=field.ToLower().Contains("deform")||field.Contains("变形");
                int w=Math.Max(800,_picture.ClientSize.Width);
                int h=Math.Max(500,_picture.ClientSize.Height);
                Bitmap bmp=ContourRenderer.Render(_mesh,w,h,field,showDeformed);
                Image old=_picture.Image;_picture.Image=bmp;if(old!=null)old.Dispose();}
            catch(Exception ex){Log("绘图失败："+ex.Message);}
        }

        private void SetButtonsEnabled(bool enabled)
        {
            if(IsDisposed)return;
            if(InvokeRequired){try{Invoke(new Action<bool>(SetButtonsEnabled),enabled);}catch{}return;}
            _miImport.Enabled=enabled;_miExport.Enabled=enabled;
            _miGenCST.Enabled=enabled;_miGenLT6.Enabled=enabled;_miGenQ4.Enabled=enabled;
            _miStructCG.Enabled=enabled;_miStructGMRES.Enabled=enabled;_miStructBiCGSTAB.Enabled=enabled;_miStructAMG.Enabled=enabled;
            _miNS.Enabled=enabled;_miNSUnsteady.Enabled=enabled;
            Cursor=enabled?Cursors.Default:Cursors.WaitCursor;
        }

        protected override void OnResize(EventArgs e)
        {base.OnResize(e);if(_picture!=null&&_mesh!=null)RefreshPlot();}
    }

    public sealed class FEMParams
    {
        [System.ComponentModel.Category("几何")]
        [System.ComponentModel.Description("梁长度/方腔边长 L")]
        public double Length{get;set;}=10.0;
        [System.ComponentModel.Category("几何")]
        [System.ComponentModel.Description("梁高度 H")]
        public double Height{get;set;}=2.0;
        [System.ComponentModel.Category("几何")]public int Nx{get;set;}=20;
        [System.ComponentModel.Category("几何")]public int Ny{get;set;}=4;

        [System.ComponentModel.Category("结构材料")][System.ComponentModel.Description("杨氏模量 E [Pa]")]public double Young{get;set;}=2.1e11;
        [System.ComponentModel.Category("结构材料")]public double Poisson{get;set;}=0.3;
        [System.ComponentModel.Category("结构材料")]public double Thickness{get;set;}=1.0;
        [System.ComponentModel.Category("结构材料")]public double Density{get;set;}=7850.0;
        [System.ComponentModel.Category("荷载")]public double LoadY{get;set;}=-1000.0;

        [System.ComponentModel.Category("流体")]public double FluidDensity{get;set;}=1.0;
        [System.ComponentModel.Category("流体")][System.ComponentModel.Description("动力粘度 μ (Re=ρUL/μ)")]public double Viscosity{get;set;}=0.01;
        [System.ComponentModel.Category("流体")]public double LidVelocity{get;set;}=1.0;
        [System.ComponentModel.Category("流体")]public double TimeStep{get;set;}=0.01;
        [System.ComponentModel.Category("流体")]public int NumTimeSteps{get;set;}=20;

        [System.ComponentModel.Category("SIMPLE")]public int MaxSIMPLE{get;set;}=100;
        [System.ComponentModel.Category("SIMPLE")]public double NSTol{get;set;}=1e-4;
        [System.ComponentModel.Category("SIMPLE")]public double AlphaU{get;set;}=0.7;
        [System.ComponentModel.Category("SIMPLE")]public double AlphaP{get;set;}=0.3;
    }
}
