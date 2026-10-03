# FEM2D — 二维有限元计算桌面程序 (C# / .NET 10.0 / WinForms)

本项目是一套从零搭建的二维有限元计算桌面程序，包含结构刚性分析与不可压缩流体 Navier–Stokes 求解两大模块，使用 Visual Studio 2026 的工程格式（`.csproj` + `.slnx`），面向 .NET 10.0 Windows 桌面（WinForms）运行。

## ✨ 主要特性

- **线弹性结构静力分析**（平面应变假设），提供悬臂梁等计算实例，含位移、应力、Von Mises 等结果
- **三种单元**：
  - **CST3**：三节点常应变三角单元
  - **LT6**：六节点二次三角单元（面积坐标 + 3 点 Hammer 积分）
  - **Q4**：四节点双线性等参四边形单元（2×2 高斯积分）
- **网格来源**：
  - 从类 Abaqus INP 文本文件导入（`*NODE`、`*ELEMENT`、`*BOUNDARY`、`*CLOAD`、`*DSLOAD`）
  - 矩形区域结构化自动剖分（三种单元一键生成），自动施加悬臂梁/顶盖驱动方腔边界
- **高性能数值代数**：
  - **CSR 格式稀疏矩阵**（COO 装配器 + 多线程 `Parallel.For` 单元组装、多线程 SpMV）
  - **共轭梯度法 (CG)** 求解对称正定结构刚度方程
  - **重启 GMRES(m)**（Givens 旋转 + 雅可比对角预条件）求解非对称系统
- **Navier–Stokes 求解器**：
  - 不可压缩 N–S 方程，**SIMPLE 算法**（动量预测 → 压力泊松 → 速度/压力修正）
  - 对流项 **Picard 线性化**（`u·∇u`）
  - 向后欧拉 **非稳态时间项** `∂u/∂t`（支持多时间步推进）
  - 速度/压力欠松弛（αu、αp），雅可比预条件 GMRES 内迭代
- **结果输出**：
  - 数据文件：节点/单元 CSV、VTK Legacy 非结构化网格（ParaView 可直接打开）
  - 彩色云图：Jet 色阶渲染网格、位移、σxx/σyy/σxy、Von Mises、速度大小、压力、变形图，右侧含图例色条
  - 位图 PNG 1600×1200 导出
- **WinForm 图形界面**：左侧 PropertyGrid 参数面板、右侧 PictureBox 云图 + 实时文本日志、顶部菜单（文件/网格/求解/帮助）
- **编码规范**：全部 C# 源代码显式类型声明，**不使用 `var` 关键字**

## 📁 目录结构

```
FEM2D_project/
├── FEM2D.slnx                    # Visual Studio 2026 解决方案
├── README.md                     # 本文档
└── FEM2D/
    ├── FEM2D.csproj              # .NET 10.0 WinForms 项目
    ├── Program.cs                # 程序入口
    ├── MainForm.cs               # WinForm 主界面
    ├── Core/                     # 线性代数与求解器核心
    │   ├── DenseVector.cs        # 稠密向量 + BLAS-1
    │   ├── SparseMatrixCSR.cs    # CSR 稀疏矩阵 + COO 并行装配器
    │   ├── Solvers.cs            # GMRES(m) + CG 迭代求解器
    │   └── ParallelTool.cs       # 多线程工具
    ├── Elements/                 # 单元库
    │   ├── FiniteElement.cs      # 单元抽象基类
    │   ├── CST3Element.cs        # 三节点三角单元
    │   ├── LT6Element.cs         # 六节点二次三角单元
    │   └── Q4Element.cs          # 四节点四边形单元
    ├── Mesh/                     # 网格数据与 I/O
    │   ├── FEMesh.cs             # 网格数据容器
    │   ├── MeshReader.cs         # INP 文件导入
    │   └── StructuredMeshGenerator.cs  # 矩形区域自动剖分（含 NS 算例）
    ├── Solvers/                  # 物理求解器
    │   ├── StructuralSolver.cs   # 线弹性结构静力（并行组装+CG/GMRES）
    │   └── NavierStokesSolver.cs # NS SIMPLE（对流+粘性+压力耦合+非稳态）
    ├── IO/
    │   └── ResultWriter.cs       # CSV + VTK + PNG 输出
    ├── Rendering/
    │   └── ContourRenderer.cs    # GDI+ Jet 色阶云图（带图例）
    └── Examples/                 # 计算实例
        ├── cantilever_cst.inp    # 悬臂梁 CST3 算例
        └── lid_driven_q4.inp     # 顶盖驱动方腔 Q4 算例
```

## 🔨 编译与运行

### 环境要求
- **Visual Studio 2026**（或 .NET SDK 10.0.100+）
- Windows 操作系统（WinForms 依赖 `Microsoft.WindowsDesktop.App`）
- .NET Desktop Runtime 10.0

### 用 Visual Studio 2026
1. 打开 `FEM2D.slnx`
2. 选择 `Any CPU / Debug` 或 `Release`，按 `F5` 即可编译运行

### 命令行
```powershell
cd FEM2D_project
dotnet build -c Release
dotnet run --project FEM2D/FEM2D.csproj -c Release
```

## 🖱️ 使用指南

1. **生成或导入网格**
   - 点击菜单 **[网格]** → 选择单元类型（CST3 / LT6 / Q4）自动生成矩形网格；
   - 或点击 **[文件]→[导入网格 INP...]** 导入 `Examples/cantilever_cst.inp` 等文件。
2. **在左侧参数面板调整参数**
   - 几何（LengthX/Y）、网格密度（Nx/Ny）
   - 材料（E、ν、Thickness）、荷载（FixLeft、TopPressure、EndLoad）
   - 流体参数（Rho、Mu、LidVelocity）、SIMPLE 松弛因子（αu、αp）、非稳态设置（Dt、NTimeSteps）
3. **运行求解**
   - **[求解]→[结构静力分析 (CG)]**：悬臂梁等算例；可选 GMRES 非对称求解器
   - **[求解]→[NS 定常 (SIMPLE)]**：顶盖驱动方腔
   - **[求解]→[NS 非稳态]**：带时间项的瞬态 NS
4. **查看结果**：右下角文本框实时显示迭代日志；右上云图可通过 **显示字段** 下拉框切换：
   - 结构：位移大小 / 变形图 / σxx / σyy / σxy / Von Mises
   - 流体：速度大小 / 压力
5. **导出结果**：**[文件]→[导出结果...]**，选择目录后将写出：
   - `nodes_displacement.csv` / `nodes_ns.csv`、`elements_stress.csv`
   - `structural.vtk` / `ns.vtk`（ParaView 可视化）
   - `vonmises.png` / `deformed.png` / `velocity.png` / `pressure.png`（1600×1200 云图）

## 📄 输入文件格式（类 Abaqus INP）

```
*NODE
<id>, x, y
...
*ELEMENT, TYPE=CPS3|CPS6|CPS4
<id>, n1, n2, ...
*BOUNDARY
<node>, <dof(1=ux/u,2=uy/v,3=p)>, <value>
*CLOAD
<node>, <dof>, <value>
*DSLOAD
<elem>, <face>, <pressure>
```

编号均从 **1** 开始；行以 `**` 开头视为注释。

## 🧮 数值方法摘要

| 模块 | 方法 |
|------|------|
| 单元积分 | CST 1 点；LT6 3 点 Hammer；Q4 2×2 高斯 |
| 全局矩阵存储 | CSR（Compressed Sparse Row） |
| 装配方式 | `Parallel.For` 单元级别并行 + COO 三元组合并去重 |
| Dirichlet 边界 | 大罚函数法（penalty = 1e20 × 对角元） |
| 结构求解 | CG + 雅可比预条件（对称正定）；或 GMRES(m=80) |
| NS 对流 | Picard 定点线性化 `C(u_old) u_new` |
| NS 时间 | 向后 Euler，集总质量矩阵 M/Δt |
| 压力-速度耦合 | SIMPLE：动量预测 → 压力泊松（标量 Laplacian 近似 + 对角 A⁻¹）→ 修正 |
| 欠松弛 | `A/αu`，右端加 `(1/αu - 1) diag(A) u_prev`；压力 `p += αp p'` |

## 📊 计算实例

1. **悬臂梁（CST3）**：`Examples/cantilever_cst.inp`，L=10 m, H=1 m, E=210 GPa, ν=0.3，左端固支，右端 Fy=-10 kN。可验证端部挠度近似解析解 `δ ≈ FL³/(3EI)`。
2. **顶盖驱动方腔（Q4）**：1×1 方形，u_top=1，Re=ρUL/μ=100，SIMPLE 迭代 200 步内收敛，速度云图呈现典型主涡结构。

## 🔧 扩展方向

- 支持平面应力、轴对称、平面梁/桁架单元
- 添加 Petrov-Galerkin (SUPG/PSPG) 稳定化以提升高 Re 数精度
- 支持更多求解器：BiCGSTAB、AMG 预条件、多核 SIMD SpMV
- 引入 Newton-Raphson 求解几何/材料非线性
- 支持瞬态结构动力学（Newmark / HHT-α）

## 📜 许可证

本代码仅供学习与研究使用，欢迎在此基础上扩展。
