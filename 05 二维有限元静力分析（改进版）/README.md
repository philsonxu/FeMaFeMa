# FEM2D — 二维有限元计算桌面程序

基于 C# / .NET 10.0 / WinForm 的二维有限元分析平台，支持**结构刚性分析**与**不可压缩 Navier-Stokes 流动**两类问题。

## ✨ 功能特性

| 模块 | 能力 |
|------|------|
| **单元类型** | 3 节点常应变三角单元 (CST3)、6 节点二次三角单元 (LT6/P2)、4 节点双线性四边形 (Q4) |
| **网格** | 类 Abaqus INP 文件导入；矩形区域结构化网格自动剖分（悬臂梁 / 顶盖驱动方腔预设） |
| **线性代数** | CSR 稀疏矩阵存储 + 线程安全 COO 装配器；多线程（ThreadPool + Parallel）单元组装；Jacobi 预条件重启 GMRES(m)；CG 共轭梯度法 |
| **结构分析** | 平面应变线弹性；penalty 法 Dirichlet 边界；集中力荷载；节点应力平均外推 (σxx/σyy/σxy/Von Mises) |
| **流体分析** | 不可压缩 Navier-Stokes 方程；SIMPLE 算法（动量预测→压力泊松→速度/压力修正）；Picard 线性化对流项；向后欧拉非稳态时间项；欠松弛；人工压缩压力稳定化 |
| **可视化** | GDI+ 实时云图（Jet 色阶）；网格线框、位移、变形前后对比（虚线轮廓）、应力、速度、压力；右侧色条图例 + 数值标注 |
| **结果输出** | CSV 节点/单元数据表；VTK Legacy 非结构化网格（ParaView 可读）；1600×1200 PNG 彩色云图 |
| **界面** | WinForm 桌面程序；PropertyGrid 参数面板；字段切换下拉框；后台 Task 求解不阻塞 UI；实时彩色日志 |

## 📁 目录结构

```
FEM2D/
├── FEM2D.slnx                        # Visual Studio 2026 解决方案
├── README.md
└── FEM2D/
    ├── FEM2D.csproj                  # .NET 10.0 WinForm 项目
    ├── Program.cs                    # 程序入口
    ├── MainForm.cs                   # WinForm 主界面
    ├── Core/
    │   ├── DenseVector.cs            # 稠密向量 + BLAS1（并行 Dot/AddScaled）
    │   ├── SparseMatrixCSR.cs        # CSR 稀疏矩阵 + 线程安全 COO 装配器（支持矩形）
    │   ├── GMRESSolver.cs            # GMRES(m) + Givens 旋转 + Jacobi 左预条件
    │   ├── CGSolver.cs               # 共轭梯度法
    │   └── ParallelTool.cs           # 多线程并行封装
    ├── Elements/
    │   ├── FiniteElement.cs          # 单元抽象基类（刚度/质量/对流/粘性/压力梯度/应力）
    │   ├── CST3Element.cs            # 3 节点三角单元（解析 B 矩阵）
    │   ├── LT6Element.cs             # 6 节点二次三角单元（3 点面积坐标高斯积分）
    │   └── Q4Element.cs              # 4 节点四边形单元（2×2 高斯积分，等参）
    ├── Mesh/
    │   ├── FEMesh.cs                 # 网格容器
    │   ├── MeshReader.cs             # INP 文件导入器
    │   └── StructuredMeshGenerator.cs  # 结构化网格自动剖分
    ├── Solvers/
    │   ├── StructuralSolver.cs       # 结构静力求解器（CG/GMRES）
    │   └── NavierStokesSolver.cs     # NS 方程 SIMPLE 求解器（定常/非稳态）
    ├── IO/
    │   └── ResultWriter.cs           # CSV / VTK 输出
    ├── Rendering/
    │   └── ContourRenderer.cs        # GDI+ 云图/色条绘制
    └── Examples/
        ├── cantilever_cst.inp        # 悬臂梁 CST3 算例
        └── lid_driven_q4.inp         # 顶盖驱动方腔 Q4 算例
```

## 🔨 编译与运行

### 环境要求
- **Visual Studio 2026**（或 .NET SDK 10.0）
- Windows 操作系统（WinForm）

### 命令行编译
```bash
cd FEM2D
dotnet build -c Release
dotnet run --project FEM2D/FEM2D.csproj
```

### Visual Studio 打开
直接双击 `FEM2D.slnx`，按 **F5** 运行即可。

## 📖 使用指南

### 1. 准备网格
- **自动剖分**：在 [网格] 菜单中选择单元类型（CST3 / LT6 / Q4）和场景（悬臂梁 / 顶盖驱动方腔），程序会按左侧参数面板的几何尺寸和网格密度自动生成。
- **文件导入**：[文件] → [导入网格 INP...]，打开 `.inp` 文件（示例见 `Examples/` 目录）。

### 2. 参数设置
左侧 PropertyGrid 面板可调整：
- **几何**：Length（长度）、Height（高度）、Nx / Ny（网格分段数）
- **材料**：YoungModulus（杨氏模量）、PoissonRatio（泊松比）、Density（密度）、Viscosity（动力粘度 μ）、Thickness（厚度）
- **荷载**：TipLoad（悬臂梁端力）、LidVelocity（顶盖速度）
- **求解**：DeltaT（时间步长）、TimeSteps（时间步数）、Unsteady（是否非稳态）、DeformScale（变形放大系数）

### 3. 执行求解
- **结构静力 (CG / GMRES)**：线弹性平面应变问题，适用于悬臂梁等固体力学算例。
- **NS 定常 (SIMPLE)**：不可压缩定常流动，适用于顶盖驱动方腔等经典算例。
- **NS 非稳态 (SIMPLE+Δt)**：向后欧拉时间离散的瞬态流动。

### 4. 查看结果
- 右侧"显示字段"下拉框可切换 9 种物理量：网格、位移大小、变形图、σxx、σyy、σxy、Von Mises、速度大小、压力。
- 云图右侧显示 Jet 色阶图例（最小/中间/最大值）。
- 下方黑色日志窗口实时输出求解过程（迭代步、残差）。

### 5. 导出结果
[文件] → [导出结果...]，选择目录后将输出：
- `nodes.csv` — 节点坐标 + 所有物理量（逗号分隔，可直接用 Excel/Pandas 打开）
- `elements.csv` — 单元拓扑（类型 + 节点编号）
- `result.vtk` — VTK Legacy 非结构化网格（可直接拖入 ParaView 做三维后处理）
- `contour_*.png` — 各字段 1600×1200 高清彩色云图

## 📄 INP 文件格式

采用类 Abaqus 关键字格式，节点编号从 1 开始：

```
*MATERIAL, E=2.1e11, NU=0.3, RHO=7800, MU=1e-3, T=0.01
*NODE
id, x, y
...
*ELEMENT, TYPE=CST3|LT6|Q4
id, n1, n2, ...
...
*BOUNDARY          ! 结构位移边界（dof: 1=x, 2=y）
node, dof, value
...
*CLOAD             ! 集中力
node, dof, value
...
*VELOCITY          ! 流体速度边界
node, dof, value
...
```

示例见 `Examples/cantilever_cst.inp` 与 `Examples/lid_driven_q4.inp`。

## 🧮 数值方法概要

### 结构静力
- **平面应变本构** D 矩阵：CST 采用解析 B 矩阵；LT6 与 Q4 分别采用 3 点 Hammer 与 2×2 Gauss 积分。
- **并行组装**：每个单元在 ThreadPool 线程中独立计算刚度矩阵，通过线程安全的 COO 装配器合并三元组，最后转 CSR。
- **Dirichlet 边界**：penalty 法（penalty = max|diag| × 10¹²）。
- **求解器**：CG（默认，对称正定）或 GMRES(m=50) + Jacobi 对角预条件。
- **应力恢复**：单元应力在节点处求值后，按相邻单元数做节点平均。

### Navier-Stokes (SIMPLE)
- **动量方程**：(M/Δt + μK + C(u*)) u = -Gp + (M/Δt)u_n
  - M：一致质量矩阵
  - K：粘性扩散矩阵（μ∇²）
  - C(u*)：Picard 线性化对流项（ρ(u·∇)u，用当前速度线性化）
  - G：压力梯度矩阵（-∫N_p ∇·N_u dΩ）
- **压力泊松方程**：Ap p' = Gᵀu^，其中 Ap ≈ Gᵀ diag(A)⁻¹ G + εMp（ε=0.1 人工压缩稳定化，避免同阶单元的压力棋盘）
- **速度修正**：u ← u − α_u · diag(A)⁻¹ G p'（α_u = 0.7）
- **压力修正**：p ← p + α_p p'（α_p = 0.3）
- **非稳态**：向后欧拉（一阶隐式），Δt 可设；定常时 Δt→∞（10³⁰）近似稳态。
- **压力参考点**：固定节点 0 的压力为 0，避免压力奇异性。

## 📊 内置算例

1. **悬臂梁**（CST3/LT6/Q4）：长 10m，高 1m，厚 0.01m，E=210GPa，ν=0.3；左端固支，右端中点受竖向力 -1000N。可观察到位移（向下弯曲）与 σxx 弯曲应力分布。
2. **顶盖驱动方腔**（Q4/CST3）：1×1 方腔，ρ=1，μ=0.01，顶盖 u=1 m/s（Re=100）。可观察到中心主涡与角涡结构。

## 🔧 扩展方向

- 增加轴对称、平面应力选项
- 增加 8/9 节点二次四边形 Q8/Q9
- 增加 Neumann 边界面力积分
- 增加材料非线性（弹塑性）与几何非线性
- 集成 PETSc / Intel MKL / cuSPARSE 等外部高性能求解器
- 增加截面曲线绘制、积分点查询等后处理

## ⚠️ 注意事项

- **代码规范**：所有 C# 源码严格遵循"不使用 `var`"约定，所有局部变量均显式声明类型。
- **运行环境**：面向 .NET 10.0，需 Visual Studio 2026 或最新 .NET SDK 10。
- **SIMPLE 收敛性**：粗网格或高雷诺数下可能需要更多外迭代或更小的欠松弛因子；建议网格数 ≥ 20×20 以获得合理流场。
- **时间步长**：非稳态计算应满足 CFL 条件：Δt < Δx / U_max。

---
© FEM2D 项目 — 仅供学习与研究使用。
