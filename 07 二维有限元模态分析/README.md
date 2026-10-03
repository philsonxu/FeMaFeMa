# ModalFEM2D — 二维有限元桌面程序（模态分析 / 静力 / Navier-Stokes）

ModalFEM2D 是一套从零构建的 Windows 桌面二维有限元计算程序，使用 **C# 13 + .NET 10.0 + WinForms**。核心特色是 **二维结构模态分析**（子空间逆迭代 + Rayleigh-Ritz + AMG 预条件），同时提供线弹性静力求解与不可压缩 Navier-Stokes SIMPLE 求解，支持 CST3/LT6/Q4 三种二维单元、多核 SIMD SpMV、VTK/CSV/PNG 输出。

## 功能特性

| 模块 | 实现要点 |
|---|---|
| **单元类型** | 3 节点常应变三角（CST3）、6 节点二次三角（LT6，3 点 Hammer 积分）、4 节点双线性四边形（Q4，2×2 高斯）|
| **网格** | 类 Abaqus INP 文件导入；内置悬臂梁 / 简支梁 / 顶盖驱动方腔自动剖分 |
| **线性代数** | CSR 稀疏矩阵 + COO 并行装配；`System.Numerics.Vector<T>` SIMD SpMV（自动派发 SSE/AVX/AVX-512）+ 多线程按行并行 |
| **Krylov 求解器** | CG（对称正定）、GMRES(m)（重启 + Givens 正交化）、BiCGSTAB（非对称） |
| **预条件器** | Jacobi 对角；**AMG 代数多重网格（Smoothed Aggregation + Galerkin RAP + V-cycle + Gauss-Jordan 粗层直接求逆）** |
| **结构静力** | 平面应变（Plane Strain）本构；多线程单元组装；penalty 法施加 Dirichlet；CG/GMRES/BiCGSTAB + Jacobi/AMG；单元应力恢复 + 节点平均 |
| **模态分析 ⭐** | 子空间逆迭代 + Rayleigh-Ritz 投影（子空间维数 `m=2p+4`）；M-正交 Gram-Schmidt；内迭代 CG + AMG；小投影广义特征问题用 Cholesky→标准对称阵→Jacobi 旋转求特征对；输出前 N 阶固有频率（Hz）与归一化振型 |
| **Navier–Stokes** | 不可压缩牛顿流体；SIMPLE 分块迭代（动量预测→压力泊松→速度/压力修正）；Picard 线性化对流项；向后欧拉非稳态时间项；αu/αp 欠松弛；等阶插值 + 压力 Laplacian 稳定化 |
| **结果输出** | CSV（节点/单元表）、VTK Legacy 非结构化网格（ParaView 可读，含各阶模态）、PNG Jet 色阶云图（含色条、变形前后对比、虚线未变形轮廓）、`modes.txt` 频率报告 |
| **WinForm 界面** | 菜单 + PropertyGrid 参数面板 + 字段下拉（网格/位移/变形/σxx/σyy/τxy/Von Mises/速度/压力/模态振型）+ 模态阶次选择 + PictureBox 自适应云图 + 彩色实时日志；后台 Task 求解不阻塞 UI |

## 目录结构

```
ModalFEM2D/
├── ModalFEM2D.sln              # VS 经典 sln（兼容旧版 VS）
├── ModalFEM2D.slnx             # Visual Studio 2026 解决方案
├── README.md
└── ModalFEM2D/
    ├── ModalFEM2D.csproj       # net10.0-windows, UseWindowsForms, AllowUnsafeBlocks
    ├── Program.cs              # WinForm 入口
    ├── MainForm.cs             # 主界面（菜单/参数/云图/日志）
    ├── Core/                   # 稀疏线性代数
    │   ├── DenseVector.cs      # 稠密向量 + BLAS1（并行 Dot/Axpy）
    │   ├── SparseMatrixCSR.cs  # CSR + SIMD SpMV + COO 装配器
    │   ├── IPreconditioner.cs  # 预条件接口 + Jacobi
    │   ├── AMGPreconditioner.cs # Smoothed Aggregation AMG
    │   └── KrylovSolvers.cs    # CG / GMRES(m) / BiCGSTAB
    ├── Elements/
    │   ├── FiniteElement.cs    # 单元抽象基类（刚度/质量/粘性/对流/压力梯度/应力/面积）
    │   ├── CST3Element.cs      # 3 节点三角
    │   ├── LT6Element.cs       # 6 节点二次三角
    │   ├── Q4Element.cs        # 4 节点四边形
    │   └── ElementFactory.cs   # 单元工厂
    ├── Mesh/
    │   ├── FEMesh.cs           # 网格数据容器（节点/单元/BC/结果场）
    │   ├── MeshReader.cs       # INP 读入
    │   └── StructuredMeshGenerator.cs  # 矩形自动剖分
    ├── Solvers/
    │   ├── StructuralSolver.cs # 线弹性静力
    │   ├── ModalSolver.cs      # 子空间逆迭代模态分析 ⭐
    │   └── NavierStokesSolver.cs  # SIMPLE + 向后欧拉
    ├── IO/ResultWriter.cs      # CSV / VTK / PNG
    ├── Rendering/ContourRenderer.cs  # GDI+ Jet 云图 + 色条 + 变形
    └── Examples/
        ├── cantilever_cst.inp  # 悬臂梁 CST3 算例
        └── lid_driven_q4.inp   # 顶盖驱动方腔 Q4 算例
```

## 编译与运行

### 前置条件
- **.NET 10 SDK**（Visual Studio 2026 自带，或从 https://dotnet.microsoft.com/ 下载）
- Windows 10/11（使用 WinForms）

### 命令行
```bash
cd ModalFEM2D
dotnet build -c Release
dotnet run --project ModalFEM2D/ModalFEM2D.csproj -c Release
```

### Visual Studio 2026
打开 `ModalFEM2D.slnx`（推荐）或 `ModalFEM2D.sln`，按 F5 运行。

## 使用指南

1. **建立网格**：
   - 菜单 `[网格] → 悬臂梁 CST3 / 悬臂梁 LT6 / 简支梁 CST3 (模态) / 顶盖驱动方腔 Q4` 一键自动生成；
   - 或 `[文件] → 导入 INP...` 打开 `Examples/` 下的现成算例；
   - 左侧 PropertyGrid 可修改几何尺寸、网格密度、材料参数、荷载、模态阶数、SIMPLE 松弛因子等。
2. **求解**：
   - `[求解] → 结构静力 (CG/GMRES/BiCGSTAB+Jacobi, CG+AMG)` 做线性弹性静力；
   - `[求解] → 模态分析 (子空间迭代 + AMG)` ⭐ 做模态分析；
   - `[求解] → NS 定常/NS 非稳态` 运行 SIMPLE。
3. **后处理**：
   - 左下"显示"下拉切换 10 种物理量（网格/位移/变形/σxx/σyy/τxy/Von Mises/速度/压力/模态振型）；
   - "模态阶次"下拉在提取完模态后可逐阶切换；
   - 右下彩色日志窗口实时显示迭代收敛过程。
4. **导出**：`[文件] → 导出结果...` 选择目录，自动写出 `nodes.csv`、`elements.csv`、`modes.txt`、`result.vtk`、`contour.png`（1600×1200）。

## 输入文件格式（INP）

类 Abaqus 的关键字文本（节点、单元 1-based，导入时自动转 0-based）：

```
*NODE
id, x, y
...
*ELEMENT, TYPE=CST3 | LT6 | Q4 | S4
eid, n1, n2, n3 [, n4, n5, n6]
...
*BOUNDARY
nid, first_dof, last_dof, value     # 1=ux, 2=uy
*CLOAD
nid, dir, value
*MATERIAL
E, nu, rho, mu, thickness
```

## 数值方法要点

- **模态分析**：采用经典子空间逆迭代（Bathe）——初始子空间由随机向量张成；每步解 `K Y = M X`（CG+AMG 内迭代精度 1e-10）；在 Y 张成的子空间内通过 Rayleigh-Ritz 投影求解小广义特征问题，利用 Mp 的 Cholesky 分解化为对称特征问题，再用经典 Jacobi 旋转法求特征对；M-正交 Gram-Schmidt 保证数值稳定。Dirichlet 约束通过刚度对角加 `1e18` 惩罚项施加，提取时自动跳过惩罚模式（λ≈1e18）。

- **AMG 预条件**：Smoothed Aggregation——按对角权值选种子点，通过强耦合阈值（θ=0.08）BFS 聚合；tentative 延拓算子（piecewise constant）经一次阻尼 Jacobi 平滑；粗层矩阵 `Ac = P^T A P` 通过行式 CSR×CSR 三积计算；V(2,2)-cycle（前后各 2 步加权 Jacobi 光滑，ω=0.6）；最粗层 n<80 切换为 Gauss-Jordan 稠密直接求逆。

- **SIMPLE 算法**：
  1. 装配动量算子 `A_u = (ν/ρ)∇² + (u·∇) + (1/Δt)M_lump`（隐式欠松弛：A_u/αu + (1/αu−1)·diag(A_u)）；
  2. 解 `A_u u* = rhs`（BiCGSTAB+Jacobi）；
  3. 压力 Poisson 方程 `−D·A_u⁻¹·G p' = −D·u*`，加等阶插值稳定化项 `ε∫∇p·∇q dΩ`；
  4. 速度修正 `u = u* − αu·A_u⁻¹·G p'`，压力修正 `p += αp·p'`；
  5. 迭代至 `max|div(u)| < tol`。

- **多核 SIMD SpMV**：按行分块分配到 `Environment.ProcessorCount` 个线程；每个线程内部非零元按 `Vector<double>.Count` 分块，`xv` 通过 `stackalloc Span<double>` gather 后构造向量，`Vector.Dot` 一次完成一个向量块的乘加。运行时自动适配 CPU 的 SIMD 宽度（SSE2/AVX/AVX2/AVX-512/SVE）。

## 代码规范

- 所有 C# 源码**不使用 `var` 关键字**（全部显式类型声明）；
- 所有模块按职责分层（`Core`/`Elements`/`Mesh`/`Solvers`/`IO`/`Rendering`）；
- 单元通过抽象基类 `FiniteElement` 统一接口，新增单元类型只需继承并实现 `ShapeAndDerivs`/`GetIntegrationPoints` 即可；
- 并行装配使用线程局部 COO 装配器 + WaitHandle 同步，无全局锁竞争。

## 注意事项

- 模态分析推荐在 1000+ 自由度网格上使用 `CG+AMG` 组合；对于粗网格，AMG 可能仅构造少量级别，但依然比 Jacobi 收敛快；
- NS SIMPLE 为教学实现，顶盖驱动方腔在 Re≤400（μ=0.01, L=1, U=1）范围内可稳定收敛；高 Re 需要更精细的网格与更小的 αp；
- 等阶速度-压力插值使用了 Laplacian 稳定化（系数 0.3），若出现棋盘格压力振荡可适当增大该系数（见 `NavierStokesSolver.cs` 中 `stabCoef`）。

## 许可

教学/科研演示代码，作者不承担任何使用责任；可自由修改扩展。
