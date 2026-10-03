# 二维多物理场耦合有限元 (MultiPhysics FEM 2D)

一套基于 **C# / WinForms / .NET 10.0** 的桌面型二维有限元计算程序，覆盖结构、热传导、Navier-Stokes 流动，以及热-结构 / Boussinesq 流-热两种多物理场耦合分析。

## 主要特性

### 物理场与耦合

| 模块 | 求解器 | 说明 |
|------|--------|------|
| 结构静力 | `StructuralSolver` | 平面应力线弹性，多线程单元装配，penalty Dirichlet，集中力/均布荷载，单元应力→节点面积加权平均 |
| 稳态/瞬态热传导 | `ThermalSolver` | Fourier 导热模型；集中热容 + 向后欧拉时间积分；支持 Dirichlet 温度 / Neumann 热流 / 体热源 |
| **★ 热-结构单向耦合** | `ThermalStructuralSolver` | T → ε_th = αΔT → F_th = B^T D ε_th → 结构位移/应力 |
| 等温 NS | `NavierStokesSolver` | SIMPLE 算法：动量预测 → 压力 Poisson → 速度/压力修正；Picard 线性化对流 + α_u/α_p 欠松弛 |
| NS 非稳态 | 同上 | 向后欧拉 + 集中质量时间项 |
| **★ Boussinesq 自然对流** | NavierStokesSolver + Thermal | 三场（u,p,T）固定点耦合，浮力体积力 f_y = −ρgβ(T−T_ref) 注入动量方程，温度含对流-扩散 |

### 有限单元
| 单元 | 节点 | 积分 | 支持物理场 |
|------|------|------|------------|
| **CST3**（3 节点常应变三角） | 3 | 1 点形心 | 结构刚/质、热导/热容、NS 粘性/对流/压力梯度/压力质量、热载荷 |
| **LT6**（6 节点二次三角） | 6 | 3 点 Hammer 面积坐标 | 同上 |
| **Q4**（4 节点双线性四边形） | 4 | 2×2 Gauss-Legendre | 同上 |

每个单元严格实现 `FiniteElement` 抽象基类的全部接口（刚度/质量/导/容/粘性/对流/压力梯度/压力质量/能量对流/应力/热载荷/形函数）。

### 稀疏线性代数
- **CSR 稀疏矩阵**：压缩行存储，`MultiplyAdd` 支持 α·A·x + β·y 标准 BLAS-2 接口；
- **COO 三元组装配器**：线程安全（内部 lock），支持每线程批量合并；
- **多核 SIMD SpMV**：`System.Numerics.Vector<double>` 自动适配 SSE2 / AVX / AVX2 / AVX-512 等向量宽度，行块 ThreadPool 并行，小矩阵自动串行；
- **Krylov 子空间求解器**：
  - CG（对称正定）
  - GMRES(m) + Givens 旋转重正交
  - BiCGSTAB（非对称对流系统，内存 O(n)）
- **预条件子**：
  - Jacobi 对角（快速、通用）
  - **AMG 代数多重网格**（Smoothed Aggregation）：贪心种子 + 强耦合阈值 θ 聚合 → damped-Jacobi 平滑 prolongation → Galerkin 粗化 `A_c = P^T A P`（自研多线程 CSR 三积）→ V(2,2)-cycle；最粗层 Gauss-Jordan 稠密求逆

### 网格
- **INP 导入**：类 Abaqus 文本格式（关键字 `*NODE / *ELEMENT,TYPE=CPS3|CPS4|CPS6 / *BOUNDARY`）；节点 1-based 自动转换；
- **自动剖分**：`StructuredMeshGenerator` 一键生成 5 种预设算例：
  - 悬臂梁 CST3（端载弯曲）
  - 悬臂梁 Q4
  - **★ 热-结构耦合梁** CST3（顶热/底冷/左固支）
  - 顶盖驱动方腔 Q4（NS 验证算例）
  - **★ 自然对流方腔 Q4**（Boussinesq：左热/右冷/4 壁无滑）

### 结果输出
- **CSV**：节点表（x, y, ux, uy, vx, vy, p, T, σxx, σyy, τxy, VM, ε_p）+ 单元表
- **VTK Legacy**：非结构化网格，ParaView 可读，含位移/速度/压力/温度/应力/等效应力/塑性应变向量与标量
- **PNG 云图**：Jet 色阶（蓝-青-绿-黄-红）+ 右侧数值色条 + 黑色细线网格 + 灰色虚线未变形轮廓；自动适配画布；支持 12 种物理量实时切换

### WinForm 界面
- 顶栏菜单：文件 / 网格 / 求解 / 帮助
- 左侧：`PropertyGrid` 参数面板（几何、结构材料、热、流体、时间/迭代参数）
- Krylov / 预条件 / 显示字段下拉
- 右上：自适应云图预览
- 右下：黑色实时日志窗口（绿色文字）
- 所有求解在后台 ThreadPool 运行，UI 不阻塞

## 工程文件

- `MultiPhysicsFEM2D.slnx` / `MultiPhysicsFEM2D.sln`：Visual Studio 2026 解决方案
- `MultiPhysicsFEM2D/MultiPhysicsFEM2D.csproj`：`net10.0-windows`，`UseWindowsForms=true`，`AllowUnsafeBlocks=true`（SIMD SpMV stackalloc）
- **零 `var` 关键字**：全项目显式类型声明，可通过 `grep -rn "\bvar\b" *.cs` 验证
- 示例输入文件：`Examples/cantilever_cst.inp`（悬臂梁）、`Examples/lid_driven_q4.inp`（顶盖驱动方腔）

## 目录结构

```
MultiPhysicsFEM2D/
├── MultiPhysicsFEM2D.sln / .slnx / README.md
└── MultiPhysicsFEM2D/
    ├── MultiPhysicsFEM2D.csproj
    ├── Program.cs · MainForm.cs
    ├── Core/        DenseVector · SparseMatrixCSR(CooAssembler) ·
    │               Preconditioners (Identity/Jacobi) · AMGPreconditioner ·
    │               KrylovSolvers (CG/GMRES/BiCGSTAB)
    ├── Elements/    FiniteElement (抽象基类) · CST3Element · LT6Element ·
    │               Q4Element · ElementFactory
    ├── Mesh/        FEMesh · MeshReader (INP 导入) · StructuredMeshGenerator
    ├── Solvers/     StructuralSolver · ThermalSolver ·
    │               ThermalStructuralSolver ⭐ · NavierStokesSolver (Boussinesq) ⭐
    ├── IO/ResultWriter.cs (CSV/VTK/PNG)
    ├── Rendering/ContourRenderer.cs (GDI+ Jet 色阶云图)
    └── Examples/    cantilever_cst.inp · lid_driven_q4.inp
```

## 编译与运行

### Visual Studio 2026
1. 打开 `MultiPhysicsFEM2D.slnx`（或 `.sln`）
2. 选择 **Release | Any CPU**，按 **Ctrl+F5** 启动

### 命令行（需要 .NET 10 SDK）
```bash
cd MultiPhysicsFEM2D
dotnet build -c Release
dotnet run --project MultiPhysicsFEM2D/MultiPhysicsFEM2D.csproj -c Release
```

## 使用指南

### 快速体验

| 步骤 | 算例 | 推荐参数 |
|------|------|----------|
| 1 | 热-结构耦合梁 | 默认参数；`[网格]→★ 热-结构耦合梁 CST3` → `[求解]→★ 热-结构耦合(单向)`；查看"温度 T"→"等效应力 Von Mises" |
| 2 | 顶盖驱动 NS | 设 `Nx=Ny=30`, `Mu=0.01`, `UTop=1.0`, `AlphaU=0.7`, `AlphaP=0.3`, `NSOuter=500` → `[网格]→顶盖驱动方腔 Q4` → `[求解]→NS 定常 SIMPLE` |
| 3 | Boussinesq 自然对流 | 设 `Nx=Ny=30`, `Mu=0.001`, `K=0.001`, `RhoFluid=1`, `Beta=0.00034`, `TLeft=1,TRight=0,Gravity=10` → `[网格]→★ 自然对流方腔 Q4` → `[求解]→★ 流-热 Boussinesq 自然对流`；查看"温度 T"（S 形等温线） |
| 4 | 静力分析 | `[网格]→悬臂梁 CST3` → `[求解]→结构静力`（可在左侧切换 Krylov=CG/BiCGSTAB，预条件=Jacobi/AMG） |

### 参数说明
- **Length / Height / Nx / Ny**：矩形几何与网格分辨率
- **E / Nu / Density / Thickness**：结构材料（弹性模量、泊松比、密度、厚度）
- **EndLoad**：悬臂梁右端中点 Y 向集中力
- **K / Cp / Alpha / T0 / Thot / Tcold**：热导率、比热、热膨胀系数、参考温度、壁面温度
- **RhoFluid / MuFluid / Beta**：流体密度、动力粘度、体积膨胀系数（Boussinesq）
- **UTop**：顶盖速度（NS）
- **AlphaU / AlphaP**：SIMPLE 欠松弛因子（典型 0.3~0.8 / 0.1~0.3）
- **DT / TimeSteps**：非稳态时间步长与步数

### Krylov × 预条件推荐组合

| 规模 | 推荐组合 |
|------|----------|
| 小（< 5k DOF） | 任意 Krylov + Jacobi |
| 结构大网格（> 10k DOF） | **CG + AMG**（椭圆算子 h-无关收敛） |
| NS 动量（非对称对流） | **BiCGSTAB + Jacobi**（内存低、稳定） |

## 输入文件格式（INP）

```
*NODE
<nodeId>, <x>, <y>                  # 节点编号 1-based
*ELEMENT, TYPE=CPS3 | CPS4 | CPS6
<elemId>, <n1>, <n2>, ...           # 单元连接（1-based）
*BOUNDARY
<nodeId>, <BCType>, <value>         # FixX/FixY/FixBoth/TractionX/TractionY
                                    # Vx/Vy/P/Temp/HeatFlux
```

支持的边界条件关键字：`FixX, FixY, FixBoth, ForceX(Fx), ForceY(Fy), Vx, Vy, P, Temp(T), HeatFlux(Q), Encastre`。

## 数值方法要点

1. **SIMPLE 压力修正**：
   - 动量方程使用 H（粘性 + 对流 + 瞬态质量）组装，以 Jacobi 迭代近似求 u\*
   - 压力 Poisson 用 Laplacian + 压力质量稳定化（避免棋盘格）
   - 速度修正采用逐点近似 u' = −H_ii^{−1} ∇p'，α_u/α_p 欠松弛
2. **AMG 粗化**：强耦合阈值 θ·√|A_ii·A_jj|，damped-Jacobi(ω=0.7) 平滑 prolongation；每层前后各 2 步光滑（V(2,2)）；最粗层 ≤64 个未知量时 Gauss-Jordan 直接求逆
3. **SIMD SpMV**：`Vector<double>` 自动向量化，对非零元按 `vlen` 分块 gather-x，`Vector.Dot(vv, xv)` 完成点积
4. **多线程并行**：单元矩阵组装按行块分配给 ThreadPool（WaitHandle 同步），每线程局部数组累加后合并到 CooAssembler，避免逐 entry lock

## 扩展方向

- 单元层：8 节点 Serendipity 四边形（Q8）、9 节点 Lagrangian 四边形（Q9）、4 节点四面体（扩展三维）
- 本构：J2 von Mises 弹塑性（径向返回映射 + 一致切线）、超弹性 Neo-Hookean
- 求解：弧长法（Riks/Crisfield）、瞬态动力学 Newmark/HHT、FGMRES、LU/ILU(0) 预条件
- 耦合：FSI 双向耦合（移动网格/ALE）、热-电耦合（Joule 热）、压电耦合
- 可视化：变形动画导出（GIF/MP4）、切片、矢量箭头 glyph

## 许可与注意事项

本项目为教学与快速原型验证用途，结果正确性取决于网格质量与参数；SIMPLE 压力方程采用简化处理（Jacobi 内解 + 稳定化），高 Re 流动可能需要更精确的 Rhie-Chow 插值或更严格的压力离散；学术/工程使用建议与商用软件（Abaqus/ANSYS/OpenFOAM）做算例校核。
