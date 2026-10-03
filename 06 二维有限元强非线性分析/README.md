# 二维强非线性有限元桌面程序 NonlinearFEM2D

基于 C# / WinForms / .NET 10 开发的二维平面有限元桌面程序，核心特色为 **J2 von Mises 弹塑性强非线性分析**（Newton-Raphson 增量迭代 + 径向返回映射），同时保留线弹性静力与不可压缩 Navier–Stokes 求解能力。

---

## 主要功能

### 结构分析
- **强非线性结构分析**：J2 von Mises 弹塑性本构（各向同性线性硬化、Simo–Hughes 向后 Euler 径向返回、一致弹塑性切线）；
- **Newton–Raphson 增量求解**：支持力加载与位移加载两种控制模式，多线程并行装配一致切线 Kt 与内力 Rint；
- **几何非线性**：总拉格朗日（TL）初应力刚度 Kσ，可与 J2 塑性叠加实现材料+几何双重非线性；
- **线弹性静力**：平面应变线弹性刚度，支持 CG / GMRES / BiCGSTAB + Jacobi / AMG 求解。

### 单元库
| 单元 | 类型 | 积分点 | 非线性 |
|------|------|--------|--------|
| CST3 | 3 节点常应变三角 | 1 点（中心） | ✅ 完整 J2 一致切线 |
| LT6  | 6 节点二次三角（P2） | 3 点 Hammer | 线性切线（后处理兼容） |
| Q4   | 4 节点双线性四边形 | 2×2 Gauss | 线性切线（后处理兼容） |

> 非线性求解器以 CST3 为主要算例单元（每个单元一个 Gauss 点，径向返回精确），LT6/Q4 的线性切线可用于线弹性分析。

### 稀疏线性代数
- **CSR 稀疏存储** + **COO 多线程装配器**（每线程局部 builder，合并阶段加锁）；
- **多核 SIMD SpMV**：使用 `System.Numerics.Vector<double>` 自动派发 SSE / AVX / AVX-512，按行块 `Parallel.For` 并行；
- **Krylov 三兄弟**：
  - CG（对称正定）
  - GMRES(m)（Givens 旋转、重启 m=80，通用非对称）
  - BiCGSTAB（稳定双共轭梯度，NS/动量方程首选）
- **AMG 代数多重网格预条件**（Smoothed Aggregation）：
  - 强耦合阈值 BFS 聚合；
  - 阻尼 Jacobi 平滑 prolongation；
  - Galerkin 粗化 `Ac = RAP`（自研多线程 CSR×CSR 三积）；
  - V(2,2)-cycle，最粗层 Gauss-Jordan 稠密求逆；
  - 对椭圆算子收敛率近似 h-无关。

### 流体 Navier–Stokes
- **SIMPLE 算法**：
  1. 动量预测（隐式 u*，含 Picard 线性化对流、粘性、质量项）；
  2. 压力 Poisson 方程；
  3. 速度/压力修正（欠松弛 αu、αp）；
- **Picard 线性化对流项** `(u*·∇)u`；
- **向后 Euler 非稳态时间项** `(ρ/Δt)M(u-uⁿ)`；
- 支持定常（时间项关闭）与非稳态两种模式。

### 网格 I/O
- **INP 文件导入**：类 Abaqus 文本格式（*NODE/*ELEMENT/*BOUNDARY/*MATERIAL）；
- **自动剖分**：一键生成悬臂梁（CST3/Q4）、简支梁、顶盖驱动方腔（Q4）矩形网格，自动施加固支/荷载/壁面边界条件。

### 结果输出
- **数据文件**：
  - `nodes.csv`（节点坐标+位移+速度+压力+应力+等效塑性应变）；
  - `elements.csv`（单元连接）；
  - `result.vtk`（VTK Legacy 非结构化网格，ParaView 直接打开，含位移/速度/压力/VM/塑性应变张量）；
- **彩色图**：Jet 色阶 PNG 云图 + 右侧数值色条 + 变形前后虚线对比（黑色虚线为未变形轮廓）；
- 支持字段：网格、位移大小、变形图、σxx、σyy、τxy、Von Mises、等效塑性应变 εp、速度大小、压力。

---

## 目录结构

```
NonlinearFEM2D/
├── NonlinearFEM2D.sln / .slnx        # Visual Studio 2026 解决方案
├── README.md                          # 本文档
└── NonlinearFEM2D/                    # 项目目录
    ├── NonlinearFEM2D.csproj          # net10.0-windows 项目
    ├── Program.cs                     # 入口
    ├── MainForm.cs                    # WinForm 主界面
    ├── Core/
    │   ├── DenseVector.cs             # 稠密向量（BLAS1, 并行 Dot）
    │   ├── SparseMatrixCSR.cs         # CSR 稀疏矩阵 + COO 装配器 + SIMD 并行 SpMV
    │   ├── ILinearSolvers.cs          # Krylov 接口 + Identity/Jacobi 预条件
    │   ├── KrylovSolvers.cs           # CG / GMRES(m) / BiCGSTAB
    │   └── AMGPreconditioner.cs       # AMG Smoothed Aggregation + V-cycle
    ├── Materials/
    │   └── J2PlasticityMaterial.cs    # J2 von Mises 弹塑性本构（径向返回）
    ├── Elements/
    │   ├── FiniteElement.cs           # 单元基类（含 GPHistory）
    │   ├── CST3Element.cs             # 3 节点三角单元（完整 J2 切线）
    │   ├── LT6Element.cs              # 6 节点二次三角
    │   ├── Q4Element.cs               # 4 节点四边形
    │   └── ElementFactory.cs
    ├── Mesh/
    │   ├── FEMesh.cs                  # 网格数据容器
    │   ├── MeshReader.cs              # INP 文件导入
    │   └── StructuredMeshGenerator.cs # 自动剖分
    ├── Solvers/
    │   ├── StructuralSolver.cs        # 线弹性静力（多求解器/预条件组合）
    │   ├── NonlinearStaticSolver.cs   # J2 Newton-Raphson 强非线性
    │   └── NavierStokesSolver.cs      # SIMPLE 非稳态/定常 NS
    ├── IO/ResultWriter.cs             # CSV / VTK / PNG 导出
    ├── Rendering/ContourRenderer.cs   # GDI+ Jet 云图 + 色条 + 变形图
    └── Examples/
        └── cantilever_cst.inp         # 示例 INP（2×1 悬臂梁）
```

---

## 编译与运行

### 环境要求
- **.NET 10 SDK**（SDK 10.0.xxx）
- Windows 10/11（WinForms，net10.0-windows）
- 建议 Visual Studio 2022 17.12+（即 VS2026）

### 命令行
```bash
cd NonlinearFEM2D
dotnet build -c Release
dotnet run --project NonlinearFEM2D/NonlinearFEM2D.csproj -c Release
```

### Visual Studio
1. 启动 Visual Studio 2022；
2. 打开 `NonlinearFEM2D.slnx`（若 .slnx 打不开，使用 `NonlinearFEM2D.sln`）；
3. Ctrl+F5 运行。

---

## 使用指南

1. 启动程序后，左侧 PropertyGrid 设置几何、材料、荷载、流体参数：
   - `Length / Height / Nx / Ny`：悬臂梁尺寸与网格密度；
   - `E / Nu / Rho / Thick / SigmaY0 / Harden`：弹性模量、泊松比、密度、厚度、初始屈服应力、硬化模量；
   - `EndLoad / LoadSteps`：右端集中力、载荷步数；
   - `Viscosity / LidVelocity / CavityN / Dt / TimeSteps / AlphaU / AlphaP`：NS 参数（粘度、顶盖速度、方腔网格、时间步、欠松弛）。
2. **[网格]** 菜单：
   - `悬臂梁 CST3 (J2 弹塑性)` → 生成 CST3 悬臂梁（左固支，右顶集中力）；
   - `悬臂梁 Q4` → Q4 悬臂梁（线弹性验证）；
   - `顶盖驱动方腔 Q4 (NS)` → 1m×1m 方腔，顶盖匀速拖动。
3. **[求解]** 菜单：
   - 线性静力四种组合（CG/GMRES/BiCGSTAB+Jacobi、CG+AMG）；
   - ★ **强非线性 J2 塑性（力加载-Newton+AMG）**：推荐使用 BiCGSTAB+AMG，大网格收敛最快；
   - ★ 强非线性 J2 塑性（位移加载）：位移控加载（右端顶部竖向位移 10% 梁高），适合后屈服路径；
   - ★ 几何+材料双重非线性：开启初应力刚度 Kσ；
   - NS 定常 SIMPLE / NS 非稳态。
4. 右上"显示字段"下拉可实时切换 10 种物理量云图（Jet 色阶含数值色条）；Krylov/预条件下拉供线性与非线性求解器选择组合。
5. **[文件]→[导出结果...]** 将 CSV、VTK、PNG 输出到指定目录。

---

## 数值方法要点

### J2 弹塑性径向返回（平面应变）
1. 弹性试探 σ^tr = σ_n + D_e : Δε；
2. 偏应力 s^tr 与 von Mises 等效应力；
3. 屈服函数 f^tr = ||s^tr|| - √(2/3)σ_y(ε̄_p)；若 f^tr ≤ 0 则保持弹性；
4. 否则 Δγ = f^tr / (2μ + 2H/3)，σ_{n+1} = σ^tr - 2μ Δγ n，ε̄_p += √(2/3) Δγ；
5. 一致切线 D_alg = 2μ [(1-β)P_dev + γ_f n⊗n] + K 1⊗1。

### AMG Smoothed Aggregation
- 按 |A_ii| 降序选种子，强连接阈值 θ√(|A_ii·A_jj|) BFS 聚合；
- Tentative P（聚合指示）→ 阻尼 Jacobi 平滑 P ← (I - ω D⁻¹A) P；
- Galerkin 粗化 Ac = P^T A P（多线程三积）；
- V-cycle：前/后各 2 次加权 Jacobi 光滑；最粗层 (n<64) 稠密 Gauss-Jordan 直接求逆。

### SIMPLE
- 动量方程 (M/Δt + C(u*) + K) u* = f - G p^n；
- 压力 Poisson：L p' = -D·u* / Δt；
- 修正 u = u* - (Δt/ρ)∇p'；p^{n+1} = p^n + α_p p'；
- 速度欠松弛 u ← α_u u* + (1-α_u)u_old。

---

## 输入文件格式（.inp 示例）

```
** 注释
*NODE
1, x, y
2, x, y
...
*ELEMENT, TYPE=CPS3
1, n1, n2, n3          # 1-based 节点号，三角形 3 个角点
*MATERIAL
E 2.1e11 nu 0.3 rho 7850 t 0.01 sy 2.5e8 h 2.1e9
*BOUNDARY
nodeId, dof(1=u,2=v), value
```
CST3 单元类型关键字：`CPS3` / `TRI3` / `CST`；Q4：`CPE4` / `QUAD` / `Q4`；LT6：`CPS6` / `TRI6` / `LT6`。

---

## 典型算例建议

| 算例 | 推荐菜单 | 备注 |
|------|---------|------|
| 悬臂梁弹性验证 | 悬臂梁 Q4 → 线性静力 CG+AMG | 挠度与材料力学公式对比 |
| 悬臂梁塑性扩展 | 悬臂梁 CST3 → 强非线性 J2 塑性（力加载） | 观察 Von Mises 与 ε_p 云图；Nx=40,Ny=10 效果较好 |
| 悬臂梁后屈服路径 | 悬臂梁 CST3 → 位移加载 | 避免力加载下软化段发散 |
| 顶盖驱动方腔 | 顶盖驱动方腔 Q4 → NS 定常 | Re=ρUL/μ≈1000 时可见主涡（默认 ρ=7850,U=1,L=1 时 Re 过大，建议把 Rho 调为 100、Viscosity 0.01） |

> NS 默认参数适合高粘度 creeping 流；观察高 Re 涡结构请适当调整 Rho/Viscosity 并使用 NS 非稳态（时间步更小）。

---

## 扩展方向
- CST3 已实现完整 J2 一致切线，将径向返回逻辑复制到 Q4 的 4 个 Gauss 点即可支持四边形弹塑性；
- 弧长法（Riks/Crisfield）替换 Newton-Raphson 以跟踪极限点后屈曲路径；
- 接触/摩擦边界条件；
- 模态分析（已在早期 ModalFEM2D 版本中实现，可直接搬过来 + J2 切线模态）；
- 瞬态动力学（Newmark/HHT-α）+ J2 本构；
- 用 ILGPU/CUDA 加速 SpMV/RAP 以支撑十万级自由度。

---

## 代码规范
- 所有 C# 源码**不使用 `var`**（全量类型显式声明）；
- 命名空间按 `NonlinearFEM2D.{Core,Materials,Elements,Mesh,Solvers,IO,Rendering}` 分层；
- 项目目标框架：`net10.0-windows`，UseWindowsForms=true，AllowUnsafeBlocks=true（stackalloc 用于 SIMD gather）。

---

## 许可
本程序仅用于教学/科研示例，用户可自由修改与扩展。
