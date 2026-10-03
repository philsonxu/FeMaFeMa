# Fem2DFluid — 二维有限元流体计算程序（C# / .NET 10 / WinForm）

一个用 C# 3.0 手写的二维有限元流体计算桌面程序，适用于势流 / 渗流 / 低雷诺数 Navier-Stokes 流动分析，支持 CST 线性三角单元与 LT6 二次三角单元，稀疏矩阵存储，GMRES+ILU0 迭代求解器，SIMPLE 压力修正算法，非稳态时间推进，DXF/文本导入，彩色云图可视化。

## ✨ 功能特性

| 模块 | 功能 |
|------|------|
| 单元类型 | 3 节点 CST 线性三角单元（1 点高斯积分）；6 节点 LT6 二次三角单元（3 点高斯积分，O(h³)精度） |
| 网格处理 | TXT 节点/单元文本导入；netDxf 读取 DXF（Line/Circle/LwPolyline）；Bowyer-Watson Delaunay 自动三角剖分；边界边等距细分；多连通域（外框+内孔）多边形裁剪；LT6 边中点自动生成；节点连续重编号 |
| 物理模型 | 势流/渗流（∇·(k∇φ)=S）；Navier-Stokes 不可压缩流动（SIMPLE 算法+Rhie-Chow 动量插值+压力修正+欠松弛）；非稳态隐式欧拉时间推进；集中质量矩阵 |
| 线性求解器 | CSR 压缩稀疏行存储；COO 三元组 → CSR 转换；Parallel.For 行级锁并行装配；GMRES(30) 重启广义最小残差（Arnoldi+Givens 旋转）；ILU0 不完全 LU 预处理 |
| 输出 | TXT 求解报告（节点/单元明细+统计）；CSV 节点/单元表；VTK 非结构网格（ParaView 直接打开）；Jet 彩色云图 PNG（带右侧色条图例+网格+速度矢量+标题+坐标范围） |
| 界面 | WinForm 纯代码构建：菜单栏 / 工具栏 / 左侧模型树 / 云图画布 / 黑色日志面板 / 状态栏 |

## 📁 项目结构

```
Fem2DFluid/
├── Fem2DFluid.slnx          Visual Studio 2026 解决方案（简练格式）
├── Fem2DFluid.csproj        SDK 风格项目文件（net10.0-windows）
├── Program.cs               程序入口
├── Properties/AssemblyInfo.cs
├── Models/                  数据模型
│   ├── Node.cs              节点（坐标/BC/速度/压力）
│   ├── Element.cs           单元（CST/LT6节点、面积、形心、速度）
│   ├── Edge.cs              无序边结构体
│   ├── Material.cs          材料（渗透系数/密度/粘度）
│   ├── BoundaryEdge.cs      有向边界边
│   ├── FemModel.cs          整体模型容器
│   └── FemResult.cs         求解结果统计
├── Services/                业务逻辑
│   ├── SparseMatrix.cs      CSR/COO + ILU0 + GMRES
│   ├── FemSolver.cs         CST/LT6装配/SIMPLE/非稳态/导出
│   ├── Triangulator.cs      Delaunay剖分/边界细分/裁剪/LT6中点
│   ├── MeshImporter.cs      TXT/DXF 导入
│   ├── SampleBuilder.cs     三个内置示例
│   └── ResultExporter.cs    CSV/Tecplot/TXT/PNG 导出
└── Forms/
    └── MainForm.cs          WinForm 主界面
```

## 🚀 快速开始

1. Visual Studio 2026 打开 `Fem2DFluid.slnx`；
2. 右键解决方案 → **还原 NuGet 程序包**（自动下载 netDxf 3.0.0）；
3. **Ctrl + F5** 运行；默认加载「示例1：方腔势流验证」；
4. 工具栏：
   - 选择示例 / 打开 TXT 或 DXF 文件；
   - 设置网格尺寸 `h`（0.20 左右推荐）；
   - 勾选「LT6」使用二次单元；勾选「NS流动」启用 SIMPLE 求解 N-S 方程；
   - 点「🔺 剖分」生成网格；点「🧮 求解」开始计算；
   - 非稳态模式下点「▶ 时间步」推进 Δt；
   - 点「📄 导出」输出 TXT/CSV/VTK/PNG 结果文件。

## 📐 内置计算实例

| 示例 | 说明 | 解析/预期结果 |
|------|------|---------------|
| **示例1：方腔势流验证** | 21×21 结构化网格，上下壁 φ=0，左右 Neumann | φ=x，Vx=1, Vy=0，验证求解器正确性 |
| **示例2：圆柱绕流（NS）** | 矩形域(0,0)-(4,1) + 中心(2,0.5)半径0.2圆柱，左入口φ=0右出口φ=4，Re≈20 | 圆柱两侧速度加速，流线绕圆柱偏转 |
| **示例3：达西坝体渗流** | 21×11 矩形坝体，上游水头10，下游水头2，k=0.001 | 势从上游向下游近似线性递减 |

三个示例运行即出结果，PNG 云图均带右侧 Jet 色条图例、网格线、速度矢量、坐标标题。

## 🧮 理论说明

### CST 线性三角单元（势流）
刚度矩阵：`K_ij = (b_i·b_j + c_i·c_j)/(4Δ)·k`，其中 `b_i = y_j - y_k, c_i = x_k - x_j`（循环 i→j→k→i）。

### LT6 二次三角单元
6 个节点：3 角点 + 3 边中点；形函数 `N_1 = (2L_1-1)L_1, N_4 = 4L_1L_2`；3 点高斯积分（面积坐标权重 1/3 或 1/6）。

### SIMPLE 算法
1. 由猜测压力场 p* 解动量方程得 u*,v*；
2. Rhie-Chow 插值消除棋盘格压力振荡；
3. 解压力修正方程 p'；
4. 修正速度 u=u*+u', p=p*+α_p·p'；
5. 欠松弛 α_u=0.7, α_p=0.3，迭代至连续性残差 < 1e-5。

### GMRES(30) + ILU0
- ILU0：在原稀疏结构内做不完全 LU 分解作为预条件子；
- GMRES(m)：Arnoldi 过程构造 Krylov 子空间正交基，Givens 旋转求解最小二乘残差，每 m=30 步重启一次。

## 📄 导出格式

- **TXT 报告**：含节点数/单元数/装配时间/GMRES迭代次数/残差/各场最值，附节点 φ/u/v/p/V 值与单元形心值表；
- **CSV**：`_nodes.csv`（Id,X,Y,BC,U,V,P,Phi,Vmag）、`_elements.csv`（Id,A,B,C,Area,Cx,Cy,Vx,Vy,Vmag,P）；
- **VTK**：`.vtu` 非结构网格文件，ParaView 可直接打开做流线/等值面/动画；
- **PNG 云图**：24 级 Jet 色标填充三角单元 + 黑色网格线 + 红色速度矢量 + 右侧垂直色条（Min/Max/中值标注）+ 标题（场名+最值）+ 坐标边框与 Xmin/Xmax/Ymin/Ymax 标注。

## 🔧 依赖

- **.NET 10.0**（`net10.0-windows`）
- **netDxf 3.0.0**（NuGet 自动还原，DXF 文件读写）
- 求解器/可视化/导出全部使用 .NET 基础类库手写，无其他第三方依赖

## ⚙️ 代码规范

- 全程显式类型声明，**不使用 `var`** 关键字；
- 不使用 C# 4+ 的 `async/await`、字符串插值 `$""`、null 条件运算符 `?.`；
- 允许使用 C# 3.0 特性：`List<>`、`Dictionary<>`、匿名方法 `delegate(){}`、对象初始化器；
- WinForm 控件全部由代码动态构建，无 Designer.cs 依赖；
- 所有文本输出使用 `Encoding.UTF8`，中文界面正常显示。

## 📝 License

仅供学习与研究使用，可自由修改分发。
