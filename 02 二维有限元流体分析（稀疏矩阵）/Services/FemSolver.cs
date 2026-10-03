using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Threading;
using System.Threading.Tasks;
using System.IO;
using Fem2DFluid.Models;

namespace Fem2DFluid.Services
{
    /// <summary>
    /// 二维势流有限元求解器（稀疏矩阵版本，支持LT6二次三角形 + 稳态/非稳态）
    /// 控制方程：∇·(k ∇φ) = f  在Ω内   （稳态）
    ///          ∂φ/∂t - ∇·(k ∇φ) = f   （非稳态，隐式欧拉）
    /// Dirichlet: φ = φ0          在Γ_D
    /// Neumann:  k ∂φ/∂n = q      在Γ_N
    ///
    /// 单元：线性CST(3节点, 1点高斯) 或 二次LT6(6节点, 3点高斯)；
    /// 刚度矩阵使用 CSR 压缩稀疏行格式存储，COO三元组并行装配；
    /// ILU0预处理的 GMRES(m) 重启广义最小残差法求解。
    /// </summary>
    public class FemSolver
    {
        /// <summary>非稳态时间步长（0 = 稳态求解）</summary>
        public double Dt = 0.0;

        /// <summary>上次求解结果（用于时间推进作为初值）</summary>
        private double[] _prevPhi;

        /// <summary>求解耗时（秒）</summary>
        public double AssemblyTimeMs;
        public double SolveTimeMs;
        public int Iterations;
        public double FinalResidual;

        public static FemResult Solve(FemModel model)
        {
            FemSolver solver = new FemSolver();
            return solver.SolveSteady(model);
        }

        /// <summary>稳态求解</summary>
        public FemResult SolveSteady(FemModel model)
        {
            Dt = 0.0;
            _prevPhi = null;
            return RunSolve(model);
        }

        /// <summary>非稳态时间步进：从上一步phi出发推进Dt秒</summary>
        public FemResult AdvanceTime(FemModel model, double dt)
        {
            Dt = dt;
            // 保存当前phi作为初值
            int n = model.Nodes.Count;
            _prevPhi = new double[n];
            Dictionary<int, int> idxMap = BuildIndexMap(model);
            for (int i = 0; i < n; i++) _prevPhi[i] = model.Nodes[i].Phi;
            return RunSolve(model);
        }

        private FemResult RunSolve(FemModel model)
        {
            FemResult result = new FemResult();
            Stopwatch swTotal = Stopwatch.StartNew();
            int n = model.Nodes.Count;
            if (n == 0) { result.Success = false; result.Message = "没有节点"; return result; }

            Dictionary<int, int> idxMap = BuildIndexMap(model);
            bool useLT6 = model.Elements.Count > 0 && model.Elements[0].Order >= 1;

            // ====== COO 三元组并行装配稀疏刚度矩阵K与质量矩阵M ======
            Stopwatch swAsm = Stopwatch.StartNew();
            CooMatrix cooK = new CooMatrix(n);
            CooMatrix cooM = new CooMatrix(n);
            double[] F = new double[n];
            object lockF = new object();

            int elemCount = model.Elements.Count;
            Parallel.For(0, elemCount, ei =>
            {
                Element e = model.Elements[ei];
                Material mat = model.GetMaterial(e.MatId) ?? model.Materials[0];
                int[] ids = e.NodeIds6();
                Node[] nodes = new Node[6];
                int[] idxs = new int[6];
                bool ok = true;
                for (int k = 0; k < 6; k++)
                {
                    if (ids[k] <= 0) { nodes[k] = null; idxs[k] = -1; }
                    else
                    {
                        nodes[k] = model.GetNode(ids[k]);
                        if (nodes[k] == null) { ok = false; break; }
                        idxs[k] = idxMap[ids[k]];
                    }
                }
                if (!ok) return;

                double[,] Ke, Me;
                double[] Fe;
                double area;
                if (e.Order >= 1 && nodes[3] != null && nodes[4] != null && nodes[5] != null)
                    AssembleLT6(nodes, mat, out Ke, out Me, out Fe, out area);
                else
                    AssembleCST(nodes, mat, out Ke, out Me, out Fe, out area);
                if (area <= 1e-14) return;

                // 写入COO
                int nloc = (e.Order >= 1) ? 6 : 3;
                for (int r = 0; r < nloc; r++)
                {
                    int ir = idxs[r];
                    if (ir < 0) continue;
                    lock (lockF) { F[ir] += Fe[r]; }
                    for (int c = 0; c < nloc; c++)
                    {
                        int ic = idxs[c];
                        if (ic < 0) continue;
                        cooK.Add(ir, ic, Ke[r, c]);
                        cooM.Add(ir, ic, Me[r, c]);
                    }
                }
            });

            CsrMatrix K = cooK.ToCsr();
            CsrMatrix M = Dt > 1e-14 ? cooM.ToCsr() : null;

            // 非稳态：形成 A = M/Dt + K；右端 F = M/Dt * phi_old + F
            CsrMatrix A;
            double[] rhs = new double[n];
            if (Dt > 1e-14 && M != null && _prevPhi != null)
            {
                CsrMatrix Md = M.Scale(1.0 / Dt);
                A = CsrMatrix.Add(Md, K);
                double[] Mphi = new double[n];
                Md.Mult(_prevPhi, Mphi);
                for (int i = 0; i < n; i++) rhs[i] = Mphi[i] + F[i];
            }
            else
            {
                A = K;
                for (int i = 0; i < n; i++) rhs[i] = F[i];
            }
            swAsm.Stop();
            AssemblyTimeMs = swAsm.Elapsed.TotalMilliseconds;

            // ====== 施加Dirichlet边界条件（置1法，稀疏版） ======
            List<int> dirichletNodes = new List<int>();
            bool[] isD = new bool[n];
            for (int i = 0; i < n; i++)
            {
                if (model.Nodes[i].BCType == 1)
                {
                    dirichletNodes.Add(i);
                    isD[i] = true;
                }
            }
            if (dirichletNodes.Count == 0)
            {
                result.Success = false;
                result.Message = "未指定Dirichlet边界条件，问题无唯一解";
                return result;
            }
            for (int di = 0; di < dirichletNodes.Count; di++)
            {
                int i = dirichletNodes[di];
                double val = model.Nodes[i].BCValue;
                if (Dt > 1e-14 && _prevPhi != null)
                {
                    // 非稳态Dirichlet固定为目标值
                }
                for (int j = 0; j < n; j++)
                {
                    if (j == i || isD[j]) continue;
                    double v = A.Get(j, i);
                    if (Math.Abs(v) > 1e-18)
                    {
                        rhs[j] -= v * val;
                        A.SetZero(j, i);
                    }
                }
                A.ZeroRow(i);
                A.Set(i, i, 1.0);
                rhs[i] = val;
            }

            // ====== Neumann边界（第二类）：按相邻Neumann节点距离近似线积分 ======
            for (int i = 0; i < n; i++)
            {
                Node ni = model.Nodes[i];
                if (ni.BCType == 2)
                {
                    double sumLen = 0.0;
                    int cnt = 0;
                    for (int j = 0; j < n; j++)
                    {
                        if (i == j) continue;
                        Node nj = model.Nodes[j];
                        if (nj.BCType == 2)
                        {
                            double dx = nj.X - ni.X, dy = nj.Y - ni.Y;
                            double dist = Math.Sqrt(dx * dx + dy * dy);
                            if (dist < 1e-8) continue;
                            sumLen += dist;
                            cnt++;
                            if (cnt >= 6) break;
                        }
                    }
                    double halfL = sumLen * 0.5;
                    rhs[i] += ni.BCValue * halfL;
                }
            }

            // ====== ILU0 预处理 + GMRES(m) 求解 ======
            Stopwatch swSol = Stopwatch.StartNew();
            Ilu0 precond = new Ilu0(A);
            double[] phi = new double[n];
            int iter;
            double resNorm;
            bool conv = Gmres.Solve(A, rhs, phi, precond, 1e-10, 2000, 30, out iter, out resNorm);
            swSol.Stop();
            SolveTimeMs = swSol.Elapsed.TotalMilliseconds;
            Iterations = iter;
            FinalResidual = resNorm;

            if (!conv)
            {
                result.Success = false;
                result.Message = "求解失败：GMRES迭代未收敛（" + iter + "步，残差" + resNorm.ToString("E3") + "）";
                return result;
            }
            for (int i = 0; i < n; i++) model.Nodes[i].Phi = phi[i];

            // ====== 计算速度（单元常梯度 或 LT6 节点梯度） ======
            double phiMin = double.MaxValue, phiMax = double.MinValue;
            double vMin = double.MaxValue, vMax = double.MinValue;
            if (useLT6)
                ComputeVelocitiesLT6(model, idxMap);
            else
                ComputeVelocitiesCST(model);

            for (int i = 0; i < n; i++)
            {
                if (model.Nodes[i].Phi < phiMin) phiMin = model.Nodes[i].Phi;
                if (model.Nodes[i].Phi > phiMax) phiMax = model.Nodes[i].Phi;
            }
            for (int ei = 0; ei < model.Elements.Count; ei++)
            {
                double vm = model.Elements[ei].Vmag;
                if (vm < vMin) vMin = vm;
                if (vm > vMax) vMax = vm;
            }
            swTotal.Stop();
            result.Success = true;
            result.Message = "求解成功（" + (useLT6 ? "LT6" : "CST") + "单元, GMRES迭代" + iter + "步, 残差" + resNorm.ToString("E2") + "）";
            result.NodeCount = n;
            result.ElementCount = model.Elements.Count;
            result.PhiMin = phiMin;
            result.PhiMax = phiMax;
            result.VmagMin = vMin;
            result.VmagMax = vMax;
            result.SolveTime = swTotal.Elapsed;
            return result;
        }

        private Dictionary<int, int> BuildIndexMap(FemModel model)
        {
            int n = model.Nodes.Count;
            Dictionary<int, int> idxMap = new Dictionary<int, int>();
            for (int i = 0; i < n; i++) idxMap[model.Nodes[i].Id] = i;
            return idxMap;
        }

        // ========= CST线性单元装配：1点高斯积分 =========
        private static void AssembleCST(Node[] nodes, Material mat,
            out double[,] Ke, out double[,] Me, out double[] Fe, out double area)
        {
            Node n1 = nodes[0], n2 = nodes[1], n3 = nodes[2];
            double x1 = n1.X, y1 = n1.Y;
            double x2 = n2.X, y2 = n2.Y;
            double x3 = n3.X, y3 = n3.Y;
            double ar2 = (x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1);
            area = 0.5 * Math.Abs(ar2);
            double A = Math.Abs(ar2) * 0.5;
            double b1 = y2 - y3, b2 = y3 - y1, b3 = y1 - y2;
            double c1 = x3 - x2, c2 = x1 - x3, c3 = x2 - x1;
            double coef = mat.K / (4 * A);
            Ke = new double[3, 3];
            Me = new double[3, 3];
            Fe = new double[3];
            int[] idx = new int[] { 0, 1, 2 };
            double[] bb = new double[] { b1, b2, b3 };
            double[] cc = new double[] { c1, c2, c3 };
            for (int r = 0; r < 3; r++)
            {
                for (int c = 0; c < 3; c++)
                    Ke[r, c] = coef * (bb[r] * bb[c] + cc[r] * cc[c]);
            }
            // 集中质量矩阵：M_ii = Δ/3（行和法）
            double mDiag = A / 3.0;
            for (int r = 0; r < 3; r++) Me[r, r] = mDiag;
            // 源项
            double fLoad = mat.Source * A / 3.0;
            for (int r = 0; r < 3; r++) Fe[r] = fLoad;
        }

        // ========= LT6二次三角形装配：3点高斯积分（面积坐标） =========
        // 6节点顺序：角点(1,2,3) 对应 L1,L2,L3，边中点：4(边1-2), 5(边2-3), 6(边3-1)
        // 二次形函数：N_i = (2L_i - 1)L_i (i=1,2,3角点) ; N_{i+3} = 4L_i L_j
        // 3点高斯积分点在面积坐标下：(2/3,1/6,1/6),(1/6,2/3,1/6),(1/6,1/6,2/3)，权均为 1/3
        private static void AssembleLT6(Node[] nodes, Material mat,
            out double[,] Ke, out double[,] Me, out double[] Fe, out double area)
        {
            double x1 = nodes[0].X, y1 = nodes[0].Y;
            double x2 = nodes[1].X, y2 = nodes[1].Y;
            double x3 = nodes[2].X, y3 = nodes[2].Y;
            double ar2 = (x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1);
            area = 0.5 * Math.Abs(ar2);
            double A = Math.Abs(ar2) * 0.5;
            // 形函数对面积坐标的导数（常数）
            // ∂N/∂x 与 ∂N/∂y 需通过 Jacobian 转换。简化：在每个高斯点计算。
            Ke = new double[6, 6];
            Me = new double[6, 6];
            Fe = new double[6];
            // 3点高斯
            double[,] gp = new double[,] { { 2.0 / 3, 1.0 / 6, 1.0 / 6 }, { 1.0 / 6, 2.0 / 3, 1.0 / 6 }, { 1.0 / 6, 1.0 / 6, 2.0 / 3 } };
            double[] wt = new double[] { 1.0 / 3, 1.0 / 3, 1.0 / 3 };
            // 角点坐标数组
            double[] px = new double[] { x1, x2, x3 };
            double[] py = new double[] { y1, y2, y3 };
            // 在每个高斯点计算形函数导数/值
            for (int ip = 0; ip < 3; ip++)
            {
                double L1 = gp[ip, 0], L2 = gp[ip, 1], L3 = gp[ip, 2];
                // 形函数值
                double[] Nv = new double[6];
                Nv[0] = (2 * L1 - 1) * L1;
                Nv[1] = (2 * L2 - 1) * L2;
                Nv[2] = (2 * L3 - 1) * L3;
                Nv[3] = 4 * L1 * L2;
                Nv[4] = 4 * L2 * L3;
                Nv[5] = 4 * L3 * L1;
                // 形函数对面积坐标的偏导：dN/dL1, dN/dL2, dN/dL3
                // dN1/dL1=4L1-1, dN1/dL2=0, dN1/dL3=0 等
                double[,] dNdL = new double[6, 3];
                dNdL[0, 0] = 4 * L1 - 1; dNdL[0, 1] = 0;          dNdL[0, 2] = 0;
                dNdL[1, 0] = 0;          dNdL[1, 1] = 4 * L2 - 1; dNdL[1, 2] = 0;
                dNdL[2, 0] = 0;          dNdL[2, 1] = 0;          dNdL[2, 2] = 4 * L3 - 1;
                dNdL[3, 0] = 4 * L2;     dNdL[3, 1] = 4 * L1;     dNdL[3, 2] = 0;
                dNdL[4, 0] = 0;          dNdL[4, 1] = 4 * L3;     dNdL[4, 2] = 4 * L2;
                dNdL[5, 0] = 4 * L3;     dNdL[5, 1] = 0;          dNdL[5, 2] = 4 * L1;
                // Jacobian矩阵：∂(x,y)/∂(L1,L2) （L3=1-L1-L2）
                // x = Σ N_i x_i ; y = Σ N_i y_i; 但对线性几何映射（3角点即可）J是常数
                // ∂x/∂L1 = x1 - x3; ∂x/∂L2 = x2 - x3
                // ∂y/∂L1 = y1 - y3; ∂y/∂L2 = y2 - y3
                double dxdl1 = x1 - x3, dxdl2 = x2 - x3;
                double dydl1 = y1 - y3, dydl2 = y2 - y3;
                double detJ = dxdl1 * dydl2 - dxdl2 * dydl1; // = 2*Δ
                double invDet = 1.0 / detJ;
                // J^{-1}: ∂L/∂x = 1/detJ * (dydl2, -dxdl2; -dydl1, dxdl1)
                double dL1dx = dydl2 * invDet, dL1dy = -dxdl2 * invDet;
                double dL2dx = -dydl1 * invDet, dL2dy = dxdl1 * invDet;
                double dL3dx = -dL1dx - dL2dx, dL3dy = -dL1dy - dL2dy;
                double w = wt[ip] * detJ; // detJ 已是2Δ，乘权后dΩ
                for (int a = 0; a < 6; a++)
                {
                    double dNadx = dNdL[a, 0] * dL1dx + dNdL[a, 1] * dL2dx + dNdL[a, 2] * dL3dx;
                    double dNady = dNdL[a, 0] * dL1dy + dNdL[a, 1] * dL2dy + dNdL[a, 2] * dL3dy;
                    for (int b = 0; b < 6; b++)
                    {
                        double dNbdx = dNdL[b, 0] * dL1dx + dNdL[b, 1] * dL2dx + dNdL[b, 2] * dL3dx;
                        double dNbdy = dNdL[b, 0] * dL1dy + dNdL[b, 1] * dL2dy + dNdL[b, 2] * dL3dy;
                        Ke[a, b] += mat.K * (dNadx * dNbdx + dNady * dNbdy) * w;
                        Me[a, b] += Nv[a] * Nv[b] * w;
                    }
                    Fe[a] += mat.Source * Nv[a] * w;
                }
            }
            // 集中质量（行和法）：M_ii = Σ_j M_ij
            double[,] MLump = new double[6, 6];
            for (int a = 0; a < 6; a++)
            {
                double rowSum = 0.0;
                for (int b = 0; b < 6; b++) rowSum += Me[a, b];
                MLump[a, a] = rowSum;
            }
            Me = MLump;
        }

        // ========= 速度计算：CST（单元常梯度） =========
        private static void ComputeVelocitiesCST(FemModel model)
        {
            for (int ei = 0; ei < model.Elements.Count; ei++)
            {
                Element e = model.Elements[ei];
                Node n1 = model.GetNode(e.N1);
                Node n2 = model.GetNode(e.N2);
                Node n3 = model.GetNode(e.N3);
                if (n1 == null || n2 == null || n3 == null) continue;
                double x1 = n1.X, y1 = n1.Y;
                double x2 = n2.X, y2 = n2.Y;
                double x3 = n3.X, y3 = n3.Y;
                double area = 0.5 * ((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1));
                if (Math.Abs(area) < 1e-12) area = 1e-12;
                double b1 = y2 - y3, b2 = y3 - y1, b3 = y1 - y2;
                double c1 = x3 - x2, c2 = x1 - x3, c3 = x2 - x1;
                double vx = (b1 * n1.Phi + b2 * n2.Phi + b3 * n3.Phi) / (2 * area);
                double vy = (c1 * n1.Phi + c2 * n2.Phi + c3 * n3.Phi) / (2 * area);
                e.Vx = vx;
                e.Vy = vy;
                e.Vmag = Math.Sqrt(vx * vx + vy * vy);
                Triangulator.ComputeElementGeom(e, model);
            }
        }

        // ========= 速度计算：LT6 单元节点梯度（在3个内部高斯点平均后赋给节点） =========
        private static void ComputeVelocitiesLT6(FemModel model, Dictionary<int, int> idxMap)
        {
            int n = model.Nodes.Count;
            double[] vxSum = new double[n];
            double[] vySum = new double[n];
            int[] count = new int[n];

            for (int ei = 0; ei < model.Elements.Count; ei++)
            {
                Element e = model.Elements[ei];
                if (e.Order < 1) { ComputeVelocitiesCST(model); continue; }
                int[] ids = e.NodeIds6();
                Node[] nds = new Node[6];
                for (int k = 0; k < 6; k++) nds[k] = model.GetNode(ids[k]);
                if (nds[3] == null || nds[4] == null || nds[5] == null) { ComputeVelocitiesCST(model); continue; }

                double x1 = nds[0].X, y1 = nds[0].Y;
                double x2 = nds[1].X, y2 = nds[1].Y;
                double x3 = nds[2].X, y3 = nds[2].Y;
                double dxdl1 = x1 - x3, dxdl2 = x2 - x3;
                double dydl1 = y1 - y3, dydl2 = y2 - y3;
                double detJ = dxdl1 * dydl2 - dxdl2 * dydl1;
                double invDet = 1.0 / detJ;
                double dL1dx = dydl2 * invDet, dL1dy = -dxdl2 * invDet;
                double dL2dx = -dydl1 * invDet, dL2dy = dxdl1 * invDet;
                double dL3dx = -dL1dx - dL2dx, dL3dy = -dL1dy - dL2dy;

                // 3点高斯点梯度
                double[,] gp = new double[,] { { 2.0 / 3, 1.0 / 6, 1.0 / 6 }, { 1.0 / 6, 2.0 / 3, 1.0 / 6 }, { 1.0 / 6, 1.0 / 6, 2.0 / 3 } };
                double gvx = 0, gvy = 0;
                for (int ip = 0; ip < 3; ip++)
                {
                    double L1 = gp[ip, 0], L2 = gp[ip, 1], L3 = gp[ip, 2];
                    double[,] dNdL = new double[6, 3];
                    dNdL[0, 0] = 4 * L1 - 1; dNdL[1, 1] = 4 * L2 - 1; dNdL[2, 2] = 4 * L3 - 1;
                    dNdL[3, 0] = 4 * L2; dNdL[3, 1] = 4 * L1;
                    dNdL[4, 1] = 4 * L3; dNdL[4, 2] = 4 * L2;
                    dNdL[5, 0] = 4 * L3; dNdL[5, 2] = 4 * L1;
                    double vxp = 0, vyp = 0;
                    for (int a = 0; a < 6; a++)
                    {
                        double dNadx = dNdL[a, 0] * dL1dx + dNdL[a, 1] * dL2dx + dNdL[a, 2] * dL3dx;
                        double dNady = dNdL[a, 0] * dL1dy + dNdL[a, 1] * dL2dy + dNdL[a, 2] * dL3dy;
                        vxp += dNadx * nds[a].Phi;
                        vyp += dNady * nds[a].Phi;
                    }
                    gvx += vxp; gvy += vyp;
                }
                gvx /= 3.0; gvy /= 3.0;
                e.Vx = gvx; e.Vy = gvy; e.Vmag = Math.Sqrt(gvx * gvx + gvy * gvy);
                Triangulator.ComputeElementGeom(e, model);

                // 累加到6个节点做平均
                for (int k = 0; k < 6; k++)
                {
                    int gi = idxMap[ids[k]];
                    vxSum[gi] += gvx; vySum[gi] += gvy; count[gi]++;
                }
            }
            // 节点速度 = 周围单元平均
            for (int i = 0; i < n; i++)
            {
                if (count[i] > 0)
                {
                    model.Nodes[i].Vx = vxSum[i] / count[i];
                    model.Nodes[i].Vy = vySum[i] / count[i];
                }
            }
        }

        /// <summary>
        /// 将计算结果导出为文本报告文件（.txt）+ 节点/单元数据CSV + VTK（ParaView可视化）
        /// </summary>
        public static void ExportTextResults(FemModel model, FemResult result, string outDir)
        {
            if (!Directory.Exists(outDir)) Directory.CreateDirectory(outDir);
            string timestamp = DateTime.Now.ToString("yyyyMMdd_HHmmss");

            // 1. 文本报告
            string reportPath = Path.Combine(outDir, "Report_" + timestamp + ".txt");
            using (StreamWriter sw = new StreamWriter(reportPath))
            {
                sw.WriteLine("===========================================================");
                sw.WriteLine("   FEM 2D Fluid Analysis Report  二维势流有限元计算报告");
                sw.WriteLine("===========================================================");
                sw.WriteLine("生成时间: " + DateTime.Now.ToString("yyyy-MM-dd HH:mm:ss"));
                sw.WriteLine("节点数:   " + result.NodeCount);
                sw.WriteLine("单元数:   " + result.ElementCount);
                sw.WriteLine("单元类型: " + (model.Elements.Count > 0 && model.Elements[0].Order >= 1 ? "LT6二次三角形(6节点)" : "CST线性三角形(3节点)"));
                sw.WriteLine("装配耗时: {0:F2} ms", result.SolveTime.TotalMilliseconds - 0); // 使用总耗时近似
                sw.WriteLine("求解迭代: " + (model.Elements.Count > 0 ? "" : ""));
                sw.WriteLine("势 φ 范围: [{0:F6}, {1:F6}]", result.PhiMin, result.PhiMax);
                sw.WriteLine("速度|V| 范围: [{0:F6}, {1:F6}]", result.VmagMin, result.VmagMax);
                sw.WriteLine("求解信息: " + result.Message);
                sw.WriteLine("-----------------------------------------------------------");
                sw.WriteLine("节点列表 (Id, X, Y, BCType, BCValue, Phi, Vx, Vy, Vmag):");
                sw.WriteLine("-----------------------------------------------------------");
                for (int i = 0; i < model.Nodes.Count; i++)
                {
                    Node nd = model.Nodes[i];
                    double vx = nd.Vx, vy = nd.Vy, vm = Math.Sqrt(vx * vx + vy * vy);
                    // 若节点无速度（CST模式），取周围单元平均
                    if (nd.NodeType == 0 && Math.Abs(vx) < 1e-12 && Math.Abs(vy) < 1e-12)
                    {
                        double sx = 0, sy = 0; int cnt = 0;
                        for (int e = 0; e < model.Elements.Count; e++)
                        {
                            Element el = model.Elements[e];
                            if (el.N1 == nd.Id || el.N2 == nd.Id || el.N3 == nd.Id)
                            { sx += el.Vx; sy += el.Vy; cnt++; }
                        }
                        if (cnt > 0) { vx = sx / cnt; vy = sy / cnt; vm = Math.Sqrt(vx * vx + vy * vy); }
                    }
                    sw.WriteLine("{0,6}  {1,10:F6}  {2,10:F6}  BC={3}  BCv={4,8:F4}  φ={5,10:F6}  Vx={6,10:F6}  Vy={7,10:F6}  |V|={8,10:F6}",
                        nd.Id, nd.X, nd.Y, nd.BCType, nd.BCValue, nd.Phi, vx, vy, vm);
                }
                sw.WriteLine("-----------------------------------------------------------");
                sw.WriteLine("单元列表 (Id, N1,N2,N3, 面积, Vx, Vy, Vmag, 中心):");
                sw.WriteLine("-----------------------------------------------------------");
                for (int e = 0; e < model.Elements.Count; e++)
                {
                    Element el = model.Elements[e];
                    sw.WriteLine("{0,5}  ({1,5},{2,5},{3,5})  Area={4,10:F6}  Vx={5,10:F6}  Vy={6,10:F6}  |V|={7,10:F6}  C=({8:F4},{9:F4})",
                        el.Id, el.N1, el.N2, el.N3, el.Area, el.Vx, el.Vy, el.Vmag, el.Cx, el.Cy);
                }
                sw.WriteLine("===========================================================");
                sw.WriteLine("报告结束");
            }

            // 2. 节点CSV
            string nodeCsv = Path.Combine(outDir, "Nodes_" + timestamp + ".csv");
            using (StreamWriter sw = new StreamWriter(nodeCsv))
            {
                sw.WriteLine("Id,X,Y,BCType,BCValue,Phi");
                for (int i = 0; i < model.Nodes.Count; i++)
                {
                    Node nd = model.Nodes[i];
                    sw.WriteLine("{0},{1},{2},{3},{4},{5}", nd.Id, nd.X, nd.Y, nd.BCType, nd.BCValue, nd.Phi);
                }
            }

            // 3. 单元CSV
            string elemCsv = Path.Combine(outDir, "Elements_" + timestamp + ".csv");
            using (StreamWriter sw = new StreamWriter(elemCsv))
            {
                sw.WriteLine("Id,N1,N2,N3,N4,N5,N6,Area,Cx,Cy,Vx,Vy,Vmag");
                for (int e = 0; e < model.Elements.Count; e++)
                {
                    Element el = model.Elements[e];
                    sw.WriteLine("{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12}",
                        el.Id, el.N1, el.N2, el.N3, el.N4, el.N5, el.N6, el.Area, el.Cx, el.Cy, el.Vx, el.Vy, el.Vmag);
                }
            }

            // 4. VTK（非结构网格，ParaView直接打开）
            string vtkPath = Path.Combine(outDir, "Result_" + timestamp + ".vtk");
            Dictionary<int, int> idx = new Dictionary<int, int>();
            for (int i = 0; i < model.Nodes.Count; i++) idx[model.Nodes[i].Id] = i;
            using (StreamWriter sw = new StreamWriter(vtkPath))
            {
                sw.WriteLine("# vtk DataFile Version 3.0");
                sw.WriteLine("Fem2D Fluid Result");
                sw.WriteLine("ASCII");
                sw.WriteLine("DATASET UNSTRUCTURED_GRID");
                sw.WriteLine("POINTS " + model.Nodes.Count + " double");
                for (int i = 0; i < model.Nodes.Count; i++)
                    sw.WriteLine("{0} {1} 0.0", model.Nodes[i].X, model.Nodes[i].Y);
                int connLen = 0;
                for (int e = 0; e < model.Elements.Count; e++)
                    connLen += 1 + (model.Elements[e].Order >= 1 && model.Elements[e].N4 > 0 ? 6 : 3);
                sw.WriteLine("CELLS " + model.Elements.Count + " " + connLen);
                for (int e = 0; e < model.Elements.Count; e++)
                {
                    Element el = model.Elements[e];
                    if (el.Order >= 1 && el.N4 > 0 && el.N5 > 0 && el.N6 > 0)
                        sw.WriteLine("6 {0} {1} {2} {3} {4} {5}",
                            idx[el.N1], idx[el.N2], idx[el.N3], idx[el.N4], idx[el.N5], idx[el.N6]);
                    else
                        sw.WriteLine("3 {0} {1} {2}", idx[el.N1], idx[el.N2], idx[el.N3]);
                }
                sw.WriteLine("CELL_TYPES " + model.Elements.Count);
                for (int e = 0; e < model.Elements.Count; e++)
                    sw.WriteLine(model.Elements[e].Order >= 1 && model.Elements[e].N4 > 0 ? 22 : 5);
                sw.WriteLine("POINT_DATA " + model.Nodes.Count);
                sw.WriteLine("SCALARS Phi double 1");
                sw.WriteLine("LOOKUP_TABLE default");
                for (int i = 0; i < model.Nodes.Count; i++) sw.WriteLine(model.Nodes[i].Phi);
                sw.WriteLine("VECTORS Velocity double");
                for (int i = 0; i < model.Nodes.Count; i++)
                    sw.WriteLine("{0} {1} 0.0", model.Nodes[i].Vx, model.Nodes[i].Vy);
            }
        }
    }

    // ====================================================================
    // COO 三元组（带线程安全Add，原子累加）
    // ====================================================================
    internal class CooMatrix
    {
        private int _n;
        private List<int> _ri, _ci;
        private List<double> _vv;
        private object _lock = new object();

        public CooMatrix(int n)
        {
            _n = n;
            _ri = new List<int>();
            _ci = new List<int>();
            _vv = new List<double>();
        }

        public void Add(int i, int j, double v)
        {
            if (Math.Abs(v) < 1e-30) return;
            lock (_lock)
            {
                _ri.Add(i);
                _ci.Add(j);
                _vv.Add(v);
            }
        }

        public CsrMatrix ToCsr()
        {
            int nnz = _ri.Count;
            int[] rowCount = new int[_n + 1];
            for (int k = 0; k < nnz; k++) rowCount[_ri[k] + 1]++;
            int[] rowPtr = new int[_n + 1];
            for (int i = 0; i < _n; i++) rowPtr[i + 1] = rowPtr[i] + rowCount[i + 1];
            int[] cols = new int[nnz];
            double[] vals = new double[nnz];
            int[] cursor = new int[_n];
            for (int i = 0; i < _n; i++) cursor[i] = rowPtr[i];
            for (int k = 0; k < nnz; k++)
            {
                int r = _ri[k];
                int pos = cursor[r]++;
                cols[pos] = _ci[k];
                vals[pos] = _vv[k];
            }
            for (int i = 0; i < _n; i++)
            {
                int p0 = rowPtr[i];
                int p1 = rowPtr[i + 1];
                if (p1 - p0 <= 1) continue;
                for (int a = p0 + 1; a < p1; a++)
                {
                    int tc = cols[a]; double tv = vals[a];
                    int b = a - 1;
                    while (b >= p0 && cols[b] > tc)
                    { cols[b + 1] = cols[b]; vals[b + 1] = vals[b]; b--; }
                    cols[b + 1] = tc; vals[b + 1] = tv;
                }
                int write = p0;
                for (int a = p0; a < p1; a++)
                {
                    if (write > p0 && cols[write - 1] == cols[a])
                        vals[write - 1] += vals[a];
                    else
                    { cols[write] = cols[a]; vals[write] = vals[a]; write++; }
                }
                for (int a = write; a < p1; a++) { vals[a] = 0.0; cols[a] = -1; }
                rowPtr[i + 1] = write + (rowPtr[i + 1] - p1);
            }
            int newNnz = rowPtr[_n];
            int[] newCols = new int[newNnz];
            double[] newVals = new double[newNnz];
            int[] newRowPtr = new int[_n + 1];
            int wp = 0;
            for (int i = 0; i < _n; i++)
            {
                newRowPtr[i] = wp;
                for (int k = rowPtr[i]; k < rowPtr[i + 1]; k++)
                {
                    if (cols[k] >= 0)
                    { newCols[wp] = cols[k]; newVals[wp] = vals[k]; wp++; }
                }
            }
            newRowPtr[_n] = wp;
            return new CsrMatrix(_n, newRowPtr, newCols, newVals);
        }
    }

    // ====================================================================
    // CSR 压缩稀疏行
    // ====================================================================
    internal class CsrMatrix
    {
        public int N;
        public int[] RowPtr;
        public int[] Cols;
        public double[] Vals;

        public CsrMatrix(int n, int[] rp, int[] c, double[] v)
        { N = n; RowPtr = rp; Cols = c; Vals = v; }

        public int FindPos(int i, int j)
        {
            int p0 = RowPtr[i], p1 = RowPtr[i + 1];
            int lo = p0, hi = p1 - 1;
            while (lo <= hi)
            {
                int mid = (lo + hi) >> 1;
                if (Cols[mid] == j) return mid;
                if (Cols[mid] < j) lo = mid + 1; else hi = mid - 1;
            }
            return -1;
        }
        public double Get(int i, int j) { int p = FindPos(i, j); return p < 0 ? 0.0 : Vals[p]; }
        public void Set(int i, int j, double v) { int p = FindPos(i, j); if (p >= 0) Vals[p] = v; }
        public void SetZero(int i, int j) { int p = FindPos(i, j); if (p >= 0) Vals[p] = 0.0; }
        public void ZeroRow(int i)
        { int p0 = RowPtr[i], p1 = RowPtr[i + 1]; for (int k = p0; k < p1; k++) Vals[k] = 0.0; }

        public void Mult(double[] x, double[] y)
        {
            for (int i = 0; i < N; i++)
            {
                double s = 0.0;
                for (int k = RowPtr[i]; k < RowPtr[i + 1]; k++) s += Vals[k] * x[Cols[k]];
                y[i] = s;
            }
        }

        public CsrMatrix Scale(double alpha)
        {
            double[] nv = new double[Vals.Length];
            for (int k = 0; k < Vals.Length; k++) nv[k] = Vals[k] * alpha;
            return new CsrMatrix(N, (int[])RowPtr.Clone(), (int[])Cols.Clone(), nv);
        }

        public static CsrMatrix Add(CsrMatrix A, CsrMatrix B)
        {
            // 简单实现：把B加到A（输出新矩阵），保持A结构+扩展B列
            CooMatrix coo = new CooMatrix(A.N);
            for (int i = 0; i < A.N; i++)
            {
                for (int k = A.RowPtr[i]; k < A.RowPtr[i + 1]; k++) coo.Add(i, A.Cols[k], A.Vals[k]);
                for (int k = B.RowPtr[i]; k < B.RowPtr[i + 1]; k++) coo.Add(i, B.Cols[k], B.Vals[k]);
            }
            return coo.ToCsr();
        }
    }

    // ====================================================================
    // ILU0 不完全LU预处理（标准Saad算法）
    // ====================================================================
    internal class Ilu0
    {
        private int _n;
        private int[] _rowPtr, _cols;
        private double[] _lu;
        private int[] _diagPos;

        public Ilu0(CsrMatrix A)
        {
            _n = A.N;
            _rowPtr = (int[])A.RowPtr.Clone();
            _cols = (int[])A.Cols.Clone();
            _lu = (double[])A.Vals.Clone();
            _diagPos = new int[_n];
            for (int i = 0; i < _n; i++)
            {
                _diagPos[i] = -1;
                for (int k = _rowPtr[i]; k < _rowPtr[i + 1]; k++)
                    if (_cols[k] == i) { _diagPos[i] = k; break; }
            }
            for (int i = 0; i < _n; i++)
            {
                int d = _diagPos[i];
                if (d < 0 || Math.Abs(_lu[d]) < 1e-30) continue;
                for (int kk = _rowPtr[i]; kk < d; kk++)
                {
                    int k = _cols[kk];
                    if (k < 0 || k >= i) continue;
                    double lik = _lu[kk] / _lu[_diagPos[k]];
                    _lu[kk] = lik;
                    int kd = _diagPos[k];
                    for (int kj = kd + 1; kj < _rowPtr[k + 1]; kj++)
                    {
                        int j = _cols[kj];
                        int ip = FindInRow(i, j);
                        if (ip >= 0) _lu[ip] -= lik * _lu[kj];
                    }
                }
            }
        }
        private int FindInRow(int i, int j)
        {
            for (int k = _rowPtr[i]; k < _rowPtr[i + 1]; k++)
            {
                if (_cols[k] == j) return k;
                if (_cols[k] > j) return -1;
            }
            return -1;
        }
        public void Solve(double[] r, double[] z)
        {
            double[] y = new double[_n];
            for (int i = 0; i < _n; i++)
            {
                double s = r[i];
                for (int k = _rowPtr[i]; k < _diagPos[i]; k++) s -= _lu[k] * y[_cols[k]];
                y[i] = s;
            }
            for (int i = _n - 1; i >= 0; i--)
            {
                double s = y[i];
                for (int k = _diagPos[i] + 1; k < _rowPtr[i + 1]; k++) s -= _lu[k] * z[_cols[k]];
                z[i] = s / _lu[_diagPos[i]];
            }
        }
    }

    // ====================================================================
    // GMRES(m) 重启广义最小残差法（右预处理）
    // 参考 Saad, Algorithm 6.9
    // ====================================================================
    internal static class Gmres
    {
        public static bool Solve(CsrMatrix A, double[] b, double[] x, Ilu0 M,
            double tol, int maxIter, int mRestart, out int iter, out double resNorm)
        {
            int n = A.N;
            iter = 0;
            resNorm = 0.0;
            double[] r = new double[n];
            double[] w = new double[n];
            double bnorm = 0.0;
            for (int i = 0; i < n; i++) bnorm += b[i] * b[i];
            bnorm = Math.Sqrt(bnorm);
            if (bnorm < 1e-30) bnorm = 1.0;

            // Arnoldi 过程所需：V 是 (m+1) 个 n-向量；H 是 (m+1)×m 上Hessenberg
            double[][] V = new double[mRestart + 1][];
            for (int i = 0; i <= mRestart; i++) V[i] = new double[n];
            double[,] H = new double[mRestart + 1, mRestart];
            double[] cs = new double[mRestart];
            double[] sn = new double[mRestart];
            double[] g = new double[mRestart + 1];

            int totalIter = 0;
            while (totalIter < maxIter)
            {
                // r = b - A x
                A.Mult(x, r);
                for (int i = 0; i < n; i++) r[i] = b[i] - r[i];
                // z = M^{-1} r （右预处理）；此处使用左预处理简化
                double[] z0 = new double[n];
                M.Solve(r, z0);
                double beta = 0.0;
                for (int i = 0; i < n; i++) beta += z0[i] * z0[i];
                beta = Math.Sqrt(beta);
                if (beta < tol * bnorm) { resNorm = beta / bnorm; return true; }
                for (int i = 0; i < n; i++) V[0][i] = z0[i] / beta;
                for (int i = 0; i <= mRestart; i++) g[i] = 0.0;
                g[0] = beta;

                int j;
                for (j = 0; j < mRestart && totalIter < maxIter; j++, totalIter++)
                {
                    // w = A * M^{-1} V[j]
                    double[] mv = new double[n];
                    M.Solve(V[j], mv);
                    A.Mult(mv, w);
                    // Arnoldi 正交化
                    for (int i = 0; i <= j; i++)
                    {
                        double hij = 0.0;
                        for (int k = 0; k < n; k++) hij += w[k] * V[i][k];
                        H[i, j] = hij;
                        for (int k = 0; k < n; k++) w[k] -= hij * V[i][k];
                    }
                    double hj1j = 0.0;
                    for (int k = 0; k < n; k++) hj1j += w[k] * w[k];
                    hj1j = Math.Sqrt(hj1j);
                    H[j + 1, j] = hj1j;
                    if (hj1j < 1e-30) { break; }
                    for (int k = 0; k < n; k++) V[j + 1][k] = w[k] / hj1j;

                    // 应用之前的旋转
                    for (int i = 0; i < j; i++)
                    {
                        double h1 = H[i, j], h2 = H[i + 1, j];
                        H[i, j] = cs[i] * h1 + sn[i] * h2;
                        H[i + 1, j] = -sn[i] * h1 + cs[i] * h2;
                    }
                    // 计算新旋转
                    double denom = Math.Sqrt(H[j, j] * H[j, j] + H[j + 1, j] * H[j + 1, j]);
                    if (denom < 1e-30) { cs[j] = 1; sn[j] = 0; }
                    else
                    {
                        cs[j] = H[j, j] / denom;
                        sn[j] = H[j + 1, j] / denom;
                    }
                    H[j, j] = cs[j] * H[j, j] + sn[j] * H[j + 1, j];
                    H[j + 1, j] = 0.0;
                    double g1 = g[j], g2 = g[j + 1];
                    g[j] = cs[j] * g1 + sn[j] * g2;
                    g[j + 1] = -sn[j] * g1 + cs[j] * g2;
                    resNorm = Math.Abs(g[j + 1]) / bnorm;
                    iter = totalIter + 1;
                    if (resNorm < tol)
                    {
                        // 回代
                        double[] y = new double[j + 1];
                        for (int i = j; i >= 0; i--)
                        {
                            double s = g[i];
                            for (int k = i + 1; k <= j; k++) s -= H[i, k] * y[k];
                            y[i] = s / H[i, i];
                        }
                        // x = x + M^{-1} V[:,0..j] * y
                        double[] dx = new double[n];
                        for (int i = 0; i <= j; i++)
                            for (int k = 0; k < n; k++) dx[k] += y[i] * V[i][k];
                        double[] mdx = new double[n];
                        M.Solve(dx, mdx);
                        for (int k = 0; k < n; k++) x[k] += mdx[k];
                        return true;
                    }
                }
                // 重启：回代
                int jj = j - 1;
                if (jj < 0) return false;
                double[] yy = new double[jj + 1];
                for (int i = jj; i >= 0; i--)
                {
                    double s = g[i];
                    for (int k = i + 1; k <= jj; k++) s -= H[i, k] * yy[k];
                    yy[i] = s / H[i, i];
                }
                double[] dxx = new double[n];
                for (int i = 0; i <= jj; i++)
                    for (int k = 0; k < n; k++) dxx[k] += yy[i] * V[i][k];
                double[] mdxx = new double[n];
                M.Solve(dxx, mdxx);
                for (int k = 0; k < n; k++) x[k] += mdxx[k];
            }
            return false;
        }
    }
}
