using System;
using System.Collections.Generic;
using System.Diagnostics;
using Fem2DFluid.Models;

namespace Fem2DFluid.Services
{
    /// <summary>
    /// 二维势流有限元求解器（稀疏矩阵版本）
    /// 控制方程：∇·(k ∇φ) = f  在Ω内
    /// Dirichlet: φ = φ0          在Γ_D
    /// Neumann:  k ∂φ/∂n = q      在Γ_N
    ///
    /// 线性三角单元（CST），1点高斯积分；
    /// 刚度矩阵使用 CSR 压缩稀疏行格式存储；
    /// 边界条件施加后使用 ILU0 预处理的 BiCGStab 迭代求解，可支撑大规模问题。
    /// </summary>
    public class FemSolver
    {
        /// <summary>
        /// 求解模型：组装稀疏刚度矩阵、施加边界条件、迭代求解
        /// </summary>
        public static FemResult Solve(FemModel model)
        {
            FemResult result = new FemResult();
            Stopwatch sw = Stopwatch.StartNew();
            int n = model.Nodes.Count;
            if (n == 0) { result.Success = false; result.Message = "没有节点"; return result; }

            // 建立 id -> 连续索引映射（节点ID不一定连续）
            Dictionary<int, int> idxMap = new Dictionary<int, int>();
            for (int i = 0; i < n; i++) idxMap[model.Nodes[i].Id] = i;

            // ====== 用 COO 三元组装配稀疏刚度矩阵 ======
            CooMatrix coo = new CooMatrix(n);
            double[] F = new double[n];

            for (int ei = 0; ei < model.Elements.Count; ei++)
            {
                Element e = model.Elements[ei];
                Material mat = model.GetMaterial(e.MatId) ?? model.Materials[0];
                Node n1 = model.GetNode(e.N1);
                Node n2 = model.GetNode(e.N2);
                Node n3 = model.GetNode(e.N3);
                if (n1 == null || n2 == null || n3 == null) continue;
                int i1 = idxMap[e.N1], i2 = idxMap[e.N2], i3 = idxMap[e.N3];

                double x1 = n1.X, y1 = n1.Y;
                double x2 = n2.X, y2 = n2.Y;
                double x3 = n3.X, y3 = n3.Y;
                double area = 0.5 * ((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1));
                if (Math.Abs(area) < 1e-12) continue;
                double A = Math.Abs(area);

                double b1 = y2 - y3, b2 = y3 - y1, b3 = y1 - y2;
                double c1 = x3 - x2, c2 = x1 - x3, c3 = x2 - x1;
                double coef = mat.K / (4 * A);

                int[] idx = new int[] { i1, i2, i3 };
                double[] bb = new double[] { b1, b2, b3 };
                double[] cc = new double[] { c1, c2, c3 };
                for (int r = 0; r < 3; r++)
                {
                    for (int cc2 = 0; cc2 < 3; cc2++)
                    {
                        coo.Add(idx[r], idx[cc2], coef * (bb[r] * bb[cc2] + cc[r] * cc[cc2]));
                    }
                }
                // 源项载荷：∫ f N_i dΩ = f*Δ/3
                double fSrc = mat.Source;
                double fLoad = fSrc * A / 3.0;
                F[i1] += fLoad; F[i2] += fLoad; F[i3] += fLoad;
            }

            // ====== COO -> CSR ======
            CsrMatrix K = coo.ToCsr();

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

            // 置1法：对Dirichlet行i，先把其他行的 K[j,i] 贡献扣除到F并清0，
            //       再对Dirichlet行i，令 K[i,j]=0 (j!=i), K[i,i]=1, F[i]=val
            for (int di = 0; di < dirichletNodes.Count; di++)
            {
                int i = dirichletNodes[di];
                double val = model.Nodes[i].BCValue;
                // 第一步：遍历第i列，对所有非Dirichlet行j，把K[j,i]*val加到F[j]并清K[j,i]=0
                // 利用CSR逐行查找列i（稀疏结构内线性查找，代价O(nnz/n·n)≈O(n)，小规模足够）
                for (int j = 0; j < n; j++)
                {
                    if (j == i || isD[j]) continue;
                    double v = K.Get(j, i);
                    if (Math.Abs(v) > 1e-18)
                    {
                        F[j] -= v * val;
                        K.SetZero(j, i);
                    }
                }
                // 第二步：Dirichlet行i清零非对角、对角设为1、右端设val
                K.ZeroRow(i);
                K.Set(i, i, 1.0);
                F[i] = val;
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
                            if (cnt >= 4) break;
                        }
                    }
                    double halfL = sumLen * 0.5;
                    F[i] += ni.BCValue * halfL;
                }
            }

            // ====== ILU0 预处理 + BiCGStab 求解 ======
            Ilu0 precond = new Ilu0(K);
            double[] phi = new double[n];
            int iter;
            bool conv = BiCgStab.Solve(K, F, phi, precond, 1e-10, 2000, out iter);
            if (!conv)
            {
                result.Success = false;
                result.Message = "求解失败：迭代未收敛（" + iter + "步），矩阵可能奇异或病态";
                return result;
            }
            for (int i = 0; i < n; i++) model.Nodes[i].Phi = phi[i];

            // ====== 计算单元内速度（势梯度 u=∂φ/∂x, v=∂φ/∂y） ======
            double phiMin = double.MaxValue, phiMax = double.MinValue;
            double vMin = double.MaxValue, vMax = double.MinValue;
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
            sw.Stop();
            result.Success = true;
            result.Message = "求解成功（迭代" + iter + "步，稀疏ILU0-BiCGStab）";
            result.NodeCount = n;
            result.ElementCount = model.Elements.Count;
            result.PhiMin = phiMin;
            result.PhiMax = phiMax;
            result.VmagMin = vMin;
            result.VmagMax = vMax;
            result.SolveTime = sw.Elapsed;
            return result;
        }
    }

    // ====================================================================
    // 稀疏矩阵：COO 三元组（装配阶段用）
    // ====================================================================
    internal class CooMatrix
    {
        private int _n;
        private List<int> _ri, _ci;
        private List<double> _vv;

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
            _ri.Add(i);
            _ci.Add(j);
            _vv.Add(v);
        }

        public CsrMatrix ToCsr()
        {
            int nnz = _ri.Count;
            // 统计每行非零数
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
            // 每行内按列号排序、合并相同列号
            for (int i = 0; i < _n; i++)
            {
                int p0 = rowPtr[i];
                int p1 = rowPtr[i + 1];
                if (p1 - p0 <= 1) continue;
                // 简单插入排序（每行非零数通常很小：<10）
                for (int a = p0 + 1; a < p1; a++)
                {
                    int tc = cols[a]; double tv = vals[a];
                    int b = a - 1;
                    while (b >= p0 && cols[b] > tc)
                    {
                        cols[b + 1] = cols[b];
                        vals[b + 1] = vals[b];
                        b--;
                    }
                    cols[b + 1] = tc;
                    vals[b + 1] = tv;
                }
                // 合并相同列
                int write = p0;
                for (int a = p0; a < p1; a++)
                {
                    if (write > p0 && cols[write - 1] == cols[a])
                    {
                        vals[write - 1] += vals[a];
                    }
                    else
                    {
                        cols[write] = cols[a];
                        vals[write] = vals[a];
                        write++;
                    }
                }
                // 收缩行尾（实际CSR不回收空间，将多余位置置零并调整rowPtr）
                // 为简化，重新构造rowPtr/cols/vals
                // 这里直接写回并把write之后的内容清零
                for (int a = write; a < p1; a++) { vals[a] = 0.0; cols[a] = -1; }
                rowPtr[i + 1] = write + (rowPtr[i + 1] - p1); // 收缩多余长度
            }
            // 第二次压缩：去掉cols/vals中 -1 的项
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
                    {
                        newCols[wp] = cols[k];
                        newVals[wp] = vals[k];
                        wp++;
                    }
                }
            }
            newRowPtr[_n] = wp;
            return new CsrMatrix(_n, newRowPtr, newCols, newVals);
        }
    }

    // ====================================================================
    // 稀疏矩阵：CSR 压缩稀疏行格式
    // ====================================================================
    internal class CsrMatrix
    {
        public int N;
        public int[] RowPtr;
        public int[] Cols;
        public double[] Vals;

        public CsrMatrix(int n, int[] rp, int[] c, double[] v)
        {
            N = n; RowPtr = rp; Cols = c; Vals = v;
        }

        // 定位 (i,j) 的存储位置，-1表示不存在
        public int FindPos(int i, int j)
        {
            int p0 = RowPtr[i], p1 = RowPtr[i + 1];
            // 二分查找
            int lo = p0, hi = p1 - 1;
            while (lo <= hi)
            {
                int mid = (lo + hi) >> 1;
                if (Cols[mid] == j) return mid;
                if (Cols[mid] < j) lo = mid + 1; else hi = mid - 1;
            }
            return -1;
        }

        public double Get(int i, int j)
        {
            int p = FindPos(i, j);
            return p < 0 ? 0.0 : Vals[p];
        }

        public void Set(int i, int j, double v)
        {
            int p = FindPos(i, j);
            if (p >= 0) Vals[p] = v;
            // 若不存在，需要插入；CSR插入代价大，这里仅用于Dirichlet清列，不存在则忽略
            // （因为Dirichlet列清0时目标列i的值已经被遍历行内清除，列i上的非零项是K[i,i]=1通过ZeroRow+Set写入）
        }

        /// <summary>将(i,j)位置置为0（不破坏稀疏结构，仅将值写0）</summary>
        public void SetZero(int i, int j)
        {
            int p = FindPos(i, j);
            if (p >= 0) Vals[p] = 0.0;
        }

        public void ZeroRow(int i)
        {
            int p0 = RowPtr[i], p1 = RowPtr[i + 1];
            for (int k = p0; k < p1; k++) Vals[k] = 0.0;
        }

        // y = A*x
        public void Mult(double[] x, double[] y)
        {
            for (int i = 0; i < N; i++)
            {
                double s = 0.0;
                for (int k = RowPtr[i]; k < RowPtr[i + 1]; k++)
                {
                    s += Vals[k] * x[Cols[k]];
                }
                y[i] = s;
            }
        }
    }

    // ====================================================================
    // ILU0 不完全LU预处理（基于CSR稀疏结构）——标准Saad算法
    // 参考：Yousef Saad, "Iterative Methods for Sparse Linear Systems", 算法10.1
    // 在原矩阵A的稀疏结构内分解：A ≈ LU
    // L 为单位下三角（对角元为1，存于 lu 的严格下三角部分）
    // U 为上三角（对角元存于 lu 的对角位置）
    // ====================================================================
    internal class Ilu0
    {
        private int _n;
        private int[] _rowPtr, _cols;
        private double[] _lu;       // L/U重叠存储
        private int[] _diagPos;    // 每行对角元在 lu 中的位置

        public Ilu0(CsrMatrix A)
        {
            _n = A.N;
            // 拷贝A的稀疏结构与数值
            _rowPtr = (int[])A.RowPtr.Clone();
            _cols = (int[])A.Cols.Clone();
            _lu = (double[])A.Vals.Clone();
            _diagPos = new int[_n];
            // 找到每行对角元的位置
            for (int i = 0; i < _n; i++)
            {
                _diagPos[i] = -1;
                for (int k = _rowPtr[i]; k < _rowPtr[i + 1]; k++)
                {
                    if (_cols[k] == i) { _diagPos[i] = k; break; }
                }
            }
            // ===== 标准 ILU0 分解（Saad 算法 10.1）=====
            // for i = 1..n-1
            //   for k = 1..i-1 and (i,k) in Sparsity(A)
            //     aik = aik / akk
            //     for j = k+1..n-1 and (i,j) in Sparsity(A)
            //       if (k,j) in Sparsity(A) then aij = aij - aik * akj
            for (int i = 0; i < _n; i++)
            {
                int d = _diagPos[i];
                if (d < 0 || Math.Abs(_lu[d]) < 1e-30) continue;
                // 遍历i行中的严格下三角元素 k < i（在lu中位置 < d）
                for (int kk = _rowPtr[i]; kk < d; kk++)
                {
                    int k = _cols[kk];
                    if (k < 0 || k >= i) continue;
                    double lik = _lu[kk] / _lu[_diagPos[k]];   // aik = aik / akk
                    _lu[kk] = lik;
                    int kd = _diagPos[k];
                    // 对k行中列号 j > k 的元素（上三角部分），若(i,j)存在则更新 aij -= lik * akj
                    for (int kj = kd + 1; kj < _rowPtr[k + 1]; kj++)
                    {
                        int j = _cols[kj];
                        int ip = FindInRow(i, j);
                        if (ip >= 0)
                        {
                            _lu[ip] -= lik * _lu[kj];
                        }
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

        // z = M^{-1} r，即求解 LU z = r（前向/回代）
        public void Solve(double[] r, double[] z)
        {
            // L 是单位下三角（对角1），U 是上三角
            // 前向 L y = r
            double[] y = new double[_n];
            for (int i = 0; i < _n; i++)
            {
                double s = r[i];
                for (int k = _rowPtr[i]; k < _diagPos[i]; k++)
                {
                    s -= _lu[k] * y[_cols[k]];
                }
                y[i] = s;
            }
            // 回代 U z = y
            for (int i = _n - 1; i >= 0; i--)
            {
                double s = y[i];
                for (int k = _diagPos[i] + 1; k < _rowPtr[i + 1]; k++)
                {
                    s -= _lu[k] * z[_cols[k]];
                }
                z[i] = s / _lu[_diagPos[i]];
            }
        }
    }

    // ====================================================================
    // BiCGStab（稳定双共轭梯度法）带左预处理
    // ====================================================================
    internal static class BiCgStab
    {
        public static bool Solve(CsrMatrix A, double[] b, double[] x, Ilu0 M, double tol, int maxIter, out int iter)
        {
            int n = A.N;
            double[] r = new double[n];
            double[] r0 = new double[n];
            double[] p = new double[n];
            double[] Ap = new double[n];
            double[] s = new double[n];
            double[] As = new double[n];
            double[] Mp = new double[n];
            double[] Ms = new double[n];

            // r = b - A x
            A.Mult(x, r);
            for (int i = 0; i < n; i++) r[i] = b[i] - r[i];
            for (int i = 0; i < n; i++) r0[i] = r[i];

            double rho = 1.0, alpha = 1.0, omega = 1.0;
            double bnorm = 0.0;
            for (int i = 0; i < n; i++) bnorm += b[i] * b[i];
            bnorm = Math.Sqrt(bnorm);
            if (bnorm < 1e-30) bnorm = 1.0;

            iter = 0;
            for (;;)
            {
                double rhoNew = 0.0;
                for (int i = 0; i < n; i++) rhoNew += r0[i] * r[i];
                if (Math.Abs(rhoNew) < 1e-30) { return false; }
                double beta = (rhoNew / rho) * (alpha / omega);
                rho = rhoNew;
                // p = r + beta*(p - omega*Ap)
                for (int i = 0; i < n; i++) p[i] = r[i] + beta * (p[i] - omega * Ap[i]);
                // Mp = M^{-1} p
                M.Solve(p, Mp);
                // Ap = A*Mp
                A.Mult(Mp, Ap);
                // alpha = rho / (r0·Ap)
                double r0Ap = 0.0;
                for (int i = 0; i < n; i++) r0Ap += r0[i] * Ap[i];
                if (Math.Abs(r0Ap) < 1e-30) { return false; }
                alpha = rho / r0Ap;
                // s = r - alpha*Ap
                for (int i = 0; i < n; i++) s[i] = r[i] - alpha * Ap[i];
                // 收敛检查：|s| < tol*|b|
                double snorm = 0.0;
                for (int i = 0; i < n; i++) snorm += s[i] * s[i];
                snorm = Math.Sqrt(snorm);
                if (snorm < tol * bnorm)
                {
                    for (int i = 0; i < n; i++) x[i] += alpha * Mp[i];
                    iter++;
                    return true;
                }
                // Ms = M^{-1} s
                M.Solve(s, Ms);
                // As = A*Ms
                A.Mult(Ms, As);
                // omega = (As·s)/(As·As)
                double AsS = 0.0, AsAs = 0.0;
                for (int i = 0; i < n; i++) { AsS += As[i] * s[i]; AsAs += As[i] * As[i]; }
                if (AsAs < 1e-30) { return false; }
                omega = AsS / AsAs;
                // x = x + alpha*Mp + omega*Ms
                for (int i = 0; i < n; i++) x[i] += alpha * Mp[i] + omega * Ms[i];
                // r = s - omega*As
                for (int i = 0; i < n; i++) r[i] = s[i] - omega * As[i];
                double rnorm = 0.0;
                for (int i = 0; i < n; i++) rnorm += r[i] * r[i];
                rnorm = Math.Sqrt(rnorm);
                iter++;
                if (rnorm < tol * bnorm) return true;
                if (iter >= maxIter) return false;
                if (Math.Abs(omega) < 1e-30) return false;
            }
        }
    }
}
