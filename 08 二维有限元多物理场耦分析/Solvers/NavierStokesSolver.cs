using System;
using System.Collections;
using System.Collections.Generic;
using System.Threading;
using MultiPhysicsFEM2D.Core;
using MultiPhysicsFEM2D.Elements;
using MultiPhysicsFEM2D.Mesh;

namespace MultiPhysicsFEM2D.Solvers
{
    public class NSResult
    {
        public int OuterIterations;
        public double FinalResU, FinalResP, FinalResT;
        public bool Transient;
    }

    public class NavierStokesSolver
    {
        public bool Transient = false;
        public double TimeStep = 0.01;
        public int TimeSteps = 20;
        public int OuterIterations = 200;
        public double AlphaU = 0.7;
        public double AlphaP = 0.3;
        public double Tol = 1e-4;
        public double Gravity = 9.81;
        public double BetaExpansion = 0.00034;
        public double TRef = 0.0;
        public bool IncludeThermal = false; // Boussinesq 耦合
        public Action<string> Log;

        private CooAssembler AsmU, AsmT;
        private SparseMatrixCSR _M, _Hdiag, _AinvDiag;
        private DenseVector _rhs, _pRhs, _rhsT;
        private int _nu, _np, _nT;
        private int _nNodes;
        private FEMesh _mesh;

        public NSResult Solve(FEMesh mesh)
        {
            _mesh = mesh;
            _nNodes = mesh.NumNodes;
            _nu = 2 * _nNodes; _np = _nNodes; _nT = _nNodes;
            if (Transient)
            {
                // 时间循环：每步调稳态 SIMPLE
                bool oldT = Transient;
                Transient = false;
                // 初始 u=p=0，T=Tref
                Array.Clear(mesh.VelocityU, 0, _nNodes);
                Array.Clear(mesh.VelocityV, 0, _nNodes);
                Array.Clear(mesh.Pressure, 0, _nNodes);
                if (IncludeThermal)
                    for (int i = 0; i < _nNodes; i++) mesh.Temperature[i] = TRef;
                NSResult last = null;
                for (int s = 0; s < TimeSteps; s++)
                {
                    if (Log != null) Log(string.Format("NS time step {0}/{1} dt={2}", s + 1, TimeSteps, TimeStep));
                    last = RunSIMPLE(mesh);
                }
                Transient = oldT;
                return last;
            }
            else return RunSIMPLE(mesh);
        }

        private NSResult RunSIMPLE(FEMesh mesh)
        {
            _nNodes = mesh.NumNodes;
            _nu = 2 * _nNodes; _np = _nNodes; _nT = _nNodes;
            double[] uOld = new double[_nNodes], vOld = new double[_nNodes];
            Array.Copy(mesh.VelocityU, uOld, _nNodes); Array.Copy(mesh.VelocityV, vOld, _nNodes);
            int maxOut = OuterIterations;
            double resU = 1, resP = 1, resT = 1;
            for (int it = 1; it <= maxOut; it++)
            {
                BuildMomentumSystem(mesh, out SparseMatrixCSR H, out DenseVector b);
                // Jacobi 预条件的简单近似 A = diag(H)
                DenseVector diagInv = new DenseVector(_nu);
                for (int i = 0; i < _nu; i++)
                {
                    double d = 1.0;
                    for (int k = H.RowPtr[i]; k < H.RowPtr[i + 1]; k++)
                        if (H.ColIdx[k] == i) { d = H.Values[k]; break; }
                    diagInv.Values[i] = Math.Abs(d) > 1e-30 ? 1.0 / d : 1.0;
                }
                DenseVector uStar = new DenseVector(_nu);
                for (int i = 0; i < _nNodes; i++) { uStar.Values[2 * i] = mesh.VelocityU[i]; uStar.Values[2 * i + 1] = mesh.VelocityV[i]; }
                JacobiSmooth(H, b, uStar, diagInv, 30); // 近似内解
                // 写 u*
                double[] uStarX = new double[_nNodes], uStarY = new double[_nNodes];
                for (int i = 0; i < _nNodes; i++) { uStarX[i] = uStar.Values[2 * i]; uStarY[i] = uStar.Values[2 * i + 1]; }
                // 压力泊松：Ap p' = -div(u*) · (diag(H)^{-1})
                SparseMatrixCSR Ap;
                DenseVector bp;
                BuildPressurePoisson(mesh, diagInv, uStarX, uStarY, out Ap, out bp);
                DenseVector pPrime = new DenseVector(_np);
                ApplyDirichletScalar(Ap, bp, BCType.Pressure, mesh, 0.0);
                JacobiSmoothScalar(Ap, bp, pPrime, 60);
                // 压力修正 + 速度修正
                for (int i = 0; i < _nNodes; i++) mesh.Pressure[i] += AlphaP * pPrime.Values[i];
                // 速度修正 u' = -Hdiag^-1 G p'  （逐点近似）
                DenseVector dpdx = new DenseVector(_nNodes), dpdy = new DenseVector(_nNodes);
                ComputeNodalGradient(mesh, pPrime.Values, dpdx.Values, dpdy.Values);
                for (int i = 0; i < _nNodes; i++)
                {
                    mesh.VelocityU[i] = uStarX[i] - AlphaU * diagInv.Values[2 * i] * dpdx.Values[i];
                    mesh.VelocityV[i] = uStarY[i] - AlphaU * diagInv.Values[2 * i + 1] * dpdy.Values[i];
                }
                ApplyVelocityDirichlet(mesh);
                // 残差
                double du = 0, duT = 0, dp = 0;
                for (int i = 0; i < _nNodes; i++)
                {
                    du += (mesh.VelocityU[i] - uOld[i]) * (mesh.VelocityU[i] - uOld[i])
                        + (mesh.VelocityV[i] - vOld[i]) * (mesh.VelocityV[i] - vOld[i]);
                    duT += mesh.VelocityU[i] * mesh.VelocityU[i] + mesh.VelocityV[i] * mesh.VelocityV[i];
                }
                // 热方程（Boussinesq）
                if (IncludeThermal)
                {
                    SolveTemperature(mesh);
                }
                resU = Math.Sqrt(du) / Math.Sqrt(Math.Max(duT, 1e-30));
                if (Log != null && (it % 10 == 0 || it < 5 || resU < Tol))
                    Log(string.Format("  SIMPLE outer {0}: resU = {1:E3}", it, resU));
                Array.Copy(mesh.VelocityU, uOld, _nNodes); Array.Copy(mesh.VelocityV, vOld, _nNodes);
                if (resU < Tol)
                    return new NSResult { OuterIterations = it, FinalResU = resU, FinalResP = resP, FinalResT = resT, Transient = Transient };
            }
            return new NSResult { OuterIterations = maxOut, FinalResU = resU, FinalResP = resP, FinalResT = resT, Transient = Transient };
        }

        private void BuildMomentumSystem(FEMesh mesh, out SparseMatrixCSR H, out DenseVector b)
        {
            CooAssembler coo = new CooAssembler(_nu, _nu);
            DenseVector rhs = new DenseVector(_nu);
            int nE = mesh.Elements.Count;
            int threads = Environment.ProcessorCount;
            int block = (nE + threads - 1) / threads;
            WaitHandle[] whs = new WaitHandle[threads];
            List<Tuple<int[], double[,]>>[] perB = new List<Tuple<int[], double[,]>>[threads];
            double[][] perF = new double[threads][];
            for (int t = 0; t < threads; t++) { perB[t] = new List<Tuple<int[], double[,]>>(); perF[t] = new double[_nu]; }
            for (int t = 0; t < threads; t++)
            {
                int s = t * block, e = Math.Min(s + block, nE);
                int tt = t;
                ManualResetEvent mre = new ManualResetEvent(false); whs[t] = mre;
                ThreadPool.QueueUserWorkItem(delegate
                {
                    for (int ei = s; ei < e; ei++)
                    {
                        FiniteElement el = mesh.Elements[ei];
                        int nn = el.NodesPerElement;
                        int ldof = 2 * nn;
                        double[,] Hloc = new double[ldof, ldof];
                        double[,] Visc = new double[ldof, ldof];
                        double[,] Conv = new double[ldof, ldof];
                        el.BuildViscousMatrix(Visc);
                        double[] uN = new double[nn], vN = new double[nn];
                        for (int k = 0; k < nn; k++) { uN[k] = mesh.VelocityU[el.NodeIds[k]]; vN[k] = mesh.VelocityV[el.NodeIds[k]]; }
                        el.BuildConvectionMatrix(uN, vN, Conv);
                        double[,] Mass = new double[ldof, ldof];
                        el.BuildMassMatrix(Mass, true);
                        for (int i = 0; i < ldof; i++) for (int j = 0; j < ldof; j++)
                        {
                            double v = Visc[i, j] + Conv[i, j];
                            if (Transient) v += Mass[i, j] / TimeStep;
                            Hloc[i, j] = v;
                        }
                        // 压力梯度项加入 rhs: -∫ ∇p·N → -Gx·p, -Gy·p
                        int pn = el.PressureNodesPerElement;
                        double[,] Gx = new double[pn, ldof];
                        double[,] Gy = new double[pn, ldof];
                        el.BuildGradientMatrix(Gx, Gy);
                        double[] pN = new double[nn];
                        for (int k = 0; k < nn; k++) pN[k] = mesh.Pressure[el.NodeIds[k]];
                        // Gx^T * pN  (rhs_u), Gy^T * pN (rhs_v)
                        int[] dofMap = new int[ldof];
                        for (int k = 0; k < nn; k++) { dofMap[2 * k] = 2 * el.NodeIds[k]; dofMap[2 * k + 1] = 2 * el.NodeIds[k] + 1; }
                        for (int i = 0; i < ldof; i++)
                        {
                            double fu = 0, fv = 0;
                            for (int q = 0; q < pn; q++)
                            {
                                int pi = (q < nn) ? el.NodeIds[q] : el.NodeIds[0]; // Q4/LT6 pNodeMap == same
                                fu -= Gx[q, i] * pN[q];
                                fv -= Gy[q, i] * pN[q];
                            }
                            // Boussinesq y 向浮力
                            if (IncludeThermal)
                            {
                                double[] Tloc = new double[nn];
                                for (int k = 0; k < nn; k++) Tloc[k] = mesh.Temperature[el.NodeIds[k]];
                                // ∫ N_i ρ g β (T-Tref) dy
                                // 用形函数近似，平均T
                                double Tav = 0;
                                for (int k = 0; k < nn; k++) Tav += Tloc[k] / nn;
                                double area = el.ComputeArea(), th = el.Thickness;
                                double Ni = 1.0 / nn;
                                if (i % 2 == 1) fv += -el.Density * Gravity * BetaExpansion * (Tav - TRef) * area * th * Ni;
                            }
                            perF[tt][dofMap[i]] += fu + (i % 2 == 1 ? fv : 0); // careful: fu 是 x 分量
                            if (i % 2 == 0) perF[tt][dofMap[i]] += fu;
                            else perF[tt][dofMap[i]] += fv;
                        }
                        // 瞬态质量项旧值
                        if (Transient)
                        {
                            for (int i = 0; i < nn; i++)
                            {
                                double mu = Mass[2 * i, 2 * i] / TimeStep * uN[i];
                                double mv = Mass[2 * i + 1, 2 * i + 1] / TimeStep * vN[i];
                                perF[tt][2 * el.NodeIds[i]] += mu;
                                perF[tt][2 * el.NodeIds[i] + 1] += mv;
                            }
                        }
                        perB[tt].Add(Tuple.Create(dofMap, Hloc));
                    }
                    mre.Set();
                });
            }
            WaitHandle.WaitAll(whs);
            for (int t = 0; t < threads; t++) { coo.AddLocalBulk(perB[t]); for (int i = 0; i < _nu; i++) rhs.Values[i] += perF[t][i]; }
            // Dirichlet 速度
            H = coo.BuildCSR();
            b = rhs;
            ApplyDirichletVelocity(H, b, mesh);
        }

        private void ApplyDirichletVelocity(SparseMatrixCSR H, DenseVector b, FEMesh mesh)
        {
            foreach (BoundaryCondition bc in mesh.BCs)
            {
                if (bc.Type == BCType.VelocityX) H.ApplyDirichlet(2 * bc.NodeId, 1e16, bc.Value, b);
                else if (bc.Type == BCType.VelocityY) H.ApplyDirichlet(2 * bc.NodeId + 1, 1e16, bc.Value, b);
            }
        }

        private void ApplyVelocityDirichlet(FEMesh mesh)
        {
            foreach (BoundaryCondition bc in mesh.BCs)
            {
                if (bc.Type == BCType.VelocityX) mesh.VelocityU[bc.NodeId] = bc.Value;
                else if (bc.Type == BCType.VelocityY) mesh.VelocityV[bc.NodeId] = bc.Value;
            }
        }

        private void ApplyDirichletScalar(SparseMatrixCSR A, DenseVector b, BCType bct, FEMesh mesh, double penalty)
        {
            foreach (BoundaryCondition bc in mesh.BCs)
                if (bc.Type == bct) A.ApplyDirichlet(bc.NodeId, penalty, bc.Value, b);
        }

        private void BuildPressurePoisson(FEMesh mesh, DenseVector diagInvU, double[] uStarX, double[] uStarY,
            out SparseMatrixCSR Ap, out DenseVector bp)
        {
            // Ap ≈ G · diag(H)^{-1} · G^T   逐点近似：Ap ~ -Laplacian（用 Laplacian 近似）
            CooAssembler coo = new CooAssembler(_np, _np);
            bp = new DenseVector(_np);
            int nE = mesh.Elements.Count;
            int threads = Environment.ProcessorCount;
            int block = (nE + threads - 1) / threads;
            WaitHandle[] whs = new WaitHandle[threads];
            List<Tuple<int[], double[,]>>[] perB = new List<Tuple<int[], double[,]>>[threads];
            double[][] perF = new double[threads][];
            for (int t = 0; t < threads; t++) { perB[t] = new List<Tuple<int[], double[,]>>(); perF[t] = new double[_np]; }
            for (int t = 0; t < threads; t++)
            {
                int s = t * block, e = Math.Min(s + block, nE);
                int tt = t;
                ManualResetEvent mre = new ManualResetEvent(false); whs[t] = mre;
                ThreadPool.QueueUserWorkItem(delegate
                {
                    for (int ei = s; ei < e; ei++)
                    {
                        FiniteElement el = mesh.Elements[ei];
                        int nn = el.NodesPerElement;
                        // 用单元压力质量/Laplacian 组装 Ap：简单取 Laplacian 近似
                        // Kp ≈ ∫ ∇N_p · ∇N_p dΩ  （压力 Poisson 近似）
                        double[,] dummy = new double[nn, nn];
                        // 借助热传导矩阵构造 ∇N·∇N
                        double oldK = el.Conductivity; el.Conductivity = 1.0;
                        double[,] Kt = new double[nn, nn];
                        el.BuildConductivityMatrix(Kt);
                        el.Conductivity = oldK;
                        // Mp
                        double[,] Mp = new double[nn, nn];
                        el.BuildPressureMassMatrix(Mp);
                        // Ap += (Kt + 稳定化系数·Mp)
                        double[,] Aploc = new double[nn, nn];
                        double stab = 0.01; // 稳定化
                        for (int i = 0; i < nn; i++)
                            for (int j = 0; j < nn; j++)
                                Aploc[i, j] = Kt[i, j] + stab * Mp[i, j] / Math.Max(el.ComputeArea(), 1e-10);
                        int[] dofMap = new int[nn];
                        for (int k = 0; k < nn; k++) dofMap[k] = el.NodeIds[k];
                        perB[tt].Add(Tuple.Create(dofMap, Aploc));
                        // bp = -div(u*) —— 用 ∫ N div(u*) ≈ -∫ ∇N · u* dΩ 分部积分
                        double[,] Gx = new double[nn, 2 * nn];
                        double[,] Gy = new double[nn, 2 * nn];
                        el.BuildGradientMatrix(Gx, Gy);
                        for (int i = 0; i < nn; i++)
                        {
                            double div = 0;
                            for (int j = 0; j < nn; j++)
                            {
                                div -= Gx[i, 2 * j] * uStarX[el.NodeIds[j]];
                                div -= Gy[i, 2 * j + 1] * uStarY[el.NodeIds[j]];
                            }
                            perF[tt][el.NodeIds[i]] += div;
                        }
                    }
                    mre.Set();
                });
            }
            WaitHandle.WaitAll(whs);
            for (int t = 0; t < threads; t++) { coo.AddLocalBulk(perB[t]); for (int i = 0; i < _np; i++) bp.Values[i] += perF[t][i]; }
            Ap = coo.BuildCSR();
            // 压力参考点：在 Dirichlet pressure 处已由 ApplyDirichlet 处理
        }

        private void ComputeNodalGradient(FEMesh mesh, double[] f, double[] dfdx, double[] dfdy)
        {
            Array.Clear(dfdx, 0, _nNodes); Array.Clear(dfdy, 0, _nNodes);
            double[] wsum = new double[_nNodes];
            foreach (FiniteElement el in mesh.Elements)
            {
                int nn = el.NodesPerElement;
                // 计算单元中心梯度
                double[] fLoc = new double[nn];
                for (int k = 0; k < nn; k++) fLoc[k] = f[el.NodeIds[k]];
                double gx = 0, gy = 0;
                // 用简单方法：循环 Gauss 点求平均
                if (el is CST3Element)
                {
                    // 直接 B 矩阵梯度（常数）
                    double x0 = el.X[0], y0 = el.Y[0], x1 = el.X[1], y1 = el.Y[1], x2 = el.X[2], y2 = el.Y[2];
                    double b0 = y1 - y2, b1 = y2 - y0, b2 = y0 - y1;
                    double c0 = x2 - x1, c1 = x0 - x2, c2 = x1 - x0;
                    double area = el.ComputeArea();
                    double inv2A = 1.0 / (2 * area);
                    gx = inv2A * (b0 * fLoc[0] + b1 * fLoc[1] + b2 * fLoc[2]);
                    gy = inv2A * (c0 * fLoc[0] + c1 * fLoc[1] + c2 * fLoc[2]);
                    double w = area * el.Thickness / nn;
                    for (int k = 0; k < nn; k++)
                    {
                        dfdx[el.NodeIds[k]] += w * gx;
                        dfdy[el.NodeIds[k]] += w * gy;
                        wsum[el.NodeIds[k]] += w;
                    }
                }
                else
                {
                    // 用形函数导数数值积分后按面积平均到节点
                    double area = el.ComputeArea();
                    // 对 Q4/LT6 使用形心平均
                    // 简化：在各单元形心位置用 2×2/3 点积分
                    double aW = area * el.Thickness;
                    // 调用内部方法：临时用 BuildConductivity 方式获得 ∫∇N dΩ 难以；改以单元常数近似
                    // 简单用 ShapeAt 0,0 求导数
                    double[] Nc = el.ShapeAt(0, 0);
                    double[,] dNloc = el.ShapeDerivLocal(0, 0);
                    // 需要 Jacobi; 为简洁，按面积平均分配梯度到节点，用最粗的有限差
                    for (int k = 0; k < nn; k++)
                    {
                        dfdx[el.NodeIds[k]] += 0;
                        dfdy[el.NodeIds[k]] += 0;
                        wsum[el.NodeIds[k]] += aW / nn;
                    }
                    // 使用 Green-Gauss: ∫∂f/∂x dΩ = ∮ f n_x dS ≈ ∑e f_e * n_x * l_e/2
                    // 近似：用邻接边平均。由于这是简化实现，改用单元中心周围节点有限差分（为避免过复杂，用简单两点差分）
                    // 直接采用 ShapeDerivLocal+Jacobian 通用代码：
                    ComputeElementAverageGradient(el, fLoc, out gx, out gy);
                    double w = aW / nn;
                    for (int k = 0; k < nn; k++)
                    {
                        dfdx[el.NodeIds[k]] += w * gx;
                        dfdy[el.NodeIds[k]] += w * gy;
                        wsum[el.NodeIds[k]] += w - aW / nn; // 抵消之前的零加
                    }
                }
            }
            for (int i = 0; i < _nNodes; i++)
                if (wsum[i] > 1e-30) { dfdx[i] /= wsum[i]; dfdy[i] /= wsum[i]; }
        }

        private void ComputeElementAverageGradient(FiniteElement el, double[] fLoc, out double gx, out double gy)
        {
            // 简单实现：用中心差分式（对角）
            int nn = el.NodesPerElement;
            gx = 0; gy = 0;
            // 调用单元的形函数导数（仅适合 CST3 精度）；对 Q4 用 2×2 点平均
            if (el is Q4Element)
            {
                // 2×2 gauss: 直接重构 dNdx,dNdy
                // 简化处理：复制 CST 处理过简单有限差分
                double fx0 = 0, fx1 = 0, fy0 = 0, fy1 = 0;
                for (int k = 0; k < nn; k++) { fx0 += el.X[k] * fLoc[k]; fx1 += el.X[k]; fy0 += el.Y[k] * fLoc[k]; fy1 += el.Y[k]; }
                double xAvg = fx1 / nn, yAvg = fy1 / nn, fAvg = 0;
                for (int k = 0; k < nn; k++) fAvg += fLoc[k] / nn;
                double numx = 0, denx = 0;
                for (int k = 0; k < nn; k++) { numx += (el.X[k] - xAvg) * (fLoc[k] - fAvg); denx += (el.X[k] - xAvg) * (el.X[k] - xAvg); }
                gx = denx > 1e-30 ? numx / denx : 0;
                double numy = 0, deny = 0;
                for (int k = 0; k < nn; k++) { numy += (el.Y[k] - yAvg) * (fLoc[k] - fAvg); deny += (el.Y[k] - yAvg) * (el.Y[k] - yAvg); }
                gy = deny > 1e-30 ? numy / deny : 0;
            }
            else
            {
                double x0 = el.X[0], y0 = el.Y[0], x1 = el.X[1], y1 = el.Y[1], x2 = el.X[2], y2 = el.Y[2];
                double b0 = y1 - y2, b1 = y2 - y0, b2 = y0 - y1;
                double c0 = x2 - x1, c1 = x0 - x2, c2 = x1 - x0;
                double area = Math.Abs(b0 * c1 - b1 * c0) * 0.5;
                if (area < 1e-30) { gx = 0; gy = 0; return; }
                double inv2A = 1.0 / (2 * area);
                gx = inv2A * (b0 * fLoc[0] + b1 * fLoc[1] + b2 * fLoc[2]);
                gy = inv2A * (c0 * fLoc[0] + c1 * fLoc[1] + c2 * fLoc[2]);
            }
        }

        private void JacobiSmooth(SparseMatrixCSR A, DenseVector b, DenseVector x, DenseVector diagInv, int iters)
        {
            int n = A.Nrows;
            DenseVector Ax = new DenseVector(n);
            for (int it = 0; it < iters; it++)
            {
                A.Multiply(x, Ax);
                for (int i = 0; i < n; i++)
                {
                    double r = b.Values[i] - Ax.Values[i];
                    x.Values[i] += 0.7 * diagInv.Values[i] * r;
                }
            }
        }

        private void JacobiSmoothScalar(SparseMatrixCSR A, DenseVector b, DenseVector x, int iters)
        {
            int n = A.Nrows;
            double[] di = new double[n];
            for (int i = 0; i < n; i++)
            {
                double d = 1.0;
                for (int k = A.RowPtr[i]; k < A.RowPtr[i + 1]; k++)
                    if (A.ColIdx[k] == i) { d = A.Values[k]; break; }
                di[i] = Math.Abs(d) > 1e-30 ? 1.0 / d : 1.0;
            }
            DenseVector Ax = new DenseVector(n);
            for (int it = 0; it < iters; it++)
            {
                A.Multiply(x, Ax);
                for (int i = 0; i < n; i++)
                {
                    double r = b.Values[i] - Ax.Values[i];
                    x.Values[i] += 0.7 * di[i] * r;
                }
            }
        }

        private void SolveTemperature(FEMesh mesh)
        {
            int ndof = _nNodes;
            CooAssembler coo = new CooAssembler(ndof, ndof);
            DenseVector rhs = new DenseVector(ndof);
            int nE = mesh.Elements.Count;
            int threads = Environment.ProcessorCount;
            int block = (nE + threads - 1) / threads;
            WaitHandle[] whs = new WaitHandle[threads];
            List<Tuple<int[], double[,]>>[] perB = new List<Tuple<int[], double[,]>>[threads];
            double[][] perF = new double[threads][];
            for (int t = 0; t < threads; t++) { perB[t] = new List<Tuple<int[], double[,]>>(); perF[t] = new double[ndof]; }
            for (int t = 0; t < threads; t++)
            {
                int s = t * block, e = Math.Min(s + block, nE);
                int tt = t;
                ManualResetEvent mre = new ManualResetEvent(false); whs[t] = mre;
                ThreadPool.QueueUserWorkItem(delegate
                {
                    for (int ei = s; ei < e; ei++)
                    {
                        FiniteElement el = mesh.Elements[ei];
                        int nn = el.NodesPerElement;
                        double[,] Kt = new double[nn, nn];
                        double[,] Ct = new double[nn, nn];
                        double[,] Ec = new double[nn, nn];
                        el.BuildConductivityMatrix(Kt);
                        el.BuildCapacityMatrix(Ct, true);
                        double[] uN = new double[nn], vN = new double[nn];
                        for (int k = 0; k < nn; k++) { uN[k] = mesh.VelocityU[el.NodeIds[k]]; vN[k] = mesh.VelocityV[el.NodeIds[k]]; }
                        el.BuildEnergyConvection(uN, vN, Ec);
                        double[,] A = new double[nn, nn];
                        for (int i = 0; i < nn; i++) for (int j = 0; j < nn; j++) A[i, j] = Kt[i, j] + Ec[i, j];
                        if (Transient)
                            for (int i = 0; i < nn; i++) A[i, i] += Ct[i, i] / TimeStep;
                        int[] dofMap = new int[nn];
                        for (int k = 0; k < nn; k++) dofMap[k] = el.NodeIds[k];
                        perB[tt].Add(Tuple.Create(dofMap, A));
                        if (Transient)
                        {
                            for (int i = 0; i < nn; i++) perF[tt][el.NodeIds[i]] += Ct[i, i] / TimeStep * mesh.Temperature[el.NodeIds[i]];
                        }
                    }
                    mre.Set();
                });
            }
            WaitHandle.WaitAll(whs);
            for (int t = 0; t < threads; t++) { coo.AddLocalBulk(perB[t]); for (int i = 0; i < ndof; i++) rhs.Values[i] += perF[t][i]; }
            SparseMatrixCSR A = coo.BuildCSR();
            DenseVector T = new DenseVector(ndof);
            for (int i = 0; i < ndof; i++) T.Values[i] = mesh.Temperature[i];
            ApplyDirichletScalar(A, rhs, BCType.Temperature, mesh, 1e16);
            JacobiSmoothScalar(A, rhs, T, 20);
            for (int i = 0; i < ndof; i++) mesh.Temperature[i] = T.Values[i];
        }
    }
}
