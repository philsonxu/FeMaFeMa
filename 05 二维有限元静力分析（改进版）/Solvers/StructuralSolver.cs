// StructuralSolver.cs - 二维线弹性静力求解器（平面应变）
// 多线程装配 + Krylov 求解（CG / GMRES / BiCGSTAB）+ Jacobi / AMG 预条件
// 结果直接写入 mesh 的 DisplacementU/V、StressXX/YY/XY/VonMises（单元级）
using System;
using System.Collections.Generic;
using FEM2D.Core;
using FEM2D.Elements;
using FEM2D.Mesh;

namespace FEM2D.Solvers
{
    public enum KrylovSolverKind { CG, GMRES, BiCGSTAB }
    public enum PreconditionerKind { None, Jacobi, AMG }

    public sealed class StructuralSolver
    {
        public event Action<string> OnLog;
        public event Action<int, double> OnIteration;
        public KrylovSolverKind SolverKind = KrylovSolverKind.CG;
        public PreconditionerKind PreconditionerKind = PreconditionerKind.Jacobi;
        public double Tolerance = 1.0e-8;
        public int MaxIterations = 5000;

        public void Solve(FEMesh mesh)
        {
            int ndof = 2 * mesh.NumNodes;
            Log(string.Format("[Struct] 装配整体刚度矩阵, ndof={0}, elements={1}", ndof, mesh.NumElements));

            CooAssembler coo = new CooAssembler();
            DenseVector rhs = new DenseVector(ndof);

            ParallelTool.For(0, mesh.NumElements, ie =>
            {
                FiniteElement e = mesh.Elements[ie];
                double[,] ke = e.BuildStiffnessMatrix();
                int nen = e.NodesPerElement;
                int[] dofMap = new int[2 * nen];
                for (int a = 0; a < nen; a++)
                {
                    dofMap[2 * a] = 2 * e.NodeIds[a];
                    dofMap[2 * a + 1] = 2 * e.NodeIds[a] + 1;
                }
                coo.AddBlock(dofMap, ke);
            });

            foreach (Tuple<int, int, double> nf in mesh.NodalForces)
            {
                int dof = 2 * nf.Item1 + nf.Item2;
                rhs.Values[dof] += nf.Item3;
            }

            SparseMatrixCSR K = coo.ToCSR(ndof);

            double[] diag = K.GetDiagonals();
            double maxDiag = 0.0;
            for (int i = 0; i < ndof; i++) { double d = Math.Abs(diag[i]); if (d > maxDiag) maxDiag = d; }
            double penalty = maxDiag * 1.0e12;
            foreach (Tuple<int, int, double> bc in mesh.DirichletBCs)
            {
                int dof = 2 * bc.Item1 + bc.Item2;
                K.ApplyDirichlet(dof, bc.Item3, rhs, penalty);
            }

            // 选择预条件子
            IPreconditioner M;
            string precondName;
            if (PreconditionerKind == PreconditionerKind.AMG)
            {
                Log("[Struct] 构造 AMG (Smoothed Aggregation) 预条件 ...");
                DateTime tA = DateTime.Now;
                M = new AMGPreconditioner(K);
                Log(string.Format("[Struct] AMG 构造完成, 用时 {0:F2}s", (DateTime.Now - tA).TotalSeconds));
                precondName = "AMG";
            }
            else if (PreconditionerKind == PreconditionerKind.Jacobi)
            { M = new JacobiPreconditioner(K); precondName = "Jacobi"; }
            else { M = new IdentityPreconditioner(); precondName = "None"; }

            string solverName = SolverKind.ToString();
            Log(string.Format("[Struct] 求解：{0} + {1}", solverName, precondName));

            DenseVector u = new DenseVector(ndof);
            bool conv = false; int iters = 0; double res = 0.0;

            DateTime t0 = DateTime.Now;
            Action<int, double> progress = delegate(int it, double r)
            {
                if (it % 50 == 0 || it < 3)
                    Log(string.Format("  [{0}+{1}] iter={2}, relRes={3:E3}", solverName, precondName, it, r));
                OnIteration?.Invoke(it, r);
            };

            if (SolverKind == KrylovSolverKind.GMRES)
            {
                GMRESSolver solver = new GMRESSolver();
                solver.Tolerance = Tolerance;
                solver.MaxIterations = MaxIterations;
                solver.Restart = 80;
                solver.Preconditioner = M;
                solver.OnProgress += progress;
                conv = solver.Solve(K, rhs, u);
                iters = solver.LastIterations; res = solver.LastResidual;
            }
            else if (SolverKind == KrylovSolverKind.BiCGSTAB)
            {
                BiCGSTABSolver solver = new BiCGSTABSolver();
                solver.Tolerance = Tolerance;
                solver.MaxIterations = MaxIterations;
                solver.Preconditioner = M;
                solver.OnProgress += progress;
                conv = solver.Solve(K, rhs, u);
                iters = solver.LastIterations; res = solver.LastResidual;
            }
            else
            {
                CGSolver solver = new CGSolver();
                solver.Tolerance = Tolerance;
                solver.MaxIterations = MaxIterations;
                solver.Preconditioner = M;
                solver.OnProgress += progress;
                conv = solver.Solve(K, rhs, u);
                iters = solver.LastIterations; res = solver.LastResidual;
            }
            Log(string.Format("[Struct] 求解完成, solver={0}+{1}, iter={2}, res={3:E3}, conv={4}, t={5:F2}s",
                solverName, precondName, iters, res, conv, (DateTime.Now - t0).TotalSeconds));

            // 写回位移
            Array.Clear(mesh.DisplacementU, 0, mesh.NumNodes);
            Array.Clear(mesh.DisplacementV, 0, mesh.NumNodes);
            for (int i = 0; i < mesh.NumNodes; i++)
            {
                mesh.DisplacementU[i] = u.Values[2 * i];
                mesh.DisplacementV[i] = u.Values[2 * i + 1];
            }

            // 单元级应力（单元中心）
            foreach (FiniteElement e in mesh.Elements)
            {
                int nen = e.NodesPerElement;
                double[] ue = new double[2 * nen];
                for (int a = 0; a < nen; a++)
                {
                    ue[2 * a] = u.Values[2 * e.NodeIds[a]];
                    ue[2 * a + 1] = u.Values[2 * e.NodeIds[a] + 1];
                }
                double sxx, syy, sxy, vm;
                double[,] c = e.GetCoordsMatrix();
                e.ComputeStress(c, ue, e.YoungModulus, e.PoissonRatio, out sxx, out syy, out sxy, out vm);
                e.StressXX = sxx; e.StressYY = syy; e.StressXY = sxy; e.VonMises = vm;
                mesh.StressXX[e.Id] = sxx;
                mesh.StressYY[e.Id] = syy;
                mesh.StressXY[e.Id] = sxy;
                mesh.VonMises[e.Id] = vm;
            }
        }

        private void Log(string msg) { OnLog?.Invoke(msg); }
    }
}
