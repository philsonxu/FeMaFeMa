using System;
using System.Collections.Generic;
using System.Threading;
using MultiPhysicsFEM2D.Core;
using MultiPhysicsFEM2D.Elements;
using MultiPhysicsFEM2D.Mesh;

namespace MultiPhysicsFEM2D.Solvers
{
    public enum KrylovKind { CG, GMRES, BiCGSTAB }
    public enum PrecondKind { None, Jacobi, AMG }

    public class StructuralResult
    {
        public int Iterations;
        public double Residual;
        public string SolverUsed;
    }

    public class StructuralSolver
    {
        public KrylovKind Krylov = KrylovKind.CG;
        public PrecondKind Precond = PrecondKind.Jacobi;
        public double Penalty = 1e20;
        public double Tol = 1e-8;
        public int MaxIt = 5000;
        public Action<string> Log;

        public StructuralResult Solve(FEMesh mesh)
        {
            int nNodes = mesh.NumNodes;
            int ndof = 2 * nNodes;
            CooAssembler coo = new CooAssembler(ndof, ndof);
            DenseVector f = new DenseVector(ndof);
            int threads = Environment.ProcessorCount;
            int block = (mesh.Elements.Count + threads - 1) / threads;
            WaitHandle[] whs = new WaitHandle[threads];
            List<Tuple<int[], double[,]>>[] perBatch = new List<Tuple<int[], double[,]>>[threads];
            double[][] perF = new double[threads][];
            for (int t = 0; t < threads; t++) { perBatch[t] = new List<Tuple<int[], double[,]>>(); perF[t] = new double[ndof]; }
            for (int t = 0; t < threads; t++)
            {
                int start = t * block;
                int end = Math.Min(start + block, mesh.Elements.Count);
                int tt = t;
                ManualResetEvent mre = new ManualResetEvent(false);
                whs[t] = mre;
                ThreadPool.QueueUserWorkItem(delegate
                {
                    double[] fLocal = perF[tt];
                    double[,] Ke = new double[6, 6]; // 最大是 LT6 12
                    for (int ei = start; ei < end; ei++)
                    {
                        FiniteElement e = mesh.Elements[ei];
                        int nn = e.NodesPerElement;
                        int ldof = 2 * nn;
                        double[,] Klocal;
                        if (nn == 6) Klocal = new double[12, 12];
                        else if (nn == 4) Klocal = new double[8, 8];
                        else Klocal = new double[6, 6];
                        for (int i = 0; i < ldof; i++) for (int j = 0; j < ldof; j++) Klocal[i, j] = 0;
                        e.BuildStiffnessMatrix(Klocal);
                        int[] dofMap = new int[ldof];
                        for (int k = 0; k < nn; k++) { dofMap[2 * k] = 2 * e.NodeIds[k]; dofMap[2 * k + 1] = 2 * e.NodeIds[k] + 1; }
                        perBatch[tt].Add(Tuple.Create(dofMap, Klocal));
                    }
                    mre.Set();
                });
            }
            WaitHandle.WaitAll(whs);
            for (int t = 0; t < threads; t++)
            {
                coo.AddLocalBulk(perBatch[t]);
                for (int i = 0; i < ndof; i++) f.Values[i] += perF[t][i];
            }
            // 施加集中力
            foreach (BoundaryCondition bc in mesh.BCs)
            {
                if (bc.Type == BCType.TractionX) f.Values[2 * bc.NodeId] += bc.Value;
                if (bc.Type == BCType.TractionY) f.Values[2 * bc.NodeId + 1] += bc.Value;
            }
            SparseMatrixCSR K = coo.BuildCSR();
            DenseVector u = new DenseVector(ndof);
            // Dirichlet
            foreach (BoundaryCondition bc in mesh.BCs)
            {
                if (bc.Type == BCType.FixX) K.ApplyDirichlet(2 * bc.NodeId, Penalty, bc.Value, f);
                else if (bc.Type == BCType.FixY) K.ApplyDirichlet(2 * bc.NodeId + 1, Penalty, bc.Value, f);
                else if (bc.Type == BCType.FixBoth)
                {
                    K.ApplyDirichlet(2 * bc.NodeId, Penalty, bc.Value, f);
                    K.ApplyDirichlet(2 * bc.NodeId + 1, Penalty, bc.Value, f);
                }
            }
            IPreconditioner pc;
            if (Precond == PrecondKind.None) pc = new IdentityPreconditioner();
            else if (Precond == PrecondKind.AMG) { AMGPreconditioner ap = new AMGPreconditioner(); ap.Build(K); pc = ap; }
            else pc = new JacobiPreconditioner(); pc.Build(K);
            int iters = -1; string solverName = Krylov.ToString();
            IterationCallback cb = null;
            if (Log != null) cb = delegate (int it, double r)
            {
                if (it % 20 == 0 || it < 5 || r < Tol)
                    Log(string.Format("  {0} iter {1}: rel.res = {2:E4}", solverName, it, r));
            };
            if (Krylov == KrylovKind.CG) iters = KrylovSolvers.CG(K, f, u, pc, Tol, MaxIt, cb);
            else if (Krylov == KrylovKind.GMRES) iters = KrylovSolvers.GMRES(K, f, u, pc, Tol, KrylovSolvers.GMRES_MAXIT, KrylovSolvers.GMRES_M, cb);
            else iters = KrylovSolvers.BiCGSTAB(K, f, u, pc, Tol, KrylovSolvers.BICG_MAXIT, cb);
            // 写回位移
            for (int i = 0; i < nNodes; i++)
            {
                mesh.DisplacementU[i] = u.Values[2 * i];
                mesh.DisplacementV[i] = u.Values[2 * i + 1];
            }
            // 计算单元应力
            foreach (FiniteElement e in mesh.Elements)
            {
                int nn = e.NodesPerElement;
                double[] eu = new double[nn]; double[] ev = new double[nn];
                for (int k = 0; k < nn; k++) { eu[k] = mesh.DisplacementU[e.NodeIds[k]]; ev[k] = mesh.DisplacementV[e.NodeIds[k]]; }
                double sxx, syy, sxy, vm;
                e.ComputeStress(eu, ev, out sxx, out syy, out sxy, out vm);
                e.StressXX = sxx; e.StressYY = syy; e.StressXY = sxy; e.VonMises = vm;
            }
            mesh.ProjectStressFromElements();
            // 残差
            DenseVector rr = new DenseVector(ndof);
            K.Multiply(u, rr);
            double residual = 0;
            for (int i = 0; i < ndof; i++) residual += (f.Values[i] - rr.Values[i]) * (f.Values[i] - rr.Values[i]);
            residual = Math.Sqrt(residual);
            return new StructuralResult { Iterations = iters, Residual = residual, SolverUsed = solverName + "+" + pc.GetType().Name };
        }
    }
}
