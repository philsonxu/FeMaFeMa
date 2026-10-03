using System;
using System.Threading;
using NonlinearFEM2D.Core;
using NonlinearFEM2D.Elements;
using NonlinearFEM2D.Mesh;

namespace NonlinearFEM2D.Solvers
{
    public enum KrylovKind { CG, GMRES, BiCGSTAB }
    public enum PrecondKind { None, Jacobi, AMG }

    public sealed class StructuralSolver
    {
        public KrylovKind Krylov = KrylovKind.CG;
        public PrecondKind Precond = PrecondKind.Jacobi;
        public double Tolerance = 1.0e-8;
        public int MaxIterations = 5000;
        public Action<string> Log;
        public Action<int, double> OnIter;

        public void Solve(FEMesh mesh)
        {
            int ndof = 2 * mesh.NumNodes;
            CooBuilder Kcoo = new CooBuilder(ndof, ndof);
            DenseVector f = new DenseVector(ndof);
            int nElem = mesh.Elements.Count;
            int cores = Math.Max(1, Environment.ProcessorCount);
            int chunk = (nElem + cores - 1) / cores;
            ManualResetEvent[] evts = new ManualResetEvent[cores];
            for (int t = 0; t < cores; t++) evts[t] = new ManualResetEvent(false);
            for (int tid = 0; tid < cores; tid++)
            {
                int begin = tid * chunk;
                int end = Math.Min(begin + chunk, nElem);
                ThreadPool.QueueUserWorkItem(delegate (object state)
                {
                    CooBuilder local = new CooBuilder(ndof, ndof);
                    for (int ie = begin; ie < end; ie++)
                    {
                        FiniteElement e = mesh.Elements[ie];
                        int npe = e.NodeIds.Length;
                        double[,] Ke = new double[npe * 2, npe * 2];
                        e.BuildLinearStiffness(Ke);
                        int[] dofs = new int[npe * 2];
                        for (int k = 0; k < npe; k++) { dofs[2 * k] = 2 * e.NodeIds[k]; dofs[2 * k + 1] = 2 * e.NodeIds[k] + 1; }
                        local.AddLocal(dofs, Ke);
                    }
                    lock (Kcoo)
                    {
                        // merge local -> Kcoo
                        SparseMatrixCSR tmp = local.Build();
                        for (int i = 0; i < tmp.Rows; i++)
                            for (int kk = tmp.RowPtr[i]; kk < tmp.RowPtr[i + 1]; kk++)
                                Kcoo.Add(i, tmp.ColIdx[kk], tmp.Values[kk]);
                    }
                    evts[(int)state].Set();
                }, tid);
            }
            WaitHandle.WaitAll(evts);
            for (int i = 0; i < mesh.NumNodes; i++)
            {
                f[2 * i] = mesh.NodeForceX[i];
                f[2 * i + 1] = mesh.NodeForceY[i];
            }
            double penalty = 1.0e30;
            foreach (BoundaryCondition bc in mesh.Dirichlet)
            {
                int row = 2 * bc.NodeId + bc.Dof;
                Kcoo.Add(row, row, penalty);
                f[row] = bc.Value * penalty;
            }
            SparseMatrixCSR K = Kcoo.Build();
            IPreconditioner M = BuildPrecond(this.Precond, K);
            ILinearSolver solver = BuildSolver(this.Krylov);
            solver.Tolerance = this.Tolerance; solver.MaxIterations = this.MaxIterations;
            if (this.OnIter != null) solver.OnIteration = this.OnIter;
            DenseVector u = new DenseVector(ndof);
            int iters = solver.Solve(K, f, u, M);
            if (this.Log != null) this.Log(string.Format("[线性静力] {0}+{1} iters={2} relRes={3:E4}",
                this.Krylov, this.Precond, iters, solver.FinalResidual));
            for (int i = 0; i < mesh.NumNodes; i++)
            {
                mesh.DisplacementU[i] = u[2 * i];
                mesh.DisplacementV[i] = u[2 * i + 1];
            }
            // 应力恢复（按单元平均→节点）
            int[] cnt = new int[mesh.NumNodes];
            Array.Clear(mesh.StressXX, 0, mesh.NumNodes);
            Array.Clear(mesh.StressYY, 0, mesh.NumNodes);
            Array.Clear(mesh.StressXY, 0, mesh.NumNodes);
            Array.Clear(mesh.VonMises, 0, mesh.NumNodes);
            foreach (FiniteElement e in mesh.Elements)
            {
                int npe = e.NodeIds.Length;
                double[] eu = new double[npe], ev = new double[npe];
                for (int k = 0; k < npe; k++) { eu[k] = mesh.DisplacementU[e.NodeIds[k]]; ev[k] = mesh.DisplacementV[e.NodeIds[k]]; }
                double sx, sy, sxy;
                e.ComputeStress(eu, ev, out sx, out sy, out sxy);
                double vm = Math.Sqrt(Math.Max(0.0, sx * sx - sx * sy + sy * sy + 3.0 * sxy * sxy));
                for (int k = 0; k < npe; k++)
                {
                    int nid = e.NodeIds[k];
                    mesh.StressXX[nid] += sx; mesh.StressYY[nid] += sy; mesh.StressXY[nid] += sxy; mesh.VonMises[nid] += vm;
                    cnt[nid]++;
                }
            }
            for (int i = 0; i < mesh.NumNodes; i++) if (cnt[i] > 0)
                {
                    mesh.StressXX[i] /= cnt[i]; mesh.StressYY[i] /= cnt[i]; mesh.StressXY[i] /= cnt[i]; mesh.VonMises[i] /= cnt[i];
                }
        }

        internal static IPreconditioner BuildPrecond(PrecondKind pk, SparseMatrixCSR K)
        {
            IPreconditioner M;
            if (pk == PrecondKind.Jacobi) M = new JacobiPreconditioner();
            else if (pk == PrecondKind.AMG) M = new AMGPreconditioner();
            else M = new IdentityPreconditioner();
            M.Build(K);
            return M;
        }
        internal static ILinearSolver BuildSolver(KrylovKind ks)
        {
            if (ks == KrylovKind.GMRES) return new GMRESSolver();
            if (ks == KrylovKind.BiCGSTAB) return new BiCGSTABSolver();
            return new CGSolver();
        }
    }
}
