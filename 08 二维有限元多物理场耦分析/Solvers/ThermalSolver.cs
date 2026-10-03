using System;
using System.Collections;
using System.Collections.Generic;
using System.Threading;
using MultiPhysicsFEM2D.Core;
using MultiPhysicsFEM2D.Elements;
using MultiPhysicsFEM2D.Mesh;

namespace MultiPhysicsFEM2D.Solvers
{
    public class ThermalResult
    {
        public int Iterations;
        public double Residual;
        public bool Steady;
        public int TimeSteps;
    }

    public class ThermalSolver
    {
        public bool Steady = true;
        public double TimeStep = 0.01;
        public int TimeSteps = 10;
        public double Penalty = 1e20;
        public double Tol = 1e-8;
        public int MaxIt = 3000;
        public KrylovKind Krylov = KrylovKind.CG;
        public PrecondKind Precond = PrecondKind.Jacobi;
        public Action<string> Log;
        public double BodySource = 0.0; // Q [W/m^3] 体热源

        public ThermalResult Solve(FEMesh mesh)
        {
            int nNodes = mesh.NumNodes;
            int ndof = nNodes;
            CooAssembler coo = new CooAssembler(ndof, ndof);
            DenseVector rhs = new DenseVector(ndof);
            int nElem = mesh.Elements.Count;
            int threads = Environment.ProcessorCount;
            int block = (nElem + threads - 1) / threads;
            WaitHandle[] whs = new WaitHandle[threads];
            List<Tuple<int[], double[,]>>[] perB = new List<Tuple<int[], double[,]>>[threads];
            double[][] perF = new double[threads][];
            for (int t = 0; t < threads; t++) { perB[t] = new List<Tuple<int[], double[,]>>(); perF[t] = new double[ndof]; }
            for (int t = 0; t < threads; t++)
            {
                int start = t * block, end = Math.Min(start + block, nElem);
                int tt = t;
                ManualResetEvent mre = new ManualResetEvent(false); whs[t] = mre;
                ThreadPool.QueueUserWorkItem(delegate
                {
                    for (int ei = start; ei < end; ei++)
                    {
                        FiniteElement e = mesh.Elements[ei];
                        int nn = e.NodesPerElement;
                        double[,] Kt = new double[nn, nn];
                        double[,] Ct = new double[nn, nn];
                        e.BuildConductivityMatrix(Kt);
                        e.BuildCapacityMatrix(Ct, true);
                        double[,] Mat;
                        if (Steady) { Mat = Kt; }
                        else
                        {
                            Mat = new double[nn, nn];
                            for (int i = 0; i < nn; i++) for (int j = 0; j < nn; j++) Mat[i, j] = Ct[i, j] / TimeStep + Kt[i, j];
                        }
                        int[] dofMap = new int[nn];
                        for (int k = 0; k < nn; k++) dofMap[k] = e.NodeIds[k];
                        perB[tt].Add(Tuple.Create(dofMap, Mat));
                        if (!Steady)
                        {
                            for (int i = 0; i < nn; i++)
                            {
                                double s = 0;
                                for (int j = 0; j < nn; j++) s += Ct[i, j] / TimeStep * mesh.Temperature[e.NodeIds[j]];
                                perF[tt][e.NodeIds[i]] += s;
                            }
                        }
                        // 体热源 Q：f_i += ∫ N_i Q = Q * area * thickness / nn
                        if (Math.Abs(BodySource) > 1e-12)
                        {
                            double a = e.ComputeArea() * e.Thickness / nn;
                            for (int i = 0; i < nn; i++) perF[tt][e.NodeIds[i]] += BodySource * a;
                        }
                    }
                    mre.Set();
                });
            }
            WaitHandle.WaitAll(whs);
            for (int t = 0; t < threads; t++)
            {
                coo.AddLocalBulk(perB[t]);
                for (int i = 0; i < ndof; i++) rhs.Values[i] += perF[t][i];
            }
            SparseMatrixCSR A = coo.BuildCSR();
            DenseVector T = new DenseVector(ndof);
            // Dirichlet
            foreach (BoundaryCondition bc in mesh.BCs)
            {
                if (bc.Type == BCType.Temperature) { A.ApplyDirichlet(bc.NodeId, Penalty, bc.Value, rhs); T.Values[bc.NodeId] = bc.Value; }
                else if (bc.Type == BCType.HeatFlux) { rhs.Values[bc.NodeId] += bc.Value; }
            }
            IPreconditioner pc = new JacobiPreconditioner();
            if (Precond == PrecondKind.AMG) { AMGPreconditioner ap = new AMGPreconditioner(); ap.Build(A); pc = ap; }
            else pc.Build(A);
            IterationCallback cb = null;
            if (Log != null) cb = delegate (int it, double r) { if (it % 10 == 0 || r < Tol) Log(string.Format("  Thermal {0} iter {1}: res={2:E4}", Krylov, it, r)); };
            int iters = -1;
            if (Krylov == KrylovKind.CG) iters = KrylovSolvers.CG(A, rhs, T, pc, Tol, MaxIt, cb);
            else if (Krylov == KrylovKind.BiCGSTAB) iters = KrylovSolvers.BiCGSTAB(A, rhs, T, pc, Tol, MaxIt, cb);
            else iters = KrylovSolvers.GMRES(A, rhs, T, pc, Tol, KrylovSolvers.GMRES_MAXIT, KrylovSolvers.GMRES_M, cb);
            for (int i = 0; i < nNodes; i++) mesh.Temperature[i] = T.Values[i];
            DenseVector rr = new DenseVector(ndof); A.Multiply(T, rr);
            double res = 0;
            for (int i = 0; i < ndof; i++) res += (rhs.Values[i] - rr.Values[i]) * (rhs.Values[i] - rr.Values[i]);
            return new ThermalResult { Iterations = iters, Residual = Math.Sqrt(res), Steady = Steady, TimeSteps = Steady ? 1 : TimeSteps };
        }

        public ThermalResult SolveTransient(FEMesh mesh)
        {
            // 初始温度
            for (int i = 0; i < mesh.NumNodes; i++) mesh.Temperature[i] = mesh.ReferenceTemperature;
            ThermalResult last = null;
            bool steady0 = Steady;
            Steady = false;
            for (int s = 0; s < TimeSteps; s++)
            {
                if (Log != null) Log(string.Format("Thermal step {0}/{1}  dt={2}", s + 1, TimeSteps, TimeStep));
                last = Solve(mesh);
            }
            Steady = steady0;
            return last;
        }
    }
}
