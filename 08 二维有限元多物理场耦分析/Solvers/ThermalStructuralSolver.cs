using System;
using System.Collections;
using System.Collections.Generic;
using System.Threading;
using MultiPhysicsFEM2D.Core;
using MultiPhysicsFEM2D.Elements;
using MultiPhysicsFEM2D.Mesh;

namespace MultiPhysicsFEM2D.Solvers
{
    public class ThermalStructuralResult
    {
        public ThermalResult Thermal;
        public StructuralResult Structural;
    }

    public class ThermalStructuralSolver
    {
        public bool SteadyThermal = true;
        public double TimeStep = 0.01;
        public int TimeSteps = 10;
        public KrylovKind Krylov = KrylovKind.CG;
        public PrecondKind Precond = PrecondKind.Jacobi;
        public double BodySource = 0.0;
        public Action<string> Log;

        public ThermalStructuralResult Solve(FEMesh mesh)
        {
            ThermalStructuralResult res = new ThermalStructuralResult();
            // 先解热传导
            ThermalSolver ts = new ThermalSolver();
            ts.Steady = SteadyThermal;
            ts.TimeStep = TimeStep; ts.TimeSteps = TimeSteps;
            ts.Krylov = Krylov; ts.Precond = Precond;
            ts.BodySource = BodySource;
            ts.Log = Log;
            if (Log != null) Log("=== Phase 1: Thermal ===");
            if (SteadyThermal) res.Thermal = ts.Solve(mesh);
            else res.Thermal = ts.SolveTransient(mesh);

            // 装配结构矩阵，把热应变载荷叠加到 f
            if (Log != null) Log("=== Phase 2: Structural with thermal load ===");
            int nNodes = mesh.NumNodes;
            int ndof = 2 * nNodes;
            CooAssembler coo = new CooAssembler(ndof, ndof);
            DenseVector f = new DenseVector(ndof);
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
                        int ldof = 2 * nn;
                        double[,] K = new double[ldof, ldof];
                        el.BuildStiffnessMatrix(K);
                        double[] floc = new double[ldof];
                        double[] TN = new double[nn];
                        for (int k = 0; k < nn; k++) TN[k] = mesh.Temperature[el.NodeIds[k]];
                        el.BuildThermalLoadVector(TN, el.ReferenceTemperature, floc);
                        int[] dofMap = new int[ldof];
                        for (int k = 0; k < nn; k++) { dofMap[2 * k] = 2 * el.NodeIds[k]; dofMap[2 * k + 1] = 2 * el.NodeIds[k] + 1; }
                        perB[tt].Add(Tuple.Create(dofMap, K));
                        for (int i = 0; i < ldof; i++) perF[tt][dofMap[i]] += floc[i];
                    }
                    mre.Set();
                });
            }
            WaitHandle.WaitAll(whs);
            for (int t = 0; t < threads; t++)
            {
                coo.AddLocalBulk(perB[t]);
                for (int i = 0; i < ndof; i++) f.Values[i] += perF[t][i];
            }
            // 集中力
            foreach (BoundaryCondition bc in mesh.BCs)
            {
                if (bc.Type == BCType.TractionX) f.Values[2 * bc.NodeId] += bc.Value;
                if (bc.Type == BCType.TractionY) f.Values[2 * bc.NodeId + 1] += bc.Value;
            }
            SparseMatrixCSR K = coo.BuildCSR();
            DenseVector u = new DenseVector(ndof);
            foreach (BoundaryCondition bc in mesh.BCs)
            {
                if (bc.Type == BCType.FixX) K.ApplyDirichlet(2 * bc.NodeId, 1e20, bc.Value, f);
                else if (bc.Type == BCType.FixY) K.ApplyDirichlet(2 * bc.NodeId + 1, 1e20, bc.Value, f);
                else if (bc.Type == BCType.FixBoth)
                {
                    K.ApplyDirichlet(2 * bc.NodeId, 1e20, bc.Value, f);
                    K.ApplyDirichlet(2 * bc.NodeId + 1, 1e20, bc.Value, f);
                }
            }
            IPreconditioner pc = new JacobiPreconditioner();
            if (Precond == PrecondKind.AMG) { AMGPreconditioner ap = new AMGPreconditioner(); ap.Build(K); pc = ap; }
            else pc.Build(K);
            IterationCallback cb = null;
            if (Log != null) cb = delegate (int it, double r) { if (it % 20 == 0 || r < 1e-8) Log(string.Format("  Structure {0} iter {1}: res={2:E4}", Krylov, it, r)); };
            int iters = -1;
            if (Krylov == KrylovKind.CG) iters = KrylovSolvers.CG(K, f, u, pc, 1e-8, 5000, cb);
            else if (Krylov == KrylovKind.GMRES) iters = KrylovSolvers.GMRES(K, f, u, pc, 1e-8, KrylovSolvers.GMRES_MAXIT, KrylovSolvers.GMRES_M, cb);
            else iters = KrylovSolvers.BiCGSTAB(K, f, u, pc, 1e-8, KrylovSolvers.BICG_MAXIT, cb);
            for (int i = 0; i < nNodes; i++)
            {
                mesh.DisplacementU[i] = u.Values[2 * i];
                mesh.DisplacementV[i] = u.Values[2 * i + 1];
            }
            foreach (FiniteElement el in mesh.Elements)
            {
                int nn = el.NodesPerElement;
                double[] eu = new double[nn], ev = new double[nn];
                for (int k = 0; k < nn; k++) { eu[k] = mesh.DisplacementU[el.NodeIds[k]]; ev[k] = mesh.DisplacementV[el.NodeIds[k]]; }
                double sxx, syy, sxy, vm;
                el.ComputeStress(eu, ev, out sxx, out syy, out sxy, out vm);
                el.StressXX = sxx; el.StressYY = syy; el.StressXY = sxy; el.VonMises = vm;
            }
            mesh.ProjectStressFromElements();
            DenseVector rr = new DenseVector(ndof); K.Multiply(u, rr);
            double resid = 0;
            for (int i = 0; i < ndof; i++) resid += (f.Values[i] - rr.Values[i]) * (f.Values[i] - rr.Values[i]);
            res.Structural = new StructuralResult { Iterations = iters, Residual = Math.Sqrt(resid), SolverUsed = Krylov.ToString() };
            return res;
        }
    }
}
