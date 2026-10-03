using System;
using System.Threading.Tasks;
using FEM2D.Core;
using FEM2D.Elements;
using FEM2D.Mesh;

namespace FEM2D.Solvers
{
    /// <summary>
    /// 结构线弹性静力分析结果
    /// </summary>
    public sealed class StructuralResult
    {
        public DenseVector Displacement;  // 2*N
        public DenseVector Reaction;      // 2*N
        public double[] StressXX;
        public double[] StressYY;
        public double[] StressXY;
        public double[] VonMises;
        public double MaxVM;
        public double MaxDisp;
    }

    /// <summary>
    /// 线弹性平面应变结构求解器。使用多线程单元组装，CG 或 GMRES 求解。
    /// </summary>
    public sealed class StructuralSolver
    {
        public double YoungModulus { get; set; } = 210e9;
        public double PoissonRatio { get; set; } = 0.3;
        public double Thickness { get; set; } = 1.0;
        public bool UseGMRES { get; set; } = false;
        public int Iterations { get; private set; }

        public StructuralResult Solve(FEMesh mesh, Action<string> log = null)
        {
            int ndof = 2 * mesh.NumNodes;
            int nElem = mesh.NumElements;
            FiniteElement[] elems = new FiniteElement[nElem];
            for (int e = 0; e < nElem; e++) elems[e] = ElementFactory.Create(mesh.ElementTypes[e]);

            SparseMatrixCSR.CooBuilder builder = new SparseMatrixCSR.CooBuilder(ndof);
            DenseVector rhs = new DenseVector(ndof);

            double[,] ke = null;
            int maxNpe = 6;
            // 单元刚度矩阵临时数组（按该单元 NPE）
            object lockKe = new object();

            if (log != null) log($"开始组装刚度矩阵，共 {nElem} 个单元，{ndof} 个自由度…");

            Parallel.For(0, nElem, e =>
            {
                FiniteElement fe = elems[e];
                int npe = mesh.Elements[e].Length;
                int edof = 2 * npe;
                double[,] kel = new double[edof, edof];
                fe.ComputeStiffness(mesh, e, YoungModulus, PoissonRatio, Thickness, kel);

                int[] conn = mesh.Elements[e];
                int[] dofMap = new int[edof];
                for (int k = 0; k < npe; k++)
                {
                    dofMap[2 * k] = 2 * conn[k];
                    dofMap[2 * k + 1] = 2 * conn[k] + 1;
                }
                lock (lockKe) { builder.AddLocal(dofMap, kel); }
            });

            // 集中力
            foreach ((int node, int dof, double val) in mesh.ConcentratedForces)
            {
                rhs[2 * node + dof] += val;
            }

            SparseMatrixCSR K = builder.Build(true);
            if (log != null) log("刚度矩阵装配完成，非零元：" + K.Values.Length);

            // 罚函数施加 Dirichlet 边界条件
            double penalty = 1e20;
            double[] diagInv = K.DiagonalInverse();
            DenseVector u = new DenseVector(ndof);
            foreach ((int node, int dof, double val) in mesh.DirichletBCs)
            {
                int row = 2 * node + dof;
                double kv = 1.0 / diagInv[row] * penalty;
                // 重建 COO 增加对角项
                // 这里直接通过 COO builder 重建一个新矩阵较繁；改为右端项修改+简单方式：
                rhs[row] = kv * val;
                // 向 builder 追加
                builder.Add(row, row, kv - (1.0 / diagInv[row]));
                u[row] = val;
            }
            K = builder.Build(true);
            diagInv = K.DiagonalInverse();

            if (log != null) log("边界条件施加完成，开始求解线性系统…");

            bool ok;
            if (UseGMRES)
            {
                GMRESSolver gmres = new GMRESSolver { Tolerance = 1e-8, KrylovDim = 80, MaxRestart = 200, PrintLevel = 0 };
                ok = gmres.Solve(K, rhs, u, diagInv);
                Iterations = gmres.FinalIterations;
                if (log != null) log($"GMRES 收敛，迭代 = {Iterations}, 残差 = {gmres.FinalResidual:E4}");
            }
            else
            {
                CGSolver cg = new CGSolver { Tolerance = 1e-8, MaxIterations = 10000 };
                ok = cg.Solve(K, rhs, u, diagInv);
                Iterations = cg.FinalIterations;
                if (log != null) log($"CG 收敛，迭代 = {Iterations}, 残差 = {cg.FinalResidual:E4}");
            }

            // 反力 R = K u - f
            DenseVector R = new DenseVector(ndof);
            K.MatVec(u, R);
            for (int i = 0; i < ndof; i++) R[i] -= rhs[i];

            // 单元应力恢复
            double[] sx = new double[nElem];
            double[] sy = new double[nElem];
            double[] sxy = new double[nElem];
            double[] vm = new double[nElem];
            double maxVM = 0;
            double maxDisp = 0;

            Parallel.For(0, nElem, e =>
            {
                int[] conn = mesh.Elements[e];
                int npe = conn.Length;
                double[] ue = new double[2 * npe];
                for (int k = 0; k < npe; k++)
                {
                    ue[2 * k] = u[2 * conn[k]];
                    ue[2 * k + 1] = u[2 * conn[k] + 1];
                }
                double[] s = elems[e].ComputeStress(mesh, e, YoungModulus, PoissonRatio, ue);
                sx[e] = s[0]; sy[e] = s[1]; sxy[e] = s[2]; vm[e] = s[3];
            });
            for (int e = 0; e < nElem; e++) if (Math.Abs(vm[e]) > maxVM) maxVM = Math.Abs(vm[e]);
            for (int i = 0; i < ndof; i++) { double d = Math.Abs(u[i]); if (d > maxDisp) maxDisp = d; }

            return new StructuralResult
            {
                Displacement = u,
                Reaction = R,
                StressXX = sx,
                StressYY = sy,
                StressXY = sxy,
                VonMises = vm,
                MaxVM = maxVM,
                MaxDisp = maxDisp
            };
        }
    }
}
