using System;

namespace ModalFEM2D.Core
{

    public delegate void IterationCallback(int iter, double relRes);

    public static class KrylovSolvers
    {
        public static int CG(SparseMatrixCSR A, DenseVector b, DenseVector x,
                             IPreconditioner M, int maxIter, double tol, IterationCallback cb)
        {
            int n = A.NumRows;
            DenseVector r = new DenseVector(n);
            DenseVector z = new DenseVector(n);
            DenseVector p = new DenseVector(n);
            DenseVector Ap = new DenseVector(n);
            A.MatVec(x, r);
            for (int i = 0; i < n; i++) r.Values[i] = b.Values[i] - r.Values[i];
            double bnorm = b.Norm2(); if (bnorm < 1e-30) bnorm = 1.0;
            M.Apply(r, z);
            p.CopyFrom(z);
            double rz = r.Dot(z);
            int it = 0;
            for (; it < maxIter; it++)
            {
                A.MatVec(p, Ap);
                double pAp = p.Dot(Ap);
                if (Math.Abs(pAp) < 1e-30) break;
                double alpha = rz / pAp;
                x.Axpy(alpha, p);
                r.Axpy(-alpha, Ap);
                double rn = r.Norm2();
                if (cb != null) cb(it + 1, rn / bnorm);
                if (rn / bnorm < tol) return it + 1;
                M.Apply(r, z);
                double rzNew = r.Dot(z);
                double beta = rzNew / rz;
                rz = rzNew;
                for (int i = 0; i < n; i++) p.Values[i] = z.Values[i] + beta * p.Values[i];
            }
            return it;
        }

        public static int GMRES(SparseMatrixCSR A, DenseVector b, DenseVector x,
                                IPreconditioner M, int m, int maxOuter, double tol, IterationCallback cb)
        {
            int n = A.NumRows;
            double bnrm = b.Norm2(); if (bnrm < 1e-30) bnrm = 1.0;
            int totalIt = 0;
            DenseVector r = new DenseVector(n);
            DenseVector w = new DenseVector(n);
            DenseVector[] V = new DenseVector[m + 1];
            for (int i = 0; i <= m; i++) V[i] = new DenseVector(n);
            double[,] H = new double[m + 1, m];
            double[] cs = new double[m];
            double[] sn = new double[m];
            double[] g = new double[m + 1];

            for (int outer = 0; outer < maxOuter; outer++)
            {
                A.MatVec(x, r);
                for (int i = 0; i < n; i++) r.Values[i] = b.Values[i] - r.Values[i];
                M.Apply(r, V[0]);
                double beta = V[0].Norm2();
                if (beta < 1e-30) return totalIt;
                V[0].Scale(1.0 / beta);
                for (int i = 0; i <= m; i++) g[i] = 0.0;
                g[0] = beta;
                int j;
                for (j = 0; j < m; j++)
                {
                    A.MatVec(V[j], w);
                    DenseVector ww = new DenseVector(n);
                    M.Apply(w, ww);
                    w = ww;
                    for (int i = 0; i <= j; i++)
                    {
                        H[i, j] = V[i].Dot(w);
                        w.Axpy(-H[i, j], V[i]);
                    }
                    H[j + 1, j] = w.Norm2();
                    if (H[j + 1, j] < 1e-30) { j++; break; }
                    V[j + 1].CopyFrom(w); V[j + 1].Scale(1.0 / H[j + 1, j]);
                    // Apply previous Givens
                    for (int i = 0; i < j; i++)
                    {
                        double t = cs[i] * H[i, j] + sn[i] * H[i + 1, j];
                        H[i + 1, j] = -sn[i] * H[i, j] + cs[i] * H[i + 1, j];
                        H[i, j] = t;
                    }
                    double d = Math.Sqrt(H[j, j] * H[j, j] + H[j + 1, j] * H[j + 1, j]);
                    cs[j] = H[j, j] / d; sn[j] = H[j + 1, j] / d;
                    H[j, j] = cs[j] * H[j, j] + sn[j] * H[j + 1, j];
                    H[j + 1, j] = 0.0;
                    double gtmp = cs[j] * g[j];
                    g[j + 1] = -sn[j] * g[j];
                    g[j] = gtmp;
                    double res = Math.Abs(g[j + 1]) / bnrm;
                    totalIt++;
                    if (cb != null) cb(totalIt, res);
                    if (res < tol) { j++; break; }
                }
                // Backsolve H y = g (top j x j)
                int kk = j - 1;
                double[] y = new double[j];
                for (int i = kk; i >= 0; i--)
                {
                    double s = g[i];
                    for (int l = i + 1; l < j; l++) s -= H[i, l] * y[l];
                    y[i] = s / H[i, i];
                }
                DenseVector dx = new DenseVector(n);
                for (int i = 0; i < j; i++) dx.Axpy(y[i], V[i]);
                // Apply M^-1 on the Krylov combination? We already preconditioned the search directions;
                // V columns are already in preconditioned space. Add directly.
                x.Axpy(1.0, dx);
                A.MatVec(x, r);
                for (int i = 0; i < n; i++) r.Values[i] = b.Values[i] - r.Values[i];
                double rn = r.Norm2();
                if (rn / bnrm < tol) return totalIt;
                if (j < m) break;
            }
            return totalIt;
        }

        public static int BiCGSTAB(SparseMatrixCSR A, DenseVector b, DenseVector x,
                                   IPreconditioner M, int maxIter, double tol, IterationCallback cb)
        {
            int n = A.NumRows;
            DenseVector r = new DenseVector(n);
            DenseVector r0 = new DenseVector(n);
            DenseVector p = new DenseVector(n);
            DenseVector v = new DenseVector(n);
            DenseVector s = new DenseVector(n);
            DenseVector t = new DenseVector(n);
            DenseVector phat = new DenseVector(n);
            DenseVector shat = new DenseVector(n);
            A.MatVec(x, r);
            for (int i = 0; i < n; i++) r.Values[i] = b.Values[i] - r.Values[i];
            r0.CopyFrom(r);
            p.CopyFrom(r);
            double rho = 1.0, alpha = 1.0, omega = 1.0;
            double bn = b.Norm2(); if (bn < 1e-30) bn = 1.0;
            int it = 0;
            for (; it < maxIter; it++)
            {
                double rhoNew = r0.Dot(r);
                if (Math.Abs(rhoNew) < 1e-30) break;
                double beta = (rhoNew / rho) * (alpha / omega);
                rho = rhoNew;
                for (int i = 0; i < n; i++) p.Values[i] = r.Values[i] + beta * (p.Values[i] - omega * v.Values[i]);
                M.Apply(p, phat);
                A.MatVec(phat, v);
                double tv = r0.Dot(v);
                if (Math.Abs(tv) < 1e-30) break;
                alpha = rho / tv;
                s.CopyFrom(r); s.Axpy(-alpha, v);
                double sn = s.Norm2();
                if (sn / bn < tol) { x.Axpy(alpha, phat); if (cb != null) cb(it + 1, sn / bn); return it + 1; }
                M.Apply(s, shat);
                A.MatVec(shat, t);
                double tt = t.Dot(t);
                omega = tt > 1e-30 ? t.Dot(s) / tt : 0.0;
                x.Axpy(alpha, phat);
                x.Axpy(omega, shat);
                r.CopyFrom(s); r.Axpy(-omega, t);
                double rn = r.Norm2();
                if (cb != null) cb(it + 1, rn / bn);
                if (rn / bn < tol) return it + 1;
                if (Math.Abs(omega) < 1e-30) break;
            }
            return it;
        }
    }
}
