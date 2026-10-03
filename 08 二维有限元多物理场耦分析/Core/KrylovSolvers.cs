using System;

namespace MultiPhysicsFEM2D.Core
{
    public delegate void IterationCallback(int iter, double res);

    public static class KrylovSolvers
    {
        public const int CG_MAXIT = 5000;
        public const int GMRES_M = 60;
        public const int GMRES_MAXIT = 200;
        public const int BICG_MAXIT = 5000;

        public static int CG(SparseMatrixCSR A, DenseVector b, DenseVector x,
            IPreconditioner M, double tol, int maxit, IterationCallback cb)
        {
            int n = A.Nrows;
            DenseVector r = new DenseVector(n);
            DenseVector z = new DenseVector(n);
            DenseVector p = new DenseVector(n);
            DenseVector Ap = new DenseVector(n);
            A.Multiply(x, r);
            for (int i = 0; i < n; i++) r.Values[i] = b.Values[i] - r.Values[i];
            double bnorm = b.Norm2(); if (bnorm < 1e-30) bnorm = 1.0;
            M.Apply(r, z);
            p.CopyFrom(z);
            double rz = r.Dot(z);
            for (int it = 1; it <= maxit; it++)
            {
                A.Multiply(p, Ap);
                double pAp = p.Dot(Ap);
                if (Math.Abs(pAp) < 1e-30) return -1;
                double alpha = rz / pAp;
                x.Axpy(alpha, p);
                r.Axpy(-alpha, Ap);
                double res = r.Norm2() / bnorm;
                if (cb != null) cb(it, res);
                if (res < tol) return it;
                M.Apply(r, z);
                double rzNew = r.Dot(z);
                double beta = rzNew / rz;
                rz = rzNew;
                for (int i = 0; i < n; i++) p.Values[i] = z.Values[i] + beta * p.Values[i];
            }
            return -1;
        }

        public static int GMRES(SparseMatrixCSR A, DenseVector b, DenseVector x,
            IPreconditioner M, double tol, int maxit, int m, IterationCallback cb)
        {
            int n = A.Nrows;
            double bnorm = b.Norm2(); if (bnorm < 1e-30) bnorm = 1.0;
            DenseVector r = new DenseVector(n);
            DenseVector w = new DenseVector(n);
            DenseVector[] V = new DenseVector[m + 1];
            for (int i = 0; i <= m; i++) V[i] = new DenseVector(n);
            double[,] H = new double[m + 1, m];
            double[] cs = new double[m], sn = new double[m], g = new double[m + 1];
            DenseVector Mw = new DenseVector(n);
            int totalIt = 0;
            while (totalIt < maxit)
            {
                A.Multiply(x, r);
                for (int i = 0; i < n; i++) r.Values[i] = b.Values[i] - r.Values[i];
                M.Apply(r, Mw);
                double beta = Mw.Norm2();
                if (beta / bnorm < tol) return totalIt;
                for (int i = 0; i < n; i++) V[0].Values[i] = Mw.Values[i] / beta;
                Array.Clear(g, 0, g.Length); g[0] = beta;
                Array.Clear(H, 0, H.Length); Array.Clear(cs, 0, cs.Length); Array.Clear(sn, 0, sn.Length);
                int j = 0; bool conv = false;
                for (; j < m && totalIt < maxit; j++, totalIt++)
                {
                    A.Multiply(V[j], w); M.Apply(w, Mw); w.CopyFrom(Mw);
                    for (int i = 0; i <= j; i++)
                    {
                        H[i, j] = w.Dot(V[i]);
                        w.Axpy(-H[i, j], V[i]);
                    }
                    H[j + 1, j] = w.Norm2();
                    if (H[j + 1, j] != 0.0)
                        for (int i = 0; i < n; i++) V[j + 1].Values[i] = w.Values[i] / H[j + 1, j];
                    for (int i = 0; i < j; i++)
                    {
                        double h1 = cs[i] * H[i, j] + sn[i] * H[i + 1, j];
                        double h2 = -sn[i] * H[i, j] + cs[i] * H[i + 1, j];
                        H[i, j] = h1; H[i + 1, j] = h2;
                    }
                    double denom = Math.Sqrt(H[j, j] * H[j, j] + H[j + 1, j] * H[j + 1, j]);
                    if (denom < 1e-30) { cs[j] = 1; sn[j] = 0; }
                    else { cs[j] = H[j, j] / denom; sn[j] = H[j + 1, j] / denom; }
                    H[j, j] = cs[j] * H[j, j] + sn[j] * H[j + 1, j];
                    H[j + 1, j] = 0.0;
                    double g0 = g[j];
                    g[j] = cs[j] * g0;
                    g[j + 1] = -sn[j] * g0;
                    double res = Math.Abs(g[j + 1]) / bnorm;
                    if (cb != null) cb(totalIt + 1, res);
                    if (res < tol) { conv = true; j++; break; }
                }
                double[] y = new double[j];
                for (int i = j - 1; i >= 0; i--)
                {
                    double s = g[i];
                    for (int k = i + 1; k < j; k++) s -= H[i, k] * y[k];
                    y[i] = s / H[i, i];
                }
                for (int i = 0; i < j; i++) x.Axpy(y[i], V[i]);
                if (conv) return totalIt;
            }
            return -1;
        }

        public static int BiCGSTAB(SparseMatrixCSR A, DenseVector b, DenseVector x,
            IPreconditioner M, double tol, int maxit, IterationCallback cb)
        {
            int n = A.Nrows;
            double bnorm = b.Norm2(); if (bnorm < 1e-30) bnorm = 1.0;
            DenseVector r = new DenseVector(n);
            DenseVector rHat = new DenseVector(n);
            DenseVector p = new DenseVector(n);
            DenseVector v = new DenseVector(n);
            DenseVector s = new DenseVector(n);
            DenseVector t = new DenseVector(n);
            DenseVector Mp = new DenseVector(n);
            DenseVector Ms = new DenseVector(n);
            A.Multiply(x, r);
            for (int i = 0; i < n; i++) { r.Values[i] = b.Values[i] - r.Values[i]; rHat.Values[i] = r.Values[i]; }
            double rho = 1.0, alpha = 1.0, omega = 1.0;
            p.SetZero(); v.SetZero();
            for (int it = 1; it <= maxit; it++)
            {
                double rhoNew = rHat.Dot(r);
                if (Math.Abs(rhoNew) < 1e-30) return -1;
                double beta = (rhoNew / rho) * (alpha / omega);
                rho = rhoNew;
                for (int i = 0; i < n; i++) p.Values[i] = r.Values[i] + beta * (p.Values[i] - omega * v.Values[i]);
                M.Apply(p, Mp);
                A.Multiply(Mp, v);
                double rHv = rHat.Dot(v);
                if (Math.Abs(rHv) < 1e-30) return -1;
                alpha = rho / rHv;
                for (int i = 0; i < n; i++) s.Values[i] = r.Values[i] - alpha * v.Values[i];
                double sNorm = s.Norm2();
                if (sNorm / bnorm < tol) { x.Axpy(alpha, Mp); if (cb != null) cb(it, sNorm / bnorm); return it; }
                M.Apply(s, Ms);
                A.Multiply(Ms, t);
                double tt = t.Dot(t);
                if (tt < 1e-30) omega = 0; else omega = t.Dot(s) / tt;
                for (int i = 0; i < n; i++) x.Values[i] += alpha * Mp.Values[i] + omega * Ms.Values[i];
                for (int i = 0; i < n; i++) r.Values[i] = s.Values[i] - omega * t.Values[i];
                double res = r.Norm2() / bnorm;
                if (cb != null) cb(it, res);
                if (res < tol) return it;
                if (Math.Abs(omega) < 1e-30) return -1;
            }
            return -1;
        }
    }
}
