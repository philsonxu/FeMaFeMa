using System;

namespace NonlinearFEM2D.Materials
{
    /// <summary>
    /// J2 von Mises 弹塑性本构：各向同性线性硬化 + 关联流动 + 平面应变。
    /// 应力应变约定 (3 分量, 平面应变): [σxx, σyy, σzz, τxy], 应变同理。
    /// 向后 Euler 返回映射 (Simo-Hughes)。
    /// </summary>
    public sealed class J2PlasticityMaterial
    {
        public double E;            // Young
        public double Nu;           // Poisson
        public double SigmaY0;      // 初始屈服应力
        public double Hiso;         // 等向硬化模量
        public double Lambda;
        public double Mu;
        public double Kbulk;

        public J2PlasticityMaterial(double e, double nu, double sy0, double h)
        {
            this.E = e; this.Nu = nu; this.SigmaY0 = sy0; this.Hiso = h;
            this.Mu = e / (2.0 * (1.0 + nu));
            this.Lambda = e * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
            this.Kbulk = this.Lambda + 2.0 * this.Mu / 3.0;
        }

        // 将总应变增量 deps 更新到应力 sig 与累积塑性应变 eqp；返回一致弹塑性切线 Dalg(4x4)
        // 调用方在 sig, eqp 中提供当前已知历史，返回后就地更新。
        public void UpdateStressAndTangent(double[] sig, double[] deps, ref double eqp, double[,] Dalg)
        {
            double[,] De = new double[4, 4];
            double mu = this.Mu;
            double lam = this.Lambda;
            De[0, 0] = lam + 2.0 * mu; De[0, 1] = lam; De[0, 2] = lam;
            De[1, 0] = lam; De[1, 1] = lam + 2.0 * mu; De[1, 2] = lam;
            De[2, 0] = lam; De[2, 1] = lam; De[2, 2] = lam + 2.0 * mu;
            De[3, 3] = mu;
            // 弹性试探
            double[] sigT = new double[4];
            for (int i = 0; i < 4; i++)
            {
                double s = sig[i];
                for (int j = 0; j < 4; j++) s += De[i, j] * deps[j];
                sigT[i] = s;
            }
            double sigm = (sigT[0] + sigT[1] + sigT[2]) / 3.0;
            double[] sT = new double[4] { sigT[0] - sigm, sigT[1] - sigm, sigT[2] - sigm, sigT[3] };
            double sTnorm = Math.Sqrt(sT[0] * sT[0] + sT[1] * sT[1] + sT[2] * sT[2] + 2.0 * sT[3] * sT[3]);
            double sigmaY = this.SigmaY0 + this.Hiso * eqp;
            double fTrial = sTnorm - Math.Sqrt(2.0 / 3.0) * sigmaY;
            // 默认 Dalg = De
            for (int i = 0; i < 4; i++) for (int j = 0; j < 4; j++) Dalg[i, j] = De[i, j];
            if (fTrial <= 1.0e-12)
            {
                // 弹性
                for (int i = 0; i < 4; i++) sig[i] = sigT[i];
                return;
            }
            // 塑性：径向返回
            double dgamma = fTrial / (2.0 * mu + 2.0 * this.Hiso / 3.0);
            double nf = (sTnorm > 1.0e-30) ? sTnorm : 1.0;
            double[] n = new double[4];
            for (int i = 0; i < 4; i++) n[i] = sT[i] / nf;
            for (int i = 0; i < 4; i++) sig[i] = sigT[i] - 2.0 * mu * dgamma * n[i];
            eqp += Math.Sqrt(2.0 / 3.0) * dgamma;
            // 一致切线
            double beta = 2.0 * mu * dgamma / sTnorm;
            double gamma_fac = 1.0 / (1.0 + this.Hiso / (3.0 * mu));
            for (int i = 0; i < 4; i++) for (int j = 0; j < 4; j++)
            {
                double I4 = (i == j) ? 1.0 : 0.0;
                double I4dev = I4 - ((i < 3 && j < 3) ? 1.0 / 3.0 : 0.0);
                // 注意：对 4 分量 I4dev 对剪分量(3)应为 0.5 形式的偏量算子? 对 4 分量 (xx,yy,zz,xy)，
                // 偏量投影 P[i,j] = I4 - (i<3 && j<3 ? 1/3 : 0) 正好给出偏量（剪分量 2*τ 已在 sTnorm 以 2*s3^2 考虑）
                // 修正剪分量偏差投影：i==3 且 j==3 时是 0.5；跨角元剪-正为0；正-剪为0；正-正为δ_ij-1/3。
                double Pij;
                if (i < 3 && j < 3) Pij = ((i == j) ? 1.0 : 0.0) - 1.0 / 3.0;
                else if (i == 3 && j == 3) Pij = 0.5;
                else Pij = 0.0;
                double nnij = n[i] * n[j];
                Dalg[i, j] = 2.0 * mu * ((1.0 - beta) * Pij + gamma_fac * nnij)
                           + this.Kbulk * ((i < 3 && j < 3) ? 1.0 : 0.0);
            }
        }

        public static double VonMises(double[] sig)
        {
            double sxx = sig[0], syy = sig[1], szz = sig[2], sxy = sig[3];
            double s1 = sxx - syy, s2 = syy - szz, s3 = szz - sxx;
            double v = 0.5 * (s1 * s1 + s2 * s2 + s3 * s3) + 3.0 * sxy * sxy;
            return Math.Sqrt(Math.Max(0.0, v));
        }
    }
}
