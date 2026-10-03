using System;

namespace FEM2D.Elements
{
    public sealed class LT6Element : FiniteElement
    {
        public override ElementType Type => ElementType.LT6;
        public override int NodesPerElement => 6;
        public override int DofPerNode => 2;
        public override int NumGaussPoints => 3;

        private static readonly double[] _L1 = new double[] { 0.5, 0.5, 0.0 };
        private static readonly double[] _L2 = new double[] { 0.0, 0.5, 0.5 };
        private static readonly double[] _L3 = new double[] { 0.5, 0.0, 0.5 };
        private static readonly double[] _W = new double[] { 1.0/6.0, 1.0/6.0, 1.0/6.0 };

        private static double[,] DMatrix(double E, double nu)
        {
            double c = E/((1.0+nu)*(1.0-2.0*nu));
            double[,] D = new double[3,3];
            D[0,0]=c*(1.0-nu); D[0,1]=c*nu; D[0,2]=0.0;
            D[1,0]=c*nu; D[1,1]=c*(1.0-nu); D[1,2]=0.0;
            D[2,0]=0.0; D[2,1]=0.0; D[2,2]=c*(0.5-nu);
            return D;
        }

        private static void ShapeLocal(double L1,double L2,double L3,double[] N,double[] dNdL1,double[] dNdL2,double[] dNdL3)
        {
            N[0]=L1*(2.0*L1-1.0); N[1]=L2*(2.0*L2-1.0); N[2]=L3*(2.0*L3-1.0);
            N[3]=4.0*L1*L2; N[4]=4.0*L2*L3; N[5]=4.0*L3*L1;
            dNdL1[0]=4.0*L1-1.0; dNdL2[0]=0.0; dNdL3[0]=0.0;
            dNdL1[1]=0.0; dNdL2[1]=4.0*L2-1.0; dNdL3[1]=0.0;
            dNdL1[2]=0.0; dNdL2[2]=0.0; dNdL3[2]=4.0*L3-1.0;
            dNdL1[3]=4.0*L2; dNdL2[3]=4.0*L1; dNdL3[3]=0.0;
            dNdL1[4]=0.0; dNdL2[4]=4.0*L3; dNdL3[4]=4.0*L2;
            dNdL1[5]=4.0*L3; dNdL2[5]=0.0; dNdL3[5]=4.0*L1;
        }

        private double Jacobian(double[,] coords, double[] dNdL1, double[] dNdL2, double[] dNdL3, double[] dNdx, double[] dNdy)
        {
            double x1=coords[0,0],y1=coords[0,1],x2=coords[1,0],y2=coords[1,1],x3=coords[2,0],y3=coords[2,1];
            double dxdl1=x1-x3, dxdl2=x2-x3, dydl1=y1-y3, dydl2=y2-y3;
            double detJ=dxdl1*dydl2-dxdl2*dydl1;
            double inv=1.0/detJ;
            double dldx1=dydl2*inv, dldx2=-dydl1*inv, dldy1=-dxdl2*inv, dldy2=dxdl1*inv;
            for(int i=0;i<6;i++)
            {
                dNdx[i]=dNdL1[i]*dldx1+dNdL2[i]*dldx2+dNdL3[i]*(-dldx1-dldx2);
                dNdy[i]=dNdL1[i]*dldy1+dNdL2[i]*dldy2+dNdL3[i]*(-dldy1-dldy2);
            }
            return detJ;
        }

        public override double ComputeArea(double[,] coords)
        {
            double x1=coords[0,0],y1=coords[0,1],x2=coords[1,0],y2=coords[1,1],x3=coords[2,0],y3=coords[2,1];
            return 0.5*Math.Abs((x2-x1)*(y3-y1)-(x3-x1)*(y2-y1));
        }

        public override void ComputeStiffness(double[,] coords, double E, double nu, double thickness, double[,] ke)
        {
            double[,] D=DMatrix(E,nu); Array.Clear(ke,0,ke.Length);
            double[] dNdL1=new double[6],dNdL2=new double[6],dNdL3=new double[6],dNdx=new double[6],dNdy=new double[6],N=new double[6];
            double[,] B=new double[3,12];
            for(int gp=0;gp<3;gp++)
            {
                double L1=_L1[gp],L2=_L2[gp],L3=_L3[gp];
                ShapeLocal(L1,L2,L3,N,dNdL1,dNdL2,dNdL3);
                double detJ=Jacobian(coords,dNdL1,dNdL2,dNdL3,dNdx,dNdy);
                double w=_W[gp]*detJ*2.0*thickness;
                Array.Clear(B,0,B.Length);
                for(int i=0;i<6;i++){B[0,2*i]=dNdx[i];B[1,2*i+1]=dNdy[i];B[2,2*i]=dNdy[i];B[2,2*i+1]=dNdx[i];}
                for(int i=0;i<12;i++)for(int j=0;j<12;j++){double s=0.0;for(int k=0;k<3;k++)for(int m=0;m<3;m++)s+=B[k,i]*D[k,m]*B[m,j];ke[i,j]+=s*w;}
            }
        }

        public override void ComputeMass(double[,] coords, double rho, double thickness, double[,] me)
        {
            Array.Clear(me,0,me.Length);
            double[] N=new double[6],dNdL1=new double[6],dNdL2=new double[6],dNdL3=new double[6],dNdx=new double[6],dNdy=new double[6];
            double[] qL1=new double[]{1.0/3.0,0.6,0.2,0.2}; double[] qL2=new double[]{1.0/3.0,0.2,0.6,0.2};
            double[] qL3=new double[]{1.0/3.0,0.2,0.2,0.6}; double[] qW=new double[]{-27.0/48.0,25.0/48.0,25.0/48.0,25.0/48.0};
            for(int gp=0;gp<4;gp++)
            {
                ShapeLocal(qL1[gp],qL2[gp],qL3[gp],N,dNdL1,dNdL2,dNdL3);
                double detJ=Jacobian(coords,dNdL1,dNdL2,dNdL3,dNdx,dNdy);
                double w=rho*thickness*qW[gp]*detJ*2.0;
                for(int a=0;a<6;a++)for(int b=0;b<6;b++){double v=N[a]*N[b]*w;me[2*a,2*b]+=v;me[2*a+1,2*b+1]+=v;}
            }
        }

        public override void ComputeConvection(double[,] coords, double[] uNodal, double thickness, double[,] ce)
        {
            Array.Clear(ce,0,ce.Length);
            double[] N=new double[6],dNdL1=new double[6],dNdL2=new double[6],dNdL3=new double[6],dNdx=new double[6],dNdy=new double[6];
            for(int gp=0;gp<3;gp++)
            {
                ShapeLocal(_L1[gp],_L2[gp],_L3[gp],N,dNdL1,dNdL2,dNdL3);
                double detJ=Jacobian(coords,dNdL1,dNdL2,dNdL3,dNdx,dNdy);
                double w=_W[gp]*detJ*2.0*thickness;
                double u=0.0,v=0.0;
                for(int i=0;i<6;i++){u+=N[i]*uNodal[2*i];v+=N[i]*uNodal[2*i+1];}
                for(int a=0;a<6;a++)for(int b=0;b<6;b++){double val=N[a]*(u*dNdx[b]+v*dNdy[b])*w;ce[2*a,2*b]+=val;ce[2*a+1,2*b+1]+=val;}
            }
        }

        public override void ComputeViscous(double[,] coords, double mu, double thickness, double[,] ve)
        {
            Array.Clear(ve,0,ve.Length);
            double[] N=new double[6],dNdL1=new double[6],dNdL2=new double[6],dNdL3=new double[6],dNdx=new double[6],dNdy=new double[6];
            for(int gp=0;gp<3;gp++)
            {
                ShapeLocal(_L1[gp],_L2[gp],_L3[gp],N,dNdL1,dNdL2,dNdL3);
                double detJ=Jacobian(coords,dNdL1,dNdL2,dNdL3,dNdx,dNdy);
                double w=mu*thickness*_W[gp]*detJ*2.0;
                for(int a=0;a<6;a++)for(int b=0;b<6;b++){double v=(dNdx[a]*dNdx[b]+dNdy[a]*dNdy[b])*w;ve[2*a,2*b]+=v;ve[2*a+1,2*b+1]+=v;}
            }
        }

        public override void ComputePressureGradient(double[,] coords, double thickness, double[,] gx, double[,] gy)
        {
            Array.Clear(gx,0,gx.Length); Array.Clear(gy,0,gy.Length);
            double[] Nv=new double[6],dNdL1=new double[6],dNdL2=new double[6],dNdL3=new double[6],dNdx=new double[6],dNdy=new double[6];
            double[] Np=new double[3];
            for(int gp=0;gp<3;gp++)
            {
                double L1=_L1[gp],L2=_L2[gp],L3=_L3[gp];
                ShapeLocal(L1,L2,L3,Nv,dNdL1,dNdL2,dNdL3);
                double detJ=Jacobian(coords,dNdL1,dNdL2,dNdL3,dNdx,dNdy);
                double w=_W[gp]*detJ*2.0*thickness;
                Np[0]=L1;Np[1]=L2;Np[2]=L3;
                for(int a=0;a<6;a++)for(int p=0;p<3;p++){gx[2*a,p]+=dNdx[a]*Np[p]*w;gy[2*a+1,p]+=dNdy[a]*Np[p]*w;}
            }
        }

        public override void ComputeStress(double[,] coords, double[] ue, double E, double nu, out double sxx, out double syy, out double sxy, out double von)
        {
            double[,] D=DMatrix(E,nu);
            double[] N=new double[6],dNdL1=new double[6],dNdL2=new double[6],dNdL3=new double[6],dNdx=new double[6],dNdy=new double[6];
            ShapeLocal(1.0/3.0,1.0/3.0,1.0/3.0,N,dNdL1,dNdL2,dNdL3);
            Jacobian(coords,dNdL1,dNdL2,dNdL3,dNdx,dNdy);
            double exx=0.0,eyy=0.0,gxy=0.0;
            for(int i=0;i<6;i++){exx+=dNdx[i]*ue[2*i];eyy+=dNdy[i]*ue[2*i+1];gxy+=dNdy[i]*ue[2*i]+dNdx[i]*ue[2*i+1];}
            sxx=D[0,0]*exx+D[0,1]*eyy; syy=D[1,0]*exx+D[1,1]*eyy; sxy=D[2,2]*gxy;
            von=Math.Sqrt(sxx*sxx-sxx*syy+syy*syy+3.0*sxy*sxy);
        }

        public override void ShapeFunctions(double[,] coords, double xi, double eta, double[] N, double[] dNdx, double[] dNdy, out double detJ)
        {
            double L1=xi,L2=eta,L3=1.0-xi-eta;
            double[] dNdL1=new double[6],dNdL2=new double[6],dNdL3=new double[6];
            ShapeLocal(L1,L2,L3,N,dNdL1,dNdL2,dNdL3);
            detJ=Jacobian(coords,dNdL1,dNdL2,dNdL3,dNdx,dNdy);
        }
    }
}
