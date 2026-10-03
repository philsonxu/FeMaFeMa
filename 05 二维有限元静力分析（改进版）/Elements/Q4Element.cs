using System;

namespace FEM2D.Elements
{
    public sealed class Q4Element : FiniteElement
    {
        public override ElementType Type => ElementType.Q4;
        public override int NodesPerElement => 4;
        public override int DofPerNode => 2;
        public override int NumGaussPoints => 4;

        private static readonly double[] _GP = new double[] { -1.0/Math.Sqrt(3.0), 1.0/Math.Sqrt(3.0) };

        private static double[,] DMatrix(double E, double nu)
        {
            double c = E/((1.0+nu)*(1.0-2.0*nu));
            double[,] D = new double[3,3];
            D[0,0]=c*(1.0-nu); D[0,1]=c*nu; D[0,2]=0.0;
            D[1,0]=c*nu; D[1,1]=c*(1.0-nu); D[1,2]=0.0;
            D[2,0]=0.0; D[2,1]=0.0; D[2,2]=c*(0.5-nu);
            return D;
        }

        private static void ShapeNat(double xi,double eta,double[] N,double[] dNdxi,double[] dNdeta)
        {
            N[0]=0.25*(1.0-xi)*(1.0-eta); N[1]=0.25*(1.0+xi)*(1.0-eta);
            N[2]=0.25*(1.0+xi)*(1.0+eta); N[3]=0.25*(1.0-xi)*(1.0+eta);
            dNdxi[0]=-0.25*(1.0-eta); dNdeta[0]=-0.25*(1.0-xi);
            dNdxi[1]=0.25*(1.0-eta); dNdeta[1]=-0.25*(1.0+xi);
            dNdxi[2]=0.25*(1.0+eta); dNdeta[2]=0.25*(1.0+xi);
            dNdxi[3]=-0.25*(1.0+eta); dNdeta[3]=0.25*(1.0-xi);
        }

        private static double Jacobian(double[,] coords, double[] dNdxi, double[] dNdeta, double[] dNdx, double[] dNdy)
        {
            double dxdxi=0.0,dxdeta=0.0,dydxi=0.0,dydeta=0.0;
            for(int i=0;i<4;i++){dxdxi+=dNdxi[i]*coords[i,0];dydxi+=dNdxi[i]*coords[i,1];dxdeta+=dNdeta[i]*coords[i,0];dydeta+=dNdeta[i]*coords[i,1];}
            double detJ=dxdxi*dydeta-dxdeta*dydxi; double inv=1.0/detJ;
            for(int i=0;i<4;i++){dNdx[i]=(dydeta*dNdxi[i]-dydxi*dNdeta[i])*inv;dNdy[i]=(-dxdeta*dNdxi[i]+dxdxi*dNdeta[i])*inv;}
            return detJ;
        }

        public override double ComputeArea(double[,] coords)
        {
            double area=0.0; double[] N=new double[4],dNdxi=new double[4],dNdeta=new double[4],dNdx=new double[4],dNdy=new double[4];
            for(int i=0;i<2;i++)for(int j=0;j<2;j++){ShapeNat(_GP[i],_GP[j],N,dNdxi,dNdeta);double detJ=Jacobian(coords,dNdxi,dNdeta,dNdx,dNdy);area+=detJ;}
            return area;
        }

        public override void ComputeStiffness(double[,] coords, double E, double nu, double thickness, double[,] ke)
        {
            double[,] D=DMatrix(E,nu); Array.Clear(ke,0,ke.Length);
            double[] N=new double[4],dNdxi=new double[4],dNdeta=new double[4],dNdx=new double[4],dNdy=new double[4];
            double[,] B=new double[3,8];
            for(int i=0;i<2;i++)for(int j=0;j<2;j++)
            {
                ShapeNat(_GP[i],_GP[j],N,dNdxi,dNdeta);
                double detJ=Jacobian(coords,dNdxi,dNdeta,dNdx,dNdy);
                double w=detJ*thickness;
                Array.Clear(B,0,B.Length);
                for(int k=0;k<4;k++){B[0,2*k]=dNdx[k];B[1,2*k+1]=dNdy[k];B[2,2*k]=dNdy[k];B[2,2*k+1]=dNdx[k];}
                for(int a=0;a<8;a++)for(int b=0;b<8;b++){double s=0.0;for(int k=0;k<3;k++)for(int m=0;m<3;m++)s+=B[k,a]*D[k,m]*B[m,b];ke[a,b]+=s*w;}
            }
        }

        public override void ComputeMass(double[,] coords, double rho, double thickness, double[,] me)
        {
            Array.Clear(me,0,me.Length);
            double[] N=new double[4],dNdxi=new double[4],dNdeta=new double[4],dNdx=new double[4],dNdy=new double[4];
            for(int i=0;i<2;i++)for(int j=0;j<2;j++)
            {
                ShapeNat(_GP[i],_GP[j],N,dNdxi,dNdeta);
                double detJ=Jacobian(coords,dNdxi,dNdeta,dNdx,dNdy);
                double w=rho*thickness*detJ;
                for(int a=0;a<4;a++)for(int b=0;b<4;b++){double v=N[a]*N[b]*w;me[2*a,2*b]+=v;me[2*a+1,2*b+1]+=v;}
            }
        }

        public override void ComputeConvection(double[,] coords, double[] uNodal, double thickness, double[,] ce)
        {
            Array.Clear(ce,0,ce.Length);
            double[] N=new double[4],dNdxi=new double[4],dNdeta=new double[4],dNdx=new double[4],dNdy=new double[4];
            for(int i=0;i<2;i++)for(int j=0;j<2;j++)
            {
                ShapeNat(_GP[i],_GP[j],N,dNdxi,dNdeta);
                double detJ=Jacobian(coords,dNdxi,dNdeta,dNdx,dNdy);
                double w=thickness*detJ;
                double u=0.0,v=0.0;
                for(int k=0;k<4;k++){u+=N[k]*uNodal[2*k];v+=N[k]*uNodal[2*k+1];}
                for(int a=0;a<4;a++)for(int b=0;b<4;b++){double val=N[a]*(u*dNdx[b]+v*dNdy[b])*w;ce[2*a,2*b]+=val;ce[2*a+1,2*b+1]+=val;}
            }
        }

        public override void ComputeViscous(double[,] coords, double mu, double thickness, double[,] ve)
        {
            Array.Clear(ve,0,ve.Length);
            double[] N=new double[4],dNdxi=new double[4],dNdeta=new double[4],dNdx=new double[4],dNdy=new double[4];
            for(int i=0;i<2;i++)for(int j=0;j<2;j++)
            {
                ShapeNat(_GP[i],_GP[j],N,dNdxi,dNdeta);
                double detJ=Jacobian(coords,dNdxi,dNdeta,dNdx,dNdy);
                double w=mu*thickness*detJ;
                for(int a=0;a<4;a++)for(int b=0;b<4;b++){double v=(dNdx[a]*dNdx[b]+dNdy[a]*dNdy[b])*w;ve[2*a,2*b]+=v;ve[2*a+1,2*b+1]+=v;}
            }
        }

        public override void ComputePressureGradient(double[,] coords, double thickness, double[,] gx, double[,] gy)
        {
            Array.Clear(gx,0,gx.Length); Array.Clear(gy,0,gy.Length);
            double[] N=new double[4],dNdxi=new double[4],dNdeta=new double[4],dNdx=new double[4],dNdy=new double[4];
            for(int i=0;i<2;i++)for(int j=0;j<2;j++)
            {
                ShapeNat(_GP[i],_GP[j],N,dNdxi,dNdeta);
                double detJ=Jacobian(coords,dNdxi,dNdeta,dNdx,dNdy);
                double w=thickness*detJ;
                for(int a=0;a<4;a++)for(int p=0;p<4;p++){gx[2*a,p]+=dNdx[a]*N[p]*w;gy[2*a+1,p]+=dNdy[a]*N[p]*w;}
            }
        }

        public override void ComputeStress(double[,] coords, double[] ue, double E, double nu, out double sxx, out double syy, out double sxy, out double von)
        {
            double[,] D=DMatrix(E,nu);
            double[] N=new double[4],dNdxi=new double[4],dNdeta=new double[4],dNdx=new double[4],dNdy=new double[4];
            ShapeNat(0.0,0.0,N,dNdxi,dNdeta);
            Jacobian(coords,dNdxi,dNdeta,dNdx,dNdy);
            double exx=0.0,eyy=0.0,gxy=0.0;
            for(int i=0;i<4;i++){exx+=dNdx[i]*ue[2*i];eyy+=dNdy[i]*ue[2*i+1];gxy+=dNdy[i]*ue[2*i]+dNdx[i]*ue[2*i+1];}
            sxx=D[0,0]*exx+D[0,1]*eyy; syy=D[1,0]*exx+D[1,1]*eyy; sxy=D[2,2]*gxy;
            von=Math.Sqrt(sxx*sxx-sxx*syy+syy*syy+3.0*sxy*sxy);
        }

        public override void ShapeFunctions(double[,] coords, double xi, double eta, double[] N, double[] dNdx, double[] dNdy, out double detJ)
        {
            double[] dNdxi=new double[4],dNdeta=new double[4];
            ShapeNat(xi,eta,N,dNdxi,dNdeta);
            detJ=Jacobian(coords,dNdxi,dNdeta,dNdx,dNdy);
        }
    }
}
