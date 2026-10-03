using System;

namespace FEM2D.Elements
{
    public sealed class CST3Element : FiniteElement
    {
        public override ElementType Type => ElementType.CST3;
        public override int NodesPerElement => 3;
        public override int DofPerNode => 2;
        public override int NumGaussPoints => 1;

        public override double ComputeArea(double[,] coords)
        {
            double x1 = coords[0,0], y1 = coords[0,1];
            double x2 = coords[1,0], y2 = coords[1,1];
            double x3 = coords[2,0], y3 = coords[2,1];
            return 0.5 * Math.Abs((x2-x1)*(y3-y1) - (x3-x1)*(y2-y1));
        }

        private void BMatrix(double[,] coords, double[,] B, out double A)
        {
            double x1=coords[0,0],y1=coords[0,1],x2=coords[1,0],y2=coords[1,1],x3=coords[2,0],y3=coords[2,1];
            double b1=y2-y3,b2=y3-y1,b3=y1-y2;
            double c1=x3-x2,c2=x1-x3,c3=x2-x1;
            A=0.5*Math.Abs(b1*c2-b2*c1);
            double i2A=1.0/(2.0*A);
            Array.Clear(B,0,B.Length);
            B[0,0]=b1*i2A; B[1,1]=c1*i2A; B[2,0]=c1*i2A; B[2,1]=b1*i2A;
            B[0,2]=b2*i2A; B[1,3]=c2*i2A; B[2,2]=c2*i2A; B[2,3]=b2*i2A;
            B[0,4]=b3*i2A; B[1,5]=c3*i2A; B[2,4]=c3*i2A; B[2,5]=b3*i2A;
        }

        private static void DMatrix(double E, double nu, double[,] D)
        {
            double c = E/((1.0+nu)*(1.0-2.0*nu));
            D[0,0]=c*(1.0-nu); D[0,1]=c*nu; D[0,2]=0.0;
            D[1,0]=c*nu; D[1,1]=c*(1.0-nu); D[1,2]=0.0;
            D[2,0]=0.0; D[2,1]=0.0; D[2,2]=c*(0.5-nu);
        }

        public override void ComputeStiffness(double[,] coords, double E, double nu, double thickness, double[,] ke)
        {
            double[,] B=new double[3,6], DM=new double[3,3];
            double A; BMatrix(coords,B,out A); DMatrix(E,nu,DM);
            double[,] DB=new double[3,6];
            for(int i=0;i<3;i++)for(int j=0;j<6;j++){double s=0.0;for(int k=0;k<3;k++)s+=DM[i,k]*B[k,j];DB[i,j]=s;}
            double tA=thickness*A;
            Array.Clear(ke,0,ke.Length);
            for(int i=0;i<6;i++)for(int j=0;j<6;j++){double s=0.0;for(int k=0;k<3;k++)s+=B[k,i]*DB[k,j];ke[i,j]=s*tA;}
        }

        public override void ComputeMass(double[,] coords, double rho, double thickness, double[,] me)
        {
            double A=ComputeArea(coords); double v=rho*thickness*A/12.0;
            Array.Clear(me,0,me.Length);
            for(int i=0;i<3;i++)for(int j=0;j<3;j++){double n=(i==j)?2.0:1.0;me[2*i,2*j]=v*n;me[2*i+1,2*j+1]=v*n;}
        }

        public override void ComputeConvection(double[,] coords, double[] uNodal, double thickness, double[,] ce)
        {
            double A=ComputeArea(coords); double i2A=1.0/(2.0*A);
            double x1=coords[0,0],y1=coords[0,1],x2=coords[1,0],y2=coords[1,1],x3=coords[2,0],y3=coords[2,1];
            double[] b=new double[]{y2-y3,y3-y1,y1-y2}; double[] c=new double[]{x3-x2,x1-x3,x2-x1};
            double[] dNdx=new double[3],dNdy=new double[3];
            for(int i=0;i<3;i++){dNdx[i]=b[i]*i2A;dNdy[i]=c[i]*i2A;}
            double uB=(uNodal[0]+uNodal[2]+uNodal[4])/3.0; double vB=(uNodal[1]+uNodal[3]+uNodal[5])/3.0;
            double w=thickness*A; Array.Clear(ce,0,ce.Length);
            for(int a=0;a<3;a++)for(int bb=0;bb<3;bb++){double v=(1.0/3.0)*(uB*dNdx[bb]+vB*dNdy[bb])*w;ce[2*a,2*bb]=v;ce[2*a+1,2*bb+1]=v;}
        }

        public override void ComputeViscous(double[,] coords, double mu, double thickness, double[,] ve)
        {
            double A=ComputeArea(coords); double i2A=1.0/(2.0*A);
            double x1=coords[0,0],y1=coords[0,1],x2=coords[1,0],y2=coords[1,1],x3=coords[2,0],y3=coords[2,1];
            double[] b=new double[]{y2-y3,y3-y1,y1-y2}; double[] c=new double[]{x3-x2,x1-x3,x2-x1};
            Array.Clear(ve,0,ve.Length);
            for(int a=0;a<3;a++)for(int bb=0;bb<3;bb++){double v=mu*thickness*0.25*(b[a]*b[bb]+c[a]*c[bb])/A;ve[2*a,2*bb]=v;ve[2*a+1,2*bb+1]=v;}
        }

        public override void ComputePressureGradient(double[,] coords, double thickness, double[,] gx, double[,] gy)
        {
            double A=ComputeArea(coords); double i2A=1.0/(2.0*A);
            double x1=coords[0,0],y1=coords[0,1],x2=coords[1,0],y2=coords[1,1],x3=coords[2,0],y3=coords[2,1];
            double[] b=new double[]{y2-y3,y3-y1,y1-y2}; double[] c=new double[]{x3-x2,x1-x3,x2-x1};
            double w=thickness*A/3.0; Array.Clear(gx,0,gx.Length); Array.Clear(gy,0,gy.Length);
            for(int a=0;a<3;a++)for(int p=0;p<3;p++)if(a==p){gx[2*a,p]=w*b[a]*i2A;gy[2*a+1,p]=w*c[a]*i2A;}
        }

        public override void ComputeStress(double[,] coords, double[] ue, double E, double nu, out double sxx, out double syy, out double sxy, out double von)
        {
            double[,] B=new double[3,6]; double A; BMatrix(coords,B,out A);
            double[,] D=new double[3,3]; DMatrix(E,nu,D);
            double[] eps=new double[3];
            for(int i=0;i<3;i++){double s=0.0;for(int j=0;j<6;j++)s+=B[i,j]*ue[j];eps[i]=s;}
            sxx=D[0,0]*eps[0]+D[0,1]*eps[1]; syy=D[1,0]*eps[0]+D[1,1]*eps[1]; sxy=D[2,2]*eps[2];
            von=Math.Sqrt(sxx*sxx-sxx*syy+syy*syy+3.0*sxy*sxy);
        }

        public override void ShapeFunctions(double[,] coords, double xi, double eta, double[] N, double[] dNdx, double[] dNdy, out double detJ)
        {
            N[0]=xi; N[1]=eta; N[2]=1.0-xi-eta;
            double A=ComputeArea(coords); double i2A=1.0/(2.0*A);
            double x1=coords[0,0],y1=coords[0,1],x2=coords[1,0],y2=coords[1,1],x3=coords[2,0],y3=coords[2,1];
            dNdx[0]=(y2-y3)*i2A; dNdy[0]=(x3-x2)*i2A;
            dNdx[1]=(y3-y1)*i2A; dNdy[1]=(x1-x3)*i2A;
            dNdx[2]=(y1-y2)*i2A; dNdy[2]=(x2-x1)*i2A;
            detJ=2.0*A;
        }
    }
}
