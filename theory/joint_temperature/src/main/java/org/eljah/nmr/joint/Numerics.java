package org.eljah.nmr.joint;

import java.util.Arrays;
import java.util.function.DoubleUnaryOperator;

/** Small deterministic numerical routines; no third-party dependencies. */
final class Numerics {
    private Numerics() {}
    record Eigen(double[] values, double[][] vectors) {}
    record RankFit(double[] singular, double[][] vectors, double[][] fitted,
                   double sse, double maximum, double retainedFraction) {}
    record Regression(double intercept, double slope, double sse, int n) {}
    record LeastSquares(double[] coefficients, double condition) {}
    record Minimum(double x, double value, boolean boundary) {}

    static double dot(double[] a, double[] b) {
        double s=0; for (int i=0;i<a.length;i++) s+=a[i]*b[i]; return s;
    }
    static double sumSquares(double[] a) { return dot(a,a); }
    static double[][] copy(double[][] a) {
        double[][] r=new double[a.length][]; for(int i=0;i<a.length;i++) r[i]=a[i].clone(); return r;
    }

    /** Symmetric Jacobi rotations; eigenvectors are columns, sorted descending. */
    static Eigen eigen(double[][] input) {
        int n=input.length; double[][] a=copy(input), v=new double[n][n];
        for(int i=0;i<n;i++)v[i][i]=1;
        for(int iter=0;iter<200*n*n;iter++) {
            int p=0,q=Math.min(1,n-1); double off=0,scale=0;
            for(int i=0;i<n;i++) {
                scale=Math.max(scale,Math.abs(a[i][i]));
                for(int j=i+1;j<n;j++)if(Math.abs(a[i][j])>off){off=Math.abs(a[i][j]);p=i;q=j;}
            }
            if(off <= 2e-15*Math.max(scale,1e-300)) break;
            double tau=(a[q][q]-a[p][p])/(2*a[p][q]);
            double t=Math.copySign(1.0,tau)/(Math.abs(tau)+Math.hypot(1,tau));
            double c=1/Math.hypot(1,t),s=t*c,app=a[p][p],aqq=a[q][q],apq=a[p][q];
            a[p][p]=app-t*apq; a[q][q]=aqq+t*apq; a[p][q]=a[q][p]=0;
            for(int k=0;k<n;k++)if(k!=p&&k!=q){
                double kp=a[k][p],kq=a[k][q];
                a[k][p]=a[p][k]=c*kp-s*kq; a[k][q]=a[q][k]=s*kp+c*kq;
            }
            for(int k=0;k<n;k++){double kp=v[k][p],kq=v[k][q];v[k][p]=c*kp-s*kq;v[k][q]=s*kp+c*kq;}
        }
        Integer[] idx=new Integer[n];for(int i=0;i<n;i++)idx[i]=i;
        Arrays.sort(idx,(i,j)->Double.compare(a[j][j],a[i][i]));
        double[] w=new double[n];double[][] out=new double[n][n];
        for(int k=0;k<n;k++){
            w[k]=a[idx[k]][idx[k]];
            for(int i=0;i<n;i++)out[i][k]=v[i][idx[k]];
            int imax=0;for(int i=1;i<n;i++)if(Math.abs(out[i][k])>Math.abs(out[imax][k]))imax=i;
            if(out[imax][k]<0)for(int i=0;i<n;i++)out[i][k]*=-1;
        }
        return new Eigen(w,out);
    }

    static RankFit rank(double[][] a,int rank) {
        int n=a.length,p=a[0].length;double[][] gram=new double[p][p];
        for(double[] row:a)for(int j=0;j<p;j++)for(int k=0;k<p;k++)gram[j][k]+=row[j]*row[k];
        Eigen e=eigen(gram);double[] singular=new double[p];
        for(int j=0;j<p;j++)singular[j]=Math.sqrt(Math.max(0,e.values[j]));
        double[][] fitted=new double[n][p];double sse=0,max=0,total=0,kept=0;
        for(int j=0;j<p;j++){total+=Math.max(0,e.values[j]);if(j<rank)kept+=Math.max(0,e.values[j]);}
        for(int i=0;i<n;i++){
            for(int k=0;k<rank;k++){
                double score=0;for(int j=0;j<p;j++)score+=a[i][j]*e.vectors[j][k];
                for(int j=0;j<p;j++)fitted[i][j]+=score*e.vectors[j][k];
            }
            for(int j=0;j<p;j++){double r=fitted[i][j]-a[i][j];sse+=r*r;max=Math.max(max,Math.abs(r));}
        }
        return new RankFit(singular,e.vectors,fitted,sse,max,total==0?1:kept/total);
    }

    /** Stable intercept/slope fit using centered predictors. */
    static Regression regress(double[] x,double[] y) {
        if(x.length!=y.length||x.length<2)throw new IllegalArgumentException("Regression requires at least two pairs");
        double mx=Arrays.stream(x).average().orElseThrow(),my=Arrays.stream(y).average().orElseThrow();
        double xx=0,xy=0;for(int i=0;i<x.length;i++){xx+=(x[i]-mx)*(x[i]-mx);xy+=(x[i]-mx)*(y[i]-my);}
        if(!(xx>1e-280))throw new IllegalArgumentException("Degenerate predictor");
        double b=xy/xx,a=my-b*mx,sse=0;
        for(int i=0;i<x.length;i++){double r=a+b*x[i]-y[i];sse+=r*r;}
        return new Regression(a,b,sse,x.length);
    }

    /** Two predictors, no intercept; modified Gram-Schmidt QR after column scaling. */
    static LeastSquares leastSquares2(double[] x,double[] z,double[] y) {
        double nx=Math.sqrt(sumSquares(x)),nz=Math.sqrt(sumSquares(z));
        if(!(nx>0&&nz>0))throw new IllegalArgumentException("Zero predictor");
        double[] q1=new double[x.length],q2=new double[x.length];
        for(int i=0;i<x.length;i++)q1[i]=x[i]/nx;
        double r12=dot(q1,z)/nz;
        for(int i=0;i<x.length;i++)q2[i]=z[i]/nz-r12*q1[i];
        double r22=Math.sqrt(sumSquares(q2));
        if(r22<1e-13)throw new IllegalArgumentException("Collinear identity predictors");
        for(int i=0;i<x.length;i++)q2[i]/=r22;
        double b=dot(q2,y)/r22, a=dot(q1,y)-r12*b;
        double abs=Math.min(1,Math.abs(r12));
        return new LeastSquares(new double[]{a/nx,b/nz},Math.sqrt((1+abs)/(1-abs)));
    }

    /** Global coarse scan plus local golden search; includes both endpoints. */
    static Minimum gridGolden(DoubleUnaryOperator fn,double low,double high,int grid) {
        double best=low,bestValue=fn.applyAsDouble(low);int bestIndex=0;
        for(int i=1;i<=grid;i++){
            double x=low+(high-low)*i/grid,y=fn.applyAsDouble(x);
            if(y<bestValue){best=x;bestValue=y;bestIndex=i;}
        }
        double l=low+(high-low)*Math.max(0,bestIndex-1)/grid;
        double r=low+(high-low)*Math.min(grid,bestIndex+1)/grid;
        double golden=(Math.sqrt(5)-1)/2;
        double c=r-golden*(r-l),d=l+golden*(r-l),fc=fn.applyAsDouble(c),fd=fn.applyAsDouble(d);
        for(int it=0;it<120&&r-l>2e-13*(1+Math.abs(l)+Math.abs(r));it++){
            if(fc<fd){r=d;d=c;fd=fc;c=r-golden*(r-l);fc=fn.applyAsDouble(c);}
            else{l=c;c=d;fc=fd;d=l+golden*(r-l);fd=fn.applyAsDouble(d);}
        }
        double mid=(l+r)/2,val=fn.applyAsDouble(mid);
        if(val<bestValue){best=mid;bestValue=val;}
        return new Minimum(best,bestValue,Math.abs(best-low)<1e-6*(high-low)||Math.abs(best-high)<1e-6*(high-low));
    }
}
