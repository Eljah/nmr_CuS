package org.eljah.nmr.relaxation;

import java.util.*;
import java.util.function.DoubleUnaryOperator;

/** Deterministic low-dimensional variable projection and constrained least squares. */
final class Numerics {
    private Numerics() {}
    record Linear(double[] coefficients,double objective) {}
    record Minimum(double x,double value) {}

    /** Exact nonnegative two-column least squares, including both boundary faces. */
    static Linear nnls2(double[] y,double[] scale,double[] first,double[] second) {
        int n=y.length;
        double g00=0,g01=0,g11=0,b0=0,b1=0;
        for(int i=0;i<n;i++) {
            double u=first[i]/scale[i],v=second[i]/scale[i],z=y[i]/scale[i];
            g00+=u*u;g01+=u*v;g11+=v*v;b0+=u*z;b1+=v*z;
        }
        double c0=g00>0?Math.max(0,b0/g00):0,c1=0;
        double best=objective(y,scale,first,second,c0,c1);
        double test1=g11>0?Math.max(0,b1/g11):0;
        double test=objective(y,scale,first,second,0,test1);
        if(test<best){best=test;c0=0;c1=test1;}
        double determinant=g00*g11-g01*g01;
        if(determinant>1e-13*Math.max(g00*g11,1e-300)) {
            double u=(b0*g11-b1*g01)/determinant,v=(b1*g00-b0*g01)/determinant;
            if(u>=0 && v>=0) {
                test=objective(y,scale,first,second,u,v);
                if(test<best){best=test;c0=u;c1=v;}
            }
        }
        return new Linear(new double[]{c0,c1},best);
    }
    static double objective(double[] y,double[] s,double[] x0,double[] x1,double c0,double c1) {
        double sum=0;
        for(int j=0;j<y.length;j++){double r=(c0*x0[j]+c1*x1[j]-y[j])/s[j];sum+=r*r;}
        return sum;
    }

    /** Enumerate active sets. Used only for <=9 regression coefficients. */
    static Linear nnls(double[][] x,double[] y,double[] scale) {
        int n=y.length,p=x[0].length;
        if(p>12)throw new IllegalArgumentException("Active-set enumeration is limited to 12 columns.");
        double[][] gram=new double[p][p];double[] rhs=new double[p];
        for(int i=0;i<n;i++)for(int j=0;j<p;j++) {
            rhs[j]+=x[i][j]*y[i]/(scale[i]*scale[i]);
            for(int k=0;k<p;k++)gram[j][k]+=x[i][j]*x[i][k]/(scale[i]*scale[i]);
        }
        double[] best=new double[p];double bestObjective=objective(x,y,scale,best);
        for(int mask=1;mask<(1<<p);mask++) {
            int k=Integer.bitCount(mask);int[] selected=new int[k];
            for(int j=0,r=0;j<p;j++)if((mask&(1<<j))!=0)selected[r++]=j;
            double[][] a=new double[k][k];double[] b=new double[k];double[] norms=new double[k];
            boolean singular=false;
            for(int j=0;j<k;j++){norms[j]=Math.sqrt(gram[selected[j]][selected[j]]);if(norms[j]<1e-100)singular=true;}
            if(singular)continue;
            for(int j=0;j<k;j++){
                b[j]=rhs[selected[j]]/norms[j];
                for(int l=0;l<k;l++)a[j][l]=gram[selected[j]][selected[l]]/(norms[j]*norms[l]);
            }
            double[] solved=solve(a,b);
            if(solved==null)continue;
            double[] candidate=new double[p];boolean feasible=true;
            for(int j=0;j<k;j++) {
                candidate[selected[j]]=solved[j]/norms[j];
                if(candidate[selected[j]] < -1e-9)feasible=false;
                candidate[selected[j]]=Math.max(0,candidate[selected[j]]);
            }
            if(!feasible)continue;
            double objective=objective(x,y,scale,candidate);
            if(objective<bestObjective){bestObjective=objective;best=candidate;}
        }
        return new Linear(best,bestObjective);
    }
    static double objective(double[][] x,double[] y,double[] scale,double[] c) {
        double sum=0;
        for(int i=0;i<y.length;i++){double z=0;for(int j=0;j<c.length;j++)z+=c[j]*x[i][j];double r=(z-y[i])/scale[i];sum+=r*r;}
        return sum;
    }
    static double[] solve(double[][] a,double[] b) {
        int n=b.length;
        for(int j=0;j<n;j++) {
            int pivot=j;for(int k=j+1;k<n;k++)if(Math.abs(a[k][j])>Math.abs(a[pivot][j]))pivot=k;
            if(Math.abs(a[pivot][j])<1e-11)return null;
            double[] row=a[pivot];a[pivot]=a[j];a[j]=row;double val=b[pivot];b[pivot]=b[j];b[j]=val;
            for(int k=j+1;k<n;k++) {
                double ratio=a[k][j]/a[j][j];
                for(int l=j;l<n;l++)a[k][l]-=ratio*a[j][l];b[k]-=ratio*b[j];
            }
        }
        double[] result=new double[n];
        for(int j=n-1;j>=0;j--){double z=b[j];for(int k=j+1;k<n;k++)z-=a[j][k]*result[k];result[j]=z/a[j][j];}
        return result;
    }
    /** Scan the entire bounded domain, then refine every grid-local minimum. */
    static Minimum minimize(DoubleUnaryOperator f,double low,double high,double step) {
        int n=(int)Math.ceil((high-low)/step);
        double[] values=new double[n+1];
        Minimum best=new Minimum(low,Double.POSITIVE_INFINITY);
        for(int j=0;j<=n;j++) {
            double x=Math.min(high,low+j*step);values[j]=f.applyAsDouble(x);
            if(values[j]<best.value())best=new Minimum(x,values[j]);
        }
        for(int j=1;j<n;j++)if(values[j]<=values[j-1] && values[j]<=values[j+1]
                && (values[j]<values[j-1] || values[j]<values[j+1])) {
            Minimum m=golden(f,low+(j-1)*step,Math.min(high,low+(j+1)*step));
            if(m.value()<best.value())best=m;
        }
        return best;
    }
    static Minimum golden(DoubleUnaryOperator f,double a,double b) {
        double q=(Math.sqrt(5)-1)/2,c=b-q*(b-a),d=a+q*(b-a),fc=f.applyAsDouble(c),fd=f.applyAsDouble(d);
        for(int k=0;k<75 && b-a>1e-6;k++){
            if(fc<fd){b=d;d=c;fd=fc;c=b-q*(b-a);fc=f.applyAsDouble(c);}
            else{a=c;c=d;fc=fd;d=a+q*(b-a);fd=f.applyAsDouble(d);}
        }
        return fc<fd?new Minimum(c,fc):new Minimum(d,fd);
    }
    static double rms(double sum,int n){return n>0?Math.sqrt(sum/n):Double.NaN;}
}
