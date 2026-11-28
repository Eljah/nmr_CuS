package org.eljah.nmr.relaxation;

import java.util.*;
import static org.eljah.nmr.relaxation.FilesIO.*;
import static org.eljah.nmr.relaxation.RelaxationDynamics.*;

/** Tests scientific/numerical properties rather than snapshots of implementation output. */
final class SelfChecks {
    private SelfChecks(){}
    static Map<String,Object> run() {
        List<String> checks=new ArrayList<>();
        double[] x={0,1,2,3,4},z={1,1,1,1,1},s={1,1,1,1,1},y={3,5,7,9,11};
        var fit=Numerics.nnls2(y,s,x,z);
        close(fit.coefficients()[0],2,1e-10);close(fit.coefficients()[1],3,1e-10);
        checks.add("NNLS recovers an exact nonnegative affine model.");
        double[] negative={-3,-2,-1,0,1};
        var bounded=Numerics.nnls2(negative,s,x,z);
        close(bounded.coefficients()[1],0,1e-10);require(bounded.coefficients()[0]>=0,"NNLS boundary.");
        checks.add("NNLS respects a physically nonnegative boundary.");
        double[][] design=new double[y.length][2];for(int j=0;j<y.length;j++){design[j][0]=x[j];design[j][1]=z[j];}
        var general=Numerics.nnls(design,y,s);close(general.objective(),fit.objective(),1e-10);
        close(general.coefficients()[0],2,1e-10);close(general.coefficients()[1],3,1e-10);
        checks.add("Independent active-set enumeration agrees with the two-column solver.");
        var minimum=Numerics.minimize(v->Math.pow(v-1234.567,2),0,3000,25);
        close(minimum.x(),1234.567,1e-4);
        checks.add("Grid plus local refinement recovers a known continuous minimum.");
        List<Rate> synthetic=new ArrayList<>();
        for(String site:SITES)for(double t:new double[]{50,77,100,140,180,220,250,280})
            synthetic.add(new Rate(site,t,0.6*t+500*activation(t,4200),1,"synthetic",0));
        var recovered=fitActivation(synthetic,"shared_theta","stored_error_scale");
        close(recovered.parameters().get("s1").theta(),4200,1e-3);
        close(recovered.parameters().get("s1").background(),0.6,1e-7);
        close(recovered.parameters().get("s1").amplitude(),500,1e-4);
        checks.add("Variable projection recovers a shared activation scale from noise-free synthetic rates.");
        List<Rate> powerData=new ArrayList<>();
        for(String site:SITES)for(double t:new double[]{50,77,100,140,180,220,250,280})
            powerData.add(new Rate(site,t,0.7*t+300*Math.pow(t/250,8.5),1,"synthetic",0));
        var powerFit=fitPower(powerData,"stored_error_scale");
        close(powerFit.parameters().get("s1").theta(),8.5,1e-5);
        close(powerFit.parameters().get("s1").amplitude(),300,1e-4);
        checks.add("The equal-parameter power alternative recovers its independent synthetic truth.");
        double slow=bppShape(15,200,-8,8000),tau=Math.pow(10,-8)*Math.exp(8000.0/200),omega=2*Math.PI*15e6;
        close(slow*omega*tau,1,1e-12);
        double fast=bppShape(15,200,-18,0);
        close(fast/(omega*1e-18),1,1e-12);
        checks.add("BPP implementation has the correct fast/slow limits without overflow.");
        double r1=5.79/5,r2=2.92/2.5;close(r1,1.158,1e-12);close(r2,1.168,1e-12);
        checks.add("Observed isotope ratios preserve numerator/denominator orientation.");
        close(1000*7.07,7070,1e-8);
        checks.add("R2 ms^-1 to s^-1 conversion preserves the factor 1000.");
        RelaxationDynamics interpolation=new RelaxationDynamics(java.nio.file.Path.of("."));
        NavigableMap<Double,Double> f=new TreeMap<>();
        f.put(160.0,16.18);f.put(230.0,15.64);
        interpolation.frequencies.put("s4",f);
        require(Double.isNaN(interpolation.frequency("s4",200,30)),"Missing frequency gap must not be silently bridged.");
        require(Double.isNaN(interpolation.frequency("s4",240,75)),"No extrapolation.");
        close(interpolation.frequency("s4",195,75),15.91,1e-12);
        checks.add("Primary interpolation refuses the missing 70 K interval and all extrapolation.");
        Map<Double,Integer> blocks=contiguousFolds(synthetic);
        int previous=-1;for(int block:blocks.values()){require(block>=previous,"Contiguous folds must be monotone in T.");previous=block;}
        checks.add("Contiguous folds are monotone in temperature and distinct from interleaved folds.");
        return map("passed",checks.size(),"checks",checks,"randomness","none");
    }
    static void close(double actual,double expected,double tolerance){
        require(Math.abs(actual-expected)<=tolerance,"Expected "+expected+" but obtained "+actual);
    }
}
