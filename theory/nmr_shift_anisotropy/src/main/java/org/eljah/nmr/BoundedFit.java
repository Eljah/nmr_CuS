package org.eljah.nmr;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/** Damped finite-difference Gauss-Newton; bounded LOCAL search, no global claim. */
final class BoundedFit {
    enum Kind { FREE_ETA_ISOTROPIC_K, AXIAL_ETA_AXIAL_K, FREE_ETA_AXIAL_K }
    static final double[] DIFFERENCE_STEPS = {.002,.002,.002,.002,.0015,.0015,.0015,.0015,.001,.002,.002,.00012,.00012,.00012,.00012};
    record Bounds(double[] low, double[] high, int[] active) { }
    record Iteration(int iteration, double sse, double lambda, double projectedGradient,
                     double largestScaledStep, boolean accepted, long evaluations) { }
    record Result(double[] p, double initialSse, double sse, int iterations, long evaluations,
                  boolean converged, String reason, List<Iteration> trace, double seconds) { }

    private BoundedFit() { }

    static Bounds bounds(Kind kind, double kappaLimit) {
        double[] lo = new double[15], hi = new double[15];
        Arrays.fill(hi,0,4,1);
        Arrays.fill(lo,4,8,74.8); Arrays.fill(hi,4,8,76.8);
        lo[8]=.025;hi[8]=.25;lo[9]=lo[10]=.05;hi[9]=hi[10]=.8;
        Arrays.fill(lo,11,15,-kappaLimit);Arrays.fill(hi,11,15,kappaLimit);
        if(kind==Kind.AXIAL_ETA_AXIAL_K)Arrays.fill(hi,0,4,0);
        if(kind==Kind.FREE_ETA_ISOTROPIC_K){Arrays.fill(lo,11,15,0);Arrays.fill(hi,11,15,0);}
        List<Integer> indices=new ArrayList<>();
        for(int i=0;i<15;i++)if(hi[i]>lo[i])indices.add(i);
        return new Bounds(lo,hi,indices.stream().mapToInt(Integer::intValue).toArray());
    }

    static double[] clamp(double[] input, Bounds b) {
        double[] result=input.clone();
        for(int i=0;i<result.length;i++)result[i]=Math.max(b.low()[i],Math.min(b.high()[i],result[i]));
        return result;
    }

    static Result fit(PowderModel model,double[] start,Bounds bounds,int maximumIterations) {
        long startTime=System.nanoTime(),startEvaluations=model.evaluations;
        double[] p=clamp(start,bounds);
        PowderModel.Evaluation current=model.evaluate(p);
        double initial=current.sse(),lambda=.01;
        List<Iteration> trace=new ArrayList<>();
        String reason="iteration_limit";boolean converged=false;
        int stagnant=0,iteration=0;
        for(;iteration<maximumIterations;iteration++) {
            double[][] jac=jacobian(model,p,bounds);
            int n=bounds.active().length;
            double[][] gram=new double[n][n];double[] gradient=new double[n];
            for(int i=0;i<n;i++)for(int row=0;row<current.residual().length;row++) {
                gradient[i]+=jac[i][row]*current.residual()[row];
                for(int j=0;j<n;j++)gram[i][j]+=jac[i][row]*jac[j][row];
            }
            double pg=0;
            boolean[] blocked=new boolean[n];
            for(int i=0;i<n;i++) {
                int id=bounds.active()[i];
                blocked[i]=(p[id]<=bounds.low()[id]+1e-12&&gradient[i]>0)
                    ||(p[id]>=bounds.high()[id]-1e-12&&gradient[i]<0);
                if(!blocked[i])pg=Math.max(pg,Math.abs(gradient[i])*(bounds.high()[id]-bounds.low()[id]));
            }
            if(pg<2e-6){reason="small_projected_gradient";converged=true;break;}
            boolean accepted=false;double scaledStep=0;
            double before=current.sse();
            for(int attempt=0;attempt<12;attempt++) {
                double[][] damped=new double[n][n];double[] rhs=new double[n];
                for(int i=0;i<n;i++) {
                    for(int j=0;j<n;j++)damped[i][j]=blocked[i]||blocked[j]?0:gram[i][j];
                    damped[i][i]+=blocked[i]?1:lambda*Math.max(gram[i][i],1e-8);
                    rhs[i]=blocked[i]?0:-gradient[i];
                }
                double[] step=PowderModel.solve(damped,rhs);
                if(step==null){lambda*=10;continue;}
                scaledStep=0;
                for(int i=0;i<n;i++) {
                    int id=bounds.active()[i];
                    scaledStep=Math.max(scaledStep,Math.abs(step[i])/(bounds.high()[id]-bounds.low()[id]));
                }
                double fraction=scaledStep>.15?.15/scaledStep:1;
                double[] trial=p.clone();scaledStep=0;
                for(int i=0;i<n;i++) {
                    int id=bounds.active()[i];
                    trial[id]=Math.max(bounds.low()[id],Math.min(bounds.high()[id],p[id]+fraction*step[i]));
                    scaledStep=Math.max(scaledStep,Math.abs(trial[id]-p[id])/(bounds.high()[id]-bounds.low()[id]));
                }
                PowderModel.Evaluation candidate=model.evaluate(trial);
                if(candidate.sse()<current.sse()-1e-12) {
                    p=trial;current=candidate;accepted=true;lambda=Math.max(1e-9,lambda/3);break;
                }
                lambda=Math.min(1e14,lambda*5);
            }
            trace.add(new Iteration(iteration,current.sse(),lambda,pg,scaledStep,accepted,model.evaluations-startEvaluations));
            if(!accepted){reason="no_decreasing_step_at_finite_difference_resolution";break;}
            double improvement=before-current.sse();
            if(improvement<2e-7*Math.max(1,before))stagnant++;else stagnant=0;
            if(stagnant>=3||scaledStep<1e-5){reason="small_local_improvements";converged=true;break;}
        }
        return new Result(p,initial,current.sse(),trace.size(),model.evaluations-startEvaluations,
            converged,reason,trace,(System.nanoTime()-startTime)/1e9);
    }

    /** Parameter-major residual derivatives, after refitting the linear nuisance parameters. */
    static double[][] jacobian(PowderModel model,double[] p,Bounds bounds) {
        double[][] jac=new double[bounds.active().length][];
        for(int col=0;col<bounds.active().length;col++) {
            int index=bounds.active()[col];double step=DIFFERENCE_STEPS[index];
            double[] lo=p.clone(),hi=p.clone();
            lo[index]=Math.max(bounds.low()[index],p[index]-step);
            hi[index]=Math.min(bounds.high()[index],p[index]+step);
            double[] rl=model.evaluate(lo).residual(),rh=model.evaluate(hi).residual();
            jac[col]=new double[rl.length];
            for(int row=0;row<rl.length;row++)jac[col][row]=(rh[row]-rl[row])/(hi[index]-lo[index]);
        }
        return jac;
    }
}
