package org.eljah.nmr;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.Locale;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;

/** Reproducible conditional comparison of EFG asymmetry and axial magnetic shift. */
public final class NmrShiftAnisotropy {
    private NmrShiftAnisotropy() { }
    public static void main(String[] args) throws Exception {
        String mode = "all";
        for (String arg : args) {
            if (arg.startsWith("--mode=")) mode = arg.substring(7);
            else if (arg.equals("--help")) { System.out.println("Use --mode=all|baseline|fit|polish|finish|gridcheck. all runs baseline, local multistart, dense local polish and conditional diagnostics. finish/gridcheck need results/final_parameters.csv."); return; }
            else throw new IllegalArgumentException("Unknown option: " + arg);
        }
        if (!List.of("all","baseline","fit","polish","finish","gridcheck").contains(mode)) throw new IllegalArgumentException("Unknown mode: " + mode);
        PowderModel.Data[] data=baseline();
        if(mode.equals("fit")||mode.equals("all"))search(data);
        if(mode.equals("polish")||mode.equals("all"))NmrAnalysis.polish(data);
        if(mode.equals("finish")||mode.equals("all"))NmrAnalysis.finish(data);
        if(mode.equals("gridcheck"))NmrAnalysis.gridCheck(data);
    }

    record Candidate(BoundedFit.Kind kind,String label,double kappaLimit,PowderModel.Grid fitGrid,BoundedFit.Result result,
                     PowderModel.Grid evaluationGrid,double evaluationSse) { }

    private static void search(PowderModel.Data[] data) throws IOException {
        List<NmrFiles.Reference> refs=NmrFiles.references(Path.of("inputs/reference_solutions.csv"));
        List<Candidate> selected=new ArrayList<>();
        try(BufferedWriter trials=NmrFiles.writer(Path.of("results/search_trials.csv"));
            BufferedWriter iterations=NmrFiles.writer(Path.of("results/optimizer_iterations.csv"))) {
            List<Object> header=new ArrayList<>(List.of("model","start_label","stage","n_u","n_phi","kappa_limit","initial_SSE","final_SSE","iterations","evaluations","local_converged","termination","elapsed_seconds"));
            header.addAll(Arrays.asList(PowderModel.PARAMETER_NAMES));NmrFiles.csv(trials,header.toArray());
            NmrFiles.csv(iterations,"model","start_label","stage","iteration","SSE","damping","projected_gradient","largest_scaled_step","accepted","evaluations");
            for(BoundedFit.Kind kind:BoundedFit.Kind.values()) {
                System.out.println("SEARCH "+kind+" (local multistart; no global guarantee)");
                boolean axial=kind==BoundedFit.Kind.AXIAL_ETA_AXIAL_K;
                boolean isotropic=kind==BoundedFit.Kind.FREE_ETA_ISOTROPIC_K;
                PowderModel.Grid coarse=axial?new PowderModel.Grid(640,1):isotropic?new PowderModel.Grid(80,160):new PowderModel.Grid(40,80);
                PowderModel.Grid refine=axial?new PowderModel.Grid(1280,1):new PowderModel.Grid(80,160);
                PowderModel.Grid dense=axial?new PowderModel.Grid(2560,1):new PowderModel.Grid(160,320);
                List<double[]> starts=new ArrayList<>();List<String> labels=new ArrayList<>();
                if(axial){starts.add(refs.get(3).p().clone());labels.add("archived_axial_isotropic");}
                starts.add(refs.get(0).p().clone());labels.add("archived_free_seed29");
                starts.add(refs.get(1).p().clone());labels.add("archived_free_seed67");
                if(!isotropic) {
                    double[][] seeds={{-.004,.003,.005,-.006},{.006,-.004,-.003,.005}};
                    for(int k=0;k<2;k++){
                        double[] p=refs.get(k).p().clone();System.arraycopy(seeds[k],0,p,11,4);
                        starts.add(p);labels.add("deterministic_kappa_"+k);
                    }
                }
                PowderModel coarseModel=new PowderModel(data,coarse,false);
                List<Candidate> candidates=new ArrayList<>();
                for(int i=0;i<starts.size();i++) {
                    BoundedFit.Result fit=BoundedFit.fit(coarseModel,starts.get(i),BoundedFit.bounds(kind,.02),isotropic?38:axial?70:48);
                    logTrial(trials,iterations,kind,labels.get(i),"coarse",coarse,.02,fit);
                    candidates.add(new Candidate(kind,labels.get(i),.02,coarse,fit,coarse,fit.sse()));
                }
                PowderModel refinementModel=new PowderModel(data,refine,false);
                List<Candidate> ranked=new ArrayList<>();
                for(Candidate c:candidates)ranked.add(new Candidate(kind,c.label(),.02,coarse,c.result(),refine,refinementModel.evaluate(c.result().p()).sse()));
                ranked.sort(Comparator.comparingDouble(Candidate::evaluationSse));
                List<Candidate> refined=new ArrayList<>();
                for(int i=0;i<Math.min(2,ranked.size());i++) {
                    Candidate c=ranked.get(i);
                    BoundedFit.Result fit=BoundedFit.fit(refinementModel,c.result().p(),BoundedFit.bounds(kind,.02),isotropic?24:axial?55:45);
                    logTrial(trials,iterations,kind,c.label(),"refined",refine,.02,fit);
                    refined.add(new Candidate(kind,c.label(),.02,refine,fit,refine,fit.sse()));
                }
                PowderModel denseModel=new PowderModel(data,dense,false);
                Candidate best=null;
                for(Candidate c:refined) {
                    double e=denseModel.evaluate(c.result().p()).sse();
                    Candidate value=new Candidate(kind,c.label(),.02,refine,c.result(),dense,e);
                    if(best==null||e<best.evaluationSse())best=value;
                }
                selected.add(best);writeBest(selected);
                saveCandidateSpectra(data,denseModel,best,"best_"+kind.name().toLowerCase(Locale.ROOT));
                System.out.printf(Locale.ROOT,"SELECTED %s: dense E=%.12f, eta=%s, kappa=%s%n",kind,best.evaluationSse(),
                    Arrays.toString(Arrays.copyOfRange(best.result().p(),0,4)),Arrays.toString(Arrays.copyOfRange(best.result().p(),11,15)));
                boolean boundary=false;
                for(int i=11;i<15;i++)if(Math.abs(best.result().p()[i])>.0199)boundary=true;
                if(boundary&&!isotropic) {
                    BoundedFit.Result fit=BoundedFit.fit(refinementModel,best.result().p(),BoundedFit.bounds(kind,.04),55);
                    logTrial(trials,iterations,kind,best.label()+"_expanded_kappa","expanded_bound",refine,.04,fit);
                    Candidate expanded=new Candidate(kind,best.label()+"_expanded_kappa",.04,refine,fit,dense,denseModel.evaluate(fit.p()).sse());
                    selected.add(expanded);writeBest(selected);
                    saveCandidateSpectra(data,denseModel,expanded,"expanded_"+kind.name().toLowerCase(Locale.ROOT));
                }
            }
        }
        System.out.println("SEARCH FINISHED. Selected feasible local solutions saved; no global optimum or statistical confidence is claimed.");
    }

    private static void saveCandidateSpectra(PowderModel.Data[] data,PowderModel model,Candidate c,String prefix)throws IOException {
        PowderModel.Evaluation e=model.evaluate(c.result().p());
        for(int i=0;i<2;i++)NmrFiles.spectrum(Path.of("results/"+prefix+"_"+data[i].series+"ms.csv"),data[i],e.series()[i]);
    }

    private static void logTrial(BufferedWriter trials,BufferedWriter iterations,BoundedFit.Kind kind,String label,
                                 String stage,PowderModel.Grid grid,double bound,BoundedFit.Result fit)throws IOException {
        List<Object> row=new ArrayList<>(List.of(kind,label,stage,grid.nu(),grid.nphi(),bound,fit.initialSse(),fit.sse(),
            fit.iterations(),fit.evaluations(),fit.converged(),fit.reason(),fit.seconds()));
        for(double p:fit.p())row.add(p);NmrFiles.csv(trials,row.toArray());trials.flush();
        for(BoundedFit.Iteration i:fit.trace())NmrFiles.csv(iterations,kind,label,stage,i.iteration(),i.sse(),i.lambda(),i.projectedGradient(),i.largestScaledStep(),i.accepted(),i.evaluations());
        iterations.flush();
        System.out.printf(Locale.ROOT,"  %s %s: %.9f -> %.9f, %.1f s, %s%n",stage,label,fit.initialSse(),fit.sse(),fit.seconds(),fit.reason());
    }

    private static void writeBest(List<Candidate> candidates)throws IOException {
        try(BufferedWriter w=NmrFiles.writer(Path.of("results/best_parameters.csv"))) {
            List<Object> header=new ArrayList<>(List.of("model","start_label","fit_n_u","fit_n_phi","kappa_limit","fit_SSE",
                "evaluation_n_u","evaluation_n_phi","evaluation_SSE","local_converged","termination"));
            header.addAll(Arrays.asList(PowderModel.PARAMETER_NAMES));NmrFiles.csv(w,header.toArray());
            for(Candidate c:candidates) {
                List<Object> row=new ArrayList<>(List.of(c.kind(),c.label(),c.fitGrid().nu(),c.fitGrid().nphi(),c.kappaLimit(),c.result().sse(),
                    c.evaluationGrid().nu(),c.evaluationGrid().nphi(),c.evaluationSse(),c.result().converged(),c.result().reason()));
                for(double v:c.result().p())row.add(v);NmrFiles.csv(w,row.toArray());
            }
        }
    }

    static PowderModel.Data[] baseline() throws IOException {
        Files.createDirectories(Path.of("results"));
        PowderModel.Data[] data = PowderModel.readData(Path.of("inputs/spectrum_integrals.csv"));
        if (data[0].n != 245 || data[1].n != 245) throw new IllegalStateException("Archive mask/count mismatch.");
        double[] spinErrors = checkSpinReference();
        System.out.printf(Locale.ROOT,"Independent 4x4 Hamiltonian check: max transition error %.3g MHz, RF error %.3g, tensor-RF error %.3g.%n",
            spinErrors[0],spinErrors[1],spinErrors[2]);
        List<NmrFiles.Reference> references = NmrFiles.references(Path.of("inputs/reference_solutions.csv"));
        try (BufferedWriter w = NmrFiles.writer(Path.of("results/baseline_reproduction.csv"))) {
            NmrFiles.csv(w,"label","n_u","n_phi","archived_SSE","java_SSE","absolute_error","SSE_500ms","SSE_20ms","elapsed_seconds","status");
            for (NmrFiles.Reference ref : references) {
                long start = System.nanoTime();
                PowderModel model = new PowderModel(data, ref.grid(), false);
                PowderModel.Evaluation result = model.evaluate(ref.p());
                double error = Math.abs(result.sse() - ref.expected());
                boolean pass = error < 2e-8;
                NmrFiles.csv(w,ref.label(),ref.grid().nu(),ref.grid().nphi(),ref.expected(),result.sse(),error,
                    result.series()[0].sse(),result.series()[1].sse(),(System.nanoTime()-start)/1e9,pass?"PASS":"FAIL");
                w.flush();
                System.out.printf(Locale.ROOT,"Baseline %s %dx%d: E=%.15f, archive error=%.3g, %.2f s, %s%n",
                    ref.label(),ref.grid().nu(),ref.grid().nphi(),result.sse(),error,(System.nanoTime()-start)/1e9,pass?"PASS":"FAIL");
                for (int i = 0; i < 2; i++) NmrFiles.spectrum(Path.of("results/baseline_"+ref.label()+"_"+data[i].series+"ms.csv"),data[i],result.series()[i]);
                if (!pass) throw new IllegalStateException("Baseline mismatch: " + ref.label());
            }
        }
        return data;
    }

    private static double[] checkSpinReference() throws IOException {
        List<String> rows = Files.readAllLines(Path.of("inputs/hamiltonian_reference_numpy.csv"));
        SpinHamiltonian spin = new SpinHamiltonian();
        double[] maximum = new double[4];
        try (BufferedWriter w = NmrFiles.writer(Path.of("results/hamiltonian_checks.csv"))) {
            NmrFiles.csv(w,"case","level_a","level_b","transition_error_MHz","archive_RF_error","tensor_RF_error","orthogonality_error");
            for (int row = 1; row < rows.size(); row++) {
                String[] f = rows.get(row).split(",");
                double eta=Double.parseDouble(f[1]),nqr=Double.parseDouble(f[2]),centre=Double.parseDouble(f[3]),kappa=Double.parseDouble(f[4]);
                double nx=Double.parseDouble(f[5]),ny=Double.parseDouble(f[6]),nz=Double.parseDouble(f[7]);
                int a=Integer.parseInt(f[8]),b=Integer.parseInt(f[9]);
                double gp=1-kappa/3,gz=1+2*kappa/3;
                spin.solve(centre*gp*nx,centre*gp*ny,centre*gz*nz,nqr/Math.sqrt(1+eta*eta/3),eta);
                double[] errors={Math.abs(spin.eigenvalues[b]-spin.eigenvalues[a]-Double.parseDouble(f[10])),
                    Math.abs(spin.weight(a,b,nx,ny,nz,1,1)-Double.parseDouble(f[11])),
                    Math.abs(spin.weight(a,b,nx,ny,nz,gp,gz)-Double.parseDouble(f[12])),spin.orthogonalityError()};
                for(int i=0;i<4;i++)maximum[i]=Math.max(maximum[i],errors[i]);
                NmrFiles.csv(w,f[0],a,b,errors[0],errors[1],errors[2],errors[3]);
            }
        }
        if(maximum[0]>1e-9||maximum[1]>1e-10||maximum[2]>1e-10||maximum[3]>1e-11)
            throw new IllegalStateException("Independent Hermitian eigensolver/RF check failed.");
        return maximum;
    }
}
