package org.eljah.nmr;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;

/** Fixed-parameter numerical checks and explicitly conditional diagnostics. */
final class NmrAnalysis {
    private NmrAnalysis() { }
    record Solution(String name,BoundedFit.Kind kind,String label,double limit,double[] p,
                    int fitNu,int fitPhi,boolean converged,String termination) { }

    static List<Solution> readSolutions(Path path) throws IOException {
        List<String> lines=Files.readAllLines(path);Map<String,Integer> index=new HashMap<>();
        String[] headings=lines.get(0).split(",");for(int i=0;i<headings.length;i++)index.put(headings[i],i);
        List<Solution> result=new ArrayList<>();
        for(int row=1;row<lines.size();row++) {
            if(lines.get(row).isBlank())continue;
            String[] values=lines.get(row).split(",",-1);double[] p=new double[15];
            for(int i=0;i<p.length;i++)p[i]=Double.parseDouble(values[index.get(PowderModel.PARAMETER_NAMES[i])]);
            BoundedFit.Kind kind=BoundedFit.Kind.valueOf(values[index.get("model")]);
            double limit=Double.parseDouble(values[index.get("kappa_limit")]);
            String name=kind.name().toLowerCase(Locale.ROOT)+(limit>.020001?"_expanded":"");
            result.add(new Solution(name,kind,values[index.get("start_label")],limit,p,
                Integer.parseInt(values[index.get("fit_n_u")]),Integer.parseInt(values[index.get("fit_n_phi")]),
                Boolean.parseBoolean(values[index.get("local_converged")]),values[index.get("termination")]));
        }
        return result;
    }

    static void polish(PowderModel.Data[] data) throws IOException {
        List<Solution> selected=readSolutions(Path.of("results/best_parameters.csv"));
        try(BufferedWriter w=NmrFiles.writer(Path.of("results/final_parameters.csv"));
            BufferedWriter trace=NmrFiles.writer(Path.of("results/dense_polish_iterations.csv"));
            BufferedWriter summary=NmrFiles.writer(Path.of("results/dense_polish_summary.csv"))) {
            List<Object> h=new ArrayList<>(List.of("model","start_label","fit_n_u","fit_n_phi","kappa_limit","fit_SSE","local_converged","termination"));
            h.addAll(Arrays.asList(PowderModel.PARAMETER_NAMES));NmrFiles.csv(w,h.toArray());
            NmrFiles.csv(trace,"model","iteration","SSE","damping","projected_gradient","largest_scaled_step","accepted","evaluations");
            NmrFiles.csv(summary,"model","n_u","n_phi","initial_SSE","final_SSE","iterations","evaluations","local_converged","termination","elapsed_seconds");
            for(Solution s:selected) {
                boolean axial=s.kind()==BoundedFit.Kind.AXIAL_ETA_AXIAL_K;
                PowderModel.Grid grid=axial?new PowderModel.Grid(2560,1):new PowderModel.Grid(160,320);
                PowderModel model=new PowderModel(data,grid,false);
                BoundedFit.Result fit=BoundedFit.fit(model,s.p(),BoundedFit.bounds(s.kind(),s.limit()),s.kind()==BoundedFit.Kind.FREE_ETA_AXIAL_K?32:18);
                List<Object> row=new ArrayList<>(List.of(s.kind(),s.label()+"_dense_polish",grid.nu(),grid.nphi(),s.limit(),fit.sse(),fit.converged(),fit.reason()));
                for(double value:fit.p())row.add(value);NmrFiles.csv(w,row.toArray());w.flush();
                for(BoundedFit.Iteration it:fit.trace())NmrFiles.csv(trace,s.name(),it.iteration(),it.sse(),it.lambda(),it.projectedGradient(),it.largestScaledStep(),it.accepted(),it.evaluations());trace.flush();
                NmrFiles.csv(summary,s.name(),grid.nu(),grid.nphi(),fit.initialSse(),fit.sse(),fit.iterations(),fit.evaluations(),fit.converged(),fit.reason(),fit.seconds());summary.flush();
                System.out.printf(Locale.ROOT,"DENSE POLISH %s: %.12f -> %.12f (%s; %.1f s)%n",s.name(),fit.initialSse(),fit.sse(),fit.reason(),fit.seconds());
            }
        }
    }

    static void finish(PowderModel.Data[] data) throws IOException {
        List<Solution> solutions=readSolutions(Path.of("results/final_parameters.csv"));
        try(BufferedWriter comparison=NmrFiles.writer(Path.of("results/final_model_comparison.csv"));
            BufferedWriter convergence=NmrFiles.writer(Path.of("results/angular_convergence.csv"));
            BufferedWriter rf=NmrFiles.writer(Path.of("results/rf_tensor_sensitivity.csv"));
            BufferedWriter parameters=NmrFiles.writer(Path.of("results/site_parameters.csv"));
            BufferedWriter nuisance=NmrFiles.writer(Path.of("results/linear_nuisance.csv"));
            BufferedWriter boundaries=NmrFiles.writer(Path.of("results/parameter_bounds.csv"));
            BufferedWriter scan=NmrFiles.writer(Path.of("results/conditional_kappa_scan.csv"))) {
            NmrFiles.csv(comparison,"model","n_u","n_phi","joint_SSE","SSE_500ms","SSE_20ms","joint_RMS_normalized","n_points","n_nonlinear","n_linear","sigma_MHz","rho500","rho20","maximum_abs_kappa","active_parameter_bounds","local_converged","termination");
            NmrFiles.csv(convergence,"model","n_u","n_phi","joint_SSE","SSE_500ms","SSE_20ms","change_from_previous_SSE","RMS_prediction_change_normalized","max_prediction_change_normalized","nonlinear_parameters_fixed");
            NmrFiles.csv(rf,"model","n_u","n_phi","archive_RF_SSE","tensor_RF_SSE","SSE_change","RMS_prediction_change_normalized","max_prediction_change_normalized","linear_nuisance_refitted","nonlinear_parameters_fixed");
            NmrFiles.csv(parameters,"model","site","NQR63_MHz","NQR65_MHz","eta","nuQ63_MHz","nuQ65_MHz","nu_iso63_MHz","kappa","nu_perp63_MHz","nu_parallel63_MHz","nu_parallel_minus_perp63_MHz","field_reference_T");
            NmrFiles.csv(nuisance,"model","series_ms","normalization","n_used","SSE","amplitude_s1_normalized","amplitude_s2_normalized","amplitude_s3_normalized","amplitude_s4_normalized","baseline_b0_normalized","baseline_b1_normalized","minimum_predicted_normalized","negative_predictions_used");
            NmrFiles.csv(boundaries,"model","parameter","value","lower","upper","active_in_search","fraction_distance_to_nearest_bound","boundary_at_tolerance_1e-4_range");
            NmrFiles.csv(scan,"model","clipped_absolute_kappa_limit","joint_SSE","SSE_minus_final","nonlinear_others_fixed","linear_nuisance_refitted","interpretation");
            for(Solution s:solutions) {
                boolean axial=s.kind()==BoundedFit.Kind.AXIAL_ETA_AXIAL_K;
                int[][] grids=axial?new int[][]{{640,1},{1280,1},{2560,1},{5120,1}}:new int[][]{{80,160},{160,320},{320,640}};
                PowderModel.Evaluation previous=null,finalEvaluation=null;PowderModel finalModel=null;
                for(int[] g:grids) {
                    PowderModel model=new PowderModel(data,new PowderModel.Grid(g[0],g[1]),false);
                    PowderModel.Evaluation e=model.evaluate(s.p());double[] changes=previous==null?new double[]{Double.NaN,Double.NaN}:predictionChanges(previous,e,data);
                    NmrFiles.csv(convergence,s.name(),g[0],g[1],e.sse(),e.series()[0].sse(),e.series()[1].sse(),previous==null?Double.NaN:e.sse()-previous.sse(),changes[0],changes[1],true);
                    convergence.flush();previous=e;finalEvaluation=e;finalModel=model;
                }
                PowderModel.Evaluation e=finalEvaluation;PowderModel.Grid grid=finalModel.grid;
                BoundedFit.Bounds b=BoundedFit.bounds(s.kind(),s.limit());int nBoundary=0;
                for(int i=0;i<s.p().length;i++) {
                    boolean active=b.high()[i]>b.low()[i];double distance=active?Math.min(s.p()[i]-b.low()[i],b.high()[i]-s.p()[i])/(b.high()[i]-b.low()[i]):Double.NaN;
                    boolean boundary=active&&distance<1e-4;if(boundary)nBoundary++;
                    NmrFiles.csv(boundaries,s.name(),PowderModel.PARAMETER_NAMES[i],s.p()[i],b.low()[i],b.high()[i],active,distance,boundary);
                }
                double maximumK=0;for(int i=11;i<15;i++)maximumK=Math.max(maximumK,Math.abs(s.p()[i]));
                NmrFiles.csv(comparison,s.name(),grid.nu(),grid.nphi(),e.sse(),e.series()[0].sse(),e.series()[1].sse(),Math.sqrt(e.sse()/(data[0].n+data[1].n)),data[0].n+data[1].n,b.active().length,12,s.p()[8],s.p()[9],s.p()[10],maximumK,nBoundary,s.converged(),s.termination());comparison.flush();
                for(int series=0;series<2;series++) {
                    PowderModel.SeriesResult r=e.series()[series];NmrFiles.spectrum(Path.of("results/final_"+s.name()+"_"+data[series].series+"ms.csv"),data[series],r);
                    double minimum=Double.POSITIVE_INFINITY;int negative=0;
                    for(int i=0;i<data[series].x.length;i++)if(data[series].used[i]){minimum=Math.min(minimum,r.predictedNormalized()[i]);if(r.predictedNormalized()[i]<0)negative++;}
                    NmrFiles.csv(nuisance,s.name(),data[series].series,data[series].scale,data[series].n,r.sse(),r.amplitudesNormalized()[0],r.amplitudesNormalized()[1],r.amplitudesNormalized()[2],r.amplitudesNormalized()[3],r.baselineNormalized()[0],r.baselineNormalized()[1],minimum,negative);
                }
                for(int site=0;site<4;site++) {
                    double eta=s.p()[site],k=s.p()[11+site],centre=s.p()[4+site],factor=Math.sqrt(1+eta*eta/3);
                    NmrFiles.csv(parameters,s.name(),"s"+(site+1),PowderModel.NQR[site][0],PowderModel.NQR[site][1],eta,PowderModel.NQR[site][0]/factor,PowderModel.NQR[site][1]/factor,centre,k,centre*(1-k/3),centre*(1+2*k/3),centre*k,PowderModel.FIELD_REFERENCE_T);
                }
                PowderModel tensorModel=new PowderModel(data,grid,true);PowderModel.Evaluation tensor=tensorModel.evaluate(s.p());double[] change=predictionChanges(e,tensor,data);
                NmrFiles.csv(rf,s.name(),grid.nu(),grid.nphi(),e.sse(),tensor.sse(),tensor.sse()-e.sse(),change[0],change[1],true,true);rf.flush();
                for(int series=0;series<2;series++)NmrFiles.spectrum(Path.of("results/rf_tensor_"+s.name()+"_"+data[series].series+"ms.csv"),data[series],tensor.series()[series]);
                if(s.kind()!=BoundedFit.Kind.FREE_ETA_ISOTROPIC_K)for(double limit:new double[]{0,.001,.002,.003,.004,.005,.0075,.01,.02,.04}) {
                    if(limit>s.limit()+1e-9)continue;
                    double[] p=s.p().clone();for(int i=11;i<15;i++)p[i]=Math.max(-limit,Math.min(limit,p[i]));
                    double conditional=finalModel.evaluate(p).sse();
                    NmrFiles.csv(scan,s.name(),limit,conditional,conditional-e.sse(),true,true,"feasible_curve_not_reoptimized_profile");scan.flush();
                }
                identify(data,s);
                forecasts(data,s,e);
                System.out.printf(Locale.ROOT,"FINAL %s: E=%.12f (%dx%d); tensor-RF change=%+.8f; active bounds=%d%n",s.name(),e.sse(),grid.nu(),grid.nphi(),tensor.sse()-e.sse(),nBoundary);
            }
        }
        writeMetadata(data,solutions);
        gridCheck(data);
    }

    static void gridCheck(PowderModel.Data[] data)throws IOException {
        try(BufferedWriter w=NmrFiles.writer(Path.of("results/ultradense_grid_check.csv"))) {
            NmrFiles.csv(w,"model","previous_n_u","previous_n_phi","check_n_u","check_n_phi","previous_SSE","check_SSE","SSE_change","RMS_prediction_change_normalized","max_prediction_change_normalized","nonlinear_parameters_fixed","linear_nuisance_refitted");
            for(Solution s:readSolutions(Path.of("results/final_parameters.csv"))) {
                boolean axial=s.kind()==BoundedFit.Kind.AXIAL_ETA_AXIAL_K;
                PowderModel.Grid before=axial?new PowderModel.Grid(5120,1):new PowderModel.Grid(320,640);
                PowderModel.Grid after=axial?new PowderModel.Grid(10240,1):new PowderModel.Grid(640,1280);
                PowderModel.Evaluation first=new PowderModel(data,before,false).evaluate(s.p());
                PowderModel.Evaluation second=new PowderModel(data,after,false).evaluate(s.p());
                double[] change=predictionChanges(first,second,data);
                NmrFiles.csv(w,s.name(),before.nu(),before.nphi(),after.nu(),after.nphi(),first.sse(),second.sse(),second.sse()-first.sse(),change[0],change[1],true,true);w.flush();
                for(int series=0;series<2;series++)NmrFiles.spectrum(Path.of("results/ultradense_"+s.name()+"_"+data[series].series+"ms.csv"),data[series],second.series()[series]);
                System.out.printf(Locale.ROOT,"ULTRADENSE CHECK %s: %.12f -> %.12f; prediction RMS change=%.6g%n",s.name(),first.sse(),second.sse(),change[0]);
            }
        }
    }

    private static double[] predictionChanges(PowderModel.Evaluation first,PowderModel.Evaluation second,PowderModel.Data[] data) {
        double sum=0,max=0;int n=0;
        for(int series=0;series<2;series++)for(int row=0;row<data[series].x.length;row++)if(data[series].used[row]) {
            double difference=second.series()[series].predictedNormalized()[row]-first.series()[series].predictedNormalized()[row];
            sum+=difference*difference;max=Math.max(max,Math.abs(difference));n++;
        }
        return new double[]{Math.sqrt(sum/n),max};
    }

    private static void identify(PowderModel.Data[] data,Solution s) throws IOException {
        // A local derivative diagnostic, not a likelihood Hessian or confidence interval.
        PowderModel.Grid grid=s.kind()==BoundedFit.Kind.AXIAL_ETA_AXIAL_K?new PowderModel.Grid(2560,1):new PowderModel.Grid(160,320);
        PowderModel model=new PowderModel(data,grid,false);BoundedFit.Bounds bounds=BoundedFit.bounds(s.kind(),s.limit());
        double[][] jac=BoundedFit.jacobian(model,s.p(),bounds);int n=jac.length;double[][] gram=new double[n][n];
        for(int i=0;i<n;i++)for(int j=0;j<n;j++)for(int row=0;row<jac[i].length;row++)gram[i][j]+=jac[i][row]*jac[j][row];
        double[][] scaled=new double[n][n];double[] scales=new double[n];
        for(int i=0;i<n;i++)scales[i]=bounds.high()[bounds.active()[i]]-bounds.low()[bounds.active()[i]];
        for(int i=0;i<n;i++)for(int j=0;j<n;j++)scaled[i][j]=gram[i][j]*scales[i]*scales[j];
        RealEigen eig=RealEigen.of(scaled);double[][] inverse=new double[n][n];double largest=Math.max(0,eig.values[n-1]);
        for(int k=0;k<n;k++)if(eig.values[k]>largest*1e-12)for(int i=0;i<n;i++)for(int j=0;j<n;j++)inverse[i][j]+=eig.vectors[i][k]*eig.vectors[j][k]/eig.values[k];
        double eigenResidual=0,orthogonalError=0,gramNorm=0;
        for(int i=0;i<n;i++)for(int j=0;j<n;j++)gramNorm=Math.max(gramNorm,Math.abs(scaled[i][j]));
        for(int k=0;k<n;k++)for(int i=0;i<n;i++) {
            double residual=-eig.values[k]*eig.vectors[i][k];for(int j=0;j<n;j++)residual+=scaled[i][j]*eig.vectors[j][k];
            eigenResidual=Math.max(eigenResidual,Math.abs(residual)/Math.max(1,gramNorm));
            double dot=0;for(int j=0;j<n;j++)dot+=eig.vectors[j][i]*eig.vectors[j][k];
            orthogonalError=Math.max(orthogonalError,Math.abs(dot-(i==k?1:0)));
        }
        if(eigenResidual>1e-9||orthogonalError>1e-10)throw new IllegalStateException("Local Gram eigensystem failed consistency check.");
        try(BufferedWriter modes=NmrFiles.writer(Path.of("results/local_identifiability_"+s.name()+".csv"));
            BufferedWriter correlation=NmrFiles.writer(Path.of("results/local_parameter_correlation_"+s.name()+".csv"));
            BufferedWriter derivatives=NmrFiles.writer(Path.of("results/local_derivative_norms_"+s.name()+".csv"))) {
            List<Object> header=new ArrayList<>(List.of("mode_weak_to_strong","eigenvalue_scaled_JTJ","singular_value_scaled_J","relative_singular_value","n_u","n_phi","pseudoinverse_included","eigensystem_relative_residual","eigenvectors_orthogonality_error"));
            for(int i:bounds.active())header.add(PowderModel.PARAMETER_NAMES[i]+"_scaled_direction");NmrFiles.csv(modes,header.toArray());
            for(int k=0;k<n;k++) {
                List<Object> row=new ArrayList<>(List.of(k+1,eig.values[k],Math.sqrt(Math.max(0,eig.values[k])),Math.sqrt(Math.max(0,eig.values[k])/largest),grid.nu(),grid.nphi(),eig.values[k]>largest*1e-12,eigenResidual,orthogonalError));
                for(int i=0;i<n;i++)row.add(eig.vectors[i][k]);NmrFiles.csv(modes,row.toArray());
            }
            NmrFiles.csv(correlation,"parameter_a","parameter_b","local_inverse_gram_correlation","not_a_statistical_confidence_result");
            NmrFiles.csv(derivatives,"parameter","finite_difference_step","search_range_scale","residual_derivative_norm","scaled_derivative_norm","n_u","n_phi");
            for(int i=0;i<n;i++) {
                int id=bounds.active()[i];NmrFiles.csv(derivatives,PowderModel.PARAMETER_NAMES[id],BoundedFit.DIFFERENCE_STEPS[id],scales[i],Math.sqrt(gram[i][i]),Math.sqrt(scaled[i][i]),grid.nu(),grid.nphi());
                for(int j=i+1;j<n;j++)NmrFiles.csv(correlation,PowderModel.PARAMETER_NAMES[id],PowderModel.PARAMETER_NAMES[bounds.active()[j]],inverse[i][i]*inverse[j][j]>0?inverse[i][j]/Math.sqrt(inverse[i][i]*inverse[j][j]):Double.NaN,true);
            }
        }
    }

    private static void forecasts(PowderModel.Data[] data,Solution s,PowderModel.Evaluation reference) throws IOException {
        PowderModel.Grid grid=s.kind()==BoundedFit.Kind.AXIAL_ETA_AXIAL_K?new PowderModel.Grid(5120,1):new PowderModel.Grid(320,640);
        PowderModel model=new PowderModel(data,grid,false);
        try(BufferedWriter w=NmrFiles.writer(Path.of("results/field_forecasts_"+s.name()+".csv"));
            BufferedWriter moments=NmrFiles.writer(Path.of("results/field_forecast_summary_"+s.name()+".csv"))) {
            NmrFiles.csv(w,"field_T","frequency_MHz","s1_Cu63_kernel","s1_Cu65_kernel","s2_Cu63_kernel","s2_Cu65_kernel","s3_Cu63_kernel","s3_Cu65_kernel","s4_Cu63_kernel","s4_Cu65_kernel","total500_unit_area","total20_unit_area","baseline_excluded");
            NmrFiles.csv(moments,"field_T","series_ms","area_before_normalizing","mean_MHz","standard_deviation_MHz","assumption");
            for(double field:new double[]{4.5,PowderModel.FIELD_REFERENCE_T,9}) {
                double ratio=field/PowderModel.FIELD_REFERENCE_T;double[][][] profiles=new double[4][][];
                double min=Double.POSITIVE_INFINITY,max=Double.NEGATIVE_INFINITY;
                for(int site=0;site<4;site++) {
                    profiles[site]=model.profiles(site,s.p()[site],s.p()[4+site],s.p()[11+site],s.p()[8],ratio);
                    min=Math.min(min,s.p()[4+site]*ratio-18);max=Math.max(max,s.p()[4+site]*ratio*PowderModel.GAMMA_RATIO+18);
                }
                double step=.02;int count=(int)Math.ceil((max-min)/step)+1;double[][] values=new double[count][10];double[] area=new double[2];
                for(int i=0;i<count;i++) {
                    double frequency=min+i*step;
                    for(int site=0;site<4;site++)for(int isotope=0;isotope<2;isotope++)values[i][2*site+isotope]=PowderModel.interpolate(profiles[site][isotope],frequency-s.p()[4+site]*ratio*(isotope==0?1:PowderModel.GAMMA_RATIO));
                    for(int series=0;series<2;series++)for(int site=0;site<4;site++)values[i][8+series]+=reference.series()[series].amplitudesNormalized()[site]*(values[i][2*site]+s.p()[9+series]*values[i][2*site+1]);
                    for(int series=0;series<2;series++)area[series]+=values[i][8+series]*step;
                }
                double[] mean=new double[2],second=new double[2];
                for(int i=0;i<count;i++) {
                    double frequency=min+i*step;List<Object> row=new ArrayList<>(List.of(field,frequency));
                    for(int k=0;k<8;k++)row.add(values[i][k]);
                    for(int series=0;series<2;series++) {
                        double density=values[i][8+series]/area[series];row.add(density);mean[series]+=frequency*density*step;second[series]+=frequency*frequency*density*step;
                    }
                    row.add(true);NmrFiles.csv(w,row.toArray());
                }
                for(int series=0;series<2;series++)NmrFiles.csv(moments,field,data[series].series,area[series],mean[series],Math.sqrt(Math.max(0,second[series]-mean[series]*mean[series])),"fixed_NQR_eta_kappa_relative_amplitudes_rho_and_sigma_MHz_centres_scale_with_B_no_background");
            }
        }
    }

    private static void writeMetadata(PowderModel.Data[] data,List<Solution> solutions)throws IOException {
        StringBuilder out=new StringBuilder("{\n  \"model\": \"exact I=3/2 Hermitian Hamiltonian; axial effective Zeeman anisotropy coaxial with EFG\",\n");
        out.append("  \"field_reference_T\": ").append(NmrFiles.number(PowderModel.FIELD_REFERENCE_T)).append(",\n");
        out.append("  \"effective_kappa_definition\": \"(nu_parallel - nu_perp)/nu_iso; absolute Knight shift is not calibrated\",\n");
        out.append("  \"objective\": \"unweighted SSE after separate per-series maximum normalization; descriptive, not chi-square\",\n");
        out.append("  \"n_used\": [").append(data[0].n).append(", ").append(data[1].n).append("],\n");
        out.append("  \"normalization_500_20\": ").append(NmrFiles.array(new double[]{data[0].scale,data[1].scale})).append(",\n");
        out.append("  \"primary_RF\": \"archive isotropic RF operator; anisotropic static Zeeman Hamiltonian\",\n");
        out.append("  \"RF_sensitivity\": \"full axial RF tensor at fixed nonlinear parameters; linear nuisance reoptimized\",\n");
        out.append("  \"optimizer\": \"bounded finite-difference damped Gauss-Newton; deterministic local multistart and one dense warm polish; no global guarantee\",\n");
        out.append("  \"kappa_scan\": \"clipped fixed-nonlinear feasible curve with linear nuisance refitted; not a profile minimum or confidence bound\",\n");
        out.append("  \"identifiability\": \"local variable-projection finite-difference Jacobian; inverse Gram correlations are geometric, not statistical confidence\",\n");
        out.append("  \"forecasts\": \"fixed spectroscopic parameters at other fields; no absolute Knight calibration, no extrapolated linear background, no relaxation prediction\",\n");
        out.append("  \"models\": [\n");
        for(int i=0;i<solutions.size();i++){Solution s=solutions.get(i);if(i>0)out.append(",\n");out.append("    {\"name\": ").append(NmrFiles.quote(s.name())).append(", \"kappa_search_limit\": ").append(NmrFiles.number(s.limit())).append(", \"local_termination\": ").append(NmrFiles.quote(s.termination())).append(", \"p\": ").append(NmrFiles.array(s.p())).append("}");}
        out.append("\n  ]\n}\n");Files.writeString(Path.of("results/analysis_metadata.json"),out);
    }

    /** Jacobi eigensystem for a small real symmetric Gram matrix; eigenvectors in columns. */
    static final class RealEigen {
        final double[] values;final double[][] vectors;
        private RealEigen(double[] values,double[][] vectors){this.values=values;this.vectors=vectors;}
        static RealEigen of(double[][] matrix) {
            int n=matrix.length;double[][] a=new double[n][n],v=new double[n][n];
            for(int i=0;i<n;i++){a[i]=matrix[i].clone();v[i][i]=1;}
            for(int iteration=0;iteration<100*n*n;iteration++) {
                int p=0,q=1;double maximum=0,diagonal=0;
                for(int i=0;i<n;i++){diagonal=Math.max(diagonal,Math.abs(a[i][i]));for(int j=i+1;j<n;j++)if(Math.abs(a[i][j])>maximum){maximum=Math.abs(a[i][j]);p=i;q=j;}}
                if(maximum<1e-14*Math.max(1,diagonal))break;
                double tau=(a[q][q]-a[p][p])/(2*a[p][q]);double t=Math.copySign(1,tau)/(Math.abs(tau)+Math.sqrt(1+tau*tau));double c=1/Math.sqrt(1+t*t),ss=t*c;
                double ap=a[p][p],aq=a[q][q],off=a[p][q];a[p][p]=ap-t*off;a[q][q]=aq+t*off;a[p][q]=a[q][p]=0;
                for(int i=0;i<n;i++)if(i!=p&&i!=q){double ip=a[i][p],iq=a[i][q];a[i][p]=a[p][i]=c*ip-ss*iq;a[i][q]=a[q][i]=ss*ip+c*iq;}
                for(int i=0;i<n;i++){double ip=v[i][p],iq=v[i][q];v[i][p]=c*ip-ss*iq;v[i][q]=ss*ip+c*iq;}
            }
            Integer[] order=new Integer[n];for(int i=0;i<n;i++)order[i]=i;Arrays.sort(order,Comparator.comparingDouble(i->a[i][i]));
            double[] values=new double[n];double[][] vectors=new double[n][n];
            for(int k=0;k<n;k++){values[k]=a[order[k]][order[k]];for(int i=0;i<n;i++)vectors[i][k]=v[i][order[k]];}
            return new RealEigen(values,vectors);
        }
    }
}
