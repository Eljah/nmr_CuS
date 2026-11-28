package org.eljah.nmr.relaxation;

import java.io.*;
import java.nio.charset.StandardCharsets;
import java.nio.file.*;
import java.util.*;
import java.util.function.Predicate;
import static org.eljah.nmr.relaxation.FilesIO.*;

/**
 * Reproducible descriptive analysis of archived 63Cu NQR rates.
 * Java 17 API only. No raw recovery curves or new experimental data are invented.
 */
public final class RelaxationDynamics {
    static final String[] SITES={"s1","s2","s3","s4"};
    static final String[] ALL_SITES={"s1","s2","s3","s4","c"};
    static final String[] SUFFIXES={"1195","1526","1620","1696","1485"};
    static final double Q63_Q65=1.081, GAMMA65_GAMMA63=1.071210656566704;
    static final double T_REF=250.0, THETA_MAX=16000.0, K_B_EV=8.617333262145e-5;
    static final List<String> WEIGHTS=List.of("stored_error_scale","relative_floor_5pct");

    record Rate(String site,double t,double y,double err,String source,int sourceRow) {}
    record Params(double theta,double background,double amplitude) {}
    record Fit(String model,String weighting,Map<String,Params> parameters,double objective) {}
    record Nearest(String label,double distance,boolean complete) {}
    record BppPoint(double logTau,double theta,double objective,Map<String,Params> parameters) {}

    final Path base,inputs,output;
    final List<Rate> r1=new ArrayList<>(),r2=new ArrayList<>(),r2Summary=new ArrayList<>();
    final Map<String,NavigableMap<Double,Double>> frequencies=new LinkedHashMap<>();
    final List<List<Object>> conflicts=new ArrayList<>(),quarantine=new ArrayList<>();
    final Set<Double> ambiguousR2=new TreeSet<>();
    final Map<String,Object> summary=new LinkedHashMap<>();
    final Map<String,Object> audit=new LinkedHashMap<>();

    RelaxationDynamics(Path base){this.base=base;inputs=base.resolve("inputs");output=base.resolve("outputs");}

    public static void main(String[] args)throws Exception {
        Locale.setDefault(Locale.ROOT);
        Path cwd=Path.of("").toAbsolutePath().normalize();
        Path base=Files.exists(cwd.resolve("pom.xml")) && cwd.getFileName().toString().equals("relaxation_dynamics")
                ?cwd:cwd.resolve("theory/relaxation_dynamics");
        for(int j=0;j<args.length;j++)if(args[j].equals("--base"))base=Path.of(args[++j]).toAbsolutePath();
        RelaxationDynamics app=new RelaxationDynamics(base);
        if(Arrays.asList(args).contains("--self-test")){
            System.out.println(json(SelfChecks.run(),0));return;
        }
        if(Arrays.asList(args).contains("--materialize")){
            app.materialize();System.out.println("Materialized exact Git inputs under "+app.inputs);return;
        }
        app.run();
    }

    void materialize()throws Exception {
        Path repository=base.getParent().getParent();
        String revision=new String(git(repository,"rev-parse","HEAD"),StandardCharsets.UTF_8).trim();
        List<Object> manifest=new ArrayList<>();
        for(String object:List.of("02","03","04","05")){
            List<String> names=new ArrayList<>();
            if(object.equals("02"))names.add("Data1.csv");
            if(object.equals("03") || object.equals("04")){
                for(String suffix:SUFFIXES)names.add("T1S4F"+suffix+".csv");
                if(object.equals("03"))names.add("T1S1.csv");
            }
            if(object.equals("05"))names.addAll(List.of("T2s1.csv","T2s2.csv","T2s3.csv","T2s4.csv","T2datd.csv"));
            names.add("column_provenance.json");
            for(String name:names){
                String source="data_analysis/thesis_embedded/numeric/diploma10_object_"+object+"_01/"+name;
                Path dest=inputs.resolve("raw/object_"+object+"/"+name);
                byte[] bytes=git(repository,"show",revision+":"+source);
                Files.createDirectories(dest.getParent());Files.write(dest,bytes);
                manifest.add(map("git_path",source,"local_path",base.relativize(dest).toString().replace('\\','/'),
                        "bytes",bytes.length,"sha256",sha256(bytes)));
            }
        }
        for(String extension:List.of("pdf","txt")){
            String source="new_evidence/thesis/diploma_magister_geerite_2009_5."+extension;
            byte[] bytes=git(repository,"show",revision+":"+source);
            Path dest=inputs.resolve("reference/diploma_magister_geerite_2009_5."+extension);
            Files.createDirectories(dest.getParent());Files.write(dest,bytes);
            manifest.add(map("git_path",source,"local_path",base.relativize(dest).toString().replace('\\','/'),
                    "bytes",bytes.length,"sha256",sha256(bytes)));
        }
        writeJson(inputs.resolve("manifest.json"),map("git_revision",revision,"raw_files",manifest,
                "manual_transcription","isotope_observations.csv: visually checked thesis p.27; 65Cu errors unavailable"));
        StringBuilder digests=new StringBuilder();
        for(Object item:manifest){
            Map<?,?> entry=(Map<?,?>)item;
            digests.append(entry.get("sha256")).append("  ").append(entry.get("local_path")).append('\n');
        }
        digests.append(sha256(Files.readAllBytes(inputs.resolve("isotope_observations.csv"))))
                .append("  inputs/isotope_observations.csv\n");
        Files.writeString(inputs.resolve("manifest.sha256"),digests,StandardCharsets.UTF_8);
    }
    static byte[] git(Path cwd,String...arguments)throws Exception {
        List<String> command=new ArrayList<>();command.add("git");command.addAll(Arrays.asList(arguments));
        Process process=new ProcessBuilder(command).directory(cwd.toFile()).redirectError(ProcessBuilder.Redirect.INHERIT).start();
        byte[] result=process.getInputStream().readAllBytes();
        if(process.waitFor()!=0)throw new IOException("Git command failed: "+command);
        return result;
    }
    Path raw(String object,String filename){return inputs.resolve("raw/object_"+object+"/"+filename);}

    void run()throws Exception {
        if(!Files.exists(raw("04","T1S4F1195.csv")))
            throw new IllegalStateException("Run --materialize once, or provide the preserved inputs directory.");
        Files.createDirectories(output);
        Map<String,Object> checks=SelfChecks.run();writeJson(output.resolve("checks.json"),checks);
        verifyManifest();
        readAndAudit();
        isotopeRatios();
        activationAnalysis();
        contiguousValidation();
        powerAnalysis();
        bppAnalysis();
        r2Analysis();
        writeCsv(output.resolve("provenance_conflicts.csv"),
                List.of("kind","site","T_K","primary_value","alternative_value","primary_error","alternative_error","decision"),conflicts);
        writeCsv(output.resolve("quarantine.csv"),List.of("source","site","T_K","quantity","value","reason"),quarantine);
        writeJson(output.resolve("audit.json"),audit);
        summary.put("audit",audit);
        summary.put("statistical_status","Descriptive objectives and temperature-held-out prediction only. Stored errors are not established standard deviations. No p-values, confidence intervals or unique microscopic barrier claims.");
        summary.put("raw_recovery_available",false);
        writeJson(output.resolve("SUMMARY.json"),summary);
        writeReport();
        System.out.println("Finished. Main results: "+output.resolve("SUMMARY.json"));
    }
    void verifyManifest()throws IOException {
        int checked=0;
        for(String line:Files.readAllLines(inputs.resolve("manifest.sha256"),StandardCharsets.UTF_8)){
            if(line.isBlank())continue;
            String expected=line.substring(0,64),relative=line.substring(66);
            String actual=sha256(Files.readAllBytes(base.resolve(relative)));
            require(actual.equals(expected),"Input SHA-256 mismatch: "+relative);checked++;
        }
        audit.put("input_integrity",map("sha256_files_checked",checked,"status","all exact bytes match the materialized manifest"));
    }

    void readAndAudit()throws Exception {
        List<List<Object>> normalized=new ArrayList<>(),ttChecks=new ArrayList<>();
        int matched=0;double maxDifference=0;
        for(int index=0;index<ALL_SITES.length;index++){
            String site=ALL_SITES[index],name="T1S4F"+SUFFIXES[index]+".csv";
            Map<Double,Rate> byTemperature=new TreeMap<>();
            for(var row:csv(raw("04",name))){
                double t=number(row,"T"),rate=number(row,"W1"),err=number(row,"ErrW1");
                if(!Double.isFinite(t)||!Double.isFinite(rate))continue;
                require(t>0 && rate>0 && err>0,"Invalid R1 rate/error.");
                Rate point=new Rate(site,t,rate,err,"raw/object_04/"+name,(int)number(row,"origin_row"));
                r1.add(point);byTemperature.put(t,point);
                normalized.add(row(site,t,rate,err,rate/t,err/t,1/rate,63,point.source(),point.sourceRow()));
            }
            for(var row:csv(raw("03",name))){
                double t=number(row,"T"),value=number(row,"T1T");
                if(!Double.isFinite(t)||!Double.isFinite(value))continue;
                Rate point=byTemperature.get(t);
                if(point==null){
                    quarantine.add(row("raw/object_03/"+name,site,t,"R1_over_T",value,"No matching obj04 rate or error; implied R1="+value*t));
                }else{
                    double difference=value-point.y()/t;
                    maxDifference=Math.max(maxDifference,Math.abs(difference));matched++;
                    ttChecks.add(row(site,t,value,point.y()/t,difference));
                }
            }
        }
        require(matched==144,"Unexpected T1T matching count: "+matched);
        require(maxDifference<1e-12,"T1T is not equal to R1/T in this source revision.");
        writeCsv(output.resolve("r1_normalized.csv"),
                List.of("site","T_K","R1_s-1","error_scale_s-1","R1_over_T_s-1_K-1","error_over_T_s-1_K-1","T1_s","isotope","source","origin_row"),normalized);
        writeCsv(output.resolve("r1_t1t_audit.csv"),List.of("site","T_K","object03_T1T","object04_R1_over_T","difference"),ttChecks);
        audit.put("T1T",map("meaning","R1/T, not T1*T","matching_points",matched,"max_abs_difference",maxDifference,
                "independent_measurements",false,"extra_s4_290K","quarantined"));
        audit.put("R1",map("isotope",63,"rate_unit","s^-1","stored_errors","Graph17 paired error bars; statistical interpretation unknown",
                "recovery_fit_in_thesis","A*(1-B*exp(-t/T1)), thesis p.16; no individual delay/intensity samples in these tables"));
        audit.put("unassigned_T1S1","Copied for provenance but not assigned to spectral group s1; the worksheet name is not a site identification.");
        audit.put("phase_labels","Thesis uses the historical name geerite for s1-s4. This module preserves spectral labels and does not independently verify the later anilite phase reassignment.");
        audit.put("legacy_activation_sign",map("source_page",28,"printed","R1=alpha*T+beta*exp(+Ea/T)",
                "decision","For positive Ea and beta, this term decreases when heated, contrary to the data; likely missing minus. Implemented negative sign, not literal transcription."));

        for(String site:ALL_SITES)frequencies.put(site,new TreeMap<>());
        int synthetic=0;double maxRatioError=0;
        for(var row:csv(raw("02","Data1.csv"))){
            double t=number(row,"freq");
            if(!Double.isFinite(t)){
                double extra=number(row,"si");if(Double.isFinite(extra))quarantine.add(row("raw/object_02/Data1.csv","unknown",null,"si",extra,"No temperature; not an independent resonance observation"));
                continue;
            }
            for(String site:ALL_SITES){
                double nu=number(row,site),nu65=number(row,site.equals("c")?"ci65":site+"i65");
                if(Double.isFinite(nu))frequencies.get(site).put(t,nu);
                if(Double.isFinite(nu)&&Double.isFinite(nu65)){
                    maxRatioError=Math.max(maxRatioError,Math.abs(nu65-nu/Q63_Q65));synthetic++;
                }
            }
        }
        require(maxRatioError<1e-12,"Unexpected 65Cu frequency generation.");
        audit.put("frequency_table",map("temperature_column_named","freq","unit_T","K","unit_nu","MHz",
                "65Cu_columns","exactly 63Cu/1.081; predictions, not independent measurements",
                "checked_pairs",synthetic,"max_abs_difference_MHz",maxRatioError,
                "default_interpolation","linear only inside <=30 K bracketing interval; never extrapolate",
                "s4_gap","160 to 230 K; 75 K bridge is labelled sensitivity only"));
        for(String site:SITES)for(var row:csv(raw("05","T2"+site+".csv"))){
            double t=number(row,"T"),y=number(row,"B"),err=number(row,"C");
            if(Double.isFinite(t)&&Double.isFinite(y))r2.add(new Rate(site,t,y,err,"raw/object_05/T2"+site+".csv",(int)number(row,"origin_row")));
        }
        for(var row:csv(raw("05","T2datd.csv"))){
            double t=number(row,"A");
            for(String site:ALL_SITES){
                double y=number(row,site.equals("c")?"C":site),err=number(row,site.equals("c")?"cerr":site+"err");
                if(Double.isFinite(t)&&Double.isFinite(y)){
                    Rate p=new Rate(site,t,y,err,"raw/object_05/T2datd.csv",(int)number(row,"origin_row"));
                    r2Summary.add(p);if(site.equals("c"))r2.add(p);
                }
            }
        }
        for(Rate a:r2)if(a.site().equals("s3"))for(Rate b:r2)
            if(b.site().equals("s4")&&a.t()==b.t()&&a.y()==b.y()&&a.err()==b.err()){
                ambiguousR2.add(a.t());
                conflicts.add(row("identical_s3_s4_rate_and_error","s3/s4",a.t(),a.y(),b.y(),a.err(),b.err(),
                        "At most one independent unresolved/ambiguous observation; exclude both or retain s3 only"));
            }
        for(Rate alternate:r2Summary)if(!alternate.site().equals("c")){
            Optional<Rate> primary=r2.stream().filter(p->p.site().equals(alternate.site())&&p.t()==alternate.t()).findFirst();
            if(primary.isEmpty())conflicts.add(row("summary_only",alternate.site(),alternate.t(),null,alternate.y(),null,alternate.err(),"Preserve; only summary_s2 sensitivity uses s2 alternatives"));
            else if(primary.get().y()!=alternate.y()||primary.get().err()!=alternate.err())
                conflicts.add(row("conflicting_value",alternate.site(),alternate.t(),primary.get().y(),alternate.y(),primary.get().err(),alternate.err(),"Separate sensitivity branch; never average"));
        }
        List<List<Object>> normalizedR2=new ArrayList<>();
        for(Rate p:r2)normalizedR2.add(row(p.site(),p.t(),p.y(),p.err(),1000*p.y(),1000*p.err(),1/p.y(),
                ambiguousR2.contains(p.t())&&(p.site().equals("s3")||p.site().equals("s4")),63,p.source(),p.sourceRow()));
        writeCsv(output.resolve("r2_normalized.csv"),List.of("site","T_K","R2_ms-1","error_scale_ms-1","R2_s-1","error_scale_s-1","T2_ms","shared_s3_s4_ambiguous","isotope","source","origin_row"),normalizedR2);
        audit.put("R2",map("unit","ms^-1 (multiply by 1000 for s^-1)","primary_source","Graph6 legend points to T2s1-T2s4; CuS from T2datd",
                "identical_s3_s4_temperatures_K",ambiguousR2,"s2_260K_individual",18.076,"s2_260K_summary",12.55667,
                "sensitivity_policy","Exclude both ambiguous branches or retain s3 only; no test treats them as independent",
                "no_low_T_s1_s2_s4_series",true));
        Map<String,Object> counts=new LinkedHashMap<>();
        for(String site:ALL_SITES)counts.put(site,r1.stream().filter(p->p.site().equals(site)).count());
        audit.put("R1_counts",counts);
    }
    static void require(boolean condition,String text){if(!condition)throw new IllegalStateException(text);}

    void isotopeRatios()throws IOException {
        List<List<Object>> rows=new ArrayList<>();
        double magnetic=GAMMA65_GAMMA63*GAMMA65_GAMMA63,quadrupolar=1/(Q63_Q65*Q63_Q65);
        var observations=csv(inputs.resolve("isotope_observations.csv"));
        for(String site:List.of("s2","s4")){
            double rate63=Double.NaN,rate65=Double.NaN;
            for(var item:observations)if(item.get("site").equals(site)){
                if(number(item,"isotope")==63)rate63=number(item,"R1_s-1");else rate65=number(item,"R1_s-1");
            }
            double ratio=rate65/rate63;
            rows.add(row(site,4.2,rate63,rate65,ratio,magnetic,quadrupolar,100*(ratio/magnetic-1),
                    "No 65Cu error; no uncertainty or mixture fraction inferred"));
        }
        writeCsv(output.resolve("isotope_ratios.csv"),List.of("site","T_K","R1_63_s-1","R1_65_s-1","R65_over_R63","flat_J_magnetic","flat_J_quadrupolar","percent_above_flat_J_magnetic","limitation"),rows);
        List<List<Object>> spectralRows=new ArrayList<>();
        double omega63=2*Math.PI*15e6,omega65=omega63/Q63_Q65;
        for(double log=-13;log<=-3;log+=0.25){
            double tau=Math.pow(10,log);
            double densityRatio=(1+omega63*omega63*tau*tau)/(1+omega65*omega65*tau*tau);
            spectralRows.add(row(log,tau,magnetic*densityRatio,quadrupolar*densityRatio));
        }
        writeCsv(output.resolve("isotope_J_sensitivity.csv"),List.of("log10_tau_s","tau_s","magnetic_R65_over_R63_at_15MHz","quadrupolar_R65_over_R63_at_15MHz"),spectralRows);
        summary.put("isotope",map("s2_R65_R63",5.79/5.0,"s4_R65_R63",2.92/2.5,
                "flat_J_magnetic",magnetic,"flat_J_quadrupolar",quadrupolar,
                "interpretation","Compatible with predominantly magnetic relaxation under a flat-J assumption; 65Cu errors and independent J constraints are absent.",
                "J_sensitivity","single Lorentzian at 15 MHz is an illustration, not a fit or complete I=3/2 recovery model"));
    }

    static double scale(Rate point,String weighting){
        return weighting.equals("stored_error_scale")?point.err():Math.max(point.err(),0.05*point.y());
    }
    static double activation(double t,double theta){return Math.exp(-theta*(1/t-1/T_REF));}
    static double predict(Rate point,Params p){return p.background()*point.t()+p.amplitude()*activation(point.t(),p.theta());}
    static Params fitAt(List<Rate> points,double theta,String weights){
        double[] y=new double[points.size()],s=new double[y.length],x=new double[y.length],z=new double[y.length];
        for(int j=0;j<y.length;j++){Rate p=points.get(j);y[j]=p.y();s[j]=scale(p,weights);x[j]=p.t()/T_REF;z[j]=activation(p.t(),theta);}
        Numerics.Linear fit=Numerics.nnls2(y,s,x,z);
        return new Params(theta,fit.coefficients()[0]/T_REF,fit.coefficients()[1]);
    }
    static double objective(List<Rate> points,Map<String,Params> parameters,String weights){
        double sum=0;for(Rate p:points){double r=(predict(p,parameters.get(p.site()))-p.y())/scale(p,weights);sum+=r*r;}return sum;
    }
    static Map<String,List<Rate>> grouped(List<Rate> points){
        Map<String,List<Rate>> result=new LinkedHashMap<>();
        for(String site:SITES){List<Rate> rows=points.stream().filter(p->p.site().equals(site)).toList();if(!rows.isEmpty())result.put(site,rows);}
        return result;
    }
    static Fit fitActivation(List<Rate> points,String model,String weighting){
        var groups=grouped(points);Map<String,Params> params=new LinkedHashMap<>();
        if(model.equals("linear_T")){
            for(var entry:groups.entrySet()){
                double num=0,den=0;
                for(Rate p:entry.getValue()){double w=1/Math.pow(scale(p,weighting),2);num+=w*p.t()*p.y();den+=w*p.t()*p.t();}
                params.put(entry.getKey(),new Params(0,Math.max(0,num/den),0));
            }
        }else if(model.equals("shared_theta")){
            Numerics.Minimum optimum=Numerics.minimize(theta->{
                Map<String,Params> current=new LinkedHashMap<>();
                for(var e:groups.entrySet())current.put(e.getKey(),fitAt(e.getValue(),theta,weighting));
                return objective(points,current,weighting);
            },0,THETA_MAX,25);
            for(var e:groups.entrySet())params.put(e.getKey(),fitAt(e.getValue(),optimum.x(),weighting));
        }else{
            for(var e:groups.entrySet()){
                Numerics.Minimum optimum=Numerics.minimize(theta->{
                    Params p=fitAt(e.getValue(),theta,weighting);
                    return objective(e.getValue(),Map.of(e.getKey(),p),weighting);
                },0,THETA_MAX,25);
                params.put(e.getKey(),fitAt(e.getValue(),optimum.x(),weighting));
            }
        }
        return new Fit(model,weighting,params,objective(points,params,weighting));
    }
    static Map<Double,Integer> folds(List<Rate> points){
        Set<Double> temperatures=new TreeSet<>();for(Rate p:points)temperatures.add(p.t());
        Map<Double,Integer> result=new TreeMap<>();int k=0;for(double t:temperatures)result.put(t,k++%4);return result;
    }
    static Map<Double,Integer> contiguousFolds(List<Rate> points){
        Set<Double> temperatures=new TreeSet<>();for(Rate p:points)temperatures.add(p.t());
        Map<Double,Integer> result=new TreeMap<>();int k=0,n=temperatures.size();
        for(double t:temperatures)result.put(t,(k++)*4/n);return result;
    }
    static double[] metrics(List<Rate> points,Map<Rate,Double> predictions,String weights){
        double raw=0,relative=0,scaled=0;
        for(Rate p:points){double d=predictions.get(p)-p.y();raw+=d*d;relative+=Math.pow(d/p.y(),2);scaled+=Math.pow(d/scale(p,weights),2);}
        return new double[]{raw,relative,scaled,Numerics.rms(raw,points.size()),Numerics.rms(relative,points.size()),Numerics.rms(scaled,points.size())};
    }
    void activationAnalysis()throws IOException {
        List<Rate> data=r1.stream().filter(p->!p.site().equals("c")&&p.t()>=50).toList();
        var folds=folds(data);
        List<List<Object>> comparison=new ArrayList<>(),parameters=new ArrayList<>(),predictions=new ArrayList<>(),profiles=new ArrayList<>();
        List<Object> summaryFits=new ArrayList<>();
        for(String weighting:WEIGHTS){
            for(String model:List.of("linear_T","shared_theta","separate_theta")){
                Fit fit=fitActivation(data,model,weighting);
                Map<Rate,Double> fitted=new HashMap<>(),heldout=new HashMap<>();
                for(Rate p:data)fitted.put(p,predict(p,fit.parameters().get(p.site())));
                for(int fold=0;fold<4;fold++){
                    final int f=fold;
                    List<Rate> train=data.stream().filter(p->folds.get(p.t())!=f).toList();
                    Fit holdFit=fitActivation(train,model,weighting);
                    for(Rate p:data)if(folds.get(p.t())==fold)heldout.put(p,predict(p,holdFit.parameters().get(p.site())));
                }
                double[] m=metrics(data,fitted,weighting),cv=metrics(data,heldout,weighting);
                int count=model.equals("linear_T")?4:model.equals("shared_theta")?9:12;
                comparison.add(row(weighting,model,data.size(),count,m[2],m[3],m[4],m[5],cv[2],cv[3],cv[4],cv[5]));
                summaryFits.add(map("weighting",weighting,"model",model,"n",data.size(),"parameters",count,
                        "scaled_SSE",m[2],"relative_RMS",m[4],"heldout_relative_RMS",cv[4],"heldout_scaled_RMS",cv[5]));
                for(var e:fit.parameters().entrySet()){
                    Params p=e.getValue();long n=data.stream().filter(v->v.site().equals(e.getKey())).count();
                    double logB=p.amplitude()>0?Math.log(p.amplitude())+p.theta()/T_REF:Double.NaN;
                    parameters.add(row(weighting,model,e.getKey(),n,p.theta(),p.theta()*K_B_EV,p.background(),p.amplitude(),logB,
                            Double.isFinite(logB)?Math.exp(logB):0,p.theta()==0||p.theta()>=THETA_MAX-1));
                }
                for(Rate p:data)predictions.add(row(weighting,model,p.site(),p.t(),p.y(),p.err(),scale(p,weighting),
                        fitted.get(p),heldout.get(p),folds.get(p.t())));
            }
            var groups=grouped(data);
            for(double theta=0;theta<=THETA_MAX;theta+=25){
                double total=0;List<Double> objectives=new ArrayList<>();
                for(String site:SITES){Params p=fitAt(groups.get(site),theta,weighting);double o=objective(groups.get(site),Map.of(site,p),weighting);objectives.add(o);total+=o;}
                profiles.add(row(weighting,theta,total,objectives.get(0),objectives.get(1),objectives.get(2),objectives.get(3)));
            }
        }
        writeCsv(output.resolve("activation_comparison.csv"),List.of("weighting","model","n","parameter_count","scaled_SSE","RMS_s-1","relative_RMS","scaled_RMS","heldout_scaled_SSE","heldout_RMS_s-1","heldout_relative_RMS","heldout_scaled_RMS"),comparison);
        writeCsv(output.resolve("activation_parameters.csv"),List.of("weighting","model","site","n","theta_K","theta_kB_eV","a_s-1_K-1","C_at_250K_s-1","log_B_s-1","B_s-1","theta_at_bound"),parameters);
        writeCsv(output.resolve("activation_predictions.csv"),List.of("weighting","model","site","T_K","observed","stored_error","effective_scale","predicted","heldout_prediction","heldout_fold"),predictions);
        writeCsv(output.resolve("activation_profile.csv"),List.of("weighting","theta_K","shared_scaled_SSE","s1_profile_SSE","s2_profile_SSE","s3_profile_SSE","s4_profile_SSE"),profiles);
        summary.put("activation",map("T_min_K",50,"parameterization","a_i*T + C_i*exp[-theta*(1/T-1/250)]; B_i=C_i*exp(theta/250)",
                "theta_search_K",List.of(0,THETA_MAX),"a_and_C_constraints",">=0","cross_validation","4 folds; sorted unique temperatures assigned index modulo 4; every site at that temperature held out jointly",
                "heldout_columns_mean","interleaved temperature folds, not contiguous temperature blocks",
                "fit_comparison",summaryFits));
        System.out.println("Activation comparison complete: "+comparison.size()+" descriptive models.");
    }
    void contiguousValidation()throws IOException{
        List<Rate> data=r1.stream().filter(p->!p.site().equals("c")&&p.t()>=50).toList();
        Map<Double,Integer> groups=contiguousFolds(data);
        List<List<Object>> comparisons=new ArrayList<>(),predictions=new ArrayList<>(),parameters=new ArrayList<>();
        List<Object> entries=new ArrayList<>(),foldRanges=new ArrayList<>();
        for(int fold=0;fold<4;fold++){
            final int f=fold;
            DoubleSummaryStatistics range=data.stream().filter(p->groups.get(p.t())==f).mapToDouble(Rate::t).summaryStatistics();
            foldRanges.add(map("fold",fold,"heldout_min_T",range.getMin(),"heldout_max_T",range.getMax(),"heldout_n",range.getCount()));
        }
        for(String weighting:WEIGHTS)for(String model:List.of("linear_T","shared_theta","separate_theta")){
            Map<Rate,Double> heldout=new HashMap<>();int atBound=0;
            for(int fold=0;fold<4;fold++){
                final int f=fold;
                List<Rate> train=data.stream().filter(p->groups.get(p.t())!=f).toList();
                List<Rate> test=data.stream().filter(p->groups.get(p.t())==f).toList();
                Fit fit=fitActivation(train,model,weighting);
                double minHeld=test.stream().mapToDouble(Rate::t).min().orElseThrow(),maxHeld=test.stream().mapToDouble(Rate::t).max().orElseThrow();
                for(var e:fit.parameters().entrySet()){
                    Params p=e.getValue();boolean bound=!model.equals("linear_T")&&(p.theta()<=1||p.theta()>=THETA_MAX-1);
                    if(bound)atBound++;
                    parameters.add(row(weighting,model,fold,e.getKey(),minHeld,maxHeld,p.theta(),p.background(),p.amplitude(),bound));
                }
                for(Rate p:test){
                    double prediction=predict(p,fit.parameters().get(p.site()));heldout.put(p,prediction);
                    predictions.add(row(weighting,model,fold,p.site(),p.t(),p.y(),prediction,scale(p,weighting)));
                }
            }
            double[] cv=metrics(data,heldout,weighting);
            comparisons.add(row(weighting,model,data.size(),cv[2],cv[3],cv[4],cv[5],atBound));
            entries.add(map("weighting",weighting,"model",model,"n",data.size(),"contiguous_heldout_scaled_SSE",cv[2],
                    "contiguous_heldout_relative_RMS",cv[4],"contiguous_heldout_scaled_RMS",cv[5],"site_fold_theta_bound_count",atBound));
        }
        writeCsv(output.resolve("activation_contiguous_cv_comparison.csv"),List.of("weighting","model","n","contiguous_heldout_scaled_SSE","contiguous_heldout_RMS_s-1","contiguous_heldout_relative_RMS","contiguous_heldout_scaled_RMS","site_fold_theta_bound_count"),comparisons);
        writeCsv(output.resolve("activation_contiguous_cv_predictions.csv"),List.of("weighting","model","fold","site","T_K","observed","heldout_prediction","effective_scale"),predictions);
        writeCsv(output.resolve("activation_contiguous_cv_parameters.csv"),List.of("weighting","model","fold","site","heldout_T_min_K","heldout_T_max_K","theta_K","a_s-1_K-1","C_at_250K_s-1","theta_at_bound"),parameters);
        summary.put("activation_contiguous_validation",map("scheme","sorted unique temperatures assigned floor(4*index/N); all sites at a held-out temperature excluded together",
                "fold_ranges",foldRanges,"results",entries,"interpretation","Includes extrapolation outside training-temperature support; report separately from interleaved validation."));
        System.out.println("Contiguous-temperature validation complete.");
    }
    static Params powerAt(List<Rate> points,String weights,double exponent){
        int n=points.size();double[] y=new double[n],s=new double[n],first=new double[n],second=new double[n];
        for(int j=0;j<n;j++){Rate p=points.get(j);y[j]=p.y();s[j]=scale(p,weights);first[j]=p.t()/T_REF;second[j]=Math.pow(p.t()/T_REF,exponent);}
        var fit=Numerics.nnls2(y,s,first,second);
        return new Params(exponent,fit.coefficients()[0]/T_REF,fit.coefficients()[1]);
    }
    static double predictPower(Rate point,Params parameter){
        return parameter.background()*point.t()+parameter.amplitude()*Math.pow(point.t()/T_REF,parameter.theta());
    }
    static Fit fitPower(List<Rate> data,String weights){
        var grouped=grouped(data);
        Numerics.Minimum best=Numerics.minimize(exponent->{
            double sum=0;
            for(var entry:grouped.entrySet()){
                Params p=powerAt(entry.getValue(),weights,exponent);
                for(Rate point:entry.getValue())sum+=Math.pow((predictPower(point,p)-point.y())/scale(point,weights),2);
            }
            return sum;
        },1.2,24,0.05);
        Map<String,Params> parameters=new LinkedHashMap<>();
        for(var entry:grouped.entrySet())parameters.put(entry.getKey(),powerAt(entry.getValue(),weights,best.x()));
        return new Fit("shared_power",weights,parameters,best.value());
    }
    void powerAnalysis()throws IOException{
        List<Rate> data=r1.stream().filter(p->!p.site().equals("c")&&p.t()>=50).toList();
        List<List<Object>> comparisons=new ArrayList<>(),parameters=new ArrayList<>(),predictions=new ArrayList<>(),profiles=new ArrayList<>(),cvParameters=new ArrayList<>();
        List<Object> summaries=new ArrayList<>();
        for(String weighting:WEIGHTS){
            Fit fit=fitPower(data,weighting);
            Map<Rate,Double> fitted=new HashMap<>(),interleaved=new HashMap<>(),contiguous=new HashMap<>();
            for(Rate p:data)fitted.put(p,predictPower(p,fit.parameters().get(p.site())));
            for(String scheme:List.of("interleaved","contiguous")){
                Map<Double,Integer> folds=scheme.equals("interleaved")?folds(data):contiguousFolds(data);
                for(int fold=0;fold<4;fold++){
                    final int f=fold;Fit train=fitPower(data.stream().filter(p->folds.get(p.t())!=f).toList(),weighting);
                    for(var entry:train.parameters().entrySet()){
                        Params p=entry.getValue();
                        cvParameters.add(row(weighting,scheme,fold,entry.getKey(),p.theta(),p.background(),p.amplitude(),p.theta()<=1.200001||p.theta()>=23.999999));
                    }
                    for(Rate p:data)if(folds.get(p.t())==fold)
                        (scheme.equals("interleaved")?interleaved:contiguous).put(p,predictPower(p,train.parameters().get(p.site())));
                }
            }
            double[] train=metrics(data,fitted,weighting),inter=metrics(data,interleaved,weighting),contig=metrics(data,contiguous,weighting);
            double exponent=fit.parameters().get("s1").theta();
            comparisons.add(row(weighting,"shared_power",data.size(),9,exponent,train[2],train[3],train[4],inter[2],inter[3],inter[4],contig[2],contig[3],contig[4]));
            summaries.add(map("weighting",weighting,"model","shared_power","parameters",9,"exponent",exponent,"scaled_SSE",train[2],
                    "relative_RMS",train[4],"interleaved_heldout_relative_RMS",inter[4],"contiguous_heldout_relative_RMS",contig[4]));
            for(var entry:fit.parameters().entrySet()){
                Params p=entry.getValue();parameters.add(row(weighting,entry.getKey(),p.theta(),p.background(),p.amplitude(),p.theta()<=1.200001||p.theta()>=23.999999));
            }
            for(Rate p:data)predictions.add(row(weighting,p.site(),p.t(),p.y(),fitted.get(p),interleaved.get(p),contiguous.get(p)));
            for(int j=0;j<=456;j++){
                double power=1.2+0.05*j,total=0;
                for(var entry:grouped(data).entrySet()){
                    Params p=powerAt(entry.getValue(),weighting,power);
                    for(Rate point:entry.getValue())total+=Math.pow((predictPower(point,p)-point.y())/scale(point,weighting),2);
                }
                profiles.add(row(weighting,power,total));
            }
        }
        writeCsv(output.resolve("power_comparison.csv"),List.of("weighting","model","n","parameter_count","power_p","scaled_SSE","RMS_s-1","relative_RMS","interleaved_heldout_scaled_SSE","interleaved_heldout_RMS_s-1","interleaved_heldout_relative_RMS","contiguous_heldout_scaled_SSE","contiguous_heldout_RMS_s-1","contiguous_heldout_relative_RMS"),comparisons);
        writeCsv(output.resolve("power_parameters.csv"),List.of("weighting","site","power_p","a_s-1_K-1","C_at_250K_s-1","power_at_bound"),parameters);
        writeCsv(output.resolve("power_predictions.csv"),List.of("weighting","site","T_K","observed","predicted","interleaved_heldout_prediction","contiguous_heldout_prediction"),predictions);
        writeCsv(output.resolve("power_profile.csv"),List.of("weighting","power_p","shared_scaled_SSE"),profiles);
        writeCsv(output.resolve("power_cv_parameters.csv"),List.of("weighting","scheme","fold","site","power_p","a_s-1_K-1","C_at_250K_s-1","power_at_bound"),cvParameters);
        summary.put("power_alternative",map("model","R1_i=a_i*T+C_i*(T/250)^p; common p in [1.2,24], a_i,C_i>=0",
                "parameter_count",9,"interpretation","Phenomenological curvature control, not an established Raman process","results",summaries));
        System.out.println("Equal-parameter power-law alternative complete.");
    }

    double frequency(String site,double temperature,double maxGap){
        var series=frequencies.get(site);
        if(series.containsKey(temperature))return series.get(temperature);
        var left=series.floorEntry(temperature);var right=series.ceilingEntry(temperature);
        if(left==null||right==null||right.getKey()-left.getKey()>maxGap)return Double.NaN;
        return left.getValue()+(right.getValue()-left.getValue())*(temperature-left.getKey())/(right.getKey()-left.getKey());
    }
    static double bppShape(double nuMHz,double t,double logTau,double theta){
        double omega=2*Math.PI*nuMHz*1e6,omegaRef=2*Math.PI*15e6;
        double u=Math.log(omega)+logTau*Math.log(10)+theta/t;
        double a=Math.exp(-Math.abs(u));
        return (omegaRef/omega)*a/(1+a*a);
    }
    BppPoint fitBppAt(List<Rate> data,String weighting,double logTau,double theta){
        Map<String,Params> parameters=new LinkedHashMap<>();double total=0;
        for(var e:grouped(data).entrySet()){
            List<Rate> points=e.getValue();int n=points.size();
            double[] y=new double[n],scale=new double[n],x=new double[n],z=new double[n];
            for(int j=0;j<n;j++){Rate p=points.get(j);y[j]=p.y();scale[j]=scale(p,weighting);x[j]=p.t()/T_REF;z[j]=bppShape(frequency(p.site(),p.t(),30),p.t(),logTau,theta);}
            var fit=Numerics.nnls2(y,scale,x,z);total+=fit.objective();
            parameters.put(e.getKey(),new Params(theta,fit.coefficients()[0]/T_REF,fit.coefficients()[1]));
        }
        return new BppPoint(logTau,theta,total,parameters);
    }
    void bppAnalysis()throws IOException {
        List<Rate> data=r1.stream().filter(p->!p.site().equals("c")&&p.t()>=50&&Double.isFinite(frequency(p.site(),p.t(),30))).toList();
        List<List<Object>> surface=new ArrayList<>(),tauProfiles=new ArrayList<>(),paramsRows=new ArrayList<>(),predictionRows=new ArrayList<>();
        List<Object> summaries=new ArrayList<>();
        for(String weights:WEIGHTS){
            List<BppPoint> profile=new ArrayList<>();
            for(int j=0;j<=52;j++){
                double logTau=-18+j*0.25;
                if(weights.equals("stored_error_scale"))for(double theta=0;theta<=THETA_MAX;theta+=100){
                    BppPoint p=fitBppAt(data,weights,logTau,theta);surface.add(row(logTau,theta,p.objective()));
                }
                Numerics.Minimum minimum=Numerics.minimize(theta->fitBppAt(data,weights,logTau,theta).objective(),0,THETA_MAX,50);
                profile.add(fitBppAt(data,weights,logTau,minimum.x()));
            }
            BppPoint best=profile.stream().min(Comparator.comparingDouble(BppPoint::objective)).orElseThrow();
            double low=Double.POSITIVE_INFINITY,high=Double.NEGATIVE_INFINITY;
            for(BppPoint p:profile){
                boolean near=p.objective()<=best.objective()*1.01;
                if(near){low=Math.min(low,p.logTau());high=Math.max(high,p.logTau());}
                double minU=Double.POSITIVE_INFINITY,maxU=Double.NEGATIVE_INFINITY;
                for(Rate r:data){double u=Math.log(2*Math.PI*frequency(r.site(),r.t(),30)*1e6)+p.logTau()*Math.log(10)+p.theta()/r.t();minU=Math.min(minU,u);maxU=Math.max(maxU,u);}
                tauProfiles.add(row(weights,p.logTau(),p.theta(),p.objective(),p.objective()/best.objective(),near,minU,maxU));
            }
            Fit activationReference=fitActivation(data,"shared_theta",weights);
            summaries.add(map("weighting",weights,"n",data.size(),"best_log10_tau0_s",best.logTau(),"best_theta_K",best.theta(),
                    "best_scaled_SSE",best.objective(),"activation_shared_same_subset_SSE",activationReference.objective(),
                    "within_1pct_objective_logtau_range",List.of(low,high),
                    "one_percent_is","arbitrary descriptive tolerance, not a confidence region",
                    "grid_touches_boundary",low<=-18||high>=-5));
            for(var e:best.parameters().entrySet())paramsRows.add(row(weights,e.getKey(),best.logTau(),best.theta(),e.getValue().background(),e.getValue().amplitude()));
            for(Rate p:data){Params a=best.parameters().get(p.site());double prediction=a.background()*p.t()+a.amplitude()*bppShape(frequency(p.site(),p.t(),30),p.t(),best.logTau(),best.theta());predictionRows.add(row(weights,p.site(),p.t(),p.y(),prediction));}
        }
        writeCsv(output.resolve("bpp_surface.csv"),List.of("log10_tau0_s","theta_K","scaled_SSE_stored_error_scale"),surface);
        writeCsv(output.resolve("bpp_tau0_profile.csv"),List.of("weighting","log10_tau0_s","profiled_theta_K","profiled_scaled_SSE","objective_over_min","within_one_percent_descriptive","min_log_omega_tau","max_log_omega_tau"),tauProfiles);
        writeCsv(output.resolve("bpp_parameters.csv"),List.of("weighting","site","log10_tau0_s","theta_K","a_s-1_K-1","dimensionless_J_amplitude_s-1"),paramsRows);
        writeCsv(output.resolve("bpp_predictions.csv"),List.of("weighting","site","T_K","observed_R1_s-1","predicted_R1_s-1"),predictionRows);
        summary.put("BPP",map("model","R1_i=a_i*T+A_i*omega_ref*tau/(1+(omega_i(T)*tau)^2); tau=10^logtau0*exp(theta/T); omega_ref=2*pi*15MHz",
                "interpretation","One-frequency Lorentzian spectral-density phenomenology. Not an exact I=3/2 recovery model or identification of magnetic versus quadrupolar motion.",
                "range_log10_tau0_s",List.of(-18,-5),"theta_range_K",List.of(0,THETA_MAX),"profiles",summaries,
                "excluded_s3_290K","No bracketing measured frequency; not extrapolated",
                "slow_motion_degeneracy","For omega*tau >>1, R1_fluct ~ A*omega_ref/(omega^2*tau0)*exp(-theta/T); A and tau0 cannot be separated."));
        System.out.println("BPP profiles complete.");
    }

    Nearest nearest(Rate point,boolean covellitePartners,double gap){
        double target=frequency(point.site(),point.t(),gap),minimum=Double.POSITIVE_INFINITY;
        String label="";boolean complete=Double.isFinite(target);
        for(String partner:ALL_SITES){
            if(partner.equals("c")!=covellitePartners)continue;
            double nu=frequency(partner,point.t(),gap);
            if(!Double.isFinite(nu)){complete=false;continue;}
            for(int isotope:List.of(63,65)){
                if(partner.equals(point.site())&&isotope==63)continue;
                double d=Math.abs(target-nu/(isotope==63?1:Q63_Q65));
                if(d<minimum){minimum=d;label=partner+"_"+isotope;}
            }
        }
        return new Nearest(label,Double.isFinite(minimum)?minimum:Double.NaN,complete);
    }
    static double overlap(double distance,double width){return 1/(1+Math.pow(distance/width,2));}
    static double hinge(double t){return Math.max(0,t-210)/70;}

    void r2Analysis()throws IOException {
        List<List<Object>> nearestRows=new ArrayList<>(),pairs=new ArrayList<>();
        for(Rate p:r2){
            Nearest host=nearest(p,false,30),cov=nearest(p,true,30);
            nearestRows.add(row(p.site(),p.t(),p.y(),p.err(),frequency(p.site(),p.t(),30),host.label(),host.distance(),host.complete(),
                    cov.label(),cov.distance(),cov.complete(),ambiguousR2.contains(p.t())&&(p.site().equals("s3")||p.site().equals("s4"))));
            double target=frequency(p.site(),p.t(),30);
            for(String partner:ALL_SITES)for(int isotope:List.of(63,65)){
                if(partner.equals(p.site())&&isotope==63)continue;
                double nu=frequency(partner,p.t(),30)/(isotope==63?1:Q63_Q65);
                pairs.add(row(p.site(),p.t(),63,target,partner,isotope,nu,Math.abs(target-nu),
                        p.site().equals("c")!=partner.equals("c"),isotope==65?"generated_divide_1.081":"measured_or_local_interpolation"));
            }
        }
        writeCsv(output.resolve("r2_nearest_lines.csv"),List.of("site","T_K","R2_ms-1","error_scale_ms-1","observed_nu63_MHz","nearest_host_line","host_detuning_MHz","all_host_frequencies_known","nearest_CuS_line","CuS_detuning_MHz","CuS_frequency_known","ambiguous_shared_s3_s4"),nearestRows);
        writeCsv(output.resolve("r2_all_line_pairs.csv"),List.of("observed_site","T_K","observed_isotope","observed_nu_MHz","partner_site","partner_isotope","partner_nu_MHz","detuning_MHz","interphase_comparison","partner_provenance"),pairs);

        List<List<Object>> comparisons=new ArrayList<>(),predictions=new ArrayList<>(),parameters=new ArrayList<>();
        List<Object> mainResults=new ArrayList<>();
        for(String source:List.of("individual","summary_s2")){
            List<Rate> sourceData=new ArrayList<>();
            for(Rate p:r2)if(!p.site().equals("c")&&(!source.equals("summary_s2")||!p.site().equals("s2")))sourceData.add(p);
            if(source.equals("summary_s2"))for(Rate p:r2Summary)if(p.site().equals("s2"))sourceData.add(p);
            for(String ambiguity:List.of("exclude_both","keep_s3")){
                for(double gap:List.of(30.0,75.0)){
                    List<Rate> data=sourceData.stream().filter(p->p.t()>=77)
                            .filter(p->!ambiguousR2.contains(p.t())||!p.site().equals("s4"))
                            .filter(p->!ambiguousR2.contains(p.t())||!p.site().equals("s3")||ambiguity.equals("keep_s3"))
                            .filter(p->nearest(p,false,gap).complete()).sorted(Comparator.comparing(Rate::site).thenComparingDouble(Rate::t)).toList();
                    Map<Double,Integer> folds=folds(data);
                    for(double width:List.of(0.05,0.10,0.20))for(String model:List.of("site_offsets","common_hinge","proximity","hinge_plus_proximity","site_hinges")){
                        Numerics.Linear fit=fitR2(data,model,width,gap);
                        Map<Rate,Double> full=new HashMap<>(),heldout=new HashMap<>();
                        for(Rate p:data)full.put(p,predictR2(p,fit.coefficients(),model,width,gap));
                        for(int fold=0;fold<4;fold++){
                            final int f=fold;List<Rate> train=data.stream().filter(p->folds.get(p.t())!=f).toList();
                            Numerics.Linear trainingFit=fitR2(train,model,width,gap);
                            for(Rate p:data)if(folds.get(p.t())==fold)heldout.put(p,predictR2(p,trainingFit.coefficients(),model,width,gap));
                        }
                        double[] m=metrics(data,full,"relative_floor_5pct"),cv=metrics(data,heldout,"relative_floor_5pct");
                        comparisons.add(row(source,ambiguity,gap,width,model,data.size(),fit.coefficients().length,m[2],m[3],m[4],cv[2],cv[3],cv[4]));
                        for(int j=0;j<fit.coefficients().length;j++)parameters.add(row(source,ambiguity,gap,width,model,j,fit.coefficients()[j]));
                        for(Rate p:data)predictions.add(row(source,ambiguity,gap,width,model,p.site(),p.t(),p.y(),full.get(p),heldout.get(p),folds.get(p.t()),nearest(p,false,gap).distance()));
                        if(source.equals("individual")&&gap==30&&width==0.10)
                            mainResults.add(map("ambiguity",ambiguity,"model",model,"n",data.size(),"scaled_SSE",m[2],
                                    "relative_RMS",m[4],"heldout_relative_RMS",cv[4]));
                    }
                }
            }
        }
        writeCsv(output.resolve("r2_model_comparison.csv"),List.of("source_variant","ambiguity_policy","max_frequency_interpolation_gap_K","overlap_width_MHz","model","n","parameter_count","scaled_SSE","RMS_ms-1","relative_RMS","heldout_scaled_SSE","heldout_RMS_ms-1","heldout_relative_RMS"),comparisons);
        writeCsv(output.resolve("r2_predictions.csv"),List.of("source_variant","ambiguity_policy","max_frequency_interpolation_gap_K","overlap_width_MHz","model","site","T_K","observed_R2_ms-1","predicted_R2_ms-1","heldout_prediction_R2_ms-1","heldout_fold","nearest_host_detuning_MHz"),predictions);
        writeCsv(output.resolve("r2_parameters.csv"),List.of("source_variant","ambiguity_policy","max_frequency_interpolation_gap_K","overlap_width_MHz","model","coefficient_index","value"),parameters);
        List<Rate> cov=r2.stream().filter(p->p.site().equals("c")&&p.t()>=77).toList();
        double[] y=new double[cov.size()],sc=new double[y.length],one=new double[y.length],h=new double[y.length];
        for(int j=0;j<y.length;j++){Rate p=cov.get(j);y[j]=p.y();sc[j]=scale(p,"relative_floor_5pct");one[j]=1;h[j]=hinge(p.t());}
        Numerics.Linear covFit=Numerics.nnls2(y,sc,one,h);
        summary.put("R2",map("primary_model_domain","T>=77K; all four host frequency branches known within a <=30 K interpolation bracket",
                "models","site offsets; offsets+common high-T hinge; offsets+proximity; both; site-specific hinges",
                "hinge","max(0,T-210)/70; phenomenological, fixed onset; no transition temperature fitted",
                "proximity","1/(1+(nearest within-host 63/65 detuning / width)^2), widths 0.05,0.10,0.20 MHz",
                "cross_phase_CuS","Listed separately, not fitted as microscopic interphase coupling",
                "weighting","max(stored error, 0.05*R2), used as descriptive scale",
                "ambiguity_policies",List.of("exclude_both","keep_s3"),
                "main_results",mainResults,"covellite_control",map("n",cov.size(),"offset_ms-1",covFit.coefficients()[0],"nonnegative_hinge_ms-1",covFit.coefficients()[1],"scaled_SSE",covFit.objective()),
                "causal_limit","Proximity is not proof of cross-relaxation; unresolved lines, static disorder and common temperature dependence remain alternatives."));
        System.out.println("R2 provenance sensitivities complete: "+comparisons.size()+" model comparisons.");
    }
    double[] design(Rate p,String model,double width,double gap){
        int count=model.equals("site_offsets")?4:model.equals("hinge_plus_proximity")?6:model.equals("site_hinges")?8:5;
        double[] x=new double[count];int site=Arrays.asList(SITES).indexOf(p.site());x[site]=1;
        if(model.equals("common_hinge")||model.equals("hinge_plus_proximity"))x[4]=hinge(p.t());
        if(model.equals("proximity"))x[4]=overlap(nearest(p,false,gap).distance(),width);
        if(model.equals("hinge_plus_proximity"))x[5]=overlap(nearest(p,false,gap).distance(),width);
        if(model.equals("site_hinges"))x[4+site]=hinge(p.t());
        return x;
    }
    Numerics.Linear fitR2(List<Rate> points,String model,double width,double gap){
        double[][] x=new double[points.size()][];double[] y=new double[x.length],s=new double[x.length];
        for(int j=0;j<x.length;j++){Rate p=points.get(j);x[j]=design(p,model,width,gap);y[j]=p.y();s[j]=scale(p,"relative_floor_5pct");}
        return Numerics.nnls(x,y,s);
    }
    double predictR2(Rate p,double[] coefficients,String model,double width,double gap){
        double[] x=design(p,model,width,gap);double y=0;for(int j=0;j<x.length;j++)y+=x[j]*coefficients[j];return y;
    }
    void writeReport()throws IOException {
        String text="""
                # Проверка релаксационных и динамических гипотез

                Расчёт полностью выполнен Java 17. Численные результаты и подробный
                комментарий находятся в ANALYSIS_RU.md, SUMMARY.json и CSV этого каталога.
                Исходные CSV/JSON/PDF сохранены без изменения в ../inputs; manifest.json
                фиксирует Git-ревизию и SHA-256 каждого исходника.

                Главные ограничения: сохранены параметры прежних одноэкспоненциальных
                подгонок, а не исходные M(t); ошибки не установлены как стандартные
                отклонения; 65Cu частоты получены делением на 1.081. Два фактических
                измерения 65Cu R1 при 4.2 K транскрибированы отдельно со стр.27 диплома.
                Совпадающий блок s3/s4 R2 никогда не считается двумя независимыми
                сериями. Никаких p-values, доверительных интервалов или доказанного
                единственного барьера расчёт не выдаёт.
                """;
        Files.writeString(output.resolve("README_RU.md"),text,StandardCharsets.UTF_8);
    }
}
