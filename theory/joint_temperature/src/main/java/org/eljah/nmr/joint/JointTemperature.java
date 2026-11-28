package org.eljah.nmr.joint;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.security.MessageDigest;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HexFormat;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.function.Predicate;

/**
 * Joint analysis of archived 63Cu NQR frequencies, Java 17, no external libraries.
 *
 * No statistical error bars are supplied: every residual is descriptive and no
 * p-value, chi-square significance, or confidence interval is inferred.
 * The period-four calculation tests a fixed tensor pattern and one real scalar
 * modulation amplitude. Its constrained construction is a CONDITIONAL
 * projection, not a global nonlinear least-squares fit of all tensor elements.
 */
public final class JointTemperature {
    static final String EXPECTED_INPUT_SHA="6b735e3026e15b680811b126e0ea22fd5f38b2c78f1fb3d942626162e1a9f415";
    static final String EXPECTED_SOURCE_SHA="cd331631cb3f904dfa38ef6f61e403f39c021678019a1103763cad36851855ec";
    static final String SOURCE_COMMIT="cdbbda92a2d924daff508f912fb23d116fd70728";
    static final String[] SITES={"s1","s2","s3","s4"};
    static final String[] SCENARIOS={"original","swap_s3_s4_T_ge_230K"};
    static final String[] MODELS={"bose_common_theta","bose_separate_theta","T_squared"};
    static final int[][] ORDERS={{0,1,3,2},{0,1,2,3},{0,2,1,3}};
    static final String[] ORDER_NAMES={"opposite_14_23","opposite_13_24","opposite_12_34"};
    static final double LOG_THETA_LOW=Math.log(1),LOG_THETA_HIGH=Math.log(3000);
    static final int OPTIMIZER_GRID=400;
    static final List<Block> BLOCKS=List.of(
            new Block("cold_4p2_60K",4.2,60), new Block("middle_77_160K",77,160),
            new Block("gap_180_220K",180,220),new Block("hot_230_290K",230,290));

    record Data(double t,double[] nu,int sourceRow) { boolean complete(){return Arrays.stream(nu).allMatch(Double::isFinite);} }
    record Block(String name,double low,double high) { boolean contains(double t){return t>=low-1e-9&&t<=high+1e-9;} }
    record Score(int n,double sse,double max) { double rms(){return n==0?Double.NaN:Math.sqrt(sse/n);} }
    record Thermal(String model,double[] theta,double[] a,double[] b,boolean[] boundary) {
        double predict(double t,int site){return a[site]+b[site]*basis(model,t,theta[site]);}
    }
    record Harmonic(Data data,double mean,double b,double c,double d) {
        double power(){return b*b+c*c;}
    }
    record TensorFit(double ux,double uy,double[] singular,double firstRmsMHz2,
                     double alpha,double beta,double identityCondition,double identityRmsMHz2,
                     int dSign,boolean dSignConstant,boolean rawGramFeasible,
                     double physicalA,double physicalE,boolean physicalBoundary,
                     double physicalObjective) {}
    record TensorScore(int n,double firstSseMHz,double firstMaxMHz,double identitySseMHz2,
                       double conditionalSseMHz,double conditionalMaxMHz,int signMismatch) {
        double firstRms(){return n==0?Double.NaN:Math.sqrt(firstSseMHz/(4*n));}
        double identityRms(){return n==0?Double.NaN:Math.sqrt(identitySseMHz2/n);}
        double conditionalRms(){return n==0?Double.NaN:Math.sqrt(conditionalSseMHz/(4*n));}
    }
    record Centered(double[][] matrix,double[][] removed) {}

    static final class Csv implements AutoCloseable {
        final BufferedWriter writer;
        Csv(Path path,String...header)throws IOException{writer=Files.newBufferedWriter(path,StandardCharsets.UTF_8);row((Object[])header);}
        void row(Object...cells)throws IOException{
            for(int i=0;i<cells.length;i++){
                if(i>0)writer.write(',');
                Object o=cells[i];String s;
                if(o==null)s="";
                else if(o instanceof Double d)s=Double.isFinite(d)?Double.toString(d):"";
                else s=o.toString();
                if(s.indexOf(',')>=0||s.indexOf('"')>=0||s.indexOf('\n')>=0)writer.write('"'+s.replace("\"","\"\"")+'"');
                else writer.write(s);
            }writer.newLine();
        }
        public void close()throws IOException{writer.close();}
    }
    final Path input,out;
    final Map<String,Csv> csv=new LinkedHashMap<>();
    final List<String> failed=new ArrayList<>();
    JointTemperature(Path input,Path out){this.input=input;this.out=out;}

    public static void main(String[] args)throws Exception{
        if(args.length>2)throw new IllegalArgumentException("Usage: JointTemperature [input.csv [output-directory]]");
        Path input=Path.of(args.length>0?args[0]:"inputs/temperature_series.csv");
        Path output=Path.of(args.length>1?args[1]:"results");
        new JointTemperature(input,output).run();
    }

    void run()throws Exception{
        Locale.setDefault(Locale.ROOT);Files.createDirectories(out);
        openOutputs();
        try{
            List<Data> source=read(input);
            validate(source);
            for(String scenario:SCENARIOS){
                List<Data> data=scenario(source,scenario);
                writeTransformed(scenario,data);
                rankAnalyses(scenario,data);
                thermalAnalyses(scenario,data);
                tensorAnalyses(scenario,data);
            }
            writeManifest(source);
        }finally{for(Csv c:csv.values())c.close();}
        if(!failed.isEmpty())throw new IllegalStateException("Validation failed: "+failed);
        System.out.println("Joint-temperature analysis completed: "+out.toAbsolutePath());
        System.out.println("24 temperatures; 89 finite observations; 18 complete quartets; two fixed label scenarios.");
        System.out.println("No experimental uncertainties; no p-values. See ANALYSIS_RU.md and the output CSV files.");
    }

    void add(String name,String...header)throws IOException{csv.put(name,new Csv(out.resolve(name+".csv"),header));}
    void openOutputs()throws IOException{
        add("transformed_series","labels","T_K","site","frequency_MHz","squared_MHz2","row_centered_MHz","row_centered_MHz2","complete_quartet","source_row");
        add("rank_summary","labels","subset","quantity","centering","rank","n_temperatures","n_cells","s1","s2","s3","s4","retained_fraction","residual_rms","max_abs_residual","frequency_equivalent_rms_MHz");
        add("rank_modes","labels","quantity","centering","mode","site","right_vector","singular_value");
        add("rank_reconstructions","labels","quantity","centering","rank","T_K","site","observed","fitted","residual","frequency_fitted_MHz","source_row");
        add("rank_omit_one","labels","omitted_T_K","nu_row_rank1_rms_MHz","nu_column_rank1_rms_MHz","nu2_row_rank2_rms_MHz2");
        add("thermal_parameters","labels","model","fold","site","n_train","theta_K","a_MHz","b_MHz","theta_boundary");
        add("thermal_scores","labels","model","fold","n_train","n_test","parameters","train_rms_MHz","test_rms_MHz","train_max_abs_MHz","test_max_abs_MHz");
        add("thermal_predictions","labels","model","fold","partition","T_K","site","observed_MHz","fitted_MHz","residual_MHz","source_row");
        add("tensor_coefficients","labels","order","cyclic_sites","T_K","xmean_MHz2","B_MHz2","F_MHz2","D_MHz2","F_over_B","P_MHz4","source_row");
        add("tensor_summary","labels","order","cyclic_sites","subset","n_temperatures","theta_direction_deg","p_cos_squared","first_singular_1","first_singular_2","first_rms_MHz2","first_equivalent_rms_MHz","ratio_min","ratio_max","direction_spread_deg","D_min_MHz2","D_max_MHz2","D_constant_sign","identity_alpha","identity_beta","identity_design_condition","identity_rms_MHz2","raw_Gram_feasible_at_TLS_p","physical_A","physical_E","physical_boundary","conditional_physical_rms_MHz","conditional_physical_max_MHz","Gram_margin_cos","Gram_margin_sin","Gram_boundary_active");
        add("tensor_predictions","labels","order","fit_subset","partition","T_K","site","observed_MHz","first_only_fitted_MHz","conditional_physical_fitted_MHz","conditional_residual_MHz","predicted_mean_MHz2","observed_D_MHz2","latent_v_gauge","latent_psi_gauge","source_row");
        add("tensor_block_validation","labels","order","held_out_block","n_train","n_test","theta_direction_deg","identity_alpha","identity_beta","physical_A","physical_E","test_first_rms_MHz","test_identity_rms_MHz2","test_conditional_physical_rms_MHz","test_D_sign_mismatch");
        add("tensor_omit_one","labels","order","omitted_T_K","theta_direction_deg","identity_alpha","identity_beta","physical_A","physical_E","train_first_rms_MHz","train_identity_rms_MHz2","train_conditional_physical_rms_MHz","omitted_first_rms_MHz","omitted_identity_rms_MHz2","omitted_conditional_physical_rms_MHz");
        add("validation_checks","check","value","expected_or_limit","pass","detail");
    }

    static List<Data> read(Path path)throws IOException{
        List<String> lines=Files.readAllLines(path,StandardCharsets.UTF_8);
        if(lines.isEmpty()||!lines.get(0).equals("T_K,s1_MHz,s2_MHz,s3_MHz,s4_MHz,source_row"))throw new IllegalArgumentException("Unexpected input header");
        List<Data> data=new ArrayList<>();double previous=-1;
        for(int i=1;i<lines.size();i++){
            if(lines.get(i).isBlank())continue;
            String[] f=lines.get(i).split(",",-1);if(f.length!=6)throw new IllegalArgumentException("Input line "+(i+1));
            double t=Double.parseDouble(f[0]);double[] v=new double[4];
            if(!(t>previous))throw new IllegalArgumentException("Temperatures must be strictly increasing");previous=t;
            for(int j=0;j<4;j++){
                v[j]=f[j+1].isEmpty()?Double.NaN:Double.parseDouble(f[j+1]);
                if(Double.isFinite(v[j])&&v[j]<=0)throw new IllegalArgumentException("Positive frequencies required");
            }
            data.add(new Data(t,v,Integer.parseInt(f[5])));
        }return data;
    }
    static List<Data> scenario(List<Data> input,String name){
        List<Data> r=new ArrayList<>();
        for(Data d:input){double[] nu=d.nu.clone();if(name.equals(SCENARIOS[1])&&d.t>=230){double x=nu[2];nu[2]=nu[3];nu[3]=x;}r.add(new Data(d.t,nu,d.sourceRow));}
        return r;
    }
    static List<Data> select(List<Data> data,Predicate<Data> predicate){return data.stream().filter(predicate).toList();}
    static String orderText(int[] order){return String.join(" ",Arrays.stream(order).mapToObj(i->SITES[i]).toList());}
    static double mean(double[] a){return Arrays.stream(a).average().orElseThrow();}
    static int count(List<Data> data){int n=0;for(Data d:data)for(double v:d.nu)if(Double.isFinite(v))n++;return n;}

    void writeTransformed(String scenario,List<Data> data)throws IOException{
        for(Data d:data){double m=d.complete()?mean(d.nu):Double.NaN;double mx=d.complete()?Arrays.stream(d.nu).map(v->v*v).average().orElseThrow():Double.NaN;
            for(int s=0;s<4;s++){double v=d.nu[s];csv.get("transformed_series").row(scenario,d.t,SITES[s],v,v*v,v-m,v*v-mx,d.complete(),d.sourceRow);}}
    }

    static double[][] matrix(List<Data> data,boolean square){
        double[][] y=new double[data.size()][4];for(int i=0;i<data.size();i++)for(int s=0;s<4;s++){double v=data.get(i).nu[s];y[i][s]=square?v*v:v;}return y;
    }
    static Centered center(double[][] y,String mode){
        int n=y.length;double[] rm=new double[n],cm=new double[4];double grand=0;
        for(int i=0;i<n;i++)for(int s=0;s<4;s++){rm[i]+=y[i][s]/4;cm[s]+=y[i][s]/n;grand+=y[i][s]/(4*n);}
        double[][] a=new double[n][4],removed=new double[n][4];
        for(int i=0;i<n;i++)for(int s=0;s<4;s++){
            removed[i][s]=switch(mode){case "row"->rm[i];case "column"->cm[s];case "both"->rm[i]+cm[s]-grand;default->throw new IllegalArgumentException(mode);};
            a[i][s]=y[i][s]-removed[i][s];
        }return new Centered(a,removed);
    }
    void rankAnalyses(String scenario,List<Data> data)throws IOException{
        List<Data> complete=select(data,Data::complete);
        Map<String,List<Data>> subsets=new LinkedHashMap<>();
        subsets.put("all_complete",complete);
        subsets.put("through_160K",select(complete,d->d.t<=160));
        subsets.put("from_77K",select(complete,d->d.t>=77));
        for(var entry:subsets.entrySet())for(boolean square:new boolean[]{false,true})for(String centering:new String[]{"row","column","both"}){
            List<Data> rows=entry.getValue();double[][] y=matrix(rows,square);Centered centered=center(y,centering);String quantity=square?"nu_squared_MHz2":"nu_MHz";
            for(int rank=1;rank<=3;rank++){
                Numerics.RankFit fit=Numerics.rank(centered.matrix,rank);double frequencySse=0;
                for(int i=0;i<rows.size();i++)for(int s=0;s<4;s++){
                    double prediction=fit.fitted()[i][s]+centered.removed[i][s];
                    double freq=square?(prediction>=0?Math.sqrt(prediction):Double.NaN):prediction;
                    frequencySse+=(freq-rows.get(i).nu[s])*(freq-rows.get(i).nu[s]);
                    if(entry.getKey().equals("all_complete")&&rank<=2)csv.get("rank_reconstructions").row(scenario,quantity,centering,rank,rows.get(i).t,SITES[s],y[i][s],prediction,prediction-y[i][s],freq,rows.get(i).sourceRow);
                }
                csv.get("rank_summary").row(scenario,entry.getKey(),quantity,centering,rank,rows.size(),rows.size()*4,
                        fit.singular()[0],fit.singular()[1],fit.singular()[2],fit.singular()[3],fit.retainedFraction(),Math.sqrt(fit.sse()/(rows.size()*4)),fit.maximum(),Math.sqrt(frequencySse/(rows.size()*4)));
                if(rank==1&&entry.getKey().equals("all_complete"))for(int mode=0;mode<4;mode++)for(int s=0;s<4;s++)csv.get("rank_modes").row(scenario,quantity,centering,mode+1,SITES[s],fit.vectors()[s][mode],fit.singular()[mode]);
            }
        }
        for(Data omitted:complete){
            List<Data> remaining=select(complete,d->d!=omitted);int n=remaining.size()*4;
            double row=Math.sqrt(Numerics.rank(center(matrix(remaining,false),"row").matrix,1).sse()/n);
            double col=Math.sqrt(Numerics.rank(center(matrix(remaining,false),"column").matrix,1).sse()/n);
            double sq=Math.sqrt(Numerics.rank(center(matrix(remaining,true),"row").matrix,2).sse()/n);
            csv.get("rank_omit_one").row(scenario,omitted.t,row,col,sq);
        }
    }

    static double basis(String model,double t,double theta){
        if(model.equals("T_squared"))return (t/100)*(t/100);
        double z=theta/t;return z>700?0:1/Math.expm1(z);
    }
    static Numerics.Regression siteRegression(List<Data> data,int site,String model,double theta){
        List<Data> finite=select(data,d->Double.isFinite(d.nu[site]));double[] x=new double[finite.size()],y=new double[finite.size()];
        for(int i=0;i<finite.size();i++){x[i]=basis(model,finite.get(i).t,theta);y[i]=finite.get(i).nu[site];}
        return Numerics.regress(x,y);
    }
    static Thermal fitThermal(List<Data> train,String model){
        double[] theta=new double[4],a=new double[4],b=new double[4];boolean[] boundary=new boolean[4];
        Arrays.fill(theta,Double.NaN);
        if(model.equals("bose_common_theta")){
            Numerics.Minimum m=Numerics.gridGolden(logtheta->{double sse=0;for(int s=0;s<4;s++)sse+=siteRegression(train,s,model,Math.exp(logtheta)).sse();return sse;},LOG_THETA_LOW,LOG_THETA_HIGH,OPTIMIZER_GRID);
            Arrays.fill(theta,Math.exp(m.x()));Arrays.fill(boundary,m.boundary());
        }
        for(int s=0;s<4;s++){
            final int site=s;
            if(model.equals("bose_separate_theta")){
                Numerics.Minimum m=Numerics.gridGolden(logtheta->siteRegression(train,site,model,Math.exp(logtheta)).sse(),LOG_THETA_LOW,LOG_THETA_HIGH,OPTIMIZER_GRID);
                theta[s]=Math.exp(m.x());boundary[s]=m.boundary();
            }
            Numerics.Regression r=siteRegression(train,s,model,theta[s]);a[s]=r.intercept();b[s]=r.slope();
        }return new Thermal(model,theta,a,b,boundary);
    }
    static Score scoreThermal(List<Data> data,Thermal fit){
        int n=0;double sse=0,max=0;for(Data d:data)for(int s=0;s<4;s++)if(Double.isFinite(d.nu[s])){double r=fit.predict(d.t,s)-d.nu[s];sse+=r*r;max=Math.max(max,Math.abs(r));n++;}return new Score(n,sse,max);
    }
    void thermalAnalyses(String scenario,List<Data> data)throws IOException{
        for(String model:MODELS){
            double cvSse=0,cvMax=0;int cvN=0;
            for(int fold=-1;fold<BLOCKS.size();fold++){
                Block block=fold<0?null:BLOCKS.get(fold);String foldName=fold<0?"all":block.name;
                List<Data> train=fold<0?data:select(data,d->!block.contains(d.t));
                List<Data> test=fold<0?List.of():select(data,d->block.contains(d.t));
                Thermal fit=fitThermal(train,model);Score training=scoreThermal(train,fit),testing=scoreThermal(test,fit);
                int parameters=model.equals("bose_common_theta")?9:model.equals("bose_separate_theta")?12:8;
                csv.get("thermal_scores").row(scenario,model,foldName,training.n,testing.n,parameters,training.rms(),testing.rms(),training.max,testing.max);
                for(int s=0;s<4;s++){final int site=s;int n=select(train,d->Double.isFinite(d.nu[site])).size();csv.get("thermal_parameters").row(scenario,model,foldName,SITES[s],n,fit.theta[s],fit.a[s],fit.b[s],fit.boundary[s]);}
                for(Data d:data)for(int s=0;s<4;s++)if(Double.isFinite(d.nu[s])){
                    boolean heldOut=block!=null&&block.contains(d.t);double pred=fit.predict(d.t,s);
                    csv.get("thermal_predictions").row(scenario,model,foldName,heldOut?"test":"train",d.t,SITES[s],d.nu[s],pred,pred-d.nu[s],d.sourceRow);
                }
                if(fold>=0){cvSse+=testing.sse;cvN+=testing.n;cvMax=Math.max(cvMax,testing.max);}
            }
            csv.get("thermal_scores").row(scenario,model,"pooled_four_blocks",0,cvN,model.equals("bose_common_theta")?9:model.equals("bose_separate_theta")?12:8,Double.NaN,Math.sqrt(cvSse/cvN),Double.NaN,cvMax);
        }
    }

    static List<Harmonic> harmonics(List<Data> data,int[] order){
        List<Harmonic> out=new ArrayList<>();for(Data d:data){double[] x=new double[4];for(int m=0;m<4;m++)x[m]=d.nu[order[m]]*d.nu[order[m]];
            out.add(new Harmonic(d,mean(x),(x[0]-x[2])/2,(x[1]-x[3])/2,(x[0]-x[1]+x[2]-x[3])/4));}return out;
    }
    static double gramLower(double e,double p,int sign){
        double den1=4*(e+sign),den2=4*(e-sign);
        double l1=p==0?0:den1>0?p/den1:Double.POSITIVE_INFINITY;
        double l2=p==1?0:den2>0?(1-p)/den2:Double.POSITIVE_INFINITY;
        return Math.max(l1,l2);
    }
    static TensorFit fitTensor(List<Harmonic> h){
        int n=h.size();double[][] bc=new double[n][2];double[] z=new double[n],d=new double[n],y=new double[n];
        int sign=(int)Math.signum(h.get(0).d);boolean same=sign!=0;
        for(int i=0;i<n;i++){
            Harmonic r=h.get(i);bc[i][0]=r.b;bc[i][1]=r.c;z[i]=r.power()/r.d;d[i]=r.d;y[i]=r.mean;
            if(Math.signum(r.d)!=sign||Math.abs(r.d)<1e-12)same=false;
        }
        Numerics.RankFit first=Numerics.rank(bc,1);double ux=first.vectors()[0][0],uy=first.vectors()[1][0];
        if(ux<0){ux=-ux;uy=-uy;}
        Numerics.LeastSquares identity=Numerics.leastSquares2(z,d,y);
        double alpha=identity.coefficients()[0],beta=identity.coefficients()[1],identitySse=0;
        for(int i=0;i<n;i++){double r=alpha*z[i]+beta*d[i]-y[i];identitySse+=r*r;}
        double p=ux*ux;
        boolean gram=same&&Math.signum(alpha)==sign&&Math.signum(beta)==sign&&Math.abs(beta)>=1-1e-12&&alpha*(beta+1)>=p/4-1e-12&&alpha*(beta-1)>=(1-p)/4-1e-12;
        double physicalA=Double.NaN,physicalE=Double.NaN,objective=Double.NaN;boolean boundary=false;
        if(same){
            double[] zz=new double[n],dd=new double[n];
            for(int i=0;i<n;i++){double projection=ux*h.get(i).b+uy*h.get(i).c;dd[i]=Math.abs(d[i]);zz[i]=projection*projection/dd[i];}
            final int dsign=sign;final double pp=p;
            Numerics.Minimum minimum=Numerics.gridGolden(loge->{double e=Math.exp(loge),a=Math.max(gramLower(e,pp,dsign),(Numerics.dot(zz,y)-e*Numerics.dot(zz,dd))/Numerics.sumSquares(zz));double sse=0;
                for(int i=0;i<n;i++){double residual=a*zz[i]+e*dd[i]-y[i];sse+=residual*residual;}return sse;
            },Math.log(1+1e-10),Math.log(1e6),OPTIMIZER_GRID);
            physicalE=Math.exp(minimum.x());physicalA=Math.max(gramLower(physicalE,p,sign),(Numerics.dot(zz,y)-physicalE*Numerics.dot(zz,dd))/Numerics.sumSquares(zz));
            boundary=minimum.boundary();objective=minimum.value();
        }
        return new TensorFit(ux,uy,first.singular(),Math.sqrt(first.sse()/(2*n)),alpha,beta,identity.condition(),Math.sqrt(identitySse/n),sign,same,gram,physicalA,physicalE,boundary,objective);
    }
    static double[] reconstruct(double mean,double b,double c,double d){return new double[]{mean+d+b,mean-d+c,mean+d-b,mean-d-c};}
    static TensorScore scoreTensor(List<Harmonic> h,int[] order,TensorFit fit){
        double firstSse=0,firstMax=0,identitySse=0,physicalSse=0,physicalMax=0;int mismatch=0;
        for(Harmonic row:h){
            double score=fit.ux*row.b+fit.uy*row.c,bhat=fit.ux*score,chat=fit.uy*score;
            double[] first=reconstruct(row.mean,bhat,chat,row.d);
            double meanhat=Math.signum(row.d)==fit.dSign?fit.physicalA*score*score/Math.abs(row.d)+fit.physicalE*Math.abs(row.d):Double.NaN;
            double[] physical=reconstruct(meanhat,bhat,chat,row.d);
            double identity=fit.alpha*row.power()/row.d+fit.beta*row.d-row.mean;identitySse+=identity*identity;
            if(Math.signum(row.d)!=fit.dSign)mismatch++;
            for(int m=0;m<4;m++){
                double f=first[m]>=0?Math.sqrt(first[m]):Double.NaN;
                double p=physical[m]>=0?Math.sqrt(physical[m]):Double.NaN;
                double rf=f-row.data.nu[order[m]],rp=p-row.data.nu[order[m]];
                firstSse+=rf*rf;physicalSse+=rp*rp;firstMax=Math.max(firstMax,Math.abs(rf));physicalMax=Math.max(physicalMax,Math.abs(rp));
            }
        }
        return new TensorScore(h.size(),firstSse,firstMax,identitySse,physicalSse,physicalMax,mismatch);
    }
    void tensorAnalyses(String scenario,List<Data> data)throws IOException{
        List<Data> complete=select(data,Data::complete);
        for(int oi=0;oi<ORDERS.length;oi++){
            int[] order=ORDERS[oi];String name=ORDER_NAMES[oi],sites=orderText(order);
            List<Harmonic> all=harmonics(complete,order);
            for(Harmonic r:all)csv.get("tensor_coefficients").row(scenario,name,sites,r.data.t,r.mean,r.b,r.c,r.d,r.c/r.b,r.power(),r.data.sourceRow);
            Map<String,List<Harmonic>> subsets=new LinkedHashMap<>();subsets.put("all_complete",all);
            subsets.put("through_160K",all.stream().filter(r->r.data.t<=160).toList());
            subsets.put("from_77K",all.stream().filter(r->r.data.t>=77).toList());
            for(var entry:subsets.entrySet()){
                List<Harmonic> rows=entry.getValue();TensorFit fit=fitTensor(rows);TensorScore score=scoreTensor(rows,order,fit);
                double ratioMin=Double.POSITIVE_INFINITY,ratioMax=Double.NEGATIVE_INFINITY,amin=Double.POSITIVE_INFINITY,amax=Double.NEGATIVE_INFINITY,dmin=Double.POSITIVE_INFINITY,dmax=Double.NEGATIVE_INFINITY;
                double angle=Math.atan2(fit.uy,fit.ux);
                for(Harmonic r:rows){double ratio=r.c/r.b;ratioMin=Math.min(ratioMin,ratio);ratioMax=Math.max(ratioMax,ratio);double da=Math.atan2(r.c,r.b)-angle;while(da>Math.PI/2)da-=Math.PI;while(da< -Math.PI/2)da+=Math.PI;amin=Math.min(amin,da);amax=Math.max(amax,da);dmin=Math.min(dmin,r.d);dmax=Math.max(dmax,r.d);}
                double p=fit.ux*fit.ux,A=fit.physicalA,E=fit.physicalE;
                double marginCos=A*(E+fit.dSign)-p/4,marginSin=A*(E-fit.dSign)-(1-p)/4;
                csv.get("tensor_summary").row(scenario,name,sites,entry.getKey(),rows.size(),Math.toDegrees(angle),p,fit.singular[0],fit.singular[1],fit.firstRmsMHz2,score.firstRms(),ratioMin,ratioMax,Math.toDegrees(amax-amin),dmin,dmax,fit.dSignConstant,fit.alpha,fit.beta,fit.identityCondition,fit.identityRmsMHz2,fit.rawGramFeasible,A,E,fit.physicalBoundary,score.conditionalRms(),score.conditionalMaxMHz,marginCos,marginSin,Math.min(marginCos,marginSin)<1e-8);
                if(entry.getKey().equals("all_complete"))writeTensorPredictions(scenario,name,"all_complete","train",rows,order,fit);
            }
            for(Block block:BLOCKS){
                List<Harmonic> train=all.stream().filter(r->!block.contains(r.data.t)).toList(),test=all.stream().filter(r->block.contains(r.data.t)).toList();
                if(test.isEmpty())continue;
                TensorFit fit=fitTensor(train);TensorScore score=scoreTensor(test,order,fit);
                csv.get("tensor_block_validation").row(scenario,name,block.name,train.size(),test.size(),Math.toDegrees(Math.atan2(fit.uy,fit.ux)),fit.alpha,fit.beta,fit.physicalA,fit.physicalE,score.firstRms(),score.identityRms(),score.conditionalRms(),score.signMismatch);
                writeTensorPredictions(scenario,name,"held_out_"+block.name,"test",test,order,fit);
            }
            for(Harmonic omitted:all){
                List<Harmonic> train=all.stream().filter(r->r!=omitted).toList();TensorFit fit=fitTensor(train);
                TensorScore training=scoreTensor(train,order,fit),testing=scoreTensor(List.of(omitted),order,fit);
                csv.get("tensor_omit_one").row(scenario,name,omitted.data.t,Math.toDegrees(Math.atan2(fit.uy,fit.ux)),fit.alpha,fit.beta,fit.physicalA,fit.physicalE,training.firstRms(),training.identityRms(),training.conditionalRms(),testing.firstRms(),testing.identityRms(),testing.conditionalRms());
            }
        }
    }
    void writeTensorPredictions(String scenario,String name,String subset,String partition,List<Harmonic> rows,int[] order,TensorFit fit)throws IOException{
        for(Harmonic row:rows){
            double z=fit.ux*row.b+fit.uy*row.c,bhat=fit.ux*z,chat=fit.uy*z;
            double predictedMean=Math.signum(row.d)==fit.dSign?fit.physicalA*z*z/Math.abs(row.d)+fit.physicalE*Math.abs(row.d):Double.NaN;
            double[] first=reconstruct(row.mean,bhat,chat,row.d),physical=reconstruct(predictedMean,bhat,chat,row.d);
            double psi=Math.sqrt(Math.abs(row.d)),v=Math.sqrt(fit.physicalA)*z/psi;
            for(int m=0;m<4;m++){
                double f=first[m]>=0?Math.sqrt(first[m]):Double.NaN,p=physical[m]>=0?Math.sqrt(physical[m]):Double.NaN;
                csv.get("tensor_predictions").row(scenario,name,subset,partition,row.data.t,SITES[order[m]],row.data.nu[order[m]],f,p,p-row.data.nu[order[m]],predictedMean,row.d,v,psi,row.data.sourceRow);
            }
        }
    }

    void check(String name,double value,double expected,double tolerance,String detail)throws IOException{
        boolean pass=Double.isFinite(value)&&Math.abs(value-expected)<=tolerance;
        csv.get("validation_checks").row(name,value,expected+" +/- "+tolerance,pass,detail);if(!pass)failed.add(name);
    }
    void validate(List<Data> data)throws Exception{
        String hash=sha(input);boolean known=hash.equals(EXPECTED_INPUT_SHA);
        csv.get("validation_checks").row("input_sha256",hash,EXPECTED_INPUT_SHA,known,"Authoritative archived adapter, no interpolation");
        if(!known)throw new IllegalArgumentException("Input SHA-256 differs; use the archived data for this reproducible analysis");
        check("n_temperatures",data.size(),24,0,"Finite measured temperature rows");
        check("n_observations",count(data),89,0,"Missing values stay missing");
        check("n_complete_quartets",select(data,Data::complete).size(),18,0,"No inferred quartet values");
        check("label_scenario_preserves_count",count(scenario(data,SCENARIOS[1])),89,0,"Swap includes missing cells, not frequency sorting");
        Path source=input.resolveSibling("temperature_source_Data1.csv");
        if(Files.isRegularFile(source)){
            String rawsha=sha(source);boolean pass=rawsha.equals(EXPECTED_SOURCE_SHA);csv.get("validation_checks").row("original_Data1_sha256",rawsha,EXPECTED_SOURCE_SHA,pass,"Original source bytes including line endings");if(!pass)failed.add("source_sha256");
        }
        double[][] knownSvd={{3,0},{0,4},{0,0}};Numerics.RankFit rank=Numerics.rank(knownSvd,1);
        check("svd_known_largest",rank.singular()[0],4,1e-12,"Independent diagonal-matrix singular value");
        check("svd_known_rank1_sse",rank.sse(),9,1e-12,"Eckart-Young residual");
        Numerics.Minimum min=Numerics.gridGolden(x->(x-1.23456789)*(x-1.23456789)+2,-3,4,80);
        check("bounded_optimizer_quadratic",min.x(),1.23456789,2e-7,"Known interior minimum");
        List<Data> full=select(data,Data::complete);
        double originalColumn=Math.sqrt(Numerics.rank(center(matrix(full,false),"column").matrix,1).sse()/(full.size()*4));
        check("independent_pilot_nu_column_rank1",originalColumn,0.017603359696195503,2e-10,"Previously independent NumPy pilot, not code-generated expectation");
        double originalRow=Math.sqrt(Numerics.rank(center(matrix(full,false),"row").matrix,1).sse()/(full.size()*4));
        check("independent_pilot_nu_row_rank1",originalRow,0.1430771846174922,2e-10,"Previously independent NumPy pilot");
        Thermal thermal=fitThermal(data,"bose_common_theta");
        check("independent_pilot_common_Bose_theta",thermal.theta[0],188.66454957462528,0.002,"Independent nonlinear fit with 89 observations");
        check("independent_pilot_common_Bose_rms",scoreThermal(data,thermal).rms(),0.024880116398679177,2e-10,"Independent nonlinear fit");
        // Synthetic data from explicit independent tensor-vector components in 5D.
        // U=(10,0,0,0,0), C=(0.8,1.2,0.2,0,0), S=(-0.5,0,0.3,0.4,0).
        double[] u={10,0,0,0,0},c={0.8,1.2,0.2,0,0},s={-0.5,0,0.3,0.4,0};
        List<Data> synthetic=new ArrayList<>();
        for(int j=0;j<12;j++){
            double v=1+0.015*j,psi=0.2+0.025*j;double[] nu=new double[4];
            for(int m=0;m<4;m++){double q=0;for(int k=0;k<5;k++){double component=v*u[k]+psi*(c[k]*Math.cos(Math.PI*m/2)+s[k]*Math.sin(Math.PI*m/2));q+=component*component;}nu[m]=Math.sqrt(q);}
            synthetic.add(new Data(10+10*j,nu,j+1));
        }
        int[] identityOrder={0,1,2,3};List<Harmonic> h=harmonics(synthetic,identityOrder);TensorFit fit=fitTensor(h);TensorScore score=scoreTensor(h,identityOrder,fit);
        check("explicit_tensor_first_harmonic_rms",fit.firstRmsMHz2,0,2e-12,"Generated by vector tensors, not by the fitted harmonic identity");
        check("explicit_tensor_identity_rms",fit.identityRmsMHz2,0,1e-8,"Two-amplitude tensor construction must satisfy second identity");
        check("explicit_tensor_conditional_frequency_rms",score.conditionalRms(),0,2e-7,"Constrained recovery of physically feasible synthetic tensor sequence");
        check("explicit_tensor_Gram_feasible",fit.rawGramFeasible?1:0,1,0,"The true Gram matrix is positive semidefinite");
        // A rotating first-harmonic direction is not a fixed C,S tensor pattern.
        List<Harmonic> rotating=new ArrayList<>();for(int j=0;j<12;j++){double angle=0.05*j,amplitude=3+0.1*j;rotating.add(new Harmonic(synthetic.get(j),100,amplitude*Math.cos(angle),amplitude*Math.sin(angle),1));}
        TensorFit rejected=fitTensor(rotating);boolean detect=rejected.firstRmsMHz2>0.1;
        csv.get("validation_checks").row("rotating_direction_counterexample",rejected.firstRmsMHz2,"> 0.1 MHz2",detect,"Necessary first-harmonic test must detect a changing direction");if(!detect)failed.add("rotating_direction_counterexample");
    }
    static String sha(Path file)throws Exception{return HexFormat.of().formatHex(MessageDigest.getInstance("SHA-256").digest(Files.readAllBytes(file)));}
    void writeManifest(List<Data> data)throws Exception{
        String text="{\n"+
                "  \"module\": \"joint-temperature\",\n"+
                "  \"source_language\": \"Java 17\",\n"+
                "  \"runtime\": \""+System.getProperty("java.runtime.version")+"\",\n"+
                "  \"input_sha256\": \""+sha(input)+"\",\n"+
                "  \"source_Data1_sha256\": \""+EXPECTED_SOURCE_SHA+"\",\n"+
                "  \"source_commit\": \""+SOURCE_COMMIT+"\",\n"+
                "  \"n_temperatures\": "+data.size()+",\n"+
                "  \"n_finite_observations\": "+count(data)+",\n"+
                "  \"n_complete_quartets\": "+select(data,Data::complete).size()+",\n"+
                "  \"bose_theta_bounds_K\": [1, 3000],\n"+
                "  \"physical_E_bounds\": [1.0000000001, 1000000],\n"+
                "  \"scalar_optimizer_grid_intervals\": "+OPTIMIZER_GRID+",\n"+
                "  \"label_scenarios\": [\"original\", \"swap_s3_s4_T_ge_230K\"],\n"+
                "  \"period4_global_classes\": 3,\n"+
                "  \"experimental_sigma_supplied\": false,\n"+
                "  \"formal_significance_calculated\": false,\n"+
                "  \"random_sampling\": false,\n"+
                "  \"physical_construction\": \"conditional TLS direction plus Gram-constrained identity; observed D retained; not global tensor least squares\"\n"+
                "}\n";
        Files.writeString(out.resolve("calculation_manifest.json"),text,StandardCharsets.UTF_8);
    }
}
