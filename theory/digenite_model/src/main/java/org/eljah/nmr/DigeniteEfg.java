package org.eljah.nmr;

import java.io.*;
import java.nio.charset.StandardCharsets;
import java.nio.file.*;
import java.util.*;
import java.util.zip.GZIPOutputStream;

/**
 * Periodic point-charge electric-field gradients for an explicitly specified
 * rhombohedral digenite model. No external dependencies. Java 17 or newer.
 *
 * Charges are dimensionless z=q/e, lengths in Angstrom. Raw Hessians G have
 * units Angstrom^-3; the elementary charge and Coulomb constant are factored out.
 * V = electrostaticFactor * G in V/m^2. No Sternheimer correction is assumed.
 * Inherited and fitted scales are kept separate from bare predictions.
 */
public final class DigeniteEfg {
    private static final double SQRT_PI = Math.sqrt(Math.PI);
    private static final double TWO_PI = 2 * Math.PI;
    private static final double ELEMENTARY_CHARGE = 1.602176634e-19;
    private static final double PLANCK = 6.62607015e-34;
    private static final String[] ORBITS = {"DgA", "DgB", "DgC", "DgD", "DgE"};
    private static final int[] MULTIPLICITIES = {1, 2, 2, 2, 2};
    private static final String[] COMPONENTS = {"xx", "yy", "zz", "xy", "xz", "yz"};

    private record Vec(double x, double y, double z) {
        Vec add(Vec b) { return new Vec(x+b.x, y+b.y, z+b.z); }
        Vec sub(Vec b) { return new Vec(x-b.x, y-b.y, z-b.z); }
        Vec mul(double c) { return new Vec(c*x, c*y, c*z); }
        double dot(Vec b) { return x*b.x+y*b.y+z*b.z; }
        double norm2() { return dot(this); }
        double norm() { return Math.sqrt(norm2()); }
        Vec cross(Vec b) { return new Vec(y*b.z-z*b.y, z*b.x-x*b.z, x*b.y-y*b.x); }
        double at(int i) { return i == 0 ? x : i == 1 ? y : z; }
    }

    private static final class Cell {
        final Vec[] a, inverseRows, reciprocal;
        final double volume;
        Cell(Vec a1, Vec a2, Vec a3) {
            a = new Vec[]{a1,a2,a3};
            volume = a1.dot(a2.cross(a3));
            if (!(volume > 0)) throw new IllegalArgumentException("Cell volume must be positive");
            inverseRows = new Vec[]{a2.cross(a3).mul(1/volume),
                    a3.cross(a1).mul(1/volume), a1.cross(a2).mul(1/volume)};
            reciprocal = Arrays.stream(inverseRows).map(v -> v.mul(TWO_PI)).toArray(Vec[]::new);
        }
        static Cell rhombohedral(double length, double angleDeg) {
            double c = Math.cos(Math.toRadians(angleDeg));
            double h = length * Math.sqrt((1+2*c)/3);
            double rho = length * Math.sqrt(2*(1-c)/3);
            return new Cell(new Vec(rho,0,h), new Vec(-rho/2,Math.sqrt(3)*rho/2,h),
                    new Vec(-rho/2,-Math.sqrt(3)*rho/2,h));
        }
        Vec cart(Vec f) { return a[0].mul(f.x).add(a[1].mul(f.y)).add(a[2].mul(f.z)); }
        Vec fractional(Vec r) { return new Vec(inverseRows[0].dot(r),inverseRows[1].dot(r),inverseRows[2].dot(r)); }
        Cell scale(double s) { return new Cell(a[0].mul(s),a[1].mul(s),a[2].mul(s)); }
    }

    private record Atom(String name, String orbit, String element, Vec frac, int partner) {}
    private record Model(Cell cell, List<Atom> atoms, int[] targets) {
        Model scale(double s) { return new Model(cell.scale(s),atoms,targets); }
        Model displacePair(int orbitIndex, Vec positiveDisplacement) {
            if (orbitIndex == 0) throw new IllegalArgumentException("DgA is an inversion-centre site and must not be displaced by this operation");
            Vec df = cell.fractional(positiveDisplacement);
            List<Atom> moved = new ArrayList<>();
            for (Atom at: atoms) moved.add(at.orbit.equals(ORBITS[orbitIndex])
                    ? new Atom(at.name,at.orbit,at.element,at.frac.add(df.mul(at.partner)),at.partner) : at);
            return new Model(cell,List.copyOf(moved),targets);
        }
    }

    private static final class Tensor {
        final double[] v = new double[6];
        Tensor() {}
        Tensor(double[] a) { System.arraycopy(a,0,v,0,6); }
        Tensor add(Tensor b) { Tensor t = new Tensor(); for(int i=0;i<6;i++) t.v[i]=v[i]+b.v[i]; return t; }
        Tensor scale(double s) { Tensor t=new Tensor(); for(int i=0;i<6;i++) t.v[i]=s*v[i]; return t; }
        double trace() { return v[0]+v[1]+v[2]; }
        double maxDifference(Tensor b) { double d=0; for(int i=0;i<6;i++) d=Math.max(d,Math.abs(v[i]-b.v[i])); return d; }
        double[][] matrix() { return new double[][]{{v[0],v[3],v[4]},{v[3],v[1],v[5]},{v[4],v[5],v[2]}}; }
    }

    /** Kahan sums avoid systematic accumulation error in long real/reciprocal sums. */
    private static final class TensorSum {
        final double[] s=new double[6], c=new double[6];
        void add(double[] values) {
            for (int i=0;i<6;i++) { double y=values[i]-c[i], t=s[i]+y; c[i]=(t-s[i])-y; s[i]=t; }
        }
        Tensor tensor() { return new Tensor(s); }
    }

    private record Eigen(double xx, double yy, double zz, double eta, Vec zAxis) {
        double p() { return Math.abs(zz)*Math.sqrt(1+eta*eta/3); }
    }
    private record Result(Tensor real, Tensor reciprocal, Tensor self, Tensor total,
                          Eigen eigen, long realTerms, long reciprocalVectors) {}
    private record Scenario(String name, double[] charges) {}
    private record Assignment(String scenario, String[] mapping, int omitted,
                              double[] bare, double[] inherited, double[] fitted,
                              double bareScale, double inheritedScale, double omittedBare,
                              double omittedInherited, double omittedFitted, double rmse, double maxError) {}

    private static final class Config {
        final Properties p = new Properties();
        Config(String[] args) throws IOException {
            String[][] defaults = {
                    {"a","16.16"},{"angleDeg","13.933333333333334"},
                    {"x.DgA","0.5"},{"x.DgB","0.06"},{"x.DgC","0.133"},{"x.DgD","0.25"},{"x.DgE","0.35"},
                    {"ewaldAlphas","0.15,0.25,0.4"},{"mainAlpha","0.25"},
                    {"realFactor","6"},{"reciprocalFactor","12"},
                    {"q63Barn","-0.211"},{"q65Barn","-0.195"},
                    {"electrostaticFactor","1.43996454687e21"},{"inheritedScaleMHz","38.2916669626624"},
                    {"observedMHz","11.97,15.27,16.18,16.96"},{"output","results"},
                    {"rawTerms","true"},{"sensitivities","true"},{"directSums","true"},
                    {"directRadii","20,40,60,80,120"},{"uniformLambdas","0.99,0.995,1,1.005,1.01"},
                    {"angleAlternativeDeg","13.93"},{"displacementAngstrom","0.01"},
                    {"transverseExample","true"},{"verify","true"},{"convergenceTol","1e-9"}};
            for(String[] row: defaults) p.setProperty(row[0],row[1]);
            for(String arg:args) if(arg.startsWith("--config=")) {
                Properties read=new Properties();
                try(Reader r=Files.newBufferedReader(Path.of(arg.substring(9)),StandardCharsets.UTF_8)) { read.load(r); }
                for(String key:read.stringPropertyNames()) {
                    if(!p.containsKey(key)) throw new IllegalArgumentException("Unknown property: "+key);
                    p.setProperty(key,read.getProperty(key));
                }
            }
            for(String arg:args) {
                if(arg.startsWith("--config=")) continue;
                if(!arg.startsWith("--") || !arg.contains("=")) throw new IllegalArgumentException("Expected --key=value: "+arg);
                int eq=arg.indexOf('='); String key=arg.substring(2,eq);
                if(!p.containsKey(key)) throw new IllegalArgumentException("Unknown option: "+key);
                p.setProperty(key,arg.substring(eq+1));
            }
            if(array("observedMHz").length!=4) throw new IllegalArgumentException("Exactly four observed frequencies are required");
            if(!(d("a")>0) || !(d("angleDeg")>0 && d("angleDeg")<120)) throw new IllegalArgumentException("Invalid rhombohedral geometry");
            if(d("realFactor")<=0 || d("reciprocalFactor")<=0 || d("mainAlpha")<=0) throw new IllegalArgumentException("Ewald cutoffs and alpha must be positive");
            for(double a:array("ewaldAlphas")) if(a<=0) throw new IllegalArgumentException("Ewald alpha must be positive");
        }
        double d(String key) { return Double.parseDouble(p.getProperty(key)); }
        boolean b(String key) { return Boolean.parseBoolean(p.getProperty(key)); }
        double[] array(String key) { return Arrays.stream(p.getProperty(key).split(",")).map(String::trim).mapToDouble(Double::parseDouble).toArray(); }
        Path output() { return Path.of(p.getProperty("output")); }
    }

    private static final class Csv implements AutoCloseable {
        final BufferedWriter out;
        Csv(Path file, boolean gzip) throws IOException {
            OutputStream stream=Files.newOutputStream(file);
            if(gzip) stream=new GZIPOutputStream(stream,65536);
            out=new BufferedWriter(new OutputStreamWriter(stream,StandardCharsets.UTF_8),65536);
        }
        void row(Object... values) throws IOException {
            for(int i=0;i<values.length;i++) {
                if(i>0) out.write(',');
                String s=String.valueOf(values[i]);
                if(s.indexOf(',')>=0 || s.indexOf('"')>=0 || s.indexOf('\n')>=0) {
                    out.write('"'); out.write(s.replace("\"","\"\"")); out.write('"');
                } else out.write(s);
            }
            out.newLine();
        }
        void row(List<?> values) throws IOException { row(values.toArray()); }
        @Override public void close() throws IOException { out.close(); }
    }

    private static Model model(Config cfg, double angle) {
        Cell c=Cell.rhombohedral(cfg.d("a"),angle);
        List<Atom> atoms=new ArrayList<>(); int[] targets=new int[5];
        for(int o=0;o<5;o++) {
            double x=cfg.d("x."+ORBITS[o]); targets[o]=atoms.size();
            atoms.add(new Atom(ORBITS[o]+"+",ORBITS[o],"Cu",new Vec(x,x,x),1));
            if(MULTIPLICITIES[o]==2) atoms.add(new Atom(ORBITS[o]+"-",ORBITS[o],"Cu",new Vec(-x,-x,-x),-1));
        }
        atoms.add(new Atom("S0","S0","S",new Vec(0,0,0),0));
        for(int o=1;o<=2;o++) for(int sign:new int[]{1,-1}) {
            double x=sign*o*0.2; atoms.add(new Atom("S"+o+(sign>0?"+":"-"),"S"+o,"S",new Vec(x,x,x),sign));
        }
        return new Model(c,List.copyOf(atoms),targets);
    }

    private static List<Scenario> scenarios(Model model) {
        List<Scenario> scenarios=new ArrayList<>();
        for(int extra=-1;extra<5;extra++) {
            double[] q=new double[model.atoms.size()];
            for(int j=0;j<q.length;j++) {
                Atom at=model.atoms.get(j);
                q[j]=at.element.equals("S") ? -2 : extra<0 ? 10.0/9 : 1;
                if(extra>=0 && at.orbit.equals(ORBITS[extra])) q[j]+=1.0/MULTIPLICITIES[extra];
            }
            double net=Arrays.stream(q).sum();
            if(Math.abs(net)>1e-12) throw new IllegalStateException("Non-neutral charges: "+net);
            scenarios.add(new Scenario(extra<0?"uniform_Cu10over9":"extra_charge_"+ORBITS[extra],q));
        }
        return List.copyOf(scenarios);
    }

    /** erf Taylor series for x<=1.5; regularized incomplete gamma Q continued fraction above it. */
    static double erfc(double x) {
        if(x<0) return 2-erfc(-x);
        if(x==0) return 1;
        if(x<=1.5) {
            double term=x, sum=x, compensation=0;
            for(int n=1;n<300;n++) {
                term *= -x*x/n;
                double add=term/(2*n+1), y=add-compensation, t=sum+y;
                compensation=(t-sum)-y; sum=t;
                if(Math.abs(add)<2e-17*Math.abs(sum)) return 1-2*sum/SQRT_PI;
            }
            throw new ArithmeticException("erf Taylor series did not converge");
        }
        double t=x*x;
        if(t>745) return 0;
        double b=t+0.5, c=1e300, d=1/b, h=d;
        for(int n=1;n<10000;n++) {
            double an=-n*(n-0.5); b+=2;
            d=an*d+b; if(Math.abs(d)<1e-300) d=Math.copySign(1e-300,d);
            c=b+an/c; if(Math.abs(c)<1e-300) c=Math.copySign(1e-300,c);
            d=1/d; double delta=d*c; h*=delta;
            if(Math.abs(delta-1)<3e-16) return Math.exp(-t+0.5*Math.log(t)-0.5*Math.log(Math.PI))*h;
        }
        throw new ArithmeticException("Gamma Q continued fraction did not converge");
    }

    private static double[] hessian(Vec r, double a, double b, double q) {
        return new double[]{q*(a*r.x*r.x+b),q*(a*r.y*r.y+b),q*(a*r.z*r.z+b),
                q*a*r.x*r.y,q*a*r.x*r.z,q*a*r.y*r.z};
    }

    private static Result ewald(Model m, Scenario scenario, int targetIndex, double alpha,
                                Config cfg, Csv rawReal, Csv rawReciprocal) throws IOException {
        Atom target=m.atoms.get(targetIndex); Cell c=m.cell;
        double rcut=cfg.d("realFactor")/alpha, kcut=cfg.d("reciprocalFactor")*alpha;
        double rcut2=rcut*rcut, kcut2=kcut*kcut;
        TensorSum real=new TensorSum(), reciprocal=new TensorSum(); long nr=0, nk=0;
        for(int j=0;j<m.atoms.size();j++) {
            Atom source=m.atoms.get(j); Vec delta=source.frac.sub(target.frac);
            int[] lo=new int[3],hi=new int[3];
            for(int i=0;i<3;i++) {
                double bound=c.inverseRows[i].norm()*rcut;
                lo[i]=(int)Math.ceil(-delta.at(i)-bound-1e-12);
                hi[i]=(int)Math.floor(-delta.at(i)+bound+1e-12);
            }
            for(int n1=lo[0];n1<=hi[0];n1++) for(int n2=lo[1];n2<=hi[1];n2++) for(int n3=lo[2];n3<=hi[2];n3++) {
                Vec r=c.cart(delta.add(new Vec(n1,n2,n3))); double r2=r.norm2();
                if(r2<1e-24) {
                    if(j!=targetIndex || n1!=0 || n2!=0 || n3!=0) throw new IllegalArgumentException("Overlapping sites in the input model");
                    continue;
                }
                if(r2>rcut2) continue;
                double radius=Math.sqrt(r2), ex=Math.exp(-alpha*alpha*r2), ec=erfc(alpha*radius);
                double a=3*ec/(r2*r2*radius)+6*alpha*ex/(SQRT_PI*r2*r2)+4*alpha*alpha*alpha*ex/(SQRT_PI*r2);
                double b=-ec/(r2*radius)-2*alpha*ex/(SQRT_PI*r2);
                double[] term=hessian(r,a,b,scenario.charges[j]); real.add(term); nr++;
                if(rawReal!=null) {
                    List<Object> row=new ArrayList<>(List.of(target.name,j,source.name,source.orbit,scenario.charges[j],n1,n2,n3,r.x,r.y,r.z,radius));
                    append(row,term); rawReal.row(row);
                }
            }
        }
        int[] bounds=new int[3]; for(int i=0;i<3;i++) bounds[i]=(int)Math.ceil(c.a[i].norm()*kcut/TWO_PI);
        for(int n1=-bounds[0];n1<=bounds[0];n1++) for(int n2=-bounds[1];n2<=bounds[1];n2++) for(int n3=-bounds[2];n3<=bounds[2];n3++) {
            if(n1==0 && n2==0 && n3==0) continue;
            Vec k=c.reciprocal[0].mul(n1).add(c.reciprocal[1].mul(n2)).add(c.reciprocal[2].mul(n3));
            double k2=k.norm2(); if(k2>kcut2) continue; nk++;
            double prefactor=-4*Math.PI/c.volume*Math.exp(-k2/(4*alpha*alpha))/k2;
            for(int j=0;j<m.atoms.size();j++) {
                Atom source=m.atoms.get(j); Vec df=target.frac.sub(source.frac);
                double phase=TWO_PI*(n1*df.x+n2*df.y+n3*df.z), qcos=scenario.charges[j]*Math.cos(phase);
                double[] term=hessian(k,prefactor,0,qcos); reciprocal.add(term);
                if(rawReciprocal!=null) {
                    List<Object> row=new ArrayList<>(List.of(target.name,n1,n2,n3,j,source.name,source.orbit,scenario.charges[j],k.x,k.y,k.z,Math.sqrt(k2),phase,qcos));
                    append(row,term); rawReciprocal.row(row);
                }
            }
        }
        double s=4*alpha*alpha*alpha*scenario.charges[targetIndex]/(3*SQRT_PI);
        Tensor self=new Tensor(new double[]{s,s,s,0,0,0});
        Tensor rt=real.tensor(), kt=reciprocal.tensor(), total=rt.add(kt).add(self);
        return new Result(rt,kt,self,total,eigen(total),nr,nk);
    }

    /** Jacobi rotations for a real symmetric 3x3 matrix; eigenvalues retain sign. */
    private static Eigen eigen(Tensor t) {
        double[][] a=t.matrix(), v={{1,0,0},{0,1,0},{0,0,1}};
        for(int iter=0;iter<100;iter++) {
            int p=0,q=1; double largest=Math.abs(a[0][1]);
            for(int i=0;i<3;i++) for(int j=i+1;j<3;j++) if(Math.abs(a[i][j])>largest) {p=i;q=j;largest=Math.abs(a[i][j]);}
            double norm=Math.max(Math.abs(a[0][0]),Math.max(Math.abs(a[1][1]),Math.abs(a[2][2])));
            if(largest<=1e-16*Math.max(norm,1e-12)) break;
            double tau=(a[q][q]-a[p][p])/(2*a[p][q]);
            double rot=Math.copySign(1.0,tau)/(Math.abs(tau)+Math.sqrt(1+tau*tau));
            double cs=1/Math.sqrt(1+rot*rot), sn=rot*cs, apq=a[p][q];
            a[p][p]-=rot*apq; a[q][q]+=rot*apq; a[p][q]=a[q][p]=0;
            for(int j=0;j<3;j++) if(j!=p && j!=q) {
                double jp=a[j][p],jq=a[j][q];
                a[j][p]=a[p][j]=cs*jp-sn*jq; a[j][q]=a[q][j]=sn*jp+cs*jq;
            }
            for(int j=0;j<3;j++) {double jp=v[j][p],jq=v[j][q];v[j][p]=cs*jp-sn*jq;v[j][q]=sn*jp+cs*jq;}
        }
        Integer[] order={0,1,2}; Arrays.sort(order,Comparator.comparingDouble(i->Math.abs(a[i][i])));
        double xx=a[order[0]][order[0]], yy=a[order[1]][order[1]], zz=a[order[2]][order[2]];
        double eta=Math.abs(zz)<1e-25?0:(xx-yy)/zz;
        if(Math.abs(eta)<1e-12) eta=0;
        return new Eigen(xx,yy,zz,eta,new Vec(v[0][order[2]],v[1][order[2]],v[2][order[2]]));
    }

    private static double frequency(Eigen e, double qBarn, Config cfg) {
        return ELEMENTARY_CHARGE*Math.abs(qBarn)*1e-28*cfg.d("electrostaticFactor")*e.p()/(2*PLANCK)/1e6;
    }
    private static double inherited(Eigen e, Config cfg) { return cfg.d("inheritedScaleMHz")*e.p(); }
    private static void append(List<Object> row,double[] a) { for(double v:a) row.add(v); }
    private static void appendTensorHeader(List<Object> row,String prefix) { for(String part:COMPONENTS) row.add(prefix+part); }
    private static Csv csv(Path out,String name) throws IOException { return new Csv(out.resolve(name),name.endsWith(".gz")); }

    private static void writeInputs(Path out, Model m, List<Scenario> scenarios, Config cfg) throws IOException {
        try(Writer w=Files.newBufferedWriter(out.resolve("effective.properties"),StandardCharsets.UTF_8)) {
            for(String key:new TreeSet<>(cfg.p.stringPropertyNames())) w.write(key+"="+cfg.p.getProperty(key)+"\n");
        }
        try(Csv c=csv(out,"cell_vectors.csv")) {
            c.row("basis","index","x","y","z","unit","volume_A3");
            for(int i=0;i<3;i++) { Vec a=m.cell.a[i],b=m.cell.reciprocal[i];
                c.row("direct",i+1,a.x,a.y,a.z,"Angstrom",m.cell.volume);
                c.row("reciprocal",i+1,b.x,b.y,b.z,"inverse_Angstrom",m.cell.volume);
            }
        }
        try(Csv c=csv(out,"expanded_atoms.csv")) {
            c.row("atom_index","atom","orbit","element","partner","fx_unwrapped","fy_unwrapped","fz_unwrapped","x_A","y_A","z_A");
            for(int j=0;j<m.atoms.size();j++) {Atom a=m.atoms.get(j);Vec r=m.cell.cart(a.frac);c.row(j,a.name,a.orbit,a.element,a.partner,a.frac.x,a.frac.y,a.frac.z,r.x,r.y,r.z);}
        }
        try(Csv c=csv(out,"charges.csv")) {
            c.row("scenario","atom_index","atom","orbit","charge_e","net_cell_charge_e");
            for(Scenario s:scenarios) for(int j=0;j<m.atoms.size();j++) c.row(s.name,j,m.atoms.get(j).name,m.atoms.get(j).orbit,s.charges[j],Arrays.stream(s.charges).sum());
        }
        try(Csv c=csv(out,"near_neighbors_lt5A.csv")) {
            c.row("target","source","element","n1","n2","n3","dx_A","dy_A","dz_A","r_A");
            for(int ti:m.targets) {
                Atom target=m.atoms.get(ti); List<Object[]> rows=new ArrayList<>();
                for(Atom source:m.atoms) {
                    Vec d=source.frac.sub(target.frac);int[] bound=new int[3];
                    for(int i=0;i<3;i++) bound[i]=(int)Math.ceil(5*m.cell.inverseRows[i].norm()+Math.abs(d.at(i)))+1;
                    for(int n1=-bound[0];n1<=bound[0];n1++) for(int n2=-bound[1];n2<=bound[1];n2++) for(int n3=-bound[2];n3<=bound[2];n3++) {
                        Vec r=m.cell.cart(d.add(new Vec(n1,n2,n3)));double radius=r.norm();
                        if(radius>1e-10 && radius<5) rows.add(new Object[]{target.name,source.name,source.element,n1,n2,n3,r.x,r.y,r.z,radius});
                    }
                }
                rows.sort(Comparator.comparingDouble(r->(Double)r[9]));for(Object[] row:rows)c.row(row);
            }
        }
    }

    private static void writeResult(Csv csv, String scenario, String site, double alpha, Result r, Config cfg, boolean subtotals) throws IOException {
        List<Object> row=new ArrayList<>(List.of(scenario,site,alpha,cfg.d("realFactor")/alpha,cfg.d("reciprocalFactor")*alpha));
        if(subtotals) {append(row,r.real.v);append(row,r.reciprocal.v);append(row,r.self.v);}
        append(row,r.total.v); row.add(r.total.trace());row.add(r.eigen.xx);row.add(r.eigen.yy);row.add(r.eigen.zz);row.add(r.eigen.eta);row.add(r.eigen.p());
        row.add(r.eigen.zz*cfg.d("electrostaticFactor"));
        row.add(frequency(r.eigen,cfg.d("q63Barn"),cfg));row.add(frequency(r.eigen,cfg.d("q65Barn"),cfg));
        row.add(inherited(r.eigen,cfg));row.add(inherited(r.eigen,cfg)*Math.abs(cfg.d("q65Barn")/cfg.d("q63Barn")));
        row.add(r.eigen.zAxis.x);row.add(r.eigen.zAxis.y);row.add(r.eigen.zAxis.z);row.add(r.realTerms);row.add(r.reciprocalVectors);csv.row(row);
    }
    private static void resultHeader(Csv csv,boolean subtotals) throws IOException {
        List<Object> row=new ArrayList<>(List.of("scenario","site","ewald_alpha_invA","Rcut_A","Kcut_invA"));
        if(subtotals) {appendTensorHeader(row,"real_G_");appendTensorHeader(row,"reciprocal_G_");appendTensorHeader(row,"self_G_");}
        appendTensorHeader(row,"total_G_");row.addAll(List.of("trace_G","principal_Gxx","principal_Gyy","principal_Gzz_signed","eta","P_absGzz_sqrt1pluseta2over3",
                "Vzz_signed_V_per_m2","bare_63_MHz","bare_65_MHz","inherited_anilite_63_MHz_NOT_absolute","inherited_anilite_65_MHz_NOT_absolute","Vzz_axis_x","Vzz_axis_y","Vzz_axis_z","real_term_count","reciprocal_vector_count"));csv.row(row);
    }

    private static void assignments(Path out,List<Scenario> scenarios,Result[][] baseline,Config cfg) throws IOException {
        double[] observed=cfg.array("observedMHz");
        try(Csv all=csv(out,"assignments_all720.csv");Csv best=csv(out,"assignments_best6.csv")) {
            List<Object> header=new ArrayList<>(List.of("scenario","rank","omitted_site","s1_assigned_site","s2_assigned_site","s3_assigned_site","s4_assigned_site"));
            for(String p:new String[]{"observed","bare","inherited_anilite_NOT_absolute","fitted"}) for(int i=1;i<=4;i++) header.add(p+"_s"+i+"_MHz");
            header.addAll(List.of("fitted_multiplier_of_bare","fitted_multiplier_of_inherited","omitted_bare_MHz","omitted_inherited_NOT_absolute_MHz","omitted_fitted_MHz","fitted_RMSE_MHz","fitted_max_abs_residual_MHz"));all.row(header);best.row(header);
            for(int s=0;s<scenarios.size();s++) {
                List<Assignment> list=new ArrayList<>();
                for(int omitted=0;omitted<5;omitted++) {
                    int[] choices=new int[4];int j=0;for(int i=0;i<5;i++)if(i!=omitted)choices[j++]=i;
                    List<int[]> permutations=new ArrayList<>();permute(choices,0,permutations);
                    for(int[] p:permutations) {
                        double[] f=new double[4],inh=new double[4],fitted=new double[4];String[] map=new String[4];double dot=0,den=0,idot=0,iden=0;
                        for(int i=0;i<4;i++) {Eigen e=baseline[s][p[i]].eigen;f[i]=frequency(e,cfg.d("q63Barn"),cfg);inh[i]=inherited(e,cfg);map[i]=ORBITS[p[i]];dot+=f[i]*observed[i];den+=f[i]*f[i];idot+=inh[i]*observed[i];iden+=inh[i]*inh[i];}
                        double scale=dot/den,is=idot/iden,sum=0,max=0;
                        for(int i=0;i<4;i++) {fitted[i]=scale*f[i];double d=fitted[i]-observed[i];sum+=d*d;max=Math.max(max,Math.abs(d));}
                        Eigen missing=baseline[s][omitted].eigen;double missingBare=frequency(missing,cfg.d("q63Barn"),cfg);
                        list.add(new Assignment(scenarios.get(s).name,map,omitted,f,inh,fitted,scale,is,missingBare,inherited(missing,cfg),scale*missingBare,Math.sqrt(sum/4),max));
                    }
                }
                list.sort(Comparator.comparingDouble(Assignment::rmse));int rank=0;
                for(Assignment a:list) {
                    List<Object> row=new ArrayList<>(List.of(a.scenario,++rank,ORBITS[a.omitted]));row.addAll(Arrays.asList(a.mapping));append(row,observed);append(row,a.bare);append(row,a.inherited);append(row,a.fitted);
                    row.addAll(List.of(a.bareScale,a.inheritedScale,a.omittedBare,a.omittedInherited,a.omittedFitted,a.rmse,a.maxError));all.row(row);if(rank==1)best.row(row);
                }
                if(list.size()!=120)throw new IllegalStateException("Assignment enumeration error");
            }
        }
    }
    private static void permute(int[] p,int from,List<int[]> out) {
        if(from==p.length) {out.add(p.clone());return;}
        for(int i=from;i<p.length;i++) {int t=p[from];p[from]=p[i];p[i]=t;permute(p,from+1,out);t=p[from];p[from]=p[i];p[i]=t;}
    }

    private static Tensor direct(Model m,Scenario s,int ti,double cutoff) {
        TensorSum sum=new TensorSum();Vec target=m.atoms.get(ti).frac;
        for(int j=0;j<m.atoms.size();j++) {
            Vec delta=m.atoms.get(j).frac.sub(target);int[] lo=new int[3],hi=new int[3];
            for(int i=0;i<3;i++) {double b=cutoff*m.cell.inverseRows[i].norm();lo[i]=(int)Math.ceil(-delta.at(i)-b-1e-12);hi[i]=(int)Math.floor(-delta.at(i)+b+1e-12);}
            for(int n1=lo[0];n1<=hi[0];n1++)for(int n2=lo[1];n2<=hi[1];n2++)for(int n3=lo[2];n3<=hi[2];n3++) {
                Vec r=m.cell.cart(delta.add(new Vec(n1,n2,n3)));double r2=r.norm2();if(r2<1e-24||r2>cutoff*cutoff)continue;
                double r3=r2*Math.sqrt(r2);sum.add(hessian(r,3/(r3*r2),-1/r3,s.charges[j]));
            }
        }
        return sum.tensor();
    }

    private static void directSums(Path out,Model m,Scenario uniform,Result[] reference,Config cfg) throws IOException {
        try(Csv c=csv(out,"direct_spherical_crosscheck.csv")) {
            c.row("site","sphere_radius_A","direct_Gzz_signed","direct_eta","direct_bare_63_MHz","ewald_Gzz_signed","ewald_bare_63_MHz","max_tensor_difference_G","note");
            for(double radius:cfg.array("directRadii")) for(int o=0;o<5;o++) {
                Tensor t=direct(m,uniform,m.targets[o],radius);Eigen e=eigen(t);
                c.row(ORBITS[o],radius,e.zz,e.eta,frequency(e,cfg.d("q63Barn"),cfg),reference[o].eigen.zz,frequency(reference[o].eigen,cfg.d("q63Barn"),cfg),t.maxDifference(reference[o].total),"spherical_truncation_crosscheck_not_primary_result");
            }
        }
    }

    private static void sensitivityRow(Csv c,String kind,String varied,double value,String evaluation,int o,Result r,Result ref,Config cfg) throws IOException {
        double f=frequency(r.eigen,cfg.d("q63Barn"),cfg), f0=frequency(ref.eigen,cfg.d("q63Barn"),cfg);
        c.row(kind,varied,value,evaluation,ORBITS[o],r.eigen.zz,r.eigen.eta,r.eigen.p(),f,inherited(r.eigen,cfg),f0,f/f0-1,r.total.maxDifference(ref.total),"model_sensitivity_NOT_measured_temperature_dependence");
    }
    private static void sensitivities(Path out,Model m,Scenario uniform,Result[] baseline,Config cfg,Csv checks) throws IOException {
        try(Csv c=csv(out,"sensitivities_uniform_charges.csv")) {
            c.row("variation","varied_parameter_or_orbit","value","evaluation","site","Gzz_signed","eta","P","bare_63_MHz","inherited_NOT_absolute_63_MHz","reference_bare_63_MHz","relative_frequency_change","max_tensor_difference_from_reference_G","interpretation");
            for(double lambda:cfg.array("uniformLambdas")) for(int o=0;o<5;o++) {
                double factor=Math.pow(lambda,-3);Tensor t=baseline[o].total.scale(factor);
                Result r=new Result(baseline[o].real.scale(factor),baseline[o].reciprocal.scale(factor),baseline[o].self.scale(factor),t,eigen(t),0,0);
                sensitivityRow(c,"uniform_linear_scale","lambda",lambda,"analytic_lambda_minus3",o,r,baseline[o],cfg);
            }
            Model scaled=m.scale(0.99);
            for(int o=0;o<5;o++) {
                Result r=ewald(scaled,uniform,scaled.targets[o],cfg.d("mainAlpha"),cfg,null,null);
                sensitivityRow(c,"uniform_linear_scale","lambda",0.99,"full_ewald_check",o,r,baseline[o],cfg);
                check(checks,"uniform_scale_0.99_"+ORBITS[o],r.total.maxDifference(baseline[o].total.scale(Math.pow(0.99,-3))),cfg.d("convergenceTol"));
            }
            Model alternate=model(cfg,cfg.d("angleAlternativeDeg"));
            for(int o=0;o<5;o++) sensitivityRow(c,"rhombohedral_angle","angleDeg",cfg.d("angleAlternativeDeg"),"full_ewald",o,
                    ewald(alternate,uniform,alternate.targets[o],cfg.d("mainAlpha"),cfg,null,null),baseline[o],cfg);
            for(int moved=1;moved<5;moved++) for(int sign:new int[]{-1,1}) {
                double delta=sign*cfg.d("displacementAngstrom");Model changed=m.displacePair(moved,new Vec(0,0,delta));
                for(int o=0;o<5;o++) sensitivityRow(c,"inversion_paired_axial_displacement_A",ORBITS[moved],delta,"full_ewald",o,
                        ewald(changed,uniform,changed.targets[o],cfg.d("mainAlpha"),cfg,null,null),baseline[o],cfg);
            }
            if(cfg.b("transverseExample")) for(int sign:new int[]{-1,1}) {
                double delta=sign*cfg.d("displacementAngstrom");Model changed=m.displacePair(1,new Vec(delta,0,0));
                for(int o=0;o<5;o++) sensitivityRow(c,"illustrative_transverse_pair_displacement_A","DgB",delta,"full_ewald",o,
                        ewald(changed,uniform,changed.targets[o],cfg.d("mainAlpha"),cfg,null,null),baseline[o],cfg);
            }
        }
    }

    private static int failedChecks=0;
    private static void check(Csv c,String name,double absoluteError,double tolerance) throws IOException {
        boolean ok=Double.isFinite(absoluteError)&&absoluteError<=tolerance;if(!ok)failedChecks++;
        c.row(name,absoluteError,tolerance,ok?"PASS":"FAIL");
    }
    private static void independentChecks(Csv c,Model m,List<Scenario> scenarios,Result[][] results,Config cfg) throws IOException {
        double[][] erfcReferences={{0,1},{0.5,0.4795001221869535},{1,0.15729920705028513},{1.5,0.033894853524689274},{2,0.004677734981047265},{3,0.000022090496998585438},{5,1.5374597944280351e-12},{6,2.1519736712498916e-17}};
        for(double[] pair:erfcReferences) check(c,"erfc_relative_x="+pair[0],Math.abs(erfc(pair[0])/pair[1]-1),2e-14);
        for(int s=0;s<scenarios.size();s++) for(int o=0;o<5;o++) {
            Result r=results[s][o];check(c,"trace_"+scenarios.get(s).name+"_"+ORBITS[o],Math.abs(r.total.trace()),cfg.d("convergenceTol"));
            check(c,"axial_symmetry_"+scenarios.get(s).name+"_"+ORBITS[o],Math.max(Math.abs(r.total.v[0]-r.total.v[1]),Math.max(Math.abs(r.total.v[3]),Math.max(Math.abs(r.total.v[4]),Math.abs(r.total.v[5])))),cfg.d("convergenceTol"));
        }
        for(int o=1;o<5;o++) {
            Result partner=ewald(m,scenarios.get(0),m.targets[o]+1,cfg.d("mainAlpha"),cfg,null,null);
            check(c,"inversion_partner_"+ORBITS[o],partner.total.maxDifference(results[0][o].total),cfg.d("convergenceTol"));
        }
        // Independent reference values are ONLY assertions; they are never used to generate result tables.
        boolean nominal=Math.abs(cfg.d("a")-16.16)<1e-12&&Math.abs(cfg.d("angleDeg")-13.933333333333334)<1e-12;
        double[] positions={0.5,0.06,0.133,0.25,0.35};for(int o=0;o<5;o++)nominal&=Math.abs(cfg.d("x."+ORBITS[o])-positions[o])<1e-12;
        if(nominal) {
            double[][] ref={
                    {-0.0290586668021,0.2174141780431,0.4137980909365,-0.0584572731234,0.0175286755961},
                    {-0.1699618089853,0.2630588405144,0.3510351492399,-0.0148580376443,0.0551120979850},
                    {0.0129979021174,0.1717406274693,0.4519117266279,-0.1071258548199,0.0301952742837},
                    {-0.0836112579880,0.2673262577970,0.3530826142616,-0.0273484961180,0.0550589759541},
                    {0.0116237389859,0.1694168602073,0.4337796877401,-0.1048645602283,0.0520975547807},
                    {0.0132065207677,0.2383506354634,0.4477998059648,-0.0162897990668,-0.0860288138287}};
            for(int s=0;s<6;s++)for(int o=0;o<5;o++)check(c,"independent_Gzz_"+scenarios.get(s).name+"_"+ORBITS[o],Math.abs(results[s][o].eigen.zz-ref[s][o]),5e-10);
        } else c.row("independent_nominal_reference", "not_applicable_to_modified_geometry","", "SKIP");
    }

    public static void main(String[] args) throws Exception {
        if(Arrays.asList(args).contains("--help")) {
            System.out.println("Usage: java -jar digenite-efg-1.0.0.jar [--config=model.properties] [--output=results] [--key=value]\nAll settings are listed in model.properties. No libraries or network access are needed.\nRaw G units: Angstrom^-3 for dimensionless charges z=q/e; bare frequencies are separate from inherited/fitted values.");return;
        }
        Locale.setDefault(Locale.ROOT);Config cfg=new Config(args);Path out=cfg.output();Files.createDirectories(out);
        Model m=model(cfg,cfg.d("angleDeg"));List<Scenario> scenarios=scenarios(m);writeInputs(out,m,scenarios,cfg);
        Result[][] baseline=new Result[6][5];double mainAlpha=cfg.d("mainAlpha");
        try(Csv tensors=csv(out,"tensors_and_frequencies_30.csv");Csv convergence=csv(out,"alpha_convergence_subtotals.csv");Csv checks=csv(out,"validation_checks.csv")) {
            resultHeader(tensors,false);resultHeader(convergence,true);checks.row("check","absolute_or_relative_error_as_named","tolerance","status");
            Csv rawReal=cfg.b("rawTerms")?csv(out,"raw_real_uniform_mainalpha.csv.gz"):null;
            Csv rawReciprocal=cfg.b("rawTerms")?csv(out,"raw_reciprocal_uniform_mainalpha.csv.gz"):null;
            try {
                if(rawReal!=null) {List<Object> h=new ArrayList<>(List.of("target","source_index","source","orbit","charge_e","n1","n2","n3","rx_A","ry_A","rz_A","r_A"));appendTensorHeader(h,"term_G_");rawReal.row(h);}
                if(rawReciprocal!=null) {List<Object> h=new ArrayList<>(List.of("target","h","k","l","source_index","source","orbit","charge_e","kx_invA","ky_invA","kz_invA","k_invA","phase_rad","q_cos_phase"));appendTensorHeader(h,"term_G_");rawReciprocal.row(h);}
                for(int s=0;s<6;s++)for(int o=0;o<5;o++) {
                    baseline[s][o]=ewald(m,scenarios.get(s),m.targets[o],mainAlpha,cfg,s==0?rawReal:null,s==0?rawReciprocal:null);
                    writeResult(tensors,scenarios.get(s).name,ORBITS[o],mainAlpha,baseline[s][o],cfg,false);
                    for(double alpha:cfg.array("ewaldAlphas")) {
                        Result r=Math.abs(alpha-mainAlpha)<1e-14?baseline[s][o]:ewald(m,scenarios.get(s),m.targets[o],alpha,cfg,null,null);
                        writeResult(convergence,scenarios.get(s).name,ORBITS[o],alpha,r,cfg,true);
                        if(cfg.b("verify"))check(checks,"alpha_"+alpha+"_"+scenarios.get(s).name+"_"+ORBITS[o],r.total.maxDifference(baseline[s][o].total),cfg.d("convergenceTol"));
                    }
                }
            } finally {if(rawReal!=null)rawReal.close();if(rawReciprocal!=null)rawReciprocal.close();}
            if(cfg.b("verify"))independentChecks(checks,m,scenarios,baseline,cfg);
            if(cfg.b("sensitivities"))sensitivities(out,m,scenarios.get(0),baseline[0],cfg,checks);
            assignments(out,scenarios,baseline,cfg);
            if(cfg.b("directSums"))directSums(out,m,scenarios.get(0),baseline[0],cfg);
        }
        System.out.printf(Locale.ROOT,"Computed %d atoms, %d charge scenarios, 30 EFG tensors, 720 assignments. Cell volume %.12f A^3.%n",m.atoms.size(),scenarios.size(),m.cell.volume);
        for(int o=0;o<5;o++)System.out.printf(Locale.ROOT,"%s: Gzz=%+.12g eta=%.6g bare63=%.12g MHz%n",ORBITS[o],baseline[0][o].eigen.zz,baseline[0][o].eigen.eta,frequency(baseline[0][o].eigen,cfg.d("q63Barn"),cfg));
        System.out.println("Results: "+out.toAbsolutePath()+"; failed validation checks: "+failedChecks);
        if(failedChecks>0)throw new IllegalStateException("Numerical validation failed; inspect validation_checks.csv");
    }
}
