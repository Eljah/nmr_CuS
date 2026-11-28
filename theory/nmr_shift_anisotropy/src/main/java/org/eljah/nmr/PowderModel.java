package org.eljah.nmr;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/** Archive-compatible powder spectrum, convolution and variable-projection NNLS. */
final class PowderModel {
    static final double GAMMA_RATIO = 1.071210656566704;
    static final double FIELD_REFERENCE_T = 6.720566;
    static final double[][] NQR = {{11.97,11.07},{15.27,14.13},{16.18,14.97},{16.96,15.69}};
    static final double STEP = .01, REL_MIN = -18;
    static final int BINS = 3601;
    static final String[] PARAMETER_NAMES = {"eta_s1","eta_s2","eta_s3","eta_s4",
        "nuZ63_s1_MHz","nuZ63_s2_MHz","nuZ63_s3_MHz","nuZ63_s4_MHz",
        "sigma_MHz","rho500","rho20","kappa_s1","kappa_s2","kappa_s3","kappa_s4"};

    record Grid(int nu, int nphi) { }
    record Key(int site, long eta, long centre, long kappa, long sigma, long field, boolean rfTensor) { }

    static final class Data {
        final String series;
        final double[] x, observed, normalized, t;
        final boolean[] used;
        final int[] sourceRow;
        final double scale, sumT, sumTT;
        final int n;
        Data(String series, List<double[]> rows) {
            this.series = series;
            int count = rows.size();
            x = new double[count]; observed = new double[count]; normalized = new double[count];
            t = new double[count]; used = new boolean[count]; sourceRow = new int[count];
            double maximum = Double.NEGATIVE_INFINITY, st = 0, stt = 0;
            int included = 0;
            for (int i = 0; i < count; i++) {
                double[] r = rows.get(i);
                x[i] = r[0]; observed[i] = r[1]; sourceRow[i] = (int)r[2];
                t[i] = (x[i] - 78) / 6;
                used[i] = Math.abs(x[i] - 75.941) > .12 && Math.abs(x[i] - 81.346) > .12;
                if (used[i]) {
                    included++; maximum = Math.max(maximum, observed[i]);
                    st += t[i]; stt += t[i] * t[i];
                }
            }
            if (!(maximum > 0) || included < 6) throw new IllegalArgumentException("Insufficient positive data.");
            scale = maximum; n = included; sumT = st; sumTT = stt;
            for (int i = 0; i < count; i++) normalized[i] = observed[i] / scale;
        }
    }

    record SeriesResult(double[] predictedNormalized, double[] residualUsed, double[] amplitudesNormalized,
                        double[] baselineNormalized, double sse, double[][] design) { }
    record Evaluation(double sse, double[] residual, SeriesResult[] series) { }

    final Data[] data;
    final Grid grid;
    final boolean rfTensor;
    private final Map<Key,double[][]> rawCache = lru(320);
    private final Map<Key,double[][]> smoothCache = lru(640);
    long evaluations, profileCalculations;

    PowderModel(Data[] data, Grid grid, boolean rfTensor) {
        this.data = data; this.grid = grid; this.rfTensor = rfTensor;
    }

    private static <K,V> Map<K,V> lru(int maximum) {
        return new LinkedHashMap<>(maximum, .75f, true) {
            @Override protected boolean removeEldestEntry(Map.Entry<K,V> eldest) { return size() > maximum; }
        };
    }

    static Data[] readData(Path path) throws IOException {
        List<String> lines = Files.readAllLines(path);
        String[] headers = lines.get(0).split(",", -1);
        Map<String,Integer> index = new LinkedHashMap<>();
        for (int i = 0; i < headers.length; i++) index.put(headers[i].trim(), i);
        List<double[]> longRows = new ArrayList<>(), shortRows = new ArrayList<>();
        for (int row = 1; row < lines.size(); row++) {
            if (lines.get(row).isBlank()) continue;
            String[] fields = lines.get(row).split(",", -1);
            double x = parse(fields[index.get("Freq")]);
            if (!Double.isFinite(x) || Math.abs(x - 81.4) <= 1e-7 + 1e-5 * 81.4) continue;
            int origin = (int)parse(fields[index.get("origin_row")]);
            double y500 = parse(fields[index.get("IntQ1")]), y20 = parse(fields[index.get("IntQ")]);
            if (Double.isFinite(y500)) longRows.add(new double[]{x, y500, origin});
            if (Double.isFinite(y20)) shortRows.add(new double[]{x, y20, origin});
        }
        return new Data[]{new Data("500", longRows), new Data("20", shortRows)};
    }

    private static double parse(String value) {
        try { return Double.parseDouble(value.trim()); } catch (NumberFormatException e) { return Double.NaN; }
    }

    Evaluation evaluate(double[] p) {
        evaluations++;
        double[][][] shape = new double[4][][];
        for (int site = 0; site < 4; site++) shape[site] = profiles(site, p[site], p[4 + site], p[11 + site], p[8], 1);
        SeriesResult[] result = new SeriesResult[2];
        double[] residual = new double[data[0].n + data[1].n];
        int offset = 0;
        double total = 0;
        for (int series = 0; series < 2; series++) {
            Data d = data[series];
            double[][] design = new double[d.x.length][4];
            for (int row = 0; row < d.x.length; row++) for (int site = 0; site < 4; site++) {
                design[row][site] = interpolate(shape[site][0], d.x[row] - p[4 + site])
                    + p[9 + series] * interpolate(shape[site][1], d.x[row] - p[4 + site] * GAMMA_RATIO);
            }
            result[series] = solveLinear(d, design);
            System.arraycopy(result[series].residualUsed(), 0, residual, offset, d.n);
            offset += d.n; total += result[series].sse();
        }
        return new Evaluation(total, residual, result);
    }

    double[][] profiles(int site, double eta, double centre, double kappa, double sigma, double fieldRatio) {
        Key smoothKey = new Key(site, Double.doubleToLongBits(eta), Double.doubleToLongBits(centre),
            Double.doubleToLongBits(kappa), Double.doubleToLongBits(sigma), Double.doubleToLongBits(fieldRatio), rfTensor);
        double[][] cached = smoothCache.get(smoothKey);
        if (cached != null) return cached;
        Key rawKey = new Key(site, smoothKey.eta(), smoothKey.centre(), smoothKey.kappa(), 0, smoothKey.field(), rfTensor);
        double[][] raw = rawCache.get(rawKey);
        if (raw == null) {
            raw = new double[2][];
            for (int isotope = 0; isotope < 2; isotope++) raw[isotope] = sticks(eta, NQR[site][isotope],
                centre * fieldRatio * (isotope == 0 ? 1 : GAMMA_RATIO), kappa, grid, rfTensor);
            rawCache.put(rawKey, raw); profileCalculations++;
        }
        double[][] smoothed = {gaussian(raw[0], sigma), gaussian(raw[1], sigma)};
        smoothCache.put(smoothKey, smoothed);
        return smoothed;
    }

    static double[] sticks(double eta, double nqr, double centre, double kappa, Grid grid, boolean rfTensor) {
        double nuQ = nqr / Math.sqrt(1 + eta * eta / 3);
        double gp = 1 - kappa / 3, gz = 1 + 2 * kappa / 3;
        double[] histogram = new double[BINS];
        SpinHamiltonian spin = new SpinHamiltonian();
        // Axial EFG and coaxial axial shift are azimuth-independent exactly.
        int nphi = eta == 0 ? 1 : grid.nphi();
        double norm = 1.0 / (grid.nu() * nphi * STEP);
        for (int iu = 0; iu < grid.nu(); iu++) {
            double nz = (iu + .5) / grid.nu();
            double radial = Math.sqrt(1 - nz * nz);
            for (int ip = 0; ip < nphi; ip++) {
                double phi = (ip + .5) * Math.PI / (2 * nphi);
                double nx = radial * Math.cos(phi), ny = radial * Math.sin(phi);
                spin.solve(centre * gp * nx, centre * gp * ny, centre * gz * nz, nuQ, eta);
                for (int a = 0; a < 3; a++) for (int b = a + 1; b < 4; b++) {
                    double z = (spin.eigenvalues[b] - spin.eigenvalues[a] - centre - REL_MIN) / STEP;
                    int lower = (int)Math.floor(z);
                    if (lower < 0 || lower >= BINS - 1) continue;
                    double fraction = z - lower;
                    double weight = spin.weight(a, b, nx, ny, nz, rfTensor ? gp : 1, rfTensor ? gz : 1) * norm;
                    histogram[lower] += weight * (1 - fraction);
                    histogram[lower + 1] += weight * fraction;
                }
            }
        }
        return histogram;
    }

    static double[] gaussian(double[] input, double sigma) {
        double binsSigma = sigma / STEP;
        int radius = (int)(4 * binsSigma + .5);
        double[] kernel = new double[2 * radius + 1];
        double sum = 0;
        for (int i = -radius; i <= radius; i++) { kernel[i + radius] = Math.exp(-.5 * i * i / (binsSigma * binsSigma)); sum += kernel[i + radius]; }
        for (int i = 0; i < kernel.length; i++) kernel[i] /= sum;
        double[] output = new double[input.length];
        // Scatter is equivalent to symmetric constant-zero convolution, but avoids zero bins.
        for (int i = 0; i < input.length; i++) if (input[i] != 0) {
            int from = Math.max(-radius, -i), to = Math.min(radius, input.length - i - 1);
            for (int k = from; k <= to; k++) output[i + k] += input[i] * kernel[k + radius];
        }
        return output;
    }

    static double interpolate(double[] profile, double relativeFrequency) {
        double coordinate = (relativeFrequency - REL_MIN) / STEP;
        int index = (int)Math.floor(coordinate);
        if (index < 0 || index >= profile.length - 1) return 0;
        double fraction = coordinate - index;
        return profile[index] * (1 - fraction) + profile[index + 1] * fraction;
    }

    private static SeriesResult solveLinear(Data d, double[][] design) {
        double[] sumA = new double[4], sumTA = new double[4];
        double sumY = 0, sumTY = 0;
        for (int row = 0; row < d.x.length; row++) if (d.used[row]) {
            sumY += d.normalized[row]; sumTY += d.t[row] * d.normalized[row];
            for (int col = 0; col < 4; col++) { sumA[col] += design[row][col]; sumTA[col] += d.t[row] * design[row][col]; }
        }
        double denominator = d.n * d.sumTT - d.sumT * d.sumT;
        double by0 = (sumY * d.sumTT - sumTY * d.sumT) / denominator;
        double by1 = (d.n * sumTY - d.sumT * sumY) / denominator;
        double[] ba0 = new double[4], ba1 = new double[4];
        for (int c = 0; c < 4; c++) {
            ba0[c] = (sumA[c] * d.sumTT - sumTA[c] * d.sumT) / denominator;
            ba1[c] = (d.n * sumTA[c] - d.sumT * sumA[c]) / denominator;
        }
        double[][] gram = new double[4][4];
        double[] rhs = new double[4], a = new double[4];
        for (int row = 0; row < d.x.length; row++) if (d.used[row]) {
            double y = d.normalized[row] - by0 - by1 * d.t[row];
            for (int c = 0; c < 4; c++) a[c] = design[row][c] - ba0[c] - ba1[c] * d.t[row];
            for (int c = 0; c < 4; c++) {
                rhs[c] += a[c] * y;
                for (int k = 0; k < 4; k++) gram[c][k] += a[c] * a[k];
            }
        }
        // Exact active-set enumeration of a four-variable nonnegative LS problem.
        double[] best = new double[4];
        double bestObjective = 0;
        for (int subset = 1; subset < 16; subset++) {
            int size = Integer.bitCount(subset);
            int[] columns = new int[size];
            for (int c = 0, k = 0; c < 4; c++) if ((subset & (1 << c)) != 0) columns[k++] = c;
            double[][] g = new double[size][size]; double[] y = new double[size];
            for (int c = 0; c < size; c++) {
                y[c] = rhs[columns[c]];
                for (int k = 0; k < size; k++) g[c][k] = gram[columns[c]][columns[k]];
            }
            double[] solution = solve(g, y);
            if (solution == null) continue;
            boolean feasible = true;
            for (double value : solution) if (value < -1e-10) feasible = false;
            if (!feasible) continue;
            double objective = 0;
            for (int c = 0; c < size; c++) {
                solution[c] = Math.max(0, solution[c]);
                objective -= 2 * rhs[columns[c]] * solution[c];
                for (int k = 0; k < size; k++) objective += solution[c] * gram[columns[c]][columns[k]] * solution[k];
            }
            if (objective < bestObjective) {
                bestObjective = objective; Arrays.fill(best, 0);
                for (int c = 0; c < size; c++) best[columns[c]] = solution[c];
            }
        }
        double b0 = by0, b1 = by1;
        for (int c = 0; c < 4; c++) { b0 -= best[c] * ba0[c]; b1 -= best[c] * ba1[c]; }
        double[] prediction = new double[d.x.length], residual = new double[d.n];
        double sse = 0;
        for (int row = 0, usedIndex = 0; row < d.x.length; row++) {
            prediction[row] = b0 + b1 * d.t[row];
            for (int c = 0; c < 4; c++) prediction[row] += design[row][c] * best[c];
            if (d.used[row]) { double r = prediction[row] - d.normalized[row]; residual[usedIndex++] = r; sse += r * r; }
        }
        return new SeriesResult(prediction, residual, best, new double[]{b0,b1}, sse, design);
    }

    /** Partial-pivot dense solve; null means a numerically singular active set. */
    static double[] solve(double[][] matrix, double[] rhs) {
        int n = rhs.length; double[][] a = new double[n][n + 1];
        double maximum = 0;
        for (int i = 0; i < n; i++) { for (int j = 0; j < n; j++) { a[i][j] = matrix[i][j]; maximum = Math.max(maximum, Math.abs(a[i][j])); } a[i][n] = rhs[i]; }
        for (int col = 0; col < n; col++) {
            int pivot = col;
            for (int i = col + 1; i < n; i++) if (Math.abs(a[i][col]) > Math.abs(a[pivot][col])) pivot = i;
            if (Math.abs(a[pivot][col]) < 1e-13 * Math.max(1e-20, maximum)) return null;
            double[] swap = a[col]; a[col] = a[pivot]; a[pivot] = swap;
            for (int row = col + 1; row < n; row++) {
                double factor = a[row][col] / a[col][col];
                for (int j = col; j <= n; j++) a[row][j] -= factor * a[col][j];
            }
        }
        double[] x = new double[n];
        for (int i = n - 1; i >= 0; i--) { double value = a[i][n]; for (int j = i + 1; j < n; j++) value -= a[i][j] * x[j]; x[i] = value / a[i][i]; }
        return x;
    }
}
