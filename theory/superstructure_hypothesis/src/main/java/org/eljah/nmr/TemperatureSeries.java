package org.eljah.nmr;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HexFormat;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;

/**
 * Computes temperature-dependent scalar summaries of archived NQR frequencies.
 * Missing cells remain missing. Metrics use complete quartets only; the fixed
 * hypothetical cyclic order is [s1,s2,s4,s3] at every temperature. This class
 * never optimizes or exchanges branch assignments at crossings.
 */
public final class TemperatureSeries {
    private static final String[] LABELS = {"s1", "s2", "s3", "s4"};
    private static final String[] COLORS = {"#00858d", "#c66b13", "#2864ac", "#934a9d"};
    private static final int[] CYCLIC_ORDER = {0, 1, 3, 2};

    private TemperatureSeries() { }

    private record Measurement(double temperature, double[] frequencies, String sourceRow) {
        boolean complete() {
            return Arrays.stream(frequencies).allMatch(Double::isFinite);
        }
    }
    private record Metrics(Measurement measurement, double mean, double span, double rms,
                           double normalizedRms, double maxMinRatio,
                           double b, double c, double d2) {
        double rmse() { return Math.abs(d2); }
    }
    private record Check(String name, double error, double tolerance) {
        boolean pass() { return Double.isFinite(error) && Math.abs(error) <= tolerance; }
    }

    public static void main(String[] args) throws IOException {
        Path input = Path.of("inputs/temperature_series.csv");
        Path output = Path.of("results");
        Path figure = Path.of("figures/temperature_branches.svg");
        double maxGap = 40;
        for (String arg : args) {
            if (arg.equals("--help") || arg.equals("-h")) {
                System.out.println("Usage: java -cp target/classes org.eljah.nmr.TemperatureSeries "
                    + "[--input=inputs/temperature_series.csv] [--output=results] "
                    + "[--figure=figures/temperature_branches.svg] [--max-gap-k=40]");
                System.out.println("CSV columns: T_K,s1_MHz,s2_MHz,s3_MHz,s4_MHz,source_row. "
                    + "Empty frequency cells are missing. Only complete quartets enter metrics.");
                System.out.println("Fixed cyclic order [s1,s2,s4,s3]; no reassignment at any temperature.");
                return;
            } else if (arg.startsWith("--input=")) {
                input = Path.of(arg.substring("--input=".length()));
            } else if (arg.startsWith("--output=")) {
                output = Path.of(arg.substring("--output=".length()));
            } else if (arg.startsWith("--figure=")) {
                figure = Path.of(arg.substring("--figure=".length()));
            } else if (arg.startsWith("--max-gap-k=")) {
                maxGap = Double.parseDouble(arg.substring("--max-gap-k=".length()));
            } else throw new IllegalArgumentException("Unknown argument: " + arg + "; use --help.");
        }
        if (!Double.isFinite(maxGap) || maxGap <= 0) {
            throw new IllegalArgumentException("The graphical connection gap limit must be positive in K.");
        }
        byte[] bytes = Files.readAllBytes(input);
        List<Measurement> rows = readMeasurements(new String(bytes, StandardCharsets.UTF_8));
        List<Metrics> metrics = new ArrayList<>();
        for (Measurement row : rows) if (row.complete()) metrics.add(metrics(row));
        List<Check> checks = validate(rows, metrics);
        Files.createDirectories(output);
        writeMetrics(output.resolve("temperature_metrics.csv"), metrics);
        writeSummary(output.resolve("temperature_summary.json"), input, sha256(bytes), rows, metrics, maxGap, checks);
        writeFigure(figure, rows, metrics.size(), maxGap);
        long failures = checks.stream().filter(c -> !c.pass()).count();
        System.out.printf(Locale.ROOT,
            "Temperature series: %d rows, %d complete quartets, %d incomplete rows; T=%.6g to %.6g K.%n",
            rows.size(), metrics.size(), rows.size() - metrics.size(),
            rows.get(0).temperature(), rows.get(rows.size() - 1).temperature());
        if (!metrics.isEmpty()) {
            Metrics first = metrics.get(0), last = metrics.get(metrics.size() - 1);
            System.out.printf(Locale.ROOT,
                "Complete-quartet endpoints %.6g -> %.6g K:%n"
                + "  mean %.12f -> %.12f MHz; span %.12f -> %.12f MHz%n"
                + "  rms splitting %.12f -> %.12f MHz; normalized rms %.12f -> %.12f%n"
                + "  fixed-order first-harmonic RMSE %.12f -> %.12f MHz%n",
                first.measurement().temperature(), last.measurement().temperature(),
                first.mean(), last.mean(), first.span(), last.span(), first.rms(), last.rms(),
                first.normalizedRms(), last.normalizedRms(), first.rmse(), last.rmse());
        }
        System.out.printf(Locale.ROOT, "Temperature consistency checks: %d PASS, %d FAIL.%n",
            checks.size() - failures, failures);
        System.out.println("Fixed hypothetical order [s1,s2,s4,s3]; missing values are not imputed or interpreted as merging.");
        System.out.println("Results: " + output.toAbsolutePath().normalize());
        System.out.println("Figure: " + figure.toAbsolutePath().normalize());
        if (failures != 0) throw new IllegalStateException("Temperature consistency checks failed.");
    }

    private static List<Measurement> readMeasurements(String text) {
        if (text.startsWith("\uFEFF")) text = text.substring(1);
        List<List<String>> parsed = parseCsv(text);
        if (parsed.isEmpty()) throw new IllegalArgumentException("The temperature CSV is empty.");
        Map<String, Integer> columns = new LinkedHashMap<>();
        for (int i = 0; i < parsed.get(0).size(); i++) {
            String key = parsed.get(0).get(i).trim();
            if (columns.put(key, i) != null) throw new IllegalArgumentException("Duplicate CSV column: " + key);
        }
        for (String key : new String[]{"T_K", "s1_MHz", "s2_MHz", "s3_MHz", "s4_MHz", "source_row"}) {
            if (!columns.containsKey(key)) throw new IllegalArgumentException("Missing CSV column: " + key);
        }
        List<Measurement> result = new ArrayList<>();
        for (int i = 1; i < parsed.size(); i++) {
            List<String> cells = parsed.get(i);
            if (cells.size() == 1 && cells.get(0).isBlank()) continue;
            if (cells.size() != parsed.get(0).size()) {
                throw new IllegalArgumentException("CSV record " + (i + 1) + " has the wrong field count.");
            }
            double temperature = finiteNumber(cells.get(columns.get("T_K")), "T_K", i + 1);
            if (temperature < 0) throw new IllegalArgumentException("Temperature cannot be negative in K.");
            double[] f = new double[4];
            for (int branch = 0; branch < 4; branch++) {
                String key = LABELS[branch] + "_MHz";
                String value = cells.get(columns.get(key)).trim();
                f[branch] = value.isEmpty() ? Double.NaN : finiteNumber(value, key, i + 1);
                if (Double.isFinite(f[branch]) && f[branch] <= 0) {
                    throw new IllegalArgumentException("Frequency must be positive at CSV record " + (i + 1));
                }
            }
            result.add(new Measurement(temperature, f, cells.get(columns.get("source_row"))));
        }
        if (result.isEmpty()) throw new IllegalArgumentException("No temperature records found.");
        if (result.stream().noneMatch(r -> Arrays.stream(r.frequencies()).anyMatch(Double::isFinite))) {
            throw new IllegalArgumentException("No measured frequencies found.");
        }
        result.sort(Comparator.comparingDouble(Measurement::temperature));
        return result;
    }

    private static double finiteNumber(String value, String column, int row) {
        try {
            double result = Double.parseDouble(value.trim());
            if (Double.isFinite(result)) return result;
        } catch (NumberFormatException ignored) { }
        throw new IllegalArgumentException("Invalid " + column + " at CSV record " + row + ": " + value);
    }

    // Handles empty cells, quoted commas, doubled quotes and CRLF/LF records.
    private static List<List<String>> parseCsv(String text) {
        List<List<String>> records = new ArrayList<>();
        List<String> row = new ArrayList<>();
        StringBuilder cell = new StringBuilder();
        boolean quoted = false;
        for (int i = 0; i < text.length(); i++) {
            char ch = text.charAt(i);
            if (quoted) {
                if (ch == '"') {
                    if (i + 1 < text.length() && text.charAt(i + 1) == '"') {
                        cell.append('"');
                        i++;
                    } else quoted = false;
                } else cell.append(ch);
            } else if (ch == '"' && cell.length() == 0) {
                quoted = true;
            } else if (ch == ',') {
                row.add(cell.toString());
                cell.setLength(0);
            } else if (ch == '\n' || ch == '\r') {
                if (ch == '\r' && i + 1 < text.length() && text.charAt(i + 1) == '\n') i++;
                row.add(cell.toString());
                records.add(row);
                row = new ArrayList<>();
                cell.setLength(0);
            } else cell.append(ch);
        }
        if (quoted) throw new IllegalArgumentException("Unclosed quote in temperature CSV.");
        if (cell.length() > 0 || !row.isEmpty()) {
            row.add(cell.toString());
            records.add(row);
        }
        return records;
    }

    private static Metrics metrics(Measurement row) {
        if (!row.complete()) throw new IllegalArgumentException("Metrics require a complete quartet.");
        double[] f = row.frequencies();
        double mean = Arrays.stream(f).sum() / 4;
        double min = Arrays.stream(f).min().orElseThrow();
        double max = Arrays.stream(f).max().orElseThrow();
        double variance = 0;
        for (double value : f) variance += (value - mean) * (value - mean) / 4;
        double rms = Math.sqrt(variance);
        double b = (f[0] - f[3]) / 2;
        double c = (f[1] - f[2]) / 2;
        double d2 = (f[0] - f[1] + f[3] - f[2]) / 4;
        return new Metrics(row, mean, max - min, rms, rms / mean, max / min, b, c, d2);
    }

    private static List<Check> validate(List<Measurement> rows, List<Metrics> metrics) {
        List<Check> checks = new ArrayList<>();
        checks.add(new Check("metrics_count_equals_complete_input_rows",
            metrics.size() - rows.stream().filter(Measurement::complete).count(), 0));
        double parseval = 0, rmse = 0, interpolation = 0, mean = 0, normalized = 0;
        double scale = 1;
        for (Metrics m : metrics) {
            double[] observed = m.measurement().frequencies();
            scale = Math.max(scale, Arrays.stream(observed).max().orElseThrow());
            double varianceFromFourier = (m.b() * m.b() + m.c() * m.c()) / 2 + m.d2() * m.d2();
            parseval = Math.max(parseval, Math.abs(m.rms() * m.rms() - varianceFromFourier));
            double sse = 0;
            double reconstructedMean = 0;
            for (int i = 0; i < 4; i++) {
                double theta = Math.PI * i / 2;
                double predicted = m.mean() + m.b() * Math.cos(theta) + m.c() * Math.sin(theta);
                double actual = observed[CYCLIC_ORDER[i]];
                sse += (predicted - actual) * (predicted - actual);
                interpolation = Math.max(interpolation,
                    Math.abs(predicted + m.d2() * Math.cos(2 * theta) - actual));
                reconstructedMean += predicted / 4;
            }
            mean = Math.max(mean, Math.abs(reconstructedMean - m.mean()));
            rmse = Math.max(rmse, Math.abs(Math.sqrt(sse / 4) - m.rmse()));
            normalized = Math.max(normalized, Math.abs(m.normalizedRms() * m.mean() - m.rms()));
        }
        double tolerance = 2e-12 * scale;
        checks.add(new Check("max_variance_vs_fixed_order_Parseval_MHz2", parseval, tolerance * scale));
        checks.add(new Check("max_direct_first_harmonic_RMSE_vs_abs_d2_MHz", rmse, tolerance));
        checks.add(new Check("max_two_harmonic_reconstruction_error_MHz", interpolation, tolerance));
        checks.add(new Check("max_mean_vs_nu0_MHz", mean, tolerance));
        checks.add(new Check("max_normalized_rms_reconstruction_error_MHz", normalized, tolerance));
        List<List<String>> probe = parseCsv("a,b,c\r\n1,,\"x,y\"\r\n2,3,\"quote\"\"word\"\r\n");
        checks.add(new Check("CSV_empty_frequency_cell_remains_empty", probe.get(1).get(1).isEmpty() ? 0 : 1, 0));
        checks.add(new Check("CSV_quoted_source_fields_are_preserved",
            probe.get(1).get(2).equals("x,y") && probe.get(2).get(2).equals("quote\"word") ? 0 : 1, 0));
        return checks;
    }

    private static void writeMetrics(Path path, List<Metrics> metrics) throws IOException {
        try (BufferedWriter w = Files.newBufferedWriter(path, StandardCharsets.UTF_8)) {
            csv(w, "T_K", "source_row", "s1_MHz", "s2_MHz", "s3_MHz", "s4_MHz",
                "mean_MHz", "span_MHz", "rms_splitting_MHz", "normalized_rms", "max_to_min_ratio",
                "nu0_MHz", "b_MHz", "c_MHz", "d2_MHz", "first_harmonic_RMSE_MHz",
                "fixed_cyclic_order", "order_optimized_at_this_temperature");
            for (Metrics m : metrics) {
                double[] f = m.measurement().frequencies();
                csv(w, m.measurement().temperature(), m.measurement().sourceRow(), f[0], f[1], f[2], f[3],
                    m.mean(), m.span(), m.rms(), m.normalizedRms(), m.maxMinRatio(), m.mean(), m.b(), m.c(), m.d2(),
                    m.rmse(), "s1;s2;s4;s3", false);
            }
        }
    }

    private static void writeSummary(Path path, Path input, String hash, List<Measurement> rows,
                                     List<Metrics> metrics, double maxGap, List<Check> checks) throws IOException {
        StringBuilder j = new StringBuilder("{\n");
        j.append("  \"schema_version\": 1,\n  \"input_csv\": ").append(quote(input.toString().replace('\\', '/'))).append(",\n");
        j.append("  \"input_sha256\": ").append(quote(hash)).append(",\n");
        j.append("  \"input_rows\": ").append(rows.size()).append(",\n");
        j.append("  \"complete_quartets\": ").append(metrics.size()).append(",\n");
        j.append("  \"incomplete_rows\": ").append(rows.size() - metrics.size()).append(",\n");
        j.append("  \"rows_sorted_by_temperature_for_display\": true,\n");
        j.append("  \"all_recorded_temperature_min_K\": ").append(number(rows.get(0).temperature())).append(",\n");
        j.append("  \"all_recorded_temperature_max_K\": ").append(number(rows.get(rows.size() - 1).temperature())).append(",\n");
        j.append("  \"metrics_use_only_complete_quartets\": true,\n");
        j.append("  \"missing_values_imputed\": false,\n");
        j.append("  \"fixed_cyclic_order\": [\"s1\", \"s2\", \"s4\", \"s3\"],\n");
        j.append("  \"cyclic_order_optimized_at_each_temperature\": false,\n");
        j.append("  \"spatial_order_measured\": false,\n");
        j.append("  \"coefficient_formulas\": {\"nu0\": \"(s1+s2+s3+s4)/4\", \"b\": \"(s1-s4)/2\", "
            + "\"c\": \"(s2-s3)/2\", \"d2\": \"(s1-s2+s4-s3)/4\", \"first_harmonic_RMSE\": \"abs(d2)\"},\n");
        j.append("  \"rms_splitting_definition\": \"sqrt(sum_i((nu_i-mean)^2)/4); population RMS\",\n");
        j.append("  \"normalized_rms_definition\": \"rms_splitting/mean; dimensionless fraction\",\n");
        j.append("  \"incomplete_quartet_temperatures_K\": ")
            .append(temperatureArray(rows.stream().filter(r -> !r.complete()).toList())).append(",\n");
        j.append("  \"branches\": {\n");
        for (int branch = 0; branch < 4; branch++) {
            int index = branch;
            List<Measurement> missing = rows.stream().filter(r -> !Double.isFinite(r.frequencies()[index])).toList();
            j.append("    ").append(quote(LABELS[branch])).append(": {\"observed_points\": ")
                .append(rows.size() - missing.size()).append(", \"missing_temperatures_K\": ")
                .append(temperatureArray(missing)).append('}').append(branch < 3 ? ",\n" : "\n");
        }
        j.append("  },\n  \"first_complete_quartet\": ")
            .append(metrics.isEmpty() ? "null" : metricsJson(metrics.get(0))).append(",\n");
        j.append("  \"last_complete_quartet\": ")
            .append(metrics.isEmpty() ? "null" : metricsJson(metrics.get(metrics.size() - 1))).append(",\n");
        j.append("  \"endpoint_comparison\": ");
        if (metrics.isEmpty()) j.append("null,\n");
        else {
            Metrics first = metrics.get(0), last = metrics.get(metrics.size() - 1);
            j.append("{\n");
            jsonNumber(j, "last_minus_first_mean_MHz", last.mean() - first.mean());
            jsonNumber(j, "last_minus_first_span_MHz", last.span() - first.span());
            jsonNumber(j, "last_minus_first_rms_splitting_MHz", last.rms() - first.rms());
            jsonNumber(j, "last_minus_first_normalized_rms", last.normalizedRms() - first.normalizedRms());
            jsonNumber(j, "last_over_first_span", ratio(last.span(), first.span()));
            jsonNumber(j, "last_over_first_rms_splitting", ratio(last.rms(), first.rms()));
            jsonNumber(j, "last_over_first_normalized_rms", ratio(last.normalizedRms(), first.normalizedRms()));
            j.append("    \"last_minus_first_fixed_order_RMSE_MHz\": ").append(number(last.rmse() - first.rmse())).append("\n  },\n");
        }
        j.append("  \"graph\": {\"all_input_temperatures_retained\": true, \"segments_break_at_missing_cells\": true, "
            + "\"maximum_connection_gap_K\": ").append(number(maxGap))
            .append(", \"curves_are_guides_to_eye_not_fits\": true},\n");
        j.append("  \"interpretation_cautions\": [\n"
            + "    \"The mean is the average of a complete quartet at that temperature, not a predicted high-temperature line.\",\n"
            + "    \"Missing frequencies are not evidence of coalescence; their cause is not inferred.\",\n"
            + "    \"Branch labels are retained from the source; no branch exchange is made near crossings.\",\n"
            + "    \"The fixed cyclic assignment is hypothetical and does not determine a crystal structure or a propagation vector.\",\n"
            + "    \"No statistical significance is assigned without frequency uncertainties and their covariance.\"\n  ],\n");
        j.append("  \"validation\": {\"passed\": ").append(checks.stream().filter(Check::pass).count())
            .append(", \"failed\": ").append(checks.stream().filter(c -> !c.pass()).count()).append("},\n");
        j.append("  \"validation_checks\": [\n");
        for (int i = 0; i < checks.size(); i++) {
            Check c = checks.get(i);
            j.append("    {\"name\": ").append(quote(c.name())).append(", \"absolute_error\": ").append(number(c.error()))
                .append(", \"absolute_tolerance\": ").append(number(c.tolerance())).append(", \"passed\": ").append(c.pass())
                .append('}').append(i < checks.size() - 1 ? ",\n" : "\n");
        }
        j.append("  ]\n}\n");
        Files.writeString(path, j.toString(), StandardCharsets.UTF_8);
    }

    private static String metricsJson(Metrics m) {
        return "{\"T_K\": " + number(m.measurement().temperature()) + ", \"source_row\": " + quote(m.measurement().sourceRow())
            + ", \"mean_MHz\": " + number(m.mean()) + ", \"span_MHz\": " + number(m.span())
            + ", \"rms_splitting_MHz\": " + number(m.rms()) + ", \"normalized_rms\": " + number(m.normalizedRms())
            + ", \"max_to_min_ratio\": " + number(m.maxMinRatio()) + ", \"nu0_MHz\": " + number(m.mean())
            + ", \"b_MHz\": " + number(m.b()) + ", \"c_MHz\": " + number(m.c()) + ", \"d2_MHz\": " + number(m.d2())
            + ", \"first_harmonic_RMSE_MHz\": " + number(m.rmse()) + "}";
    }

    private static void writeFigure(Path path, List<Measurement> rows, int completeCount, double maxGap) throws IOException {
        if (path.getParent() != null) Files.createDirectories(path.getParent());
        double minT = rows.get(0).temperature(), maxT = rows.get(rows.size() - 1).temperature();
        double minF = Double.POSITIVE_INFINITY, maxF = Double.NEGATIVE_INFINITY;
        for (Measurement row : rows) for (double f : row.frequencies()) if (Double.isFinite(f)) {
            minF = Math.min(minF, f);
            maxF = Math.max(maxF, f);
        }
        double tStep = niceStep(Math.max(maxT - minT, 1) / 6);
        double xMin = Math.max(0, Math.floor((minT - 0.02 * (maxT - minT)) / tStep) * tStep);
        double xMax = Math.ceil(maxT / tStep) * tStep;
        if (xMax <= xMin) xMax = xMin + tStep;
        double fSpan = Math.max(maxF - minF, Math.max(maxF * 0.05, 0.1));
        double fStep = niceStep(fSpan / 6);
        double yMin = Math.floor((minF - 0.07 * fSpan) / fStep) * fStep;
        double yMax = Math.ceil((maxF + 0.07 * fSpan) / fStep) * fStep;
        double left = 105, right = 1110, top = 155, bottom = 555;
        StringBuilder s = new StringBuilder("""
            <?xml version="1.0" encoding="UTF-8"?>
            <svg xmlns="http://www.w3.org/2000/svg" width="1200" height="735" viewBox="0 0 1200 735" role="img" aria-labelledby="title desc">
            <title id="title">Archived NQR branches versus temperature</title>
            <desc id="desc">All available input frequencies are plotted using their original branch labels. Lines stop at missing values and large temperature gaps. Missing values are not treated as merged lines.</desc>
            <rect width="1200" height="735" fill="#fff"/>
            <g font-family="DejaVu Sans,Arial,sans-serif" fill="#192433">
            <text x="55" y="47" font-size="27" font-weight="700">Archived NQR branches versus temperature</text>
            """);
        svgText(s, 55, 79, 16, rows.size() + " recorded temperatures · " + completeCount
            + " complete quartets · original branch labels retained", "start", "#475569");
        for (int branch = 0; branch < 4; branch++) {
            double x = 130 + branch * 235;
            svgLine(s, x, 120, x + 35, 120, COLORS[branch], 2.3);
            marker(s, x + 17.5, 120, branch, 4.5);
            svgText(s, x + 49, 125, 16, LABELS[branch], "start", "#192433");
        }
        for (int k = 0; k <= Math.round((yMax - yMin) / fStep); k++) {
            double value = yMin + k * fStep;
            double y = map(value, yMin, yMax, bottom, top);
            svgLine(s, left, y, right, y, "#e2e8f0", 1);
            svgText(s, left - 16, y + 5, 14, compact(value, 3), "end", "#475569");
        }
        for (int k = 0; k <= Math.round((xMax - xMin) / tStep); k++) {
            double value = xMin + k * tStep;
            double x = map(value, xMin, xMax, left, right);
            svgLine(s, x, top, x, bottom, "#e2e8f0", 1);
            svgText(s, x, bottom + 28, 14, compact(value, 2), "middle", "#475569");
        }
        svgLine(s, left, top, left, bottom, "#64748b", 1.3);
        svgLine(s, left, bottom, right, bottom, "#64748b", 1.3);
        for (int branch = 0; branch < 4; branch++) {
            s.append("<g id=\"").append(LABELS[branch]).append("\">\n");
            Measurement previous = null;
            for (Measurement row : rows) {
                double frequency = row.frequencies()[branch];
                if (!Double.isFinite(frequency)) {
                    previous = null;
                    continue;
                }
                if (previous != null && row.temperature() - previous.temperature() <= maxGap) {
                    svgLine(s, map(previous.temperature(), xMin, xMax, left, right),
                        map(previous.frequencies()[branch], yMin, yMax, bottom, top),
                        map(row.temperature(), xMin, xMax, left, right), map(frequency, yMin, yMax, bottom, top),
                        COLORS[branch], 2.2);
                }
                previous = row;
            }
            for (Measurement row : rows) if (Double.isFinite(row.frequencies()[branch])) {
                marker(s, map(row.temperature(), xMin, xMax, left, right),
                    map(row.frequencies()[branch], yMin, yMax, bottom, top), branch, 3.8);
            }
            s.append("</g>\n");
        }
        svgText(s, (left + right) / 2, bottom + 63, 16, "Temperature (K)", "middle", "#192433");
        s.append("<text x=\"31\" y=\"355\" font-size=\"16\" text-anchor=\"middle\" transform=\"rotate(-90 31 355)\">NQR frequency (MHz)</text>\n");
        svgText(s, 55, 658, 15,
            "Markers: input centres; lines: guides to the eye. Missing values are not interpreted as merged lines.", "start", "#475569");
        svgText(s, 55, 691, 14, "Segments break at empty cells or a temperature gap > " + compact(maxGap, 2)
            + " K. Scalar metrics use only complete quartets.", "start", "#475569");
        s.append("</g>\n</svg>\n");
        Files.writeString(path, s.toString(), StandardCharsets.UTF_8);
    }

    private static double niceStep(double raw) {
        double base = Math.pow(10, Math.floor(Math.log10(raw)));
        double normalized = raw / base;
        return (normalized <= 1 ? 1 : normalized <= 2 ? 2 : normalized <= 5 ? 5 : 10) * base;
    }

    private static double map(double value, double min, double max, double from, double to) {
        return from + (value - min) * (to - from) / (max - min);
    }

    private static void marker(StringBuilder s, double x, double y, int branch, double r) {
        String style = " fill=\"#fff\" stroke=\"" + COLORS[branch] + "\" stroke-width=\"1.7\"";
        if (branch == 0) {
            s.append(String.format(Locale.ROOT, "<circle cx=\"%.3f\" cy=\"%.3f\" r=\"%.3f\"%s/>%n", x, y, r, style));
        } else if (branch == 1) {
            s.append(String.format(Locale.ROOT, "<rect x=\"%.3f\" y=\"%.3f\" width=\"%.3f\" height=\"%.3f\"%s/>%n",
                x - r, y - r, 2 * r, 2 * r, style));
        } else if (branch == 2) {
            s.append(String.format(Locale.ROOT, "<path d=\"M %.3f %.3f L %.3f %.3f L %.3f %.3f Z\"%s/>%n",
                x, y - 1.2 * r, x + 1.1 * r, y + r, x - 1.1 * r, y + r, style));
        } else {
            s.append(String.format(Locale.ROOT, "<path d=\"M %.3f %.3f L %.3f %.3f L %.3f %.3f L %.3f %.3f Z\"%s/>%n",
                x, y - 1.3 * r, x + 1.1 * r, y, x, y + 1.3 * r, x - 1.1 * r, y, style));
        }
    }

    private static void svgLine(StringBuilder s, double x1, double y1, double x2, double y2, String color, double width) {
        s.append(String.format(Locale.ROOT,
            "<line x1=\"%.3f\" y1=\"%.3f\" x2=\"%.3f\" y2=\"%.3f\" stroke=\"%s\" stroke-width=\"%.2f\"/>%n",
            x1, y1, x2, y2, color, width));
    }

    private static void svgText(StringBuilder s, double x, double y, int size, String text, String anchor, String fill) {
        s.append(String.format(Locale.ROOT,
            "<text x=\"%.3f\" y=\"%.3f\" font-size=\"%d\" text-anchor=\"%s\" fill=\"%s\">%s</text>%n",
            x, y, size, anchor, fill, text.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")));
    }

    private static void csv(BufferedWriter w, Object... fields) throws IOException {
        for (int i = 0; i < fields.length; i++) {
            if (i > 0) w.write(',');
            String value = fields[i] instanceof Double d ? number(d) : fields[i].toString();
            if (value.contains(",") || value.contains("\"") || value.contains("\n") || value.contains("\r")) {
                w.write('"');
                w.write(value.replace("\"", "\"\""));
                w.write('"');
            } else w.write(value);
        }
        w.write('\n');
    }

    private static String temperatureArray(List<Measurement> rows) {
        StringBuilder result = new StringBuilder("[");
        for (int i = 0; i < rows.size(); i++) {
            if (i > 0) result.append(", ");
            result.append(number(rows.get(i).temperature()));
        }
        return result.append(']').toString();
    }

    private static void jsonNumber(StringBuilder j, String key, double value) {
        j.append("    ").append(quote(key)).append(": ").append(number(value)).append(",\n");
    }

    private static String number(double value) {
        if (!Double.isFinite(value)) return "null";
        if (value == 0) return "0";
        return String.format(Locale.ROOT, "%.15g", value);
    }

    private static String compact(double value, int decimals) {
        String result = String.format(Locale.ROOT, "%1$." + decimals + "f", value);
        return result.indexOf('.') < 0 ? result : result.replaceAll("0+$", "").replaceAll("\\.$", "");
    }

    private static double ratio(double numerator, double denominator) {
        return denominator == 0 ? Double.NaN : numerator / denominator;
    }

    private static String quote(String text) {
        return '"' + text.replace("\\", "\\\\").replace("\"", "\\\"")
            .replace("\n", "\\n").replace("\r", "\\r").replace("\t", "\\t") + '"';
    }

    private static String sha256(byte[] bytes) {
        try {
            return HexFormat.of().formatHex(MessageDigest.getInstance("SHA-256").digest(bytes));
        } catch (NoSuchAlgorithmException e) {
            throw new IllegalStateException("SHA-256 is required by the Java platform.", e);
        }
    }
}
