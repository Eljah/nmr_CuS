package org.eljah.nmr.review;

import org.knowm.xchart.BitmapEncoder;
import org.knowm.xchart.VectorGraphicsEncoder;
import org.knowm.xchart.XYChart;
import org.knowm.xchart.XYChartBuilder;
import org.knowm.xchart.XYSeries;
import org.knowm.xchart.style.Styler;
import org.knowm.xchart.style.lines.SeriesLines;
import org.knowm.xchart.style.markers.SeriesMarkers;

import java.awt.Color;
import java.awt.Font;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.security.MessageDigest;
import java.util.ArrayList;
import java.util.HexFormat;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;

/** Publication figures from committed numerical tables. No model fitting is done here. */
public final class ReviewFigures {
    private static final String[] SITES = {"s1", "s2", "s3", "s4", "CuS"};
    private static final Color[] COLORS = {
        new Color(0x2166ac), new Color(0xd66a00), new Color(0x16836b),
        new Color(0x8b4caf), new Color(0x5e6975)
    };
    private final Path repository;
    private final Path output;
    private final List<Path> sources = new ArrayList<>();

    private ReviewFigures(Path repository) throws IOException {
        this.repository = repository;
        this.output = repository.resolve("theory/combined_review/figures");
        Files.createDirectories(output);
    }

    public static void main(String[] args) throws Exception {
        Locale.setDefault(Locale.ROOT);
        System.setProperty("java.awt.headless", "true");
        Path root = args.length > 0 ? Path.of(args[0]).toAbsolutePath().normalize() : findRepository();
        ReviewFigures plots = new ReviewFigures(root);
        plots.nqrTemperature();
        Path tensor = root.resolve("theory/joint_temperature/results/tensor_coefficients.csv");
        if (Files.isRegularFile(tensor)) plots.tensorDirection();
        Path relaxation = root.resolve("theory/relaxation_dynamics/outputs/r1_normalized.csv");
        if (Files.isRegularFile(relaxation)) {
            plots.relaxation(false);
            plots.relaxation(true);
        }
        if (Files.isRegularFile(root.resolve("theory/relaxation_dynamics/outputs/bpp_tau0_profile.csv")))
            plots.bppProfile();
        if (Files.isRegularFile(root.resolve("theory/nmr_shift_anisotropy/results/final_free_eta_axial_k_500ms.csv"))) {
            plots.nmr(500);
            plots.nmr(20);
        }
        plots.writeProvenance();
        System.out.println("Figures written to " + plots.output);
    }

    private static Path findRepository() {
        for (Path p = Path.of("").toAbsolutePath(); p != null; p = p.getParent()) {
            if (Files.isDirectory(p.resolve("theory/superstructure_hypothesis"))) return p;
        }
        throw new IllegalArgumentException("Run inside the repository, or pass its root path.");
    }

    private XYChart chart(String title, String x, String y) {
        XYChart chart = new XYChartBuilder().width(1180).height(720)
                .title(title).xAxisTitle(x).yAxisTitle(y).build();
        var s = chart.getStyler();
        s.setChartBackgroundColor(Color.WHITE);
        s.setPlotBackgroundColor(Color.WHITE);
        s.setPlotBorderVisible(false);
        s.setLegendPosition(Styler.LegendPosition.OutsideE);
        s.setLegendBackgroundColor(Color.WHITE);
        s.setLegendBorderColor(Color.WHITE);
        s.setChartTitleFont(new Font(Font.SANS_SERIF, Font.BOLD, 22));
        s.setAxisTitleFont(new Font(Font.SANS_SERIF, Font.PLAIN, 19));
        s.setAxisTickLabelsFont(new Font(Font.SANS_SERIF, Font.PLAIN, 16));
        s.setLegendFont(new Font(Font.SANS_SERIF, Font.PLAIN, 17));
        s.setMarkerSize(7);
        s.setSeriesColors(COLORS);
        s.setPlotGridLinesColor(new Color(0xe2e6ea));
        s.setXAxisDecimalPattern("0");
        return chart;
    }

    private void nqrTemperature() throws IOException {
        Path input = repository.resolve("theory/superstructure_hypothesis/inputs/temperature_series.csv");
        List<Map<String,String>> rows = read(input);
        XYChart chart = chart("Четыре архивные ветви ЯКР ⁶³Cu", "Температура, K", "Частота, МГц");
        chart.getStyler().setXAxisMin(0.0);
        chart.getStyler().setXAxisMax(300.0);
        chart.getStyler().setYAxisMin(11.0);
        chart.getStyler().setYAxisMax(17.4);
        chart.getStyler().setYAxisDecimalPattern("0.0");
        for (int i = 0; i < 4; i++) {
            List<Double> x = new ArrayList<>(), y = new ArrayList<>();
            for (Map<String,String> row : rows) {
                double value = number(row, SITES[i] + "_MHz");
                if (Double.isFinite(value)) {x.add(number(row,"T_K")); y.add(value);}
            }
            XYSeries series = chart.addSeries(SITES[i], x, y);
            series.setXYSeriesRenderStyle(XYSeries.XYSeriesRenderStyle.Scatter);
            series.setMarker(SeriesMarkers.CIRCLE);
            series.setMarkerColor(COLORS[i]);
        }
        Path predictions = repository.resolve("theory/joint_temperature/results/thermal_predictions.csv");
        if (Files.isRegularFile(predictions)) {
            List<Map<String,String>> fitted = read(predictions);
            for (int i = 0; i < 4; i++) {
                List<Double> x = new ArrayList<>(), y = new ArrayList<>();
                for (Map<String,String> row : fitted) {
                    if ("original".equals(row.get("labels")) && "bose_common_theta".equals(row.get("model"))
                            && "all".equals(row.get("fold")) && SITES[i].equals(row.get("site"))) {
                        x.add(number(row,"T_K")); y.add(number(row,"fitted_MHz"));
                    }
                }
                XYSeries series = chart.addSeries("Общий тепловой ход " + SITES[i], x, y);
                series.setMarker(SeriesMarkers.NONE);
                series.setLineColor(COLORS[i]);
                series.setShowInLegend(false);
            }
        }
        save(chart, "nqr_temperature");
    }

    private void tensorDirection() throws IOException {
        List<Map<String,String>> rows = read(repository.resolve("theory/joint_temperature/results/tensor_coefficients.csv"));
        XYChart chart = chart("Меняется ли рисунок первой тензорной гармоники?", "Температура, K", "Изменение угла пары (B, F), градусы");
        String[] classes = {"opposite_14_23", "opposite_13_24", "opposite_12_34"};
        String[] labels = {"Пары 1–4 / 2–3", "Пары 1–3 / 2–4", "Пары 1–2 / 3–4"};
        chart.getStyler().setXAxisMin(0.0);
        chart.getStyler().setXAxisMax(280.0);
        chart.getStyler().setYAxisDecimalPattern("0");
        for (int i=0; i<classes.length; i++) {
            List<Double> x = new ArrayList<>(), y = new ArrayList<>();
            double reference = Double.NaN;
            for (Map<String,String> row : rows) {
                if (!"original".equals(row.get("labels")) || !classes[i].equals(row.get("order"))) continue;
                double angle = Math.atan2(number(row,"F_MHz2"),number(row,"B_MHz2"));
                if (!Double.isFinite(reference)) reference = angle;
                x.add(number(row,"T_K"));
                y.add(Math.toDegrees(Math.IEEEremainder(angle-reference,2*Math.PI)));
            }
            if (x.isEmpty()) throw new IOException("Missing tensor class " + classes[i]);
            XYSeries series = chart.addSeries(labels[i],x,y);
            series.setMarker(SeriesMarkers.CIRCLE);
            series.setLineColor(COLORS[i]);
            series.setMarkerColor(COLORS[i]);
        }
        XYSeries zero=chart.addSeries("Фиксированная форма",new double[]{0,280},new double[]{0,0});
        zero.setMarker(SeriesMarkers.NONE);
        zero.setLineColor(Color.GRAY);
        zero.setLineStyle(SeriesLines.DASH_DASH);
        save(chart,"tensor_harmonic_direction");
    }

    private void relaxation(boolean lowTemperature) throws IOException {
        List<Map<String,String>> rows = read(repository.resolve("theory/relaxation_dynamics/outputs/r1_normalized.csv"));
        XYChart chart = chart(lowTemperature ? "Низкотемпературная релаксация ⁶³Cu"
                : "Температурная зависимость R₁/T: четыре группы и CuS", "Температура, K", "R₁/T, с⁻¹ K⁻¹");
        chart.getStyler().setXAxisMin(0.0);
        chart.getStyler().setXAxisMax(lowTemperature ? 80.0 : 305.0);
        chart.getStyler().setYAxisLogarithmic(!lowTemperature);
        chart.getStyler().setYAxisDecimalPattern(lowTemperature ? "0.0" : "0.##");
        for (int i = 0; i < SITES.length; i++) {
            List<Double> x = new ArrayList<>(), y = new ArrayList<>(), error = new ArrayList<>();
            for (Map<String,String> row : rows) {
                String site = normalizeSite(row.get("site"));
                double t = number(row, "T_K");
                if (!SITES[i].equals(site) || (lowTemperature && t > 80)) continue;
                double r = number(row, "R1_s-1"), e = number(row, "error_scale_s-1");
                if (!Double.isFinite(r) || !(t > 0)) continue;
                x.add(t); y.add(r / t); error.add(e / t);
            }
            if (x.isEmpty()) continue;
            XYSeries series = chart.addSeries(SITES[i].equals("CuS") ? "CuS (контроль)" : SITES[i], x, y, error);
            series.setMarker(i == 4 ? SeriesMarkers.DIAMOND : SeriesMarkers.CIRCLE);
            series.setMarkerColor(COLORS[i]);
            series.setLineColor(COLORS[i]);
            series.setLineStyle(i == 4 ? SeriesLines.DASH_DASH : SeriesLines.SOLID);
        }
        save(chart, lowTemperature ? "r1_over_t_low_temperature" : "r1_over_t");
    }

    private static String normalizeSite(String value) {
        if (value == null) return "";
        String s = value.toLowerCase(Locale.ROOT);
        return s.equals("c") || s.equals("cus") || s.contains("covellite") || s.contains("ковел") ? "CuS" : s;
    }

    private void bppProfile() throws IOException {
        List<Map<String,String>> rows = read(repository.resolve("theory/relaxation_dynamics/outputs/bpp_tau0_profile.csv"));
        XYChart chart = chart("Плато параметра τ₀ в однолоренцевой модели",
                "log₁₀(τ₀ / с)", "Критерий S / минимальное S");
        chart.getStyler().setXAxisMin(-18.0);
        chart.getStyler().setXAxisMax(-5.0);
        chart.getStyler().setYAxisMin(0.98);
        chart.getStyler().setYAxisDecimalPattern("0.0");
        String[] weightings = {"stored_error_scale", "relative_floor_5pct"};
        String[] names = {"Архивные Err", "Нижний порог 5%"};
        for (int i = 0; i < weightings.length; i++) {
            List<Double> x = new ArrayList<>(), y = new ArrayList<>();
            for (Map<String,String> row : rows) {
                if (!weightings[i].equals(row.get("weighting"))) continue;
                x.add(number(row, "log10_tau0_s"));
                y.add(number(row, "objective_over_min"));
            }
            XYSeries series = chart.addSeries(names[i], x, y);
            series.setLineColor(COLORS[i]);
            series.setMarkerColor(COLORS[i]);
            series.setMarker(SeriesMarkers.CIRCLE);
        }
        XYSeries threshold = chart.addSeries("Описательный допуск 1%",
                new double[]{-18, -5}, new double[]{1.01, 1.01});
        threshold.setMarker(SeriesMarkers.NONE);
        threshold.setLineColor(Color.GRAY);
        threshold.setLineStyle(SeriesLines.DASH_DASH);
        save(chart, "bpp_parameter_profile");
    }

    private void nmr(int repetitionMs) throws IOException {
        Path folder = repository.resolve("theory/nmr_shift_anisotropy/results");
        String[] models = {"free_eta_isotropic_k", "axial_eta_axial_k", "free_eta_axial_k"};
        String[] names = {"A: η свободна, K изотр.",
                "B: η = 0, K аксиален", "C: η свободна, K аксиален"};
        XYChart chart = chart("Порошковый ЯМР: период повторения " + repetitionMs + " мс",
                "Частота, МГц", "Нормированная интенсивность");
        chart.getStyler().setXAxisMin(72.0);
        chart.getStyler().setXAxisMax(85.0);
        chart.getStyler().setYAxisDecimalPattern("0.0");
        chart.getStyler().setMarkerSize(4);
        XYChart residuals = chart("Остатки ЯМР: период повторения " + repetitionMs + " мс",
                "Частота, МГц", "Модель − измерение, норм. ед.");
        residuals.getStyler().setXAxisMin(72.0);
        residuals.getStyler().setXAxisMax(85.0);
        residuals.getStyler().setYAxisDecimalPattern("0.00");
        for (int i = 0; i < models.length; i++) {
            List<Map<String,String>> rows = read(folder.resolve(
                    "final_" + models[i] + "_" + repetitionMs + "ms.csv"));
            double scale = rows.stream().filter(r -> Boolean.parseBoolean(r.get("used_in_fit")))
                    .mapToDouble(r -> number(r, "observed")).max().orElseThrow();
            if (!(scale > 0)) throw new IOException("Invalid NMR normalization");
            if (i == 0) {
                List<Double> usedX = new ArrayList<>(), usedY = new ArrayList<>();
                List<Double> excludedX = new ArrayList<>(), excludedY = new ArrayList<>();
                for (Map<String,String> row : rows) {
                    boolean used = Boolean.parseBoolean(row.get("used_in_fit"));
                    (used ? usedX : excludedX).add(number(row, "frequency_MHz"));
                    (used ? usedY : excludedY).add(number(row, "observed") / scale);
                }
                XYSeries observed = chart.addSeries("Измерение", usedX, usedY);
                observed.setXYSeriesRenderStyle(XYSeries.XYSeriesRenderStyle.Scatter);
                observed.setMarker(SeriesMarkers.CIRCLE);
                observed.setMarkerColor(new Color(0x42464a));
                if (!excludedX.isEmpty()) {
                    XYSeries excluded = chart.addSeries("Исключено маской", excludedX, excludedY);
                    excluded.setXYSeriesRenderStyle(XYSeries.XYSeriesRenderStyle.Scatter);
                    excluded.setMarker(SeriesMarkers.CROSS);
                    excluded.setMarkerColor(new Color(0xaaaaaa));
                }
            }
            List<Double> x = new ArrayList<>(), y = new ArrayList<>();
            List<Double> rx = new ArrayList<>(), ry = new ArrayList<>();
            for (Map<String,String> row : rows) {
                x.add(number(row, "frequency_MHz"));
                y.add(number(row, "predicted") / scale);
                if (Boolean.parseBoolean(row.get("used_in_fit"))) {
                    rx.add(number(row, "frequency_MHz"));
                    ry.add(number(row, "residual") / scale);
                }
            }
            XYSeries series = chart.addSeries(names[i], x, y);
            series.setMarker(SeriesMarkers.NONE);
            series.setLineColor(COLORS[i]);
            series.setLineWidth(2.0f);
            XYSeries residual = residuals.addSeries(names[i], rx, ry);
            residual.setMarker(SeriesMarkers.NONE);
            residual.setLineColor(COLORS[i]);
        }
        XYSeries zero = residuals.addSeries("Нулевой остаток", new double[]{72,85}, new double[]{0,0});
        zero.setMarker(SeriesMarkers.NONE);
        zero.setLineColor(Color.GRAY);
        zero.setLineStyle(SeriesLines.DASH_DASH);
        save(chart, "nmr_comparison_fullscale_" + repetitionMs + "ms");
        chart.getStyler().setYAxisMin(-0.02);
        chart.getStyler().setYAxisMax(1.06);
        save(chart, "nmr_comparison_" + repetitionMs + "ms");
        save(residuals, "nmr_residuals_" + repetitionMs + "ms");
    }

    private List<Map<String,String>> read(Path path) throws IOException {
        if (!sources.contains(path)) sources.add(path);
        List<String> lines = Files.readAllLines(path, StandardCharsets.UTF_8);
        if (lines.isEmpty()) throw new IOException("Empty CSV: " + path);
        List<String> header = csvLine(lines.get(0));
        List<Map<String,String>> rows = new ArrayList<>();
        for (int i = 1; i < lines.size(); i++) {
            if (lines.get(i).isBlank()) continue;
            List<String> fields = csvLine(lines.get(i));
            if (fields.size() != header.size()) throw new IOException("CSV column count at " + path + ":" + (i + 1));
            Map<String,String> row = new LinkedHashMap<>();
            for (int j = 0; j < header.size(); j++) row.put(header.get(j), fields.get(j));
            rows.add(row);
        }
        return rows;
    }

    private static List<String> csvLine(String line) {
        List<String> fields = new ArrayList<>();
        StringBuilder value = new StringBuilder();
        boolean quoted = false;
        for (int i = 0; i < line.length(); i++) {
            char c = line.charAt(i);
            if (c == '"') {
                if (quoted && i + 1 < line.length() && line.charAt(i + 1) == '"') {value.append(c); i++;}
                else quoted = !quoted;
            } else if (c == ',' && !quoted) {fields.add(value.toString()); value.setLength(0);}
            else value.append(c);
        }
        if (quoted) throw new IllegalArgumentException("Unterminated CSV quote");
        fields.add(value.toString());
        return fields;
    }

    private static double number(Map<String,String> row, String name) {
        String v = row.get(name);
        return v == null || v.isBlank() ? Double.NaN : Double.parseDouble(v);
    }

    private void save(XYChart chart, String name) throws IOException {
        String base = output.resolve(name).toString();
        BitmapEncoder.saveBitmap(chart, base, BitmapEncoder.BitmapFormat.PNG);
        VectorGraphicsEncoder.saveVectorGraphic(chart, base, VectorGraphicsEncoder.VectorGraphicsFormat.SVG);
    }

    private void writeProvenance() throws Exception {
        StringBuilder s = new StringBuilder("{\n  \"generator\": \"Java 17 / XChart 3.8.8\",\n  \"fit_performed\": false,\n  \"source_files\": [\n");
        for (int i = 0; i < sources.size(); i++) {
            Path p = sources.get(i);
            String hash = HexFormat.of().formatHex(MessageDigest.getInstance("SHA-256").digest(Files.readAllBytes(p)));
            s.append("    {\"path\": \"").append(repository.relativize(p).toString().replace('\\','/'))
                    .append("\", \"sha256\": \"").append(hash).append("\"}")
                    .append(i + 1 == sources.size() ? "\n" : ",\n");
        }
        s.append("  ],\n  \"notes\": [\n")
                .append("    \"NQR points are archived measurements; optional lines are the shared-Bose fit, not additional measurements.\",\n")
                .append("    \"The tensor-direction angle is in the (B,F) Fourier-coefficient plane; it is not an atomic or EFG-axis rotation angle.\",\n")
                .append("    \"R1 error bars preserve archived error scales; their confidence convention is unknown.\",\n")
                .append("    \"Lines joining R1 points are visual guides, not fitted models.\",\n")
                .append("    \"BPP one-percent tolerance is descriptive, not a confidence interval.\",\n")
                .append("    \"NMR data and predictions use the maximum retained observed intensity of each repetition dataset; masked points remain marked.\",\n")
                .append("    \"NMR comparison panels zoom to normalized intensity -0.02 to 1.06; fullscale panels retain the taller masked instrumental peaks.\"\n  ]\n}\n");
        Files.writeString(output.resolve("PROVENANCE.json"), s, StandardCharsets.UTF_8);
    }
}
