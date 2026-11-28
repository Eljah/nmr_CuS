package org.eljah.nmr;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;

final class NmrFiles {
    private NmrFiles() { }
    record Reference(String label, PowderModel.Grid grid, double expected, double[] p) { }
    static List<Reference> references(Path path) throws IOException {
        List<Reference> result = new ArrayList<>();
        List<String> lines = Files.readAllLines(path);
        for (int line = 1; line < lines.size(); line++) {
            if (lines.get(line).isBlank()) continue;
            String[] c = lines.get(line).split(",");
            double[] p = new double[15];
            for (int i = 0; i < 11; i++) p[i] = Double.parseDouble(c[i + 4]);
            result.add(new Reference(c[0], new PowderModel.Grid(Integer.parseInt(c[1]), Integer.parseInt(c[2])),
                Double.parseDouble(c[3]), p));
        }
        return result;
    }
    static BufferedWriter writer(Path path) throws IOException {
        if (path.getParent() != null) Files.createDirectories(path.getParent());
        return Files.newBufferedWriter(path, StandardCharsets.UTF_8);
    }
    static void csv(BufferedWriter w, Object... fields) throws IOException {
        for (int i = 0; i < fields.length; i++) {
            if (i > 0) w.write(',');
            String value = fields[i] instanceof Double d ? number(d) : fields[i].toString();
            if (value.contains(",") || value.contains("\"") || value.contains("\n")) {
                w.write('"'); w.write(value.replace("\"", "\"\"")); w.write('"');
            } else w.write(value);
        }
        w.write('\n');
    }
    static String number(double value) {
        if (!Double.isFinite(value)) return "null";
        if (value == 0) return "0";
        return String.format(Locale.ROOT, "%.16g", value);
    }
    static String quote(String text) {
        return '"' + text.replace("\\", "\\\\").replace("\"", "\\\"").replace("\n", "\\n").replace("\r", "\\r") + '"';
    }
    static String array(double[] values) {
        StringBuilder s = new StringBuilder("[");
        for (int i = 0; i < values.length; i++) { if (i > 0) s.append(", "); s.append(number(values[i])); }
        return s.append(']').toString();
    }
    static void spectrum(Path path, PowderModel.Data d, PowderModel.SeriesResult r) throws IOException {
        try (BufferedWriter w = writer(path)) {
            csv(w,"source_row","frequency_MHz","observed","predicted","residual","used_in_fit","predicted_normalized",
                "linear_baseline","s1_component","s2_component","s3_component","s4_component");
            for (int i = 0; i < d.x.length; i++) {
                double pred = r.predictedNormalized()[i] * d.scale;
                csv(w,d.sourceRow[i],d.x[i],d.observed[i],pred,pred - d.observed[i],d.used[i],r.predictedNormalized()[i],
                    (r.baselineNormalized()[0] + r.baselineNormalized()[1] * d.t[i]) * d.scale,
                    r.design()[i][0] * r.amplitudesNormalized()[0] * d.scale,
                    r.design()[i][1] * r.amplitudesNormalized()[1] * d.scale,
                    r.design()[i][2] * r.amplitudesNormalized()[2] * d.scale,
                    r.design()[i][3] * r.amplitudesNormalized()[3] * d.scale);
            }
        }
    }
}
