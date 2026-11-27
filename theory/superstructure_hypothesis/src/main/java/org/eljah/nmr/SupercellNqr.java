package org.eljah.nmr;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Locale;

/**
 * Four labelled frequencies on the four phases theta_m = pi*m/2.
 *
 * This is a scalar-frequency interpolation, not a crystal-structure solution.
 * The input labels contain no information about the spatial cyclic order.
 * No uncertainty, statistical rejection threshold, atomic displacement or
 * temperature dependence is inferred from the supplied line centres.
 *
 * Convention: f_m = nu0 + b*cos(theta_m) + c*sin(theta_m) + d2*cos(2*theta_m).
 * The first harmonic is A1*cos(theta_m + phi1), so phi1 = atan2(-c, b).
 * On this four-point grid sin(2*theta_m) vanishes and its coefficient cannot
 * be identified. The second-harmonic interpolation has zero residual degrees
 * of freedom. A harmonic EFG tensor need not produce a harmonic frequency.
 */
public final class SupercellNqr {
    private static final double[] DEFAULT_FREQUENCIES = {11.97, 15.27, 16.18, 16.96};
    private static final double[][] BASIS = {
        {1, 1, 1, 1}, {1, 0, -1, 0}, {0, 1, 0, -1}, {1, -1, 1, -1}
    };
    private static final int[][][] PAIRS = {
        {{0, 1}, {2, 3}}, {{0, 2}, {1, 3}}, {{0, 3}, {1, 2}}
    };
    private static final String ORDER_NOTE =
        "Spatial order is not measured; all 24 labelled assignments are enumerated.";
    private static final String INTERPOLATION_NOTE =
        "Four coefficients interpolate four samples exactly; zero residual degrees of freedom."
        + " No structure, propagation vector or charge-density-wave mechanism is identified.";

    private SupercellNqr() { }

    private record Fit(int id, int[] labels, double[] observed, int pairId,
                       double nu0, double b, double c, double d2,
                       double amplitude1, double phase1,
                       double[] first, double[] exact, double rmse, double exactRmse) { }
    private record Check(String name, double value, double expected, double tolerance) {
        boolean pass() {
            return Double.isFinite(value) && Math.abs(value - expected) <= tolerance;
        }
    }

    public static void main(String[] args) throws IOException {
        double[] frequencies = DEFAULT_FREQUENCIES.clone();
        Path output = Path.of("results");
        Path figure = Path.of("figures/four_subcells.svg");
        for (String arg : args) {
            if (arg.equals("--help") || arg.equals("-h")) {
                System.out.println("Usage: java -jar superstructure-hypothesis-1.0.0.jar "
                    + "[--frequencies=11.97,15.27,16.18,16.96] [--output=results] "
                    + "[--figure=figures/four_subcells.svg]");
                System.out.println("Four positive frequencies in MHz, in input-label order. "
                    + "All 24 spatial assignments are evaluated with equal weights.");
                System.out.println("Phase convention: A1*cos(theta+phi1). "
                    + "No physical spatial order or experimental errors are assumed.");
                return;
            } else if (arg.startsWith("--frequencies=")) {
                String[] parts = arg.substring("--frequencies=".length()).split(",", -1);
                if (parts.length != 4) {
                    throw new IllegalArgumentException("Exactly four frequencies are required.");
                }
                frequencies = new double[4];
                for (int i = 0; i < 4; i++) {
                    frequencies[i] = Double.parseDouble(parts[i].trim());
                }
            } else if (arg.startsWith("--output=")) {
                output = Path.of(arg.substring("--output=".length()));
            } else if (arg.startsWith("--figure=")) {
                figure = Path.of(arg.substring("--figure=".length()));
            } else {
                throw new IllegalArgumentException("Unknown argument: " + arg + "; use --help.");
            }
        }
        for (double f : frequencies) {
            if (!Double.isFinite(f) || f <= 0 || f > 1e100) {
                throw new IllegalArgumentException("Frequencies must be finite and positive in MHz.");
            }
        }

        List<int[]> permutations = new ArrayList<>();
        enumerate(new int[4], new boolean[4], 0, permutations);
        List<Fit> fits = new ArrayList<>();
        for (int[] labels : permutations) {
            fits.add(fit(fits.size() + 1, labels, frequencies));
        }
        double scale = Arrays.stream(frequencies).max().orElseThrow();
        double tolerance = 2e-12 * Math.max(1.0, scale);
        Fit best = fits.get(0);
        for (Fit f : fits) {
            // Resolve numerical ties by the lexicographic input-label order.
            if (f.rmse() < best.rmse() - tolerance) best = f;
        }
        List<Check> checks = validate(frequencies, fits, best, tolerance);
        Files.createDirectories(output);
        writeOrders(output.resolve("all_orders.csv"), fits);
        writePairs(output.resolve("opposite_pairs.csv"), frequencies, fits);
        writeBest(output.resolve("best_order_samples.csv"), best);
        writeInvariants(output.resolve("tensor_invariants.csv"), frequencies);
        writeChecks(output.resolve("validation_checks.csv"), checks);
        writeSummary(output.resolve("summary.json"), frequencies, fits, best, checks, tolerance);
        writeFigure(figure, best);

        long failures = checks.stream().filter(c -> !c.pass()).count();
        double min = Arrays.stream(frequencies).min().orElseThrow();
        System.out.printf(Locale.ROOT,
            "24 labelled orders; 3 opposite-pair classes. Best order (MHz): %s%n"
            + "nu0=%.12f, b=%.12f, c=%.12f, d2=%.12f MHz%n"
            + "A1=%.12f MHz, phi1=%.12f degrees (cos(theta+phi1))%n"
            + "First-harmonic RMSE=%.12f MHz; exact-interpolation RMSE=%.3g MHz%n"
            + "Observed max/min=%.12f; eta-only bound=%.12f; minimum |Vzz| ratio=%.12f%n"
            + "Consistency checks: %d PASS, %d FAIL. Results: %s%n",
            Arrays.toString(best.observed()), best.nu0(), best.b(), best.c(), best.d2(),
            best.amplitude1(), Math.toDegrees(best.phase1()), best.rmse(), best.exactRmse(),
            scale / min, Math.sqrt(4.0 / 3.0), Math.max(1, scale / min / Math.sqrt(4.0 / 3.0)),
            checks.size() - failures, failures, output.toAbsolutePath().normalize());
        System.out.println(ORDER_NOTE);
        System.out.println(INTERPOLATION_NOTE);
        System.out.println("Figure: " + figure.toAbsolutePath().normalize());
        if (failures != 0) throw new IllegalStateException(failures + " consistency checks failed.");
    }

    private static void enumerate(int[] order, boolean[] used, int n, List<int[]> result) {
        if (n == order.length) {
            result.add(order.clone());
            return;
        }
        for (int label = 0; label < order.length; label++) {
            if (!used[label]) {
                used[label] = true;
                order[n] = label;
                enumerate(order, used, n + 1, result);
                used[label] = false;
            }
        }
    }

    private static Fit fit(int id, int[] labels, double[] frequencies) {
        double[] observed = new double[4];
        for (int i = 0; i < 4; i++) observed[i] = frequencies[labels[i]];
        // Real DFT of a length-four sequence, with the Nyquist term undoubled.
        double nu0 = (observed[0] + observed[1] + observed[2] + observed[3]) / 4.0;
        double b = (observed[0] - observed[2]) / 2.0;
        double c = (observed[1] - observed[3]) / 2.0;
        double d2 = (observed[0] - observed[1] + observed[2] - observed[3]) / 4.0;
        double[] first = new double[4];
        double[] exact = new double[4];
        for (int m = 0; m < 4; m++) {
            first[m] = nu0 + b * BASIS[1][m] + c * BASIS[2][m];
            exact[m] = first[m] + d2 * BASIS[3][m];
        }
        int zeroPosition = 0;
        while (labels[zeroPosition] != 0) zeroPosition++;
        int oppositeLabel = labels[(zeroPosition + 2) % 4];
        return new Fit(id, labels.clone(), observed, oppositeLabel,
            nu0, b, c, d2, Math.hypot(b, c), Math.atan2(-c, b),
            first, exact, rmse(first, observed), rmse(exact, observed));
    }

    private static double rmse(double[] prediction, double[] observed) {
        double sum = 0;
        for (int i = 0; i < prediction.length; i++) sum += square(prediction[i] - observed[i]);
        return Math.sqrt(sum / prediction.length);
    }

    private static double square(double x) { return x * x; }
    private static double dot(double[] a, double[] b) {
        double result = 0;
        for (int i = 0; i < a.length; i++) result += a[i] * b[i];
        return result;
    }

    private static List<Check> validate(double[] input, List<Fit> fits, Fit best, double tol) {
        List<Check> checks = new ArrayList<>();
        checks.add(new Check("number_of_labelled_permutations", fits.size(), 24, 0));
        HashSet<String> unique = new HashSet<>();
        for (Fit f : fits) unique.add(Arrays.toString(f.labels()));
        checks.add(new Check("number_of_unique_label_orders", unique.size(), 24, 0));
        for (int pair = 1; pair <= 3; pair++) {
            int current = pair;
            long count = fits.stream().filter(f -> f.pairId() == current).count();
            checks.add(new Check("orders_in_opposite_pair_class_" + pair, count, 8, 0));
        }
        double projectionError = 0, exactError = 0, normalError = 0;
        double rmseError = 0, phaseError = 0, classError = 0;
        for (Fit f : fits) {
            double[] coefficients = {f.nu0(), f.b(), f.c(), f.d2()};
            for (int k = 0; k < 4; k++) {
                double independentProjection = dot(f.observed(), BASIS[k]) / dot(BASIS[k], BASIS[k]);
                projectionError = Math.max(projectionError, Math.abs(coefficients[k] - independentProjection));
            }
            double[] residual = new double[4];
            for (int m = 0; m < 4; m++) {
                residual[m] = f.first()[m] - f.observed()[m];
                exactError = Math.max(exactError, Math.abs(f.exact()[m] - f.observed()[m]));
                double byPhase = f.nu0() + f.amplitude1() * Math.cos(Math.PI * m / 2 + f.phase1());
                phaseError = Math.max(phaseError, Math.abs(byPhase - f.first()[m]));
            }
            for (int k = 0; k < 3; k++) normalError = Math.max(normalError, Math.abs(dot(residual, BASIS[k])));
            rmseError = Math.max(rmseError, Math.abs(f.rmse() - Math.abs(f.d2())));
            int[][] pair = PAIRS[f.pairId() - 1];
            double pairRmse = Math.abs(input[pair[0][0]] + input[pair[0][1]]
                - input[pair[1][0]] - input[pair[1][1]]) / 4;
            classError = Math.max(classError, Math.abs(f.rmse() - pairRmse));
        }
        checks.add(new Check("max_closed_form_vs_orthogonal_DFT_projection_MHz", projectionError, 0, tol));
        checks.add(new Check("max_exact_interpolation_error_MHz", exactError, 0, tol));
        checks.add(new Check("max_first_harmonic_normal_equation_residual_MHz", normalError, 0, 4 * tol));
        checks.add(new Check("max_RMSE_minus_abs_d2_MHz", rmseError, 0, tol));
        checks.add(new Check("max_amplitude_phase_reconstruction_error_MHz", phaseError, 0, tol));
        checks.add(new Check("max_order_RMSE_vs_opposite_pair_contrast_MHz", classError, 0, tol));
        double pairMinimum = Double.POSITIVE_INFINITY;
        for (int[][] pair : PAIRS) {
            pairMinimum = Math.min(pairMinimum, Math.abs(input[pair[0][0]] + input[pair[0][1]]
                - input[pair[1][0]] - input[pair[1][1]]) / 4);
        }
        checks.add(new Check("best_RMSE_vs_independent_minimum_of_three_pairs_MHz", best.rmse(), pairMinimum, tol));

        // Independent tests of the EFG invariant: diagonal and rotated tensors.
        // W = (eQ/h)*V in MHz; Tr(W^2) must equal 6*nu_NQR^2 for I=3/2.
        double invariantError = 0, rotatedError = 0, traceError = 0;
        double angle = 0.731;
        double[][] rotation = {
            {Math.cos(angle), -Math.sin(angle), 0},
            {Math.sin(angle), Math.cos(angle), 0}, {0, 0, 1}
        };
        for (double nu : input) {
            for (double eta : new double[]{0, 0.37, 1}) {
                double wzz = 2 * nu / Math.sqrt(1 + eta * eta / 3);
                double[] eigenvalues = {-(1 - eta) * wzz / 2, -(1 + eta) * wzz / 2, wzz};
                double invariant = dot(eigenvalues, eigenvalues);
                double target = 6 * nu * nu;
                invariantError = Math.max(invariantError, Math.abs(invariant / target - 1));
                double[][] rotated = new double[3][3];
                for (int i = 0; i < 3; i++) {
                    for (int j = 0; j < 3; j++) {
                        for (int k = 0; k < 3; k++) {
                            rotated[i][j] += rotation[i][k] * eigenvalues[k] * rotation[j][k];
                        }
                    }
                }
                double traceSquare = 0;
                for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) {
                    traceSquare += rotated[i][j] * rotated[j][i];
                }
                rotatedError = Math.max(rotatedError, Math.abs(traceSquare / target - 1));
                traceError = Math.max(traceError,
                    Math.abs((rotated[0][0] + rotated[1][1] + rotated[2][2]) / wzz));
            }
        }
        checks.add(new Check("max_relative_invariant_error_12_diagonal_tensors", invariantError, 0, 2e-14));
        checks.add(new Check("max_relative_invariant_error_12_rotated_tensors", rotatedError, 0, 2e-14));
        checks.add(new Check("max_normalized_trace_of_rotated_tensors", traceError, 0, 2e-14));
        checks.add(new Check("eta_only_bound_squared", square(Math.sqrt(4.0 / 3.0)), 4.0 / 3.0, 1e-14));

        if (Arrays.equals(input, DEFAULT_FREQUENCIES)) {
            checks.add(new Check("reference_best_label_order_1_2_4_3",
                Arrays.equals(best.labels(), new int[]{0, 1, 3, 2}) ? 1 : 0, 1, 0));
            checks.add(new Check("reference_nu0_MHz", best.nu0(), 15.095, tol));
            checks.add(new Check("reference_b_MHz", best.b(), -2.495, tol));
            checks.add(new Check("reference_c_MHz", best.c(), -0.455, tol));
            checks.add(new Check("reference_d2_MHz", best.d2(), -0.630, tol));
            checks.add(new Check("reference_first_harmonic_RMSE_MHz", best.rmse(), 0.630, tol));
        }
        return checks;
    }

    private static void writeOrders(Path path, List<Fit> fits) throws IOException {
        try (BufferedWriter w = writer(path)) {
            csv(w, "order_id", "opposite_pair_class", "m0_input_label", "m1_input_label", "m2_input_label", "m3_input_label",
                "f0_MHz", "f1_MHz", "f2_MHz", "f3_MHz", "nu0_MHz", "b_MHz", "c_MHz", "d2_MHz",
                "first_harmonic_amplitude_MHz", "first_harmonic_phase_rad_cos_plus", "first_harmonic_phase_deg_cos_plus",
                "first_harmonic_RMSE_MHz", "first_harmonic_SSE_MHz2", "second_harmonic_exact_RMSE_MHz",
                "abs_d2_MHz", "abs_d2_over_first_amplitude", "second_harmonic_sine_coefficient_identifiable",
                "second_harmonic_interpolation_residual_dof", "spatial_order_measured");
            for (Fit f : fits) csv(w, f.id(), f.pairId(), f.labels()[0] + 1, f.labels()[1] + 1,
                f.labels()[2] + 1, f.labels()[3] + 1, f.observed()[0], f.observed()[1], f.observed()[2], f.observed()[3],
                f.nu0(), f.b(), f.c(), f.d2(), f.amplitude1(),
                f.amplitude1() == 0 ? "" : f.phase1(), f.amplitude1() == 0 ? "" : Math.toDegrees(f.phase1()),
                f.rmse(), 4 * square(f.rmse()), f.exactRmse(), Math.abs(f.d2()),
                f.amplitude1() == 0 ? "" : Math.abs(f.d2()) / f.amplitude1(), false, 0, false);
        }
    }

    private static void writePairs(Path path, double[] frequencies, List<Fit> fits) throws IOException {
        try (BufferedWriter w = writer(path)) {
            csv(w, "opposite_pair_class", "pair_a_label_1", "pair_a_label_2", "pair_b_label_1", "pair_b_label_2",
                "pair_a_frequency_1_MHz", "pair_a_frequency_2_MHz", "pair_b_frequency_1_MHz", "pair_b_frequency_2_MHz",
                "pair_a_sum_MHz", "pair_b_sum_MHz", "signed_difference_a_minus_b_MHz", "absolute_sum_difference_MHz",
                "first_harmonic_RMSE_MHz", "labelled_order_count");
            for (int p = 0; p < 3; p++) {
                int[][] pair = PAIRS[p];
                double a = frequencies[pair[0][0]] + frequencies[pair[0][1]];
                double b = frequencies[pair[1][0]] + frequencies[pair[1][1]];
                int id = p + 1;
                long count = fits.stream().filter(f -> f.pairId() == id).count();
                csv(w, id, pair[0][0] + 1, pair[0][1] + 1, pair[1][0] + 1, pair[1][1] + 1,
                    frequencies[pair[0][0]], frequencies[pair[0][1]], frequencies[pair[1][0]], frequencies[pair[1][1]],
                    a, b, a - b, Math.abs(a - b), Math.abs(a - b) / 4, count);
            }
        }
    }

    private static void writeBest(Path path, Fit f) throws IOException {
        try (BufferedWriter w = writer(path)) {
            csv(w, "m", "theta_rad", "input_label", "observed_MHz", "first_harmonic_MHz", "first_minus_observed_MHz",
                "second_harmonic_exact_MHz", "exact_minus_observed_MHz");
            for (int m = 0; m < 4; m++) csv(w, m, Math.PI * m / 2, f.labels()[m] + 1, f.observed()[m],
                f.first()[m], f.first()[m] - f.observed()[m], f.exact()[m], f.exact()[m] - f.observed()[m]);
        }
    }

    private static void writeInvariants(Path path, double[] frequencies) throws IOException {
        try (BufferedWriter w = writer(path)) {
            csv(w, "input_label", "nu_NQR_MHz", "nu_NQR_squared_MHz2", "Tr_W_squared_MHz2",
                "Frobenius_norm_W_MHz", "abs_Wzz_min_eta1_MHz", "abs_Wzz_max_eta0_MHz",
                "W_definition", "assumption");
            for (int i = 0; i < frequencies.length; i++) {
                double nu = frequencies[i];
                csv(w, i + 1, nu, nu * nu, 6 * nu * nu, Math.sqrt(6) * nu,
                    Math.sqrt(3) * nu, 2 * nu, "W=(eQ/h)*V expressed in MHz",
                    "I=3/2; pure quadrupole Hamiltonian; symmetric traceless EFG; 0<=eta<=1");
            }
        }
    }

    private static void writeChecks(Path path, List<Check> checks) throws IOException {
        try (BufferedWriter w = writer(path)) {
            csv(w, "check", "value", "expected", "absolute_tolerance", "absolute_error", "status");
            for (Check c : checks) csv(w, c.name(), c.value(), c.expected(), c.tolerance(),
                Math.abs(c.value() - c.expected()), c.pass() ? "PASS" : "FAIL");
        }
    }

    private static void writeSummary(Path path, double[] frequencies, List<Fit> fits,
                                     Fit best, List<Check> checks, double tolerance) throws IOException {
        double min = Arrays.stream(frequencies).min().orElseThrow();
        double max = Arrays.stream(frequencies).max().orElseThrow();
        double bound = Math.sqrt(4.0 / 3.0);
        long bestCount = fits.stream().filter(f -> Math.abs(f.rmse() - best.rmse()) <= tolerance).count();
        long passes = checks.stream().filter(Check::pass).count();
        StringBuilder j = new StringBuilder();
        j.append("{\n  \"schema_version\": 1,\n  \"frequency_unit\": \"MHz\",\n");
        j.append("  \"input_frequencies_MHz\": ").append(jsonArray(frequencies)).append(",\n");
        j.append("  \"input_provenance\": \"User-provided line centres; no uncertainties or covariance supplied.\",\n");
        j.append("  \"weights\": \"equal\",\n  \"theta_m\": \"pi*m/2, m=0,1,2,3\",\n");
        j.append("  \"scalar_model\": \"nu0+b*cos(theta)+c*sin(theta)+d2*cos(2*theta)\",\n");
        j.append("  \"first_harmonic_phase_convention\": \"A1*cos(theta+phi1); phi1=atan2(-c,b)\",\n");
        j.append("  \"physical_spatial_order_measured\": false,\n");
        j.append("  \"order_caution\": ").append(quote(ORDER_NOTE)).append(",\n");
        j.append("  \"labelled_orders\": 24,\n  \"opposite_pair_classes\": 3,\n");
        j.append("  \"orders_per_pair_class\": 8,\n  \"best_tied_labelled_orders\": ").append(bestCount).append(",\n");
        j.append("  \"tie_selection\": \"first lexicographic input-label order within numerical tolerance\",\n");
        j.append("  \"best\": {\n    \"input_labels_one_based\": [");
        for (int i = 0; i < 4; i++) {
            if (i > 0) j.append(", ");
            j.append(best.labels()[i] + 1);
        }
        j.append("],\n    \"cyclic_frequencies_MHz\": ").append(jsonArray(best.observed())).append(",\n");
        jsonNumber(j, "nu0_MHz", best.nu0(), 4);
        jsonNumber(j, "b_MHz", best.b(), 4);
        jsonNumber(j, "c_MHz", best.c(), 4);
        jsonNumber(j, "d2_MHz", best.d2(), 4);
        jsonNumber(j, "first_harmonic_amplitude_MHz", best.amplitude1(), 4);
        jsonNumber(j, "first_harmonic_phase_rad_cos_plus", best.amplitude1() == 0 ? Double.NaN : best.phase1(), 4);
        jsonNumber(j, "first_harmonic_phase_deg_cos_plus", best.amplitude1() == 0 ? Double.NaN : Math.toDegrees(best.phase1()), 4);
        jsonNumber(j, "first_harmonic_RMSE_MHz", best.rmse(), 4);
        jsonNumber(j, "first_harmonic_SSE_MHz2", 4 * square(best.rmse()), 4);
        jsonNumber(j, "abs_d2_over_first_amplitude", best.amplitude1() == 0 ? Double.NaN : Math.abs(best.d2()) / best.amplitude1(), 4);
        j.append("    \"first_harmonic_predicted_MHz\": ").append(jsonArray(best.first())).append(",\n");
        j.append("    \"second_harmonic_exact_MHz\": ").append(jsonArray(best.exact())).append(",\n");
        jsonNumber(j, "second_harmonic_exact_RMSE_MHz", best.exactRmse(), 4);
        j.append("    \"second_harmonic_sine_coefficient_identifiable\": false,\n");
        j.append("    \"second_harmonic_interpolation_residual_dof\": 0\n  },\n");
        j.append("  \"tensor_invariant\": {\n");
        j.append("    \"assumptions\": \"I=3/2; pure quadrupole Hamiltonian; symmetric traceless EFG\",\n");
        j.append("    \"formula_SI\": \"nu_Hz^2=(e*Q)^2*Tr(V^2)/(6*h^2)\",\n");
        j.append("    \"principal_values_identity\": \"Tr(V^2)=(3/2)*Vzz^2*(1+eta^2/3)\",\n");
        j.append("    \"normalized_tensor\": \"W=(eQ/h)*V expressed in MHz\",\n");
        j.append("    \"normalized_identity\": \"Tr(W^2)=6*nu_NQR_MHz^2\",\n");
        j.append("    \"normalization_uses_no_selected_nuclear_Q_value\": true,\n");
        j.append("    \"full_EFG_tensor_identified\": false\n  },\n");
        j.append("  \"eta_only_test\": {\n");
        j.append("    \"assumptions\": \"same isotope; fixed abs(Vzz); 0<=eta<=1; pure quadrupole frequencies\",\n");
        jsonNumber(j, "maximum_possible_frequency_ratio", bound, 4);
        jsonNumber(j, "observed_maximum_to_minimum_frequency_ratio", max / min, 4);
        jsonNumber(j, "observed_Tr_V_squared_maximum_to_minimum_ratio", square(max / min), 4);
        j.append("    \"fixed_abs_Vzz_compatible\": ").append(max / min <= bound).append(",\n");
        jsonNumber(j, "minimum_required_abs_Vzz_maximum_to_minimum_ratio", Math.max(1, max / min / bound), 4);
        j.append("    \"bound_is_not_a_charge_or_displacement_estimate\": true\n  },\n");
        j.append("  \"interpretation_cautions\": [\n    ").append(quote(INTERPOLATION_NOTE)).append(",\n");
        j.append("    \"The scalar first-harmonic check does not exclude a first harmonic of the EFG tensor.\",\n");
        j.append("    \"No statistical rejection level is inferred without experimental frequency uncertainties.\",\n");
        j.append("    \"The fitted nu0 at the input temperature is not a predicted high-temperature resonance.\"\n  ],\n");
        j.append("  \"validation\": {\"passed\": ").append(passes).append(", \"failed\": ")
            .append(checks.size() - passes).append("}\n}\n");
        Files.writeString(path, j.toString(), StandardCharsets.UTF_8);
    }

    /** A deterministic SVG: curves evaluate the displayed scalar formula directly. */
    private static void writeFigure(Path path, Fit best) throws IOException {
        if (path.getParent() != null) Files.createDirectories(path.getParent());
        double low = Arrays.stream(best.observed()).min().orElseThrow();
        double high = Arrays.stream(best.observed()).max().orElseThrow();
        for (int n = 0; n <= 360; n++) {
            double theta = Math.PI * n / 240.0;
            double first = best.nu0() + best.b() * Math.cos(theta) + best.c() * Math.sin(theta);
            double second = first + best.d2() * Math.cos(2 * theta);
            low = Math.min(low, Math.min(first, second));
            high = Math.max(high, Math.max(first, second));
        }
        double span = Math.max(high - low, Math.max(Math.abs(high) * 0.05, 1e-10));
        double rawStep = span / 6;
        double magnitude = Math.pow(10, Math.floor(Math.log10(rawStep)));
        double normalized = rawStep / magnitude;
        double step = (normalized <= 1 ? 1 : normalized <= 2 ? 2 : normalized <= 5 ? 5 : 10) * magnitude;
        double yMin = Math.floor((low - 0.08 * span) / step) * step;
        double yMax = Math.ceil((high + 0.08 * span) / step) * step;
        final double left = 115, right = 1085, top = 170, bottom = 535;
        StringBuilder s = new StringBuilder();
        s.append("""
            <?xml version="1.0" encoding="UTF-8"?>
            <svg xmlns="http://www.w3.org/2000/svg" width="1200" height="735" viewBox="0 0 1200 735" role="img" aria-labelledby="title desc">
            <title id="title">Four NQR lines and a hypothetical period-four scalar modulation</title>
            <desc id="desc">Four measured line centres in one assumed spatial order. The first harmonic has a nonzero residual. Adding a second harmonic interpolates all four samples by construction; this does not establish a physical superstructure.</desc>
            <rect width="1200" height="735" fill="#ffffff"/>
            <g font-family="DejaVu Sans,Arial,sans-serif" fill="#192433">
            <text x="55" y="48" font-size="27" font-weight="700">Four NQR lines: a hypothetical period-four assignment</text>
            """);
        svgText(s, 55, 82, 16, "Equal weights · first-harmonic RMSE = " + String.format(Locale.ROOT, "%.3f", best.rmse())
            + " MHz · spatial order is not measured", "start", "#475569");
        s.append("<line x1=\"115\" y1=\"126\" x2=\"153\" y2=\"126\" stroke=\"#192433\" stroke-width=\"0\"/>\n");
        svgCircle(s, 134, 126, 4.5, "#192433", "#192433", 1);
        svgText(s, 167, 132, 15, "Observed centres", "start", "#192433");
        s.append("<line x1=\"365\" y1=\"126\" x2=\"403\" y2=\"126\" stroke=\"#c66b13\" stroke-width=\"3\"/>\n");
        svgText(s, 418, 132, 15, "First harmonic", "start", "#192433");
        s.append("<line x1=\"620\" y1=\"126\" x2=\"658\" y2=\"126\" stroke=\"#1d5ca8\" stroke-width=\"2.5\" stroke-dasharray=\"8 5\"/>\n");
        svgCircle(s, 639, 126, 7.5, "#ffffff", "#1d5ca8", 2);
        svgText(s, 674, 132, 15, "Two-harmonic interpolation", "start", "#192433");
        for (int k = 0; k <= Math.round((yMax - yMin) / step); k++) {
            double value = yMin + k * step;
            double y = plotY(value, yMin, yMax, top, bottom);
            s.append(String.format(Locale.ROOT,
                "<line x1=\"%.3f\" y1=\"%.3f\" x2=\"%.3f\" y2=\"%.3f\" stroke=\"#e2e8f0\"/>%n", left, y, right, y));
            svgText(s, left - 16, y + 5, 14, shortNumber(value, 4), "end", "#475569");
        }
        for (int m = 0; m < 4; m++) {
            double x = left + m * (right - left) / 3;
            s.append(String.format(Locale.ROOT,
                "<line x1=\"%.3f\" y1=\"%.3f\" x2=\"%.3f\" y2=\"%.3f\" stroke=\"#e2e8f0\"/>%n", x, top, x, bottom));
            svgText(s, x, bottom + 28, 15, "m = " + m, "middle", "#192433");
        }
        s.append(String.format(Locale.ROOT,
            "<path d=\"M %.3f %.3f L %.3f %.3f L %.3f %.3f\" fill=\"none\" stroke=\"#64748b\" stroke-width=\"1.3\"/>%n",
            left, top, left, bottom, right, bottom));
        for (int harmonic = 1; harmonic <= 2; harmonic++) {
            s.append("<path d=\"");
            for (int n = 0; n <= 360; n++) {
                double m = n / 120.0;
                double theta = Math.PI * m / 2;
                double frequency = best.nu0() + best.b() * Math.cos(theta) + best.c() * Math.sin(theta);
                if (harmonic == 2) frequency += best.d2() * Math.cos(2 * theta);
                s.append(String.format(Locale.ROOT, "%s %.3f %.3f ", n == 0 ? "M" : "L",
                    left + m * (right - left) / 3, plotY(frequency, yMin, yMax, top, bottom)));
            }
            s.append("\" fill=\"none\" stroke=\"").append(harmonic == 1 ? "#c66b13" : "#1d5ca8")
                .append("\" stroke-width=\"2.7\"");
            if (harmonic == 2) s.append(" stroke-dasharray=\"8 5\"");
            s.append("/>\n");
        }
        for (int m = 0; m < 4; m++) {
            double x = left + m * (right - left) / 3;
            double observedY = plotY(best.observed()[m], yMin, yMax, top, bottom);
            double firstY = plotY(best.first()[m], yMin, yMax, top, bottom);
            s.append(String.format(Locale.ROOT,
                "<line x1=\"%.3f\" y1=\"%.3f\" x2=\"%.3f\" y2=\"%.3f\" stroke=\"#c66b13\" stroke-width=\"1.2\" stroke-dasharray=\"3 4\"/>%n",
                x, observedY, x, firstY));
            s.append(String.format(Locale.ROOT,
                "<rect x=\"%.3f\" y=\"%.3f\" width=\"9\" height=\"9\" fill=\"#fff\" stroke=\"#c66b13\" stroke-width=\"2\"/>%n",
                x - 4.5, firstY - 4.5));
            svgCircle(s, x, plotY(best.exact()[m], yMin, yMax, top, bottom), 8, "#fff", "#1d5ca8", 2);
            svgCircle(s, x, observedY, 4, "#192433", "#192433", 1);
            svgText(s, x + (m == 3 ? -15 : 15), observedY + (firstY < observedY ? 22 : -15), 14,
                shortNumber(best.observed()[m], 4), m == 3 ? "end" : "start", "#192433");
        }
        svgText(s, (left + right) / 2, bottom + 64, 16,
            "Subcell label m (hypothetical assignment)", "middle", "#192433");
        s.append("<text x=\"32\" y=\"352\" font-size=\"16\" text-anchor=\"middle\" transform=\"rotate(-90 32 352)\">Frequency (MHz)</text>\n");
        svgText(s, 55, 642, 15,
            "Blue rings coincide with the observations by construction: four coefficients for four samples.", "start", "#475569");
        svgText(s, 55, 672, 14,
            "A1 = " + shortNumber(best.amplitude1(), 6) + " MHz; φ1 = " + shortNumber(Math.toDegrees(best.phase1()), 6)
            + "° for cos(θ + φ1); d2 = " + shortNumber(best.d2(), 6) + " MHz.", "start", "#475569");
        svgText(s, 55, 702, 14,
            "Zero interpolation residual does not establish a crystal structure, propagation vector or physical mechanism.", "start", "#475569");
        s.append("</g>\n</svg>\n");
        Files.writeString(path, s.toString(), StandardCharsets.UTF_8);
    }

    private static double plotY(double frequency, double low, double high, double top, double bottom) {
        return bottom - (frequency - low) * (bottom - top) / (high - low);
    }

    private static void svgText(StringBuilder s, double x, double y, int size, String text, String anchor, String fill) {
        s.append(String.format(Locale.ROOT,
            "<text x=\"%.3f\" y=\"%.3f\" font-size=\"%d\" text-anchor=\"%s\" fill=\"%s\">%s</text>%n",
            x, y, size, anchor, fill, text.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")));
    }

    private static void svgCircle(StringBuilder s, double x, double y, double radius, String fill, String stroke, double width) {
        s.append(String.format(Locale.ROOT,
            "<circle cx=\"%.3f\" cy=\"%.3f\" r=\"%.3f\" fill=\"%s\" stroke=\"%s\" stroke-width=\"%.2f\"/>%n",
            x, y, radius, fill, stroke, width));
    }

    private static String shortNumber(double value, int decimals) {
        String result = String.format(Locale.ROOT, "%1$." + decimals + "f", value);
        return result.indexOf('.') < 0 ? result : result.replaceAll("0+$", "").replaceAll("\\.$", "");
    }

    private static BufferedWriter writer(Path path) throws IOException {
        return Files.newBufferedWriter(path, StandardCharsets.UTF_8);
    }

    private static void csv(BufferedWriter w, Object... fields) throws IOException {
        for (int i = 0; i < fields.length; i++) {
            if (i > 0) w.write(',');
            String value = fields[i] instanceof Double d ? number(d) : fields[i].toString();
            if (value.contains(",") || value.contains("\"") || value.contains("\n")) {
                w.write('"');
                w.write(value.replace("\"", "\"\""));
                w.write('"');
            } else w.write(value);
        }
        w.write('\n');
    }

    private static String number(double x) {
        if (!Double.isFinite(x)) return "null";
        if (x == 0) return "0";
        return String.format(Locale.ROOT, "%.15g", x);
    }

    private static void jsonNumber(StringBuilder j, String key, double value, int indent) {
        j.append(" ".repeat(indent)).append(quote(key)).append(": ").append(number(value)).append(",\n");
    }

    private static String jsonArray(double[] array) {
        StringBuilder result = new StringBuilder("[");
        for (int i = 0; i < array.length; i++) {
            if (i > 0) result.append(", ");
            result.append(number(array[i]));
        }
        return result.append(']').toString();
    }

    private static String quote(String text) {
        return '"' + text.replace("\\", "\\\\").replace("\"", "\\\"")
            .replace("\n", "\\n").replace("\r", "\\r").replace("\t", "\\t") + '"';
    }
}
