#!/usr/bin/env python3
"""Export auditable anilite point-charge EFG calculations and every Ewald term.

This is an independent implementation of the archived geometrical proxy, not a
calculation of the complete electronic EFG.  Charges are trial numbers z_j in
q_j = z_j e.  Numerical G has units Angstrom^-3 in this convention; equivalently
one sometimes describes it as e/Angstrom^3 before restoring the Coulomb factor.

Run from any directory: python theory/export_efg_intermediates.py
No experimental input is modified.  The old geometry_efg_proxy.py is not run.
"""
from __future__ import annotations

import csv
import gzip
import hashlib
import io
import itertools
import json
import platform
import re
from pathlib import Path

import numpy as np
import scipy
from scipy.constants import elementary_charge, epsilon_0
from scipy.special import erfc


HERE = Path(__file__).resolve().parent
BASE = HERE.parent
OUT = HERE / "results" / "efg"
CIF = BASE / "research_literature" / "anilite_amcsd0009376_transcribed.cif"
ARCHIVED = BASE / "new_evidence" / "structure_model" / "geometry_efg_proxy.json"
FIT_TABLE = BASE / "powder_model" / "FINAL_MODEL_COMPARISON.csv"
COMPONENTS = ("aa", "bb", "cc", "ab", "ac", "bc")
COMPONENT_PAIRS = ((0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2))
ALPHA_REFERENCE = 0.3
ALPHAS = (0.2, 0.3, 0.4)
DIRECT_RADIUS = 100.0
DIRECT_CUTOFFS = (15.0, 20.0, 30.0, 40.0, 60.0, 80.0, 100.0)
SI_FACTOR = elementary_charge / (4 * np.pi * epsilon_0) * 1e30


def save_json(name, data):
    (OUT / name).write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")


def save_csv(name, rows):
    rows = iter(rows)
    first = next(rows)
    path = OUT / name
    if name.endswith(".gz"):
        # mtime and embedded filename are fixed for deterministic gzip bytes.
        with path.open("wb") as raw:
            with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as gz:
                with io.TextIOWrapper(gz, encoding="utf-8", newline="") as text:
                    writer = csv.DictWriter(text, fieldnames=list(first))
                    writer.writeheader()
                    writer.writerow(first)
                    writer.writerows(rows)
    else:
        with path.open("w", newline="", encoding="utf-8") as text:
            writer = csv.DictWriter(text, fieldnames=list(first))
            writer.writeheader()
            writer.writerow(first)
            writer.writerows(rows)


def six(t):
    return np.array([t[..., i, j] for i, j in COMPONENT_PAIRS]).T


def matrix(v):
    t = np.zeros((3, 3), dtype=float)
    for x, (i, j) in zip(v, COMPONENT_PAIRS):
        t[i, j] = t[j, i] = x
    return t


def component_dict(t, prefix="G_"):
    return {prefix + label: float(t[i, j]) for label, (i, j) in zip(COMPONENTS, COMPONENT_PAIRS)}


def properties(t):
    values, axes = np.linalg.eigh((t + t.T) / 2)
    order = np.argsort(abs(values))
    values, axes = values[order], axes[:, order]
    # Eigenvector signs have no physical meaning; canonicalize for saved output.
    for j in range(3):
        if axes[np.argmax(abs(axes[:, j])), j] < 0:
            axes[:, j] *= -1
    eta = float((values[0] - values[1]) / values[2]) if abs(values[2]) > 1e-14 else None
    assert eta is None or -1e-10 <= eta <= 1 + 1e-10
    return {
        "tensor_abc": t.tolist(),
        "eigenvalues_XX_YY_ZZ": values.tolist(),
        "principal_axes_columns_XX_YY_ZZ_in_abc": axes.tolist(),
        "eta_pc": eta,
        "abs_GZZ": float(abs(values[2])),
        "trace": float(np.trace(t)),
        "frobenius_norm": float(np.linalg.norm(t)),
        "determinant": float(np.linalg.det(t)),
        "eigen_reconstruction_max_error": float(np.max(abs(axes @ np.diag(values) @ axes.T - t))),
        "bare_lattice_EFG_SI_V_per_m2": (t * SI_FACTOR).tolist(),
        "bare_lattice_abs_VZZ_SI_V_per_m2": float(abs(values[2]) * SI_FACTOR),
        "relative_NQR_shape_proxy_P": float(abs(values[2]) * np.sqrt(1 + eta**2 / 3)) if eta is not None else 0.0,
    }


def translations(nmax):
    return np.array(list(itertools.product(*(range(-int(n), int(n) + 1) for n in nmax))), dtype=int)


def read_structure():
    text = CIF.read_text()
    cell = np.array([float(re.search(rf"_cell_length_{a}\s+(\S+)", text).group(1)) for a in "abc"])
    sites = {
        m.group(1): np.array([float(m.group(i)) for i in (2, 3, 4)])
        for m in re.finditer(r"^(Cu\d+|S\d+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)", text, re.M)
    }
    # Explicit affine operations, in the same setting as the preserved CIF.
    operators = [
        (np.diag([1, 1, 1]), [0, 0, 0]),
        (np.diag([1, -1, 1]), [0, .5, 0]),
        (np.diag([-1, 1, -1]), [0, .5, 0]),
        (np.diag([-1, 1, 1]), [.5, .5, .5]),
        (np.diag([1, -1, -1]), [.5, .5, .5]),
        (np.diag([1, 1, -1]), [.5, 0, .5]),
        (np.diag([-1, -1, 1]), [.5, 0, .5]),
        (np.diag([-1, -1, -1]), [0, 0, 0]),
    ]
    expanded, symmetry, multiplicities = [], [], {}
    for label, position in sites.items():
        orbit = {}
        for op_id, (rotation, shift) in enumerate(operators, start=1):
            frac = tuple(np.round((rotation @ position + shift) % 1, 8))
            orbit.setdefault(frac, []).append(op_id)
            symmetry.append({"orbit": label, "operator_id": op_id,
                             "fractional_a": frac[0], "fractional_b": frac[1], "fractional_c": frac[2]})
        multiplicities[label] = len(orbit)
        for orbit_index, (frac, op_ids) in enumerate(sorted(orbit.items())):
            expanded.append({
                "atom_index": len(expanded), "orbit": label, "orbit_index": orbit_index,
                "symmetry_operator_ids": ";".join(map(str, op_ids)),
                **{f"fractional_{a}": float(frac[i]) for i, a in enumerate("abc")},
                **{f"cartesian_{a}_A": float(frac[i] * cell[i]) for i, a in enumerate("abc")},
            })
    assert len(expanded) == 44 and sum(x["orbit"].startswith("Cu") for x in expanded) == 28
    labels = np.array([row["orbit"] for row in expanded])
    coordinates = np.array([[row[f"fractional_{a}"] for a in "abc"] for row in expanded])
    return cell, sites, expanded, symmetry, multiplicities, labels, coordinates


def make_charge_scenarios(sites, multiplicities):
    scenarios = {"uniform_Cu_8over7": {s: 8 / 7 if s.startswith("Cu") else -2.0 for s in sites}}
    for target in ("Cu1", "Cu2", "Cu3", "Cu4", "Cu5"):
        q = {s: 1.0 if s.startswith("Cu") else -2.0 for s in sites}
        q[target] += 4 / multiplicities[target]
        scenarios[f"excess_charge_on_{target}"] = q
    for name, q in scenarios.items():
        assert abs(sum(q[s] * multiplicities[s] for s in sites)) < 1e-12, name
    return scenarios


def ewald_basis(observation, cell, atom_cart, alpha):
    """Return uncharged real kernels and reciprocal geometry, self excluded."""
    rc, kc = 6 / alpha, 12 * alpha
    image = translations(np.ceil(rc / cell).astype(int) + 1)
    source_index = np.tile(np.arange(len(atom_cart)), len(image))
    image_expanded = np.repeat(image, len(atom_cart), axis=0)
    r = np.tile(atom_cart, (len(image), 1)) + image_expanded * cell - observation
    radius = np.linalg.norm(r, axis=1)
    keep = (radius > 1e-8) & (radius < rc)
    source_index, image_expanded, r, radius = source_index[keep], image_expanded[keep], r[keep], radius[keep]
    damp = np.exp(-(alpha * radius)**2)
    complementary_error = erfc(alpha * radius)
    radial_A = (3 * complementary_error / radius**5
                + 6 * alpha / np.sqrt(np.pi) * damp / radius**4
                + 4 * alpha**3 / np.sqrt(np.pi) * damp / radius**2)
    radial_B = -complementary_error / radius**3 - 2 * alpha / np.sqrt(np.pi) * damp / radius**2
    kernels = radial_A[:, None, None] * r[:, :, None] * r[:, None, :] + radial_B[:, None, None] * np.eye(3)
    indices = translations(np.ceil(kc * cell / (2 * np.pi)).astype(int))
    k = indices * 2 * np.pi / cell
    k2 = np.einsum("ni,ni->n", k, k)
    keep = (k2 > 0) & (k2 < kc**2)
    indices, k, k2 = indices[keep], k[keep], k2[keep]
    coefficient = -4 * np.pi / np.prod(cell) * np.exp(-k2 / (4 * alpha**2)) / k2
    return {
        "source_index": source_index, "image": image_expanded, "r": r, "radius": radius,
        "erfc": complementary_error, "exp": damp, "A": radial_A, "B": radial_B, "kernels": kernels,
        "hkl": indices, "k": k, "k2": k2, "reciprocal_coefficient": coefficient,
        "cos_atoms": np.cos(k @ atom_cart.T), "sin_atoms": np.sin(k @ atom_cart.T),
        "cos_observation": np.cos(k @ observation), "sin_observation": np.sin(k @ observation),
        "real_cutoff_A": rc, "reciprocal_cutoff_Ainv": kc,
    }


def evaluate_ewald(basis, q_atoms, self_charge, alpha):
    real_terms = basis["kernels"] * q_atoms[basis["source_index"], None, None]
    structure_re = basis["cos_atoms"] @ q_atoms
    structure_im = -basis["sin_atoms"] @ q_atoms
    charge_cos = basis["cos_observation"] * structure_re - basis["sin_observation"] * structure_im
    reciprocal_terms = (basis["reciprocal_coefficient"] * charge_cos)[:, None, None] * basis["k"][:, :, None] * basis["k"][:, None, :]
    real_sum, reciprocal_sum = real_terms.sum(axis=0), reciprocal_terms.sum(axis=0)
    self_term = np.eye(3) * 4 * alpha**3 / (3 * np.sqrt(np.pi)) * self_charge
    total_before = real_sum + reciprocal_sum + self_term
    trace_before = float(np.trace(total_before))
    assert abs(trace_before) < 1e-9
    projection = -np.eye(3) * trace_before / 3
    total = total_before + projection
    return {
        "real_terms": real_terms, "reciprocal_terms": reciprocal_terms,
        "structure_re": structure_re, "structure_im": structure_im, "charge_cos": charge_cos,
        "real_sum": real_sum, "reciprocal_sum": reciprocal_sum, "self_term": self_term,
        "roundoff_trace_projection": projection, "total_before_projection": total_before,
        "trace_before_projection": trace_before, "total": total,
    }


def empirical_frequency_comparison(tensors):
    observed = np.array([11.97, 15.27, 16.18, 16.96])
    high_sites = ("Cu1", "Cu2", "Cu4", "Cu5")
    with FIT_TABLE.open(newline="") as f:
        fit_rows = list(csv.DictReader(f))
    branch_eta = {
        ("L" if row["branch"] == "29" else "H"): np.array([float(row[f"eta_s{i}"]) for i in range(1, 5)])
        for row in fit_rows if row["model"] == "Central transitions; nonnegative linear background"
    }
    all_rows, group_rows, best = [], [], []
    for scenario, sites in tensors.items():
        rows = []
        for order in itertools.permutations(high_sites):
            proxies = np.array([sites[s]["relative_NQR_shape_proxy_P"] for s in order])
            eta = np.array([sites[s]["eta_pc"] for s in order])
            scale = float(np.dot(proxies, observed) / np.dot(proxies, proxies))
            predicted = scale * proxies
            residual = predicted - observed
            fractions = residual / observed
            cu3_prediction = scale * sites["Cu3"]["relative_NQR_shape_proxy_P"]
            row = {
                "scenario": scenario, "assignment_s1_s2_s3_s4": ";".join(order),
                "scale_MHz_per_proxy_unit": scale,
                "high_groups_SSE_MHz2": float(np.sum(residual**2)),
                "high_groups_RMSE_MHz": float(np.sqrt(np.mean(residual**2))),
                "high_groups_fractional_RMS": float(np.sqrt(np.mean(fractions**2))),
                "high_groups_max_abs_fractional_residual": float(np.max(abs(fractions))),
                "scaled_Cu3_MHz": float(cu3_prediction),
                "Cu3_within_reported_1p5_to_2p2_MHz_band": bool(1.5 <= cu3_prediction <= 2.2),
                "eta_RMS_vs_NMR_branch_L": float(np.sqrt(np.mean((eta - branch_eta["L"])**2))),
                "eta_RMS_vs_NMR_branch_H": float(np.sqrt(np.mean((eta - branch_eta["H"])**2))),
            }
            for i in range(4):
                row[f"predicted_s{i+1}_MHz"] = float(predicted[i])
                row[f"residual_s{i+1}_MHz"] = float(residual[i])
                row[f"fractional_residual_s{i+1}"] = float(fractions[i])
                row[f"model_eta_at_s{i+1}"] = float(eta[i])
            rows.append(row)
        rows.sort(key=lambda row: row["high_groups_SSE_MHz2"])
        for rank, row in enumerate(rows, start=1):
            row["frequency_SSE_rank_within_scenario"] = rank
            all_rows.append(row)
        best_row = dict(rows[0])
        best_row["all_24_assignments_Cu3_MHz_min"] = min(row["scaled_Cu3_MHz"] for row in rows)
        best_row["all_24_assignments_Cu3_MHz_max"] = max(row["scaled_Cu3_MHz"] for row in rows)
        best.append(best_row)
        for i, site in enumerate(best_row["assignment_s1_s2_s3_s4"].split(";")):
            group_rows.append({
                "scenario": scenario, "group": f"s{i+1}", "assigned_site_in_best_frequency_only_mapping": site,
                "observed_NQR_MHz": observed[i], "scaled_proxy_NQR_MHz": best_row[f"predicted_s{i+1}_MHz"],
                "residual_MHz": best_row[f"residual_s{i+1}_MHz"],
                "fractional_residual": best_row[f"fractional_residual_s{i+1}"],
                "point_charge_eta": sites[site]["eta_pc"],
                "NMR_eta_branch_L": branch_eta["L"][i], "NMR_eta_branch_H": branch_eta["H"][i],
            })
    save_csv("frequency_scaling_all_144_assignments.csv", all_rows)
    save_csv("frequency_scaling_best_by_scenario.csv", best)
    save_csv("frequency_scaling_best_assignment_group_residuals.csv", group_rows)
    save_json("frequency_scaling.json", {
        "status": "Empirical one-scale consistency exercise; neither an absolute EFG prediction nor a crystallographic assignment.",
        "observed_high_groups_MHz": observed.tolist(),
        "experimental_source": "Author thesis, 4.2 K; four 63Cu high-frequency NQR groups.",
        "reported_low_frequency_band_MHz": [1.5, 2.2],
        "proxy_definition": "P_s = |G_ZZ,s| sqrt(1 + eta_pc,s^2/3)",
        "scale_definition": "a = sum(P_i nu_i) / sum(P_i^2), for four specified high-site assignments.",
        "assignment_search": "All 24 bijections of Cu1,Cu2,Cu4,Cu5 to s1,s2,s3,s4, separately in all six charge scenarios.",
        "objective": "Unweighted frequency SSE only; eta residuals are diagnostic and not included in selecting the minimum.",
        "no_error_bars": "No noise model or statistical confidence interval is inferred; scenario/assignment spreads are conditional model sensitivity.",
        "assumption": "The same unknown positive multiplicative scale is applied to all five sites. Real electronic contributions may violate it.",
        "best_by_scenario": best,
        "best_mapping_Cu3_MHz_range_all_scenarios": [min(x["scaled_Cu3_MHz"] for x in best), max(x["scaled_Cu3_MHz"] for x in best)],
        "all_144_Cu3_MHz_range": [min(x["scaled_Cu3_MHz"] for x in all_rows), max(x["scaled_Cu3_MHz"] for x in all_rows)],
        "interpretation": "Any agreement with the 1.5--2.2 MHz band supports investigating Cu3, but cannot identify an unresolved signal in a mixed Cu1.1S sample. High-group spacing and eta mismatches remain independent checks.",
    })
    return best


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    cell, sites, expanded, symmetry, multiplicities, atom_labels, fractional = read_structure()
    atom_cart = fractional * cell
    scenarios = make_charge_scenarios(sites, multiplicities)
    save_csv("expanded_atoms.csv", expanded)
    save_csv("symmetry_expansion_before_deduplication.csv", symmetry)
    save_csv("charge_scenarios.csv", [{"scenario": name, "orbit": site, "multiplicity": multiplicities[site],
                                     "trial_charge_z_in_q_equals_z_e": q,
                                     "orbit_charge_in_e": multiplicities[site] * q}
                                    for name, values in scenarios.items() for site, q in values.items()])
    save_json("structure_and_conventions.json", {
        "source_CIF": str(CIF.relative_to(BASE)), "source_CIF_sha256": hashlib.sha256(CIF.read_bytes()).hexdigest(),
        "cell_orthorhombic_A": cell.tolist(), "volume_from_lengths_A3": float(np.prod(cell)),
        "source_note": "Preserved transcription of AMCSD0009376 / Koto--Morimoto1970; room-temperature structural reference, not a 5K refinement.",
        "asymmetric_unit": {key: value.tolist() for key, value in sites.items()}, "multiplicities": multiplicities,
        "coordinates": "Cartesian components use orthonormal crystallographic a,b,c directions; lengths in Angstrom.",
        "tensor_convention": "G is the Hessian of electrostatic potential, not the gradient of E=-grad(phi). q_j=z_j e with trial dimensionless z_j.",
        "numerical_tensor_unit": "Angstrom^-3 for dimensionless z_j; also commonly called e/Angstrom^3 when retaining the charge unit.",
        "SI_conversion": "V_lattice,ab [V/m^2] = G_ab * e/(4*pi*epsilon_0) * 10^30",
        "SI_factor_V_per_m2_per_proxy_unit": float(SI_FACTOR),
        "constants": {"elementary_charge_C": elementary_charge, "epsilon0_F_per_m": epsilon_0,
                      "source": "scipy.constants from the recorded SciPy version"},
        "alpha_reference_Ainv": ALPHA_REFERENCE, "convergence_alphas_Ainv": ALPHAS,
        "Ewald_cutoffs": "R_c=6/alpha and K_c=12*alpha; strict inequalities; the singular self image is omitted.",
        "zero_mode": "All six cells are neutral; the reciprocal k=0 term is omitted.",
        "principal_axis_order": "|G_ZZ| >= |G_YY| >= |G_XX|; eta=(G_XX-G_YY)/G_ZZ. Axis signs are canonicalized only for file reproducibility.",
        "roundoff_projection": "Saved separately: -(Tr G_raw)/3 I. The raw trace is checked before this negligible numerical correction.",
        "model_limitations": ["Trial charges are not measured oxidation states.", "No on-site orbital electronic EFG, covalence, screening or antishielding.",
                              "Bare SI numbers are unscreened lattice proxies, not the physical nuclear EFG.",
                              "Numerical convergence does not validate the physical model."],
        "software": {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__},
    })
    archived = json.loads(ARCHIVED.read_text())
    archived_sites = {s["site"]: s for s in archived["sites"]}
    real_rows, reciprocal_rows, subtotal_rows, convergence_rows, summary_rows = [], [], [], [], []
    tensors = {name: {} for name in scenarios}
    all_alpha_tensors = {}
    max_reference_difference = max_raw_trace = max_mirror_error = 0.0
    for site, frac in sites.items():
        if not site.startswith("Cu"):
            continue
        observation = frac * cell
        all_alpha_tensors[site] = {name: {} for name in scenarios}
        for alpha in ALPHAS:
            basis = ewald_basis(observation, cell, atom_cart, alpha)
            for scenario, charges in scenarios.items():
                q_atoms = np.array([charges[label] for label in atom_labels])
                result = evaluate_ewald(basis, q_atoms, charges[site], alpha)
                total = result["total"]
                all_alpha_tensors[site][scenario][str(alpha)] = total.tolist()
                max_raw_trace = max(max_raw_trace, abs(result["trace_before_projection"]))
                props = properties(total)
                reference = np.array(archived_sites[site]["periodic_ewald"][scenario]["tensor_abc"])
                difference = float(np.max(abs(total - reference)))
                max_reference_difference = max(max_reference_difference, difference)
                convergence_rows.append({"site": site, "scenario": scenario, "alpha_Ainv": alpha,
                                         "real_cutoff_A": basis["real_cutoff_A"], "reciprocal_cutoff_Ainv": basis["reciprocal_cutoff_Ainv"],
                                         "n_real_terms": len(basis["r"]), "n_reciprocal_terms": len(basis["k"]),
                                         "trace_before_projection": result["trace_before_projection"],
                                         "max_abs_difference_from_archived_reference": difference,
                                         "eta_pc": props["eta_pc"], "abs_GZZ": props["abs_GZZ"], **component_dict(total)})
                if multiplicities[site] == 4:
                    max_mirror_error = max(max_mirror_error, abs(total[0, 1]), abs(total[1, 2]))
                if alpha != ALPHA_REFERENCE:
                    continue
                tensors[scenario][site] = props
                summary_rows.append({"scenario": scenario, "site": site, "multiplicity": multiplicities[site],
                                     "eta_pc": props["eta_pc"], "abs_GZZ": props["abs_GZZ"],
                                     "bare_lattice_abs_VZZ_SI_V_per_m2": props["bare_lattice_abs_VZZ_SI_V_per_m2"],
                                     "relative_NQR_shape_proxy_P": props["relative_NQR_shape_proxy_P"],
                                     **component_dict(total)})
                for component in ("real_sum", "reciprocal_sum", "self_term", "roundoff_trace_projection", "total_before_projection", "total"):
                    subtotal_rows.append({"scenario": scenario, "site": site, "alpha_Ainv": alpha, "component": component,
                                          "trace": float(np.trace(result[component])), **component_dict(result[component])})
                for n, contribution in enumerate(result["real_terms"]):
                    atom = int(basis["source_index"][n]); r = basis["r"][n]
                    real_rows.append({"scenario": scenario, "site": site, "alpha_Ainv": alpha,
                                      "source_atom_index": atom, "source_orbit": str(atom_labels[atom]),
                                      **{f"image_{a}": int(basis["image"][n, i]) for i, a in enumerate("abc")},
                                      "charge_z": float(q_atoms[atom]),
                                      **{f"delta_{a}_A": float(r[i]) for i, a in enumerate("abc")},
                                      "distance_A": float(basis["radius"][n]), "erfc_alpha_r": float(basis["erfc"][n]),
                                      "exp_minus_alpha_r_squared": float(basis["exp"][n]),
                                      "radial_A_Ainv5": float(basis["A"][n]), "radial_B_Ainv3": float(basis["B"][n]),
                                      **component_dict(contribution)})
                for n, contribution in enumerate(result["reciprocal_terms"]):
                    reciprocal_rows.append({"scenario": scenario, "site": site, "alpha_Ainv": alpha,
                                            **{a: int(basis["hkl"][n, i]) for i, a in enumerate("hkl")},
                                            **{f"k_{a}_Ainv": float(basis["k"][n, i]) for i, a in enumerate("abc")},
                                            "k_squared_Ainv2": float(basis["k2"][n]),
                                            "structure_factor_real": float(result["structure_re"][n]),
                                            "structure_factor_imaginary": float(result["structure_im"][n]),
                                            "charge_cos_at_observation": float(result["charge_cos"][n]),
                                            "reciprocal_coefficient_before_charge_cos": float(basis["reciprocal_coefficient"][n]),
                                            **component_dict(contribution)})
        print(f"Saved tensors in memory: {site}", flush=True)
    save_csv("ewald_real_terms.csv.gz", real_rows)
    save_csv("ewald_reciprocal_terms.csv.gz", reciprocal_rows)
    term_counts = {"real": len(real_rows), "reciprocal": len(reciprocal_rows)}
    del real_rows, reciprocal_rows
    save_csv("ewald_subtotals.csv", subtotal_rows)
    save_csv("ewald_alpha_convergence.csv", convergence_rows)
    save_csv("tensor_summary.csv", summary_rows)
    save_json("tensors_and_eigensystems.json", tensors)
    save_json("ewald_convergence_tensors.json", all_alpha_tensors)

    # Independent direct Coulomb sum, reduced into per-orbit radial shells.
    images = translations(np.ceil(DIRECT_RADIUS / cell).astype(int) + 1)
    source_cart = (images[:, None, :] * cell + atom_cart[None, :, :]).reshape(-1, 3)
    source_labels = np.tile(atom_labels, len(images))
    shell_edges = np.unique(np.r_[np.arange(0, DIRECT_RADIUS + 1), 2.6, 3.3])
    shell_basis_rows, shell_scenario_rows, direct_rows = [], [], []
    max_old_direct_difference = 0.0
    for site, frac in sites.items():
        if not site.startswith("Cu"):
            continue
        delta = source_cart - frac * cell
        radius = np.linalg.norm(delta, axis=1)
        selected = (radius > 1e-8) & (radius <= DIRECT_RADIUS)
        r, radius, labels = delta[selected], radius[selected], source_labels[selected]
        r2 = np.sum(r**2, axis=1)
        kernel = (3 * r[:, :, None] * r[:, None, :] - r2[:, None, None] * np.eye(3)) / radius[:, None, None]**5
        kernel_six = six(kernel)
        shell_id = np.searchsorted(shell_edges, radius, side="left") - 1
        n_shells = len(shell_edges) - 1
        orbit_bins = {}
        for label in sites:
            selected_orbit = labels == label
            bins = np.stack([np.bincount(shell_id[selected_orbit], weights=kernel_six[selected_orbit, k], minlength=n_shells)
                             for k in range(6)], axis=1)
            counts = np.bincount(shell_id[selected_orbit], minlength=n_shells)
            orbit_bins[label] = bins
            for j in range(n_shells):
                shell_basis_rows.append({"site": site, "source_orbit": label, "r_lower_exclusive_A": float(shell_edges[j]),
                                         "r_upper_inclusive_A": float(shell_edges[j+1]), "n_images": int(counts[j]),
                                         **component_dict(matrix(bins[j]), prefix="kernel_")})
        for scenario, charges in scenarios.items():
            bins = sum(charges[label] * orbit_bins[label] for label in sites)
            cumulative = np.cumsum(bins, axis=0)
            ewald = np.array(tensors[scenario][site]["tensor_abc"])
            for j in range(n_shells):
                shell_scenario_rows.append({"scenario": scenario, "site": site,
                                            "r_lower_exclusive_A": float(shell_edges[j]), "r_upper_inclusive_A": float(shell_edges[j+1]),
                                            **component_dict(matrix(bins[j]), prefix="shell_G_"),
                                            **component_dict(matrix(cumulative[j]), prefix="cumulative_G_")})
            for cutoff in DIRECT_CUTOFFS:
                j = int(np.where(shell_edges == cutoff)[0][0]) - 1
                total = matrix(cumulative[j])
                props = properties(total)
                expected = np.array(archived_sites[site][scenario][str(cutoff)]["tensor_abc"])
                old_diff = float(np.max(abs(total - expected)))
                max_old_direct_difference = max(max_old_direct_difference, old_diff)
                direct_rows.append({"scenario": scenario, "site": site, "radius_A": cutoff,
                                    "eta_pc": props["eta_pc"], "abs_GZZ": props["abs_GZZ"],
                                    "max_abs_tensor_difference_from_Ewald": float(np.max(abs(total - ewald))),
                                    "max_abs_difference_from_archived_direct_sum": old_diff,
                                    **component_dict(total)})
        print(f"Saved direct shell comparison: {site}", flush=True)
    save_csv("direct_radial_shell_basis_by_source_orbit.csv.gz", shell_basis_rows)
    save_csv("direct_radial_shell_scenarios.csv", shell_scenario_rows)
    save_csv("direct_radius_convergence.csv", direct_rows)
    best = empirical_frequency_comparison(tensors)
    max_alpha_change = max(float(np.max(abs(np.array(t) - np.array(group[str(ALPHA_REFERENCE)]))))
                           for site in all_alpha_tensors.values() for group in site.values() for t in group.values())
    validation = {
        "source": "New independent exporter compared with archived numerical tensors; no old script was executed.",
        "charge_neutrality_max_abs_e": max(abs(sum(multiplicities[s] * q[s] for s in sites)) for q in scenarios.values()),
        "maximum_difference_from_archived_Ewald_tensor": max_reference_difference,
        "maximum_alpha_change": max_alpha_change,
        "maximum_abs_trace_before_projection": max_raw_trace,
        "maximum_4c_mirror_forbidden_component": max_mirror_error,
        "maximum_difference_from_archived_direct_spherical_sum": max_old_direct_difference,
        "Cu3_smallest_abs_GZZ_in_all_six_scenarios": all(min(data, key=lambda s: data[s]["abs_GZZ"]) == "Cu3" for data in tensors.values()),
        "maximum_eigensystem_reconstruction_error": max(x["eigen_reconstruction_max_error"] for ss in tensors.values() for x in ss.values()),
        "saved_Ewald_term_counts": term_counts,
        "tolerance_for_numerical_agreement": 1e-10,
        "physical_scope": "These checks validate summation/export consistency only, not the point-charge model as a physical description of nuclear EFGs.",
    }
    assert max_reference_difference < 1e-10
    assert max_alpha_change < 1e-10
    assert max_mirror_error < 1e-10
    assert max_old_direct_difference < 1e-10
    validation["all_checks_passed"] = True
    save_json("numerical_validation.json", validation)
    schema = """ANILITE EFG INTERMEDIATE NUMERICAL RECORDS

Rebuild: python theory/export_efg_intermediates.py
Verify saved summations: python theory/validate_efg_intermediates.py

G is the Hessian of the electrostatic potential with trial q_j=z_j e. Distances
are in Angstrom. Numeric tensors have unit Angstrom^-3 for dimensionless z_j;
multiply by the factor in structure_and_conventions.json to obtain the BARE,
UNSCREENED lattice EFG in V/m^2. This is not the full nuclear EFG. No absolute
physical NQR frequency is inferred from that SI conversion.

expanded_atoms.csv: the complete 44-atom unit cell. atom_index is zero-based and
is the join key used by the real-space term records. orbit_index is zero-based;
symmetry_operator_ids are one-based positions in the eight-operation CIF list.
symmetry_expansion_before_deduplication.csv preserves every application of all
eight operations, including coincident images on special positions.
charge_scenarios.csv: each of six trial neutral cells, by crystallographic orbit.

ewald_real_terms.csv.gz: every real-space contribution at alpha=0.3 A^-1 for all
five Cu observation sites and six charge scenarios. delta = source atom + image
translation - observation. The image integers multiply the cell lengths. Self
images at r=0 are excluded. The saved uncharged radial A and B generate the term
z_j[A r_a r_b+B delta_ab]. All six independent tensor components include z_j.

ewald_reciprocal_terms.csv.gz: all nonzero reciprocal vectors retained at the
same alpha, including both signs, in each site/scenario sum. k_i=2*pi*h_i/a_i.
The structure factor uses sum_j z_j exp(-i k dot r_j). charge_cos_at_observation
= Re[exp(+i k dot r_obs) structure_factor]. The coefficient already includes
-4*pi/V * exp[-k^2/(4 alpha^2)]/k^2; multiplying by charge_cos*k_a*k_b gives each
saved tensor component. Individual atomic inputs for the structure factor are
fully specified by expanded_atoms.csv and charge_scenarios.csv.

ewald_subtotals.csv: real sum, reciprocal sum, positive self-Gaussian correction,
total before trace projection, negligible roundoff trace projection, and total.
The self correction is +4*alpha^3*z_self/(3*sqrt(pi)) times the identity.
The raw trace is tested before applying the numerical traceless projection.
ewald_alpha_convergence.csv and ewald_convergence_tensors.json save complete
tensor results for alpha=0.2,0.3,0.4 A^-1, R_c=6/alpha, K_c=12*alpha.
tensors_and_eigensystems.json saves all eigenvalues, principal axes, invariants,
bare SI values and relative NQR proxies. Eigenvector signs are arbitrary; signs
are canonicalized by making the largest-magnitude Cartesian component positive.

direct_radial_shell_basis_by_source_orbit.csv.gz contains independent unscreened
Coulomb-kernel sums with unit source charge, grouped by observation site, source
orbit and radial shell. Multiply by the charge_scenarios.csv trial z to obtain
each scenario. Lower boundaries are excluded and upper boundaries included.
direct_radial_shell_scenarios.csv includes weighted shell and cumulative tensors.
direct_radius_convergence.csv compares spheres up to 100 A against Ewald and the
archived direct-sum results. Finite spherical sums converge slowly and need not
equal the periodic reference at any one finite cutoff.

frequency_scaling_all_144_assignments.csv retains ALL 24 high-site assignments
for EACH of six charge scenarios. No scenario is selected as experimentally
established. P=|G_ZZ|sqrt(1+eta_pc^2/3). A single positive empirical scale is fit
to four high-frequency 63Cu NQR values with unweighted least squares. The same
scale times P_Cu3 gives a conditional low-frequency estimate. Its calibration,
per-group residuals, fractional residuals, and discrepancies from two conditional
NMR eta branches are explicit. Eta discrepancies do not enter the frequency-only
ranking. frequency_scaling.json records assumptions and all scenario minima.
The resulting ranges are MODEL SENSITIVITY, not statistical uncertainty.

numerical_validation.json describes internal and archived-result checks.
saved_term_replay_validation.json is created by the separate validator, which
reads the gzip CSV files and independently re-sums the saved terms.
MANIFEST.json lists SHA-256 hashes of result files and numerical script inputs.
"""
    (OUT / "SCHEMA.txt").write_text(schema)
    manifest = {
        "inputs": [{"path": str(p.relative_to(BASE)), "sha256": hashlib.sha256(p.read_bytes()).hexdigest()} for p in [CIF, ARCHIVED, FIT_TABLE, Path(__file__)]],
        "results": [{"path": p.name, "bytes": p.stat().st_size, "sha256": hashlib.sha256(p.read_bytes()).hexdigest()}
                    for p in sorted(OUT.iterdir()) if p.is_file() and p.name not in {"MANIFEST.json", "saved_term_replay_validation.json"}],
    }
    save_json("MANIFEST.json", manifest)
    print(json.dumps({"output_directory": str(OUT), "validation": validation,
                      "frequency_best_by_scenario": [{key: row[key] for key in ["scenario", "assignment_s1_s2_s3_s4", "high_groups_RMSE_MHz", "scaled_Cu3_MHz"]} for row in best]}, indent=2))


if __name__ == "__main__":
    main()
