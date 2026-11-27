# Reproducing the anilite analysis

## Environment

The calculations were run with Python 3.12, NumPy 2.3.5, SciPy 1.17.0, Matplotlib 3.10.8, and SymPy 1.14.0. `requirements.txt` pins these versions and the optional Word-generation dependency. Use a separate virtual environment if installing them.

The numerical scripts use paths relative to their own source files. Keep the repository directory structure intact. The original files under `upload/` are inputs and are not modified by the analysis.

PDF compilation requires a normal TeX distribution with pdfLaTeX, BibTeX, Latin Modern, AMS packages, geometry, graphics, booktabs, natbib, microtype, and hyperref. Poppler's `pdftoppm` can render the resulting PDFs. Runtime-specific TeX format repair files are not part of this research archive.

## 1. Inspect original data before rerunning analyses

Verify the complete publication snapshot before regenerating any output:

```bash
python scripts/verify_snapshot.py
```

This reads the SHA-256 manifest in `provenance/SNAPSHOT_MANIFEST.json`. Regenerating outputs can intentionally change their bytes; the manifest describes the committed publication snapshot.

Read `provenance/MEASUREMENT_INVENTORY_RU.md` and `data_analysis/README.md`.

`Anilite.OPJ`, `Anilite14072016.zi`, `CuS.zi`, and the original attachment ZIP are distinct sources. The July CuS frequency series is not the later 72–85 MHz anilite integral table. The separately recorded echo-decay series has its own field and frequency.

`data_analysis/thesis_embedded/` contains recovered Origin projects, their numerical exports, source stream locations, and hashes. These tables were extracted from the supplied thesis and document files, rather than digitized from a plotted image. The extraction script is `data_analysis/extract_embedded_origin.py`; its command-line help documents the available inputs.

To repeat the original measurement-archive and OPJ extraction:

```bash
python data_analysis/analyze_measurements.py
python data_analysis/summarize_manuscript.py
```

These commands regenerate tables and descriptive plots. They do not establish phase assignments or convert uncalibrated integrals into site populations.

## 2. Analytic theory and full intermediate numerical records

```bash
python theory/symbolic_algebra.py
python theory/export_efg_intermediates.py
python theory/validate_efg_intermediates.py
```

The first command saves explicit algebraic expressions and identity checks for the Coulomb Hessian, screened Ewald Hessian, self correction, reciprocal derivative, mirror constraints, spin matrices, zero-field NQR eigenproblem, central-transition perturbation theory, and conditional least-squares scale.

The second saves the full expanded unit cell, six neutral trial-charge scenarios, every retained real-space and reciprocal-space Ewald term, self corrections, subtotal matrices, eigensystems, independent radial-shell sums, convergence checks, and all 144 frequency assignments. Details are in `theory/results/efg/SCHEMA.txt`.

The third re-reads the saved term tables and reproduces their subtotals and totals. It checks archive consistency, not the physical validity of the point-charge assumptions.

The earlier exploratory implementation is preserved separately:

```bash
python research_literature/reproduce_anilite_sites.py
python new_evidence/structure_model/geometry_efg_proxy.py
python new_evidence/structure_model/summarize_geometry_proxy.py
```

The new exporter was checked against these archived results without changing them. Exact decimal equality is not promised across arbitrary numerical-library versions; use the documented tolerances.

## 3. Powder NMR models

The final comparison is `powder_model/FINAL_MODEL_COMPARISON.csv`; the direct central-transition check is `powder_model/central_exact_validation.json`. Consult `powder_model/MODEL_ANALYSIS_RU.txt` before interpreting earlier fit files.

Rebuild the final figure from the retained solutions:

```bash
python powder_model/finalize_controls.py
```

The more expensive optimization sequence is:

```bash
python powder_model/joint_fit.py
python powder_model/refine_exact_joint.py
python powder_model/central_control.py
python powder_model/joint_literature_sensitivity.py
python powder_model/finalize_controls.py
python powder_model/extent_examples.py
```

The original profile-library NPZ files are included because this archive preserves the numerical intermediates. The two repetition-period datasets share quadrupolar parameters but have independent response amplitudes. Coil intervals are masked. Model-to-model variation is not a statistical confidence interval.

## 4. Echo decay and explanatory plots

```bash
python relaxation_analysis/analyze_t2.py
python new_evidence/eta_explanation/eta_shape_demo.py
```

The echo time is `TwoTau = 2*tau`; the numerical scale is interpreted in microseconds, with the legacy plot-unit discrepancy documented. `T2e` is an effective 1/e time for the complex-integral magnitude. The eta-shape figure is a forward-model illustration at fixed NQR frequency.

## 5. Documents

```bash
python scripts/build_pdfs.py
```

This builds the current English manuscript and full theory supplement into `papers/`. To refresh the manuscript's figures after changing an analysis, copy the generated `powder_model/powder_model_comparison.pdf` and `relaxation_analysis/export/figure_t2_reanalysis.pdf` to `latex_article/figures/`, then rebuild.

The rounded manuscript tables are explicit LaTeX text and must be reviewed if a model or dataset changes. `archive/2026-11-26/` and the older packages under `output/` preserve earlier results and wording; they do not override the revised current manuscript.
