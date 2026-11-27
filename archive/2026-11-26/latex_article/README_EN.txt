ANILITE NMR: ENGLISH LATEX MANUSCRIPT AND SUPPORTING CALCULATIONS
Working manuscript, 26 November 2026

1. MANUSCRIPT FILES

main.tex                  Main English article.
structure_priors.tex      Included section on crystallographic constraints
                          and exploratory periodic point-charge calculations.
references.bib            Bibliography.
figures/                  Two vector figures used by the article.
Anilite_NMR_article_en.pdf Compiled article, when included in the release ZIP.

Compile from the directory containing main.tex with a complete LaTeX
installation:
pdflatex -interaction=nonstopmode -halt-on-error main.tex
bibtex main
pdflatex -interaction=nonstopmode -halt-on-error main.tex
pdflatex -interaction=nonstopmode -halt-on-error main.tex

The source uses a standard two-column A4 article layout and common packages:
fontenc, inputenc, lmodern, geometry, amsmath, amssymb, bm, graphicx, booktabs,
microtype, natbib, url, and hyperref. No author list or affiliation is invented.
The PDF was rendered with pdfTeX. Temporary build logs and environment-specific
TeX format/font-map repair files are not necessary for a normal TeX installation.

2. THREE DISTINCT NUMERICAL TASKS

A. Experimental powder NMR fit (unchanged from the preceding numerical analysis)

supporting_materials/powder_model/ contains the four-group NMR model, saved
solutions, validation results, and scripts. Four frequency-labeled NQR groups
s1--s4 constrain the quadrupolar terms. They are not assignments to the five
crystallographic Cu positions. The two principal eta branches remain conditional
fits of the 20 and 500 ms spectra; adding the new structural discussion has not
replaced or silently refitted these experimental results.

The final numerical table is FINAL_MODEL_COMPARISON.csv. The direct
central-transition evaluation is central_exact_validation.json. Other relevant
files are central_control_results.json, exact_joint_refinements.json,
joint_literature_sensitivity.json, and the four central_exact_seed*.csv files.
MODEL_ANALYSIS_RU.txt explains the numerical history and limitations. The saved
joint_results.json is an earlier optimization stage retained as an input for
reproduction, not a substitute for the final comparison table.

B. Exploratory crystallographic point-charge model (new supporting calculation)

supporting_materials/new_evidence/structure_model/ contains geometry_efg_proxy.py,
summarize_geometry_proxy.py, the full tensor results, convergence checks, and
concise structural-prior tables. It reads the separately preserved transcription
of AMCSD 0009376 under supporting_materials/research_literature/.

The periodic lattice tensor is evaluated by an Ewald sum, checked against changes
of the splitting parameter and finite-radius sums. The baseline neutral charges
are q_Cu=8/7 and q_S=-2. Five additional neutral charge-distribution scenarios test
model sensitivity. These charges are assumptions, not measured Cu oxidation
states. The reported eta_pc values are ratios of eigenvalues of this electrostatic
proxy. They omit the on-site electronic charge distribution, covalency, screening,
and antishielding. They are neither measured asymmetries nor electronic-structure
predictions of the complete nuclear EFG. The scenario spans are not confidence
intervals. No absolute NQR frequency is inferred from the proxy tensor.

The structural comparison is independent of the fitted NMR parameters and is
used to formulate tentative assignment tests. It does not establish an s1--s4 to
Cu1--Cu5 mapping. The coordinate reference is near room temperature, not a 5 K
refinement of the measured material.

C. Forward-model explanation of eta sensitivity (new illustration)

supporting_materials/new_evidence/eta_explanation/ contains eta_shape_demo.py and
its CSV, JSON, PDF, and PNG outputs. This is a pedagogical calculation, not an
additional measurement or fit. It compares exact-Hamiltonian central-transition
powder profiles for eta=0.0, 0.3, 0.6, and 0.9 at the same NQR frequency, 16.18 MHz,
and Zeeman frequency, 76 MHz. For each eta, nu_Q is changed to keep the NQR
constraint satisfied. Unit-area profiles and a fixed Gaussian width illustrate
why different asymmetries can be distinguished by line shape in principle.

3. INPUTS, PROVENANCE, AND RELATIVE PATHS

The compact package begins at main.tex. All calculation folders are preserved
beneath supporting_materials/ with their original relative arrangement:
  data_analysis/export/
  powder_model/
  relaxation_analysis/
  research_literature/
  new_evidence/structure_model/
  new_evidence/eta_explanation/
  new_evidence/structure/
  new_evidence/thesis/

Do not flatten these folders. The powder scripts locate data_analysis/export as
a sibling of powder_model; the relaxation script does the same. The structural
proxy locates research_literature two directory levels above its own folder.
The eta demonstration has no external data-file dependency.

The compact inputs are extracted numeric tables and acquisition metadata. The
original user attachments and the Origin/archive extraction stage remain in the
previously supplied full Anilite_NMR_research_package.zip and are not duplicated.
Reproducing calculations from the supplied CSV files does not require Origin or
the raw attachments. Repeating the extraction itself requires the full package.

Relevant source reports are provided under new_evidence/structure/ and
new_evidence/thesis/. They distinguish retrieved primary evidence, calculations
made here, and information confirmed by the experimenter.

4. REPRODUCING THE SUPPLIED CALCULATIONS

Commands below run from the compact ZIP root. The numerical work was tested with
Python 3.12, NumPy 2.3.5, SciPy 1.17.0, and Matplotlib 3.10.8. Use these or compatible
versions. The scripts use NumPy's trapezoid function; NumPy 2.x is suitable.

Rebuild echo-decay results:
python supporting_materials/relaxation_analysis/analyze_t2.py

Rebuild the illustrative profiles:
python supporting_materials/new_evidence/eta_explanation/eta_shape_demo.py

Rebuild periodic structural proxies and their concise table:
python supporting_materials/new_evidence/structure_model/geometry_efg_proxy.py
python supporting_materials/new_evidence/structure_model/summarize_geometry_proxy.py

Rebuild the final powder figure and comparison table from saved fit parameters:
python supporting_materials/powder_model/finalize_controls.py

Repeat the powder optimization, in dependency order:
python supporting_materials/powder_model/joint_fit.py
python supporting_materials/powder_model/refine_exact_joint.py
python supporting_materials/powder_model/central_control.py
python supporting_materials/powder_model/joint_literature_sensitivity.py
python supporting_materials/powder_model/finalize_controls.py
python supporting_materials/powder_model/extent_examples.py

This last block performs substantially more work and creates intermediate NPZ
profile libraries automatically. These reproducible caches are omitted from the
compact package. Saved JSON results suffice for the shorter validation/plotting
command above. An optimization rerun may produce small numerical differences.

Optional descriptive summaries and crystallographic site expansion:
python supporting_materials/data_analysis/summarize_manuscript.py
python supporting_materials/research_literature/reproduce_anilite_sites.py

To refresh article figures after rerunning the analysis, copy the generated
powder_model/powder_model_comparison.pdf and
relaxation_analysis/export/figure_t2_reanalysis.pdf to figures/ using their current
filenames, then rebuild main.tex. The article's rounded tables are explicit
LaTeX text and should be checked against the resulting CSV files after any change
to the physical model or fitting assumptions.

5. CONVENTIONS AND REMAINING EXPERIMENTAL DOCUMENTATION

The experimenter confirms the same unoriented, ground, immobilized powder for
NQR and NMR. The historical sequence is geerite -> digenite -> anilite after repeat
X-ray diffraction before NMR. The repeat diffractogram/refinement and full
acquisition protocol should accompany a journal submission.

The coil regions are masked, not described as independently measured and
subtracted background. Model-to-model spread of eta is not a statistical
confidence interval. The T2a time variable is interpreted in microseconds and
defined as TwoTau=2*tau; the old mS plot label is retained as inconsistent metadata.
T2e is an effective 1/e envelope time, and beta>1 is a descriptive shape parameter.
The full author/affiliation list must be supplied before journal submission.
