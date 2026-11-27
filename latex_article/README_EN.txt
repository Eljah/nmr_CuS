ANILITE NMR: CURRENT ENGLISH LATEX MANUSCRIPT
Revision:27November2026

main.tex is the current article. It includes structure_priors.tex and
low_frequency_evidence.tex, uses references.bib, and reads the two vector
figures under figures/. The compiled current PDF is ../papers/Anilite_NMR_en.pdf.

Compile from this directory with a complete TeX installation:
pdflatex -interaction=nonstopmode -halt-on-error -output-directory=build main.tex
bibtex build/main
pdflatex -interaction=nonstopmode -halt-on-error -output-directory=build main.tex
pdflatex -interaction=nonstopmode -halt-on-error -output-directory=build main.tex

Alternatively, run python scripts/build_pdfs.py from the repository root.
That command also builds the full theoretical supplement in theory/.

The repository root contains the original input files, extracted measurements,
powder_model/, relaxation_analysis/, research_literature/, new_evidence/,
theory/, and provenance/. Consult ../README.md and ../REPRODUCIBILITY.md for
commands and source-to-result relationships.

Changes in this revision:
- Exact archival evidence for the1.5–2.2MHz NQR response at4.2K in nominalCu1.1S.
- Explicit distinction between observed mixed-phase intensity and a candidateCu3
  contribution, and between the thesis text and the experimenter's later account
  of resonant-circuit limitations.
- A full separate analytic derivation and all numerical Ewald intermediate terms.
- Conditional common-scale comparisons across144 assignments, including model
  discrepancies rather than selecting an established crystallographic mapping.
- Original NQR T2-rate tables recovered from embedded Origin projects; s3 values
  at6–10K replace the earlier approximate figure reading. Covellite is separated.

The main powder-NMR fits have not been silently refitted. Four spectroscopic
groups s1–s4 remain distinct from the five crystallographicCu labels. The two
eta1 branches remain comparable conditional fits. Model variation is not an
experimental confidence interval. Coil regions are masked. The echo time is
TwoTau=2*tau and the reconstructed numeric time unit is microseconds.

The preceding26November manuscript is preserved in ../archive/2026-11-26/.
The older ZIPs under ../output/ preserve the versions delivered at that date;
their supporting_materials/ layout differs from the current full Git repository.

The author/affiliation list and repeat-diffraction documentation remain to be
completed for journal submission. All supplied experimental evidence and the
complete calculations are retained for review.
