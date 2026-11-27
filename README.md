# Copper NMR/NQR: anilite Cu7S4 and the historical CuS mixtures

Research archive and working manuscript, 27 November 2026. The same complete research snapshot is published on **`main` and `master`**.

## Read the current work

| Document | PDF | Source / evidence |
|---|---|---|
| English NMR manuscript | [Anilite_NMR_en.pdf](papers/Anilite_NMR_en.pdf) | [LaTeX](latex_article/main.tex), [bibliography](latex_article/references.bib) |
| Model EFG: complete theoretical derivation | [Anilite_EFG_theory.pdf](papers/Anilite_EFG_theory.pdf) | [LaTeX](theory/efg_derivation.tex), [Russian guide](theory/README_RU.md) |
| Temperature of the structural input | — | [Why the reference structure is not a 4.2–5 K refinement](provenance/STRUCTURAL_TEMPERATURE_RU.md) |
| Archival 1.5–2.2 MHz NQR signal | — | [Source audit with exact pages and quotations](new_evidence/thesis/LOW_FREQUENCY_SIGNAL_RU.md) |
| Original measurements and recovered embedded tables | — | [Measurement inventory](provenance/MEASUREMENT_INVENTORY_RU.md), [raw attachments](upload/), [Origin extraction](data_analysis/thesis_embedded/) |

Download the complete branch using GitHub's **Code → Download ZIP**, or clone it with:

```bash
git clone --branch master https://github.com/Eljah/nmr_CuS.git
```

This is a full research archive, including original binary measurements and earlier numerical profile libraries. Its size is consequently substantially larger than the PDFs alone.

## Scientific status

The experimenter identifies the material as anilite after repeat X-ray diffraction before the NMR measurements. Earlier work called the relevant component **geerite**, and subsequently **digenite**. These historical names are preserved in the original files and explicitly distinguished from the revised interpretation. The NQR/NMR material is a ground, immobilized powder without deliberate particle alignment.

Four high-frequency NQR-labeled groups constrain the powder NMR calculation. Labels **s1–s4 are spectral groups**, not an established assignment to the five crystallographic Cu sites. The weakest group's asymmetry remains ambiguous between comparably good model branches. Metallic copper in the probe is treated as an instrumental contribution.

The thesis explicitly reports a broad, overlapping low-frequency NQR response at **1.5–2.2 MHz, 4.2 K**, for a nominal **Cu1.1S mixed sample**. Covellite has a known response in this interval. Cu3 in the reference anilite structure is a candidate for a possible additional contribution because it has the smallest EFG proxy in all six tested charge models. This is a hypothesis, not a phase-resolved resonance assignment.

The model EFG is an **unscreened periodic point-charge calculation of the near-room-temperature reference structure**, not a 4.2–5 K structural calculation. The source catalogue labels the structural entry 293 K; the transcribed CIF itself has no temperature field. The temperature of the experimenter's repeat X-ray measurement is not specified in the supplied files. It omits the on-site electronic charge density, covalency, screening, and antishielding. Its numerical convergence does not validate those physical assumptions. A conditional common frequency scale puts Cu3 at about 1.503 MHz in the baseline model, while the six scenario minima span 1.40–3.52 MHz and the models do not jointly reproduce the measured frequencies and NMR asymmetries.

## Repository map

| Directory | Contents |
|---|---|
| `upload/` | All ten original user-supplied files, preserved byte for byte |
| `data_analysis/` | Archive extraction, raw time-domain files, exported integral tables, plots, provenance, and recovered embedded Origin objects |
| `powder_model/` | NMR forward models, saved profile libraries, fits, alternative solutions, validation, and final comparison tables |
| `relaxation_analysis/` | Complex spin-echo envelope analysis and numerical exports |
| `theory/` | Analytic EFG/NQR/powder derivation, symbolic checks, independently exported Ewald terms, replay validator, and frequency-scaling comparisons |
| `new_evidence/` | Thesis copies and extracts, author-supplied papers, structural evidence, the original point-charge calculation, and the explanatory powder-shape calculation |
| `research_literature/` | Bibliographic audits, structural coordinates, reproducible site expansion, and licensed reference material |
| `latex_article/` | Current English article, bibliography, vector figures, and build records |
| `papers/` | Current compiled PDFs and rendered page previews |
| `draft/`, `output/`, `pdf_analysis/` | Earlier drafts, deliverables, archives, and document-analysis records; consult the current `papers/` PDFs for the revised results |
| `archive/` | The 26 November manuscript and earlier rendering checks |
| `provenance/` | Measurement inventory, publication scope, and checksums for the published files |

The revised NQR relaxation comparison uses actual numeric tables recovered from embedded Origin projects. For **63Cu s3 at 6–10 K**, the archived rates imply **T2 ≈ 0.121–0.141 ms**. The separate covellite entry at 4.2 K is identified explicitly. The older figure-read estimates remain in historical drafts for provenance.

## Reproduce

The full [publication manifest](provenance/SNAPSHOT_MANIFEST.json) records the size, source path, and SHA-256 of every published research file. Verify it with `python scripts/verify_snapshot.py`.

See [REPRODUCIBILITY.md](REPRODUCIBILITY.md) for dependencies and commands. A short mathematical audit is:

```bash
python theory/symbolic_algebra.py
python theory/validate_efg_intermediates.py
```

The symbolic script saves the algebraic intermediate expressions. The second command reconstructs every Ewald subtotal from the saved individual terms. Recompute the full intermediate numerical export with:

```bash
python theory/export_efg_intermediates.py
```

Build both current PDFs with a complete TeX installation:

```bash
python scripts/build_pdfs.py
```

## Provenance and publication scope

The original uploaded measurements and documents are preserved verbatim. Extracted tables retain original row indices, source object paths, and hashes. A text cell containing a comma decimal is retained as text until an explicit interpretation is made; graph digitization is not substituted for recovered numeric data.

See [provenance/PUBLISHING_SCOPE.md](provenance/PUBLISHING_SCOPE.md) for the treatment of runtime files and downloaded publication copies. Original rights and source acknowledgments remain with the supplied reference material; this archival import does not assert a new repository-wide licence.
