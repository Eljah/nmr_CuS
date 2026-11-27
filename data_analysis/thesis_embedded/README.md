# Numeric Origin data embedded in the original thesis documents

These files were extracted from the supplied DOC, RTF and DOCX containers as data. No OLE application, Word macro, Origin formula or executable attachment was run. The curves were not digitized from images.

## Sources and equivalence

- `diploma5_01.opj`–`diploma5_05.opj`: five `Contents` streams in `upload/diploma_magister_geerite_2009_5.doc`.
- `diploma10_object_01_01.opj`–`diploma10_object_05_01.opj`: the same five projects in the alternate thesis RTF from `papers_20.zip`. Each is byte-identical to its corresponding main-thesis OPJ. The OLE1 and CFB wrappers are also retained.
- `geerite_to_digenite_oleObject5_01.opj`–`geerite_to_digenite_oleObject8_01.opj`: four projects in the working DOCX. Their graph formatting differs, but every saved numeric column matches the corresponding first four thesis projects. Their CFB containers are retained.

`embedded_origin_manifest.json` records the original document SHA-256, OLE stream or DOCX member name, RTF object byte offset where applicable, and all extracted project hashes. `numeric_inventory.json` lists column sizes, finite/missing counts, ranges and CSV hashes.

Each `numeric/<project>/` folder contains worksheet CSV files and `column_provenance.json`. The latter stores full original column names, cell types and values, byte offsets, source hash and the graph's first-layer curve references. `origin_row` in a CSV is the original worksheet row, starting from one. Missing numeric values remain empty/null; text cells remain text. The full binary project is retained even where later graph layers are not interpreted by the reader.

## Scientific interpretation

The five distinct table sets are high-frequency NQR spectra of four samples; NQR line frequencies versus temperature; T1T versus temperature; spin-lattice relaxation rates; and transverse relaxation rates. These tables include processed quantities and derived columns; they are not all instrument-level raw signals.

No numeric 1.5–2.2 MHz NQR spectrum was identified. The low-frequency evidence remains the original thesis Fig. 8. `SpectraS1`–`SpectraS4` identify sample spectra, not atomic sites. In the temperature-frequency project, the misleadingly named `Data1_freq` column is temperature in kelvin according to its graph axis, and the `*i65` columns exactly equal their `63Cu` counterparts divided by 1.081. They must not be counted as independent isotope measurements.

For the transverse relaxation graph, `Graph6` stores a rate axis in ms^-1; thesis Fig. 11 identifies 63Cu. Exact low-temperature values and curve-to-column mapping are in `provenance/NQR_T2_LOW_T_EXACT.csv/.json`. Spectral group s3 is not an assignment to the crystallographic atom Cu3.

## Reproduce

From the repository root, with Python and its standard library:

```bash
python3 data_analysis/extract_embedded_origin.py
python3 data_analysis/audit_measurement_inventory.py
```

The detailed Russian audit is `provenance/MEASUREMENT_INVENTORY_RU.md`; its machine-readable counterpart includes the ten original uploads, all archive-member hashes and the checks used to establish project equivalence.
