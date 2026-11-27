# Structural input and conventions

The calculation evaluates one explicit historical structural hypothesis. It is
not a structural refinement of the NMR sample and does not convert the sample's
X-ray identification as anilite into a digenite identification.

## Sources

- Donnay G., Donnay J.D.H., Kullerud G. (1958). *Crystal and twin structure of
  digenite, Cu9S5.* American Mineralogist **43**, 228–242.
  [Publisher PDF](https://msaweb.org/AmMin/AM43/AM43_228.pdf).
- [COD 9000079](https://www.crystallography.net/cod/9000079.html), linked to
  AMCSD 0000082. The original CIF was downloaded on 2026-11-27 from
  `https://www.crystallography.net/cod/9000079.cif?download=1`.
  The downloaded revision is 294031, dated 2024-08-07.
- `COD_9000079_original.cif` is the unmodified 3505-byte download, SHA-256
  `dbba101beb32b46ab35978bb006a4f228916d6b402f7696e16d387d6d924780b`.
  Its original attribution notice is retained.

## Why the original file must not be expanded naively

The downloaded CIF specifies R-3m in the primitive rhombohedral setting and
contains both +x and -x partners in the atom loop. Expanding each row as an
independent crystallographic orbit duplicates the non-centrospecial atoms.
The file's formula field is consequently `Cu17 S9`. It is not the intended
Cu9S5 motif. Removing spatially coincident duplicates gives 9 Cu and 5 S.

`donnay1958_unique_orbits_R-3m.cif` lists the five Cu and three S orbits once,
retains the original symmetry operations, and generates 14 unique atoms.
It is explicitly labelled a derived model input. The calculation itself
constructs the same 14 atoms directly from `model.properties`; it does not
depend on a CIF parser's duplicate-handling behaviour.

The model uses the original angle 13°56′ = 13.933333333333334°, whereas COD
rounds it to 13.93°. Both choices are evaluated in the sensitivity output.
The lattice constant is 16.16 Å. No precise 4.2–5 K atomic coordinates,
displacement tensors, or temperature-dependent lattice constants are assumed.

The labels DgA–DgE are local to this calculation. In the downloaded raw CIF
they correspond to Cu1, Cu2/Cu3, Cu4/Cu5, Cu6/Cu7, Cu8/Cu9 respectively.
They do not denote Cu1–Cu5 of the anilite model or the measured s1–s4 groups.

## Experimental comparison input

The four comparison frequencies are 11.97, 15.27, 16.18 and 16.96 MHz for
63Cu at 4.2 K, from the thesis-based dataset already archived in this
repository. They are also the `thesis2009` constraints of the existing
`powder_model` analysis. These are not a newly extracted spectrum.

The inherited comparison scale 38.2916669626624 MHz per raw EFG proxy is
copied from the baseline anilite point-charge fit in
`theory/results/efg/frequency_scaling_best_by_scenario.csv`. Its use is
explicitly conditional and is not a calculation of a Sternheimer factor.

The historical moments Q(63Cu) = -0.211 b and Q(65Cu) = -0.195 b are kept
explicit in `model.properties`. See R.M. Sternheimer, Phys. Rev. **164**, 10
(1967), [DOI](https://doi.org/10.1103/PhysRev.164.10), and the discussion of
updated values in P. Pyykkö, Molecular Physics **116**, 1328–1338 (2018),
[DOI](https://doi.org/10.1080/00268976.2018.1426131). This convention is a
controllable input, not an assertion that those two rounded moments are the
unique modern standard. In particular, using Q(63Cu) = -0.220 b multiplies
all bare 63Cu frequencies by 220/211 without changing any relative EFG or eta.
