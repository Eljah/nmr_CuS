# Independent physics review of recovered spectra

Working review, 26 November 2026. Input files were read only. New descriptive calculations are in `shape_comparison.json` in this directory. No statistical significance claim is made.

## 1 What the plots support

The Origin project contains two separated spectral regions, approximately 74–78 MHz and 80–84 MHz, with sharp central peaks and weaker broad spectral intensity around them. The short- and long-repetition records retain many common spectral features but do not differ by a single scale factor.

This is consistent with a distribution of repetition-time responses among contributions to the detected signal. It does not, without pulse-program metadata, normalization and background control, prove a distribution of intrinsic sample T1 values or identify particular crystallographic sites.

The narrow maxima could contain metallic copper from the coil/leads. Remarks.pdf explicitly documents a narrow-plus-broad copper-wire/sample overlap for the earlier field/run. That note makes a hardware contribution a concrete possibility; it does not establish that every narrow peak in December comes from hardware. Conversely, a narrow peak cannot be declared a new anilite site just because it is strong.

The existence of broad intensity is directly visible. Its microscopic assignment to a quadrupolar powder pattern is a model to test. Field disorder, shift anisotropy, multiple sites, multiple phases, receiver response and processing can also contribute. A broad-plus-narrow decomposition is not unique from the line shape alone.

## 2 Export integrity issue found independently

`opj_Spec05122016.csv` contains 256 rows but 255 unique frequency values. The two rows at 81.4 MHz have the same file identifiers (425 and 426) but different intensity values:

| origin_row | IntQ (short) | IntQ1 (long) |
|---|---:|---:|
| 188 | 50.040102 | 47.737420 |
| 189 | 511.395946 | 3017.037347 |

The low first row produces an artificial-looking sharp dip in the plotted frequency-stepped data. Its cause is unknown. Do not average or silently choose one. Both rows were excluded from the descriptive shape comparison below. Keep both in the provenance export and flag the conflict in analysis.

## 3 Positions and isotope interpretation

The frequency-stepped table has its lower-band maximum at 75.90 MHz and higher-band maximum at 81.35 MHz in both repetition records. In the long-repetition series 81.30 and 81.35 MHz are nearly tied, so wording such as 81.3 MHz is adequate, but 81.30 should not be stated as the exact table maximum.

The processed Fourier spectra retained in Origin have maxima:

| Sheet | Lower region maximum | Higher region maximum |
|---|---:|---:|
| spec20ms | 75.942517 MHz | 81.345654 MHz |
| spec500ms | 75.940991 MHz | 81.345654 MHz |

These are maxima of stored processed real-channel curves, not peak-fit estimates or independently calibrated resonance frequencies. Frequency resolution of the curve does not imply matching absolute accuracy. The ratio of the two long-record peak positions is approximately 1.07117, consistent with a Cu isotope pair under the conventional resonance-reference ratio. It cannot distinguish copper metal in the probe from copper in the sample; both contain the same isotopes.

Important source-consistency check: the IUPAC 2001 tabulated bare gamma values reported by the literature agent give a ratio about 1.06926, while its reference-frequency Xi values give 1.07121. The literature agent verified a primary correction: Locher P. R. (1974), Phys. Rev. B 10, 801, DOI 10.1103/PhysRevB.10.801, reports gamma63/gamma65=0.9335232(6), consistent with the approximately 1.07121 inverse ratio. Its abstract explicitly corrects an often-quoted earlier determination. Do not combine the inconsistent tabulated old bare gamma pair for isotope calibration. Nuclear magnetic moment ratios, resonance reference shielding and absolute Knight-shift conventions must remain explicit.

## 4 Descriptive repetition-time comparison

For identical stored frequency points, excluding both ambiguous 81.4 MHz rows, ratios of summed intensities are:

| Frequency interval | Number of points | Sum(long)/sum(short) |
|---|---:|---:|
| 75.20–75.70 MHz, broad plateau left of central peak | 11 | 5.728 |
| 76.10–76.40 MHz, broad plateau right of central peak | 7 | 4.708 |
| 75.85–75.95 MHz, narrow central feature | 3 | 8.612 |
| 74.90–75.10 MHz, lower edge feature | 5 | 5.037 |
| 76.40–76.50 MHz, upper edge feature | 3 | 4.117 |
| 81.25–81.35 MHz, higher narrow feature | 3 | 8.932 |

The first three low-frequency-region ratios become approximately 5.91, 4.76 and 8.73 after subtracting a simple linearly interpolated baseline defined by median levels at 72.05–73.5 and 78–79.5 MHz. Thus their difference is not removed by this elementary baseline treatment. That baseline is only a sensitivity check; the broad tails and receiver background have not been independently established.

In the flat 72.05–73.5 MHz interval, a robust adjacent-point-difference scatter proxy is approximately 8 a.u. for the short record and 10 a.u. for the long record (median absolute deviation of successive differences, divided by 0.67449 and sqrt(2)). This is not an experimental uncertainty estimate: neighboring measurements are not replicate acquisitions and may have correlation or smooth background drift. Nevertheless, the differences between scaled spectra are far larger than the visible local point-to-point jitter. It is defensible to state that frequency-dependent repetition response is resolved in the stored records; it is not defensible to quote a p-value or a site-specific T1 distribution.

As an additional nonessential descriptive metric, raw area normalization over the full recorded interval gives an integrated long/short ratio of 5.898 and an L1 distance of 0.208 between normalized shapes (total variation 0.104). These values depend on integration bounds, baseline and unresolved components and should not be promoted to material parameters.

Even if the delays are confirmed as 20 and 500 ms, two repetition times are not a recovery curve. Under an assumed single exponential, ideal complete saturation, known normalization and no other effects, a ratio can be mapped to an apparent time constant, but those assumptions are unverified here. For quadrupolar I=3/2 selective excitation/recovery is generally multiexponential. No numerical intrinsic T1 should be published from these two curves alone.

## 5 July sweep

The audited July plot contains response near 48 MHz and 60 MHz and variation with the parameter named Satur. The series stored as 500000 is only supported by distinct trace files through 58.0 MHz; its later carried-forward values must remain excluded. Therefore the central response near 60.1 MHz has only the 100000 and 10000 records.

Do not call Satur a saturation-recovery delay with specific units until the pulse program is verified. Shape/intensity changes can be described as dependence on the stored acquisition parameter. The three stored setting values are insufficient by themselves to establish a T1 experiment. The modest response difference around 48 MHz and much stronger difference near 60 MHz is consistent with contributions having distinct acquisition responses; this is not a phase/site assignment.

The signals near 48 and 60 MHz must not be paired as the two copper isotopes (their ratio is about 1.25, not about 1.071). Nor should the 60-MHz narrow maximum be automatically assigned to anilite in view of the earlier explicit copper-wire observation near that frequency.

Combining July and December as a field-dependent sample study requires run-level fields, temperatures, the same specimen history and the same calibration/reference convention. Similar labels alone are insufficient.

## 6 Five-site Hamiltonian fit with historical NQR constraints

A five-site structural model is a valid hypothesis if the anilite structure and actual phase are established. It is not justified to equate the four historical 63Cu NQR frequencies 11.90, 15.21, 16.26 and 16.81 MHz with five uniquely identified Cu sites. Possible overlaps, weak/undetected lines, temperature differences and phase impurities remain to be tested. The 14.89 MHz CuS line must enter as a separate potential phase/background component, not as the fifth anilite frequency.

For I=3/2 a zero-field NQR line constrains

nu_NQR,j = |nu_Q,j| sqrt(1 + eta_j^2/3).

It does not independently determine nu_Q and eta. A four-line set thus gives four scalar constraints, not a full five-site EFG tensor determination, and does not establish the signs of VZZ. Historical NQR frequencies must be compared at compatible temperatures or with an explicit allowance for temperature shifts.

The useful global constraints in a new fit would be:

1. Common electronic shift tensors for both Cu isotopes to the appropriate approximation; isotope-specific nuclear gyromagnetic and quadrupole parameters from one consistent reference.
2. The same local EFG for corresponding isotope sites; nu_Q scales with Q.
3. A shared physical site population model fixed or constrained by validated crystallographic multiplicities, with separate frequency-dependent excitation/detection and relaxation factors where necessary.
4. A separately modeled narrow copper hardware contribution, preferably measured with an empty probe or a noncopper winding; hardware peak amplitudes must not be forced into the site counts.
5. A separately modeled CuS impurity contribution if supported by XRD and spectroscopy.
6. The full Hamiltonian and converged orientational averaging, with an explicit model comparison among fewer-site, five-site, and additional-phase hypotheses.

At a single field, each site can otherwise have nu_Q, eta, several shift components, orientation angles, linewidth/distribution parameters and an effective amplitude. Allowing all of these to vary creates a highly underdetermined model even though a visually excellent fit may be possible. Two isotope bands help by imposing isotope scaling but do not automatically remove the degeneracies. Two repetition times add signal contrast but not independent unconstrained site identities.

Before claiming quantitative site-resolved tensors, prefer additional fields (magnetic shifts and second-order quadrupole shifts scale differently), a temperature series, proper recovery/echo-decay experiments, and a measured hardware background. Publish profile-likelihood/bootstrap sensitivity or a comparable identifiability analysis only once a justified noise model and adequate records exist. Do not convert an unvalidated best-fit K into valence, spin susceptibility or charge disproportionation.

## Suggested honest manuscript wording

“The recovered spectra display broad spectral intensity with superposed narrow features in two frequency regions. Their relative intensities depend on the pulse-sequence repetition setting, and the two spectra cannot be superposed by a single multiplicative factor. This behavior demonstrates heterogeneous acquisition response among the detected spectral contributions. A microscopic decomposition requires explicit treatment of the copper background, pulse excitation and relaxation weighting. The present records do not by themselves determine site-resolved longitudinal relaxation times.”

“Historical NQR frequencies provide constraints for a joint isotope-dependent Hamiltonian analysis, conditional on the revised phase attribution. For spin 3/2 they constrain a combination of the quadrupole coupling and asymmetry parameter, rather than either parameter separately. A five-site structural model should therefore be treated as a constrained hypothesis whose uniqueness must be assessed.”
