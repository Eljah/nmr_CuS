# Evidence audit of supplied PDFs

Reviewed 26 November 2026. This is a working audit for the manuscript author, not a publication-ready supplement. Original files were not changed.

## Remarks.pdf

Original: `/workspace/scratch/325ebb719515/upload/Remarks.pdf`. Two scanned handwritten pages, no legible date, signature, author, or title identifying anilite. PDF producer EPSON Scan; no usable creation-date metadata. Full page renders `remarks-1.png`, `remarks-2.png` in this directory.

### Reliable readings and evidential scope

| Location | Reading | Scope and limitation |
|---|---|---|
| p. 1, middle | H = 53.200 kOe, fixed (Pers. mode) | Nominal fixed field 5.32 T; not a certified field calibration. Decimal punctuation normalized. |
| p. 1, middle | 60.015 MHz, narrow line + broad line, 63Cu, from wires + sample | Explicit warning of copper-wire background near the sample spectrum. Does not establish the origin of every peak. |
| p. 1, middle | 64.26, 65Cu? | The question mark is part of the source; isotope assignment tentative. |
| p. 1, sketch | 60.081 MHz and 60.12 MHz | Sketch annotations, not fitted tabulated resonance positions. |
| p. 1, middle | T = 5 K; f_res = 64.1 MHz | Experimental setup note, no sample/run identifier. |
| p. 1, lower | f = 59.5 MHz; tau 25 or 30 microseconds; leave 30 microseconds | Optimization of echo detection; not an independently measured T2. |
| p. 1, lower | Period 200, 400, 50 ms; leave 200 ms | Repetition-period optimization. Another lower table spans 50 ms through 1.6 s. These values are not sufficient for T1 fitting without the pulse-program definition. |
| p. 1, lower | Pulse sequence label Cu20.seq; pulses 2.5 - 30 - 5 microseconds; 200 ms; Aver = 100 | A setup condition from an optimization sheet. It must not automatically be assigned to all later July data. The filename visually appears Cu20.seq; confirm against experiment software before publication. |
| p. 1, lower | Gain 30 + 34 dB was; set 6 + 34 dB; filter 1 MHz | Gain changed during optimization. Raw amplitudes across records may need gain normalization. |
| p. 2, top | T = 5.09 K; R_therm ≈ 14.93 kOhm | Temperature note; calibration and thermometry uncertainty not supplied. |
| p. 2, middle | Spectrum; Aver = 60; 57.5 MHz -> File 3; 57.6 -> 4; 57.7 -> 5; 57.8 -> 6; 57.9 -> 7; 58.0 -> 8; 58.1 -> 9; 58.2 -> 10; 58.3 -> 11 | Strong direct provenance link to the 51-point 57.5–62.5 MHz data set in the email/CuS archive, based on file numbering. |
| p. 2, lower | 61.0 -> File 38; last point 62.5 MHz -> File 53 | Confirms 0.1-MHz spacing and 51 frequency points (inclusive). |
| p. 2, middle-right | 1024 points; 0.16 microseconds | ADC sampling settings for that record series. |
| p. 2, upper-right | Ringing ends at 80 microseconds; at 57.5 MHz seems still signal | Source contains uncertainty about ringdown versus true signal; baseline/echo windows need verification. |
| p. 2, bottom | Reading ADF is incorrect! | Processing caution; not proof that the currently supplied text export is wrong. Reader/version/affected files unknown. |
| p. 2, bottom-right | Cu1.3S1 (circled); CuS + Cu1.8S below | Chemical composition notes; no explicit XRD, phase fractions, sample serial number or anilite name. This is NOT evidence of a phase-pure anilite specimen. |

Two gamma/frequency reference calculations also occur on p. 1: 60.50 MHz with gamma/(2pi)=1.13 MHz/kOe gives 53.5398 kOe; under a 7.2 T note, 63Cu f=79.618 MHz and 65Cu f=85.288 MHz are paired with 1.1058 and 1.1846 MHz/kOe, while 1.13 and 1.21 appear nearby with a question. These are working notes and internally unresolved reference choices, not a metrological source. A nuclear gamma, a reference compound resonance, a copper-metal wire Knight shift, and a true field calibration must be kept distinct. The manuscript should quote an authoritative isotope reference and actual calibration protocol instead of deriving calibration from these scribbles.

No statement about XRD reidentification from digenite to anilite was found. No value/definition of `Satur` was found. Do not transfer the field, temperature, pulse lengths, gains, or repetition periods to the July 14 series until run correspondence is established.

## Ilya.pdf

Original: `/workspace/scratch/325ebb719515/upload/Ilya.pdf`. 53 PDF pages, actually selected printed pages 9–18, 27–38, 57–87. No title page, author or references section is included. Metadata: LaTeX with hyperref; pdfTeX-1.40.12; created 2017-02-13, modified 2017-02-15. These are file metadata dates, not experimental dates. All pages inspected as contact sheets; key equations additionally viewed at readable resolution. Extracted text is `Ilya.txt`.

This is a methodological excerpt on sodium cobaltates, with 23Na and 59Co spectra and tables. No data on copper sulfides, anilite or digenite were found. Na/Co values and experimental conditions MUST NOT be described as measurements of anilite. Existing internal reference numbers [8] etc. cannot be imported without their missing bibliography.

### Transferable theory (proposed method, not a claim already carried out on anilite)

On PDF p. 34 / printed p. 68, Eq. (4.15) is a suitable starting Hamiltonian. In frequency units, with dimensionless spin operators and common principal axes of the shift and electric-field-gradient tensors:

H/h = -gamma_bar B0 [(1+KX) IX cos(alpha)sin(beta) + (1+KY) IY sin(alpha)sin(beta) + (1+KZ) IZ cos(beta)]
      + nuQ/6 [3 IZ^2 - I(I+1) + eta(IX^2 - IY^2)].

gamma_bar = gamma/(2pi), nuQ = 3eQVZZ/[2I(2I-1)h], eta=(VXX-VYY)/VZZ, and |VZZ|>=|VYY|>=|VXX|. For both stable Cu isotopes I=3/2; the static Hamiltonian is a 4x4 Hermitian matrix. Diagonalize it for each orientation rather than assuming small quadrupole/Zeeman ratio. Consider all six possible transitions and compute frequencies from eigenvalue differences and intensities from matrix elements of the transverse RF interaction. In a general model, use a full rotated shift tensor instead of silently imposing coincident axes.

On PDF pp. 36–37 / printed pp. 70–71, the third Euler angle (orientation of B1 around B0) is eliminated analytically: the orientationally averaged transition probability is half the sum of squared matrix elements for two orthogonal transverse B1 directions. This saves one angular loop; the gamma denoting this Euler angle must not be confused with the nuclear gyromagnetic ratio.

The powder sum carries a sin(beta) weight (uniform solid angle). For the symmetry assumptions in the excerpt, alpha and beta cover the first octant and a 90x90 grid was used. Do not treat 90x90 as universally converged or first-octant symmetry as valid for arbitrary misaligned tensors. Demonstrate convergence for the actual Cu spectrum.

On PDF pp. 38–39 / printed pp. 72–73, the frequency-swept spectrum at fixed B0 is accumulated as a histogram: for each orientation, deposit transition intensities into bins within a chosen line-broadening kernel. Diagonalization is repeated by orientation, not independently in every frequency bin. The excerpt uses triangular broadening f(dnu)=max(0,1-|dnu|/BW); BW is the FWHM of this triangular kernel and its half-base width. Alternative line shapes are model choices, not established facts about anilite. Distinct sites contribute with weights, and isotope contributions require natural abundances plus appropriate excitation/detection factors.

Do not infer site occupancy directly from spin-echo spectral areas without accounting for finite excitation bandwidth, gain, receiver response, T1 saturation and T2 decay. The excerpt's transition matrix elements describe a linear-response spectral model; they do not themselves simulate a full pulse sequence or account for frequency-dependent RF sensitivity. A proposed new method should be phrased as future work until actually implemented and fitted.

### Relevant theory cautions

- Use the clear Eq. (4.15), not a literal transcription of Eq. (4.3), which appears to have a dimensional typo in its second quadrupole term.
- The excerpt's K-tensor display on printed p. 78 also has a likely typo: the final diagonal entry inside the 1/2-prefactored matrix is printed as 2, whereas it should retain KZ for a rotation about Z. The subsequent definitions KZZ=1+KZ are usable.
- Dropping alpha when eta=0 requires all other tensors to have the same axial symmetry; eta=0 alone does not eliminate an azimuthal dependence from a nonaxial magnetic shift tensor.
- The PDF contains draft cross-references (chapter ??), reinforcing that it should not be cited as a finalized paper without identifying the parent work.

## Main missing evidence for the anilite article

1. Actual repeat XRD pattern and analysis/report, with date and sample identifier.
2. Chain linking the NQR specimen previously called digenite, the repeat XRD specimen and the later NMR specimen; storage/heating/grinding history if relevant.
3. Explicit run-level metadata for the July 14 anilite measurements, particularly B0, temperature, pulse program, Satur definition/units, gain and normalization.
4. Whether nominal Cu1.3S/CuS+Cu1.8S notes describe a separate comparison specimen.
5. Author and bibliographic identity of Ilya.pdf if its methodological derivation is to be cited.

The user's account is a legitimate author-provided historical statement, but these two PDFs alone do not independently establish the claimed reidentification or phase purity.
