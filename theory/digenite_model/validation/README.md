# Independent numerical reference

`independent_ewald.py` is the separate NumPy/SciPy implementation used to
cross-check the Java lattice sums. Its original numerical outputs are retained
as `independent_ewald_saved.json` and `independent_six_charges_saved.json`.
The independent calculation uses a complex reciprocal-space structure factor
and NumPy eigensolver; the Java implementation uses explicit real cosine sums
and a Jacobi eigensolver.

Run with Python, NumPy and SciPy:

```bash
python theory/digenite_model/validation/independent_ewald.py --full-direct
```

It writes `independent_ewald_recomputed.json` next to the script. Without
`--full-direct`, only the smaller neutral-cell crosschecks are repeated.
The large n=80 direct sum uses substantially more memory than the Ewald model.
The script's main program reproduces the baseline, equal-Cu-charge result;
the Java program computes and checks all six scenarios.

The saved six-scenario JSON retains the original key
`Gzz_raw_e_per_Angstrom3`. This denotes the numerical lattice sum with charges
expressed in units of e; after e is factored out, G itself has units Å^-3.
Do not confuse G with Vzz/10^21: the conversion includes 1.43996455.

The independent calculation's electrostatic prefactor differs from the Java
input prefactor by less than 10^-9 relative. The comparison assertions use Gzz,
which is independent of this conversion and of nuclear Q.
