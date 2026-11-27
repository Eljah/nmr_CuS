"""Joint all-eta-zero control, preserving all archived inputs and fit results.

This tests local EFG axiality with measured NQR anchors. It is NOT a prediction
of the NQR frequencies or occupancies of the Donnay digenite structure.
"""
from pathlib import Path
import os
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
import sys
import json
import time
import hashlib
import functools
import numpy as np
from scipy.ndimage import gaussian_filter1d
from scipy.optimize import differential_evolution, minimize, least_squares, nnls

OUT = Path(__file__).resolve().parent
REPO = OUT.parent.parent
sys.path.insert(0, str(REPO / "powder_model"))
from powder_fit import (Fit, ROOT, NQR_SETS, GAMMA_RATIO, REL, STEP,
                        directions, stick_profile)

NQR_NAME = "thesis2009"
LOW = np.array([74.8] * 4 + [.025, .05, .05])
HIGH = np.array([76.8] * 4 + [.25, .8, .8])
BOUNDS = list(zip(LOW, HIGH))
STEPS = np.array([.0015] * 4 + [.001, .002, .002])
INPUT_PATHS = [ROOT / name for name in ["powder_fit.py", "joint_fit.py",
               "refine_exact_joint.py", "joint_results.json",
               "exact_joint_refinements.json"]]
INPUT_PATHS += [REPO / "data_analysis/export/opj_Spec05122016.csv"]


def hashes():
    return {str(p.relative_to(REPO)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in INPUT_PATHS}


def expand(q):
    return np.r_[np.zeros(4), q]


def write_json(name, data):
    (OUT / name).write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")


def boundary_flags(q):
    names = [f"nuZ_s{i}" for i in range(1, 5)] + ["sigma", "rho500", "rho20"]
    return [names[i] + (":lower" if abs(q[i]-LOW[i]) < 1e-5 else ":upper")
            for i in range(7) if min(abs(q[i]-LOW[i]), abs(q[i]-HIGH[i])) < 1e-5]


class Exact:
    def __init__(self, nu, nphi):
        self.dirs = directions(nu, nphi)
        self.fits = [Fit(np.empty(0), series, .12, False, True)
                     for series in ["500", "20"]]

    @functools.lru_cache(maxsize=512)
    def profiles(self, site, eta, centre):
        return np.array([stick_profile(eta, NQR_SETS[NQR_NAME][site, iso],
                         centre * (GAMMA_RATIO if iso else 1), self.dirs)
                         for iso in [0, 1]])

    def series(self, p, fit, rho):
        a = np.empty((len(fit.x), 4))
        for site in range(4):
            centre = p[4+site]
            raw = self.profiles(site, float(p[site]), float(centre))
            smooth = gaussian_filter1d(raw, p[8]/STEP, axis=1, mode="constant")
            a[:, site] = (np.interp(fit.x-centre, REL, smooth[0]) +
                         rho*np.interp(fit.x-centre*GAMMA_RATIO, REL, smooth[1]))
        ma = a[fit.mask]
        pa = ma - fit.Q @ (fit.Q.T @ ma)
        amps = nnls(pa, fit.py)[0]
        bg = np.linalg.lstsq(fit.base[fit.mask], fit.yy[fit.mask]-ma@amps,
                             rcond=None)[0]
        prediction = a@amps + fit.base@bg
        residual = (prediction-fit.yy)[fit.mask]
        return residual, prediction*fit.yscale, amps*fit.yscale, bg*fit.yscale

    def residual(self, q):
        p = expand(q)
        return np.concatenate([self.series(p, fit, p[9+i])[0]
                               for i, fit in enumerate(self.fits)])

    def jac(self, q):
        cols = []
        for i, step in enumerate(STEPS):
            lo, hi = q.copy(), q.copy()
            lo[i] = max(LOW[i], q[i]-step)
            hi[i] = min(HIGH[i], q[i]+step)
            cols.append((self.residual(hi)-self.residual(lo))/(hi[i]-lo[i]))
        return np.array(cols).T

    def record(self, p, prefix=None):
        r = {"parameters": p.tolist(), "eta": p[:4].tolist(),
             "nuZ63_MHz": p[4:8].tolist(), "sigma_MHz": float(p[8]),
             "rho500": float(p[9]), "rho20": float(p[10]),
             "nuQ63_MHz": (NQR_SETS[NQR_NAME][:,0]/np.sqrt(1+p[:4]**2/3)).tolist()}
        sse, n = 0., 0
        for i, fit in enumerate(self.fits):
            residual, pred, amps, bg = self.series(p, fit, p[9+i])
            e, count = float(residual@residual), int(fit.mask.sum())
            r[fit.series+"ms"] = {"normalized_SSE": e, "n": count,
                "normalized_RMS": float(np.sqrt(e/count)),
                "raw_RMS_au": float(np.sqrt(e/count)*fit.yscale),
                "normalization_au": float(fit.yscale),
                "amplitudes": amps.tolist(), "linear_correction": bg.tolist()}
            sse += e
            n += count
            if prefix:
                np.savetxt(OUT / f"{prefix}_{fit.series}ms.csv",
                    np.column_stack([fit.x, fit.y, pred, pred-fit.y, fit.mask]),
                    delimiter=",", header="frequency_MHz,observed,predicted,residual,used_in_fit",
                    comments="")
        r.update(SSE_joint=sse, n_joint=n, normalized_RMS_joint=float(np.sqrt(sse/n)))
        return r


def main():
    start = time.time()
    original_hashes = hashes()
    original = json.loads((ROOT / "exact_joint_refinements.json").read_text())
    lib = np.load(ROOT / "library_thesis2009_40x80.npz")["templates"]
    fits = [Fit(lib, s, .12, False, True) for s in ["500", "20"]]

    def objective(q):
        p = expand(q)
        return fits[0].solve(p[:10]) + fits[1].solve(np.r_[p[:9], p[10]])

    trials = []
    for i, seed in enumerate([11, 29, 67]):
        before = time.time()
        warm = np.array(original[i % 2]["parameters"])[4:]
        de = differential_evolution(objective, BOUNDS, seed=seed, x0=warm,
             maxiter=180, popsize=9, tol=5e-6, polish=False, workers=1)
        opt = minimize(objective, de.x, method="Nelder-Mead", bounds=BOUNDS,
             options={"maxiter":4000, "xatol":1e-7, "fatol":1e-10})
        rec = {"seed":seed, "q":opt.x.tolist(), "SSE_joint":float(opt.fun),
               "success":bool(opt.success), "message":str(opt.message),
               "de_SSE":float(de.fun), "de_nfev":int(de.nfev),
               "local_nfev":int(opt.nfev), "boundaries":boundary_flags(opt.x),
               "elapsed_seconds":time.time()-before}
        trials.append(rec)
        write_json("library_trials.json", trials)
        print("library", json.dumps(rec), flush=True)

    exact = Exact(80, 160)
    library_best = min(trials, key=lambda x:x["SSE_joint"])
    starts = [("best_library", np.array(library_best["q"]))]
    starts += [("projected_free_seed"+str(r["source_seed"]),
                np.array(r["parameters"])[4:]) for r in original]
    refined = []
    for label, q in starts:
        before = time.time()
        initial = float(exact.residual(q) @ exact.residual(q))
        opt = least_squares(exact.residual, q, jac=exact.jac,
            bounds=(LOW, HIGH), max_nfev=85, ftol=2e-6, xtol=2e-6,
            gtol=2e-6, x_scale="jac")
        rec = exact.record(expand(opt.x))
        rec.update(start_label=label, initial_SSE=initial,
            success=bool(opt.success), message=str(opt.message), nfev=int(opt.nfev),
            boundaries=boundary_flags(opt.x), elapsed_seconds=time.time()-before)
        refined.append(rec)
        write_json("exact_refinements_80x160.json", refined)
        print("exact80", json.dumps({k:rec[k] for k in
            ["start_label", "SSE_joint", "parameters", "boundaries", "elapsed_seconds"]}), flush=True)

    best = min(refined, key=lambda x:x["SSE_joint"])
    dense = Exact(160, 320)
    dense_control = dense.record(np.array(best["parameters"]), "axial_exact160")
    dense_free = []
    for r in original:
        free = dense.record(np.array(r["parameters"]),
                             "archived_free_eta_seed"+str(r["source_seed"])+"_exact160")
        free["archived_source_seed"] = r["source_seed"]
        free["archived_80x160_SSE"] = r["SSE_joint"]
        dense_free.append(free)
    free_best = min(dense_free, key=lambda x:x["SSE_joint"])
    result = {
        "purpose":"Local axial-EFG control with measured NQR; not a full Donnay structural prediction",
        "model":"All transitions, isotropic group shifts, common Gaussian sigma, two repetition series, NNLS group amplitudes, free signed linear correction",
        "NQR_set":NQR_NAME, "NQR_frequencies_MHz":NQR_SETS[NQR_NAME].tolist(),
        "fixed_eta":[0.,0.,0.,0.], "bounds":BOUNDS,
        "mask_halfwidth_MHz":.12, "masked_coil_centres_MHz":[75.941,81.346],
        "excluded_frequency_MHz":81.4,
        "objective":"Sum of normalized unweighted SSE; each series divided by its unmasked maximum; not a chi-square",
        "optimization":"3 independent DE+Nelder-Mead library searches; 3 direct 80x160 refinements; common 160x320 final evaluation",
        "axial_80x160_best":best, "axial_160x320":dense_control,
        "archived_free_eta_160x320":dense_free,
        "SSE_ratio_axial_over_best_free":dense_control["SSE_joint"]/free_best["SSE_joint"],
        "RMS_ratio_axial_over_best_free":dense_control["normalized_RMS_joint"]/free_best["normalized_RMS_joint"],
        "caveats":["No fitted noise model or statistical confidence/p-value",
            "Only a test of simultaneous eta=0, not separate confidence intervals for eta_s3 or eta_s4",
            "Shift anisotropy, pulse excitation, texture, additional phases and group-specific broadening remain unmodelled",
            "Free-reference parameters are archived exact 80x160 optima, evaluated here on the same dense grid without reoptimization",
            "Signed linear correction may compensate for response-model deficiencies"],
        "input_sha256_before":original_hashes,
        "input_sha256_after":hashes(),
        "archived_inputs_unchanged":original_hashes==hashes(),
        "elapsed_seconds":time.time()-start}
    write_json("RESULTS.json", result)
    print("FINAL", json.dumps({k:result[k] for k in
        ["SSE_ratio_axial_over_best_free", "RMS_ratio_axial_over_best_free",
         "archived_inputs_unchanged", "elapsed_seconds"]}), flush=True)
    print("AXIAL", dense_control["SSE_joint"], dense_control["normalized_RMS_joint"], flush=True)
    print("FREE", free_best["SSE_joint"], free_best["normalized_RMS_joint"], flush=True)


if __name__ == "__main__":
    main()
