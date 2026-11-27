"""Resolve the axial-limit polar-quadrature risk in the separate control.

At eta=0 and isotropic shift, the spectrum and transverse-RF weight are
azimuth invariant. One azimuth with a finer cos(theta) grid is equivalent
to repeating the azimuthal samples. The equality is checked below.
"""
import json
import time
import numpy as np
from scipy.optimize import least_squares, differential_evolution
from axial_control import (OUT, Exact, LOW, HIGH, BOUNDS, expand, hashes,
                           boundary_flags, write_json)


def local(exact, q, label):
    t = time.time()
    opt = least_squares(exact.residual, q, jac=exact.jac, bounds=(LOW,HIGH),
         max_nfev=100, ftol=2e-7, xtol=2e-7, gtol=2e-7, x_scale="jac")
    rec = exact.record(expand(opt.x))
    rec.update(label=label, success=bool(opt.success), message=str(opt.message),
               nfev=int(opt.nfev), boundaries=boundary_flags(opt.x),
               elapsed_seconds=time.time()-t)
    print("fine_local", label, rec["SSE_joint"], rec["elapsed_seconds"], flush=True)
    return rec


def main():
    t = time.time()
    old = json.loads((OUT/"RESULTS.json").read_text())
    p = np.array(old["axial_80x160_best"]["parameters"])
    convergence = []
    for n in [80,160,320,640,1280]:
        r = Exact(n,1).record(p)
        convergence.append({"n_cos_theta":n, "n_phi":1,
                            "SSE_joint":r["SSE_joint"]})
    write_json("axial_polar_convergence.json", convergence)
    print("polar", json.dumps(convergence), flush=True)
    exact = Exact(640,1)
    q = p[4:]
    records = [local(exact,q,"best_library_then_80x160")]

    # One additional direct-Hamiltonian global search avoids depending solely
    # on the coarser translated template library to locate axial minima.
    global_start = time.time()
    def objective(q):
        rr = exact.residual(q)
        return float(rr@rr)
    def callback(x,convergence):
        return time.time()-global_start > 110
    de = differential_evolution(objective, BOUNDS, seed=101,
        x0=np.array(records[0]["parameters"])[4:], maxiter=100, popsize=7,
        tol=2e-5, polish=False, callback=callback, workers=1)
    global_meta = {"seed":101,"SSE_joint":float(de.fun), "nfev":int(de.nfev),
        "success":bool(de.success), "message":str(de.message),
        "elapsed_seconds":time.time()-global_start}
    print("fine_global",json.dumps(global_meta),flush=True)
    records.append(local(exact,de.x,"direct_global_seed101"))
    best = min(records,key=lambda r:r["SSE_joint"])
    final = Exact(1280,1).record(np.array(best["parameters"]),"axial_exact1280")
    free = min(old["archived_free_eta_160x320"],key=lambda r:r["SSE_joint"])
    result = {"purpose":old["purpose"],
      "axial_quadrature":"1280 uniform cos(theta) midpoints, one azimuth; symmetry equivalent for axial EFG and isotropic shift",
      "azimuth_equivalence_abs_SSE_error_at80":abs(convergence[0]["SSE_joint"]-old["axial_80x160_best"]["SSE_joint"]),
      "azimuth_equivalence_abs_SSE_error_at160":abs(convergence[1]["SSE_joint"]-old["axial_160x320"]["SSE_joint"]),
      "polar_convergence_at_fixed_80grid_parameters":convergence,
      "local_refinements_640x1":records,"direct_global_search":global_meta,
      "axial_final_1280x1":final,"archived_free_eta_best_160x320":free,
      "SSE_ratio_axial_over_best_free":final["SSE_joint"]/free["SSE_joint"],
      "RMS_ratio_axial_over_best_free":final["normalized_RMS_joint"]/free["normalized_RMS_joint"],
      "caveats":old["caveats"],"input_sha256_after":hashes(),
      "archived_inputs_unchanged":old["input_sha256_before"]==hashes(),
      "elapsed_seconds":time.time()-t}
    write_json("FINAL_CONTROL.json",result)
    print("FINAL_FINE",json.dumps({k:result[k] for k in
         ["SSE_ratio_axial_over_best_free","RMS_ratio_axial_over_best_free",
          "archived_inputs_unchanged","elapsed_seconds"]}),flush=True)
    print("AXIAL",final["SSE_joint"],final["normalized_RMS_joint"],flush=True)


if __name__ == "__main__":
    main()
