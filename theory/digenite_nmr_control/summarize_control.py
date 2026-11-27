"""Produce compact comparison tables from completed control calculations."""
from pathlib import Path
import csv
import json
import platform
import numpy as np
import scipy

OUT = Path(__file__).resolve().parent
r = json.loads((OUT/"FINAL_CONTROL.json").read_text())
models = [
    ("free_eta_archived_seed67", r["archived_free_eta_best_160x320"], 11,
     "160x320 octant; archived parameters, amplitudes and baseline re-solved"),
    ("all_eta_zero", r["axial_final_1280x1"], 7,
     "1280 cos(theta) midpoints; one azimuth by axial symmetry")]
rows, parameter_rows, compact_models = [], [], []
for name, model, nonlin, grid in models:
    row = {"model":name, "nonlinear_parameters":nonlin,
           "linear_parameters":12, "nominal_total_parameters":nonlin+12,
           "n_500ms":model["500ms"]["n"], "n_20ms":model["20ms"]["n"],
           "n_joint":model["n_joint"],
           "SSE_500ms":model["500ms"]["normalized_SSE"],
           "SSE_20ms":model["20ms"]["normalized_SSE"],
           "SSE_joint":model["SSE_joint"],
           "normalized_RMS_500ms":model["500ms"]["normalized_RMS"],
           "normalized_RMS_20ms":model["20ms"]["normalized_RMS"],
           "normalized_RMS_joint":model["normalized_RMS_joint"],
           "raw_RMS_500ms_au":model["500ms"]["raw_RMS_au"],
           "raw_RMS_20ms_au":model["20ms"]["raw_RMS_au"],
           "sigma_MHz":model["sigma_MHz"], "rho500":model["rho500"],
           "rho20":model["rho20"], "evaluation_grid":grid}
    rows.append(row)
    compact_models.append({**row, "eta":model["eta"],
        "nuZ63_MHz":model["nuZ63_MHz"], "nuQ63_MHz":model["nuQ63_MHz"],
        "normalization500_au":model["500ms"]["normalization_au"],
        "normalization20_au":model["20ms"]["normalization_au"],
        "amplitudes500":model["500ms"]["amplitudes"],
        "amplitudes20":model["20ms"]["amplitudes"],
        "linear_correction500":model["500ms"]["linear_correction"],
        "linear_correction20":model["20ms"]["linear_correction"]})
    for j in range(4):
        parameter_rows.append({"model":name,"group":"s"+str(j+1),
            "eta":model["eta"][j],"nuQ63_MHz":model["nuQ63_MHz"][j],
            "nuZ63_MHz":model["nuZ63_MHz"][j],
            "amplitude500":model["500ms"]["amplitudes"][j],
            "amplitude20":model["20ms"]["amplitudes"][j]})
for filename, records in [("comparison.csv",rows),("parameters.csv",parameter_rows)]:
    with (OUT/filename).open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(records[0]),lineterminator='\n')
        w.writeheader();w.writerows(records)
summary = {"purpose":r["purpose"],
    "criterion":"E = sum_r sum_i_in_mask [(pred_ri - observed_ri)/Y_r]^2, Y_r = maximum unmasked observed intensity in series r",
    "baseline":"b0_r + b1_r*(frequency_MHz-78)/6; both coefficients unrestricted in sign",
    "linear_parameters":"Four nonnegative group amplitudes and two signed baseline coefficients per series; re-solved at every nonlinear evaluation",
    "bounds":{"eta_free":[0,1],"eta_axial":[0,0],"nuZ63_MHz":[74.8,76.8],
              "sigma_MHz":[.025,.25],"rho65_per_series":[.05,.8]},
    "parameter_count_note":"Nominal counts: four eta are fixed in the axial model. NNLS activity and identifiability are not converted into statistical degrees of freedom.",
    "NQR63_MHz":[11.97,15.27,16.18,16.96],
    "NQR65_MHz":[11.07,14.13,14.97,15.69],
    "NQR_constraint":"Hard constraint nuQ = nuNQR/sqrt(1+eta^2/3); no NQR-error penalty term",
    "isotope_Zeeman_ratio":1.071210656566704,
    "masks":{"coil_centres_MHz":[75.941,81.346],"halfwidth_MHz":.12,
             "excluded_duplicate_frequency_MHz":81.4},
    "models":compact_models,
    "SSE_increase_percent":100*(r["SSE_ratio_axial_over_best_free"]-1),
    "RMS_increase_percent":100*(r["RMS_ratio_axial_over_best_free"]-1),
    "axial_global_search_success":r["direct_global_search"]["success"],
    "archived_inputs_unchanged":r["archived_inputs_unchanged"],
    "versions":{"python":platform.python_version(),"numpy":np.__version__,"scipy":scipy.__version__},
    "limits":r["caveats"]+["A direct numerical global optimizer and several starts do not prove a mathematical global minimum",
       "The joint all-eta-zero test does not establish separate nonzero eta for s3 or s4",
       "NQR frequencies are observed anchors, not frequencies predicted from a digenite atomic structure"]}
(OUT/"summary.json").write_text(json.dumps(summary,indent=2,allow_nan=False)+"\n")
print(json.dumps({"comparison":rows,"SSE_increase_percent":summary["SSE_increase_percent"],
                  "RMS_increase_percent":summary["RMS_increase_percent"],
                  "archived_inputs_unchanged":summary["archived_inputs_unchanged"]},indent=2))
