#!/usr/bin/env python3
"""Replay saved Ewald term tables, independently of the numerical exporter."""
import csv
import gzip
import hashlib
import json
from collections import defaultdict
from pathlib import Path

import numpy as np


OUT = Path(__file__).resolve().parent / "results" / "efg"
COMPONENTS = ("aa", "bb", "cc", "ab", "ac", "bc")


def main():
    saved = {}
    with (OUT / "ewald_subtotals.csv").open(newline="") as f:
        for row in csv.DictReader(f):
            saved[row["scenario"], row["site"], row["component"]] = np.array([float(row["G_" + a]) for a in COMPONENTS])
    errors = {}
    counts = {}
    for filename, component in [("ewald_real_terms.csv.gz", "real_sum"), ("ewald_reciprocal_terms.csv.gz", "reciprocal_sum")]:
        accumulated = defaultdict(lambda: np.zeros(6))
        n = 0
        with gzip.open(OUT / filename, "rt", newline="") as f:
            for row in csv.DictReader(f):
                accumulated[row["scenario"], row["site"]] += np.array([float(row["G_" + a]) for a in COMPONENTS])
                n += 1
        counts[filename] = n
        errors[component] = max(float(np.max(abs(total - saved[scenario, site, component])))
                                for (scenario, site), total in accumulated.items())
    final_errors = []
    for scenario, site, component in saved:
        if component == "total":
            replay = sum(saved[scenario, site, p] for p in ("real_sum", "reciprocal_sum", "self_term", "roundoff_trace_projection"))
            final_errors.append(float(np.max(abs(replay - saved[scenario, site, "total"]))))
    errors["total_from_subtotals"] = max(final_errors)
    manifest = json.loads((OUT / "MANIFEST.json").read_text())
    hash_checks = {row["path"]: hashlib.sha256((OUT / row["path"]).read_bytes()).hexdigest() == row["sha256"]
                   for row in manifest["results"]}
    assert all(hash_checks.values())
    assert max(errors.values()) < 1e-11, errors
    report = {"all_checks_passed": True, "max_absolute_component_errors": errors,
              "replayed_term_counts": counts, "verified_result_file_hashes": len(hash_checks),
              "tolerance": 1e-11,
              "scope": "Numerical record replay only; this does not validate physical point charges or electronic EFG contributions."}
    (OUT / "saved_term_replay_validation.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
