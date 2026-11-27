"""Reproduce Cu-site multiplicities and nearest S distances from the saved CIF.

Uses only Python's standard library. The eight Pnma operations are written
explicitly in the same order as the transcribed AMCSD CIF. This is deliberately
a narrow calculation for this record, not a general CIF parser.
"""

from pathlib import Path
from itertools import product
import csv
import json
import math
import re

ROOT = Path(__file__).resolve().parent
data = (ROOT / "anilite_amcsd0009376_transcribed.cif").read_text()
cell = [float(re.search(rf"_cell_length_{ax}\s+(\S+)", data).group(1))
        for ax in "abc"]
atoms = {m.group(1): tuple(float(m.group(i)) for i in (2, 3, 4))
         for m in re.finditer(
             r"^(Cu\d+|S\d+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)",
             data, re.M)}


def orbit(position):
    x, y, z = position
    operations = [
        (x, y, z), (x, .5-y, z), (-x, .5+y, -z),
        (.5-x, .5+y, .5+z), (.5+x, .5-y, .5-z),
        (.5+x, y, .5-z), (.5-x, -y, .5+z), (-x, -y, -z),
    ]
    return sorted({tuple(round(v % 1.0, 8) for v in op)
                   for op in operations})


orbits = {name: orbit(pos) for name, pos in atoms.items()}
sulfur = [(label, p) for label, sites in orbits.items()
          if label.startswith("S") for p in sites]
rows = []
for name, pos in atoms.items():
    if not name.startswith("Cu"):
        continue
    distances = []
    for sulfur_label, sulfur_position in sulfur:
        for shift in product((-1, 0, 1), repeat=3):
            delta = [(q + n - p) * length
                     for q, n, p, length in zip(sulfur_position, shift, pos, cell)]
            distances.append((math.sqrt(sum(v*v for v in delta)), sulfur_label))
    nearest = sorted(distances)[:6]
    rows.append({
        "site": name,
        "fractional_coordinates": pos,
        "multiplicity": len(orbits[name]),
        "relative_population": len(orbits[name]) // 4,
        "nearest_S_distances_A": [round(d, 6) for d, _ in nearest],
        "nearest_S_labels": [label for _, label in nearest],
    })

result = {
    "source": "AMCSD 0009376, Koto & Morimoto (1970), public mirror transcription",
    "cell_A": cell,
    "Cu_atoms_per_cell": sum(len(v) for k, v in orbits.items() if k.startswith("Cu")),
    "S_atoms_per_cell": sum(len(v) for k, v in orbits.items() if k.startswith("S")),
    "sites": rows,
    "note": "Distances are a geometric calculation, not an EFG/NQR assignment. "
            "Population ratios are not automatically observed NMR intensity ratios.",
}
assert result["Cu_atoms_per_cell"] == 28
assert result["S_atoms_per_cell"] == 16
assert [r["multiplicity"] for r in rows] == [4, 4, 4, 8, 8]
(ROOT / "anilite_site_multiplicities.json").write_text(json.dumps(result, indent=2) + "\n")
with (ROOT / "anilite_site_multiplicities.csv").open("w", newline="") as stream:
    writer = csv.writer(stream)
    writer.writerow(["site", "multiplicity", "relative_population",
                     "S_distance_1_A", "S_distance_2_A", "S_distance_3_A",
                     "S_distance_4_A", "S_distance_5_A", "S_distance_6_A"])
    for row in rows:
        writer.writerow([row["site"], row["multiplicity"], row["relative_population"],
                         *row["nearest_S_distances_A"]])
print(json.dumps(result, indent=2))
