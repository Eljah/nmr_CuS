"""Build the current manuscript and theory supplement with a normal TeX installation.

Run from any directory: python scripts/build_pdfs.py
Requires pdflatex and bibtex on PATH. It does not install or download software.
"""
from pathlib import Path
import json
import shutil
import subprocess

ROOT = Path(__file__).resolve().parents[1]
DEST = ROOT / "papers"
DEST.mkdir(exist_ok=True)

def run(cwd, command):
    print(" ".join(command), flush=True)
    subprocess.run(command, cwd=cwd, check=True)

article = ROOT / "latex_article"
(article / "build").mkdir(exist_ok=True)
latex = ["pdflatex", "-interaction=nonstopmode", "-halt-on-error", "-output-directory=build"]
run(article, latex + ["main.tex"])
run(article, ["bibtex", "build/main"])
run(article, latex + ["main.tex"])
run(article, latex + ["main.tex"])

theory = ROOT / "theory"
(theory / "build").mkdir(exist_ok=True)
for _ in range(3):
    run(theory, latex + ["efg_derivation.tex"])

outputs = []
for src, name in [(article / "build/main.pdf", "Anilite_NMR_en.pdf"),
                  (theory / "build/efg_derivation.pdf", "Anilite_EFG_theory.pdf")]:
    shutil.copy2(src, DEST / name)
    outputs.append(str((DEST / name).relative_to(ROOT)))
print(json.dumps({"compiled": outputs}, indent=2))
