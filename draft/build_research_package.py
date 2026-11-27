"""Assemble the revision-2 manuscript and reproducible analysis package."""
from pathlib import Path
import hashlib
import json
import zipfile

BASE = Path(__file__).resolve().parent.parent
OUT = BASE / "output" / "Anilite_NMR_research_package.zip"

paths = [
    "output/Anilite_NMR_draft_ru.docx",
    "output/Anilite_NMR_spectra.png",
    "output/Anilite_NMR_powder_fit.png",
    "output/Anilite_NMR_T2.png",
    "output/README_PACKAGE_RU.txt",
    "draft/manuscript_v2.md",
    "draft/powder_results_v2.md",
    "draft/manuscript_v2_complete.md",
    "draft/build_document.py",
    "draft/build_research_package.py",
    "data_analysis/analyze_measurements.py",
    "data_analysis/summarize_manuscript.py",
    "data_analysis/README.md",
    "data_analysis/NOTES_DATA_RU.txt",
    "data_analysis/forensics_report.md",
    "relaxation_analysis/analyze_t2.py",
    "new_evidence/thesis/THESIS_EVIDENCE_RU.txt",
    "new_evidence/thesis/thesis_T2_figure11.png",
    "new_evidence/papers/NQR_PUBLICATIONS_AND_COIL_RU.txt",
    "new_evidence/papers/manifest.json",
    "research_literature/anilite_amcsd0009376_transcribed.cif",
    "research_literature/reproduce_anilite_sites.py",
    "research_literature/anilite_site_multiplicities.csv",
    "research_literature/anilite_site_multiplicities.json",
]
originals = sorted(p for p in (BASE / "upload").iterdir() if p.is_file())
if len(originals) != 10:
    raise RuntimeError(f"Expected ten original attachments, got {len(originals)}")
paths.extend(p.relative_to(BASE).as_posix() for p in originals)
for directory, suffixes in [
    ("data_analysis/export", {".csv", ".json", ".png", ".pdf"}),
    ("relaxation_analysis/export", {".csv", ".json", ".png", ".pdf"}),
    ("powder_model", {".py", ".txt", ".csv", ".json", ".png", ".pdf"}),
]:
    paths.extend(
        p.relative_to(BASE).as_posix()
        for p in sorted((BASE / directory).iterdir())
        if p.is_file() and p.suffix in suffixes
    )
paths = sorted(set(paths))
manifest = {
    "revision": 2,
    "date": "2026-11-26",
    "description": "Author-supplied originals, manuscript, powder calculations, and echo-decay analysis.",
    "excluded": "Reproducible NPZ caches, temporary office profiles, rendered manuscript QA pages, and separately downloaded journal PDFs.",
    "files": [],
}
for relative in paths:
    path = BASE / relative
    data = path.read_bytes()
    manifest["files"].append({
        "path": relative,
        "bytes": len(data),
        "sha256": hashlib.sha256(data).hexdigest(),
    })
with zipfile.ZipFile(OUT, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=6) as z:
    for relative in paths:
        z.write(BASE / relative, relative)
    z.writestr(
        "PACKAGE_MANIFEST_SHA256.json",
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n",
    )
with zipfile.ZipFile(OUT) as z:
    failed = z.testzip()
    if failed is not None:
        raise RuntimeError(f"Archive integrity check failed for {failed}")
print(json.dumps({
    "zip": str(OUT),
    "revision": 2,
    "file_count_including_manifest": len(paths) + 1,
    "original_attachments": len(originals),
    "size_bytes": OUT.stat().st_size,
    "crc_check": "passed",
}, ensure_ascii=False))
