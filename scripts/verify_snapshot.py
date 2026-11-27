"""Verify every archived research file against the publication snapshot.

Run from any directory: python scripts/verify_snapshot.py
Only Python's standard library is required. Files are read, never changed.
The manifest does not hash itself; its version is identified by the Git commit.
"""
from pathlib import Path
import hashlib
import json

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "provenance" / "SNAPSHOT_MANIFEST.json"


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main():
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    failures = []
    total_bytes = 0
    for entry in manifest["files"]:
        path = ROOT / entry["repository_path"]
        if not path.is_file():
            failures.append({"path": entry["repository_path"], "reason": "missing"})
            continue
        size = path.stat().st_size
        total_bytes += size
        if size != entry["bytes"] or sha256(path) != entry["sha256"]:
            failures.append({"path": entry["repository_path"], "reason": "bytes or SHA-256 differ"})
    print(json.dumps({
        "files_checked": len(manifest["files"]),
        "bytes_checked": total_bytes,
        "all_match": not failures,
        "failures": failures,
    }, ensure_ascii=False, indent=2))
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
