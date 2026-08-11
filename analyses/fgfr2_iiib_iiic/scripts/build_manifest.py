from __future__ import annotations

import hashlib
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "checks/release_manifest.tsv"


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def main() -> None:
    paths = [path for path in ROOT.rglob("*") if path.is_file() and path != OUTPUT]
    rows = ["path\tsize_bytes\tsha256"]
    for path in sorted(paths):
        rows.append(f"{path.relative_to(ROOT)}\t{path.stat().st_size}\t{digest(path)}")
    OUTPUT.write_text("\n".join(rows) + "\n")
    print(f"Wrote {len(paths)} entries to {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
