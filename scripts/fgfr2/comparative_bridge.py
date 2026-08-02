"""Feed the validated FGFR2 dataset into the shared comparative builder.

The shared comparative layer works on one protein per species: a matrix row, an
alignment column and a species-ordered axis all assume that a species contributes
one sequence. FGFR2 contributes two, so the panel has to be reduced to one model
per species before it can enter that layer — and *which* model is a scientific
choice, not a filtering detail.

The reduction here is the stated primary-reference rule from
``scripts/fgfr2/coordinate_model.py``: each species is represented by its IIIc
model where it has one, otherwise by its IIIb model. Both models remain in the
coordinate index and both remain reachable in that species' Gallery scope; only the
comparative row is a single protein, because a comparative row can only be one.

The alignment is restricted the same way, to the reference proteins, keeping the
original columns. Restricting an alignment does not move the residues that stay in
it, so the columns still mean what they meant — some become all-gap, which the
shared code already handles.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

ROOT = Path(__file__).resolve().parents[2]
for _p in (ROOT / "scripts", ROOT / "scripts" / "shared_gene_analysis"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

from fgfr2 import coordinate_model as cm  # noqa: E402
from shared_gene_analysis.msa_coordinates import read_aligned_fasta  # noqa: E402

DERIVED = ROOT / "results" / "derived" / "example"


def write_reference_alignment(model_index: Dict[str, Any], out_path: Path) -> Path:
    """The full-length alignment restricted to each species' reference protein.

    Headers are rewritten to ``<protein_id> <gene>|<species_id>``, the form the
    shared alignment parser reads, so no FGFR2-specific header handling leaks into
    the shared code.
    """
    # Keyed on species *and* isoform, not on the protein id alone. One species in
    # the panel (Pongo abelii) has a single protein sequence filed under both
    # isoform slots because its isoform could not be resolved from sequence, so a
    # protein-id match would select that species twice and the alignment would no
    # longer have one sequence per species.
    reference = {(m["species_id"], m["isoform"]): m
                 for m in cm.primary_reference_models(model_index)}
    out_path.parent.mkdir(parents=True, exist_ok=True)
    kept: List[Tuple[str, str]] = []
    for header, seq in read_aligned_fasta(cm.FULL_MSA):
        parts = [p for p in header.split("|") if p]
        if len(parts) < 3:
            continue
        model = reference.get((parts[0], parts[1]))
        if model:
            kept.append((f"{parts[2]} {cm.GENE_SYMBOL}|{model['species_id']}", seq))
    with out_path.open("w", encoding="utf-8") as fh:
        for header, seq in kept:
            fh.write(f">{header}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i + 60] + "\n")
    return out_path


def build(derived: Path = DERIVED) -> Dict[str, Any]:
    """Build the FGFR2 comparative dataset from the primary reference models."""
    from shared_gene_analysis.comparative_dataset import build_comparative_dataset

    # Resolved against the project root, because the alignment path is later reported
    # relative to it. A caller passing "runs/<id>" is as legitimate as an absolute path.
    derived = Path(derived)
    if not derived.is_absolute():
        derived = (ROOT / derived).resolve()
    index_path = derived / "website_indices" / "protein_coordinate_model.json"
    model_index = json.loads(index_path.read_text(encoding="utf-8"))

    alignment = write_reference_alignment(
        model_index,
        derived / "results" / "generic_gene_analysis" / "msa" / "primaries_msa.aln.faa")

    reference_models = cm.primary_reference_models(model_index)
    comparative_index = dict(model_index)
    comparative_index["models"] = reference_models
    comparative_index["n_models"] = len(reference_models)
    try:
        alignment_reference = str(alignment.relative_to(ROOT))
        projection_root = ROOT
    except ValueError:
        alignment_reference = "dataset:" + alignment.relative_to(derived).as_posix()
        projection_root = derived
    comparative_index["msa_coordinate_map"] = {
        **(model_index.get("msa_coordinate_map") or {}),
        "alignment_file": alignment_reference,
        "keyed_by": "species_id",
        "reason": ("Full-length FGFR2 alignment restricted to each species' primary "
                   "reference protein, so the comparative layer has one sequence per "
                   "species. Columns are the original alignment's columns."),
    }
    comparative_index["comparative_reduction"] = {
        "rule": (model_index.get("comparative_primary") or {}).get("rule", ""),
        "n_models_in_dataset": model_index.get("n_models"),
        "n_models_in_comparative_view": len(reference_models),
        "note": ("Both isoform models stay in the coordinate index and in each "
                 "species' Gallery scope. Only the comparative row is one protein, "
                 "because a comparative row can only be one."),
    }

    # The generic cross-species boundary analysis is already on the coordinate index,
    # built there by the shared dashboard over this same reduced panel. Reusing it
    # rather than rebuilding it here is what guarantees the Figure Gallery, the
    # comparative dataset and the interactive Boundary Explorer quote one set of
    # comparable-group ids, signed distances and classes.
    dataset = build_comparative_dataset(
        derived, coordinate_index=comparative_index,
        project_root=projection_root)
    dataset["comparative_reduction"] = comparative_index["comparative_reduction"]
    out = derived / "website_indices" / "comparative_dataset.json"
    out.write_text(json.dumps(dataset, indent=2, ensure_ascii=False) + "\n",
                   encoding="utf-8")

    # The reduced index is written out as well, because it is what the comparative
    # renderer must be given. Handing it the full 58-model index would put two rows
    # per species in every figure that draws one row per model — the same species
    # twice, once as IIIb and once as IIIc — which reads as thirty extra species.
    (derived / "website_indices" / "comparative_model_index.json").write_text(
        json.dumps(comparative_index, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8")
    return dataset


def main(argv: Optional[List[str]] = None) -> int:
    import argparse
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--derived", default=str(DERIVED))
    args = ap.parse_args(argv)
    dataset = build(Path(args.derived))
    print(json.dumps({
        "n_species": dataset.get("n_species"),
        "msa_columns": (dataset.get("msa") or {}).get("n_columns"),
        "msa_available": (dataset.get("msa") or {}).get("available"),
        "reduction": dataset["comparative_reduction"]["n_models_in_comparative_view"],
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
