from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Tuple

from ..common import SharedRunContext, display_species, read_json, read_tsv, to_int


def _resolve_alignment_mode(species_count: int, sequence_count: int) -> Dict[str, str]:
    if species_count >= 2:
        return {"alignment_mode": "cross_species_msa", "tab_label": "MSA"}
    if sequence_count >= 2:
        return {"alignment_mode": "isoform_alignment", "tab_label": "Isoform Alignment"}
    return {"alignment_mode": "unavailable_single_sequence", "tab_label": "Alignment"}


def _read_fasta(path: Path) -> List[Tuple[str, str, str]]:
    out: List[Tuple[str, str, str]] = []
    if not path.is_file():
        return out
    pid = ""
    sp = ""
    buf: List[str] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            if pid:
                out.append((pid, sp, "".join(buf)))
            header = line[1:].strip()
            pid = header.split()[0] if header else ""
            sp = header.split("|")[-1].strip() if "|" in header else ""
            buf = []
        else:
            buf.append(line.strip())
    if pid:
        out.append((pid, sp, "".join(buf)))
    return out


def _species_order(ctx: SharedRunContext) -> List[str]:
    slist = ctx.run_dir / "species_list.txt"
    if slist.is_file():
        return [ln.strip() for ln in slist.read_text(encoding="utf-8").splitlines() if ln.strip()]
    return []


def _build_cross_species_msa(ctx: SharedRunContext, species_count: int,
                             disclaimer: str) -> Dict[str, Any]:
    aln_path = ctx.generic_dir / "msa" / "primaries_msa.aln.faa"
    records = _read_fasta(aln_path)
    order = _species_order(ctx)
    # protein_id -> transcript_id from the core gene-model index
    tx_by_pid: Dict[str, str] = {}
    for r in read_tsv(ctx.core_dir / "gene_model_index.tsv"):
        pid = str(r.get("protein_id") or "")
        if pid and pid not in tx_by_pid:
            tx_by_pid[pid] = str(r.get("transcript_id") or "")

    def sp_rank(sp: str) -> int:
        return order.index(sp) if sp in order else len(order)

    records = sorted(records, key=lambda rec: (sp_rank(rec[1]), rec[0]))
    ref_species = order[0] if order else (records[0][1] if records else "")
    rows: List[Dict[str, Any]] = []
    for pid, sp, seq in records:
        is_primary = sp == ref_species
        rows.append({
            "species": sp,
            "display_species_name": display_species(sp),
            "isoform": "primary",
            "protein_id": pid,
            "transcript_id": tx_by_pid.get(pid, ""),
            "seq": seq,
            "is_primary": is_primary,
            "taxon_group": "",
            "is_human": False,
            "is_review": False,
        })
    n_cols = len(rows[0]["seq"]) if rows else 0
    rel_aln = str(aln_path.relative_to(ctx.run_dir)) if aln_path.is_file() else ""
    alignment = {
        "available": bool(rows),
        "label": "Cross-species primary-protein alignment",
        "file": rel_aln,
        "n_columns": n_cols,
        "rows": rows,
        "reference_sequence": next((r["protein_id"] for r in rows if r["is_primary"]),
                                   rows[0]["protein_id"] if rows else ""),
        "tool": "MAFFT",
        "sequence_count": len(rows),
    }
    return {
        "available": bool(rows),
        "mode": "cross_species_msa",
        "alignment_mode": "cross_species_msa",
        "tab_label": "MSA",
        "species_count": species_count,
        "sequence_count": len(rows),
        "reference_sequence": alignment["reference_sequence"],
        "disclaimer": disclaimer,
        "alignments": {"isoform": alignment},
        "tabs": [{"key": "isoform", "label": "MSA"}],
        "discriminating": [],
        "discriminating_columns_combined": [],
        "conservation": [],
    }


def build_msa_index(ctx: SharedRunContext) -> Dict[str, Any]:
    wi = ctx.website_indices
    iso_idx = read_json(wi / "isoform_alignment_index.json", {})
    primary_sel = read_json(wi / "primary_selection_index.json", {})
    explorer = read_json(wi / "gene_explorer_index.json", {})
    sid = str(primary_sel.get("species_id") or "")
    disp = display_species(sid)
    species_count = to_int(explorer.get("n_species")) or 1
    seq_count = to_int(iso_idx.get("sequence_count")) or len(iso_idx.get("sequences") or [])
    resolved = _resolve_alignment_mode(species_count, seq_count)

    # >=2 species: headline alignment is one primary protein per species.
    if species_count >= 2:
        cross = _build_cross_species_msa(
            ctx, species_count,
            disclaimer="Cross-species alignment of one primary protein per species "
                       "(MAFFT). Column annotations are generic.")
        if cross["available"]:
            return cross
        # fall through to isoform-level below if the cross-species alignment is missing

    if iso_idx.get("status") != "available":
        # One protein model is a property of the annotation, not a failure of this run.
        # Saying so keeps a single-isoform gene from looking broken.
        single = seq_count < 2 and species_count < 2
        return {
            "available": False,
            "mode": resolved["alignment_mode"] if seq_count >= 2 else "unavailable_single_sequence",
            "alignment_mode": resolved["alignment_mode"] if seq_count >= 2 else "unavailable_single_sequence",
            "tab_label": resolved["tab_label"],
            "species_count": species_count,
            "sequence_count": seq_count,
            "alignments": {},
            "tabs": [],
            "disclaimer": iso_idx.get("disclaimer", ""),
            "availability": {
                "state": "not_applicable" if single else "scientifically_unavailable",
                # Named after the analysis, so the page reads "Protein isoform comparison
                # not applicable" rather than a generic unavailability notice.
                "label": "Protein isoform comparison",
                "reason": (
                    "Only one distinct translated protein sequence is available. At least "
                    "two distinct protein sequences are required to detect protein-level "
                    "isoform differences." if single else
                    "No alignable protein set was recovered for this run."),
                "reason_code": "single_unique_protein_sequence" if single else "",
                "badge": "One protein sequence" if single else "",
                "prerequisite_name": "unique_protein_sequence_count",
                "prerequisite_count": seq_count,
            },
        }

    msa_rows: List[Dict[str, Any]] = []
    for seq in iso_idx.get("sequences") or []:
        if not isinstance(seq, dict):
            continue
        is_primary = bool(seq.get("is_primary"))
        pid = str(seq.get("protein_id") or "")
        msa_rows.append({
            "species": sid,
            "display_species_name": disp,
            "isoform": "primary" if is_primary else pid,
            "protein_id": pid,
            "transcript_id": seq.get("transcript_id", ""),
            "seq": seq.get("aligned_sequence", ""),
            "is_primary": is_primary,
            "taxon_group": "",
            "is_human": False,
            "is_review": False,
        })
    msa_rows.sort(key=lambda r: (not r["is_primary"], r["protein_id"]))

    aln_file = iso_idx.get("alignment_file", "")
    n_cols = iso_idx.get("alignment_length") or 0
    if msa_rows and not n_cols:
        n_cols = len(msa_rows[0].get("seq") or "")

    iso_alignment = {
        "available": bool(msa_rows),
        "label": "Within-species isoform alignment",
        "file": aln_file,
        "n_columns": n_cols,
        "rows": msa_rows,
        "reference_sequence": iso_idx.get("reference_sequence", ""),
        "tool": iso_idx.get("tool", "MAFFT"),
        "sequence_count": iso_idx.get("sequence_count", len(msa_rows)),
    }

    return {
        "available": bool(msa_rows),
        "mode": resolved["alignment_mode"],
        "alignment_mode": resolved["alignment_mode"],
        "tab_label": resolved["tab_label"],
        "species_count": species_count,
        "sequence_count": seq_count,
        "reference_sequence": iso_idx.get("reference_sequence", ""),
        "disclaimer": iso_idx.get(
            "disclaimer",
            "Within-species protein isoform alignment — not cross-species conservation.",
        ),
        "alignments": {"isoform": iso_alignment},
        "tabs": [{"key": "isoform", "label": "Within-species isoform alignment"}],
        "discriminating": [],
        "discriminating_columns_combined": [],
        "conservation": [],
    }
