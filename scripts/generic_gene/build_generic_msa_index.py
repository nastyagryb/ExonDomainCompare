"""Generic MSA handling (gene-agnostic, PART 7).

Rules:
  * multiple protein isoforms present  -> build an isoform-level MSA (MAFFT --auto)
  * multiple species present           -> (future) cross-species primary MSA
  * only one primary sequence          -> msa_status = unavailable_single_sequence

Never fabricates an alignment. Writes ``msa_index.tsv`` with an explicit status,
and (when built) the alignment file under ``msa/``. Uses the plain MAFFT call
extracted conceptually from the FGFR2 MAFFT wrapper (no cassette specifics).
"""
from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any, Dict, List

from .common import GenericContext, load_context, read_fasta, read_tsv, write_tsv
from run_fgfr2_mafft_alignments import run_mafft
from framework.portable_config import load_config

RUNTIME_CONFIG = load_config(repository_root=Path(__file__).resolve().parents[2])

COLUMNS = [
    "msa_id", "msa_kind", "n_sequences", "msa_status", "tool", "input_fasta",
    "alignment_file", "reason",
]


def _mafft_available() -> str:
    return RUNTIME_CONFIG.executable("mafft") or ""


def _run_mafft(in_fasta: Path, out_aln: Path, timeout: int = 600) -> bool:
    mafft = _mafft_available()
    if not mafft:
        return False
    ok, _reason = run_mafft(mafft, in_fasta, out_aln, timeout)
    return ok and out_aln.exists() and out_aln.stat().st_size > 0


def _fasta_species(all_faa: Path) -> Dict[str, str]:
    """Map protein_id -> species_id from headers like '>PROTEIN GENE|species_id'.

    Gene-agnostic: species come from the FASTA header, never hard-coded.
    """
    out: Dict[str, str] = {}
    if not all_faa.is_file():
        return out
    for line in all_faa.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            header = line[1:].strip()
            pid = header.split()[0] if header else ""
            sp = header.split("|")[-1].strip() if "|" in header else ""
            if pid:
                out[pid] = sp
    return out


def _write_species_fasta(all_faa: Path, pids: List[str], dest: Path) -> None:
    """Write a per-species FASTA subset preserving original records for the given pids."""
    want = set(pids)
    keep: List[str] = []
    emit = False
    for line in all_faa.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            pid = line[1:].strip().split()[0] if line[1:].strip() else ""
            emit = pid in want
        if emit:
            keep.append(line)
    dest.write_text("\n".join(keep) + ("\n" if keep else ""), encoding="utf-8")


def build(ctx: GenericContext) -> Dict[str, Any]:
    ctx.assert_not_freeze()
    iso = read_tsv(ctx.core("protein_isoform_index.tsv"))
    all_faa = ctx.core("proteins_all_isoforms.faa")
    seqs = read_fasta(all_faa)
    n_iso = len(iso) or len(seqs)
    # Species order (species_list.txt is authoritative; else distinct from the index).
    slist = ctx.run_dir / "species_list.txt"
    species_order: List[str] = []
    if slist.is_file():
        species_order = [ln.strip() for ln in slist.read_text(encoding="utf-8").splitlines() if ln.strip()]
    pid_species = _fasta_species(all_faa)
    for r in iso:
        sid = r.get("species_id")
        if sid and sid not in species_order:
            species_order.append(sid)
    n_species = len({r.get("species_id") for r in iso if r.get("species_id")}) or 1

    rows: List[Dict[str, Any]] = []
    result: Dict[str, Any] = {"msa_status": "", "alignment_file": ""}
    # Per-species within-species isoform alignments. We NEVER build one alignment
    # mixing isoforms from different species; each species is aligned on its own so
    # the hierarchy stays "species -> isoforms within that species".
    per_species: Dict[str, str] = {}
    msa_dir = ctx.generic_dir / "msa"

    # protein ids grouped by species (order preserved as encountered in the FASTA)
    by_species: Dict[str, List[str]] = {}
    for pid in seqs:
        sid = pid_species.get(pid, "")
        by_species.setdefault(sid, []).append(pid)

    ordered_species = [s for s in species_order if s in by_species] or list(by_species)
    if len(seqs) >= 2 and ordered_species:
        msa_dir.mkdir(parents=True, exist_ok=True)
        for sid in ordered_species:
            pids = by_species.get(sid, [])
            tag = sid or "unknown"
            if len(pids) < 2:
                rows.append({
                    "msa_id": f"{ctx.gene_symbol}_{tag}_isoform_msa",
                    "msa_kind": "isoform_alignment", "n_sequences": len(pids),
                    "msa_status": "unavailable_single_sequence", "tool": "none",
                    "input_fasta": "", "alignment_file": "",
                    "reason": f"Only one protein isoform for {tag}; no within-species alignment is meaningful.",
                })
                continue
            in_fasta = msa_dir / f"isoform_msa_input__{tag}.faa"
            _write_species_fasta(all_faa, pids, in_fasta)
            out_aln = msa_dir / f"isoform_msa__{tag}.aln.faa"
            ok = _run_mafft(in_fasta, out_aln)
            status = "available" if ok else "tool_unavailable"
            if ok:
                per_species[sid] = str(out_aln.relative_to(ctx.run_dir))
            rows.append({
                "msa_id": f"{ctx.gene_symbol}_{tag}_isoform_msa",
                "msa_kind": "isoform_alignment", "n_sequences": len(pids),
                "msa_status": status,
                "tool": "mafft --auto" if ok else (_mafft_available() and "mafft" or "none"),
                "input_fasta": str(in_fasta.relative_to(ctx.run_dir)),
                "alignment_file": str(out_aln.relative_to(ctx.run_dir)) if ok else "",
                "reason": (f"Within-species isoform MSA for {tag} built with MAFFT --auto."
                           if ok else "MAFFT not available or failed; alignment not built."),
            })
        # Back-compat: expose the reference (first) species' alignment under the
        # historical flat name so the stage copy + single-species readers keep working.
        ref_sid = ordered_species[0]
        ref_aln = per_species.get(ref_sid, "")
        if ref_aln:
            (msa_dir / "isoform_msa.aln.faa").write_text(
                (ctx.run_dir / ref_aln).read_text(encoding="utf-8"), encoding="utf-8")
            result["alignment_file"] = "results/generic_gene_analysis/msa/isoform_msa.aln.faa"
        result["msa_status"] = "available" if per_species else "tool_unavailable"
    else:
        rows.append({
            "msa_id": f"{ctx.gene_symbol}_isoform_msa",
            "msa_kind": "isoform_alignment",
            "n_sequences": len(seqs),
            "msa_status": "unavailable_single_sequence",
            "tool": "none",
            "input_fasta": "",
            "alignment_file": "",
            "reason": "Only one primary protein sequence; no alignment is meaningful.",
        })
        result["msa_status"] = "unavailable_single_sequence"
    result["per_species_isoform_alignments"] = per_species

    # Cross-species primary MSA: when >=2 species are analysed, align exactly one
    # primary protein per species (from proteins_primary.faa). This is the headline
    # comparative alignment for a multi-species run (generic, gene-agnostic).
    primary_faa = ctx.core("proteins_primary.faa")
    primary_seqs = read_fasta(primary_faa)
    result["cross_species_alignment_file"] = ""
    result["cross_species_status"] = "not_applicable"
    if n_species >= 2 and len(primary_seqs) >= 2:
        msa_dir = ctx.generic_dir / "msa"
        msa_dir.mkdir(parents=True, exist_ok=True)
        in_fasta = msa_dir / "primaries_msa_input.faa"
        in_fasta.write_text(primary_faa.read_text(encoding="utf-8"), encoding="utf-8")
        out_aln = msa_dir / "primaries_msa.aln.faa"
        ok = _run_mafft(in_fasta, out_aln)
        status = "available" if ok else "tool_unavailable"
        rows.append({
            "msa_id": f"{ctx.gene_symbol}_cross_species_msa",
            "msa_kind": "cross_species_primary_alignment",
            "n_sequences": len(primary_seqs),
            "msa_status": status,
            "tool": "mafft --auto" if ok else (_mafft_available() and "mafft" or "none"),
            "input_fasta": str(in_fasta.relative_to(ctx.run_dir)),
            "alignment_file": str(out_aln.relative_to(ctx.run_dir)) if ok else "",
            "reason": ("Cross-species primary MSA built with MAFFT --auto (one primary per species)."
                       if ok else "MAFFT not available or failed; cross-species alignment not built."),
        })
        result["cross_species_status"] = status
        result["cross_species_alignment_file"] = str(out_aln.relative_to(ctx.run_dir)) if ok else ""

    write_tsv(ctx.out("msa_index.tsv"), rows, COLUMNS)
    result["msa_index.tsv"] = len(rows)
    result["n_isoforms"] = n_iso
    result["n_species"] = n_species
    return result


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--run-id", required=True)
    args = ap.parse_args()
    ctx = load_context(args.run_id)
    res = build(ctx)
    print(f"OK msa_index  status={res['msa_status']}  alignment={res.get('alignment_file') or '—'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
