#!/usr/bin/env python3
"""
reconcile_fgfr2_exon_type_labels.py  (sequence-calibrated IIIb/IIIc label reconciliation)

The upstream IIIb/IIIc labels are NOT blindly trusted: several entries (incl. human & mouse)
carry swapped labels relative to the curated human FGFR2 / UniProt P21802 cassette evidence.
This module determines the VALIDATED exon type from cassette SEQUENCE evidence (local alignment
to curated UniProt-anchored IIIb/IIIc cassette references + isoform-marker residues) and sets the
final biological isoform label from sequence — never from exon order or protein QC alone, and
never from a hard-coded species list (the known species set is used only as a regression test).

Outputs:
  maps/fgfr2_exon_type_label_reconciliation.tsv
  maps/fgfr2_exon_type_label_reconciliation_summary.tsv
  maps/fgfr2_label_reconciliation_warnings.tsv
  inputs/curated_human_FGFR2_IIIb_IIIc_cassette_reference.faa

Human and mouse are positive controls: if their curated-anchor evidence fails, the pipeline stops.
"""

from __future__ import annotations

import argparse
import sys
from collections import Counter
from pathlib import Path
from typing import Dict, List, Tuple

sys.path.insert(0, str(Path(__file__).resolve().parent))
from exondomaincompare.scientific import fgfr2_msa_common as M  # noqa: E402

try:
    from Bio.Align import PairwiseAligner, substitution_matrices
    _ALN = PairwiseAligner()
    _ALN.mode = "local"
    _ALN.substitution_matrix = substitution_matrices.load("BLOSUM62")
    _ALN.open_gap_score = -11
    _ALN.extend_gap_score = -1
except Exception:  # noqa: BLE001
    _ALN = None


# Curated, UniProt P21802-anchored FGFR2 IgIII alternative-exon cassette references.
# IIIb (KGFR, epithelial) carries the SGINSSN...ICKVSNYIG exon; IIIc (BEK, mesenchymal)
# carries the GVNTTDKEI...EVLYIR exon (verified verbatim against human FGFR2 proteins).
CURATED_IIIB_REF = "HSGINSSNAEVLALFNVTEADAGEYICKVSNYIGQANQSAWLTVLP"
CURATED_IIIC_REF = "AAGVNTTDKEIEVLYIRNVTFEDAGEYTCLAGNSIGISFHSAWLTVLP"
IIIB_MARKERS = ("SGINSSNAEV", "ICKVSNYIG")
IIIC_MARKERS = ("GVNTTDKEI", "EVLYIR")

REC_COLS = ["species", "transcript_id", "protein_id", "upstream_label", "legacy_label",
            "previous_pipeline_label", "validated_exon_type", "final_isoform_label",
            "human_IIIb_identity", "human_IIIc_identity", "human_IIIb_coverage",
            "human_IIIc_coverage", "IIIb_marker_present", "IIIc_marker_present",
            "uniprot_anchor_support", "msa_discriminating_support", "sequence_evidence_source",
            "label_consistency_status", "label_reconciliation_action",
            "label_reconciliation_confidence", "label_reconciliation_warning"]
SUM_COLS = ["metric", "value"]
WARN_COLS = ["species", "isoform_context", "warning_type", "detail"]

# regression expectations (TEST ONLY — never used to drive correction)
EXPECT_SWAPPED = {"homo_sapiens", "mus_musculus", "sus_scrofa", "meleagris_gallopavo",
                  "takifugu_rubripes"}
CONTROLS = {"homo_sapiens", "mus_musculus"}


def aln_id_cov(protein: str, ref: str) -> Tuple[float, float]:
    if _ALN is None or not protein or not ref:
        return 0.0, 0.0
    try:
        a = _ALN.align(ref, protein)[0]
    except Exception:  # noqa: BLE001
        return 0.0, 0.0
    idx = a.indices
    cols = sum(1 for r, p in zip(idx[0], idx[1]) if r >= 0 and p >= 0)
    ident = sum(1 for r, p in zip(idx[0], idx[1]) if r >= 0 and p >= 0 and ref[r] == protein[p])
    return (round(ident / cols, 4) if cols else 0.0,
            round(cols / len(ref), 4) if ref else 0.0)


def build_protein_lookup(faa: Path):
    by_key, by_pid = {}, {}
    for hid, seq in M.read_fasta(faa):
        meta = {t.split("=", 1)[0]: t.split("=", 1)[1] for t in hid.split("|") if "=" in t}
        s = M.ungapped(M.clean_alignment_seq(seq))
        pid = meta.get("protein") or ""
        if pid:
            by_pid.setdefault(pid, s)
            by_key.setdefault(((meta.get("species") or "").lower(), meta.get("isoform") or "", pid), s)
    return by_key, by_pid


def protein_scores(protein: str) -> Dict[str, object]:
    iiib_id, iiib_cov = aln_id_cov(protein, CURATED_IIIB_REF)
    iiic_id, iiic_cov = aln_id_cov(protein, CURATED_IIIC_REF)
    return {"iiib_id": iiib_id, "iiic_id": iiic_id, "iiib_cov": iiib_cov, "iiic_cov": iiic_cov,
            "b_marker": any(m in protein for m in IIIB_MARKERS),
            "c_marker": any(m in protein for m in IIIC_MARKERS)}


def _resolve_member(sc: Dict[str, object], assigned: str, pairing_margin: float):
    """Given a protein's scores and its pairing-assigned type, decide validated/status/conf."""
    own_id = sc["iiib_id"] if assigned == "IIIb" else sc["iiic_id"]
    own_cov = sc["iiib_cov"] if assigned == "IIIb" else sc["iiic_cov"]
    marker_ok = sc["b_marker"] if assigned == "IIIb" else sc["c_marker"]
    conflict_marker = sc["c_marker"] if assigned == "IIIb" else sc["b_marker"]
    strong = own_id >= 0.55 and own_cov >= 0.7
    if strong and not (conflict_marker and not marker_ok):
        anchor = "strong"
        validated = assigned
    elif own_id >= 0.45 and marker_ok and not conflict_marker:
        anchor, validated = "weak", assigned
    else:
        anchor, validated = "none", "unresolved"
    msa_support = ("consistent_with_markers" if (validated in ("IIIb", "IIIc") and marker_ok)
                   else "not_corroborated" if validated in ("IIIb", "IIIc") else "n/a")
    conf = ("high" if anchor == "strong" and pairing_margin >= 0.1 else
            "medium" if anchor == "strong" else "low")
    return validated, anchor, msa_support, conf


def reconcile_species(members: List[Tuple[str, str]]):
    """members = [(upstream_label, protein_seq)] for one species (typically IIIb + IIIc).
    Returns per-member dict keyed by upstream_label with validated/final/status/conf/scores."""
    scores = {up: protein_scores(seq) for up, seq in members}
    ups = [up for up, _ in members]
    out: Dict[str, Dict[str, object]] = {}
    if set(ups) == {"IIIb", "IIIc"}:
        b, c = scores["IIIb"], scores["IIIc"]
        keep = b["iiib_id"] + c["iiic_id"]          # IIIb-protein->IIIb, IIIc-protein->IIIc
        swap = b["iiic_id"] + c["iiib_id"]          # cross assignment
        pairing_margin = round(abs(keep - swap), 4)
        assign = {"IIIb": "IIIb", "IIIc": "IIIc"} if keep >= swap else {"IIIb": "IIIc", "IIIc": "IIIb"}
        degenerate = pairing_margin < 0.05
        for up in ("IIIb", "IIIc"):
            sc = scores[up]
            assigned = assign[up]
            validated, anchor, msa_support, conf = _resolve_member(sc, assigned, pairing_margin)
            if degenerate and validated in ("IIIb", "IIIc"):
                validated, anchor, conf = "ambiguous", "ambiguous_pairing", "low"
            out[up] = _verdict(up, sc, validated, anchor, msa_support, conf)
        # enforce a valid {IIIb, IIIc} pair: resolve final-label collisions
        if out["IIIb"]["final"] == out["IIIc"]["final"]:
            strong = [up for up in ("IIIb", "IIIc") if out[up]["anchor"] == "strong"]
            if len(strong) == 1:
                keep_up = strong[0]
                other = "IIIc" if keep_up == "IIIb" else "IIIb"
                comp = "IIIc" if out[keep_up]["final"] == "IIIb" else "IIIb"
                out[other]["final"] = comp
                out[other]["validated"] = "unresolved"
                out[other]["status"] = "unresolved_no_sequence"
                out[other]["action"] = "exclude_from_primary_claim"
                out[other]["conf"] = "low"
            else:
                for up in ("IIIb", "IIIc"):
                    out[up]["final"] = up
                    out[up]["validated"] = "ambiguous"
                    out[up]["status"] = "ambiguous_label_review"
                    out[up]["action"] = "exclude_from_primary_claim"
                    out[up]["conf"] = "low"
    else:
        # non-standard member set: resolve each independently by best curated ref
        for up, _ in members:
            sc = scores[up]
            assigned = "IIIb" if sc["iiib_id"] >= sc["iiic_id"] else "IIIc"
            validated, anchor, msa_support, conf = _resolve_member(sc, assigned,
                                                                   abs(sc["iiib_id"] - sc["iiic_id"]))
            out[up] = _verdict(up, sc, validated, anchor, msa_support, conf)
    return out


def _verdict(upstream, sc, validated, anchor, msa_support, conf):
    if validated in ("IIIb", "IIIc"):
        if validated == upstream:
            status, action, final = "label_consistent", "keep_upstream_label", upstream
        else:
            status, action, final = ("swapped_relative_to_upstream",
                                     "correct_final_label_from_sequence", validated)
    elif validated == "ambiguous":
        status, action, final = "ambiguous_label_review", "manual_review_required", upstream
    else:
        status, action, final = "unresolved_no_sequence", "manual_review_required", upstream
    return {"iiib_id": sc["iiib_id"], "iiic_id": sc["iiic_id"], "iiib_cov": sc["iiib_cov"],
            "iiic_cov": sc["iiic_cov"], "b_marker": sc["b_marker"], "c_marker": sc["c_marker"],
            "validated": validated, "final": final, "anchor": anchor, "msa_support": msa_support,
            "status": status, "action": action, "conf": conf}


def main() -> int:
    ap = argparse.ArgumentParser(description="Sequence-calibrated IIIb/IIIc label reconciliation.")
    ap.add_argument("--base", type=Path, required=True)
    ap.add_argument("--no_control_gate", action="store_true",
                    help="skip the human/mouse control stop-gate (debug only)")
    args = ap.parse_args()
    base = args.base.resolve()
    dirs = M.ensure_module_dirs(base)
    maps = dirs["maps"]

    if _ALN is None:
        print("[FAIL] Biopython PairwiseAligner unavailable; cannot reconcile labels.", file=sys.stderr)
        return 3

    # write curated reference FASTA
    M.write_fasta(dirs["inputs"] / "curated_human_FGFR2_IIIb_IIIc_cassette_reference.faa",
                  [("curated|IIIb|UniProt_P21802|reference", CURATED_IIIB_REF),
                   ("curated|IIIc|UniProt_P21802|reference", CURATED_IIIC_REF)])

    coord = M.read_tsv(M.require(base, "fgfr2_current_stage_IIIb_IIIc_coordinate_audit.tsv"))
    by_key, by_pid = build_protein_lookup(M.require(base, "selected_fgfr2_proteins.faa"))

    # group cassette entries by species
    by_species: Dict[str, List[Dict[str, str]]] = {}
    for c in coord:
        sp = (c.get("species_canonical") or "").lower()
        by_species.setdefault(sp, []).append(c)

    rows: List[Dict[str, object]] = []
    warnings: List[Dict[str, object]] = []
    for sp, entries in by_species.items():
        members, meta = [], {}
        for c in entries:
            up = c.get("inferred_isoform") or ""
            pid = c.get("protein_id") or ""
            tx = c.get("transcript_id_source") or ""
            protein = by_key.get((sp, up, pid)) or by_pid.get(pid) or ""
            meta[up] = {"pid": pid, "tx": tx, "protein": protein}
            if protein:
                members.append((up, protein))
        verdicts = reconcile_species(members) if members else {}
        for c in entries:
            up = c.get("inferred_isoform") or ""
            m = meta.get(up, {})
            pid, tx, protein = m.get("pid", ""), m.get("tx", ""), m.get("protein", "")
            if not protein or up not in verdicts:
                warnings.append({"species": sp, "isoform_context": up,
                                 "warning_type": "no_protein_sequence",
                                 "detail": f"protein {pid} unavailable"})
                rows.append({"species": sp, "transcript_id": tx, "protein_id": pid,
                             "upstream_label": up, "legacy_label": up, "previous_pipeline_label": up,
                             "validated_exon_type": "unresolved", "final_isoform_label": up,
                             "human_IIIb_identity": "", "human_IIIc_identity": "",
                             "human_IIIb_coverage": "", "human_IIIc_coverage": "",
                             "IIIb_marker_present": "", "IIIc_marker_present": "",
                             "uniprot_anchor_support": "none", "msa_discriminating_support": "n/a",
                             "sequence_evidence_source": "none",
                             "label_consistency_status": "unresolved_no_sequence",
                             "label_reconciliation_action": "manual_review_required",
                             "label_reconciliation_confidence": "low",
                             "label_reconciliation_warning": "no protein sequence for evidence"})
                continue
            r = verdicts[up]
            warn = ""
            if r["status"] == "swapped_relative_to_upstream":
                warn = (f"upstream={up} but cassette matches curated {r['validated']} "
                        f"(IIIb_id={r['iiib_id']}, IIIc_id={r['iiic_id']}); corrected from sequence")
                warnings.append({"species": sp, "isoform_context": up,
                                 "warning_type": "label_swapped_relative_to_upstream", "detail": warn})
            elif r["status"] in ("ambiguous_label_review", "unresolved_no_sequence"):
                warn = f"weak/conflicting evidence (IIIb_id={r['iiib_id']}, IIIc_id={r['iiic_id']})"
                warnings.append({"species": sp, "isoform_context": up,
                                 "warning_type": r["status"], "detail": warn})
            rows.append({"species": sp, "transcript_id": tx, "protein_id": pid, "upstream_label": up,
                         "legacy_label": up, "previous_pipeline_label": up,
                         "validated_exon_type": r["validated"], "final_isoform_label": r["final"],
                         "human_IIIb_identity": r["iiib_id"], "human_IIIc_identity": r["iiic_id"],
                         "human_IIIb_coverage": r["iiib_cov"], "human_IIIc_coverage": r["iiic_cov"],
                         "IIIb_marker_present": str(r["b_marker"]).lower(),
                         "IIIc_marker_present": str(r["c_marker"]).lower(),
                         "uniprot_anchor_support": r["anchor"],
                         "msa_discriminating_support": r["msa_support"],
                         "sequence_evidence_source": "uniprot_P21802_curated_cassette_local_alignment",
                         "label_consistency_status": r["status"],
                         "label_reconciliation_action": r["action"],
                         "label_reconciliation_confidence": r["conf"],
                         "label_reconciliation_warning": warn})

    rows.sort(key=lambda d: (d["species"], d["upstream_label"]))
    M.write_tsv(maps / "fgfr2_exon_type_label_reconciliation.tsv", rows, REC_COLS)

    # ---- per-species final-pair sanity (each species must form one IIIb + one IIIc) ----
    fin_by_sp: Dict[str, List[str]] = {}
    for r in rows:
        fin_by_sp.setdefault(r["species"], []).append(str(r["final_isoform_label"]))
    for sp, labs in fin_by_sp.items():
        if sorted(labs) != ["IIIb", "IIIc"]:
            warnings.append({"species": sp, "isoform_context": "species_pair",
                             "warning_type": "final_pair_not_one_IIIb_one_IIIc",
                             "detail": f"final labels={sorted(labs)} (review/exclude from primary claim)"})

    # ---- human/mouse positive-control gate (rule 4) ----
    # The control is only meaningful when a control species is actually part of the
    # run. Distinguish two cases so custom runs are not forced to include human/mouse:
    #   * control species PRESENT but wrong labels -> real failure (stop).
    #   * control species ABSENT from this run     -> not applicable (warn, do not stop).
    # For the full-30 validated panel both controls are present, so the gate is
    # enforced exactly as before.
    _all_species = {r["species"] for r in rows}
    control_ok = True
    n_controls_present = 0
    for ctrl in CONTROLS:
        crows = [r for r in rows if r["species"] == ctrl]
        validated = {r["validated_exon_type"] for r in crows}
        if not crows:
            # control species is simply not in this run panel -> not applicable
            warnings.append({"species": ctrl, "isoform_context": "control",
                             "warning_type": "control_not_applicable",
                             "detail": "control species absent from this run panel "
                                       "(custom run); positive control not evaluated"})
            continue
        n_controls_present += 1
        if validated != {"IIIb", "IIIc"}:
            control_ok = False
            warnings.append({"species": ctrl, "isoform_context": "control",
                             "warning_type": "control_failed",
                             "detail": f"validated types={sorted(validated)} (expected IIIb & IIIc)"})

    # ---- regression test (TEST ONLY) ----
    observed_swapped = {r["species"] for r in rows
                        if r["label_consistency_status"] == "swapped_relative_to_upstream"}
    reg = {
        "expected_swapped_detected": sorted(EXPECT_SWAPPED & observed_swapped),
        "expected_swapped_missing": sorted(EXPECT_SWAPPED - observed_swapped),
        "unexpected_swapped": sorted(observed_swapped - EXPECT_SWAPPED),
    }
    for sp in reg["expected_swapped_missing"]:
        warnings.append({"species": sp, "isoform_context": "regression",
                         "warning_type": "regression_expected_swap_not_detected",
                         "detail": "known swapped species not flagged by sequence evidence"})
    M.write_tsv(maps / "fgfr2_label_reconciliation_warnings.tsv", warnings, WARN_COLS)

    sc = Counter(r["label_consistency_status"] for r in rows)
    summary = [{"metric": "n_cassettes", "value": len(rows)},
               {"metric": "n_label_consistent", "value": sc.get("label_consistent", 0)},
               {"metric": "n_swapped_relative_to_upstream", "value": sc.get("swapped_relative_to_upstream", 0)},
               {"metric": "n_ambiguous_label_review", "value": sc.get("ambiguous_label_review", 0)},
               {"metric": "n_unresolved_no_sequence", "value": sc.get("unresolved_no_sequence", 0)},
               {"metric": "species_swapped", "value": ";".join(sorted(observed_swapped))},
               {"metric": "human_mouse_control_pass", "value": str(control_ok).lower()},
               {"metric": "human_mouse_controls_present", "value": str(n_controls_present)},
               {"metric": "human_mouse_control_status",
                "value": ("not_applicable_no_control_species_in_run" if n_controls_present == 0
                          else ("pass" if control_ok else "fail"))},
               {"metric": "regression_expected_swapped", "value": ";".join(sorted(EXPECT_SWAPPED))},
               {"metric": "regression_detected", "value": ";".join(reg["expected_swapped_detected"])},
               {"metric": "regression_unexpected_swapped", "value": ";".join(reg["unexpected_swapped"])},
               {"metric": "evidence_source", "value": "uniprot_P21802_curated_cassette_local_alignment"}]
    M.write_tsv(maps / "fgfr2_exon_type_label_reconciliation_summary.tsv", summary, SUM_COLS)

    print(f"[OK] label reconciliation: {dict(sc)}")
    print(f"     swapped species: {sorted(observed_swapped)}")
    if n_controls_present == 0:
        print("     human/mouse control: not applicable (no control species in this "
              "run panel — custom run)")
    else:
        print(f"     human/mouse control pass: {control_ok} "
              f"({n_controls_present}/{len(CONTROLS)} control species present)")
    print(f"     regression detected={reg['expected_swapped_detected']} "
          f"missing={reg['expected_swapped_missing']} unexpected={reg['unexpected_swapped']}")
    if not control_ok and not args.no_control_gate:
        print("[FAIL] human/mouse positive control failed; stopping (rule 4).", file=sys.stderr)
        return 4
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
