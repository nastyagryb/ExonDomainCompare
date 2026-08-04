#!/usr/bin/env python3


from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Dict, List, Tuple

sys.path.insert(0, str(Path(__file__).resolve().parent))
from exondomaincompare.scientific import fgfr2_msa_common as M  # noqa: E402
import reconcile_fgfr2_exon_type_labels as RC  # noqa: E402


MANIFEST_COLS = [
    "msa_input_id", "species", "display_species_name", "isoform", "upstream_label",
    "legacy_label", "previous_pipeline_label", "final_isoform_label", "validated_exon_type",
    "label_consistency_status", "label_reconciliation_action", "final_label_source",
    "final_claim_status_after_rescue", "rescue_decision", "protein_id",
    "transcript_id", "sequence_type", "source_file", "native_start_aa", "native_end_aa",
    "sequence_length", "sequence_hash", "recommended_use", "final_display_class",
    "extraction_status", "extraction_warning",
]


def curated_cassette_span(seq: str, iso: str) -> Tuple[int, int]:
    ref = RC.CURATED_IIIB_REF if iso == "IIIb" else RC.CURATED_IIIC_REF
    if RC._ALN is None or not seq:
        return (None, None)
    try:
        a = RC._ALN.align(ref, seq)[0]
    except Exception:  # noqa: BLE001
        return (None, None)
    idx = a.indices
    ps = [p for r, p in zip(idx[0], idx[1]) if r >= 0 and p >= 0]
    return (min(ps) + 1, max(ps) + 1) if ps else (None, None)


def load_overrides(base: Path, dirs):
    ov, seqs = {}, {}
    p = dirs["maps"] / "fgfr2_rescue_overrides.tsv"
    for r in M.read_tsv(p):
        ov[((r["species"] or "").lower(), r.get("final_isoform_label", ""))] = r
    fa = dirs["inputs"] / "fgfr2_rescued_candidate_proteins.faa"
    for hid, seq in M.read_fasta(fa):
        parts = hid.split("|")
        if len(parts) >= 2:
            seqs[(parts[0].lower(), parts[1])] = M.ungapped(M.clean_alignment_seq(seq))
    return ov, seqs
VALIDATION_COLS = [
    "msa_input_id", "sequence_type", "check", "status", "detail",
]
# plausible cassette peptide length window (FGFR2 IIIb/IIIc alternative exon ~ 50 aa)
CASSETTE_MIN, CASSETTE_MAX = 20, 130


def build_protein_lookup(faa: Path) -> Tuple[Dict[Tuple[str, str, str], str], Dict[str, str]]:
    by_key: Dict[Tuple[str, str, str], str] = {}
    by_pid: Dict[str, str] = {}
    for hid, seq in M.read_fasta(faa):
        meta = {}
        for tok in hid.split("|"):
            if "=" in tok:
                k, v = tok.split("=", 1)
                meta[k] = v
        sp = (meta.get("species") or "").lower()
        iso = meta.get("isoform") or ""
        pid = meta.get("protein") or ""
        s = M.ungapped(M.clean_alignment_seq(seq))
        if pid:
            by_pid.setdefault(pid, s)
            by_key.setdefault((sp, iso, pid), s)
    return by_key, by_pid


def main() -> int:
    ap = argparse.ArgumentParser(description="Prepare FGFR2 MSA inputs.")
    ap.add_argument("--base", type=Path, required=True)
    args = ap.parse_args()
    base = args.base.resolve()
    dirs = M.ensure_module_dirs(base)
    inp = dirs["inputs"]

    coord = M.require(base, "fgfr2_current_stage_IIIb_IIIc_coordinate_audit.tsv")
    cmap_p = M.require(base, "fgfr2_cassette_cds_block_map.tsv")
    master_p = M.require(base, "species_qc_master.tsv", "11_pre_interpro_master")
    faa = M.require(base, "selected_fgfr2_proteins.faa")

    coord_rows = M.read_tsv(coord)
    cmap = {(r["species"].lower(), r["isoform"]): r for r in M.read_tsv(cmap_p)}
    master = {r["species"].lower(): r for r in M.read_tsv(master_p)}
    by_key, by_pid = build_protein_lookup(faa)
    recon = M.load_label_reconciliation(base)
    if not recon:
        print("[WARN] label reconciliation table missing; falling back to upstream labels.",
              file=sys.stderr)
    ov, ov_seq = load_overrides(base, dirs)
    n_override = 0

    full_items: List[Tuple[str, str]] = []
    iiib_items: List[Tuple[str, str]] = []
    iiic_items: List[Tuple[str, str]] = []
    combined_items: List[Tuple[str, str]] = []
    full_man: List[Dict[str, object]] = []
    iiib_man: List[Dict[str, object]] = []
    iiic_man: List[Dict[str, object]] = []
    combined_man: List[Dict[str, object]] = []
    validation: List[Dict[str, object]] = []
    seen_ids: Dict[str, int] = {}

    def vrow(mid, stype, check, status, detail=""):
        validation.append({"msa_input_id": mid, "sequence_type": stype, "check": check,
                           "status": status, "detail": detail})

    for c in coord_rows:
        sp = (c.get("species_canonical") or "").lower()
        up_iso = c.get("inferred_isoform") or ""
        rrow = recon.get((sp, up_iso), {})
        iso = rrow.get("final_isoform_label") or up_iso  # FINAL biological label
        pid = c.get("protein_id") or ""
        tx = c.get("transcript_id_source") or ""
        mr = master.get(sp, {})
        disp = mr.get("display_species_name", sp)
        ruse = mr.get("recommended_use", "")
        fdc = mr.get("final_display_class", "")
        tok = M.recommended_use_token(ruse)
        recon_extra = {
            "upstream_label": up_iso, "legacy_label": up_iso,
            "previous_pipeline_label": up_iso,
            "final_isoform_label": iso,
            "validated_exon_type": rrow.get("validated_exon_type", up_iso),
            "label_consistency_status": rrow.get("label_consistency_status", "no_reconciliation"),
            "label_reconciliation_action": rrow.get("label_reconciliation_action", "keep_upstream_label"),
            "final_label_source": rrow.get("final_label_source", ""),
            "final_claim_status_after_rescue": rrow.get("final_claim_status_after_rescue",
                                                        rrow.get("final_claim_status", "")),
            "rescue_decision": rrow.get("maximal_rescue_decision", rrow.get("rescue_status", "")),
        }
        _plen = M.to_int(c.get("protein_length_aa")) or M.to_int(c.get("native_protein_length_aa"))

        seq = by_key.get((sp, up_iso, pid)) or by_pid.get(pid) or ""
        src_file = faa.name
        # cassette coordinates: corrected coordinate-overlap map, else native audit
        m = cmap.get((sp, up_iso), {})
        cs = M.to_int(m.get("matched_protein_start_aa")) or M.to_int(c.get("native_protein_start_aa"))
        ce = M.to_int(m.get("matched_protein_end_aa")) or M.to_int(c.get("native_protein_end_aa"))

        # ---- Part F: maximal-rescue override (sequence-validated replacement) ----
        ovr = ov.get((sp, iso))
        if ovr:
            rseq = ov_seq.get((sp, iso), "")
            if rseq:
                seq = rseq
                pid = ovr.get("rescued_protein_id") or pid
                tx = ovr.get("rescued_transcript_id") or tx
                src_file = "fgfr2_rescued_candidate_proteins.faa"
                cs, ce = curated_cassette_span(seq, iso)
                recon_extra["validated_exon_type"] = iso
                recon_extra["label_reconciliation_action"] = "rescued_validated_candidate"
                recon_extra["final_label_source"] = ovr.get("final_label_source", "rescue_override")
                recon_extra["rescue_decision"] = ovr.get("rescue_decision", "rescued")
                n_override += 1

        base_id = f"{sp}|{iso}|{pid}|{tok}"

        # ---------- full length ----------
        fid = base_id
        if seq:
            n = seen_ids.get(fid, 0)
            seen_ids[fid] = n + 1
            full_items.append((fid, seq))
            full_man.append({
                "msa_input_id": fid, "species": sp, "display_species_name": disp,
                "isoform": iso, **recon_extra, "protein_id": pid, "transcript_id": tx,
                "sequence_type": "full_length", "source_file": src_file,
                "native_start_aa": 1, "native_end_aa": len(seq), "sequence_length": len(seq),
                "sequence_hash": M.sha256_text(seq), "recommended_use": ruse,
                "final_display_class": fdc, "extraction_status": "extracted_full_length",
                "extraction_warning": "" if n == 0 else "duplicate_input_id_resolved_by_pipeline",
            })
            vrow(fid, "full_length", "non_empty", "ok")
            bad = M.invalid_residues(seq)
            vrow(fid, "full_length", "valid_alphabet", "ok" if not bad else "fail", bad)
            vrow(fid, "full_length", "unique_id", "ok" if n == 0 else "fail",
                 "" if n == 0 else f"id_seen_{n+1}x")
        else:
            full_man.append({
                "msa_input_id": fid, "species": sp, "display_species_name": disp,
                "isoform": iso, **recon_extra, "protein_id": pid, "transcript_id": tx,
                "sequence_type": "full_length", "source_file": src_file,
                "native_start_aa": "", "native_end_aa": "", "sequence_length": 0,
                "sequence_hash": "", "recommended_use": ruse, "final_display_class": fdc,
                "extraction_status": "protein_sequence_unavailable",
                "extraction_warning": "no selected protein sequence for protein_id",
            })
            vrow(fid, "full_length", "non_empty", "fail", "protein_sequence_unavailable")

        # ---------- cassette ----------
        cas_id = base_id
        cas_seq = ""
        warn = ""
        status = "cassette_extracted"
        if not seq:
            status, warn = "cassette_not_extracted", "protein_sequence_unavailable"
        elif cs is None or ce is None or ce < cs:
            status, warn = "cassette_not_extracted", "cassette_coordinates_unavailable"
        else:
            s0 = max(1, cs)
            e0 = min(len(seq), ce)
            cas_seq = seq[s0 - 1:e0]
            if not cas_seq:
                status, warn = "cassette_not_extracted", "empty_cassette_after_clamp"
            else:
                clen = len(cas_seq)
                if clen < CASSETTE_MIN or clen > CASSETTE_MAX:
                    status = "cassette_extracted_length_outlier"
                    warn = f"cassette_length_{clen}_outside_{CASSETTE_MIN}-{CASSETTE_MAX}"
                if s0 != cs or e0 != ce:
                    warn = (warn + ";" if warn else "") + "cassette_coords_clamped_to_protein"

        man_row = {
            "msa_input_id": cas_id, "species": sp, "display_species_name": disp,
            "isoform": iso, **recon_extra, "protein_id": pid, "transcript_id": tx,
            "source_file": src_file, "native_start_aa": cs if cas_seq else "",
            "native_end_aa": ce if cas_seq else "", "sequence_length": len(cas_seq),
            "sequence_hash": M.sha256_text(cas_seq) if cas_seq else "",
            "recommended_use": ruse, "final_display_class": fdc,
            "extraction_status": status, "extraction_warning": warn,
        }
        stype = "IIIb_cassette" if iso == "IIIb" else "IIIc_cassette"
        man_row["sequence_type"] = stype
        if cas_seq:
            if iso == "IIIb":
                iiib_items.append((cas_id, cas_seq))
                iiib_man.append(man_row)
            else:
                iiic_items.append((cas_id, cas_seq))
                iiic_man.append(man_row)
            combined_items.append((cas_id, cas_seq))
            cm_row = dict(man_row)
            cm_row["sequence_type"] = "combined_cassette"
            combined_man.append(cm_row)
            vrow(cas_id, stype, "non_empty", "ok")
            bad = M.invalid_residues(cas_seq)
            vrow(cas_id, stype, "valid_alphabet", "ok" if not bad else "fail", bad)
            plausible = CASSETTE_MIN <= len(cas_seq) <= CASSETTE_MAX
            vrow(cas_id, stype, "plausible_length", "ok" if plausible else "warning",
                 f"len={len(cas_seq)}")
        else:
            # record the failed cassette in the per-isoform manifest too (transparency)
            if iso == "IIIb":
                iiib_man.append(man_row)
            else:
                iiic_man.append(man_row)
            vrow(cas_id, stype, "non_empty", "fail", warn)

    # write fasta + manifests
    M.write_fasta(inp / "fgfr2_full_length_protein_msa_input.faa", full_items)
    M.write_tsv(inp / "fgfr2_full_length_protein_msa_input_manifest.tsv", full_man, MANIFEST_COLS)
    M.write_fasta(inp / "fgfr2_IIIb_cassette_msa_input.faa", iiib_items)
    M.write_tsv(inp / "fgfr2_IIIb_cassette_msa_input_manifest.tsv", iiib_man, MANIFEST_COLS)
    M.write_fasta(inp / "fgfr2_IIIc_cassette_msa_input.faa", iiic_items)
    M.write_tsv(inp / "fgfr2_IIIc_cassette_msa_input_manifest.tsv", iiic_man, MANIFEST_COLS)
    M.write_fasta(inp / "fgfr2_IIIb_IIIc_combined_cassette_msa_input.faa", combined_items)
    M.write_tsv(inp / "fgfr2_IIIb_IIIc_combined_cassette_msa_input_manifest.tsv",
                combined_man, MANIFEST_COLS)
    M.write_tsv(inp / "fgfr2_msa_input_validation.tsv", validation, VALIDATION_COLS)

    n_fail = sum(1 for v in validation if v["status"] == "fail")
    print(f"[OK] MSA inputs: full={len(full_items)} IIIb={len(iiib_items)} "
          f"IIIc={len(iiic_items)} combined={len(combined_items)}")
    print(f"     validation rows={len(validation)} fails={n_fail} rescue_overrides_applied={n_override}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
