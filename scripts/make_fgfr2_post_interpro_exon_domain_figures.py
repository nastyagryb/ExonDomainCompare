#!/usr/bin/env python3
"""
make_fgfr2_post_interpro_exon_domain_figures.py

Post-InterProScan visualization + QC for FGFR2 exon-domain architecture.

Combines, per protein (species x isoform):
  * InterProScan protein domains (real matches only, no invented domains)
  * coding exon / CDS blocks mapped to protein amino-acid coordinates
  * IIIb / IIIc splice-exon slot (cassette) coordinates
  * protein length
  * final isoform labels + transcript / protein IDs (from the FINAL truth table)

Design contract
---------------
* This is a POST-InterPro visualization + QC step. It does NOT modify the
  pre-InterPro truth table, the FASTA files, or the primary/review membership.
* InterProScan domain calls support the *architecture*. They must NOT relabel
  IIIb/IIIc. Final IIIb/IIIc labels always come from the final truth table.
* No coordinates are invented. If exon or cassette coordinates are missing for a
  protein, that layer is marked missing and skipped for that protein.

Transmembrane helix
-------------------
InterProScan did not annotate the FGFR2 transmembrane helix. pyTMHMM predictions
(step 15) are integrated as the authoritative TM layer. FGFR2 is a single-pass
type-I receptor; where pyTMHMM reports an extra N-terminal helix it is treated as a
signal anchor, and the membrane-spanning receptor TM is the helix between the
extracellular Ig region and the kinase.

Outputs (under results/final_30_until_interpro_prepare/15_exon_domain_boundary_post_interpro/):
  tables/interpro_domain_features_normalized.tsv
  tables/pytmhmm_tm_features_normalized.tsv
  tables/exon_domain_architecture_features.tsv
  tables/fgfr2_domain_architecture_qc.tsv
  figures/per_species/{species}_{isoform}_exon_domain_architecture.{svg,pdf,png}
  figures/overview/Figure_10_all_species_FGFR2_exon_domain_architecture_primary.{svg,pdf,png}
  figures/overview/Figure_10A_IIIb_exon_domain_architecture_primary.{svg,pdf,png}
  figures/overview/Figure_10B_IIIc_exon_domain_architecture_primary.{svg,pdf,png}
  figures/overview/Figure_10C_mammals_exon_domain_architecture_primary.{svg,pdf,png}
  figures/overview/Figure_10D_nonmammals_exon_domain_architecture_primary.{svg,pdf,png}
  reports/post_interpro_exon_domain_architecture_summary.md
"""

from __future__ import annotations

import csv
import json
import sys
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

# shared paper style (also sets a writable MPLCONFIGDIR + Agg backend)
sys.path.insert(0, str(Path(__file__).resolve().parent))
from exondomaincompare.presentation import fgfr2_plot_style as st  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import FancyBboxPatch, Rectangle  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

# --------------------------------------------------------------------------- #
# paths
# --------------------------------------------------------------------------- #
REPO = Path(__file__).resolve().parent.parent


def display_path(path) -> str:
    """Repo-relative path for display/logging only; falls back to the raw path when
    BASE is a run-local relative path. Never raises and never affects outputs."""
    p = Path(path)
    try:
        return str(p.resolve().relative_to(REPO.resolve()))
    except Exception:
        return str(p)
import os as _os  # run-folder path override (RESULTS_DIR/BASE); legacy default preserved
BASE = Path(_os.environ.get("FGFR2_RESULTS_DIR") or _os.environ.get("RESULTS_DIR")
            or _os.environ.get("BASE") or (REPO / "results" / "final_30_until_interpro_prepare"))
CLOSURE = BASE / "13_final_pre_interpro_closure"
TRUTH_TABLE = CLOSURE / "final_pre_interpro_truth_table.tsv"
MANIFEST = CLOSURE / "freeze" / "final_pre_interpro_sequence_manifest.tsv"

# candidate coordinate/cassette tables (first existing wins)
COORD_CANDIDATES = [
    CLOSURE / "tables" / "figure3C_exon_to_protein_cassette_coordinate_map.tsv",
    CLOSURE / "figure3C_exon_to_protein_cassette_coordinate_map.tsv",
]

# InterProScan primary output folder (current + future preferred layout)
INTERPRO_DIR_CANDIDATES = [
    BASE / "14_interproscan_primary",
    BASE / "14_interproscan" / "primary" / "output",
]

# pyTMHMM prediction folder(s): the aggregate TSV lives in a sibling step-15 folder;
# per-protein *.summary files are the fallback. Discovery also searches recursively.
PYTMHMM_DIR_CANDIDATES = [
    BASE / "15_pytmhmm_primary",
    BASE / "15_exon_domain_boundary_post_interpro",
]

OUT_DIR = BASE / "15_exon_domain_boundary_post_interpro"
FIG_PER = OUT_DIR / "figures" / "per_species"
FIG_OVERVIEW = OUT_DIR / "figures" / "overview"
TABLES = OUT_DIR / "tables"
# native exon-block reconstruction overrides (coordinate-artifact cases), produced
# by scripts/reconstruct_exon_blocks_post_interpro.py
RECON_OVERRIDES = TABLES / "exon_block_reconstruction_overrides.json"
REPORTS = OUT_DIR / "reports"

# --------------------------------------------------------------------------- #
# colour + label constants (colour-blind-safe, consistent across the project)
# --------------------------------------------------------------------------- #
COL = {
    "signal_peptide": "#0F9B9B",   # teal
    "ig_like_domain": "#0072B2",   # blue
    "IIIb_slot": "#9C6ADE",        # purple / lavender
    "IIIc_slot": "#1B9E77",        # green
    "transmembrane": "#E69F00",    # orange (pyTMHMM receptor TM)
    "tm_anchor": "#F6C36B",        # lighter orange (N-terminal signal anchor, hatched)
    "kinase_domain": "#E6B800",    # gold / yellow
    "coding_exon": "#D9DCE1",      # light grey block
    "coding_exon_edge": "#AEB4BD",
    "fgfr_family": "#B0B6BE",      # faint receptor-family underlay
    "other_domain": "#C7CDD4",
    "warn": "#D55E00",             # amber/red outline only
    "backbone": "#8A9099",
    "ink": st.INK,
    "muted": st.MUTED,
}


# member-database preference for choosing a representative interval per domain cluster
DB_PRIORITY = {
    "CDD": 1, "Pfam": 2, "ProSiteProfiles": 3, "SMART": 4, "PRINTS": 5,
    "Gene3D": 6, "SUPERFAMILY": 7, "FunFam": 8, "PANTHER": 9, "PIRSF": 10,
    "NCBIfam": 11, "SFLD": 12, "Hamap": 12, "PIRSR": 12, "ProSitePatterns": 13,
    "MobiDBLite": 40, "Coils": 40, "AntiFam": 40,
}

# databases whose matches are family-level fingerprints / low-complexity and are
# excluded from the *drawn* representative track (still kept in the normalized table)

TAXON_GROUP_ORDER = [
    "Primates", "Other mammals", "Birds", "Reptiles", "Amphibians", "Teleost fish",
]
MAMMAL_GROUPS = {"Primates", "Other mammals"}

# cassette-boundary precision tokens that count as a genuine low-confidence flag
# ("exact" / "codon_split" are normal high-confidence phase outcomes)
LOW_PRECISION_HINTS = ("unknown", "low", "approx", "ambiguous", "uncertain")

# an N-terminal pyTMHMM helix at/under this aa position is treated as a signal
# anchor (type-I receptor signal peptide region), not the membrane-spanning TM
N_TERMINAL_ANCHOR_MAX_AA = 60


# --------------------------------------------------------------------------- #
# small IO helpers
# --------------------------------------------------------------------------- #
def read_tsv(path: Path) -> List[Dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def to_int(v, default: Optional[int] = None) -> Optional[int]:
    try:
        return int(float(str(v).strip()))
    except (ValueError, TypeError):
        return default


def truthy(v: str) -> bool:
    return str(v).strip().lower() in ("true", "1", "yes", "t")


# --------------------------------------------------------------------------- #
# data records
# --------------------------------------------------------------------------- #
@dataclass
class DomainMatch:
    db: str
    sig_acc: str
    sig_desc: str
    interpro_acc: str
    interpro_desc: str
    start: int
    end: int
    dclass: str  # simplified class


@dataclass
class ExonBlock:
    exon_id: str
    label: str
    start: int
    end: int
    is_cassette: bool
    number: int = 0   # 1-based positional exon number (N->C)


@dataclass
class Protein:
    species: str
    isoform: str            # IIIb / IIIc (pipeline key)
    display_species: str
    taxon_group: str
    final_isoform_label: str
    transcript_id: str
    protein_id: str
    protein_length: int
    claim_status: str
    interpro_seq_len: Optional[int] = None
    matches: List[DomainMatch] = field(default_factory=list)
    exons: List[ExonBlock] = field(default_factory=list)
    cassette_slot: Optional[Tuple[int, int]] = None   # (start,end) block coords used for drawing
    cassette_ref: Optional[Tuple[int, int]] = None    # reference cassette_start/end_aa (may carry native offset)
    boundary_precision: Tuple[str, str] = ("", "")
    # exon-block reconstruction (coordinate-artifact cases only)
    exon_block_source: str = "figure3C"
    exon_display_status: str = ""      # native_exon_blocks_reconstructed / cassette_only_high_confidence / ...
    recon_note: str = ""
    # pyTMHMM transmembrane predictions
    tm_segments: List[dict] = field(default_factory=list)   # all {start,end,role}
    receptor_tm: Optional[dict] = None                      # membrane-spanning TM {start,end}
    tm_anchors: List[dict] = field(default_factory=list)    # N-terminal signal-anchor helices
    tm_source: str = ""                                     # source file (relative)
    # derived / representative
    draw_domains: List[dict] = field(default_factory=list)   # {class,label,start,end}
    ig_regions: List[dict] = field(default_factory=list)
    kinase_region: Optional[dict] = None
    interpro_tm: Optional[dict] = None                      # InterProScan TM if any {start,end}
    qc: dict = field(default_factory=dict)


# --------------------------------------------------------------------------- #
# domain classification
# --------------------------------------------------------------------------- #
def classify_domain(db: str, sig_acc: str, sig_desc: str, interpro_desc: str) -> str:
    text = f"{sig_desc} {interpro_desc}".lower()
    acc = (sig_acc or "").lower()

    if "signal peptide" in text or ("signal" in text and "signalp" in db.lower()):
        return "signal_peptide"
    if "transmembrane" in text or "tm-helix" in text or "tmhmm" in db.lower():
        return "transmembrane"

    # family-level fingerprints / whole-receptor models (PANTHER, PRINTS) -> family band
    if db in ("PANTHER", "PRINTS") and (
        "receptor" in text or "fibroblast growth factor" in text or "fgfr" in text
    ):
        return "fgfr_family"

    if ("kinase" in text or "tyrosine" in text or "ptkc" in acc
            or "pk_tyr" in acc or "pkinase" in acc):
        return "kinase_domain"

    if ("immunoglob" in text or "ig-like" in text or "i-set" in text or "igi_" in acc
            or "igc2" in text or "ig subtype" in text or "v-set" in text
            or text.strip().startswith("ig ") or text.strip() == "ig"):
        return "ig_like_domain"

    if "fibroblast growth factor receptor" in text or "fgfr" in text:
        return "fgfr_family"

    return "other_domain"


def fgfr_ig_number(sig_desc: str) -> Optional[int]:
    """Canonical FGFR Ig number from CDD FGFR-specific signatures (IgI_1/2/3_FGFR)."""
    s = (sig_desc or "").lower()
    for n in (1, 2, 3):
        if f"igi_{n}_fgfr" in s:
            return n
    return None


# --------------------------------------------------------------------------- #
# representative-domain building (collapse overlapping same-class matches)
# --------------------------------------------------------------------------- #
def cluster_by_overlap(matches: List[DomainMatch]) -> List[List[DomainMatch]]:
    clusters: List[List[DomainMatch]] = []
    for m in sorted(matches, key=lambda d: (d.start, d.end)):
        placed = False
        for c in clusters:
            cs = min(x.start for x in c)
            ce = max(x.end for x in c)
            if m.start <= ce and m.end >= cs:   # any overlap
                c.append(m)
                placed = True
                break
        if not placed:
            clusters.append([m])
    return clusters


def representative(cluster: List[DomainMatch]) -> DomainMatch:
    return sorted(cluster, key=lambda d: (DB_PRIORITY.get(d.db, 30),
                                          -(d.end - d.start)))[0]


def select_receptor_tm(p: Protein) -> None:
    """Pick the biologically relevant membrane-spanning TM from pyTMHMM segments.

    FGFR2 is a single-pass type-I receptor. pyTMHMM often reports two helices: an
    N-terminal signal-anchor (~aa 20-45) and the real receptor TM between the
    extracellular Ig region and the kinase. We select the receptor TM as the last
    non-N-terminal helix upstream of (or nearest to) the kinase, and record the
    remaining N-terminal helix as a signal anchor.
    """
    segs = sorted(p.tm_segments, key=lambda s: s["start"])
    if not segs:
        p.receptor_tm, p.tm_anchors = None, []
        return
    kin = p.kinase_region["start"] if p.kinase_region else None
    non_anchor = [s for s in segs if s["start"] > N_TERMINAL_ANCHOR_MAX_AA]
    pool = non_anchor if non_anchor else segs
    if kin is not None:
        upstream = [s for s in pool if s["start"] <= kin + 25]
        receptor = max(upstream, key=lambda s: s["start"]) if upstream \
            else min(pool, key=lambda s: abs(s["start"] - kin))
    else:
        receptor = max(pool, key=lambda s: s["start"])
    p.receptor_tm = {"start": receptor["start"], "end": receptor["end"]}
    p.tm_anchors = [s for s in segs
                    if s is not receptor and s["start"] <= N_TERMINAL_ANCHOR_MAX_AA]
    # annotate roles on the stored segments
    for s in p.tm_segments:
        if s is receptor:
            s["role"] = "receptor_tm"
        elif s["start"] <= N_TERMINAL_ANCHOR_MAX_AA:
            s["role"] = "n_terminal_signal_anchor"
        else:
            s["role"] = "additional_tm"


def build_representative_domains(p: Protein) -> None:
    igs = [m for m in p.matches if m.dclass == "ig_like_domain"]
    kins = [m for m in p.matches if m.dclass == "kinase_domain"]
    itm = [m for m in p.matches if m.dclass == "transmembrane"]
    p.interpro_tm = {"start": min(m.start for m in itm),
                     "end": max(m.end for m in itm)} if itm else None

    # Ig regions, left-to-right, positional labels (Ig1/Ig2/Ig3, then Ig-like N)
    ig_regions = []
    for cluster in sorted(cluster_by_overlap(igs), key=lambda c: min(x.start for x in c)):
        rep = representative(cluster)
        canon = next((fgfr_ig_number(x.sig_desc) for x in cluster
                      if fgfr_ig_number(x.sig_desc)), None)
        ig_regions.append({"start": rep.start, "end": rep.end,
                           "canon": canon, "db": rep.db})
    for i, reg in enumerate(ig_regions, start=1):
        reg["order"] = i
        reg["label"] = f"Ig{i}" if i <= 3 else f"Ig-like {i}"
    p.ig_regions = ig_regions

    # kinase: single representative region (union of the dominant cluster)
    kinase_region = None
    if kins:
        clusters = sorted(cluster_by_overlap(kins),
                          key=lambda c: (max(x.end for x in c) - min(x.start for x in c)),
                          reverse=True)
        big = clusters[0]
        kinase_region = {"start": min(x.start for x in big),
                         "end": max(x.end for x in big),
                         "db": representative(big).db}
    p.kinase_region = kinase_region


def finalize_draw(p: Protein) -> None:
    """Assemble the drawn domain list once Ig/kinase (InterProScan) and the receptor
    TM (pyTMHMM) have been resolved."""
    draw = []
    for m in p.matches:
        if m.dclass == "signal_peptide":
            draw.append({"class": "signal_peptide", "label": "SP",
                         "start": m.start, "end": m.end})
    for reg in p.ig_regions:
        draw.append({"class": "ig_like_domain", "label": reg["label"],
                     "start": reg["start"], "end": reg["end"]})
    # pyTMHMM receptor TM (authoritative TM layer)
    if p.receptor_tm:
        draw.append({"class": "transmembrane", "label": "TM",
                     "start": p.receptor_tm["start"], "end": p.receptor_tm["end"]})
    # N-terminal signal-anchor helices (shown lighter, clearly pyTMHMM-derived)
    for a in p.tm_anchors:
        draw.append({"class": "tm_anchor", "label": "sig. anchor",
                     "start": a["start"], "end": a["end"]})
    if p.kinase_region:
        draw.append({"class": "kinase_domain", "label": "kinase",
                     "start": p.kinase_region["start"], "end": p.kinase_region["end"]})
    p.draw_domains = draw


# --------------------------------------------------------------------------- #
# biological QC
# --------------------------------------------------------------------------- #
def run_qc(p: Protein, sp_predictor_present: bool) -> None:
    warnings: List[str] = []
    ig_count = len(p.ig_regions)
    kinase_found = p.kinase_region is not None
    sp_found = any(d["class"] == "signal_peptide" for d in p.draw_domains)

    kinase_start = p.kinase_region["start"] if kinase_found else None
    kinase_end = p.kinase_region["end"] if kinase_found else None
    ig_max_end = max((r["end"] for r in p.ig_regions), default=None)

    # ---- transmembrane (pyTMHMM authoritative; InterProScan compared if present) ----
    tm = p.receptor_tm
    tm_found = tm is not None
    tm_start = tm["start"] if tm_found else None
    tm_end = tm["end"] if tm_found else None
    interpro_tm_found = p.interpro_tm is not None
    if tm_found and interpro_tm_found:
        overlap = min(tm_end, p.interpro_tm["end"]) - max(tm_start, p.interpro_tm["start"])
        tm_agreement = "pytmhmm_and_interpro_agree" if overlap > 0 else "pytmhmm_interpro_disagree"
    elif tm_found:
        tm_agreement = "pytmhmm_only"
    elif interpro_tm_found:
        tm_agreement = "interpro_only"
    else:
        tm_agreement = "none"

    # ---- domain order (Ig < TM < kinase; kinase downstream of TM) ----
    order_ok = True
    order_note = "plausible"
    if not kinase_found:
        order_ok, order_note = False, "kinase domain not detected"
    elif ig_max_end is not None and ig_max_end > kinase_start + 25:
        order_ok, order_note = False, "Ig-like region overlaps/downstream of kinase"
    if tm_found:
        if kinase_found and tm_start > kinase_start + 10:
            order_ok, order_note = False, "TM helix downstream of kinase domain"
        elif ig_max_end is not None and tm_start < ig_max_end - 40:
            order_ok, order_note = False, "TM helix upstream of/within Ig-like region"
    else:
        warnings.append("pyTMHMM transmembrane helix missing")

    # ---- cassette slot position (must be upstream of TM and kinase, not at aa1) ----
    slot_status = "unavailable"
    if p.cassette_slot:
        s0, s1 = p.cassette_slot
        if s0 <= 1:
            slot_status = "at_protein_start"
            warnings.append("cassette slot at protein start (aa<=1)")
        elif kinase_found and not (s1 < kinase_start or s0 > kinase_end):
            slot_status = "overlaps_kinase"
            warnings.append("cassette slot overlaps kinase domain")
        elif tm_found and s0 > tm_end:
            slot_status = "downstream_of_tm"
            warnings.append("cassette slot downstream of transmembrane helix")
        elif kinase_found and s0 > kinase_start:
            slot_status = "downstream_of_kinase"
            warnings.append("cassette slot downstream of kinase domain")
        else:
            slot_status = "upstream_of_tm_ok"
    else:
        warnings.append("cassette slot coordinates missing")

    # ---- exon-domain mapping ----
    if not p.exons:
        # coding exons intentionally hidden but a validated cassette slot is shown
        mapping_status = "cassette_only" if p.cassette_slot else "no_coordinate_data"
    elif p.cassette_slot is None:
        mapping_status = "cassette_missing"
    else:
        mapping_status = "mapped"

    # ---- signal-region support (pyTMHMM N-terminal anchor counts as support) ----
    sp_support = sp_found or bool(p.tm_anchors)
    if not sp_found:
        if p.tm_anchors:
            warnings.append("signal region supported by N-terminal pyTMHMM anchor "
                            "(no dedicated signal-peptide predictor)")
        else:
            warnings.append("signal peptide not annotated (no signal-peptide predictor)")

    # ---- other minor-flag bookkeeping ----
    if ig_count < 3:
        warnings.append(f"only {ig_count} Ig-like region(s) detected (FGFR2 expects ~3)")
    if not any(r["canon"] for r in p.ig_regions):
        warnings.append("Ig1/Ig2/Ig3 numbering positional only (no CDD FGFR-specific Ig numbering)")
    low_precision = False
    for prec in p.boundary_precision:
        pl = (prec or "").lower()
        if pl and any(h in pl for h in LOW_PRECISION_HINTS):
            warnings.append(f"cassette boundary precision: {prec}")
            low_precision = True
            break

    # ---- exon-block reconstruction note (coordinate-artifact cases) ----
    if p.recon_note:
        warnings.append(p.recon_note)

    # ---- final QC status ----
    hard_order_problem = (not order_ok) or slot_status in (
        "overlaps_kinase", "downstream_of_kinase", "downstream_of_tm", "at_protein_start")

    if mapping_status == "no_coordinate_data":
        final = "failed_coordinate_mapping"
    elif not kinase_found and not tm_found and ig_count == 0:
        final = "insufficient_domain_or_tm_data"
    elif hard_order_problem or not kinase_found:
        final = "review_unusual_domain_order"
    elif kinase_found and tm_found and ig_count >= 2 and order_ok:
        minor = (ig_count < 3) or (not sp_support) or low_precision
        final = "architecture_supported_with_minor_flags" if minor else "architecture_supported"
    elif tm_found and order_ok and (kinase_found or ig_count >= 1):
        # pyTMHMM supports the TM architecture while InterProScan domain data is partial
        final = "interpro_partial_tmhmm_supported"
    else:
        final = "insufficient_domain_or_tm_data"

    p.qc = {
        "expected_signal_peptide_found": sp_found,
        "signal_region_supported": sp_support,
        "expected_ig_like_domain_count": ig_count,
        "transmembrane_found": tm_found,
        "pytmhmm_tm_found": tm_found,
        "interpro_tm_found": interpro_tm_found,
        "tm_agreement": tm_agreement,
        "receptor_tm_start_aa": tm_start if tm_found else "",
        "receptor_tm_end_aa": tm_end if tm_found else "",
        "kinase_found": kinase_found,
        "domain_order_status": order_note if not order_ok else "plausible",
        "cassette_slot_position_status": slot_status,
        "exon_domain_mapping_status": mapping_status,
        "exon_block_display_status": p.exon_display_status or "figure3C_native",
        "warnings": "; ".join(warnings) if warnings else "none",
        "final_qc_status": final,
    }


# --------------------------------------------------------------------------- #
# load everything
# --------------------------------------------------------------------------- #
def find_interpro_tsv() -> Tuple[Path, Path]:
    for d in INTERPRO_DIR_CANDIDATES:
        if not d.is_dir():
            continue
        for name in ("input.fasta.tsv", "input.tsv"):
            tsv = d / name
            if tsv.exists():
                return tsv, d
        hits = sorted(d.glob("*.tsv"))
        if hits:
            return hits[0], d
    raise FileNotFoundError(
        "No InterProScan TSV found under: "
        + ", ".join(str(d) for d in INTERPRO_DIR_CANDIDATES))


def find_pytmhmm_source() -> Optional[Path]:
    """Locate a pyTMHMM prediction file. Prefer an aggregate TSV; fall back to the
    per-protein *.summary directory."""
    roots = [d for d in PYTMHMM_DIR_CANDIDATES if d.is_dir()] + [BASE]
    preferred = ["*transmembrane_hits*.tsv", "*tmhmm*summary*.tsv", "*pytmhmm*.tsv",
                 "*tmhmm*.tsv", "*topology*.tsv", "*transmembrane*.tsv",
                 "*tm_predictions*.tsv", "*tm_predictions*.csv"]
    seen: set = set()
    for root in roots:
        for pat in preferred:
            for hit in sorted(root.rglob(pat)):
                if hit.is_file() and hit not in seen:
                    seen.add(hit)
                    return hit
    # fallback: a directory of per-protein *.summary files
    for root in roots:
        summaries = sorted(root.rglob("*.summary"))
        if summaries:
            return summaries[0].parent
    return None


def _parse_topology_line(text: str) -> List[Tuple[int, int]]:
    """Parse a pyTMHMM topology line into TM (start,end) integer pairs (raw coords)."""
    tms: List[Tuple[int, int]] = []
    # supports either a single "start end label" segment or a whole "s e lab; s e lab" line
    for chunk in text.replace("\t", " ").split(";"):
        parts = chunk.split()
        if len(parts) >= 3 and parts[0].lstrip("-").isdigit() and parts[1].lstrip("-").isdigit():
            label = " ".join(parts[2:]).lower()
            if "transmembrane" in label or label in ("m", "tm", "tmhelix"):
                tms.append((int(parts[0]), int(parts[1])))
    return tms


def load_pytmhmm(proteins: Dict[Tuple[str, str], Protein]) -> Optional[Path]:
    """Attach normalized pyTMHMM TM segments (1-based aa) to each protein."""
    src = find_pytmhmm_source()
    if src is None:
        return None

    raw: Dict[Tuple[str, str], List[Tuple[int, int]]] = defaultdict(list)
    files_used = src

    def handle(seq_id: str, tms: List[Tuple[int, int]]) -> None:
        toks = seq_id.lstrip(">").split("|")
        if len(toks) < 2:
            return
        key = (toks[0], toks[1])
        if key in proteins:
            raw[key].extend(tms)

    if src.is_dir():
        for f in sorted(src.glob("*.summary")):
            acc = f.name[:-len(".summary")]
            tms: List[Tuple[int, int]] = []
            for line in f.read_text(encoding="utf-8", errors="replace").splitlines():
                tms.extend(_parse_topology_line(line))
            handle(acc, tms)
    else:
        with src.open(encoding="utf-8", errors="replace") as fh:
            header = fh.readline().rstrip("\n").split("\t")
            has_line_col = "line" in [h.strip().lower() for h in header]
            cols_lower = [h.strip().lower() for h in header]
            for line in fh:
                line = line.rstrip("\n")
                if not line:
                    continue
                cells = line.split("\t")
                seq_id = cells[0]
                if has_line_col and len(cells) >= 2:
                    handle(seq_id, _parse_topology_line(cells[1]))
                else:
                    # generic wide format: look for start/end columns
                    row = dict(zip(cols_lower, cells))
                    s = row.get("tm_start_aa") or row.get("tm_start") or row.get("start")
                    e = row.get("tm_end_aa") or row.get("tm_end") or row.get("end")
                    lab = (row.get("label") or row.get("topology") or "transmembrane").lower()
                    si, ei = to_int(s), to_int(e)
                    if si is not None and ei is not None and "transmembrane" in lab or (
                            si is not None and ei is not None and not any(
                                k in lab for k in ("inside", "outside"))):
                        handle(seq_id, [(si, ei)])

    # detect 0-based coordinates (any segment starting at 0) -> shift to 1-based aa
    all_starts = [s for pairs in raw.values() for (s, _e) in pairs]
    off = 1 if all_starts and min(all_starts) == 0 else 0

    src_rel = display_path(src) if src.is_file() else display_path(src) + "/*.summary"
    for key, pairs in raw.items():
        p = proteins[key]
        p.tm_source = src_rel
        for (s, e) in pairs:
            p.tm_segments.append({"start": s + off, "end": e + off, "role": "tm"})
        p.tm_segments.sort(key=lambda x: x["start"])
    return files_used


def apply_exon_block_overrides(proteins: Dict[Tuple[str, str], Protein]) -> None:
    """Replace template/offset exon blocks for the coordinate-artifact cases with
    native CDS-derived blocks (or hide them) and pin the cassette slot to the
    validated reference coordinate. Driven by RECON_OVERRIDES; no-op if absent."""
    if not RECON_OVERRIDES.exists():
        return
    try:
        overrides = json.loads(RECON_OVERRIDES.read_text(encoding="utf-8"))
    except (json.JSONDecodeError, OSError):
        return
    for token, ov in overrides.items():
        sp, _, iso = token.partition("|")
        p = proteins.get((sp, iso))
        if p is None:
            continue
        blocks = ov.get("exon_blocks", []) or []
        p.exons = [ExonBlock(exon_id=b.get("exon_id", ""),
                             label=b.get("label", ""),
                             start=int(b["start"]), end=int(b["end"]),
                             is_cassette=bool(b.get("is_cassette", False)),
                             number=int(b.get("number", 0)))
                   for b in blocks]
        cass = ov.get("cassette")
        if cass and cass.get("start") is not None and cass.get("end") is not None:
            p.cassette_slot = (int(cass["start"]), int(cass["end"]))
            p.cassette_ref = p.cassette_slot
        p.exon_display_status = ov.get("final_display_status", "")
        p.exon_block_source = ov.get("exon_block_source", "reconstructed")
        p.recon_note = ov.get("recon_note", "")


def clamp_exon_blocks(proteins: Dict[Tuple[str, str], Protein]) -> None:
    """Display-coordinate safety net: no displayed coding-exon block or cassette
    slot may end past the protein length. Blocks that end at most +2 aa past the
    length are clamped (codon-boundary rounding). Blocks starting entirely past
    the length are dropped. Reconstruction / hide decisions taken upstream are
    preserved; only proteins with no prior status get the minor-clamp flag."""
    for p in proteins.values():
        L = p.protein_length
        if not L:
            continue
        clamped = False
        kept: List[ExonBlock] = []
        for b in p.exons:
            if b.start is not None and b.start > L:
                clamped = True
                continue
            if b.end is not None and b.end > L:
                b.end = L
                if b.start is not None and b.start > b.end:
                    b.start = b.end
                clamped = True
            kept.append(b)
        p.exons = kept
        if p.cassette_slot:
            s0, s1 = p.cassette_slot
            if s1 > L:
                p.cassette_slot = (min(s0, L), L)
                clamped = True
        if clamped and not p.exon_display_status:
            p.exon_display_status = "minor_length_clamped"
            if not p.recon_note:
                p.recon_note = ("final coding-exon block end clamped to protein length "
                                "(+1/+2 aa codon-boundary rounding); coordinates otherwise "
                                "unchanged")


def load_proteins() -> Tuple[Dict[Tuple[str, str], Protein], dict]:
    truth = read_tsv(TRUTH_TABLE)
    manifest = read_tsv(MANIFEST)
    primary_keys = {(r["species"], r["isoform"]) for r in manifest
                    if truthy(r.get("included_in_primary_interpro", ""))}

    truth_by_key = {(r["species"], r["isoform"]): r for r in truth}

    interpro_tsv, interpro_dir = find_interpro_tsv()

    proteins: Dict[Tuple[str, str], Protein] = {}
    normalized_rows: List[dict] = []
    dbs_seen: set = set()

    with interpro_tsv.open(encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            cols = line.split("\t")
            if len(cols) < 8:
                continue
            acc = cols[0]
            toks = acc.split("|")
            if len(toks) < 2:
                continue
            species, isoform = toks[0], toks[1]
            key = (species, isoform)
            tr = truth_by_key.get(key)
            if tr is None:
                continue  # accession not in final truth table -> skip (never invent)

            seqlen = to_int(cols[2])
            db = cols[3]
            sig_acc = cols[4]
            sig_desc = cols[5]
            start = to_int(cols[6])
            end = to_int(cols[7])
            interpro_acc = cols[11] if len(cols) > 11 and cols[11] != "-" else ""
            interpro_desc = cols[12] if len(cols) > 12 and cols[12] != "-" else ""
            if start is None or end is None:
                continue
            dbs_seen.add(db)
            dclass = classify_domain(db, sig_acc, sig_desc, interpro_desc)

            p = proteins.get(key)
            if p is None:
                p = Protein(
                    species=species,
                    isoform=isoform,
                    display_species=tr.get("display_species_name", species),
                    taxon_group=tr.get("taxon_group", ""),
                    final_isoform_label=tr.get("final_isoform_label", isoform),
                    transcript_id=tr.get("transcript_id", ""),
                    protein_id=tr.get("protein_id", ""),
                    protein_length=to_int(tr.get("protein_length"), 0),
                    claim_status=tr.get("final_claim_status_after_rescue", ""),
                    interpro_seq_len=seqlen,
                )
                proteins[key] = p

            p.matches.append(DomainMatch(
                db=db, sig_acc=sig_acc, sig_desc=sig_desc,
                interpro_acc=interpro_acc, interpro_desc=interpro_desc,
                start=start, end=end, dclass=dclass))

            normalized_rows.append({
                "species": species,
                "isoform": isoform,
                "final_isoform_label": p.final_isoform_label,
                "transcript_id": p.transcript_id,
                "protein_id": p.protein_id,
                "protein_length": p.protein_length,
                "interpro_accession": interpro_acc,
                "interpro_description": interpro_desc,
                "member_database": db,
                "signature_accession": sig_acc,
                "signature_description": sig_desc,
                "domain_class_simplified": dclass,
                "domain_start_aa": start,
                "domain_end_aa": end,
                "source_file": display_path(interpro_tsv),
            })

    # attach coordinate / cassette blocks
    coord_path = next((c for c in COORD_CANDIDATES if c.exists()), None)
    coord_rows = read_tsv(coord_path) if coord_path else []
    coord_by_key: Dict[Tuple[str, str], List[dict]] = defaultdict(list)
    for r in coord_rows:
        coord_by_key[(r["species"], r["isoform"])].append(r)

    for key, p in proteins.items():
        rows = coord_by_key.get(key, [])
        for r in rows:
            s = to_int(r.get("block_start_aa"))
            e = to_int(r.get("block_end_aa"))
            if s is None or e is None:
                continue
            is_cass = truthy(r.get("is_IIIb_cassette", "")) or truthy(r.get("is_IIIc_cassette", "")) \
                or "cassette" in (r.get("feature_label", "").lower())
            p.exons.append(ExonBlock(exon_id=r.get("exon_or_cds_id", ""),
                                     label=r.get("feature_label", ""),
                                     start=s, end=e, is_cassette=is_cass))
            if p.cassette_slot is None and is_cass:
                p.cassette_slot = (s, e)
                cs = to_int(r.get("cassette_start_aa"))
                ce = to_int(r.get("cassette_end_aa"))
                if cs is not None and ce is not None:
                    p.cassette_ref = (cs, ce)
                p.boundary_precision = (r.get("boundary_left_precision", ""),
                                        r.get("boundary_right_precision", ""))
        p.exons.sort(key=lambda b: b.start)
        for i, b in enumerate(p.exons, start=1):
            b.number = i

    # apply native exon-block reconstruction overrides for coordinate-artifact cases
    apply_exon_block_overrides(proteins)
    # universal display-coordinate sanitation: clamp any exon block / cassette
    # slot that still ends past the protein length (+1/+2 codon-boundary rounding)
    clamp_exon_blocks(proteins)

    meta = {
        "interpro_tsv": interpro_tsv,
        "interpro_dir": interpro_dir,
        "coord_path": coord_path,
        "primary_keys": primary_keys,
        "normalized_rows": normalized_rows,
        "dbs_seen": dbs_seen,
        "n_truth": len(truth),
    }
    return proteins, meta


# --------------------------------------------------------------------------- #
# table writers
# --------------------------------------------------------------------------- #
def write_tsv(path: Path, rows: List[dict], columns: List[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)


def write_normalized(meta: dict) -> Path:
    cols = ["species", "isoform", "final_isoform_label", "transcript_id", "protein_id",
            "protein_length", "interpro_accession", "interpro_description",
            "member_database", "signature_accession", "signature_description",
            "domain_class_simplified", "domain_start_aa", "domain_end_aa", "source_file"]
    path = TABLES / "interpro_domain_features_normalized.tsv"
    rows = sorted(meta["normalized_rows"],
                  key=lambda r: (r["species"], r["isoform"], r["domain_start_aa"]))
    write_tsv(path, rows, cols)
    return path


def write_pytmhmm_normalized(proteins: Dict[Tuple[str, str], Protein]) -> Path:
    cols = ["species", "isoform", "final_isoform_label", "transcript_id", "protein_id",
            "protein_length", "tm_start_aa", "tm_end_aa", "prediction_source",
            "source_file", "status"]
    rows: List[dict] = []
    for key in sorted(proteins):
        p = proteins[key]
        for s in p.tm_segments:
            rows.append({
                "species": p.species, "isoform": p.isoform,
                "final_isoform_label": p.final_isoform_label,
                "transcript_id": p.transcript_id, "protein_id": p.protein_id,
                "protein_length": p.protein_length,
                "tm_start_aa": s["start"], "tm_end_aa": s["end"],
                "prediction_source": "pyTMHMM", "source_file": p.tm_source,
                "status": s.get("role", "tm"),
            })
        if not p.tm_segments:
            rows.append({
                "species": p.species, "isoform": p.isoform,
                "final_isoform_label": p.final_isoform_label,
                "transcript_id": p.transcript_id, "protein_id": p.protein_id,
                "protein_length": p.protein_length,
                "tm_start_aa": "", "tm_end_aa": "", "prediction_source": "pyTMHMM",
                "source_file": p.tm_source, "status": "no_tm_predicted",
            })
    path = TABLES / "pytmhmm_tm_features_normalized.tsv"
    write_tsv(path, rows, cols)
    return path


def write_features(proteins: Dict[Tuple[str, str], Protein]) -> Path:
    cols = ["species", "isoform", "transcript_id", "protein_id", "protein_length",
            "feature_type", "feature_label", "start_aa", "end_aa", "track", "source", "status"]
    rows: List[dict] = []
    for key in sorted(proteins):
        p = proteins[key]
        base = {"species": p.species, "isoform": p.isoform,
                "transcript_id": p.transcript_id, "protein_id": p.protein_id,
                "protein_length": p.protein_length}

        # representative domain features (Ig / kinase / SP from InterProScan; TM from pyTMHMM)
        for d in p.draw_domains:
            if d["class"] in ("transmembrane", "tm_anchor"):
                ftype = "transmembrane_pytmhmm"
                src = "pyTMHMM"
                status = "receptor_tm" if d["class"] == "transmembrane" \
                    else "n_terminal_signal_anchor"
            else:
                ftype = d["class"]
                src = "interproscan"
                status = "representative"
            rows.append({**base, "feature_type": ftype, "feature_label": d["label"],
                         "start_aa": d["start"], "end_aa": d["end"],
                         "track": "domain", "source": src, "status": status})
        # one aggregated fgfr-family band (kept but not a core drawn block)
        fam = [m for m in p.matches if m.dclass == "fgfr_family"]
        if fam:
            rows.append({**base, "feature_type": "other_domain",
                         "feature_label": "FGFR family fingerprint",
                         "start_aa": min(m.start for m in fam),
                         "end_aa": max(m.end for m in fam),
                         "track": "domain", "source": "interproscan", "status": "family_level"})

        # coding exon / CDS blocks + cassette slots (numbered N->C)
        exon_src = p.exon_block_source or "figure3C"
        exon_status = "native_reconstructed" if p.exon_display_status.startswith("native") \
            else "coordinate_mapped"
        emitted_cassette = False
        for b in p.exons:
            if b.is_cassette:
                ftype = "IIIb_slot" if p.isoform.upper() == "IIIB" else "IIIc_slot"
                rows.append({**base, "feature_type": ftype,
                             "feature_label": f"exon {b.number} ({b.label or ftype})",
                             "start_aa": b.start, "end_aa": b.end,
                             "track": "cassette", "source": "figure3C",
                             "status": "coordinate_mapped"})
                emitted_cassette = True
            else:
                rows.append({**base, "feature_type": "coding_exon",
                             "feature_label": f"exon {b.number} ({b.label or b.exon_id})",
                             "start_aa": b.start, "end_aa": b.end,
                             "track": "exon", "source": exon_src,
                             "status": exon_status})
        # reconstructed / cassette-only cases carry the cassette in cassette_slot only
        if not emitted_cassette and p.cassette_slot:
            ftype = "IIIb_slot" if p.isoform.upper() == "IIIB" else "IIIc_slot"
            s0, s1 = p.cassette_slot
            rows.append({**base, "feature_type": ftype,
                         "feature_label": f"{p.final_isoform_label} cassette (validated reference)",
                         "start_aa": s0, "end_aa": s1,
                         "track": "cassette", "source": "validated_reference_coordinate",
                         "status": p.exon_display_status or "reference_coordinate"})

        # warnings
        wtext = p.qc.get("warnings", "none")
        if wtext and wtext != "none":
            for w in wtext.split("; "):
                rows.append({**base, "feature_type": "warning", "feature_label": w,
                             "start_aa": "", "end_aa": "", "track": "warning",
                             "source": "qc", "status": p.qc.get("final_qc_status", "")})
    path = OUT_DIR / "tables" / "exon_domain_architecture_features.tsv"
    write_tsv(path, rows, cols)
    return path


def write_qc(proteins: Dict[Tuple[str, str], Protein]) -> Path:
    cols = ["species", "isoform", "protein_id", "expected_signal_peptide_found",
            "signal_region_supported", "expected_ig_like_domain_count",
            "transmembrane_found", "pytmhmm_tm_found", "interpro_tm_found",
            "tm_agreement", "receptor_tm_start_aa", "receptor_tm_end_aa",
            "kinase_found", "domain_order_status", "cassette_slot_position_status",
            "exon_domain_mapping_status", "exon_block_display_status",
            "warnings", "final_qc_status"]
    rows = []
    for key in sorted(proteins):
        p = proteins[key]
        rows.append({"species": p.species, "isoform": p.isoform, "protein_id": p.protein_id,
                     **p.qc})
    path = OUT_DIR / "tables" / "fgfr2_domain_architecture_qc.tsv"
    write_tsv(path, rows, cols)
    return path


# --------------------------------------------------------------------------- #
# per-species architecture plot
# --------------------------------------------------------------------------- #
def _round_box(ax, x0, x1, yc, h, facecolor, *, edgecolor, lw=0.8, alpha=1.0, zorder=3):
    w = max(x1 - x0, 0.1)
    pad = min(h * 0.14, w * 0.2)
    ax.add_patch(FancyBboxPatch(
        (x0, yc - h / 2), w, h,
        boxstyle=f"round,pad=0,rounding_size={pad}",
        facecolor=facecolor, edgecolor=edgecolor, linewidth=lw,
        alpha=alpha, zorder=zorder, mutation_aspect=1.0))


def slot_color(isoform: str) -> str:
    return COL["IIIb_slot"] if isoform.upper() == "IIIB" else COL["IIIc_slot"]


def plot_per_species(p: Protein) -> None:
    L = max(p.protein_length, max((b.end for b in p.exons), default=0),
            max((d["end"] for d in p.draw_domains), default=0), 1)
    fig_w = max(7.0, min(15.0, 2.4 + L / 85.0))
    fig, ax = plt.subplots(figsize=(fig_w, 3.15))
    ax.set_xlim(-L * 0.04, L * 1.05)
    ax.set_ylim(0.1, 3.05)
    ax.axis("off")

    dom_y, dom_h = 2.15, 0.52
    exon_y, exon_h = 1.15, 0.34

    # backbones
    for yc in (dom_y, exon_y):
        ax.plot([1, L], [yc, yc], color=COL["backbone"], lw=1.1, zorder=1,
                solid_capstyle="round")

    # track labels
    ax.text(-L * 0.035, dom_y, "domains", fontsize=st.FONT["small"], va="center",
            ha="right", color=COL["muted"], rotation=90)
    ax.text(-L * 0.035, exon_y, "CDS exons", fontsize=st.FONT["small"], va="center",
            ha="right", color=COL["muted"], rotation=90)

    # coding-exon blocks hidden (RefSeq protein without local native CDS coords)
    if not any(not b.is_cassette for b in p.exons):
        ax.text(L * 0.5, exon_y, "coding-exon blocks not shown (native CDS unavailable — low confidence)",
                fontsize=st.FONT["small"] - 0.5, ha="center", va="center",
                color=COL["muted"], style="italic", zorder=2)

    # exon blocks (visible rectangles) + exon numbers (N->C)
    for b in p.exons:
        if b.is_cassette:
            continue
        _round_box(ax, b.start, b.end, exon_y, exon_h, COL["coding_exon"],
                   edgecolor=COL["coding_exon_edge"], lw=0.6, zorder=2)
        w = b.end - b.start
        if w >= max(14, L * 0.018):
            ax.text((b.start + b.end) / 2, exon_y, str(b.number),
                    fontsize=st.FONT["small"] - 0.3, ha="center", va="center",
                    color="#4A4F57", zorder=5)
        else:
            # narrow exon: number just above the block, alternating to reduce overlap
            dy = 0.30 if b.number % 2 else 0.44
            ax.text((b.start + b.end) / 2, exon_y + exon_h / 2 + dy, str(b.number),
                    fontsize=st.FONT["small"] - 1.0, ha="center", va="bottom",
                    color="#7A7F87", zorder=5)

    # cassette slot (block coords) — tall highlight across both tracks
    sc = slot_color(p.isoform)
    if p.cassette_slot:
        s0, s1 = p.cassette_slot
        ax.add_patch(Rectangle((s0, exon_y - exon_h / 2 - 0.12), s1 - s0,
                               (dom_y + dom_h / 2 + 0.12) - (exon_y - exon_h / 2 - 0.12),
                               facecolor=sc, alpha=0.14, edgecolor=sc, lw=1.0,
                               zorder=1.5))
        _round_box(ax, s0, s1, exon_y, exon_h, sc, edgecolor=st.INK, lw=0.7, zorder=4)
        cass_num = next((b.number for b in p.exons if b.is_cassette), None)
        if cass_num:
            ax.text((s0 + s1) / 2, exon_y, str(cass_num),
                    fontsize=st.FONT["small"] - 0.3, ha="center", va="center",
                    color="white", fontweight="bold", zorder=5)
        ax.text((s0 + s1) / 2, exon_y - exon_h / 2 - 0.20, p.final_isoform_label,
                fontsize=st.FONT["small"], ha="center", va="top",
                color=sc, fontweight="bold")

    # domains
    for d in p.draw_domains:
        col = COL.get(d["class"], COL["other_domain"])
        _round_box(ax, d["start"], d["end"], dom_y, dom_h, col,
                   edgecolor=st.INK, lw=0.7, zorder=3)
        for xb in (d["start"], d["end"]):
            ax.plot([xb, xb], [dom_y - dom_h / 2 - 0.05, dom_y + dom_h / 2 + 0.05],
                    color=st.INK, lw=0.45, zorder=3)
        w = d["end"] - d["start"]
        lbl = d["label"]
        if w >= L * 0.06:
            ax.text((d["start"] + d["end"]) / 2, dom_y, lbl, fontsize=st.FONT["gene"],
                    ha="center", va="center", color=COL["ink"], fontweight="bold",
                    zorder=5)
        else:
            ax.text((d["start"] + d["end"]) / 2, dom_y + dom_h / 2 + 0.12, lbl,
                    fontsize=st.FONT["small"], ha="center", va="bottom",
                    color=COL["ink"], zorder=5)

    # N / C terminus
    ax.text(1, dom_y + dom_h / 2 + 0.28, "N", fontsize=st.FONT["label"], ha="center",
            va="bottom", color=COL["ink"], fontweight="bold")
    ax.text(L, dom_y + dom_h / 2 + 0.28, "C", fontsize=st.FONT["label"], ha="center",
            va="bottom", color=COL["ink"], fontweight="bold")

    # aa ruler
    ruler_y = 0.55
    ax.plot([1, L], [ruler_y, ruler_y], color=COL["muted"], lw=0.8)
    step = 100 if L <= 900 else 200
    tick = 0
    while tick <= L:
        xt = max(tick, 1)
        ax.plot([xt, xt], [ruler_y - 0.05, ruler_y + 0.05], color=COL["muted"], lw=0.7)
        ax.text(xt, ruler_y - 0.10, str(tick), fontsize=st.FONT["small"], ha="center",
                va="top", color=COL["muted"])
        tick += step
    ax.text(L, ruler_y - 0.10, "aa", fontsize=st.FONT["small"], ha="left", va="top",
            color=COL["muted"])

    # titles
    qc = p.qc.get("final_qc_status", "")
    st.title(ax, f"{p.display_species} — FGFR2 {p.final_isoform_label}",
             subtitle=f"{p.transcript_id} · {p.protein_id} · {p.protein_length} aa · QC: {qc}")

    # legend (only what is drawn / relevant)
    handles = [
        st.legend_patch(COL["ig_like_domain"], "Ig-like domain"),
    ]
    if any(d["class"] == "transmembrane" for d in p.draw_domains):
        handles.append(st.legend_patch(COL["transmembrane"], "TM helix (pyTMHMM)"))
    handles.append(st.legend_patch(COL["kinase_domain"], "kinase domain"))
    handles.append(st.legend_patch(sc, f"{p.final_isoform_label} cassette slot"))
    handles.append(st.legend_patch(COL["coding_exon"], "coding exon (CDS, numbered)"))
    if any(d["class"] == "signal_peptide" for d in p.draw_domains):
        handles.insert(0, st.legend_patch(COL["signal_peptide"], "signal peptide"))
    if any(d["class"] == "tm_anchor" for d in p.draw_domains):
        handles.append(st.legend_patch(COL["tm_anchor"], "N-term signal anchor (pyTMHMM)"))
    if qc.startswith("review") or "unusual" in qc:
        handles.append(Line2D([0], [0], marker="o", color="none",
                              markerfacecolor="none", markeredgecolor=COL["warn"],
                              markeredgewidth=1.4, markersize=8, label="QC flag"))
    ax.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, -0.16),
              ncol=len(handles), fontsize=st.FONT["legend"], frameon=False,
              handlelength=1.1, columnspacing=1.2, handletextpad=0.5)

    stem = f"{p.species}_{p.isoform}_exon_domain_architecture"
    st.savefig(fig, FIG_PER, stem)


# --------------------------------------------------------------------------- #
# overview plot
# --------------------------------------------------------------------------- #
def _order_records(records: List[Protein]) -> List[Protein]:
    def keyf(p: Protein):
        ti = TAXON_GROUP_ORDER.index(p.taxon_group) if p.taxon_group in TAXON_GROUP_ORDER \
            else len(TAXON_GROUP_ORDER)
        return (ti, p.display_species, p.isoform)
    return sorted(records, key=keyf)


def _draw_overview_row(ax, y, p: Protein, L: int, row_h: float) -> None:
    # two stacked sub-bands per row: domains (upper) and coding exons (lower)
    dom_y = y - row_h * 0.22
    exon_y = y + row_h * 0.22
    dh = row_h * 0.30
    eh = row_h * 0.24
    for yc in (dom_y, exon_y):
        ax.plot([1, p.protein_length], [yc, yc], color=COL["backbone"], lw=0.7,
                zorder=1, solid_capstyle="round")

    sc = slot_color(p.isoform)
    # ---- exon sub-band: visible grey blocks + numbers for major exons ----
    for b in p.exons:
        if b.is_cassette:
            continue
        _round_box(ax, b.start, b.end, exon_y, eh, COL["coding_exon"],
                   edgecolor=COL["coding_exon_edge"], lw=0.4, zorder=2)
        if (b.end - b.start) >= max(40, L * 0.05):
            ax.text((b.start + b.end) / 2, exon_y, str(b.number),
                    fontsize=st.FONT["small"] - 1.8, ha="center", va="center",
                    color="#5A5F67", zorder=5)
    # ---- cassette slot: highlight spanning both sub-bands + colored exon block ----
    if p.cassette_slot:
        s0, s1 = p.cassette_slot
        ax.add_patch(Rectangle((s0, dom_y - dh / 2), max(s1 - s0, 1),
                               (exon_y + eh / 2) - (dom_y - dh / 2),
                               facecolor=sc, alpha=0.16, edgecolor="none", zorder=1.4))
        _round_box(ax, s0, s1, exon_y, eh, sc, edgecolor=st.INK, lw=0.4, zorder=4)
    # ---- domain sub-band: Ig / TM / kinase ----
    for reg in p.ig_regions:
        _round_box(ax, reg["start"], reg["end"], dom_y, dh, COL["ig_like_domain"],
                   edgecolor=st.INK, lw=0.4, zorder=3)
    if p.receptor_tm:
        _round_box(ax, p.receptor_tm["start"], p.receptor_tm["end"], dom_y, dh,
                   COL["transmembrane"], edgecolor=st.INK, lw=0.4, zorder=3)
    if p.kinase_region:
        _round_box(ax, p.kinase_region["start"], p.kinase_region["end"], dom_y, dh,
                   COL["kinase_domain"], edgecolor=st.INK, lw=0.4, zorder=3)
    # QC flag marker
    qc = p.qc.get("final_qc_status", "")
    if qc.startswith("review") or "unusual" in qc or qc.startswith("insufficient") \
            or qc.startswith("failed"):
        ax.plot([-L * 0.02], [y], marker="o", markerfacecolor="none",
                markeredgecolor=COL["warn"], markeredgewidth=1.2, markersize=6,
                zorder=6, clip_on=False)
    ax.text(-L * 0.03, y, p.display_species, fontsize=st.FONT["small"] - 0.5,
            va="center", ha="right", color=COL["ink"])


def _overview_legend(fig):
    handles = [
        st.legend_patch(COL["ig_like_domain"], "Ig-like domain"),
        st.legend_patch(COL["transmembrane"], "TM helix (pyTMHMM)"),
        st.legend_patch(COL["kinase_domain"], "kinase domain"),
        st.legend_patch(COL["IIIb_slot"], "IIIb cassette slot"),
        st.legend_patch(COL["IIIc_slot"], "IIIc cassette slot"),
        st.legend_patch(COL["coding_exon"], "coding exon (CDS)"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor="none",
               markeredgecolor=COL["warn"], markeredgewidth=1.4, markersize=8,
               label="QC flag"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=7, fontsize=st.FONT["legend"],
               frameon=False, bbox_to_anchor=(0.5, 0.005))


def plot_overview(records: List[Protein], stem: str, title_text: str,
                  subtitle: str, split_by_isoform: bool = True) -> None:
    if not records:
        return
    L = max(p.protein_length for p in records)

    if split_by_isoform:
        iiib = _order_records([p for p in records if p.isoform.upper() == "IIIB"])
        iiic = _order_records([p for p in records if p.isoform.upper() == "IIIC"])
        panels = [("IIIb", iiib), ("IIIc", iiic)]
        panels = [(lbl, rows) for lbl, rows in panels if rows]
    else:
        panels = [("", _order_records(records))]

    n_rows_max = max(len(rows) for _, rows in panels)
    fig_h = max(4.5, 0.46 * n_rows_max + 1.8)
    fig, axes = plt.subplots(1, len(panels), figsize=(7.2 * len(panels), fig_h),
                             squeeze=False)
    axes = axes[0]

    for ax, (plabel, rows) in zip(axes, panels):
        ax.set_xlim(-L * 0.24, L * 1.04)
        ax.set_ylim(-1, len(rows) + 1)
        ax.invert_yaxis()
        ax.axis("off")
        prev_taxon = None
        for i, p in enumerate(rows):
            y = i
            if p.taxon_group != prev_taxon:
                ax.text(-L * 0.235, y - 0.5, p.taxon_group, fontsize=st.FONT["small"],
                        va="center", ha="left", color=COL["muted"], fontweight="bold")
                prev_taxon = p.taxon_group
            _draw_overview_row(ax, y, p, L, row_h=0.9)
        # aa ruler at top
        ax.plot([1, L], [-0.7, -0.7], color=COL["muted"], lw=0.7)
        step = 200
        t = 0
        while t <= L:
            xt = max(t, 1)
            ax.plot([xt, xt], [-0.8, -0.6], color=COL["muted"], lw=0.6)
            ax.text(xt, -0.95, str(t), fontsize=st.FONT["small"] - 0.5, ha="center",
                    va="bottom", color=COL["muted"])
            t += step
        if plabel:
            ax.set_title(f"{plabel} isoforms", fontsize=st.FONT["subtitle"],
                         fontweight="bold", loc="left", color=COL["ink"])

    fig.suptitle(title_text, fontsize=st.FONT["title"], fontweight="bold", x=0.02,
                 ha="left", y=0.995)
    fig.text(0.02, 0.965, subtitle, fontsize=st.FONT["subtitle"], color=COL["muted"],
             ha="left")
    _overview_legend(fig)
    fig.subplots_adjust(top=0.93, bottom=0.10, left=0.02, right=0.99, wspace=0.18)
    for ext in ("svg", "pdf", "png"):
        fig.savefig(FIG_OVERVIEW / f"{stem}.{ext}", facecolor="white")
    plt.close(fig)


# --------------------------------------------------------------------------- #
# report
# --------------------------------------------------------------------------- #
def write_report(proteins: Dict[Tuple[str, str], Protein], meta: dict,
                 sp_present: bool) -> Path:
    prim = [proteins[k] for k in proteins if k in meta["primary_keys"]]
    prim = sorted(prim, key=lambda p: (p.species, p.isoform))
    n = len(prim)

    def count(pred) -> int:
        return sum(1 for p in prim if pred(p))

    n_domains = count(lambda p: bool(p.matches))
    n_kinase = count(lambda p: p.kinase_region is not None)
    n_tm = count(lambda p: p.receptor_tm is not None)
    n_tm_pred = count(lambda p: bool(p.tm_segments))
    n_ig = count(lambda p: len(p.ig_regions) > 0)
    by_status: Dict[str, int] = defaultdict(int)
    for p in prim:
        by_status[p.qc.get("final_qc_status", "")] += 1

    unusual = [p for p in prim if p.qc.get("final_qc_status", "").startswith("review")
               or "unusual" in p.qc.get("final_qc_status", "")]

    lines: List[str] = []
    A = lines.append
    A("# Post-InterProScan + pyTMHMM FGFR2 exon-domain architecture — summary\n")
    A("This is a **post-InterProScan / pyTMHMM visualization + QC** step. It does **not** modify "
      "the pre-InterPro truth table, FASTA files, or primary/review membership. InterProScan "
      "provides domain annotations (Ig-like, kinase, …), pyTMHMM provides the transmembrane-helix "
      "prediction, and the final IIIb/IIIc labels still come from the final truth table. "
      "InterProScan/pyTMHMM do **not** validate or relabel IIIb/IIIc.\n")

    A("## Input files\n")
    A(f"- InterProScan TSV: `{display_path(meta['interpro_tsv'])}`")
    A(f"- InterProScan folder: `{display_path(meta['interpro_dir'])}`")
    tmsrc = meta.get("pytmhmm_source")
    A(f"- pyTMHMM predictions: `{tmsrc}`" if tmsrc else "- pyTMHMM predictions: NOT FOUND")
    A(f"- Truth table: `{display_path(TRUTH_TABLE)}`")
    A(f"- Sequence manifest: `{display_path(MANIFEST)}`")
    cp = meta["coord_path"]
    A(f"- Coordinate / cassette table: `{display_path(cp) if cp else 'NOT FOUND'}`")
    A(f"- Member databases present: {', '.join(sorted(meta['dbs_seen']))}\n")

    A("## Counts (primary proteins)\n")
    A(f"- Proteins parsed (primary, in truth table): **{n}**")
    A(f"- With a pyTMHMM TM prediction (≥1 helix): **{n_tm_pred}**")
    A(f"- With a resolved receptor TM helix (used for QC/plots): **{n_tm}**")
    A(f"- With ≥1 InterProScan domain match: **{n_domains}**")
    A(f"- With a kinase domain detected: **{n_kinase}**")
    A(f"- With ≥1 Ig-like domain detected: **{n_ig}**\n")

    A("## QC status distribution\n")
    A("| final_qc_status | count |")
    A("| --- | --- |")
    for k in ("architecture_supported", "architecture_supported_with_minor_flags",
              "interpro_partial_tmhmm_supported", "review_unusual_domain_order",
              "insufficient_domain_or_tm_data", "failed_coordinate_mapping"):
        A(f"| {k} | {by_status.get(k, 0)} |")
    A("")
    A(f"- architecture_supported: **{by_status.get('architecture_supported', 0)}**")
    A(f"- architecture_supported_with_minor_flags: "
      f"**{by_status.get('architecture_supported_with_minor_flags', 0)}**")
    A(f"- review / unusual architecture: **{len(unusual)}**\n")

    A("## Proteins with unusual / review domain order\n")
    if unusual:
        A("| species | isoform | protein_id | status | warnings |")
        A("| --- | --- | --- | --- | --- |")
        for p in unusual:
            A(f"| {p.species} | {p.final_isoform_label} | {p.protein_id} | "
              f"{p.qc.get('final_qc_status','')} | {p.qc.get('warnings','')} |")
    else:
        A("_None — no protein showed a biologically impossible domain order._")
    A("")

    A("## Notes\n")
    A("- **InterProScan** annotates protein domains (Ig-like, kinase, …); **pyTMHMM** predicts the "
      "transmembrane helix; the **final truth table** owns the IIIb/IIIc labels. Neither "
      "InterProScan nor pyTMHMM relabels IIIb/IIIc.")
    A("- pyTMHMM is used as the authoritative TM layer (InterProScan did not annotate the TM helix). "
      "Where InterProScan also reports a TM, agreement is recorded in the QC table (`tm_agreement`).")
    A("- FGFR2 is a single-pass type-I receptor. pyTMHMM often reports a second, N-terminal helix; "
      "this is treated as a signal-anchor (N-terminal signal region) and shown separately, not as "
      "the membrane-spanning receptor TM.")
    A("- Missing InterProScan TM is **not** a failure when pyTMHMM predicts the TM helix. A protein "
      "is only flagged when the pyTMHMM TM is missing or the domain order is biologically impossible "
      "(e.g. cassette slot overlapping the kinase, TM downstream of the kinase).")
    A("- Ig1/Ig2/Ig3 are numbered positionally (N→C). Where CDD FGFR-specific Ig signatures "
      "(IgI_1/2/3_FGFR) are absent, the numbering is positional only and flagged as a minor note.")
    A("- Coding exons and IIIb/IIIc cassette slots come from the existing coordinate table (numbered "
      "N→C), not from InterProScan/pyTMHMM. For major-native-offset species the reference "
      "`cassette_start/end_aa` can lie outside the species protein coordinate space; the drawn "
      "cassette slot therefore uses the protein-mapped cassette CDS block coordinates.")

    path = REPORTS / "post_interpro_exon_domain_architecture_summary.md"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")
    return path


# --------------------------------------------------------------------------- #
# main
# --------------------------------------------------------------------------- #
def main() -> int:
    st.apply_rcparams()
    for d in (FIG_PER, FIG_OVERVIEW, TABLES, REPORTS):
        d.mkdir(parents=True, exist_ok=True)

    proteins, meta = load_proteins()

    # attach pyTMHMM transmembrane predictions
    tm_src = load_pytmhmm(proteins)
    meta["pytmhmm_source"] = (display_path(tm_src) if tm_src and tm_src.is_file()
                              else (display_path(tm_src) + "/*.summary"
                                    if tm_src else None))

    # signal-peptide predictor present in the InterProScan run? (informational)
    sp_present = any(m.dclass == "signal_peptide" for p in proteins.values() for m in p.matches)

    for p in proteins.values():
        build_representative_domains(p)   # Ig / kinase / interpro-TM regions
        select_receptor_tm(p)             # pick pyTMHMM receptor TM + anchors
        finalize_draw(p)                  # assemble drawn domain list
        run_qc(p, sp_present)

    norm_path = write_normalized(meta)
    tm_path = write_pytmhmm_normalized(proteins)
    feat_path = write_features(proteins)
    qc_path = write_qc(proteins)

    # remove superseded overview variants from the earlier (InterPro-only) run
    for stale in ("Figure_10_FGFR2_exon_domain_architecture_IIIb_only",
                  "Figure_10_FGFR2_exon_domain_architecture_IIIc_only",
                  "Figure_10_FGFR2_exon_domain_architecture_mammals",
                  "Figure_10_FGFR2_exon_domain_architecture_non_mammals"):
        for ext in ("svg", "pdf", "png"):
            fp = FIG_OVERVIEW / f"{stale}.{ext}"
            if fp.exists():
                fp.unlink()

    # per-species/isoform plots (primary proteins only) — overwrite old plots
    primary = [proteins[k] for k in proteins if k in meta["primary_keys"]]
    for p in sorted(primary, key=lambda x: (x.species, x.isoform)):
        plot_per_species(p)

    # overview: main combined + split IIIb/IIIc + mammal/non-mammal panels (primary only)
    plot_overview(primary,
                  "Figure_10_all_species_FGFR2_exon_domain_architecture_primary",
                  "Figure 10 — FGFR2 exon-domain architecture across species (primary set)",
                  "InterProScan domains + pyTMHMM TM + numbered coding exons + IIIb/IIIc "
                  "cassette slot · final labels from truth table", split_by_isoform=True)
    plot_overview([p for p in primary if p.isoform.upper() == "IIIB"],
                  "Figure_10A_IIIb_exon_domain_architecture_primary",
                  "Figure 10A — FGFR2 IIIb exon-domain architecture (primary)",
                  "IIIb isoforms · InterProScan domains + pyTMHMM TM + coding exons",
                  split_by_isoform=False)
    plot_overview([p for p in primary if p.isoform.upper() == "IIIC"],
                  "Figure_10B_IIIc_exon_domain_architecture_primary",
                  "Figure 10B — FGFR2 IIIc exon-domain architecture (primary)",
                  "IIIc isoforms · InterProScan domains + pyTMHMM TM + coding exons",
                  split_by_isoform=False)
    plot_overview([p for p in primary if p.taxon_group in MAMMAL_GROUPS],
                  "Figure_10C_mammals_exon_domain_architecture_primary",
                  "Figure 10C — FGFR2 exon-domain architecture (mammals, primary)",
                  "Primates + other mammals", split_by_isoform=True)
    plot_overview([p for p in primary if p.taxon_group not in MAMMAL_GROUPS],
                  "Figure_10D_nonmammals_exon_domain_architecture_primary",
                  "Figure 10D — FGFR2 exon-domain architecture (non-mammals, primary)",
                  "Birds, reptiles, amphibians, teleost fish", split_by_isoform=True)

    report_path = write_report(proteins, meta, sp_present)

    # console summary
    n_prim = len(primary)
    print(f"[ok] parsed {len(proteins)} proteins from truth table "
          f"({n_prim} primary), {len(meta['normalized_rows'])} domain rows")
    print(f"[ok] pyTMHMM source: {meta['pytmhmm_source']}")
    print(f"[ok] tables -> {display_path(norm_path)}")
    print(f"                {display_path(tm_path)}")
    print(f"                {display_path(feat_path)}")
    print(f"                {display_path(qc_path)}")
    print(f"[ok] per-species figures -> {display_path(FIG_PER)} ({n_prim}x svg/pdf/png)")
    print(f"[ok] overview figures    -> {display_path(FIG_OVERVIEW)}")
    print(f"[ok] report -> {display_path(report_path)}")
    by = defaultdict(int)
    for p in primary:
        by[p.qc["final_qc_status"]] += 1
    for k, v in sorted(by.items()):
        print(f"     {k}: {v}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
