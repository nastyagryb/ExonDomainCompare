from __future__ import annotations

import csv
import re
import textwrap
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.ticker import FuncFormatter, MaxNLocator  # noqa: E402

from exondomaincompare.shared_gene_analysis.strand import (
    MINUS, PLUS, is_reverse, normalize_strand)
from . import fgfr2_plot_style as ps

# Exploratory-candidate overlay colour (amber; never a validated-event colour).
CANDIDATE_COLOR = "#F2B705"
CANDIDATE_EDGE = "#8A5008"
EXON_BLOCK_COLORS = ("#E6EBF2", "#C3CDDB")
AXIS_GREY = "#D2D9E6"
NON_CODING_EXON_COLOR = "#FFFFFF"

# Publication raster resolution for every figure download.
EXPORT_DPI = 300

# The two sentences every exploratory figure must carry verbatim, so a reader can
# never mistake an evidence score for a validation result.
EXPLORATORY_TAG = "Exploratory candidate"
VALIDATION_TAG = "Biological validation: not validated"

# alignment palette — kept identical to the interactive figure source
# (webapp/frontend/src/pages/viewers/alignmentFigure.js) so the static and
# interactive versions of the same figure cannot drift apart.
ALN_RESIDUE = "#8FA8BF"
ALN_RESIDUE_PRIMARY = "#0072B2"
ALN_GAP = "#EDF0F3"
ALN_VARIABLE = "#E69F00"
ALN_IDENTITY = "#117733"


# --------------------------------------------------------------------------- #
# shared style / export helpers (thin, explicit re-exports of the FGFR2 style)
# --------------------------------------------------------------------------- #
def apply_style() -> None:
    ps.apply_rcparams()


def save_figure_all_formats(fig, fig_dir: Path, stem: str, *,
                            dpi: int = EXPORT_DPI) -> None:
    fig_dir = Path(fig_dir)
    fig_dir.mkdir(parents=True, exist_ok=True)
    ps.savefig(fig, fig_dir, stem)
    fig.savefig(fig_dir / f"{stem}.png", dpi=dpi, bbox_inches="tight",
                facecolor="white")


def write_source_table(fig_dir: Path, stem: str, columns: Sequence[str],
                       rows: Sequence[Dict[str, Any]]) -> Optional[Path]:
    if not rows:
        return None
    path = Path(fig_dir) / f"{stem}.tsv"
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(columns), delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({c: row.get(c, "") for c in columns})
    return path


def figure_title(ax, text: str, subtitle: Optional[str] = None,
                 note: Optional[str] = None) -> None:
    note_lines = _wrap_for_axes(ax, note, ps.FONT["small"]) if note else []
    sub_lines = _wrap_for_axes(ax, subtitle, ps.FONT["subtitle"]) if subtitle else []
    ps.title(ax, safe_text(text))
    pad = 12.0 + 11.0 * len(sub_lines) + 10.0 * len(note_lines)
    ax.set_title(safe_text(text), fontsize=ps.FONT["title"], fontweight="bold",
                 loc="left", color=ps.INK, pad=pad)
    offset = 2.0
    for lines, size, style, step in ((note_lines, ps.FONT["small"], "italic", 10.0),
                                     (sub_lines, ps.FONT["subtitle"], "normal", 11.0)):
        for line in reversed(lines):
            ax.annotate(line, xy=(0.0, 1.0), xycoords="axes fraction",
                        xytext=(0.0, offset), textcoords="offset points",
                        fontsize=size, color=ps.MUTED, va="bottom", ha="left",
                        style=style, annotation_clip=False)
            offset += step


def shared_legend(ax, handles, *, ncol: int = 3, bbox=(0.5, -0.05),
                  loc: str = "lower center") -> None:
    ps.compact_legend(ax, handles, ncol=ncol, bbox=bbox, loc=loc)


def legend_patch(color: str, label: str):
    return ps.legend_patch(color, label)


def _legend_marker(marker: str, label: str, color: str = None):
    return plt.Line2D([0], [0], marker=marker, linestyle="none", label=label,
                      markerfacecolor=color or "white", markeredgecolor=ps.INK,
                      markeredgewidth=0.6, markersize=6)


def _int(v, default: int = 0) -> int:
    try:
        return int(float(v))
    except (TypeError, ValueError):
        return default


# --------------------------------------------------------------------------- #
# layout helpers — point-accurate placement so labels never collide
# --------------------------------------------------------------------------- #
def _axes_size_in(ax) -> Tuple[float, float]:
    fig = ax.figure
    pos = ax.get_position()
    return (max(0.2, pos.width * fig.get_figwidth()),
            max(0.2, pos.height * fig.get_figheight()))


def _below_axes(ax, points: float) -> float:
    return -points / (72.0 * _axes_size_in(ax)[1])


def _wrap_for_axes(ax, text: str, fontsize: float, factor: float = 1.15) -> List[str]:
    safe = safe_text(text)
    width_pts = _axes_size_in(ax)[0] * 72.0 * factor
    chars = max(48, int(width_pts / (0.52 * fontsize)))
    return textwrap.wrap(safe, chars) or [safe]


def _footnote(ax, text: str, points: float = 58.0) -> None:
    if not text:
        return
    lines = _wrap_for_axes(ax, text, ps.FONT["small"])
    ax.annotate("\n".join(lines), xy=(0.0, 0.0), xycoords="axes fraction",
                xytext=(0.0, -points), textcoords="offset points",
                fontsize=ps.FONT["small"], color=ps.MUTED, va="top", ha="left",
                linespacing=1.5, annotation_clip=False)


# Some house fonts (macOS Helvetica) lack arrows and stars; PNG/PDF would then
# render an empty box. Degrade to ASCII only when the active font really is
# missing the glyph, so systems with a full font keep the typographic version.
_GLYPH_FALLBACK = {"→": "->", "★": "*", "▶": ">", "◀": "<", "✓": "yes"}
_GLYPH_CACHE: Dict[str, bool] = {}


def _glyph_ok(ch: str) -> bool:
    if ch not in _GLYPH_CACHE:
        try:
            from matplotlib import font_manager
            from matplotlib.ft2font import FT2Font
            family = plt.rcParams.get("font.family") or ["sans-serif"]
            path = font_manager.findfont(font_manager.FontProperties(family=list(family)))
            _GLYPH_CACHE[ch] = FT2Font(path).get_char_index(ord(ch)) != 0
        except Exception:
            _GLYPH_CACHE[ch] = False
    return _GLYPH_CACHE[ch]


def safe_text(text: str) -> str:
    out = str(text or "")
    for ch, replacement in _GLYPH_FALLBACK.items():
        if ch in out and not _glyph_ok(ch):
            out = out.replace(ch, replacement)
    return out


def _clip(text: str, chars: int) -> str:
    out = safe_text(text)
    return out if len(out) <= chars else out[:chars - 3].rstrip() + "..."


def _readable_aa_axis(ax, length: int, label: str, *, lo: int = 0) -> None:
    ticks = MaxNLocator(nbins=10, steps=[1, 2, 5, 10], integer=True).tick_values(lo, length)
    ax.set_xticks([t for t in ticks if lo <= t <= length])
    ax.set_xlabel(label, fontsize=ps.FONT["label"])
    ax.tick_params(axis="x", labelsize=ps.FONT["tick"], length=3, pad=2)
    if ax.spines["bottom"].get_visible():
        ax.spines["bottom"].set_bounds(lo, length)


def _place_block_labels(ax, items: Sequence[Tuple[float, float, str]], *, y: float,
                        height: float, fontsize: float,
                        color: str = None, outside_color: str = None,
                        leader: bool = True, levels: int = 2,
                        above_only: bool = False) -> int:
    lo, hi = ax.get_xlim()
    span = max(1e-9, hi - lo)
    pts_per_unit = (_axes_size_in(ax)[0] * 72.0) / span
    char_w = 0.62 * fontsize
    inside_col = color or ps.INK
    outside_col = outside_color or ps.INK
    signs = (1,) if above_only else (1, -1)
    slots = [0] + [sign * lvl for lvl in range(1, levels + 1) for sign in signs]
    used = {s: -1e18 for s in slots}
    drawn = 0
    ordered = sorted((it for it in items if it[0] is not None and it[1] is not None
                      and str(it[2] or "")),
                     key=lambda it: (it[0] + it[1]) / 2.0)
    for start, end, text in ordered:
        text = str(text)
        need = len(text) * char_w + 3.0
        cx = (start + end) / 2.0
        half = (need / 2.0) / pts_per_unit
        gap = 2.5 / max(pts_per_unit, 1e-9)
        fits_inside = (end - start) * pts_per_unit >= need
        for slot in slots:
            if slot == 0 and not fits_inside:
                continue
            if cx - half < used[slot] + gap:
                continue
            if slot == 0:
                ax.text(cx, y, text, ha="center", va="center", fontsize=fontsize,
                        color=inside_col, zorder=9)
            else:
                sign = 1 if slot > 0 else -1
                ty = y + sign * height * (1.7 + 1.5 * (abs(slot) - 1))
                ax.text(cx, ty, text, ha="center", fontsize=fontsize, color=outside_col,
                        va="bottom" if sign > 0 else "top", zorder=9)
                if leader:
                    ax.plot([cx, cx], [y + sign * height, ty], color=AXIS_GREY,
                            lw=0.5, zorder=2)
            used[slot] = cx + half
            drawn += 1
            break  # placed; otherwise the label is honestly omitted
    return drawn


def _readable_nucleotide_axis(ax, lo: float, hi: float, label: str) -> None:
    ticks = MaxNLocator(nbins=8, steps=[1, 2, 5, 10]).tick_values(lo, hi)
    ax.set_xticks([t for t in ticks if lo <= t <= hi])
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _p: f"{v / 1000.0:,.0f}"))
    ax.set_xlabel(label, fontsize=ps.FONT["label"])
    ax.tick_params(axis="x", labelsize=ps.FONT["tick"], length=3, pad=2)
    if ax.spines["bottom"].get_visible():
        ax.spines["bottom"].set_bounds(lo, hi)


def _pack_rows(intervals: Sequence[Tuple[float, float]]) -> List[int]:
    ends: List[float] = []
    rows: List[int] = []
    for start, end in intervals:
        placed = False
        for r, last in enumerate(ends):
            if start >= last:
                ends[r] = end
                rows.append(r)
                placed = True
                break
        if not placed:
            ends.append(end)
            rows.append(len(ends) - 1)
    return rows


# --------------------------------------------------------------------------- #
# domain naming — coordinate-specific, numbered instances
# --------------------------------------------------------------------------- #
_NAME_EXPANSIONS = (
    (r"\bdom\b", "domain"), (r"\bsf\b", "superfamily"), (r"\bfam\b", "family"),
    (r"\bcat\b", "catalytic"), (r"\brcpt\b", "receptor"), (r"\bsub\b", "subunit"),
    (r"\bbd\b", "binding"), (r"\bTM\b", "transmembrane"),
)
_SHORT_KEYS = (
    ("ig-like", "Ig-like"), ("immunoglobulin", "Ig-like"), ("kinase", "Kinase"),
    ("p53", "p53"), ("tetramer", "Tetramer"), ("transactiv", "TAD"),
    ("fibronectin", "FN3"), ("cadherin", "Cadherin"), ("homeobox", "Homeobox"),
    ("zinc finger", "Zn finger"), ("egf", "EGF"),
)


def _pretty_domain_name(raw: str) -> str:
    name = (raw or "").replace("_", " ").strip()
    for pattern, full in _NAME_EXPANSIONS:
        name = re.sub(pattern, full, name, flags=re.IGNORECASE)
    return name or "domain"


def _short_domain_label(pretty: str) -> str:
    low = pretty.lower()
    for key, short in _SHORT_KEYS:
        if key in low:
            return short
    first = pretty.split()[0] if pretty.split() else pretty
    return first if len(first) <= 12 else first[:11] + "…"


def number_domain_instances(domains: Sequence[Dict[str, Any]]) -> List[Dict[str, Any]]:
    rows = []
    for d in domains:
        start, end = _span(d)
        rows.append({
            "start": start, "end": end,
            "accession": d.get("interpro_accession") or d.get("accession") or "",
            "instance_id": d.get("domain_instance_id") or "",
            "raw_label": (d.get("label") or d.get("interpro_name")
                          or d.get("domain_name") or d.get("domain_id") or ""),
            "source": d,
        })
    rows.sort(key=lambda r: (r["start"] or 0, r["end"] or 0))
    groups: Dict[str, List[Dict[str, Any]]] = {}
    for r in rows:
        groups.setdefault(r["accession"] or r["raw_label"], []).append(r)
    colours: Dict[str, str] = {}
    for key in sorted(groups, key=lambda k: groups[k][0]["start"] or 0):
        colours[key] = DOMAIN_PALETTE[len(colours) % len(DOMAIN_PALETTE)]
    for key, members in groups.items():
        for i, r in enumerate(members, start=1):
            src = r["source"]
            n_instances = _int(src.get("n_instances_of_accession"), len(members)) or len(members)
            number = _int(src.get("instance_number"), i) or i
            pretty = _pretty_domain_name(r["raw_label"])
            r["pretty"] = pretty
            r["instance"] = number if n_instances > 1 else None
            r["name"] = src.get("short_label") or (
                f"{pretty} {number}" if n_instances > 1 else pretty)
            short = _short_domain_label(src.get("short_label") or pretty)
            r["short"] = f"{short} {number}" if n_instances > 1 else short
            r["color"] = colours[key]
            r["legend"] = src.get("full_label") or (
                f"{r['name']} · {r['accession']}" if r["accession"] else r["name"])
            if not r["instance_id"]:
                r["instance_id"] = (f"{r['accession']}:{r['start']}-{r['end']}"
                                    if r["accession"] else f"{r['start']}-{r['end']}")
    # Resolve abbreviation collisions with the full instance label.
    seen: Dict[str, int] = {}
    for r in rows:
        seen[r["short"]] = seen.get(r["short"], 0) + 1
    for r in rows:
        if seen[r["short"]] > 1:
            r["short"] = _clip(r["source"].get("short_label") or r["name"], 24)
    rows.sort(key=lambda r: (r["start"] or 0, r["end"] or 0))
    return rows


# --------------------------------------------------------------------------- #
# gene-agnostic figures
# --------------------------------------------------------------------------- #
def plot_gene_model_overview(fig_dir: Path, stem: str, *, gene_symbol: str,
                             isoforms: Sequence[Dict[str, Any]]) -> bool:
    rows = sorted(isoforms, key=lambda r: _int(r.get("protein_length")), reverse=True)
    if not rows:
        return False
    n = len(rows)
    fig, ax = plt.subplots(figsize=(9.5, max(2.4, 0.5 * n + 1.2)))
    maxlen = max((_int(r.get("protein_length")) for r in rows), default=1) or 1
    for i, r in enumerate(rows):
        y = n - 1 - i
        length = _int(r.get("protein_length"))
        primary = str(r.get("primary_status", "")).lower() == "primary"
        color = ps.INK if primary else ps.PALETTE[i % len(ps.PALETTE)]
        ax.barh(y, length, height=0.6, color=color, alpha=0.92 if primary else 0.7,
                edgecolor=ps.INK, linewidth=0.5, zorder=3)
        label = f"{r.get('protein_id', '')} · {length} aa" + (" · primary" if primary else "")
        ax.text(length + maxlen * 0.01, y, label, va="center", ha="left",
                fontsize=ps.FONT["small"], color=ps.INK if primary else ps.MUTED)
    ax.set_yticks([])
    ax.set_xlim(0, maxlen * 1.35)
    ax.set_xlabel("Protein length (aa)", fontsize=ps.FONT["label"])
    ax.spines[["top", "right", "left"]].set_visible(False)
    figure_title(ax, f"{gene_symbol} — gene model overview",
                 subtitle=f"{n} protein isoform(s); primary selection highlighted")
    shared_legend(ax, [ps.legend_patch(ps.INK, "primary isoform"),
                       ps.legend_patch(ps.PALETTE[0], "alternative isoform")],
                  ncol=2, bbox=(0.5, -0.18))
    save_figure_all_formats(fig, fig_dir, stem)
    return True


# --------------------------------------------------------------------------- #
# transcript structure and the translated protein product
#
# Exon identity between protein models is decided on GENOMIC CDS coordinates.
# Comparing protein offsets instead makes every exon downstream of a deletion look
# altered, because the deletion shifts all later amino-acid positions — which is
# an artefact of the coordinate system, not a difference in the exon used.
# --------------------------------------------------------------------------- #
SHARED_EXON_COLOR = "#C9D3DE"
SHIFTED_EXON_COLOR = "#0072B2"
ALTERNATIVE_EXON_COLOR = "#CC79A7"
INSERTED_REGION_COLOR = "#009E73"
MISSING_MARK_COLOR = "#7C8798"
TERMINUS_COLOR = "#117733"
PRIMARY_EXON_COLOR = "#44546A"

EXON_IDENTITY_COLOR = {
    "primary": PRIMARY_EXON_COLOR,
    "shared": SHARED_EXON_COLOR,
    "shifted": SHIFTED_EXON_COLOR,
    "alternative": ALTERNATIVE_EXON_COLOR,
    "inserted": INSERTED_REGION_COLOR,
    "missing": MISSING_MARK_COLOR,
    "terminus": TERMINUS_COLOR,
}
EXON_IDENTITY_LABEL = {
    "primary": "coding exon of the primary model",
    "shared": "shared genomic exon (identical CDS interval)",
    "shifted": "shifted exon boundary (same exon, different CDS end)",
    "alternative": "alternative exon (replaces a primary exon)",
    "inserted": "inserted protein region (additional exon)",
    "missing": "missing protein region (primary exon not used)",
    "terminus": "alternative terminus",
}
EXON_IDENTITY_TABLE_COLUMNS = (
    "protein_id", "transcript_id", "model_role", "exon_label", "exon_number",
    "genomic_start", "genomic_end", "protein_start_aa", "protein_end_aa",
    "exon_identity", "reference_exon_label", "junction_aa",
)


def _genomic_span(block: Dict[str, Any]) -> Optional[Tuple[int, int]]:
    for lo_key, hi_key in (("genomic_start", "genomic_end"), ("cds_start", "cds_end")):
        lo, hi = block.get(lo_key), block.get(hi_key)
        if lo is not None and hi is not None:
            lo, hi = _int(lo), _int(hi)
            return (min(lo, hi), max(lo, hi))
    return None


def _overlaps(a: Tuple[int, int], b: Tuple[int, int]) -> bool:
    return min(a[1], b[1]) - max(a[0], b[0]) > 0


def derive_strand(blocks: Sequence[Dict[str, Any]]) -> str:
    declared = {str(b.get("strand") or "").strip() for b in blocks} - {""}
    if declared == {"+"} or declared == {"-"}:
        return declared.pop()
    spans = [s for s in (_genomic_span(b) for b in blocks) if s]
    if len(spans) < 2:
        return ""
    return "+" if spans[0][0] < spans[-1][0] else "-"


def _exon_rows(model: Dict[str, Any]) -> List[Dict[str, Any]]:
    out = []
    for i, b in enumerate(model.get("blocks") or model.get("exons") or [], start=1):
        genomic = _genomic_span(b)
        aa_start = b.get("start", b.get("protein_start_aa"))
        aa_end = b.get("end", b.get("protein_end_aa"))
        out.append({
            "label": b.get("label") or b.get("exon_label")
                     or f"E{b.get('exon_number', b.get('transcript_exon_number', i))}",
            "exon_number": _int(b.get("exon_number", b.get("transcript_exon_number", i)), i),
            "genomic": genomic,
            "aa_start": _int(aa_start) if aa_start is not None else None,
            "aa_end": _int(aa_end) if aa_end is not None else None,
            "coding": str(b.get("coding_status") or "coding").lower() == "coding",
            "source": b,
        })
    return out


def classify_transcript_models(models: Sequence[Dict[str, Any]]) -> List[Dict[str, Any]]:
    ordered = sorted(models, key=lambda m: (0 if m.get("is_primary") else 1,
                                            str(m.get("protein_id") or "")))
    primary = next((m for m in ordered if m.get("is_primary")),
                   ordered[0] if ordered else None)
    if not primary:
        return []
    p_exons = _exon_rows(primary)
    p_spans = [e["genomic"] for e in p_exons if e["genomic"]]
    p_exact = {e["genomic"]: e for e in p_exons if e["genomic"]}

    rows: List[Dict[str, Any]] = []
    for model in ordered:
        exons = _exon_rows(model)
        is_primary = bool(model.get("is_primary")) or model is primary
        spans = [e["genomic"] for e in exons if e["genomic"]]
        for e in exons:
            if is_primary:
                e["kind"] = "primary"
            elif not e["genomic"]:
                e["kind"] = "shared"
            elif e["genomic"] in p_exact:
                e["kind"] = "shared"
            elif any(_overlaps(e["genomic"], p) for p in p_spans):
                e["kind"] = "shifted"
            else:
                e["kind"] = "alternative"  # refined below once skips are known
            e["reference_label"] = ""
            if e["genomic"] and e["kind"] in ("shared", "shifted"):
                match = next((p for p in p_exons
                              if p["genomic"] and _overlaps(e["genomic"], p["genomic"])), None)
                e["reference_label"] = match["label"] if match else ""

        missing: List[Dict[str, Any]] = []
        if not is_primary:
            for idx, p in enumerate(p_exons):
                if not p["genomic"] or any(_overlaps(p["genomic"], g) for g in spans):
                    continue
                upstream = next((e for e in reversed(exons)
                                 if e["genomic"] and idx > 0
                                 and _overlaps(e["genomic"], p_exons[idx - 1]["genomic"] or (0, 0))),
                                None)
                missing.append({
                    "label": p["label"],
                    "primary_start_aa": p["aa_start"], "primary_end_aa": p["aa_end"],
                    "primary_genomic": p["genomic"],
                    "junction_aa": upstream["aa_end"] if upstream else 1,
                })
            # An alternative exon that fills the genomic window of a skipped primary
            # exon replaces it; one that fills no such window is a true insertion.
            by_genomic = sorted([e for e in exons if e["genomic"]],
                                key=lambda e: e["genomic"][0])
            for pos, e in enumerate(by_genomic):
                if e["kind"] != "alternative":
                    continue
                lo = next((by_genomic[k]["genomic"][1] for k in range(pos - 1, -1, -1)
                           if by_genomic[k]["kind"] != "alternative"), None)
                hi = next((by_genomic[k]["genomic"][0]
                           for k in range(pos + 1, len(by_genomic))
                           if by_genomic[k]["kind"] != "alternative"), None)
                replaced = [m for m in missing if m["primary_genomic"]
                            and (lo is None or m["primary_genomic"][0] >= lo)
                            and (hi is None or m["primary_genomic"][1] <= hi)]
                if replaced:
                    e["reference_label"] = ", ".join(m["label"] for m in replaced)
                else:
                    e["kind"] = "inserted"

        alt_terminus = []
        if not is_primary and exons and p_exons:
            for side, mine, theirs in (("N", exons[0], p_exons[0]),
                                       ("C", exons[-1], p_exons[-1])):
                if mine["genomic"] and theirs["genomic"] \
                        and not _overlaps(mine["genomic"], theirs["genomic"]):
                    alt_terminus.append(side)

        counts: Dict[str, int] = {}
        for e in exons:
            counts[e["kind"]] = counts.get(e["kind"], 0) + 1
        n_differing = sum(v for k, v in counts.items() if k not in ("primary", "shared"))
        rows.append({
            "protein_id": model.get("protein_id") or "",
            "transcript_id": model.get("transcript_id") or "",
            "protein_length": _int(model.get("protein_length")),
            "curation_status": model.get("curation_status") or "predicted",
            "role": model.get("role") or ("primary" if is_primary else "alternative"),
            "is_primary": is_primary,
            "strand": derive_strand([e["source"] for e in exons]),
            "exons": exons,
            "missing": missing,
            "alt_terminus": alt_terminus,
            "counts": counts,
            "n_differing": n_differing,
            "differs": bool(n_differing or missing or alt_terminus),
        })
    return rows


def exon_identity_table(rows: Sequence[Dict[str, Any]]) -> List[Dict[str, Any]]:
    out = []
    for model in rows:
        for e in model["exons"]:
            genomic = e["genomic"] or ("", "")
            out.append({
                "protein_id": model["protein_id"],
                "transcript_id": model["transcript_id"],
                "model_role": model["role"],
                "exon_label": e["label"], "exon_number": e["exon_number"],
                "genomic_start": genomic[0], "genomic_end": genomic[1],
                "protein_start_aa": e["aa_start"], "protein_end_aa": e["aa_end"],
                "exon_identity": e["kind"],
                "reference_exon_label": e.get("reference_label", ""),
                "junction_aa": "",
            })
        for m in model["missing"]:
            genomic = m["primary_genomic"] or ("", "")
            out.append({
                "protein_id": model["protein_id"],
                "transcript_id": model["transcript_id"],
                "model_role": model["role"],
                "exon_label": "", "exon_number": "",
                "genomic_start": genomic[0], "genomic_end": genomic[1],
                "protein_start_aa": "", "protein_end_aa": "",
                "exon_identity": "missing",
                "reference_exon_label": m["label"],
                "junction_aa": m["junction_aa"],
            })
    return out


def _strand_phrase(strand: str) -> str:
    normalized = normalize_strand(strand)
    if normalized == PLUS:
        return "plus strand"
    if normalized == MINUS:
        return "minus strand"
    return "strand not annotated"


def _strand_sentence(strand: str) -> str:
    normalized = normalize_strand(strand)
    if normalized == PLUS:
        return "transcribed left to right on the plus strand"
    if normalized == MINUS:
        return "transcribed right to left on the minus strand"
    return "transcription direction not annotated"


def _margin_chars(ax, x: float, fontsize: float, right: float = 0.0) -> int:
    lo, hi = ax.get_xlim()
    pts_per_unit = (_axes_size_in(ax)[0] * 72.0) / max(1e-9, hi - lo)
    return max(10, int(((right - x) * pts_per_unit - 5.0) / (0.62 * fontsize)))


def _stacked_labels(ax, x: float, y: float,
                    lines: Sequence[Tuple[str, str, str]], *, right: float = 0.0,
                    step: float = 0.215) -> None:
    rows = [(t, c, w) for t, c, w in lines if t]
    chars = _margin_chars(ax, x, ps.FONT["small"], right)
    top = y + step * (len(rows) - 1) / 2.0
    for i, (text, colour, weight) in enumerate(rows):
        ax.text(x, top - i * step, _clip(text, chars), ha="left", va="center",
                fontsize=ps.FONT["small"], color=colour, fontweight=weight)


def _model_row_labels(ax, x: float, y: float, model: Dict[str, Any],
                      *details: str, right: float = 0.0) -> None:
    marker = "  (primary)" if model["is_primary"] else ""
    lines = [(f"{model['protein_id']}{marker}", ps.INK,
              "bold" if model["is_primary"] else "normal"),
             (model["transcript_id"], ps.MUTED, "normal")]
    lines += [(d, ps.MUTED, "normal") for d in details if d]
    _stacked_labels(ax, x, y, lines, right=right)


def plot_transcript_exon_structure(fig_dir: Path, stem: str, *, gene_symbol: str,
                                   transcripts: Sequence[Dict[str, Any]],
                                   species_name: str = "",
                                   footnote: Optional[str] = None) -> bool:
    models = list(transcripts)
    if not models:
        return False
    rows = classify_transcript_models(models)
    if not rows:
        return False
    genomic = [e["genomic"] for m in rows for e in m["exons"] if e["genomic"]]
    if not genomic:
        return False
    lo = min(g[0] for g in genomic)
    hi = max(g[1] for g in genomic)
    pad = max(1, int((hi - lo) * 0.015))
    lo, hi = lo - pad, hi + pad
    max_aa = max((m["protein_length"] for m in rows), default=0) or max(
        (e["aa_end"] or 0) for m in rows for e in m["exons"])
    n = len(rows)
    primary = rows[0]
    strand = primary["strand"]
    has_non_coding = any(not e["coding"] for m in rows for e in m["exons"])

    # Each row carries a three-line label, so the panels are sized from the number
    # of rows: a row shorter than roughly 0.45 in cannot hold that label legibly.
    fig, (ax_tx, ax_pr) = plt.subplots(
        2, 1, figsize=(12.8, min(13.6, max(6.4, 0.95 * n + 4.2))),
        gridspec_kw={"height_ratios": [1, 1], "hspace": 0.30})

    # ---------------- Panel A — transcript structure (nucleotides) ------------
    left_frac = 0.32
    span = hi - lo
    tx_left = lo - span * left_frac
    ax_tx.set_xlim(tx_left, hi)
    ax_tx.set_ylim(-0.95, n - 0.05)
    label_x = tx_left + span * 0.01
    for i, model in enumerate(rows):
        y = n - 1 - i
        exons = [e for e in model["exons"] if e["genomic"]]
        if not exons:
            continue
        r_lo = min(e["genomic"][0] for e in exons)
        r_hi = max(e["genomic"][1] for e in exons)
        ax_tx.plot([r_lo, r_hi], [y, y], color=AXIS_GREY, lw=0.9, zorder=1)
        for e in exons:
            s, t = e["genomic"]
            coding = e["coding"]
            ax_tx.add_patch(plt.Rectangle(
                (s, y - 0.17), max(span * 0.0012, t - s), 0.34,
                facecolor=PRIMARY_EXON_COLOR if (model["is_primary"] and coding)
                else (SHARED_EXON_COLOR if coding else NON_CODING_EXON_COLOR),
                edgecolor=ps.INK, linewidth=0.45, zorder=4))
        arrow_x, tip = ((r_hi, r_hi + span * 0.022) if not is_reverse(strand)
                        else (r_lo, r_lo - span * 0.022))
        ax_tx.annotate("", xy=(tip, y), xytext=(arrow_x, y), annotation_clip=False,
                       arrowprops=dict(arrowstyle="-|>", color=ps.MUTED, lw=0.8,
                                       shrinkA=0, shrinkB=0))
        _model_row_labels(ax_tx, label_x, y, model,
                          f"{model['curation_status']} · {len(exons)} coding exons"
                          + (f" · {_strand_phrase(strand)}" if model["is_primary"]
                             else ""),
                          right=lo)
    _place_block_labels(ax_tx, [(e["genomic"][0], e["genomic"][1], e["label"])
                                for e in primary["exons"] if e["genomic"]],
                        y=n - 1, height=0.17, fontsize=ps.FONT["small"],
                        color="white", levels=2, above_only=True)
    ax_tx.set_yticks([])
    ax_tx.spines[["top", "right", "left"]].set_visible(False)
    _readable_nucleotide_axis(ax_tx, lo, hi,
                              "Genomic coordinate on the assembly (kb, NCBI RefSeq "
                              "annotation)")
    ps.panel_label(ax_tx, "A", x=-0.005, y=1.0)
    handles_a = [ps.legend_patch(PRIMARY_EXON_COLOR, "coding exon of the primary transcript"),
                 ps.legend_patch(SHARED_EXON_COLOR, "coding exon of an alternative transcript")]
    if has_non_coding:
        handles_a.append(ps.legend_patch(NON_CODING_EXON_COLOR, "non-coding exon"))
    handles_a.append(ps.legend_line(ps.MUTED, "intron / transcription direction", lw=0.9))
    shared_legend(ax_tx, handles_a, ncol=2, loc="upper left",
                  bbox=(0.0, _below_axes(ax_tx, 30)))

    # ---------------- Panel B — translated protein product (amino acids) -----
    pr_left = -max_aa * 0.42
    ax_pr.set_xlim(pr_left, max_aa * 1.06)
    ax_pr.set_ylim(-0.95, n - 0.05)
    label_x = pr_left + max_aa * 0.01
    kinds_seen = set()
    for i, model in enumerate(rows):
        y = n - 1 - i
        for e in model["exons"]:
            if e["aa_start"] is None or e["aa_end"] is None:
                continue
            kinds_seen.add(e["kind"])
            ax_pr.add_patch(plt.Rectangle(
                (e["aa_start"], y - 0.17), max(1, e["aa_end"] - e["aa_start"]), 0.34,
                facecolor=EXON_IDENTITY_COLOR[e["kind"]], edgecolor=ps.INK,
                linewidth=0.45, zorder=4))
        for m in model["missing"]:
            kinds_seen.add("missing")
            ax_pr.scatter([m["junction_aa"]], [y - 0.30], marker="v", s=26,
                          color=MISSING_MARK_COLOR, edgecolor=ps.INK, linewidth=0.4,
                          zorder=6, clip_on=False)
        for side in model["alt_terminus"]:
            kinds_seen.add("terminus")
            x = 1 if side == "N" else model["protein_length"]
            ax_pr.plot([x, x], [y - 0.24, y + 0.24], color=TERMINUS_COLOR, lw=1.8,
                       zorder=6)
        ax_pr.text(max_aa * 1.01, y, f"{model['protein_length']} aa", ha="left",
                   va="center", fontsize=ps.FONT["small"], color=ps.INK)
        _model_row_labels(ax_pr, label_x, y, model,
                          "reference model" if model["is_primary"]
                          else _difference_phrase(model))
    _place_block_labels(ax_pr, [(e["aa_start"], e["aa_end"], e["label"])
                                for e in primary["exons"]
                                if e["aa_start"] is not None],
                        y=n - 1, height=0.17, fontsize=ps.FONT["small"],
                        color="white", levels=2, above_only=True)
    ax_pr.set_yticks([])
    ax_pr.spines[["top", "right", "left"]].set_visible(False)
    _readable_aa_axis(ax_pr, max_aa,
                      "Amino-acid position in each translated protein model (aa)")
    ps.panel_label(ax_pr, "B", x=-0.005, y=1.0)
    handles_b = [ps.legend_patch(EXON_IDENTITY_COLOR[k], EXON_IDENTITY_LABEL[k])
                 for k in ("primary", "shared", "shifted", "alternative", "inserted")
                 if k in kinds_seen]
    if "missing" in kinds_seen:
        handles_b.append(_legend_marker("v", EXON_IDENTITY_LABEL["missing"],
                                        MISSING_MARK_COLOR))
    if "terminus" in kinds_seen:
        handles_b.append(ps.legend_line(TERMINUS_COLOR, EXON_IDENTITY_LABEL["terminus"],
                                        lw=1.8))
    ncol = 2
    shared_legend(ax_pr, handles_b, ncol=ncol, loc="upper left",
                  bbox=(0.0, _below_axes(ax_pr, 30)))

    species = f"{species_name} · " if species_name else ""
    figure_title(
        ax_tx,
        f"{gene_symbol} — transcript structure and translated protein product",
        subtitle=f"{species}{n} annotated transcript(s) · primary "
                 f"{primary['transcript_id']} → {primary['protein_id']} · "
                 f"{primary['protein_length']} aa · {_strand_sentence(strand)} · "
                 f"panel A in nucleotides, panel B in amino acids",
        note="Which coding exons does each transcript use, and which part of each "
             "translated protein do they encode? Exon identity in panel B is "
             "compared on genomic CDS coordinates; differences are descriptive "
             "annotation differences, not validated splicing events.")
    _footnote(ax_pr, footnote or "",
              points=30 + 12 * (1 + (len(handles_b) - 1) // ncol) + 8)
    save_figure_all_formats(fig, fig_dir, stem)
    return True


def _difference_phrase(model: Dict[str, Any], *, include_shared: bool = False) -> str:
    parts = []
    kinds = (("shared", "shared"),) if include_shared else ()
    for kind, word in kinds + (("shifted", "shifted"),
                               ("alternative", "alternative"),
                               ("inserted", "inserted")):
        if model["counts"].get(kind):
            parts.append(f"{model['counts'][kind]} {word}")
    if model["missing"]:
        n = len(model["missing"])
        parts.append(f"{n} exon absent" if n == 1 else f"{n} exons absent")
    if model["alt_terminus"]:
        parts.append(f"alt. {'/'.join(model['alt_terminus'])}-terminus")
    if not model["n_differing"] and not model["missing"] and not model["alt_terminus"]:
        parts.append("identical to the primary")
    return " · ".join(parts)


# --------------------------------------------------------------------------- #
# within-species protein isoform alignment
#
# The two alignment figures below are the matplotlib counterpart of the
# interactive figure source in
# ``webapp/frontend/src/pages/viewers/alignmentFigure.js`` (fullAlignmentFigureSvg
# and candidateAlignmentFigureSvg): same rows, same colour semantics, same
# derived quantities, so the static export and the on-screen figure agree.
# --------------------------------------------------------------------------- #
def alignment_rows(sequences: Sequence[Dict[str, Any]],
                   primary_id: str = "") -> List[Dict[str, Any]]:
    rows = []
    for s in sequences:
        seq = s.get("aligned_sequence") or s.get("seq") or ""
        pid = s.get("protein_id") or ""
        rows.append({
            "protein_id": pid,
            "transcript_id": s.get("transcript_id") or "",
            "is_primary": bool(s.get("is_primary")) or (bool(primary_id) and pid == primary_id),
            "seq": seq,
            "protein_length": len(seq.replace("-", "")),
            "curation_status": "curated" if pid.startswith("NP_") else "predicted",
        })
    rows.sort(key=lambda r: (0 if r["is_primary"] else 1, r["protein_id"]))
    return rows


def _column_stats(rows: Sequence[Dict[str, Any]], n_cols: int):
    variable = [False] * n_cols
    gapped = [False] * n_cols
    seqs = [r["seq"] for r in rows]
    for c in range(n_cols):
        first = None
        for seq in seqs:
            ch = seq[c] if c < len(seq) else "-"
            if ch == "-":
                gapped[c] = True
            if first is None:
                first = ch
            elif ch != first:
                variable[c] = True
    return variable, gapped


def _runs_of(flags: Sequence[bool], min_len: int = 1) -> List[Tuple[int, int]]:
    out, start = [], -1
    for c, flag in enumerate(flags):
        if flag:
            if start < 0:
                start = c
        elif start >= 0:
            if c - start >= min_len:
                out.append((start, c - 1))
            start = -1
    if start >= 0 and len(flags) - start >= min_len:
        out.append((start, len(flags) - 1))
    return out


def _identity_pct(seq: str, primary: str) -> Optional[int]:
    same = compared = 0
    for c, p in enumerate(primary):
        s = seq[c] if c < len(seq) else "-"
        if p == "-" or s == "-":
            continue
        compared += 1
        same += int(p == s)
    return round(100 * same / compared) if compared else None


def aa_to_column(primary_aligned: str, aa: Optional[int]) -> Optional[int]:
    if not primary_aligned or aa is None:
        return None
    residue = 0
    for c, ch in enumerate(primary_aligned):
        if ch != "-":
            residue += 1
            if residue == aa:
                return c
    return None


def _candidate_columns(rows, candidates) -> List[Dict[str, Any]]:
    primary = next((r for r in rows if r["is_primary"]), rows[0] if rows else None)
    if not primary:
        return []
    out = []
    for i, cand in enumerate(candidates, start=1):
        aa_start = cand.get("aa_start", cand.get("start"))
        aa_end = cand.get("aa_end", cand.get("end"))
        c0 = aa_to_column(primary["seq"], _int(aa_start) if aa_start is not None else None)
        c1 = aa_to_column(primary["seq"], _int(aa_end) if aa_end is not None else None)
        if c0 is None or c1 is None:
            continue
        out.append({
            "label": cand.get("label") or cand.get("rank_label") or f"C{i}",
            "aa_start": _int(aa_start), "aa_end": _int(aa_end),
            "col_start": c0, "col_end": c1,
            "candidate_id": cand.get("candidate_id") or cand.get("id") or "",
            "affected": candidate_affected_proteins(cand),
        })
    return out


def candidate_affected_proteins(cand: Dict[str, Any]) -> List[str]:
    listed = cand.get("affected_proteins") or cand.get("affected") or []
    if listed:
        return sorted({str(p) for p in listed})
    pairs = (cand.get("protein_isoform_evidence") or {}).get("supporting_isoform_pairs") or []
    found = set()
    for pair in pairs:
        for key in ("protein_a", "protein_b"):
            if pair.get(key):
                found.add(str(pair[key]))
    return sorted(found)


def plot_isoform_alignment_overview(fig_dir: Path, stem: str, *, gene_symbol: str,
                                    sequences: Sequence[Dict[str, Any]],
                                    candidates: Sequence[Dict[str, Any]] = (),
                                    species_name: str = "", primary_id: str = "",
                                    tool: str = "MAFFT",
                                    alignment_length: Optional[int] = None,
                                    footnote: Optional[str] = None) -> bool:
    rows = alignment_rows(sequences, primary_id)
    if len(rows) < 2:
        return False
    n_cols = alignment_length or max((len(r["seq"]) for r in rows), default=0)
    if not n_cols:
        return False
    primary = next((r for r in rows if r["is_primary"]), rows[0])
    primary_seq = primary["seq"]
    variable, gapped = _column_stats(rows, n_cols)
    min_block = max(3, round(n_cols * 0.004))
    var_blocks = _runs_of(variable, min_block)
    gap_blocks = _runs_of(gapped, min_block)
    cands = _candidate_columns(rows, candidates)

    n = len(rows)
    fig, ax = plt.subplots(figsize=(13, max(3.8, 0.62 * n + 3.2)))
    left = -0.235 * n_cols
    id_x0 = n_cols * 1.035
    id_w = n_cols * 0.085
    ax.set_xlim(left, id_x0 + id_w + n_cols * 0.03)
    ax.set_ylim(-1.0, n + 0.9)
    rh = 0.34  # half-height of an isoform row

    # global variability / conservation track
    y_var = n + 0.35
    ax.add_patch(plt.Rectangle((0, y_var - 0.18), n_cols, 0.36, facecolor="#F7F9FB",
                               edgecolor=AXIS_GREY, lw=0.5, zorder=2))
    for a, b in _runs_of(variable, 1):
        ax.add_patch(plt.Rectangle((a, y_var - 0.18), max(1, b - a + 1), 0.36,
                                   facecolor=ALN_VARIABLE, edgecolor="none",
                                   alpha=0.9, zorder=3))
    ax.text(-0.006 * n_cols, y_var, "Variable columns", ha="right", va="center",
            fontsize=ps.FONT["small"], color=ps.MUTED)

    # exploratory candidate overlays — subtle, must never hide the alignment
    for cand in cands:
        x0, x1 = cand["col_start"], cand["col_end"] + 1
        ax.axvspan(x0, x1, ymin=0.02, ymax=0.92, color=CANDIDATE_COLOR, alpha=0.12,
                   lw=0, zorder=1)
        for x in (x0, x1):
            ax.plot([x, x], [-0.55, n - 0.35], color=CANDIDATE_EDGE, lw=0.7,
                    ls=(0, (3, 2)), zorder=2)
        ax.text((x0 + x1) / 2.0, n - 0.30, cand["label"], ha="center", va="bottom",
                fontsize=ps.FONT["small"], color=CANDIDATE_EDGE)

    for i, row in enumerate(rows):
        y = n - 1 - i
        seq = row["seq"]
        ax.add_patch(plt.Rectangle((0, y - rh), n_cols, 2 * rh, facecolor=ALN_GAP,
                                   edgecolor="none", zorder=3))
        present = [(c < len(seq) and seq[c] != "-") for c in range(n_cols)]
        colour = ALN_RESIDUE_PRIMARY if row["is_primary"] else ALN_RESIDUE
        for a, b in _runs_of(present, 1):
            ax.add_patch(plt.Rectangle((a, y - rh), b - a + 1, 2 * rh, facecolor=colour,
                                       edgecolor="none",
                                       alpha=0.9 if row["is_primary"] else 0.72, zorder=4))
        differs = [variable[c] and (seq[c] if c < len(seq) else "-") !=
                   (primary_seq[c] if c < len(primary_seq) else "-") for c in range(n_cols)]
        for a, b in _runs_of(differs, 1):
            ax.add_patch(plt.Rectangle((a, y - rh), b - a + 1, 2 * rh,
                                       facecolor=ALN_VARIABLE, edgecolor="none",
                                       alpha=0.95, zorder=5))
        ax.add_patch(plt.Rectangle((0, y - rh), n_cols, 2 * rh, facecolor="none",
                                   edgecolor=AXIS_GREY, lw=0.4, zorder=6))

        marker = "  (primary)" if row["is_primary"] else ""
        ax.text(left, y + 0.16, f"{row['protein_id']}{marker}", ha="left", va="center",
                fontsize=ps.FONT["small"], color=ps.INK,
                fontweight="bold" if row["is_primary"] else "normal")
        ax.text(left, y - 0.16,
                f"{row['transcript_id']} · {row['protein_length']} aa · {row['curation_status']}",
                ha="left", va="center", fontsize=ps.FONT["small"], color=ps.MUTED)

        pct = _identity_pct(seq, primary_seq)
        if pct is None:
            ax.text(id_x0, y, "—", ha="left", va="center", fontsize=ps.FONT["small"],
                    color=ps.MUTED)
        else:
            ax.add_patch(plt.Rectangle((id_x0, y - rh * 0.7), id_w * pct / 100.0,
                                       1.4 * rh, facecolor=ALN_IDENTITY, alpha=0.30,
                                       edgecolor="none", zorder=4))
            ax.text(id_x0 + id_w * 0.04, y, f"{pct}%", ha="left", va="center",
                    fontsize=ps.FONT["small"], color=ps.INK, zorder=6)
    ax.text(id_x0, n - 0.30, "identity to primary", ha="left", va="bottom",
            fontsize=ps.FONT["small"], color=ps.MUTED)

    ax.set_yticks([])
    ax.spines[["top", "right", "left"]].set_visible(False)
    _readable_aa_axis(ax, n_cols, f"Alignment column ({tool}, 1-based)", lo=0)
    figure_title(
        ax, f"{gene_symbol} — within-species protein isoform alignment",
        subtitle=f"{species_name + ' · ' if species_name else ''}{n} protein models · "
                 f"{n_cols} alignment columns · {tool} · primary {primary['protein_id']} · "
                 f"{len(var_blocks)} major variable block(s) · {len(gap_blocks)} major gap block(s)",
        note="Within-species isoform alignment, not a cross-species conservation analysis; "
             "no isoform difference shown here is a validated splicing event.")
    handles = [ps.legend_patch(ALN_RESIDUE_PRIMARY, "primary protein (aligned residues)"),
               ps.legend_patch(ALN_RESIDUE, "alternative isoform (aligned residues)"),
               ps.legend_patch(ALN_GAP, "gap in this isoform"),
               ps.legend_patch(ALN_VARIABLE, "residue differing from primary"),
               ps.legend_patch(ALN_IDENTITY, "identity to primary")]
    if cands:
        handles.append(ps.legend_patch(CANDIDATE_COLOR, "exploratory candidate interval"))
    ncol = 3
    shared_legend(ax, handles, ncol=ncol, loc="upper center",
                  bbox=(0.5, _below_axes(ax, 44)))
    _footnote(ax, footnote or "", points=44 + 12 * (1 + (len(handles) - 1) // ncol) + 6)
    save_figure_all_formats(fig, fig_dir, stem)
    return True






def plot_protein_exon_architecture(fig_dir: Path, stem: str, *, gene_symbol: str,
                                   primary_id: str, exon_blocks: Sequence[Dict[str, Any]],
                                   candidate_regions: Sequence[Dict[str, Any]] = (),
                                   domains: Sequence[Dict[str, Any]] = (),
                                   species_name: str = "",
                                   protein_length: Optional[int] = None,
                                   footnote: Optional[str] = None) -> bool:
    blocks = []
    for b in exon_blocks:
        start = b.get("protein_start_aa", b.get("start"))
        end = b.get("protein_end_aa", b.get("end"))
        if start is None or end is None:
            continue
        label = b.get("label") or f"E{b.get('exon_number', '')}".rstrip()
        blocks.append((_int(start), _int(end), label if label != "E" else ""))
    if not blocks:
        return False
    maxaa = protein_length or max(e for _, e, _ in blocks)
    cands = []
    for c in candidate_regions:
        s = c.get("start_aa", c.get("start"))
        e = c.get("end_aa", c.get("end"))
        if s is None or e is None:
            continue
        cands.append((_int(s), _int(e), c.get("id") or c.get("candidate_id") or ""))
    dom_rows = number_domain_instances(domains) if domains else []
    boundaries = sorted({s for s, _, _ in blocks if s > 1})

    lanes = 2 + int(bool(cands)) + int(bool(dom_rows))
    fig, ax = plt.subplots(figsize=(12.5, max(3.4, 1.05 * lanes + 2.3)))
    left = -maxaa * 0.16
    ax.set_xlim(left, maxaa * 1.03)
    ax.set_ylim(-0.55, lanes - 0.45)
    h = 0.18
    y_bound = 0.0
    y_exon = 1.0
    y_cand = 2.0 if cands else None
    y_dom = (2.0 if not cands else 3.0) if dom_rows else None

    handles = []
    ax.plot([1, maxaa], [y_bound, y_bound], color=AXIS_GREY, lw=0.7, zorder=1)
    for pos in boundaries:
        ax.plot([pos, pos], [y_bound - h, y_bound + h], color=BOUNDARY_TICK_COLOR,
                lw=1.1, zorder=5)
    ax.text(left, y_bound, "Coding-exon boundaries", ha="left", va="center",
            fontsize=ps.FONT["label"], color=ps.INK)

    for i, (s, e, _lbl) in enumerate(blocks):
        ax.add_patch(plt.Rectangle((s, y_exon - h), max(1, e - s), 2 * h,
                                   facecolor=EXON_BLOCK_COLORS[i % 2], edgecolor=ps.INK,
                                   linewidth=0.4, zorder=4))
    _place_block_labels(ax, [(s, e, lbl) for s, e, lbl in blocks], y=y_exon, height=h,
                        fontsize=ps.FONT["small"])
    ax.text(left, y_exon, "Coding exons", ha="left", va="center",
            fontsize=ps.FONT["label"], color=ps.INK)
    handles += [ps.legend_patch(EXON_BLOCK_COLORS[1], "coding exon projected onto the protein"),
                ps.legend_line(BOUNDARY_TICK_COLOR, "internal coding-exon boundary", lw=1.1)]

    if cands:
        for s, e, _cid in cands:
            ax.add_patch(plt.Rectangle((s, y_cand - h), max(1, e - s), 2 * h,
                                       facecolor=CANDIDATE_COLOR, edgecolor=CANDIDATE_EDGE,
                                       linewidth=0.5, alpha=0.85, zorder=4))
            ax.axvspan(s, e, ymin=0.04, ymax=0.96, color=CANDIDATE_COLOR, alpha=0.10,
                       lw=0, zorder=1)
        _place_block_labels(ax, cands, y=y_cand, height=h, fontsize=ps.FONT["small"])
        ax.text(left, y_cand, "Exploratory candidates", ha="left", va="center",
                fontsize=ps.FONT["label"], color=ps.INK)
        handles.append(ps.legend_patch(CANDIDATE_COLOR,
                                       "exploratory candidate region (not validated)"))

    if dom_rows:
        for d in dom_rows:
            ax.add_patch(plt.Rectangle((d["start"], y_dom - h),
                                       max(1, d["end"] - d["start"]), 2 * h,
                                       facecolor=d["color"], edgecolor=ps.INK,
                                       linewidth=0.4, zorder=4))
        _place_block_labels(ax, [(d["start"], d["end"], d["short"]) for d in dom_rows],
                            y=y_dom, height=h, fontsize=ps.FONT["small"], color="white")
        ax.text(left, y_dom, "Representative domains", ha="left", va="center",
                fontsize=ps.FONT["label"], color=ps.INK)
        handles.append(ps.legend_patch(dom_rows[0]["color"], "representative InterPro domain"))

    ax.set_yticks([])
    ax.spines[["top", "right", "left"]].set_visible(False)
    _readable_aa_axis(ax, maxaa, f"Amino-acid position on {primary_id} (aa)")
    species = f"{species_name} · " if species_name else ""
    figure_title(
        ax, f"{gene_symbol} — primary exon-to-protein projection",
        subtitle=f"{species}{primary_id} · {maxaa} aa · {len(blocks)} coding exons · "
                 f"{len(boundaries)} internal coding-exon boundaries",
        note="Exploratory candidate regions are not biologically validated events.")
    ncol = min(3, max(1, len(handles)))
    shared_legend(ax, handles, ncol=ncol, loc="upper center",
                  bbox=(0.5, _below_axes(ax, 44)))
    _footnote(ax, footnote or "", points=44 + 12 * (1 + (len(handles) - 1) // ncol) + 6)
    save_figure_all_formats(fig, fig_dir, stem)
    return True


def plot_transcript_model_comparison(fig_dir: Path, stem: str, *, gene_symbol: str,
                                     models: Sequence[Dict[str, Any]],
                                     candidate: Optional[Dict[str, Any]] = None,
                                     diff_only: bool = False,
                                     species_name: str = "",
                                     footnote: Optional[str] = None) -> bool:
    rows = classify_transcript_models(models)
    if not rows:
        return False
    primary = rows[0]
    shown = [m for m in rows if m["is_primary"] or m["differs"]] if diff_only else rows
    if not shown:
        return False

    max_aa = max((m["protein_length"] for m in shown), default=1) or 1
    n = len(shown)
    fig, ax = plt.subplots(figsize=(12.8, max(4.0, min(13.6, n * 0.78 + 3.0))))
    left = -max_aa * 0.44
    ax.set_xlim(left, max_aa * 1.10)
    ax.set_ylim(-0.85, n - 0.15)
    h = 0.17
    label_x = left + max_aa * 0.008

    if candidate is not None:
        cs, ce = _int(candidate.get("start")), _int(candidate.get("end"))
        ax.axvspan(cs, ce, color=CANDIDATE_COLOR, alpha=0.13, lw=0, zorder=1)
        for x in (cs, ce):
            ax.plot([x, x], [-0.85, n - 0.55], color=CANDIDATE_EDGE, lw=0.7, ls="--",
                    zorder=2)
        ax.text((cs + ce) / 2.0, n - 0.52, candidate.get("id") or "candidate",
                ha="center", va="bottom", fontsize=ps.FONT["small"], color=CANDIDATE_EDGE)

    kinds_seen = set()
    for i, model in enumerate(shown):
        y = n - 1 - i
        for e in model["exons"]:
            if e["aa_start"] is None or e["aa_end"] is None:
                continue
            kinds_seen.add(e["kind"])
            emphasised = e["kind"] in ("shifted", "alternative", "inserted")
            ax.add_patch(plt.Rectangle(
                (e["aa_start"], y - h), max(1, e["aa_end"] - e["aa_start"]), 2 * h,
                facecolor=EXON_IDENTITY_COLOR[e["kind"]], edgecolor=ps.INK,
                linewidth=0.9 if emphasised else 0.45, zorder=4))
        for m in model["missing"]:
            kinds_seen.add("missing")
            ax.scatter([m["junction_aa"]], [y - h - 0.13], marker="v", s=30,
                       color=MISSING_MARK_COLOR, edgecolor=ps.INK, linewidth=0.4,
                       zorder=6)
        for side in model["alt_terminus"]:
            kinds_seen.add("terminus")
            x = 1 if side == "N" else model["protein_length"]
            ax.plot([x, x], [y - h - 0.06, y + h + 0.06], color=TERMINUS_COLOR,
                    lw=1.8, zorder=6)
        ax.text(max_aa * 1.02, y, f"{model['protein_length']} aa", ha="left",
                va="center", fontsize=ps.FONT["small"], color=ps.INK)
        shared = model["counts"].get("shared", 0)
        _model_row_labels(ax, label_x, y, model,
                          f"{model['curation_status']} · {len(model['exons'])} exons"
                          + (f" · {shared} shared" if shared else ""),
                          "reference model" if model["is_primary"]
                          else _difference_phrase(model))
    _place_block_labels(ax, [(e["aa_start"], e["aa_end"], e["label"])
                             for e in primary["exons"] if e["aa_start"] is not None],
                        y=n - 1, height=h, fontsize=ps.FONT["small"], color="white",
                        levels=2, above_only=True)

    ax.set_yticks([])
    ax.spines[["top", "right", "left"]].set_visible(False)
    _readable_aa_axis(ax, max_aa,
                      "Amino-acid position in each protein model (aa, common scale)")
    handles = [ps.legend_patch(EXON_IDENTITY_COLOR[k], EXON_IDENTITY_LABEL[k])
               for k in ("primary", "shared", "shifted", "alternative", "inserted")
               if k in kinds_seen]
    if "missing" in kinds_seen:
        handles.append(_legend_marker("v", EXON_IDENTITY_LABEL["missing"],
                                      MISSING_MARK_COLOR))
    if "terminus" in kinds_seen:
        handles.append(ps.legend_line(TERMINUS_COLOR, EXON_IDENTITY_LABEL["terminus"],
                                      lw=1.8))
    if candidate is not None:
        handles.append(ps.legend_patch(CANDIDATE_COLOR,
                                       "exploratory candidate region (not validated)"))
    n_diff_blocks = sum(m["n_differing"] for m in rows)
    n_blocks = sum(len(m["exons"]) for m in rows if not m["is_primary"])
    mode = "differences only" if diff_only else "all protein models"
    species = f"{species_name} · " if species_name else ""
    figure_title(
        ax, f"{gene_symbol} — transcript-model comparison ({mode})",
        subtitle=f"{species}{n} of {len(rows)} protein model(s) shown · primary "
                 f"{primary['protein_id']} first · {n_diff_blocks} of {n_blocks} "
                 f"alternative-model exons differ from the primary",
        note="Which coding exons do the alternative protein models share with the "
             "primary, and where do they diverge? Exon identity is compared on "
             "genomic CDS coordinates, so an upstream deletion does not mark "
             "genomically identical downstream exons as altered. Differences are "
             "annotation differences, not validated splicing events.")
    ncol = 2
    shared_legend(ax, handles, ncol=ncol, loc="upper left",
                  bbox=(0.0, _below_axes(ax, 34)))
    _footnote(ax, footnote or "",
              points=34 + 12 * (1 + (len(handles) - 1) // ncol) + 8)
    save_figure_all_formats(fig, fig_dir, stem)
    return True


def plot_selected_candidate_detail(fig_dir: Path, stem: str, *, gene_symbol: str,
                                   primary_id: str, exon_blocks: Sequence[Dict[str, Any]],
                                   candidate: Dict[str, Any]) -> bool:
    blocks = list(exon_blocks)
    if not blocks or not candidate:
        return False
    cs, ce = _int(candidate.get("start")), _int(candidate.get("end"))
    pad = max(20, (ce - cs))
    lo, hi = max(0, cs - pad), ce + pad
    fig, ax = plt.subplots(figsize=(11, 3.0))
    y = 0.0
    ax.add_patch(plt.Rectangle((cs, y - 0.45), max(1, ce - cs), 0.9, facecolor=CANDIDATE_COLOR,
                               edgecolor="#b8690a", linewidth=1.0, alpha=0.3, zorder=2))
    for i, b in enumerate(blocks):
        s, e = _int(b.get("protein_start_aa")), _int(b.get("protein_end_aa"))
        if e < lo or s > hi:
            continue
        col = EXON_BLOCK_COLORS[i % 2]
        ax.add_patch(plt.Rectangle((s, y - 0.25), max(1, e - s), 0.5, facecolor=col,
                                   edgecolor=ps.INK, linewidth=0.5, zorder=3))
        ax.text((s + e) / 2, y, str(b.get("exon_number", "")), va="center", ha="center",
                fontsize=ps.FONT["small"], color=ps.INK)
    ax.set_xlim(lo, hi)
    ax.set_ylim(-0.9, 0.9)
    ax.set_yticks([])
    ax.set_xlabel("Protein position (aa)", fontsize=ps.FONT["label"])
    ax.spines[["top", "right", "left"]].set_visible(False)
    ctype = candidate.get("candidate_type") or "candidate"
    figure_title(ax, f"{gene_symbol} — selected candidate detail ({candidate.get('id','C1')})",
                 subtitle=f"{primary_id}: aa {cs}–{ce} · {ctype} (exploratory, not a validated event)")
    shared_legend(ax, [ps.legend_patch(EXON_BLOCK_COLORS[1], "coding-exon block"),
                       ps.legend_patch(CANDIDATE_COLOR, "selected exploratory candidate")],
                  ncol=2, bbox=(0.5, -0.30))
    save_figure_all_formats(fig, fig_dir, stem)
    return True


# Local-neighbourhood palette. A single-species locus map has no orthology to
# assess, so there is no confidence encoding here — only what the annotation says.
TARGET_GENE_COLOR = "#44546A"
NEIGHBOUR_GENE_COLOR = "#8FA8BF"
PLACEHOLDER_LOCUS_COLOR = "#DCE1E8"
NEIGHBOURHOOD_TABLE_COLUMNS = (
    "position_in_figure", "side", "order_from_target", "locus_symbol", "locus_kind",
    "transcription_direction", "annotation_source",
)

_PLACEHOLDER_SYMBOL = re.compile(r"^(?:LOC\d+|GENE\d+|ENS\w*G\d+)$", re.IGNORECASE)


def is_placeholder_locus(symbol: str) -> bool:
    return not str(symbol or "").strip() or bool(_PLACEHOLDER_SYMBOL.match(str(symbol).strip()))


def neighbourhood_layout(gene_symbol: str, neighbours: Sequence[Dict[str, Any]],
                         target_orientation: str = "+") -> List[Dict[str, Any]]:
    def _symbol(row):
        return (row.get("neighbor_symbol") or row.get("locus_symbol")
                or row.get("gene_symbol") or "")

    up = sorted([r for r in neighbours if str(r.get("side")) == "upstream"],
                key=lambda r: _int(r.get("order")), reverse=True)
    down = sorted([r for r in neighbours if str(r.get("side")) == "downstream"],
                  key=lambda r: _int(r.get("order")))
    out = []
    for row, side in [(r, "upstream") for r in up] + [(None, "target")] \
            + [(r, "downstream") for r in down]:
        if row is None:
            out.append({"symbol": gene_symbol, "side": "target", "order": 0,
                        "kind": "target", "orientation": target_orientation or "+",
                        "source": ""})
            continue
        symbol = _symbol(row)
        out.append({
            "symbol": symbol or "unnamed locus", "side": side,
            "order": _int(row.get("order")),
            "kind": "placeholder_locus" if is_placeholder_locus(symbol)
                    else "annotated_neighbour",
            "orientation": "-" if str(row.get("orientation", "+")).startswith("-") else "+",
            "source": row.get("source") or "",
        })
    return out


def plot_synteny_neighbourhood(fig_dir: Path, stem: str, *, gene_symbol: str,
                               neighbours: Sequence[Dict[str, Any]],
                               species_name: str = "",
                               target_orientation: str = "+",
                               footnote: Optional[str] = None) -> bool:
    loci = neighbourhood_layout(gene_symbol, list(neighbours), target_orientation)
    if len(loci) < 2:
        return False
    n = len(loci)
    n_up = sum(1 for g in loci if g["side"] == "upstream")
    n_down = sum(1 for g in loci if g["side"] == "downstream")
    n_placeholder = sum(1 for g in loci if g["kind"] == "placeholder_locus")

    fig, ax = plt.subplots(figsize=(min(15.0, max(9.5, n * 1.18)), 3.4))
    # The arrow tips reach ±0.39 around each slot centre, so the axes must extend
    # past the outermost slots or the first and last gene lose their point.
    ax.set_xlim(-0.52, n - 0.48)
    ax.set_ylim(-1.25, 0.95)
    ax.plot([-0.39, n - 0.61], [0.0, 0.0], color=AXIS_GREY, lw=1.0, zorder=1)
    colours = {"target": TARGET_GENE_COLOR,
               "annotated_neighbour": NEIGHBOUR_GENE_COLOR,
               "placeholder_locus": PLACEHOLDER_LOCUS_COLOR}
    for i, g in enumerate(loci):
        x = float(i)
        strand = 1 if g["orientation"] == "+" else -1
        ps.gene_arrow(ax, x, 0.0, 0.78, 0.42, strand, colours[g["kind"]])
        is_target = g["kind"] == "target"
        ax.text(x, -0.34, g["symbol"], va="top", ha="center",
                fontsize=ps.FONT["gene"],
                color=ps.INK if g["kind"] != "placeholder_locus" else ps.MUTED,
                fontweight="bold" if is_target else "normal",
                style="italic" if g["kind"] == "placeholder_locus" else "normal")
        if not is_target:
            ax.text(x, -0.56, f"{g['side']} {g['order']}", va="top", ha="center",
                    fontsize=ps.FONT["small"], color=ps.MUTED)
        else:
            ax.text(x, 0.36, "target gene", va="bottom", ha="center",
                    fontsize=ps.FONT["small"], color=ps.INK, fontweight="bold")
    target_x = float(next(i for i, g in enumerate(loci) if g["kind"] == "target"))
    for x0, x1, text in ((-0.39, target_x - 0.5, "upstream"),
                         (target_x + 0.5, n - 0.61, "downstream")):
        if x1 <= x0:
            continue
        ax.plot([x0, x1], [-0.86, -0.86], color=AXIS_GREY, lw=0.8, zorder=1)
        ax.text((x0 + x1) / 2.0, -0.92, text, va="top", ha="center",
                fontsize=ps.FONT["small"], color=ps.MUTED)
    ax.axis("off")

    species = f"{species_name} · " if species_name else ""
    figure_title(
        ax, f"{gene_symbol} · local genomic neighbourhood",
        subtitle=f"{n_up + n_down} flanking loci shown · {n_up} upstream · "
                 f"{n_down} downstream"
                 + (f" · {n_placeholder} placeholder "
                    f"{'locus' if n_placeholder == 1 else 'loci'}"
                    if n_placeholder else ""),
        note=f"Which annotated loci flank {gene_symbol} in this assembly, and in "
             f"which direction is each transcribed? {species}Loci appear in "
             "annotated order along the assembly; the spacing is ordinal and not "
             "to scale.")
    handles = [ps.legend_patch(TARGET_GENE_COLOR, "target gene"),
               ps.legend_patch(NEIGHBOUR_GENE_COLOR, "annotated neighbouring gene")]
    if n_placeholder:
        handles.append(ps.legend_patch(PLACEHOLDER_LOCUS_COLOR,
                                       "placeholder locus (no approved symbol)"))
    handles.append(_legend_marker(">", "arrow points in the transcription direction",
                                  NEIGHBOUR_GENE_COLOR))
    ncol = 2
    shared_legend(ax, handles, ncol=ncol, loc="upper left",
                  bbox=(0.0, _below_axes(ax, 6)))
    _footnote(ax, footnote or "",
              points=6 + 12 * (1 + (len(handles) - 1) // ncol) + 8)
    save_figure_all_formats(fig, fig_dir, stem)
    return True


CANDIDATE_RANKING_TABLE_COLUMNS = (
    "rank", "candidate_id", "candidate_label", "aa_start", "aa_end", "length_aa",
    "n_affected_isoforms", "affected_isoforms", "supporting_comparisons",
    "evidence_score", "evidence_strength", "biological_validation",
)


def candidate_ranking_rows(clusters: Sequence[Dict[str, Any]]) -> List[Dict[str, Any]]:
    rows = []
    for c in clusters:
        start = c.get("aa_start", c.get("start_aa", c.get("start")))
        end = c.get("aa_end", c.get("end_aa", c.get("end")))
        if start is None or end is None:
            continue
        affected = c.get("affected_proteins") or c.get("proteins_involved") or []
        if isinstance(affected, str):
            affected = [p for p in re.split(r"[;,]\s*", affected) if p]
        score = c.get("evidence_score", c.get("overall_score", c.get("score")))
        strength = (c.get("evidence_strength") or c.get("confidence_class")
                    or c.get("confidence") or "")
        rows.append({
            "candidate_id": c.get("candidate_id") or c.get("candidate_cluster_id")
                            or c.get("id") or "",
            "candidate_label": c.get("candidate_label") or c.get("label") or "",
            "aa_start": _int(start), "aa_end": _int(end),
            "length_aa": _int(c.get("length", c.get("length_aa",
                                                    _int(end) - _int(start) + 1))),
            "affected": sorted({str(p) for p in affected}),
            "supporting_comparisons": _int(c.get("supporting_comparisons",
                                                 c.get("support_count"))),
            "evidence_score": _int(score) if score is not None else 0,
            "evidence_strength": str(strength).replace("_", " "),
            "biological_validation": "not validated",
        })
    rows.sort(key=lambda r: (-r["evidence_score"], -r["length_aa"], r["aa_start"]))
    for i, r in enumerate(rows, start=1):
        r["rank"] = i
        if not r["candidate_label"]:
            r["candidate_label"] = f"C{i}"
        r["n_affected_isoforms"] = len(r["affected"])
        r["affected_isoforms"] = ";".join(r["affected"])
    return rows


def plot_evidence_regions_on_protein(fig_dir: Path, stem: str, *, gene_symbol: str,
                                     clusters: Sequence[Dict[str, Any]],
                                     max_aa: Optional[int] = None,
                                     species_name: str = "", protein_id: str = "",
                                     n_models: Optional[int] = None,
                                     footnote: Optional[str] = None) -> bool:
    rows = candidate_ranking_rows(clusters)
    if not rows:
        return False
    n = len(rows)
    top = max((r["evidence_score"] for r in rows), default=1) or 1
    fig, ax = plt.subplots(figsize=(13.4, max(3.2, 0.46 * n + 2.9)))
    ax.set_xlim(0, 1)
    ax.set_ylim(-0.85, n + 0.55)
    cols = (0.008, 0.056, 0.135, 0.238, 0.320, 0.430, 0.560, 0.700, 0.830)
    headers = ("Rank", "Candidate", "Interval\n(aa)", "Length\n(aa)",
               "Affected\nisoforms", "Supporting\ncomparisons", "Evidence\nscore",
               "Evidence\nstrength", "Biological\nvalidation")
    for x, head in zip(cols, headers):
        ax.text(x, n - 0.22, head, ha="left", va="bottom", fontsize=ps.FONT["small"],
                color=ps.MUTED, fontweight="bold", linespacing=1.35)
    ax.plot([0, 1], [n - 0.32, n - 0.32], color=AXIS_GREY, lw=0.8, zorder=1)

    total_models = n_models or max((r["n_affected_isoforms"] for r in rows), default=0)
    for i, r in enumerate(rows):
        y = n - 1 - i
        if i % 2 == 0:
            ax.axhspan(y - 0.42, y + 0.42, color="#F6F8FA", zorder=0)
        ax.add_patch(plt.Rectangle((cols[0] - 0.006, y - 0.42), 0.004, 0.84,
                                   facecolor=CANDIDATE_COLOR, edgecolor="none",
                                   zorder=2))
        ax.text(cols[0], y, str(r["rank"]), ha="left", va="center",
                fontsize=ps.FONT["small"], color=ps.INK, zorder=3)
        ax.text(cols[1], y, r["candidate_label"], ha="left", va="center",
                fontsize=ps.FONT["small"], color=ps.INK, fontweight="bold", zorder=3)
        ax.text(cols[2], y, f"{r['aa_start']}–{r['aa_end']}", ha="left", va="center",
                fontsize=ps.FONT["small"], color=ps.INK, zorder=3)
        ax.text(cols[3], y, str(r["length_aa"]), ha="left", va="center",
                fontsize=ps.FONT["small"], color=ps.INK, zorder=3)
        affected = (f"{r['n_affected_isoforms']} of {total_models}" if total_models
                    else str(r["n_affected_isoforms"]))
        ax.text(cols[4], y, affected, ha="left", va="center",
                fontsize=ps.FONT["small"], color=ps.INK, zorder=3)
        ax.text(cols[5], y, str(r["supporting_comparisons"]), ha="left", va="center",
                fontsize=ps.FONT["small"], color=ps.INK, zorder=3)
        ax.add_patch(plt.Rectangle((cols[6], y - 0.17),
                                   0.105 * r["evidence_score"] / top, 0.34,
                                   facecolor=CANDIDATE_COLOR, alpha=0.35,
                                   edgecolor="none", zorder=2))
        ax.text(cols[6] + 0.004, y, str(r["evidence_score"]), ha="left", va="center",
                fontsize=ps.FONT["small"], color=ps.INK, zorder=3)
        ax.text(cols[7], y, r["evidence_strength"] or "—", ha="left", va="center",
                fontsize=ps.FONT["small"], color=ps.INK, zorder=3)
        ax.text(cols[8], y, r["biological_validation"], ha="left", va="center",
                fontsize=ps.FONT["small"], color=CANDIDATE_EDGE, zorder=3)

    ax.axis("off")
    species = f"{species_name} · " if species_name else ""
    reference = f"{protein_id} · " if protein_id else ""
    figure_title(
        ax, f"{gene_symbol} — exploratory candidate ranking",
        subtitle=f"{species}{reference}{n} exploratory candidate region(s) · "
                 f"amino-acid intervals on the primary protein · ranked by evidence "
                 f"score",
        note=f"Which isoform-difference regions carry the most transparent support? "
             f"{EXPLORATORY_TAG} regions only — {VALIDATION_TAG}. The evidence score "
             f"summarises isoform support, exon-boundary agreement, alignment support "
             f"and domain context; a score is a measure of support, not a validation.")
    handles = [ps.legend_patch(CANDIDATE_COLOR,
                               f"{EXPLORATORY_TAG.lower()} region · {VALIDATION_TAG}"),
               ps.legend_patch("#DCC58A", "evidence score relative to the highest-scoring "
                                          "candidate")]
    shared_legend(ax, handles, ncol=1, loc="upper left",
                  bbox=(0.0, _below_axes(ax, 6)))
    _footnote(ax, footnote or "", points=6 + 12 * len(handles) + 8)
    save_figure_all_formats(fig, fig_dir, stem)
    return True


def plot_domain_architecture(fig_dir: Path, stem: str, *, gene_symbol: str,
                             protein_id: str, protein_length: int,
                             domains: Sequence[Dict[str, Any]],
                             tm_regions: Sequence[Dict[str, Any]] = (),
                             exon_boundaries: Sequence[int] = (),
                             candidates: Sequence[Dict[str, Any]] = ()) -> bool:
    if not domains or not protein_length:
        return False
    fig, ax = plt.subplots(figsize=(12, 3.2))
    ax.plot([1, protein_length], [0, 0], color=ps.INK, lw=2, zorder=2)
    for i, domain in enumerate(domains):
        s, e = _int(domain.get("start_aa")), _int(domain.get("end_aa"))
        color = ps.PALETTE[i % len(ps.PALETTE)]
        ax.add_patch(plt.Rectangle((s, -0.22), max(1, e - s), 0.44,
                                   facecolor=color, edgecolor=ps.INK, lw=0.5, zorder=4))
        ax.text((s + e) / 2, 0.3, domain.get("domain_name") or domain.get("domain_id", ""),
                ha="center", va="bottom", fontsize=ps.FONT["small"], rotation=20)
    for tm in tm_regions:
        s, e = _int(tm.get("start_aa")), _int(tm.get("end_aa"))
        ax.add_patch(plt.Rectangle((s, -0.38), max(1, e - s), 0.76,
                                   facecolor="#E69F00", alpha=0.75, zorder=5))
    for pos in exon_boundaries:
        ax.axvline(pos, ymin=0.32, ymax=0.68, color=ps.MUTED, lw=0.5)
    for cand in candidates:
        ax.axvspan(_int(cand.get("aa_start")), _int(cand.get("aa_end")),
                   color=CANDIDATE_COLOR, alpha=0.15, lw=0)
    ax.set_xlim(0, protein_length * 1.02)
    ax.set_ylim(-0.9, 1.0)
    ax.set_yticks([])
    ax.set_xlabel("Protein position (aa)", fontsize=ps.FONT["label"])
    ax.spines[["top", "right", "left"]].set_visible(False)
    figure_title(ax, f"{gene_symbol} — domain architecture",
                 subtitle=f"{protein_id}; real InterProScan domains and pyTMHMM regions")
    shared_legend(ax, [ps.legend_patch(ps.PALETTE[0], "InterPro signature"),
                       ps.legend_patch("#E69F00", "pyTMHMM TM"),
                       ps.legend_patch(CANDIDATE_COLOR, "exploratory candidate")],
                  ncol=3, bbox=(0.5, -0.24))
    save_figure_all_formats(fig, fig_dir, stem)
    return True


def _span(f):
    return _int(f.get("start") or f.get("start_aa")), _int(f.get("end") or f.get("end_aa"))


_STACK_LANE_LABELS = {
    "domains": "Representative domains",
    "families": "Families / superfamilies",
    "signatures": "Member signatures",
    "tms": "Transmembrane topology",
    "exons": "Coding exons",
    "boundaries": "Coding-exon boundaries",
    "candidates": "Exploratory candidates",
}
FAMILY_COLOR = "#9AA7B6"
SIGNATURE_COLOR = "#6F7F92"
TM_COLOR = "#9467BD"
BOUNDARY_TICK_COLOR = "#44546A"
# Domain palette deliberately excludes the amber/orange hues reserved for
# exploratory candidate overlays, so a domain can never be read as a candidate.
DOMAIN_PALETTE = ("#0072B2", "#009E73", "#CC79A7", "#56B4E9", "#117733",
                  "#882255", "#44AA99", "#332288", "#6699CC", "#771155")


def plot_domain_feature_stack(fig_dir: Path, stem: str, *, gene_symbol: str,
                              species_name: str, protein_id: str, protein_length: int,
                              domains: Sequence[Dict[str, Any]] = (),
                              families: Sequence[Dict[str, Any]] = (),
                              signatures: Sequence[Dict[str, Any]] = (),
                              tms: Sequence[Dict[str, Any]] = (),
                              exons: Sequence[Dict[str, Any]] = (),
                              boundaries: Sequence[int] = (),
                              candidates: Sequence[Dict[str, Any]] = (),
                              lanes: Sequence[str] = ("domains", "families", "tms",
                                                      "exons", "boundaries", "candidates"),
                              subtitle: Optional[str] = None,
                              note: Optional[str] = None,
                              footnote: Optional[str] = None) -> bool:
    if not protein_length:
        return False
    order = [ln for ln in ("domains", "families", "signatures", "tms",
                           "exons", "boundaries", "candidates") if ln in lanes]
    data = {"domains": domains, "families": families, "signatures": signatures,
            "tms": tms, "exons": exons, "boundaries": boundaries,
            "candidates": candidates}
    # TM stays when domains are shown so the honest "no TM region predicted"
    # statement is made rather than silently dropping the layer.
    order = [ln for ln in order if data.get(ln) or (ln == "tms" and domains)]
    if not order:
        return False

    numbered = number_domain_instances(domains) if domains else []
    lane_items: Dict[str, List[Tuple[float, float]]] = {
        "domains": [(d["start"], d["end"]) for d in numbered],
        "families": [_span(f) for f in families],
        "signatures": [_span(f) for f in signatures],
        "tms": [_span(f) for f in tms],
        "exons": [_span(b) for b in exons],
        "boundaries": [],
        "candidates": [_span(c) for c in candidates],
    }
    packing = {ln: (_pack_rows(lane_items[ln]) if lane_items[ln] else [])
               for ln in order}
    n_sub = {ln: (max(packing[ln]) + 1 if packing[ln] else 1) for ln in order}

    block_h, sub_gap, lane_pad = 0.30, 0.42, 0.62
    lane_height = {ln: (n_sub[ln] - 1) * sub_gap + block_h + lane_pad for ln in order}
    total = sum(lane_height.values())
    fig, ax = plt.subplots(figsize=(12.5, max(3.4, 0.80 * total + 2.9)))

    lane_top: Dict[str, float] = {}
    cursor = total
    for ln in order:
        lane_top[ln] = cursor
        cursor -= lane_height[ln]

    def sub_y(ln: str, k: int) -> float:
        centre = lane_top[ln] - lane_height[ln] / 2.0
        return centre + ((n_sub[ln] - 1) / 2.0 - k) * sub_gap

    h = block_h / 2.0
    left = -protein_length * 0.235
    ax.set_xlim(left, protein_length * 1.03)
    ax.set_ylim(0.0, total)

    handles: List[Any] = []
    for i, ln in enumerate(order):
        rows = packing[ln] or [0] * max(1, len(lane_items[ln]))
        lane_centre = lane_top[ln] - lane_height[ln] / 2.0
        if i % 2 == 0:
            ax.axhspan(lane_top[ln] - lane_height[ln], lane_top[ln],
                       color="#F6F8FA", zorder=0)
        ax.text(left + protein_length * 0.012, lane_centre, _STACK_LANE_LABELS[ln],
                ha="left", va="center", fontsize=ps.FONT["label"], color=ps.INK)

        if ln == "domains":
            for d, k in zip(numbered, rows):
                y = sub_y(ln, k)
                ax.add_patch(plt.Rectangle((d["start"], y - h), max(1, d["end"] - d["start"]),
                                           2 * h, facecolor=d["color"], edgecolor=ps.INK,
                                           lw=0.5, zorder=4))
            _place_block_labels(ax, [(d["start"], d["end"], d["short"]) for d in numbered],
                                y=sub_y(ln, 0), height=h, fontsize=ps.FONT["small"],
                                color="white")
            seen = set()
            for d in numbered:
                if d["legend"] in seen:
                    continue
                seen.add(d["legend"])
                handles.append(ps.legend_patch(d["color"], d["legend"]))
        elif ln == "families":
            for f, k in zip(families, rows):
                s, e = _span(f)
                y = sub_y(ln, k)
                ax.add_patch(plt.Rectangle((s, y - h * 0.8), max(1, e - s), 1.6 * h,
                                           facecolor=FAMILY_COLOR, edgecolor=ps.INK,
                                           lw=0.4, alpha=0.85, zorder=4))
                if (e - s) > protein_length * 0.12:
                    ax.text((s + e) / 2.0, y, f.get("label") or "", ha="center",
                            va="center", fontsize=ps.FONT["small"], color=ps.INK, zorder=8)
            handles.append(ps.legend_patch(FAMILY_COLOR, "family / superfamily (InterPro)"))
        elif ln == "signatures":
            for f, k in zip(signatures, rows):
                s, e = _span(f)
                y = sub_y(ln, k)
                ax.add_patch(plt.Rectangle((s, y - h * 0.6), max(1, e - s), 1.2 * h,
                                           facecolor=SIGNATURE_COLOR, edgecolor="none",
                                           alpha=0.8, zorder=4))
            handles.append(ps.legend_patch(SIGNATURE_COLOR, "member-database signature"))
        elif ln == "tms":
            if tms:
                for f, k in zip(tms, rows):
                    s, e = _span(f)
                    y = sub_y(ln, k)
                    ax.add_patch(plt.Rectangle((s, y - h), max(1, e - s), 2 * h,
                                               facecolor=TM_COLOR, edgecolor=ps.INK,
                                               lw=0.4, zorder=4))
                _place_block_labels(
                    ax, [(_span(f)[0], _span(f)[1], f.get("label") or "TM") for f in tms],
                    y=sub_y(ln, 0), height=h, fontsize=ps.FONT["small"])
                handles.append(ps.legend_patch(TM_COLOR, "predicted transmembrane helix (pyTMHMM)"))
            else:
                ax.text(protein_length * 0.5, lane_centre,
                        "No transmembrane region predicted by pyTMHMM",
                        ha="center", va="center", fontsize=ps.FONT["small"],
                        color=ps.MUTED, style="italic")
        elif ln == "exons":
            for j, (b, k) in enumerate(zip(exons, rows)):
                s, e = _span(b)
                y = sub_y(ln, k)
                ax.add_patch(plt.Rectangle((s, y - h), max(1, e - s), 2 * h,
                                           facecolor=EXON_BLOCK_COLORS[j % 2],
                                           edgecolor=ps.INK, lw=0.3, zorder=4))
            _place_block_labels(
                ax, [(_span(b)[0], _span(b)[1],
                      b.get("label") or f"E{b.get('exon_number', '')}") for b in exons],
                y=sub_y(ln, 0), height=h, fontsize=ps.FONT["small"])
            handles.append(ps.legend_patch(EXON_BLOCK_COLORS[1],
                                           "CDS-derived coding exon"))
        elif ln == "boundaries":
            y = lane_centre
            ax.plot([1, protein_length], [y, y], color=AXIS_GREY, lw=0.7, zorder=1)
            for pos in boundaries:
                ax.plot([pos, pos], [y - h, y + h], color=BOUNDARY_TICK_COLOR,
                        lw=1.1, zorder=5)
            handles.append(ps.legend_line(BOUNDARY_TICK_COLOR,
                                          "internal coding-exon boundary", lw=1.1))
        elif ln == "candidates":
            for c, k in zip(candidates, rows):
                s, e = _span(c)
                y = sub_y(ln, k)
                ax.add_patch(plt.Rectangle((s, y - h), max(1, e - s), 2 * h,
                                           facecolor=CANDIDATE_COLOR, alpha=0.85,
                                           edgecolor=CANDIDATE_EDGE, lw=0.5, zorder=4))
            _place_block_labels(
                ax, [(_span(c)[0], _span(c)[1], c.get("id") or "") for c in candidates],
                y=sub_y(ln, 0), height=h, fontsize=ps.FONT["small"])
            handles.append(ps.legend_patch(CANDIDATE_COLOR,
                                           "exploratory candidate region (not validated)"))

    ax.set_yticks([])
    ax.spines[["top", "right", "left"]].set_visible(False)
    _readable_aa_axis(ax, protein_length,
                      f"Amino-acid position on {protein_id} (aa)")
    figure_title(ax, f"{gene_symbol} — integrated protein feature architecture",
                 subtitle=subtitle or f"{species_name} · {protein_id} · {protein_length} aa",
                 note=note)
    ncol = 3 if len(handles) > 4 else max(1, len(handles))
    n_legend_rows = 1 + (len(handles) - 1) // max(ncol, 1)
    shared_legend(ax, handles, ncol=ncol, loc="upper center",
                  bbox=(0.5, _below_axes(ax, 46)))
    _footnote(ax, footnote or "", points=46 + 12 * n_legend_rows + 6)
    save_figure_all_formats(fig, fig_dir, stem)
    return True


MEMBER_SIGNATURE_TABLE_COLUMNS = (
    "member_database", "signature_accession", "signature_name", "start_aa", "end_aa",
    "length_aa", "interpro_accession", "interpro_entry_name", "integration_status",
)
INTEGRATED_SIGNATURE_COLOR = "#6F7F92"
UNINTEGRATED_SIGNATURE_COLOR = "#FFFFFF"


def member_signature_groups(signatures: Sequence[Dict[str, Any]]
                            ) -> List[Tuple[str, List[Dict[str, Any]]]]:
    groups: Dict[str, List[Dict[str, Any]]] = {}
    for s in signatures:
        start, end = _span(s)
        if start is None or end is None:
            continue
        tip = s.get("tooltip") or {}
        database = str(s.get("source") or s.get("member_database")
                       or tip.get("member_database") or "unspecified").replace("_", " ")
        interpro = s.get("interpro_accession") or tip.get("interpro_accession") or ""
        integrated = bool(tip["is_integrated"]) if "is_integrated" in tip \
            else bool(interpro)
        groups.setdefault(database, []).append({
            "database": database,
            "accession": s.get("signature_accession") or s.get("accession") or "",
            "name": tip.get("signature_name") or s.get("signature_name")
                    or s.get("label") or "",
            "entry_name": (s.get("label") or tip.get("interpro_name") or "")
                          if integrated else "",
            "start": start, "end": end,
            "interpro_accession": interpro,
            "integrated": integrated,
        })
    for members in groups.values():
        members.sort(key=lambda r: (r["start"], r["end"], r["accession"]))
    # Databases with the most signatures first, so the reader meets the best-covered
    # evidence first; ties keep alphabetical order for a reproducible layout.
    return sorted(groups.items(), key=lambda kv: (-len(kv[1]), kv[0]))


def member_signature_table(signatures: Sequence[Dict[str, Any]]) -> List[Dict[str, Any]]:
    return [{
        "member_database": s["database"], "signature_accession": s["accession"],
        "signature_name": s["name"], "start_aa": s["start"], "end_aa": s["end"],
        "length_aa": s["end"] - s["start"] + 1,
        "interpro_accession": s["interpro_accession"] or "",
        "interpro_entry_name": s["entry_name"],
        "integration_status": "integrated" if s["integrated"] else "not integrated",
    } for _name, members in member_signature_groups(signatures) for s in members]


def plot_member_signature_supplement(fig_dir: Path, stem: str, *, gene_symbol: str,
                                     species_name: str, protein_id: str,
                                     protein_length: int,
                                     signatures: Sequence[Dict[str, Any]],
                                     footnote: Optional[str] = None) -> bool:
    groups = member_signature_groups(signatures)
    if not groups or not protein_length:
        return False
    n_sig = sum(len(members) for _n, members in groups)
    n_rows = n_sig + len(groups)
    n_integrated = sum(1 for _n, members in groups for s in members if s["integrated"])

    fig, ax = plt.subplots(figsize=(13.2, max(3.4, 0.30 * n_rows + 3.0)))
    # Three text columns share the figure with the amino-acid track: a name, an
    # accession that must never be shortened, and the InterPro entry it feeds.
    left = -protein_length * 0.58
    right = protein_length * 1.34
    ax.set_xlim(left, right)
    ax.set_ylim(-0.9, n_rows + 0.55)
    group_x = left + protein_length * 0.006
    name_x = left + protein_length * 0.030
    accession_x = left + protein_length * 0.280
    entry_x = protein_length * 1.02
    # Each text column is shortened to what its own column really holds, so a long
    # signature name cannot run into the accession beside it.
    name_chars = _margin_chars(ax, name_x, ps.FONT["small"], accession_x)
    entry_chars = _margin_chars(ax, entry_x, ps.FONT["small"], right)

    y = n_rows - 1.0
    for gi, (database, members) in enumerate(groups):
        if gi % 2 == 0:
            ax.axhspan(y - len(members) + 0.5, y + 0.5, color="#F6F8FA", zorder=0)
        ax.text(group_x, y, f"{database}  ({len(members)} signature(s))", ha="left",
                va="center", fontsize=ps.FONT["label"], color=ps.INK,
                fontweight="bold")
        ax.plot([group_x, protein_length], [y + 0.5, y + 0.5], color=AXIS_GREY,
                lw=0.5, zorder=1)
        y -= 1.0
        for s in members:
            ax.add_patch(plt.Rectangle(
                (s["start"], y - 0.24), max(1, s["end"] - s["start"]), 0.48,
                facecolor=INTEGRATED_SIGNATURE_COLOR if s["integrated"]
                else UNINTEGRATED_SIGNATURE_COLOR,
                edgecolor=INTEGRATED_SIGNATURE_COLOR,
                lw=0.5 if s["integrated"] else 0.9,
                linestyle="-" if s["integrated"] else (0, (2.5, 1.5)), zorder=4))
            ax.text(name_x, y, _clip(s["name"] or s["accession"], name_chars),
                    ha="left", va="center", fontsize=ps.FONT["small"], color=ps.INK)
            ax.text(accession_x, y, s["accession"], ha="left", va="center",
                    fontsize=ps.FONT["small"], color=ps.MUTED)
            entry = (f"{s['interpro_accession']} {s['entry_name']}".strip()
                     if s["integrated"] else "not integrated")
            ax.text(entry_x, y, _clip(entry, entry_chars), ha="left", va="center",
                    fontsize=ps.FONT["small"],
                    color=ps.MUTED if s["integrated"] else "#8A5008")
            y -= 1.0
    for x, head in ((name_x, "signature"), (accession_x, "accession"),
                    (entry_x, "InterPro entry")):
        ax.text(x, n_rows - 0.22, head, ha="left", va="bottom",
                fontsize=ps.FONT["small"], color=ps.MUTED, fontweight="bold")

    ax.set_yticks([])
    ax.spines[["top", "right", "left"]].set_visible(False)
    _readable_aa_axis(ax, protein_length, f"Amino-acid position on {protein_id} (aa)")
    figure_title(
        ax, f"{gene_symbol} — member-database signature supplement",
        subtitle=f"{species_name} · {protein_id} · {n_sig} signature(s) from "
                 f"{len(groups)} member database(s) · {n_integrated} integrated into "
                 f"an InterPro entry",
        note="Which member-database signatures underlie the representative domain "
             "annotation? Supplementary evidence: individual signatures overlap and "
             "are not independent observations, and an unintegrated signature is not "
             "an annotated domain.")
    handles = [ps.legend_patch(INTEGRATED_SIGNATURE_COLOR,
                               "signature integrated into an InterPro entry"),
               plt.Rectangle((0, 0), 1, 1, facecolor=UNINTEGRATED_SIGNATURE_COLOR,
                             edgecolor=INTEGRATED_SIGNATURE_COLOR, lw=0.9,
                             linestyle=(0, (2.5, 1.5)),
                             label="signature not integrated into an InterPro entry")]
    shared_legend(ax, handles, ncol=2, loc="upper left",
                  bbox=(0.0, _below_axes(ax, 32)))
    _footnote(ax, footnote or "", points=32 + 12 + 8)
    save_figure_all_formats(fig, fig_dir, stem)
    return True


def plot_exon_domain_boundary_distribution(fig_dir: Path, stem: str, *,
                                            gene_symbol: str,
                                            boundaries: Sequence[Dict[str, Any]]) -> bool:
    rows = [r for r in boundaries if r.get("absolute_distance_aa") not in (None, "")]
    if not rows:
        return False
    vals = [_int(r.get("absolute_distance_aa")) for r in rows]
    fig, ax = plt.subplots(figsize=(8.5, 4.3))
    ax.hist(vals, bins=min(15, max(5, len(set(vals)))), color=ps.PALETTE[0],
            edgecolor=ps.INK, linewidth=0.5)
    ax.axvline(5, color=CANDIDATE_COLOR, ls="--", lw=1, label="near-boundary threshold (5 aa)")
    ax.set_xlabel("Absolute distance to nearest domain edge (aa)", fontsize=ps.FONT["label"])
    ax.set_ylabel("Exon boundaries", fontsize=ps.FONT["label"])
    ax.spines[["top", "right"]].set_visible(False)
    figure_title(ax, f"{gene_symbol} — exon–domain boundary distances",
                 subtitle="Generic all-exon analysis on real InterProScan coordinates")
    shared_legend(ax, [ps.legend_patch(ps.PALETTE[0], "boundary count"),
                       ps.legend_patch(CANDIDATE_COLOR, "≤5 aa threshold")],
                  ncol=2, bbox=(0.5, -0.22))
    save_figure_all_formats(fig, fig_dir, stem)
    return True


# canonical boundary-class colours (color-blind-safe; align with the frontend)
BOUNDARY_CLASS_COLOR = {
    "exact_domain_edge": "#1B7837",
    "near_domain_edge": "#A6DBA0",
    "inside_domain": "#FDB863",
    "outside_annotated_domains": "#B2ABD2",
    "unavailable_or_uncertain": "#D9D9D9",
    # legacy aliases (generic core category names)
    "exact_edge": "#1B7837",
    "near_edge": "#A6DBA0",
    "outside_domain": "#B2ABD2",
    "unknown": "#D9D9D9",
}
BOUNDARY_CLASS_LABEL = {
    "exact_domain_edge": "exact domain edge (0 aa)",
    "near_domain_edge": "near domain edge (≤ threshold)",
    "inside_domain": "inside domain",
    "outside_annotated_domains": "outside annotated domains",
    "unavailable_or_uncertain": "unavailable / uncertain",
}


def _canonical_boundary_class(cls: str) -> str:
    alias = {"exact_edge": "exact_domain_edge", "near_edge": "near_domain_edge",
             "outside_domain": "outside_annotated_domains", "unknown": "unavailable_or_uncertain"}
    return alias.get(cls or "", cls or "unavailable_or_uncertain")


def _boundary_signed(b) -> int:
    v = b.get("signed_distance_aa")
    return _int(v) if v is not None else _int(b.get("signed_distance"))


def _boundary_label(b) -> str:
    raw = (b.get("label") or b.get("exon_boundary_id") or b.get("id") or "")
    if ":" in raw and "→" not in raw and "->" not in raw:
        raw = raw.split(":")[-1]
    return safe_text(raw)




def _resolve_nearest_edge(b: Dict[str, Any],
                          numbered: Sequence[Dict[str, Any]]) -> Tuple[Optional[int], str, str]:
    pos = _int(b.get("boundary_position_aa") or b.get("protein_position") or b.get("start"))
    edge_type = str(b.get("nearest_edge") or b.get("nearest_edge_type") or "").lower()
    signed = b.get("signed_distance_aa")
    if signed is None:
        signed = b.get("signed_distance")
    edge_pos = pos - _int(signed) if signed is not None else b.get("nearest_edge_position")
    edge_pos = _int(edge_pos) if edge_pos is not None else None
    name = ""
    if edge_pos is not None:
        for d in numbered:
            if d["start"] == edge_pos:
                name, edge_type = d["name"], edge_type or "start"
                break
            if d["end"] == edge_pos:
                name, edge_type = d["name"], edge_type or "end"
                break
    if not name:
        key = (_int(b.get("nearest_domain_start")), _int(b.get("nearest_domain_end")))
        name = {(d["start"], d["end"]): d["name"] for d in numbered}.get(key, "")
    if not name:
        name = _pretty_domain_name(b.get("nearest_domain_label")
                                   or b.get("nearest_domain_name") or "")
    return edge_pos, edge_type, name or "—"


def boundary_definition_footnote(threshold: int, protein_id: str,
                                 domain_source: str = "representative InterPro domains "
                                                      "(InterProScan)") -> str:
    return (
        "Distance definition: signed amino-acid distance from the internal coding-exon "
        "boundary to the nearest edge of the nearest representative domain; negative = "
        f"upstream of that edge, positive = downstream, 0 = exactly on the edge. "
        f"Near-edge threshold: ±{threshold} aa. Representative-domain source: {domain_source}. "
        f"Protein coordinate system: 1-based amino-acid positions on {protein_id}."
    )


def plot_signed_boundary_distances(fig_dir: Path, stem: str, *, gene_symbol: str,
                                   species_name: str, protein_id: str,
                                   boundaries: Sequence[Dict[str, Any]],
                                   threshold: int = 5,
                                   domains: Sequence[Dict[str, Any]] = (),
                                   footnote: Optional[str] = None) -> bool:
    rows = [b for b in boundaries if b.get("signed_distance_aa") is not None
            or b.get("signed_distance") is not None]
    if not rows:
        return False
    rows = sorted(rows, key=lambda b: _int(b.get("boundary_position_aa")
                                           or b.get("protein_position") or b.get("start") or 0))
    numbered = number_domain_instances(domains) if domains else []
    n = len(rows)
    extreme = max((abs(_boundary_signed(b)) for b in rows), default=1) or 1
    limit = extreme * 1.18 + threshold

    fig, ax = plt.subplots(figsize=(11, max(3.0, 0.42 * n + 2.8)))
    ax.set_xlim(-limit, limit)
    ax.set_ylim(-0.9, n + 0.25)
    ax.axvspan(-threshold, threshold, color=BOUNDARY_CLASS_COLOR["near_domain_edge"],
               alpha=0.30, lw=0, zorder=1)
    ax.axvline(0, color=ps.INK, lw=1.3, zorder=3)
    ax.text(0, n - 0.42, "domain edge (0 aa)", ha="center", va="bottom",
            fontsize=ps.FONT["small"], color=ps.MUTED)
    ax.text(threshold, n - 0.88, f"±{threshold} aa near-edge interval", ha="left",
            va="bottom", fontsize=ps.FONT["small"], color=ps.MUTED)

    seen, edges_seen = [], set()
    labels = []
    for i, b in enumerate(rows):
        y = n - 1 - i
        sd = _boundary_signed(b)
        cls = _canonical_boundary_class(b.get("class") or b.get("category"))
        if cls not in seen:
            seen.append(cls)
        col = BOUNDARY_CLASS_COLOR.get(cls, ps.MUTED)
        _, edge, dom = _resolve_nearest_edge(b, numbered)
        marker = ">" if edge == "start" else ("<" if edge == "end" else "o")
        edges_seen.add(edge)
        ax.plot([0, sd], [y, y], color=col, lw=1.4, zorder=4, solid_capstyle="round")
        ax.scatter([sd], [y], s=85, marker=marker, facecolor=col, edgecolor=ps.INK,
                   linewidth=0.6, zorder=5)
        labels.append(_boundary_label(b))
        ax.text(limit * 0.99, y, f"{sd:+d} aa · {dom} ({edge or 'edge'})", ha="right",
                va="center", fontsize=ps.FONT["small"], color=ps.MUTED)

    ax.set_yticks(range(n))
    ax.set_yticklabels(list(reversed(labels)), fontsize=ps.FONT["small"])
    ax.tick_params(axis="y", length=0)
    ax.set_xlabel("Signed distance to the nearest representative-domain edge (aa)",
                  fontsize=ps.FONT["label"])
    ax.tick_params(axis="x", labelsize=ps.FONT["tick"])
    ax.spines[["top", "right", "left"]].set_visible(False)
    figure_title(
        ax, f"{gene_symbol} — signed exon-boundary distances",
        subtitle=f"{species_name} · {protein_id} · {n} internal coding-exon boundaries · "
                 f"0 = representative-domain edge · shaded = ±{threshold} aa near-edge interval",
        note="Adjacent-exon transitions are labelled on the y axis; negative values lie "
             "upstream of the nearest domain edge, positive values downstream.")
    handles = [ps.legend_patch(BOUNDARY_CLASS_COLOR[c], BOUNDARY_CLASS_LABEL[c])
               for c in ("exact_domain_edge", "near_domain_edge", "inside_domain",
                         "outside_annotated_domains", "unavailable_or_uncertain")]
    if "start" in edges_seen:
        handles.append(_legend_marker(">", "nearest edge = domain start"))
    if "end" in edges_seen:
        handles.append(_legend_marker("<", "nearest edge = domain end"))
    shared_legend(ax, handles, ncol=3, loc="upper center", bbox=(0.5, _below_axes(ax, 40)))
    _footnote(ax, footnote or boundary_definition_footnote(threshold, protein_id),
              points=40 + 12 * (1 + (len(handles) - 1) // 3) + 6)
    save_figure_all_formats(fig, fig_dir, stem)
    return True


def plot_boundary_on_architecture(fig_dir: Path, stem: str, *, gene_symbol: str,
                                  species_name: str, protein_id: str, protein_length: int,
                                  domains: Sequence[Dict[str, Any]],
                                  exon_boundaries: Sequence[Dict[str, Any]],
                                  exon_blocks: Sequence[Dict[str, Any]] = (),
                                  threshold: int = 5,
                                  footnote: Optional[str] = None) -> bool:
    if not protein_length or not exon_boundaries:
        return False
    numbered = number_domain_instances(domains) if domains else []
    fig, ax = plt.subplots(figsize=(12.5, 3.9))
    left = -protein_length * 0.19
    ax.set_xlim(left, protein_length * 1.03)
    ax.set_ylim(-0.34, 1.52)

    y_exon, y_dom = 0.0, 0.72
    ax.plot([1, protein_length], [y_exon, y_exon], color=AXIS_GREY, lw=0.8, zorder=1)
    for i, b in enumerate(exon_blocks):
        s = _int(b.get("start", b.get("protein_start_aa")))
        e = _int(b.get("end", b.get("protein_end_aa")))
        ax.add_patch(plt.Rectangle((s, y_exon - 0.13), max(1, e - s), 0.26,
                                   facecolor=EXON_BLOCK_COLORS[i % 2], edgecolor=ps.INK,
                                   linewidth=0.45, zorder=3))
    ax.text(left, y_exon, "Coding exons", ha="left", va="center",
            fontsize=ps.FONT["label"], color=ps.INK)

    for d in numbered:
        ax.add_patch(plt.Rectangle((d["start"], y_dom - 0.17),
                                   max(1, d["end"] - d["start"]), 0.34,
                                   facecolor=d["color"], edgecolor=ps.INK, lw=0.5,
                                   alpha=0.95, zorder=3))
    _place_block_labels(ax, [(d["start"], d["end"], str(d["instance"] or ""))
                             for d in numbered],
                        y=y_dom, height=0.17, fontsize=ps.FONT["small"], color="white")
    ax.text(left, y_dom, "Representative domains", ha="left", va="center",
            fontsize=ps.FONT["label"], color=ps.INK)

    seen = []
    for b in exon_boundaries:
        pos = _int(b.get("boundary_position_aa") or b.get("protein_position") or b.get("start"))
        cls = _canonical_boundary_class(b.get("class") or b.get("category"))
        if cls not in seen:
            seen.append(cls)
        col = BOUNDARY_CLASS_COLOR.get(cls, ps.MUTED)
        edge_pos, _, _ = _resolve_nearest_edge(b, numbered)
        ax.plot([pos, pos], [y_exon - 0.16, 1.18], color=col, lw=1.5, zorder=5)
        ax.scatter([pos], [1.18], s=34, color=col, edgecolor=ps.INK, linewidth=0.5, zorder=6)
        if edge_pos is not None and _int(edge_pos) != pos:
            ax.plot([pos, _int(edge_pos)], [1.18, 1.18], color=col, lw=1.0, zorder=4)
            ax.scatter([_int(edge_pos)], [1.18], s=10, color=ps.INK, zorder=6)
    ax.text(left, 1.18, "Boundary vs nearest edge", ha="left", va="center",
            fontsize=ps.FONT["label"], color=ps.INK)

    ax.set_yticks([])
    ax.spines[["top", "right", "left"]].set_visible(False)
    _readable_aa_axis(ax, protein_length, f"Amino-acid position on {protein_id} (aa)")
    figure_title(
        ax, f"{gene_symbol} — coding-exon boundaries on the domain architecture",
        subtitle=f"{species_name} · {protein_id} · {len(numbered)} representative domain "
                 f"instance(s) · {len(exon_boundaries)} internal coding-exon boundaries · "
                 f"near-edge ±{threshold} aa",
        note="Connectors run from each boundary to the nearest representative-domain edge.")
    handles = [ps.legend_patch(d["color"], d["legend"]) for d in numbered[:4]]
    handles.append(ps.legend_patch(EXON_BLOCK_COLORS[1], "coding exon"))
    handles += [ps.legend_patch(BOUNDARY_CLASS_COLOR[c], BOUNDARY_CLASS_LABEL[c])
                for c in ("exact_domain_edge", "near_domain_edge", "inside_domain",
                          "outside_annotated_domains", "unavailable_or_uncertain")]
    shared_legend(ax, handles, ncol=3, loc="upper center", bbox=(0.5, _below_axes(ax, 42)))
    _footnote(ax, footnote or boundary_definition_footnote(threshold, protein_id),
              points=42 + 12 * (1 + (len(handles) - 1) // 3) + 6)
    save_figure_all_formats(fig, fig_dir, stem)
    return True


def plot_boundary_class_summary(fig_dir: Path, stem: str, *, gene_symbol: str,
                                species_name: str, protein_id: str,
                                boundaries: Sequence[Dict[str, Any]],
                                threshold: int = 5,
                                footnote: Optional[str] = None) -> bool:
    rows = list(boundaries)
    if not rows:
        return False
    order = ("exact_domain_edge", "near_domain_edge", "inside_domain",
             "outside_annotated_domains", "unavailable_or_uncertain")
    counts = {c: 0 for c in order}
    for b in rows:
        counts[_canonical_boundary_class(b.get("class") or b.get("category"))] += 1
    total = max(1, len(rows))
    top = max(counts.values(), default=1) or 1

    fig, ax = plt.subplots(figsize=(9.5, 3.9))
    ax.set_xlim(0, top * 1.18)
    ax.set_ylim(-0.7, len(order) + 0.55)

    # one stacked strip showing the composition, then the per-class counts
    x = 0.0
    for c in order:
        if not counts[c]:
            continue
        width = top * 1.18 * counts[c] / total
        ax.add_patch(plt.Rectangle((x, len(order) - 0.30), width, 0.42,
                                   facecolor=BOUNDARY_CLASS_COLOR[c], edgecolor=ps.INK,
                                   lw=0.4, zorder=3))
        share = 100.0 * counts[c] / total
        if share >= 12:
            ax.text(x + width / 2.0, len(order) - 0.09, f"{share:.0f}%", ha="center",
                    va="center", fontsize=ps.FONT["small"], color=ps.INK, zorder=5)
        x += width
    ax.text(0, len(order) + 0.20, "Composition of all internal coding-exon boundaries",
            ha="left", va="bottom", fontsize=ps.FONT["small"], color=ps.MUTED)

    for i, c in enumerate(order):
        y = len(order) - 1 - i - 0.6
        ax.barh(y, counts[c], height=0.52, color=BOUNDARY_CLASS_COLOR[c],
                edgecolor=ps.INK, linewidth=0.5, zorder=3)
        ax.text(counts[c] + top * 0.02, y,
                f"{counts[c]} ({100.0 * counts[c] / total:.0f}%)",
                va="center", ha="left", fontsize=ps.FONT["small"], color=ps.MUTED)

    ax.set_yticks([len(order) - 1 - i - 0.6 for i in range(len(order))])
    ax.set_yticklabels([BOUNDARY_CLASS_LABEL[c] for c in order],
                       fontsize=ps.FONT["small"])
    ax.tick_params(axis="y", length=0)
    ax.set_xlabel("Internal coding-exon boundaries (count)", fontsize=ps.FONT["label"])
    ax.tick_params(axis="x", labelsize=ps.FONT["tick"])
    ax.xaxis.set_major_locator(MaxNLocator(integer=True, nbins=8))
    ax.spines[["top", "right", "left"]].set_visible(False)
    figure_title(
        ax, f"{gene_symbol} — boundary-class composition",
        subtitle=f"{species_name} · {protein_id} · {len(rows)} internal coding-exon "
                 f"boundaries · mutually exclusive classes · near-edge ±{threshold} aa")
    _footnote(ax, footnote or boundary_definition_footnote(threshold, protein_id),
              points=34)
    save_figure_all_formats(fig, fig_dir, stem)
    return True


def plot_selected_boundary_detail(fig_dir: Path, stem: str, *, gene_symbol: str,
                                  species_name: str, protein_id: str,
                                  boundary: Dict[str, Any],
                                  domains: Sequence[Dict[str, Any]] = (),
                                  exon_blocks: Sequence[Dict[str, Any]] = (),
                                  threshold: int = 5) -> bool:
    if not boundary:
        return False
    pos = _int(boundary.get("protein_position") or boundary.get("boundary_position_aa") or boundary.get("start"))
    absd = boundary.get("absolute_distance")
    if absd is None:
        absd = boundary.get("absolute_distance_aa")
    win = max(40, _int(absd) + 25)
    lo, hi = max(0, pos - win), pos + win
    cls = _canonical_boundary_class(boundary.get("class") or boundary.get("category"))
    col = BOUNDARY_CLASS_COLOR.get(cls, ps.MUTED)
    fig, ax = plt.subplots(figsize=(10, 3.0))
    for i, d in enumerate(domains):
        s, e = _int(d.get("start") or d.get("start_aa")), _int(d.get("end") or d.get("end_aa"))
        if e < lo or s > hi:
            continue
        ax.add_patch(plt.Rectangle((s, 0.12), max(1, e - s), 0.3,
                                   facecolor=ps.PALETTE[i % len(ps.PALETTE)], edgecolor=ps.INK,
                                   lw=0.5, alpha=0.9, zorder=3))
        ax.text((s + e) / 2, 0.46, d.get("label") or d.get("interpro_name") or "",
                ha="center", va="bottom", fontsize=ps.FONT["small"])
    for i, b in enumerate(exon_blocks):
        s, e = _int(b.get("start") or b.get("protein_start_aa")), _int(b.get("end") or b.get("protein_end_aa"))
        if e < lo or s > hi:
            continue
        ax.add_patch(plt.Rectangle((s, -0.34), max(1, e - s), 0.26,
                                   facecolor=EXON_BLOCK_COLORS[i % 2], edgecolor=ps.INK,
                                   lw=0.4, zorder=3))
        ax.text((s + e) / 2, -0.21, b.get("label") or str(b.get("exon_number", "")),
                ha="center", va="center", fontsize=ps.FONT["small"])
    edge_pos = _int(boundary.get("nearest_edge_position")) if boundary.get("nearest_edge_position") is not None else None
    if edge_pos is not None:
        ax.axvline(edge_pos, ymin=0.34, ymax=0.72, color=ps.INK, lw=1.0, ls="--", zorder=4)
        ax.plot([pos, edge_pos], [0.02, 0.02], color=col, lw=1.4, zorder=5)
    ax.plot([pos, pos], [-0.34, 0.42], color=col, lw=2.4, zorder=6)
    ax.set_xlim(lo, hi)
    ax.set_ylim(-0.5, 0.7)
    ax.set_yticks([])
    ax.set_xlabel("Protein position (aa)", fontsize=ps.FONT["label"])
    ax.spines[["top", "right", "left"]].set_visible(False)
    label = (boundary.get("label") or boundary.get("id") or f"boundary@{pos}").replace("→", "->")
    figure_title(ax, f"{gene_symbol} — selected boundary detail ({label})",
                 subtitle=f"{species_name} · {protein_id} · signed {boundary.get('signed_distance', '?')} aa · "
                          f"{BOUNDARY_CLASS_LABEL.get(cls, cls)}")
    handles = [ps.legend_patch(col, BOUNDARY_CLASS_LABEL.get(cls, cls)),
               ps.legend_patch(ps.PALETTE[0], "representative domain"),
               ps.legend_patch(EXON_BLOCK_COLORS[1], "coding exon")]
    shared_legend(ax, handles, ncol=3, bbox=(0.5, -0.3))
    save_figure_all_formats(fig, fig_dir, stem)
    return True


def plot_boundary_evidence_supplement(fig_dir: Path, stem: str, *, gene_symbol: str,
                                      species_name: str, protein_id: str,
                                      boundaries: Sequence[Dict[str, Any]],
                                      threshold: int = 5) -> bool:
    rows = list(boundaries)
    if not rows:
        return False
    rows = sorted(rows, key=lambda b: _int(b.get("protein_position") or b.get("boundary_position_aa") or 0))
    n = len(rows)
    fig, ax = plt.subplots(figsize=(11, max(2.8, 0.44 * n + 1.6)))
    for i, b in enumerate(rows):
        y = n - 1 - i
        cls = _canonical_boundary_class(b.get("class") or b.get("category"))
        col = BOUNDARY_CLASS_COLOR.get(cls, ps.MUTED)
        pos = _int(b.get("protein_position") or b.get("boundary_position_aa"))
        sd = b.get("signed_distance")
        if sd is None:
            sd = b.get("signed_distance_aa")
        ax.scatter([0.02], [y], s=90, color=col, edgecolor=ps.INK, linewidth=0.5, zorder=4)
        txt = (f"{(b.get('label') or b.get('id','')).replace('→', '->')} · aa {pos} · nearest "
               f"{b.get('nearest_domain_label') or b.get('nearest_domain_name') or '—'} "
               f"{b.get('nearest_edge_type') or ''} · signed {sd if sd is not None else '—'} aa · "
               f"{BOUNDARY_CLASS_LABEL.get(cls, cls)}")
        ax.text(0.06, y, txt, va="center", ha="left", fontsize=ps.FONT["small"], color=ps.INK)
    ax.set_xlim(0, 1)
    ax.set_ylim(-0.8, n - 0.2)
    ax.axis("off")
    figure_title(ax, f"{gene_symbol} — boundary evidence supplement",
                 subtitle=f"{species_name} · {protein_id} · {n} internal boundaries · near-edge ±{threshold} aa")
    handles = [ps.legend_patch(BOUNDARY_CLASS_COLOR[c], BOUNDARY_CLASS_LABEL[c])
               for c in ("exact_domain_edge", "near_domain_edge", "inside_domain",
                         "outside_annotated_domains", "unavailable_or_uncertain")]
    shared_legend(ax, handles, ncol=3, bbox=(0.5, -0.16))
    save_figure_all_formats(fig, fig_dir, stem)
    return True


CANDIDATE_CONTEXT_TABLE_COLUMNS = (
    "candidate_label", "candidate_id", "aa_start", "aa_end", "length_aa", "relation",
    "overlapping_domain_instance_ids", "overlapping_domain_labels",
    "nearest_domain_instance_id", "nearest_domain_label", "nearest_edge_type",
    "nearest_edge_position_aa", "distance_to_nearest_edge_aa",
    "coincident_exon_boundary", "biological_validation",
)


def resolve_candidate_domain_context(start: int, end: int,
                                     numbered_domains: Sequence[Dict[str, Any]],
                                     boundaries: Sequence[Dict[str, Any]] = ()
                                     ) -> Dict[str, Any]:
    overlapping = [d for d in numbered_domains
                   if min(end, d["end"]) >= max(start, d["start"])]
    inside = [d for d in overlapping if d["start"] <= start and end <= d["end"]]
    if inside:
        relation = f"inside {inside[0]['name']}"
    elif overlapping:
        relation = "overlapping " + ", ".join(d["name"] for d in overlapping)
    else:
        relation = "outside annotated domains"

    nearest = None
    for d in numbered_domains:
        for edge_type, pos in (("start", d["start"]), ("end", d["end"])):
            distance = 0 if start <= pos <= end else min(abs(start - pos), abs(end - pos))
            if nearest is None or distance < nearest["distance"]:
                nearest = {"domain": d, "edge_type": edge_type, "position": pos,
                           "distance": distance}
    coincident = ""
    for b in boundaries or ():
        pos = _int(b.get("boundary_position_aa") or b.get("protein_position")
                   or b.get("start"))
        if pos in (start, end):
            coincident = _boundary_label(b) or f"aa {pos}"
            break
    return {
        "relation": relation, "overlapping": overlapping, "inside": bool(inside),
        "nearest": nearest, "coincident_exon_boundary": coincident,
    }


def candidate_context_table(contexts: Sequence[Dict[str, Any]],
                            domains: Sequence[Dict[str, Any]] = (),
                            boundaries: Sequence[Dict[str, Any]] = ()
                            ) -> List[Dict[str, Any]]:
    numbered = number_domain_instances(domains) if domains else []
    out = []
    for i, c in enumerate(contexts, start=1):
        start = _int(c.get("aa_start", c.get("start_aa", c.get("start"))))
        end = _int(c.get("aa_end", c.get("end_aa", c.get("end"))))
        ctx = resolve_candidate_domain_context(start, end, numbered, boundaries)
        nearest = ctx["nearest"] or {}
        out.append({
            "candidate_label": c.get("candidate_label") or c.get("label") or f"C{i}",
            "candidate_id": c.get("candidate_id") or c.get("id") or "",
            "aa_start": start, "aa_end": end, "length_aa": end - start + 1,
            "relation": ctx["relation"],
            "overlapping_domain_instance_ids": ";".join(
                d["instance_id"] for d in ctx["overlapping"]),
            "overlapping_domain_labels": ";".join(d["name"] for d in ctx["overlapping"]),
            "nearest_domain_instance_id": nearest.get("domain", {}).get("instance_id", ""),
            "nearest_domain_label": nearest.get("domain", {}).get("name", ""),
            "nearest_edge_type": nearest.get("edge_type", ""),
            "nearest_edge_position_aa": nearest.get("position", ""),
            "distance_to_nearest_edge_aa": nearest.get("distance", ""),
            "coincident_exon_boundary": ctx["coincident_exon_boundary"],
            "biological_validation": "not validated",
        })
    return out


def plot_candidate_domain_context(fig_dir: Path, stem: str, *, gene_symbol: str,
                                  contexts: Sequence[Dict[str, Any]],
                                  domains: Sequence[Dict[str, Any]] = (),
                                  exon_blocks: Sequence[Dict[str, Any]] = (),
                                  boundaries: Sequence[Dict[str, Any]] = (),
                                  protein_id: str = "",
                                  protein_length: Optional[int] = None,
                                  species_name: str = "",
                                  footnote: Optional[str] = None) -> bool:
    table = candidate_context_table(contexts, domains, boundaries)
    if not table:
        return False
    numbered = number_domain_instances(domains) if domains else []
    max_aa = protein_length or max(
        [r["aa_end"] for r in table]
        + [_int(b.get("end", b.get("protein_end_aa"))) for b in exon_blocks]
        + [d["end"] for d in numbered] or [1])
    n = len(table)

    # Reference lanes stack directly above the candidate rows, one unit apart, so the
    # figure carries no dead band between the annotation and the candidates.
    y_exon = (n + 0.30) if exon_blocks else None
    y_dom = (y_exon + 1.05) if (exon_blocks and numbered) else \
        ((n + 0.30) if numbered else None)
    top = max(y for y in (y_exon, y_dom, float(n - 1)) if y is not None)
    lanes = int(bool(exon_blocks)) + int(bool(numbered))
    fig, ax = plt.subplots(figsize=(12.6, max(3.4, 0.62 * n + 0.78 * lanes + 2.3)))
    left = -max_aa * 0.40
    ax.set_xlim(left, max_aa * 1.05)
    ax.set_ylim(-0.9, top + 0.42)
    label_x = left + max_aa * 0.008
    h = 0.17
    handles: List[Any] = []

    if numbered:
        for d in numbered:
            ax.add_patch(plt.Rectangle((d["start"], y_dom - h),
                                       max(1, d["end"] - d["start"]), 2 * h,
                                       facecolor=d["color"], edgecolor=ps.INK, lw=0.5,
                                       zorder=4))
        _place_block_labels(ax, [(d["start"], d["end"], d["short"]) for d in numbered],
                            y=y_dom, height=h, fontsize=ps.FONT["small"],
                            color="white", levels=1, above_only=True)
        ax.text(label_x, y_dom, "Representative domain instances", ha="left",
                va="center", fontsize=ps.FONT["label"], color=ps.INK)
        handles += [ps.legend_patch(d["color"], d["legend"]) for d in numbered]

    if exon_blocks:
        for i, b in enumerate(exon_blocks):
            s = _int(b.get("start", b.get("protein_start_aa")))
            e = _int(b.get("end", b.get("protein_end_aa")))
            ax.add_patch(plt.Rectangle((s, y_exon - h * 0.8), max(1, e - s), 1.6 * h,
                                       facecolor=EXON_BLOCK_COLORS[i % 2],
                                       edgecolor=ps.INK, lw=0.4, zorder=4))
        for b in boundaries or ():
            pos = _int(b.get("boundary_position_aa") or b.get("protein_position")
                       or b.get("start"))
            ax.plot([pos, pos], [y_exon - h * 1.5, y_exon + h * 1.5],
                    color=BOUNDARY_TICK_COLOR, lw=0.9, zorder=6)
        ax.text(label_x, y_exon, "Coding exons and internal boundaries", ha="left",
                va="center", fontsize=ps.FONT["label"], color=ps.INK)
        handles.append(ps.legend_patch(EXON_BLOCK_COLORS[1], "coding exon"))
        if boundaries:
            handles.append(ps.legend_line(BOUNDARY_TICK_COLOR,
                                          "internal coding-exon boundary", lw=0.9))

    for i, r in enumerate(table):
        y = n - 1 - i
        ax.add_patch(plt.Rectangle((r["aa_start"], y - h),
                                   max(1, r["aa_end"] - r["aa_start"]), 2 * h,
                                   facecolor=CANDIDATE_COLOR, edgecolor=CANDIDATE_EDGE,
                                   lw=0.6, zorder=5))
        edge = r["nearest_edge_position_aa"]
        if edge != "" and r["distance_to_nearest_edge_aa"] not in ("", 0):
            anchor = (r["aa_start"] if abs(r["aa_start"] - edge)
                      <= abs(r["aa_end"] - edge) else r["aa_end"])
            ax.plot([anchor, edge], [y, y], color=CANDIDATE_EDGE, lw=0.8,
                    ls=(0, (2.5, 1.5)), zorder=4)
        if edge != "":
            ax.scatter([edge], [y], s=16, marker="D", color=ps.INK, zorder=6)
        _stacked_labels(ax, label_x, y, [
            (f"{r['candidate_label']} · aa {r['aa_start']}–{r['aa_end']} · "
             f"{r['length_aa']} aa", ps.INK, "bold"),
            (r["relation"], ps.MUTED, "normal"),
            (f"at coding-exon boundary {r['coincident_exon_boundary']}"
             if r["coincident_exon_boundary"] else "", ps.MUTED, "normal"),
        ], step=0.20)
        # The nearest-edge statement belongs next to the edge it describes, where
        # there is room for it, rather than truncated in the row label.
        if r["nearest_domain_label"]:
            note = (f"nearest edge: {r['nearest_domain_label']} "
                    f"{r['nearest_edge_type']} at aa "
                    f"{r['nearest_edge_position_aa']} · "
                    f"{r['distance_to_nearest_edge_aa']} aa away")
            anchor = max(r["aa_end"], _int(r["nearest_edge_position_aa"]))
            if anchor > max_aa * 0.62:
                ax.text(min(r["aa_start"], _int(r["nearest_edge_position_aa"]))
                        - max_aa * 0.012, y, note, ha="right", va="center",
                        fontsize=ps.FONT["small"], color=ps.MUTED)
            else:
                ax.text(anchor + max_aa * 0.012, y, note, ha="left", va="center",
                        fontsize=ps.FONT["small"], color=ps.MUTED)
    ax.text(label_x, n - 1 + 0.44, f"{EXPLORATORY_TAG} regions", ha="left",
            va="bottom", fontsize=ps.FONT["label"], color=ps.INK)
    handles.append(ps.legend_patch(CANDIDATE_COLOR,
                                   f"{EXPLORATORY_TAG.lower()} region · {VALIDATION_TAG}"))
    handles.append(_legend_marker("D", "nearest representative-domain edge", ps.INK))

    ax.set_yticks([])
    ax.spines[["top", "right", "left"]].set_visible(False)
    _readable_aa_axis(ax, max_aa,
                      f"Amino-acid position on {protein_id or 'the primary protein'} (aa)")
    species = f"{species_name} · " if species_name else ""
    figure_title(
        ax, f"{gene_symbol} — exploratory candidates in their domain context",
        subtitle=f"{species}{protein_id + ' · ' if protein_id else ''}{n} candidate "
                 f"region(s) · {len(numbered)} representative domain instance(s) · "
                 f"domains resolved by instance, not by accession",
        note=f"Do the exploratory candidate regions fall inside annotated domains or "
             f"between them? {EXPLORATORY_TAG} regions only — {VALIDATION_TAG}. A "
             f"positional overlap with a domain is an observation, not a functional "
             f"claim.")
    ncol = 2
    shared_legend(ax, handles, ncol=ncol, loc="upper left",
                  bbox=(0.0, _below_axes(ax, 34)))
    _footnote(ax, footnote or "",
              points=34 + 12 * (1 + (len(handles) - 1) // ncol) + 8)
    save_figure_all_formats(fig, fig_dir, stem)
    return True
