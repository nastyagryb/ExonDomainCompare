from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence

# The seven canonical gallery categories. A card must use exactly one of them.
CATEGORIES = (
    "Exon structure",
    "Isoform analysis",
    "Domain architecture",
    "Exon–domain boundaries",
    "Genomic context",
    "Exploratory candidates",
    "Supplements",
)

COORDINATE_SYSTEM = "protein 1-based amino-acid positions"
GENOMIC_COORDINATE_SYSTEM = "assembly genomic coordinates (NCBI RefSeq annotation)"

EXPLORATORY_STATUS = "Exploratory analysis; not biologically validated."
DESCRIPTIVE_STATUS = "Descriptive annotation summary; no event is claimed as validated."


def build_caption(*, gene: str, species: str, protein_id: str, description: str,
                  coordinate_system: str = COORDINATE_SYSTEM,
                  annotation_source: str = "",
                  threshold: str = "",
                  status: str = EXPLORATORY_STATUS) -> str:
    head = f"{gene or 'Gene'}"
    if species:
        head += f" ({species})"
    if protein_id:
        head += f", primary protein {protein_id}"
    parts = [f"{head} — {description.rstrip('.')}."]
    coord = coordinate_system
    if coord and protein_id and coord == COORDINATE_SYSTEM:
        coord = f"{coord} on {protein_id}"
    if coord:
        parts.append(f"Coordinates: {coord}.")
    if annotation_source:
        parts.append(f"Annotation source: {annotation_source}.")
    if threshold:
        parts.append(f"Threshold: {threshold}.")
    if status:
        parts.append(status)
    return " ".join(parts)


def write_caption_file(fig_dir: Path, stem: str, caption: str) -> Optional[Path]:
    if not caption:
        return None
    path = Path(fig_dir) / f"{stem}.caption.txt"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(caption.strip() + "\n", encoding="utf-8")
    return path


def join_sources(sources: Sequence[str]) -> str:
    return " + ".join([s for s in sources if s])


# --------------------------------------------------------------------------- #
# Gallery cards
# --------------------------------------------------------------------------- #
GALLERY_INDEX_FILES = ("figures_index.json", "generic/figures_index.json")


def file_url(run_id: str, run_relative_path: str, *, inline: bool = False) -> str:
    suffix = "&inline=true" if inline else ""
    return f"/api/runs/{run_id}/files?path={run_relative_path}{suffix}"


def figure_card(*, figure_id: str, title: str, category: str, run_id: str,
                figure_dir: str, stem: str, scientific_question: str,
                interpretation: str, caption: str = "", kind: str = "main",
                gene_symbol: str = "", species: str = "", species_id: str = "",
                protein_id: str = "", transcript_id: str = "",
                stage: str = "post_cluster", has_table: bool = True,
                figure_type: str = "",
                source_files: Sequence[str] = ()) -> Dict[str, Any]:
    base = f"{figure_dir.rstrip('/')}/{stem}"
    card: Dict[str, Any] = {
        "figure_id": figure_id,
        "figure_type": figure_type or figure_id,
        "title": title,
        "category": category,
        "section": category,
        "kind": kind,
        "scientific_question": scientific_question,
        "interpretation": interpretation,
        "caption": caption or interpretation,
        "gene_symbol": gene_symbol,
        "species": species,
        "species_id": species_id,
        "protein_id": protein_id,
        "transcript_id": transcript_id,
        "stage": stage,
        "status": "available",
        "renderer": "matplotlib_shared_gene_plots",
        "png_url": file_url(run_id, f"{base}.png", inline=True),
        "svg_url": file_url(run_id, f"{base}.svg"),
        "pdf_url": file_url(run_id, f"{base}.pdf"),
        "source_files": list(source_files),
        "error": "",
    }
    if has_table:
        card["table_url"] = file_url(run_id, f"{base}.tsv")
        card["source_table"] = list(source_files)[0] if source_files else ""
    return card


def _card_id(entry: Any) -> str:
    if isinstance(entry, dict):
        return str(entry.get("figure_id") or entry.get("id") or "")
    return str(entry or "")


def _sort_key(card: Dict[str, Any]) -> tuple:
    category = card.get("category") or card.get("section") or ""
    index = CATEGORIES.index(category) if category in CATEGORIES else len(CATEGORIES)
    return (index, 0 if card.get("kind") == "main" else 1, _card_id(card))


def _available_entry(card: Dict[str, Any]) -> Dict[str, Any]:
    entry = {"id": _card_id(card)}
    for key in ("title", "category", "kind", "figure_type", "caption",
                "scientific_question", "interpretation", "stage", "status",
                "png_url", "svg_url", "pdf_url", "table_url", "source_files"):
        if card.get(key) not in (None, ""):
            entry[key] = card[key]
    entry.setdefault("status", "available")
    entry["error"] = card.get("error", "")
    return entry


FIGURE_SUFFIXES = (".svg", ".pdf", ".png", ".tsv", ".caption.txt")


def remove_retired_figure_files(figure_dir: Path, stems: Iterable[str]) -> int:
    removed = 0
    directory = Path(figure_dir)
    for stem in stems:
        for suffix in FIGURE_SUFFIXES:
            path = directory / f"{stem}{suffix}"
            if path.is_file():
                path.unlink()
                removed += 1
    return removed


def register_gallery_cards(run_dir: Path, cards: Sequence[Dict[str, Any]], *,
                           drop_figure_ids: Iterable[str] = ()) -> int:
    drop = {i for i in drop_figure_ids if i}
    own = {_card_id(c) for c in cards}
    written = 0
    for name in GALLERY_INDEX_FILES:
        path = Path(run_dir) / "website_indices" / name
        if not path.exists():
            continue
        try:
            doc = json.loads(path.read_text())
        except (OSError, ValueError):
            continue
        kept: List[Dict[str, Any]] = [
            f for f in (doc.get("figures") or [])
            if _card_id(f) not in drop and _card_id(f) not in own]
        doc["figures"] = sorted(kept + list(cards), key=_sort_key)
        available = doc.get("available")
        if isinstance(available, list):
            kept_available = [a for a in available if _card_id(a) not in drop
                              and _card_id(a) not in own]
            # Only the root index mirrors the figure cards in ``available``. In the
            # generic index the same key lists view names as plain strings, which are
            # not figures and must not be given figure entries.
            if all(isinstance(a, dict) for a in available):
                kept_available += [_available_entry(c) for c in cards]
            doc["available"] = kept_available
        path.write_text(json.dumps(doc, indent=2))
        written += 1
    # Registration is part of writing a card, not a step a generator may forget: a
    # card without its run and schema identity cannot be attributed to a run, and the
    # FGFR2 Gallery was the only writer that normalised afterwards.
    if written:
        figures_path = Path(run_dir) / "website_indices" / "figures_index.json"
        try:
            figures_doc = json.loads(figures_path.read_text())
            gallery_available = any(
                isinstance(card, dict) and card.get("status", "available") == "available"
                for card in figures_doc.get("figures", [])
            )
        except (OSError, ValueError):
            gallery_available = False
        if gallery_available:
            for relative in ("available_views.json", "generic/available_views.json"):
                view_path = Path(run_dir) / "website_indices" / relative
                try:
                    view_doc = json.loads(view_path.read_text())
                except (OSError, ValueError):
                    continue
                views = view_doc.get("available_views") or view_doc.get("views")
                if isinstance(views, dict):
                    views["figure_gallery"] = True
                    if "figures" in views:
                        views["figures"] = True
                    view_path.write_text(json.dumps(view_doc, indent=2) + "\n")
        try:
            from plotting.figure_registration import normalise_run
            normalise_run(Path(run_dir))
        except Exception:  # noqa: BLE001 - registration must not lose the cards
            pass
    return written
