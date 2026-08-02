"""How a run is named, identified and ordered.

Two names, deliberately separate:

* ``run_id`` is the technical identifier. It is a directory name, so it stays
  timestamped and filesystem-safe, and it never changes.
* ``run_name`` is what the user typed. It may contain spaces, capitals and
  ordinary punctuation, it is never used as a path, and two runs may share one.

When a user gives no name the interface still needs a title, so
:func:`display_label` builds a readable one from the biology of the run — the
gene and its species — rather than falling back on a placeholder like
``custom_run``, which tells a reader nothing and used to be written into the
stored metadata as if the user had chosen it.
"""
from __future__ import annotations

import re
from datetime import datetime, timezone
from typing import Any, Dict, Iterable, List, Optional, Sequence

#: A visible run name longer than this is truncated. Long enough for a sentence,
#: short enough to stay on one line in a card header.
MAX_RUN_NAME = 120

#: Characters a visible run name may contain beyond letters, digits and spaces.
_ALLOWED_PUNCTUATION = r"\-_.,:;()\[\]/+&'\u2019"
_DISALLOWED = re.compile(rf"[^\w\s{_ALLOWED_PUNCTUATION}]", re.UNICODE)
_WHITESPACE = re.compile(r"\s+")

#: Placeholders older versions stored when the user left the name empty. They
#: are names of nothing, so they are treated as "no name given".
_PLACEHOLDER_NAMES = {
    "custom_run", "custom run", "run", "new_run", "new run",
    "full30_run", "pilot_run", "untitled",
}

_RUN_ID_STAMP = re.compile(r"^(\d{4})-(\d{2})-(\d{2})[_-](\d{2})(\d{2})")


#: Suffixes the old auto-naming appended when the user gave no name.
_GENERATED_SUFFIXES = ("_core_pilot", "_pilot", "_run", "_twospecies", "_core_only_pilot")


def is_generated_name(name: str, gene_symbol: str = "") -> bool:
    """Whether a stored name was produced by the pipeline rather than by a user.

    Older runs recorded a default such as ``fgfr1_gallus_mus_core_pilot`` in the
    same field a user types into. Those names describe the run the way a
    directory does, and re-deriving the label from the gene and species reads
    better, so they are treated as "no name given". A name is taken as generated
    only when it is already a bare slug — no spaces, no capitals — and it either
    starts with the gene symbol or carries one of the old suffixes. Anything a
    person is likely to have typed keeps its exact wording.
    """
    text = (name or "").strip()
    if not text or text != slugify(text):
        return False
    gene = slugify(gene_symbol)
    return (bool(gene) and text.startswith(f"{gene}_")) \
        or text.endswith(_GENERATED_SUFFIXES)


def clean_run_name(raw: Optional[str], gene_symbol: str = "") -> str:
    """Normalise a user-entered run name for storage.

    Trims, collapses internal whitespace, drops characters that have no place in
    a label, and caps the length. Returns "" when nothing meaningful is left, so
    an empty name, a placeholder and a pipeline-generated default are the same
    thing downstream.
    """
    # Collapse whitespace after dropping characters too, so removing a dash from
    # "TP53 — mammals" does not leave a double space behind.
    text = _WHITESPACE.sub(" ", _DISALLOWED.sub("", str(raw or ""))).strip()
    if not text or text.lower() in _PLACEHOLDER_NAMES:
        return ""
    if is_generated_name(text, gene_symbol):
        return ""
    return text[:MAX_RUN_NAME].strip()


def slugify(text: str) -> str:
    """A filesystem-safe fragment for a run directory name."""
    slug = re.sub(r"[^a-z0-9]+", "_", str(text or "").lower()).strip("_")
    return re.sub(r"_{2,}", "_", slug)


def run_id_slug(run_name: str = "", *, gene_symbol: str = "",
                species: Sequence[str] = ()) -> str:
    """The descriptive part of a ``run_id``.

    A user-entered name wins, because that is what its author will look for in
    the directory listing. Otherwise the slug describes the run: gene plus one
    or two species, or gene plus a species count.
    """
    named = slugify(clean_run_name(run_name))
    if named:
        return named[:60].strip("_") or "run"
    gene = slugify(gene_symbol)
    ids = [slugify(s) for s in species if s]
    if gene and len(ids) == 1:
        return f"{gene}_{ids[0]}"[:60].strip("_")
    if gene and 2 <= len(ids) <= 2:
        return f"{gene}_{ids[0]}_{ids[1]}"[:60].strip("_")
    if gene and ids:
        return f"{gene}_{len(ids)}species"
    return gene or "run"


def short_species(species_id: str) -> str:
    """`gallus_gallus` -> `Gallus gallus`; only the genus is capitalised."""
    parts = [p for p in str(species_id or "").replace(" ", "_").split("_") if p]
    if not parts:
        return ""
    return " ".join([parts[0].capitalize(), *[p.lower() for p in parts[1:]]])


def species_summary(species: Sequence[str], limit: int = 2) -> str:
    """A concise species phrase: two names, or the first plus a remainder count."""
    names = [short_species(s) for s in species if s]
    if not names:
        return ""
    if len(names) <= limit:
        return " + ".join(names)
    return f"{names[0]} + {len(names) - 1} more"


def display_label(run_name: Optional[str] = None, *, gene_symbol: str = "",
                  species: Sequence[str] = (), species_count: int = 0,
                  run_id: str = "") -> str:
    """The title a run card shows.

    The user's name when there is one; otherwise a readable description built
    from the gene and its species. Never a placeholder, and never the bare
    ``run_id`` unless there is genuinely nothing else to say.
    """
    name = clean_run_name(run_name, gene_symbol)
    if name:
        return name
    gene = (gene_symbol or "").strip().upper()
    summary = species_summary(species)
    if gene and summary:
        return f"{gene} · {summary}"
    count = species_count or len(list(species))
    if gene and count:
        return f"{gene} · {count} species"
    if gene:
        return f"{gene} analysis"
    return run_id or "Run"


def analysis_mode(species_count: int) -> str:
    return "Single-species analysis" if species_count <= 1 else "Comparative analysis"


def run_id_timestamp(run_id: str) -> Optional[datetime]:
    """The creation time encoded in a legacy ``run_id``, if it has one."""
    match = _RUN_ID_STAMP.match(str(run_id or ""))
    if not match:
        return None
    year, month, day, hour, minute = (int(g) for g in match.groups())
    try:
        return datetime(year, month, day, hour, minute, tzinfo=timezone.utc)
    except ValueError:
        return None


def _parse_iso(value: Any) -> Optional[datetime]:
    if not value:
        return None
    text = str(value).replace("Z", "+00:00")
    try:
        parsed = datetime.fromisoformat(text)
    except ValueError:
        return None
    return parsed if parsed.tzinfo else parsed.replace(tzinfo=timezone.utc)


_EPOCH = datetime(1970, 1, 1, tzinfo=timezone.utc)


def creation_time(run: Dict[str, Any]) -> datetime:
    """The best available creation time, in the order the ordering rule requires.

    ``created_at`` first; the timestamp encoded in a legacy ``run_id`` second;
    a registry timestamp third. Filesystem modification time is never used — it
    changes whenever a run is refreshed or re-indexed and would reorder the list
    behind the user's back.
    """
    for key in ("created_at", "registry_created_at"):
        parsed = _parse_iso(run.get(key))
        if parsed:
            return parsed
    from_id = run_id_timestamp(run.get("run_id", ""))
    return from_id or _EPOCH


def sort_runs(runs: Iterable[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """Newest first, with ``run_id`` as the stable tie-breaker."""
    return sorted(runs, key=lambda r: (creation_time(r), str(r.get("run_id", ""))),
                  reverse=True)


#: Card grouping. Ordering inside a group is still newest first.
ACTIVE_STATUSES = {
    "created", "running", "pre_interpro_running", "post_interpro_running",
    "cluster_required", "cluster_running", "cluster_fetch_complete",
}
ATTENTION_STATUSES = {
    "failed", "core_model_collection_failed", "incomplete", "stopped", "partial",
}


def run_group(status: str) -> str:
    if status in ATTENTION_STATUSES:
        return "attention"
    if status == "results_ready":
        return "completed"
    return "active" if status in ACTIVE_STATUSES else "completed"


def completion_summary(*, primary_fasta_count: int = 0,
                       available_views: Optional[Dict[str, Any]] = None,
                       species_count: int = 0) -> str:
    """One sentence replacing the completed-stage checklist on a finished run."""
    views = available_views or {}
    protein = (f"{primary_fasta_count} primary protein"
               f"{'' if primary_fasta_count == 1 else 's'} analysed"
               if primary_fasta_count else
               f"{species_count} species analysed" if species_count else "Analysis complete")
    ready = [label for key, label in (
        ("domain_architecture", "Domain architecture"),
        ("boundary", "Boundary analysis"),
        ("comparative", "Comparative views"),
    ) if views.get(key)]
    if not ready:
        return f"{protein}."
    if len(ready) == 1:
        return f"{protein}. {ready[0]} is available."
    return f"{protein}. {', '.join(ready[:-1])} and {ready[-1]} are available."


def describe_failure(*, failed_stage: str = "", failed_species: str = "",
                     last_error: str = "") -> str:
    """A concise, actionable failure line for a card — not a terminal log."""
    stage_label = {
        "core_model_collection": "gene and protein model collection",
        "gene_locus_resolution": "gene locus resolution",
        "core_primary_fasta": "primary protein FASTA",
        "domain_architecture": "domain architecture",
        "boundary": "boundary analysis",
    }.get(failed_stage, failed_stage.replace("_", " "))
    who = short_species(failed_species) if failed_species else ""
    if who and stage_label:
        return f"{who}: {stage_label} could not be generated."
    if stage_label:
        return f"{stage_label.capitalize()} could not be generated."
    first_line = str(last_error or "").strip().splitlines()[:1]
    return first_line[0][:200] if first_line else "The run stopped before producing results."


def run_display_fields(cfg: Dict[str, Any], *, run_id: str, species_count: int,
                       gene_symbol: str, species: Sequence[str] = (),
                       status: str = "") -> Dict[str, Any]:
    """The label fields a run summary exposes, derived once for every caller."""
    name = clean_run_name(cfg.get("run_name"), gene_symbol)
    return {
        "run_name": name,
        "display_name": display_label(name, gene_symbol=gene_symbol, species=species,
                                      species_count=species_count, run_id=run_id),
        "species_summary": species_summary(species),
        "analysis_mode": analysis_mode(species_count),
        "group": run_group(status),
    }
