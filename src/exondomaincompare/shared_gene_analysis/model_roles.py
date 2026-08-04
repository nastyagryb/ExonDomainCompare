from __future__ import annotations

from typing import Any, Dict, List, Optional, Sequence

#: The one model of a species that carries no isoform-level claim of its own. This
#: is what a generic gene's selected primary protein is.
PRIMARY_REFERENCE = "primary_reference"

#: A model that *is* a specific validated isoform. The suffix is the isoform label
#: as the analysis validated it, so the role is readable without a lookup.
VALIDATED_ISOFORM_PREFIX = "validated_isoform_"

#: A model whose coding sequence had to be reconstructed rather than read directly.
RECONSTRUCTED_ISOFORM = "reconstructed_isoform"

#: An additional protein model kept for orientation, carrying no validated claim.
EXPLORATORY_ALTERNATIVE = "exploratory_alternative_model"

ROLE_LABEL = {
    PRIMARY_REFERENCE: "Primary reference model",
    RECONSTRUCTED_ISOFORM: "Reconstructed isoform model",
    EXPLORATORY_ALTERNATIVE: "Exploratory alternative model (not validated)",
}


def validated_isoform_role(isoform_label: str) -> str:
    return f"{VALIDATED_ISOFORM_PREFIX}{str(isoform_label).strip()}"


def is_validated_isoform_role(role: Optional[str]) -> bool:
    return str(role or "").startswith(VALIDATED_ISOFORM_PREFIX)


def isoform_of_role(role: Optional[str]) -> str:
    r = str(role or "")
    return r[len(VALIDATED_ISOFORM_PREFIX):] if is_validated_isoform_role(r) else ""


def role_label(role: Optional[str]) -> str:
    r = str(role or "")
    if is_validated_isoform_role(r):
        return f"Validated {isoform_of_role(r)} isoform model"
    return ROLE_LABEL.get(r, r or "Unspecified model role")


def known_role(role: Optional[str]) -> bool:
    r = str(role or "")
    return bool(r) and (r in ROLE_LABEL or is_validated_isoform_role(r))


def model_id(gene_symbol: str, species_id: str, discriminator: str = "") -> str:
    gene = str(gene_symbol or "gene").strip().lower()
    disc = str(discriminator or "primary").strip().replace(" ", "_")
    return f"{gene}:{species_id}:{disc}"


def role_errors(models: Sequence[Dict[str, Any]]) -> List[str]:
    errors: List[str] = []
    by_species: Dict[str, List[Dict[str, Any]]] = {}
    seen_ids: Dict[str, str] = {}

    for m in models:
        species = str(m.get("species_id") or "?")
        mid = str(m.get("model_id") or "")
        role = m.get("model_role")
        by_species.setdefault(species, []).append(m)

        if not mid:
            errors.append(f"[{species}/{m.get('protein_id')}] model_id is missing; a "
                          f"renderer must never have to infer model identity")
        elif mid in seen_ids and seen_ids[mid] != str(m.get("protein_id") or ""):
            errors.append(f"model_id {mid} is used by two different proteins "
                          f"({seen_ids[mid]} and {m.get('protein_id')})")
        else:
            seen_ids[mid] = str(m.get("protein_id") or "")

        if not known_role(role):
            errors.append(f"[{species}/{m.get('protein_id')}] model_role "
                          f"{role!r} is not one of the declared roles")

    for species, group in by_species.items():
        primaries = [m for m in group if m.get("is_primary_reference")]
        if len(primaries) != 1:
            errors.append(
                f"[{species}] exactly one model must be the primary reference, found "
                f"{len(primaries)} of {len(group)}. A comparative figure puts one row "
                f"per species, so the choice has to be recorded, not left to order.")
        roles = [str(m.get("model_role") or "") for m in group]
        if len(set(roles)) != len(roles):
            errors.append(f"[{species}] two models share the role "
                          f"{sorted(r for r in roles if roles.count(r) > 1)[0]!r}; "
                          f"distinct models need distinct roles")
    return errors
