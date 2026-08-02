"""The interactive viewers and the publication renderer are one figure.

These tests enforce the contract in `docs/architecture/figure_parity_contract.md`.
They compare the *semantic* representation — style keys, labels, class vocabulary,
feature identity — rather than DOM trees, because the two sides legitimately use
different rendering technology and only the science has to agree.

Deliberately not asserted: identical markup, identical pixel geometry. See the
contract for why, and for what remains convention.
"""

from __future__ import annotations

import json
import re
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
VIEWERS = ROOT / "webapp" / "frontend" / "src" / "pages" / "viewers"
FRONTEND = ROOT / "webapp" / "frontend" / "src"
SPEC = VIEWERS / "semanticStyles.js"
CONTRACT = ROOT / "docs" / "architecture" / "figure_parity_contract.md"

FGFR1_RUN = "2026-07-23_1100_fgfr1_gallus_core_pilot"
TP53_RUN = "2026-07-21_1436_custom_run"

REQUIRED_KEYS = [
    "coding_exon", "alternative_exon", "shared_exon", "shifted_boundary",
    "representative_domain", "family_superfamily", "tm_helix", "candidate_region",
    "boundary_exact", "boundary_near", "boundary_inside", "boundary_outside",
    "boundary_uncertain", "selected_feature",
    "primary_sequence", "alternative_sequence", "gap",
    "variable_region", "conserved_region",
]
REQUIRED_PROPS = ["fill", "stroke", "strokeWidth", "opacity", "text", "marker",
                  "labelPriority", "printFallback"]

# The frozen FGFR2 vocabularies keep their own palettes on purpose (see contract):
# boundary.js holds the FGFR2 Boundary Consistency classes, fgfr2Styles.js the FGFR2
# architecture, topology and IIIb/IIIc cassette vocabulary. Both encode features the
# generic specification does not have, so they are centralised separately rather than
# folded in.
FROZEN_VOCABULARY = {"boundary.js", "fgfr2Styles.js"}


# --------------------------------------------------------------------------- #
# Reading the shared specification from Node, so the test sees resolved values
# --------------------------------------------------------------------------- #

def _node(expr: str) -> dict:
    """Evaluate an expression against the shared spec and return its JSON."""
    script = (
        "import * as S from './semanticStyles.js';\n"
        f"process.stdout.write(JSON.stringify({expr}));\n"
    )
    tmp = VIEWERS / "_parity_probe.mjs"
    tmp.write_text(script, encoding="utf-8")
    try:
        out = subprocess.run([_node_bin(), str(tmp)], capture_output=True, text=True,
                             cwd=VIEWERS)
        if out.returncode != 0:
            raise AssertionError(f"probing the shared spec failed:\n{out.stderr}")
        return json.loads(out.stdout)
    finally:
        tmp.unlink(missing_ok=True)


def _node_bin() -> str:
    from shutil import which
    node = which("node")
    if not node:
        pytest.skip("node is required to evaluate the shared visual specification")
    return node


@pytest.fixture(scope="module")
def spec() -> dict:
    return _node("{ styles: S.FEATURE_STYLES, text: S.TEXT_ROLES, "
                 "classLabel: S.BOUNDARY_CLASS_LABEL, classStyle: S.BOUNDARY_CLASS_STYLE, "
                 "domains: S.DOMAIN_INSTANCE_COLOURS, keys: S.FEATURE_KEYS }")


# --------------------------------------------------------------------------- #
# 1. The specification itself
# --------------------------------------------------------------------------- #

def test_the_specification_covers_every_required_semantic_feature(spec):
    missing = [k for k in REQUIRED_KEYS if k not in spec["styles"]]
    assert not missing, f"semantic keys missing from the shared spec: {missing}"


@pytest.mark.parametrize("key", REQUIRED_KEYS)
def test_every_semantic_style_is_fully_specified(spec, key):
    style = spec["styles"][key]
    for prop in REQUIRED_PROPS:
        assert prop in style, f"{key} does not specify {prop}"
    assert style["fill"], f"{key} has an empty fill"
    assert isinstance(style["labelPriority"], (int, float))
    assert re.fullmatch(r"#[0-9a-fA-F]{6}", style["printFallback"]), \
        f"{key} has no print-safe grey fallback"


def test_the_declared_key_list_matches_the_implemented_styles(spec):
    assert sorted(spec["keys"]) == sorted(spec["styles"].keys()), \
        "FEATURE_KEYS and FEATURE_STYLES disagree, so a key can be silently missing"


def test_an_unknown_key_raises_instead_of_painting_a_default():
    """A silent default is how a feature ends up black in an exported SVG."""
    text = SPEC.read_text(encoding="utf-8")
    body = text.split("export function featureStyle", 1)[1].split("\n}", 1)[0]
    assert "throw new Error" in body, "featureStyle falls back instead of failing loudly"


def test_selection_is_an_outline_and_never_a_recolouring(spec):
    sel = spec["styles"]["selected_feature"]
    assert sel["fill"] == "none", \
        "a selected feature is refilled, which hides the class the reader is judging"
    assert sel["strokeWidth"] > spec["styles"]["coding_exon"]["strokeWidth"]


# --------------------------------------------------------------------------- #
# 2. No component owns a scientific colour
# --------------------------------------------------------------------------- #

def _component_files() -> list[Path]:
    files = sorted(VIEWERS.glob("*.jsx")) + sorted(VIEWERS.glob("*.js"))
    files += [FRONTEND / "pages" / "GlobalBoundaryDashboard.jsx",
              FRONTEND / "pages" / "BoundaryPage.jsx",
              FRONTEND / "components" / "shared" / "index.jsx"]
    # The spec and the figure builders legitimately name colours; the export
    # utilities legitimately name a white raster background.
    skip = {"semanticStyles.js", "figureSpec.js", "mainFigures.js", "alignmentFigure.js",
            "figureExport.js", "plotExport.js"}
    return [f for f in files
            if f.exists() and f.name not in skip and f.name not in FROZEN_VOCABULARY]


def _strip_comments(text: str) -> str:
    text = re.sub(r"/\*.*?\*/", "", text, flags=re.DOTALL)
    return "\n".join(re.sub(r"//.*$", "", line) for line in text.splitlines())


def test_no_component_hardcodes_a_colour_the_shared_spec_owns(spec):
    """Duplicated literals are how the two sides drifted apart in the first place."""
    owned = {v.lower() for s in spec["styles"].values()
             for v in (s["fill"], s["stroke"], s["text"]) if isinstance(v, str)
             and v.startswith("#")}
    owned |= {c.lower() for c in spec["domains"]}
    owned |= {r["fill"].lower() for r in spec["text"].values()}
    # Plain paper and ink encode no scientific meaning: a white panel background or
    # a near-black body text is not a feature colour, and forbidding them would be
    # noise rather than a guard.
    owned -= {"#ffffff", "#1c2433"}

    offenders: dict[str, list[str]] = {}
    for path in _component_files():
        code = _strip_comments(path.read_text(encoding="utf-8"))
        hits = sorted({m.group(0).lower()
                       for m in re.finditer(r"#[0-9a-fA-F]{6}\b", code)
                       if m.group(0).lower() in owned})
        if hits:
            offenders[path.name] = hits
    assert not offenders, (
        "these components repeat a colour the shared specification already owns; "
        f"import it instead: {offenders}")


def test_the_generic_boundary_vocabulary_has_no_palette_of_its_own():
    code = _strip_comments((VIEWERS / "boundaryClasses.js").read_text(encoding="utf-8"))
    assert not re.search(r"#[0-9a-fA-F]{6}", code), \
        "boundaryClasses.js defines colours again instead of deriving them"
    assert "semanticStyles" in code, "boundaryClasses.js is not wired to the shared spec"


def test_the_export_renderer_resolves_classes_through_the_shared_spec():
    code = (VIEWERS / "mainFigures.js").read_text(encoding="utf-8")
    assert "semanticStyles" in code, "mainFigures.js does not read the shared spec"
    assert "boundaryStyleKey" in code, "boundary classes are resolved outside the spec"
    assert "domainInstanceFill" in code, "domain instances get their own colour ramp"
    assert not re.search(r'DOMAIN_COLOURS\s*=\s*\[', code), \
        "mainFigures.js still owns a domain colour ramp"


# --------------------------------------------------------------------------- #
# 3. No scientific mark can go black without its stylesheet
# --------------------------------------------------------------------------- #

def test_no_scientific_svg_mark_depends_on_the_stylesheet():
    sys.path.insert(0, str(ROOT / "tests"))
    from audit_svg_css_dependencies import audit

    report = audit()
    assert report["total"] == 0, (
        "these marks lose their paint in a standalone SVG: "
        f"{json.dumps(report['by_file'], indent=2)}")


# --------------------------------------------------------------------------- #
# 4. Semantic parity of the shipped figures
# --------------------------------------------------------------------------- #

def _cards(run_id: str) -> list[dict]:
    p = ROOT / "runs" / run_id / "website_indices" / "figures_index.json"
    if not p.is_file():
        pytest.skip(f"reference run {run_id} is not present")
    return json.loads(p.read_text())["figures"]


def _svg_of(run_id: str, needle: str) -> str:
    for card in _cards(run_id):
        if needle in (card.get("figure_id") or "") and card.get("svg_url"):
            rel = re.sub(r"^.*?path=", "", card["svg_url"])
            from urllib.parse import unquote
            path = ROOT / "runs" / run_id / unquote(rel)
            if path.is_file():
                return path.read_text(encoding="utf-8")
    pytest.skip(f"no shipped SVG for {needle} in {run_id}")


def _primary_model(run_id: str) -> dict:
    """The primary protein model of a single-species run, as the figures read it."""
    p = (ROOT / "runs" / run_id / "website_indices" / "generic"
         / "protein_coordinate_model.json")
    if not p.is_file():
        pytest.skip("coordinate model missing")
    models = json.loads(p.read_text()).get("models") or []
    if not models:
        pytest.skip("coordinate model has no protein models")
    return models[0]


@pytest.mark.parametrize("run_id", [FGFR1_RUN, TP53_RUN])
def test_the_boundary_figure_paints_classes_in_the_shared_colours(spec, run_id):
    svg = _svg_of(run_id, "boundary_on_architecture")
    used = {m.group(1).lower() for m in re.finditer(r'fill="(#[0-9a-fA-F]{6})"', svg)}
    used |= {m.group(1).lower() for m in re.finditer(r'stroke="(#[0-9a-fA-F]{6})"', svg)}

    model = _primary_model(run_id)
    classes = {b.get("boundary_class") for b in model.get("exon_boundaries", [])}
    classes = {c for c in classes if c}
    assert classes, "the run has no classified boundaries to check"

    for cls in classes:
        key = spec["classStyle"].get(cls, "boundary_uncertain")
        colour = spec["styles"][key]["fill"].lower()
        assert colour in used, (
            f"class {cls} should be drawn in {colour} (semantic key {key}) "
            f"but that colour does not occur in the shipped figure")


@pytest.mark.parametrize("run_id", [FGFR1_RUN, TP53_RUN])
def test_the_class_legend_uses_the_shared_vocabulary(spec, run_id):
    svg = _svg_of(run_id, "boundary_class_summary")
    for cls, label in spec["classLabel"].items():
        if label.lower() in svg.lower():
            continue
        # A class absent from the data need not appear in its figure's legend.
        present = {b.get("boundary_class")
                   for b in _primary_model(run_id).get("exon_boundaries", [])}
        assert cls not in present, f"class {cls} occurs in the data but not in the legend"


@pytest.mark.parametrize("run_id", [FGFR1_RUN, TP53_RUN])
def test_candidate_labels_come_from_the_coordinate_model(run_id):
    regions = sorted(_primary_model(run_id).get("candidate_regions", []),
                     key=lambda r: r.get("start") or r.get("aa_start") or 0)
    if not regions:
        pytest.skip("the run has no candidate regions")
    expected = {r.get("id") or f"C{i}" for i, r in enumerate(regions, start=1)}

    for needle in ("primary_exon_projection", "integrated_domain_architecture"):
        svg = _svg_of(run_id, needle)
        found = set(re.findall(r"\bC(\d+)\b", svg))
        stray = {f"C{n}" for n in found} - expected
        assert not stray, f"{needle} names candidates the model does not define: {stray}"


def test_the_frozen_vocabularies_are_centralised_too():
    """A frozen palette is still a palette: it may not be scattered over components."""
    for name in FROZEN_VOCABULARY:
        assert (VIEWERS / name).is_file(), f"{name} is missing"
    fgfr2 = (VIEWERS / "fgfr2Styles.js").read_text(encoding="utf-8")
    for symbol in ("FGFR2_DOMAIN_FILL", "FGFR2_TM_FILL", "FGFR2_CASSETTE_FILL"):
        assert f"export const {symbol}" in fgfr2, f"{symbol} is not exported centrally"

    for component in ("DomainArchitecture.jsx", "BoundaryDetailTrack.jsx"):
        code = _strip_comments((VIEWERS / component).read_text(encoding="utf-8"))
        assert "fgfr2Styles" in code, f"{component} does not read the frozen vocabulary"
        assert not re.search(r'(?:iiib|iiic)\s*:\s*"#', code), \
            f"{component} redefines the IIIb/IIIc cassette colours locally"


def test_the_contract_document_exists_and_names_its_own_limits():
    text = CONTRACT.read_text(encoding="utf-8")
    for section in ("What both sides must share", "allowed to differ",
                    "still convention"):
        assert section.lower() in text.lower(), f"the contract omits '{section}'"
    assert "boundary.js" in text, "the contract does not record the frozen exception"
