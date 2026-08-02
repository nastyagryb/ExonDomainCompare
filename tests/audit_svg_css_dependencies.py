"""Audit which SVG marks in the React viewers depend on the page stylesheet.

An exported SVG is a standalone document: the application stylesheet is not
attached. Any shape that receives its `fill` only from a CSS class therefore
falls back to the SVG initial value, which is opaque black — this is the reason
exported exon blocks render as solid black rectangles.

Run:
    ./venv/bin/python tests/audit_svg_css_dependencies.py
"""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FRONTEND = ROOT / "webapp" / "frontend" / "src"
CSS = FRONTEND / "App.css"

# SVG shapes whose paint defaults to black when no fill is supplied.
SHAPES = ("rect", "circle", "path", "polygon", "ellipse", "text", "line", "polyline")
TAG_RE = re.compile(r"<(" + "|".join(SHAPES) + r")\b([^>]*?)/?>", re.DOTALL)
CLASS_RE = re.compile(r'className=(?:"([^"]*)"|\{`([^`]*)`\})')

# Helpers in semanticStyles.js that return explicit paint. A mark that spreads one
# of them carries a real fill at runtime, so counting it as CSS-dependent would
# measure the source text rather than the exported document.
PAINT_HELPERS = ("featureProps", "textProps", "boundaryProps")
SPREAD_CALL_RE = re.compile(r"\{\s*\.\.\.\s*(?:" + "|".join(PAINT_HELPERS) + r")\s*\(")
SPREAD_NAME_RE = re.compile(r"\{\s*\.\.\.\s*([A-Za-z_$][\w$]*)\s*\}")


def paint_constants(text: str) -> set[str]:
    """Names bound to a paint helper's result, e.g. `const AXIS = textProps("axis")`."""
    pattern = (r"(?:const|let|var)\s+([A-Za-z_$][\w$]*)\s*=\s*(?:"
               + "|".join(PAINT_HELPERS) + r")\s*\(")
    return set(re.findall(pattern, text))


def spreads_paint(attrs: str, constants: set[str]) -> bool:
    if SPREAD_CALL_RE.search(attrs):
        return True
    return any(m.group(1) in constants for m in SPREAD_NAME_RE.finditer(attrs))


def css_paint_rules() -> dict[str, dict[str, str]]:
    """class name -> {property: value} for fill/stroke declared in the stylesheet."""
    text = CSS.read_text(encoding="utf-8")
    rules: dict[str, dict[str, str]] = {}
    for selector, body in re.findall(r"([^{}]+)\{([^{}]*)\}", text):
        paint = {}
        for prop in ("fill", "stroke"):
            m = re.search(rf"(?<![-\w]){prop}\s*:\s*([^;]+);", body)
            if m:
                paint[prop] = m.group(1).strip()
        if not paint:
            continue
        for cls in re.findall(r"\.([\w-]+)", selector):
            rules.setdefault(cls, {}).update(paint)
    return rules


def audit() -> dict:
    paint_rules = css_paint_rules()
    findings: list[dict] = []

    for path in sorted(FRONTEND.rglob("*.jsx")):
        text = path.read_text(encoding="utf-8")
        if "<svg" not in text:
            continue
        rel = str(path.relative_to(FRONTEND))
        constants = paint_constants(text)
        for match in TAG_RE.finditer(text):
            tag, attrs = match.group(1), match.group(2)
            spread = spreads_paint(attrs, constants)
            has_fill = "fill=" in attrs or spread
            has_stroke = "stroke=" in attrs or spread
            cm = CLASS_RE.search(attrs)
            classes = " ".join(c for c in cm.groups() if c).split() if cm else []
            # Only classes that the stylesheet actually paints matter here.
            painted = [c for c in classes if c in paint_rules]
            styled_fill = [c for c in painted if "fill" in paint_rules[c]]
            styled_stroke = [c for c in painted if "stroke" in paint_rules[c]]

            problems = []
            if not has_fill and styled_fill:
                problems.append("fill only from stylesheet -> black in exported SVG")
            # `text` is included: an axis label that falls back to black is a label
            # the figure never asked for, and it is the failure that survives longest
            # unnoticed because black text still looks deliberate.
            if not has_fill and not styled_fill and tag in ("rect", "circle", "path",
                                                            "polygon", "ellipse",
                                                            "text"):
                problems.append("no fill at all -> black in exported SVG")
            if not has_stroke and styled_stroke:
                problems.append("stroke only from stylesheet -> lost in exported SVG")
            if not problems:
                continue
            findings.append({
                "file": rel,
                "line": text[:match.start()].count("\n") + 1,
                "tag": tag,
                "classes": classes,
                "problems": problems,
            })

    by_file: dict[str, int] = {}
    for f in findings:
        by_file[f["file"]] = by_file.get(f["file"], 0) + 1
    return {"total": len(findings), "by_file": by_file, "findings": findings}


if __name__ == "__main__":
    report = audit()
    print(f"CSS-dependent SVG marks: {report['total']}")
    for file, n in sorted(report["by_file"].items(), key=lambda kv: -kv[1]):
        print(f"  {n:4d}  {file}")
    out = ROOT / "tmp" / "svg_css_dependency_audit.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2))
    print(f"\nfull report: {out.relative_to(ROOT)}")
    sys.exit(0)
