// The combined exon + domain architecture figure, checked as geometry.
//
// The figure exists to make one relationship readable: where an exon boundary
// falls relative to a domain edge, per species. That only works if both tracks
// of a species share one axis, if every species gets a group, and if nothing is
// drawn outside the plotting area. Those are geometric facts, so they are
// checked on the produced marks rather than on a screenshot.

import assert from "node:assert/strict";
import {
  comparativeExonDomainArchitectureFigureSpec,
} from "../webapp/frontend/src/pages/viewers/comparativeGalleryFigures.js";

const models = [
  {
    species_id: "mus_musculus", scientific_name: "Mus musculus",
    protein_id: "NP_1", protein_length: 800, status: "available",
    exons: [
      { label: "E1", start: 1, end: 100 },
      { label: "E2", start: 100, end: 300 },
      { label: "E3", start: 300, end: 800 },
    ],
    exon_boundaries: [
      { protein_position: 100, msa_column: 105, boundary_class: "exact_domain_edge" },
      { protein_position: 300, msa_column: 310, boundary_class: "inside_domain" },
    ],
    representative_domains: [
      { interpro_accession: "IPR001", short_label: "Ig-like 1", start: 40, end: 100 },
      { interpro_accession: "IPR002", short_label: "Kinase", start: 400, end: 700 },
    ],
    tm_regions: [{ start: 350, end: 372 }],
  },
  {
    species_id: "gallus_gallus", scientific_name: "Gallus gallus",
    protein_id: "NP_2", protein_length: 760, status: "available",
    exons: [
      { label: "E1", start: 1, end: 95 },
      { label: "E2", start: 95, end: 760 },
    ],
    exon_boundaries: [
      { protein_position: 95, msa_column: 105, boundary_class: "near_domain_edge" },
    ],
    representative_domains: [
      { interpro_accession: "IPR001", short_label: "Ig-like 1", start: 38, end: 96 },
    ],
    tm_regions: [],
  },
];

const alignedExons = models.flatMap((m, i) => (m.exons || []).map((e, j) => ({
  species_id: m.species_id, scientific_name: m.scientific_name,
  exon_label: e.label, native_start: e.start, native_end: e.end,
  msa_start_column: e.start + i * 5, msa_end_column: e.end + i * 5 + j,
})));
const alignedDomains = models.flatMap((m) => (m.representative_domains || []).map((d) => ({
  species_id: m.species_id, interpro_accession: d.interpro_accession,
  label: d.short_label, msa_start_column: d.start + 4, msa_end_column: d.end + 4,
  order_along_protein: 1,
})));

function build(mode) {
  return comparativeExonDomainArchitectureFigureSpec({
    gene: "FGFR1", models, domains: alignedDomains, exons: alignedExons,
    mode, nColumns: 820,
  });
}

const checks = [];
function check(name, fn) {
  try { fn(); checks.push(`ok   ${name}`); }
  catch (err) { checks.push(`FAIL ${name}: ${err.message}`); process.exitCode = 1; }
}

for (const mode of ["native", "msa"]) {
  const fig = build(mode);

  check(`${mode}: the figure is produced`, () => {
    assert.ok(fig, "no figure returned");
    assert.ok(fig.marks.length > 20, `only ${fig.marks.length} marks`);
  });

  const texts = fig.marks.filter((m) => m.t === "text").map((m) => String(m.s));

  check(`${mode}: every species gets its own group`, () => {
    for (const m of models) {
      const binomial = m.scientific_name;
      assert.ok(texts.some((t) => t.includes(binomial.split(" ")[0])),
        `no lane label for ${binomial}`);
      assert.ok(texts.includes(m.protein_id), `no protein id for ${binomial}`);
    }
  });

  check(`${mode}: each species reports its own exon and domain counts`, () => {
    assert.ok(texts.some((t) => t === "3 exons · 2 domains"), "mouse counts missing");
    assert.ok(texts.some((t) => t === "2 exons · 1 domains"), "chicken counts missing");
  });

  check(`${mode}: both tracks of a species share one axis`, () => {
    // A connector spans from the domain track to the exon track, so its
    // existence is what proves the two tracks are on one axis.
    const connectors = fig.marks.filter((m) => m.t === "line"
      && Math.abs(m.x1 - m.x2) < 0.01 && Math.abs(m.y2 - m.y1) > 20);
    assert.ok(connectors.length >= 3,
      `expected one connector per exon boundary, found ${connectors.length}`);
  });

  check(`${mode}: nothing is drawn outside the plotting area`, () => {
    const outside = fig.marks.filter((m) => {
      const xs = [m.x, m.x1, m.x2, m.cx].filter((v) => typeof v === "number");
      return xs.some((v) => v < -1 || v > fig.width + 1);
    });
    assert.equal(outside.length, 0,
      `${outside.length} mark(s) fall outside the ${fig.width}pt canvas`);
  });

  check(`${mode}: the axis is labelled for the coordinate system in use`, () => {
    const wanted = mode === "msa" ? "MSA column" : "Amino-acid position";
    assert.ok(texts.some((t) => t.includes(wanted)),
      `axis label does not mention ${wanted}`);
  });

  check(`${mode}: boundary classes are named in the legend`, () => {
    const legend = texts.join(" | ");
    assert.ok(legend.includes("Exon boundary"), "no exon-boundary legend entry");
    assert.ok(legend.includes("Coding exon"), "no coding-exon legend entry");
  });

  check(`${mode}: the caveat about the coordinate system is stated`, () => {
    const all = texts.join(" ");
    if (mode === "msa") {
      assert.ok(all.includes("aligned, not that they are equivalent"),
        "MSA panel does not disclaim equivalence");
    } else {
      assert.ok(all.includes("not comparable position by position"),
        "native panel does not disclaim cross-species comparison");
    }
  });
}

// A transmembrane helix is a native-coordinate fact; the MSA projection carries
// no aligned coordinates for it, so it must not be invented there.
check("the transmembrane helix appears only in the native panel", () => {
  const nativeLegend = build("native").marks
    .filter((m) => m.t === "text").map((m) => String(m.s)).join(" ");
  const msaLegend = build("msa").marks
    .filter((m) => m.t === "text").map((m) => String(m.s)).join(" ");
  assert.ok(nativeLegend.includes("Transmembrane helix"));
  assert.ok(!msaLegend.includes("Transmembrane helix"));
});

check("a species with no domain annotation still gets a group", () => {
  const fig = comparativeExonDomainArchitectureFigureSpec({
    gene: "TP53",
    models: [{
      species_id: "danio_rerio", scientific_name: "Danio rerio", protein_id: "NP_3",
      protein_length: 373, status: "pending", exons: [{ label: "E1", start: 1, end: 373 }],
      exon_boundaries: [], representative_domains: [],
    }],
    domains: [], exons: [], mode: "native",
  });
  const texts = fig.marks.filter((m) => m.t === "text").map((m) => String(m.s));
  assert.ok(texts.some((t) => t.includes("domain annotation pending")),
    "a pending annotation is not stated");
  assert.ok(texts.some((t) => t === "1 exons · 0 domains"));
});

check("an empty dataset yields no figure rather than an empty one", () => {
  assert.equal(comparativeExonDomainArchitectureFigureSpec({ gene: "X", models: [] }), null);
});

// The species order is the canonical one: mammal before bird, not alphabetical.
check("species appear in the canonical taxonomic order", () => {
  const fig = build("native");
  const labels = fig.marks
    .filter((m) => m.t === "text" && /^(Mus|Gallus)/.test(String(m.s)))
    .sort((a, b) => a.y - b.y)
    .map((m) => String(m.s).split(" ")[0]);
  assert.deepEqual(labels, ["Mus", "Gallus"],
    "alphabetical order would place Gallus first");
});

console.log(checks.join("\n"));
