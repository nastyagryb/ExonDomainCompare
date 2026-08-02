// Does the interactive view really paint a feature in the shared colour?
//
// Reading the source only proves that an attribute is written. A CSS class rule
// overrides a presentation attribute, so the question can only be settled in a live
// browser, against the *computed* style. This script navigates the Gene Explorer of a
// single-species run and compares the computed fill of each scientific mark with the
// value the shared visual specification prescribes.
//
// Run against a served frontend (vite dev on 5173 or preview on 4173):
//     node scripts/verify_interactive_palette.mjs [runId] [baseUrl]

import { chromium } from "playwright";
import { mkdirSync, writeFileSync } from "node:fs";
import { artifactPath } from "./repo_paths.mjs";
import {
  FEATURE_STYLES, TEXT_ROLES, domainInstanceFill,
} from "../src/pages/viewers/semanticStyles.js";

const RUN = process.argv[2] || "2026-07-23_1100_fgfr1_gallus_core_pilot";
const BASE = process.argv[3] || "http://localhost:5173";
const OUT = artifactPath("single_species_plot_repair", "palette_verification");
mkdirSync(OUT, { recursive: true });

const rgb = (hex) => {
  const h = hex.replace("#", "");
  const n = parseInt(h.length === 3 ? h.split("").map((c) => c + c).join("") : h, 16);
  return `rgb(${(n >> 16) & 255}, ${(n >> 8) & 255}, ${n & 255})`;
};

/** selector -> the colour the shared specification prescribes for that mark. */
const EXPECT = [
  { view: "Exon Map", sel: ".exon-map .em-exon", prop: "fill",
    want: FEATURE_STYLES.coding_exon.fill, what: "coding exon block" },
  { view: "Exon Map", sel: ".exon-map .em-exon-num", prop: "fill",
    want: TEXT_ROLES.featureLabel.fill, what: "exon number on the block" },
  { view: "Exon Map", sel: ".exon-map .em-axis-lbl", prop: "fill",
    want: TEXT_ROLES.axis.fill, what: "axis tick label" },
  { view: "Exon Map", sel: ".exon-map .em-axis-end", prop: "fill",
    want: TEXT_ROLES.axisEmphasis.fill, what: "axis end label" },
  { view: "Exon Map", sel: ".exon-map .em-cand-lbl", prop: "fill",
    want: TEXT_ROLES.candidateLabel.fill, what: "candidate label" },
  { view: "Transcript comparison", sel: ".cmp-track .cmp-blk-shared", prop: "fill",
    want: FEATURE_STYLES.shared_exon.fill, what: "shared exon" },
  { view: "Transcript comparison", sel: ".cmp-track .cmp-blk-alt", prop: "fill",
    want: FEATURE_STYLES.alternative_exon.fill, what: "alternative exon" },
  { view: "Transcript comparison", sel: ".cmp-ruler .cmp-tick", prop: "fill",
    want: TEXT_ROLES.axis.fill, what: "ruler tick" },
  { view: "Domain architecture", sel: ".domain-arch .da-lane-lbl", prop: "fill",
    want: TEXT_ROLES.trackLabel.fill, what: "track label" },
  { view: "Domain architecture", sel: ".domain-arch .da-blk-lbl", prop: "fill",
    want: TEXT_ROLES.onFeatureLabel.fill, what: "label on a domain block" },
  { view: "Boundary explorer", sel: ".boundary-explorer .da-lane-lbl", prop: "fill",
    want: TEXT_ROLES.trackLabel.fill, what: "track label" },
  { view: "Boundary explorer", sel: ".boundary-explorer .em-exon", prop: "fill",
    want: FEATURE_STYLES.coding_exon.fill, what: "coding exon block" },
];

const browser = await chromium.launch();
const page = await browser.newPage({ viewport: { width: 1680, height: 1050 } });
const report = { run_id: RUN, base: BASE, checks: [], screenshots: [], notes: [] };

const shot = async (name) => {
  const file = `${OUT}/${name}.png`;
  await page.screenshot({ path: file, fullPage: true });
  report.screenshots.push(file);
};

await page.goto(BASE, { waitUntil: "networkidle" });
await page.selectOption("select.ds-select", `run:${RUN}`);
await page.waitForTimeout(1200);
await page.locator("header .nav button", { hasText: "Gene Explorer" }).first().click();
await page.waitForTimeout(2000);
await shot("01_gene_explorer");

// Walk every tab of the explorer so each viewer gets mounted at least once.
const tabs = await page.locator(".ge-tabs button, .tabs button, .subnav button")
  .allTextContents();
report.notes.push(`tabs found: ${tabs.join(" | ") || "none"}`);

// Marks only exist while their tab is mounted, so every fill has to be harvested as we
// pass through the tabs rather than once at the end.
const seenFills = new Set();

const measure = async () => {
  for (const f of await page.locator("svg rect")
    .evaluateAll((els) => els.map((el) => getComputedStyle(el).fill))) seenFills.add(f);
  for (const e of EXPECT) {
    if (report.checks.some((c) => c.selector === e.sel && c.found)) continue;
    const loc = page.locator(e.sel).first();
    const n = await page.locator(e.sel).count();
    if (n === 0) continue;
    const got = await loc.evaluate((el, prop) => getComputedStyle(el)[prop], e.prop);
    const want = rgb(e.want);
    report.checks.push({
      view: e.view, what: e.what, selector: e.sel, property: e.prop,
      expected: `${e.want} = ${want}`, computed: got, ok: got === want,
      instances: n, found: true,
    });
  }
};

await measure();
for (const [i, label] of tabs.entries()) {
  try {
    await page.locator(".ge-tabs button, .tabs button, .subnav button")
      .filter({ hasText: label }).first().click({ timeout: 4000 });
    await page.waitForTimeout(1400);
    await shot(`02_tab_${String(i + 1).padStart(2, "0")}_${label.replace(/\W+/g, "_")}`);
    await measure();
  } catch (err) {
    report.notes.push(`tab "${label}" not reachable: ${String(err).slice(0, 120)}`);
  }
}

// Anything the run does not render cannot be judged, and must be reported as such
// rather than silently counted as a pass.
for (const e of EXPECT) {
  if (!report.checks.some((c) => c.selector === e.sel)) {
    report.checks.push({ view: e.view, what: e.what, selector: e.sel,
      found: false, ok: null, note: "mark not present in this run's views" });
  }
}

// The four repeated domain instances must be four *distinguishable* colours, otherwise
// the repeated-Ig-domain repair is invisible in the interactive view.
report.domain_ramp_expected = [1, 2, 3, 4].map((i) => rgb(domainInstanceFill(i)));
report.domain_ramp_present = report.domain_ramp_expected.filter((c) => seenFills.has(c));

await browser.close();

const judged = report.checks.filter((c) => c.found);
const bad = judged.filter((c) => !c.ok);
report.summary = {
  judged: judged.length, matching: judged.length - bad.length, mismatching: bad.length,
  not_rendered: report.checks.length - judged.length,
};
writeFileSync(`${OUT}/palette_verification.json`, JSON.stringify(report, null, 2));

console.log(`run ${RUN} @ ${BASE}`);
for (const c of report.checks) {
  if (!c.found) { console.log(`  --  ${c.view} · ${c.what}: not rendered`); continue; }
  console.log(`  ${c.ok ? "ok" : "MISMATCH"}  ${c.view} · ${c.what} (${c.instances}x)`);
  if (!c.ok) console.log(`        expected ${c.expected}  got ${c.computed}`);
}
console.log(`\ndomain ramp: ${report.domain_ramp_present.length}/4 instance colours found`);
console.log(JSON.stringify(report.summary));
if (bad.length) process.exitCode = 1;
