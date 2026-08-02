// Figure Gallery screenshots + live contract check for a single-species run.
//
// Captures the Gallery overview, every scientific category, the supplement view and
// the three alignment/supplement cards, and asserts the things a screenshot alone
// cannot show: that the six categories survive frontend normalisation, that the
// member-database signature card is treated as a supplement rather than a main
// figure, and that every card carries a question and its four downloads.
//
// Usage: node scripts/gallery_single_species_screenshots.mjs [runId] [outDir]

import { chromium } from "playwright";
import { mkdirSync, writeFileSync } from "node:fs";
import { artifactPath } from "./repo_paths.mjs";

const RUN = process.argv[2] || "2026-07-23_1100_fgfr1_gallus_core_pilot";
const OUT = process.argv[3]
  || artifactPath("single_species_plot_repair", "after", "screenshots");
const BASE = "http://localhost:5173";

const CATEGORIES = [
  ["Exon structure", "09_gallery_exon_structure", 4],
  ["Isoform analysis", "10_gallery_isoform_analysis", 3],
  ["Domain architecture", "11_gallery_domain_architecture", 1],
  ["Exon–domain boundaries", "12_gallery_boundaries", 3],
  ["Genomic context", "13_gallery_genomic_context", 1],
  ["Exploratory candidates", "14_gallery_exploratory_candidates", 2],
];
const SUPPLEMENT_TITLE = "Member-database signature";

mkdirSync(OUT, { recursive: true });
const problems = [];
const report = { run_id: RUN, screenshots: [], checks: {} };

const browser = await chromium.launch({ headless: true });
const page = await browser.newPage({ viewport: { width: 1500, height: 980 } });
page.on("pageerror", (e) => problems.push(`page error: ${e.message}`));
page.on("response", (r) => {
  if (r.status() >= 400 && r.url().includes("/api/runs/")) problems.push(`${r.status()} ${r.url()}`);
});

const shot = async (name) => {
  await page.screenshot({ path: `${OUT}/${name}.png`, fullPage: true });
  report.screenshots.push(`${name}.png`);
  console.log("  saved", name);
};
const cardTitles = () => page.locator(".fig-card .fig-title-row b").allTextContents();
const categorySelect = () => page.locator(".page-head .filters select").last();
const supplementBox = () => page.locator(".page-head .filters label.check input[type=checkbox]");

await page.goto(BASE, { waitUntil: "networkidle" });
await page.selectOption("select.ds-select", `run:${RUN}`);
await page.waitForTimeout(900);
await page.locator("header .nav button", { hasText: "Figure Gallery" }).first().click();
await page.waitForTimeout(1200);

// ---- overview, with supplements hidden (the default) -----------------------
const groupsShown = async () => {
  const heads = await page.locator(".fig-group .group-title").allTextContents();
  return heads.map((h) => h.replace(/\s*·\s*\d+\s*$/, "").trim());
};
const nCards = () => page.locator(".fig-card").count();

report.checks.default_categories = await groupsShown();
report.checks.default_card_count = await nCards();
const defaultTitles = await cardTitles();
report.checks.supplement_hidden_by_default =
  !defaultTitles.some((t) => t.includes(SUPPLEMENT_TITLE));
await shot("08_gallery_full");

for (const cat of CATEGORIES.map((c) => c[0])) {
  if (!report.checks.default_categories.includes(cat)) problems.push(`category missing: ${cat}`);
}
if (!report.checks.supplement_hidden_by_default) {
  problems.push("the member-database signature card is shown as a main figure");
}

// ---- previews and per-card completeness ------------------------------------
const previews = await page.locator(".fig-card img").evaluateAll((imgs) =>
  imgs.map((i) => ({ src: i.src, ok: i.complete && i.naturalWidth > 0 })));
report.checks.previews_total = previews.length;
report.checks.previews_broken = previews.filter((p) => !p.ok).map((p) => p.src);
if (!previews.length || previews.some((p) => !p.ok)) problems.push("broken figure previews");
if (await page.locator(".fig-noimg").count()) problems.push("a card has no preview image");

const perCard = await page.locator(".fig-card").evaluateAll((cards) => cards.map((c) => {
  const t = (c.querySelector(".fig-title-row b") || {}).textContent || "";
  const links = [...c.querySelectorAll("a.btn")].map((a) => a.textContent.trim());
  return { title: t.trim(), question: Boolean(c.querySelector(".fig-q")), links };
}));
report.checks.cards = perCard;
for (const c of perCard) {
  if (!c.title) problems.push("a card has an empty title");
  if (/\.(svg|pdf|png|tsv)$/i.test(c.title)) problems.push(`filename used as title: ${c.title}`);
  if (!c.question) problems.push(`no scientific question: ${c.title}`);
  for (const f of ["SVG", "PDF", "PNG"]) {
    if (!c.links.includes(f)) problems.push(`missing ${f} download: ${c.title}`);
  }
}
const seen = new Set();
for (const c of perCard) {
  if (seen.has(c.title)) problems.push(`duplicate card: ${c.title}`);
  seen.add(c.title);
}

// ---- one screenshot per scientific category --------------------------------
for (const [cat, name, expected] of CATEGORIES) {
  await categorySelect().selectOption(cat);
  await page.waitForTimeout(500);
  const n = await nCards();
  report.checks[`count_${name}`] = n;
  if (n !== expected) problems.push(`${cat}: expected ${expected} main card(s), found ${n}`);
  await shot(name);
}
await categorySelect().selectOption("all");
await page.waitForTimeout(400);

// ---- supplements revealed --------------------------------------------------
await supplementBox().check();
await page.waitForTimeout(600);
report.checks.card_count_with_supplements = await nCards();
const withSupp = await cardTitles();
report.checks.supplement_visible_when_enabled = withSupp.some((t) => t.includes(SUPPLEMENT_TITLE));
if (!report.checks.supplement_visible_when_enabled) {
  problems.push("the supplement never appears, even with 'Show supplements' enabled");
}
if (report.checks.card_count_with_supplements <= report.checks.default_card_count) {
  problems.push("enabling supplements did not reveal an additional card");
}
await shot("15_gallery_supplements_shown");

// ---- the individual cards the report references ----------------------------
const openCard = async (needle, name) => {
  const card = page.locator(".fig-card", { hasText: needle }).first();
  if (!await card.count()) { problems.push(`card not found: ${needle}`); return; }
  await card.locator("button", { hasText: "Open" }).first().click();
  await page.waitForTimeout(700);
  const subtitle = (await page.locator(".drawer .drawer-head p").allTextContents()).join(" | ");
  report.checks[`drawer_${name}`] = subtitle;
  await shot(name);
  await page.keyboard.press("Escape");
  await page.waitForTimeout(400);
};
await openCard("Isoform alignment overview", "16_card_full_isoform_alignment");
await openCard("Candidate-associated alignment", "17_card_candidate_alignment_detail");
await openCard(SUPPLEMENT_TITLE, "18_card_member_signature_supplement");

await browser.close();
writeFileSync(`${OUT}/gallery_verification.json`, `${JSON.stringify(report, null, 2)}\n`);
console.log("\n--- checks ---");
console.log(JSON.stringify(report.checks, null, 2).slice(0, 2000));
if (problems.length) {
  console.error(`\nFAILED (${problems.length}):`);
  for (const p of problems) console.error(`  - ${p}`);
  process.exit(1);
}
console.log(`\nDONE — ${report.screenshots.length} screenshots, all Gallery checks passed`);
