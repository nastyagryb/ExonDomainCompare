import { chromium } from "playwright";
import { mkdirSync } from "node:fs";
import { artifactPath } from "./repo_paths.mjs";

const OUT = artifactPath("scientific_explorer", "screenshots");
const BASE = "http://localhost:5173";
const DATASET = "run:2026-07-11_1840_fgfr1_gallus_core_pilot";
mkdirSync(OUT, { recursive: true });

const browser = await chromium.launch({ headless: true });
const page = await browser.newPage({ viewport: { width: 1500, height: 950 } });
const errors = [];
page.on("pageerror", (error) => errors.push(`page: ${error.message}`));
page.on("response", (response) => {
  if (response.status() >= 400 && response.url().includes("/api/runs/")) {
    errors.push(`${response.status()} ${response.url()}`);
  }
});
await page.goto(BASE, { waitUntil: "networkidle" });
await page.selectOption("select.ds-select", DATASET);
await page.waitForTimeout(800);

const nav = async (label) => {
  await page.locator("header .nav button", { hasText: label }).first().click();
  await page.waitForTimeout(500);
};
const tab = async (label) => {
  const button = page.locator(".workspace .tabs button.tab", { hasText: label }).first();
  if (!await button.count()) throw new Error(`Missing tab: ${label}`);
  await button.click();
  await page.waitForTimeout(450);
};
const shot = async (name) => {
  await page.screenshot({ path: `${OUT}/${name}.png`, fullPage: true });
  console.log("saved", name);
};

await nav("Gene Explorer");
for (const [label, name] of [
  ["Summary", "01_summary"],
  ["Protein Models", "02_protein_models"],
  ["Transcript & Exon Structure", "03_transcript_exon_structure"],
  ["Protein Architecture", "04_protein_architecture"],
  ["Isoform Alignment", "05_isoform_alignment"],
  ["Exploratory Candidate Evidence", "06_candidate_evidence"],
  ["Local Gene Neighbourhood", "07_local_gene_neighbourhood"],
  ["Domain Architecture", "08_domain_pending"],
]) {
  await tab(label);
  await shot(name);
}
await nav("Exon–Domain");
await shot("09_exon_domain_pending");
await nav("Figure Gallery");
await page.waitForTimeout(1000);
const previewStatus = await page.locator(".fig-card img").evaluateAll((images) =>
  images.map((img) => ({ src: img.src, loaded: img.complete && img.naturalWidth > 0 })));
if (!previewStatus.length || previewStatus.some((item) => !item.loaded)) {
  errors.push(`Figure preview failure: ${JSON.stringify(previewStatus)}`);
}
await shot("10_figure_gallery");
await browser.close();
if (errors.length) throw new Error(errors.join("\n"));
console.log("DONE — all scientific explorer screenshots and previews verified");
