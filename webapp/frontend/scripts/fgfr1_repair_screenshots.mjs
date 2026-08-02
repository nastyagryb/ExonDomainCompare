import { chromium } from "playwright";
import { mkdirSync } from "node:fs";
import { artifactPath } from "./repo_paths.mjs";

const OUT = artifactPath("fgfr1_repair", "screenshots");
const BASE = "http://localhost:5173";
const DATASET = "run:2026-07-11_1840_fgfr1_gallus_core_pilot";
mkdirSync(OUT, { recursive: true });

const EXEC = process.env.PW_CHROME
  || "/var/folders/bg/yg2g8h9d5m349g98hm90y88h0000gn/T/cursor-sandbox-cache/efee57be99066b7de11bb33275368c15/playwright/chromium_headless_shell-1228/chrome-headless-shell-mac-arm64/chrome-headless-shell";
const browser = await chromium.launch({ headless: true, executablePath: EXEC });
const page = await browser.newPage({ viewport: { width: 1500, height: 950 } });
const errors = [];
page.on("pageerror", (error) => errors.push(`page: ${error.message}`));
page.on("response", (r) => {
  if (r.status() >= 400 && r.url().includes("/api/runs/")) errors.push(`${r.status()} ${r.url()}`);
});

await page.goto(BASE, { waitUntil: "networkidle" });
await page.selectOption("select.ds-select", DATASET);
await page.waitForTimeout(900);

const nav = async (label) => {
  await page.locator("header .nav button", { hasText: label }).first().click();
  await page.waitForTimeout(600);
};
const tab = async (label) => {
  const b = page.locator(".workspace .tabs button.tab", { hasText: label }).first();
  if (!await b.count()) throw new Error(`Missing tab: ${label}`);
  await b.click();
  await page.waitForTimeout(500);
};
const shot = async (name) => {
  await page.screenshot({ path: `${OUT}/${name}.png`, fullPage: true });
  console.log("saved", name);
};

// Overview
await nav("Overview");
await shot("00_overview");

// Gene Explorer — species panel is the sidebar; Summary is default.
await nav("Gene Explorer");
await shot("01_species_panel_and_summary");
for (const [label, name] of [
  ["Summary", "02_summary"],
  ["Protein Models", "03_protein_models"],
  ["Files", "04_files"],
  ["Transcript & Exon Structure", "05_transcript_exon_structure"],
  ["Protein Architecture", "06_protein_architecture"],
  ["Isoform Alignment", "07_isoform_alignment"],
  ["Exploratory Candidate Evidence", "08_candidate_evidence"],
  ["Local Gene Neighbourhood", "09_local_gene_neighbourhood"],
  ["Domain Architecture", "10_domain_pending"],
  ["Exon–Domain Analysis", "11_exon_domain_pending"],
]) {
  await tab(label);
  await shot(name);
}

await browser.close();
if (errors.length) throw new Error("ERRORS:\n" + errors.join("\n"));
console.log("DONE — FGFR1 repair screenshots captured with no page/API errors");
