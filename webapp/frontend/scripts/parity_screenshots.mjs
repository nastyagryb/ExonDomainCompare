import { chromium } from "playwright";
import { mkdirSync } from "node:fs";
import { artifactPath } from "./repo_paths.mjs";

const OUT = artifactPath("pipeline_parity", "screenshots");
mkdirSync(OUT, { recursive: true });
const BASE = "http://localhost:5173";
const FGFR1 = "run:2026-07-11_1840_fgfr1_gallus_core_pilot";

const sleep = (ms) => new Promise((r) => setTimeout(r, ms));

async function shot(page, name) {
  await sleep(900);
  await page.screenshot({ path: `${OUT}/${name}.png`, fullPage: true });
  console.log("saved", name);
}

async function nav(page, text) {
  const btn = page.locator("header .nav button", { hasText: text }).first();
  await btn.click();
  await sleep(700);
}

async function selectDataset(page, value) {
  await page.selectOption("select.ds-select", value);
  await sleep(1200);
}

async function tab(page, text) {
  const t = page.locator(".workspace .tabs button.tab", { hasText: text }).first();
  if (await t.count()) { await t.click(); await sleep(700); }
  else console.log("tab not found:", text);
}

// Optional explicit browser binary (some sandboxes resolve the wrong arch path).
// Set PW_CHROME to a chrome-headless-shell / chrome binary; otherwise use default.
const launchOpts = { headless: true };
if (process.env.PW_CHROME) launchOpts.executablePath = process.env.PW_CHROME;
const browser = await chromium.launch(launchOpts);
const page = await browser.newPage({ viewport: { width: 1400, height: 900 } });
page.on("pageerror", (e) => console.log("PAGEERR", e.message));

await page.goto(BASE, { waitUntil: "networkidle" });
await sleep(1500);

// FGFR2
await selectDataset(page, "example");
await nav(page, "Overview");
await shot(page, "fgfr2_overview");
await nav(page, "Gene Explorer");
await shot(page, "fgfr2_gene_explorer");
await tab(page, "Summary");
await shot(page, "fgfr2_gene_summary");
await tab(page, "Evidence");
await shot(page, "fgfr2_gene_evidence");
await nav(page, "Figure Gallery");
await shot(page, "fgfr2_figures");
await nav(page, "Boundary");
await shot(page, "fgfr2_boundary");

// FGFR1
await selectDataset(page, FGFR1);
await nav(page, "Overview");
await shot(page, "fgfr1_overview");
await nav(page, "Gene Explorer");
await shot(page, "fgfr1_gene_explorer");
await tab(page, "Summary");
await shot(page, "fgfr1_gene_summary");
await tab(page, "Evidence");
await shot(page, "fgfr1_gene_evidence");
// Exploratory candidates are evidence rows in the SAME Evidence tab, not a
// separate page/controller.
await shot(page, "fgfr1_candidate_evidence");
await tab(page, "Synteny");
await shot(page, "fgfr1_synteny");
await tab(page, "Protein architecture");
await shot(page, "fgfr1_architecture");
await nav(page, "Exon–Domain Boundaries");
await shot(page, "fgfr1_boundary_pending");
await nav(page, "Figure Gallery");
await shot(page, "fgfr1_figures");

await browser.close();
console.log("DONE");
