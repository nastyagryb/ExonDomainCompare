import { chromium } from "playwright";
import { mkdirSync } from "node:fs";
import { artifactPath, frontendPath } from "./repo_paths.mjs";

const OUT = artifactPath("fgfr1_true_reuse", "screenshots");
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

async function pickSpecies(page, name) {
  const card = page.locator(".species-card, .sp-card", { hasText: name }).first();
  if (await card.count()) { await card.click(); await sleep(600); }
}

async function tab(page, text) {
  const t = page.locator(".workspace .tabs button.tab", { hasText: text }).first();
  if (await t.count()) { await t.click(); await sleep(700); }
  else console.log("tab not found:", text);
}

const launchOpts = { headless: true };
if (process.env.PW_CHROME) launchOpts.executablePath = process.env.PW_CHROME;
else {
  launchOpts.executablePath = frontendPath(
    ".pw-browsers/chromium_headless_shell-1228/chrome-headless-shell-mac-arm64/chrome-headless-shell",
  );
}
const browser = await chromium.launch(launchOpts);
const page = await browser.newPage({ viewport: { width: 1400, height: 900 } });
page.on("pageerror", (e) => console.log("PAGEERR", e.message));

await page.goto(BASE, { waitUntil: "networkidle" });
await sleep(1500);

const TABS = ["Summary", "Isoforms", "Evidence", "Files", "Exon map",
  "Domain architecture", "Boundary", "MSA", "Synteny"];

// FGFR2 Example → Gallus
await selectDataset(page, "example");
await nav(page, "Gene Explorer");
await pickSpecies(page, "Gallus");
for (const t of TABS) {
  await tab(page, t);
  await shot(page, `fgfr2_gallus_${t.toLowerCase().replace(/\s+/g, "_")}`);
}

// FGFR1 pilot → Gallus
await selectDataset(page, FGFR1);
await nav(page, "Gene Explorer");
for (const t of TABS) {
  await tab(page, t);
  await shot(page, `fgfr1_gallus_${t.toLowerCase().replace(/\s+/g, "_")}`);
}
await tab(page, "Exploratory Candidate Evidence");
await shot(page, "fgfr1_gallus_candidate_explorer");

await browser.close();
console.log("DONE");
