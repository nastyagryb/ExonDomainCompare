import { chromium } from "playwright";
import { mkdirSync } from "node:fs";
import { artifactPath } from "./repo_paths.mjs";

const OUT = artifactPath("linked_explorer", "screenshots");
const BASE = "http://localhost:5173";
const DATASET = "run:2026-07-11_1840_fgfr1_gallus_core_pilot";
mkdirSync(OUT, { recursive: true });

const EXEC = process.env.PW_CHROME
  || "/var/folders/bg/yg2g8h9d5m349g98hm90y88h0000gn/T/cursor-sandbox-cache/efee57be99066b7de11bb33275368c15/playwright/chromium_headless_shell-1228/chrome-headless-shell-mac-arm64/chrome-headless-shell";
const browser = await chromium.launch({ headless: true, executablePath: EXEC });
const page = await browser.newPage({ viewport: { width: 1500, height: 1000 } });
const errors = [];
page.on("pageerror", (error) => errors.push(`page: ${error.message}`));
page.on("response", (r) => {
  if (r.status() >= 400 && r.url().includes("/api/")) errors.push(`${r.status()} ${r.url()}`);
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
  await page.waitForTimeout(600);
};
const shot = async (name) => {
  await page.screenshot({ path: `${OUT}/${name}.png`, fullPage: true });
  console.log("saved", name);
};

await nav("Gene Explorer");
await page.waitForTimeout(500);

// 0. Full base-module tab bar (same FGFR2 module order) with default Summary.
await shot("00_base_modules_tabbar");

// 1. Candidate ranking with C1 selected by default (top-ranked).
await tab("Exploratory Candidate Evidence");
await shot("01_candidate_ranking_c1_selected");

// 2. Linked transcript / exon highlight (C1 stays selected across tabs).
await tab("Transcript & Exon Structure");
await shot("02_linked_transcript_exon");

// 3. Linked Exon map highlight (shared FGFR2 CoordinateTrack module).
await tab("Exon map");
await shot("03_linked_exon_map");

// 4. Linked isoform-alignment highlight (auto-scrolls to C1 band).
await tab("Isoform Alignment");
await shot("04_linked_isoform_alignment");

// 5. Candidate-specific evidence.
await tab("Evidence");
await shot("05_candidate_specific_evidence");

// 6. Provenance hub.
await tab("Files");
await shot("06_provenance_hub");

// 7. Synteny base module (single-species local neighbourhood).
await tab("Synteny");
await shot("07_synteny_local_neighbourhood");

// 8. Figure gallery (pre-cluster).
await tab("Figure Gallery");
await shot("08_figure_gallery");

// 9. Domain architecture base module — honest pending state (post-cluster readiness).
await tab("Domain architecture");
await shot("09_domain_pending");

await browser.close();
if (errors.length) throw new Error("ERRORS:\n" + errors.join("\n"));
console.log("DONE — FGFR1 linked workspace screenshots captured with no page/API errors");
