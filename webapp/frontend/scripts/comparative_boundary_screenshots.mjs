// Screenshots of the Comparative Exon–Domain Boundary Explorer on the real
// two-species FGFR1 (Gallus + Mus) run.
//
// Every shot comes from the live frontend against the live API — nothing is staged and
// no fixture is substituted, so a broken view produces a broken screenshot instead of a
// reassuring one. The script fails loudly if the page reports an error or if an
// expected element is missing.
//
//   node scripts/comparative_boundary_screenshots.mjs

import { chromium } from "playwright";
import { mkdirSync } from "node:fs";
import { artifactPath } from "./repo_paths.mjs";

const OUT = artifactPath("multispecies_scientific_phase", "screenshots");
const BASE = process.env.BASE || "http://localhost:5173";
const DATASET = "run:2026-07-26_2157_fgfr1_gallus_mus_core_pilot";
mkdirSync(OUT, { recursive: true });

const browser = await chromium.launch({ headless: true });
const page = await browser.newPage({ viewport: { width: 1500, height: 1000 } });
const errors = [];
page.on("pageerror", (e) => errors.push(`page: ${e.message}`));
page.on("response", (r) => {
  if (r.status() >= 400 && r.url().includes("/api/")) {
    errors.push(`${r.status()} ${r.url()}`);
  }
});

const shot = async (name, locator = null) => {
  if (locator) {
    await locator.screenshot({ path: `${OUT}/${name}.png` });
  } else {
    await page.screenshot({ path: `${OUT}/${name}.png`, fullPage: true });
  }
  console.log("saved", name);
};
const nav = async (label) => {
  await page.locator("header .nav button", { hasText: label }).first().click();
  await page.waitForTimeout(600);
};
const need = async (locator, what) => {
  if (!await locator.count()) throw new Error(`missing: ${what}`);
  return locator;
};

await page.goto(BASE, { waitUntil: "networkidle" });
await page.selectOption("select.ds-select", DATASET);
await page.waitForTimeout(1000);

// 1. run status: both species complete
await shot("01_run_status_gallus_mus_complete");

await nav("Exon–Domain Boundaries");
await page.waitForTimeout(700);

const explorer = await need(page.locator(".cbe-explorer"),
  "the comparative explorer did not render");
const matrixCard = page.locator(".cbe-explorer .card").filter(
  { hasText: "Species × comparable-boundary groups" }).first();
await need(matrixCard, "the comparative matrix card");

// 2. the comparative matrix
await shot("02_comparative_matrix", matrixCard);

// 3. matrix hover: the full per-cell scientific record
const cells = page.locator(".cbe-matrix .bnd-heat-cell");
const nCells = await cells.count();
if (nCells < 16) throw new Error(`only ${nCells} matrix cells rendered`);
await cells.nth(2).hover();
await page.waitForTimeout(350);
await need(page.locator(".cbe-hover-card"), "the hover card");
await shot("03_matrix_hover", matrixCard);

// 4. a selected matrix cell drives every linked view
await cells.nth(2).click();
await page.waitForTimeout(500);
await shot("04_matrix_cell_selected", matrixCard);

// 5. the paired signed-distance plot
const plotCard = page.locator(".cbe-explorer .card").filter(
  { hasText: "Signed distance to nearest domain edge, per species" }).first();
await need(plotCard, "the paired plot card");
await shot("05_paired_signed_distance_plot", plotCard);

// 6. the selected comparable group: detail rows + comparative architecture
const detailCard = page.locator(".cbe-explorer .card").filter(
  { hasText: "Selected comparable boundary" }).first();
await need(detailCard, "the selected-group detail card");
await shot("06_selected_comparable_group", detailCard);

// 9. comparative local architecture on its own
await need(page.locator(".cbe-arch-svg"), "the comparative architecture");
await shot("09_comparative_architecture_detail", page.locator(".cbe-arch").first());

// 7. an active filter state
const filterCard = page.locator(".cbe-explorer .card").filter({ hasText: "Filters" }).first();
await need(filterCard, "the filter card");
await page.locator(".cbe-chip", { hasText: "Inside domain" }).first().click();
await page.waitForTimeout(400);
await page.locator(".cbe-chip", { hasText: "Exact / near edge only" }).first().click();
await page.waitForTimeout(500);
await shot("07_filter_state", explorer);

// 8. reset restores the full dataset
await page.locator("button", { hasText: "Reset all filters" }).first().click();
await page.waitForTimeout(500);
await shot("08_filters_reset", filterCard);

// 10. an inspection case, clicked, updating the linked views
const caseCard = page.locator(".cbe-explorer .card").filter(
  { hasText: "Inspection cases" }).first();
await need(caseCard, "the inspection cases card");
await page.locator(".cbe-case-btn").first().click();
await page.waitForTimeout(600);
await shot("10_inspection_case", caseCard);
await shot("10b_inspection_case_linked_detail", detailCard);

// 11. the export menu
await matrixCard.locator("button.menu-btn", { hasText: "Export" }).first().click();
await page.waitForTimeout(350);
await need(page.locator(".ui-menu.open .ui-menu-panel"), "the export menu");
await shot("11_export_menu", matrixCard);

// 12. an actually downloaded matrix PDF, proving the export path works end to end
const download = await Promise.all([
  page.waitForEvent("download", { timeout: 15000 }),
  page.locator(".ui-menu.open .menu-item", { hasText: "Matrix — PDF" }).first().click(),
]).then(([d]) => d);
const pdfPath = `${OUT}/12_exported_matrix.pdf`;
await download.saveAs(pdfPath);
console.log("saved 12_exported_matrix.pdf as", download.suggestedFilename());

// Full page for context.
await page.keyboard.press("Escape");
await page.waitForTimeout(300);
await shot("00_comparative_explorer_full_page");

await browser.close();
if (errors.length) {
  console.error(`\n${errors.length} page/API error(s):`);
  for (const e of errors) console.error(`  ${e}`);
  process.exit(1);
}
console.log("\nOK — comparative boundary explorer screenshots complete");
