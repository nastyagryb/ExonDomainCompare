/**
 * Screenshots for Part 3 (Comparative Figure Gallery) and Part 4 (Data & Downloads).
 *
 * Usage (API + Vite must be running; Playwright Chromium installed):
 *   node scripts/comparative_gallery_downloads_screenshots.mjs
 */
import { chromium } from "@playwright/test";
import { mkdirSync } from "node:fs";
import { dirname, join } from "node:path";
import { fileURLToPath } from "node:url";

const __dirname = dirname(fileURLToPath(import.meta.url));
const OUT = join(__dirname, "../../../artifacts/multispecies_scientific_phase/screenshots");
const BASE = process.env.EDC_UI || "http://localhost:5173";
const RUN = "run:2026-07-26_2157_fgfr1_gallus_mus_core_pilot";

mkdirSync(OUT, { recursive: true });

async function shot(page, name) {
  const path = join(OUT, `${name}.png`);
  await page.screenshot({ path, fullPage: false });
  console.log("saved", name);
}

const browser = await chromium.launch({ headless: true });
const page = await browser.newPage({ viewport: { width: 1400, height: 900 } });

try {
  await page.goto(`${BASE}/?dataset=${encodeURIComponent(RUN)}`, { waitUntil: "networkidle" });
  // Select dataset via UI if needed
  await page.waitForTimeout(800);

  // Figure Gallery — comparative scope
  await page.getByRole("button", { name: /Figure Gallery/i }).click();
  await page.waitForTimeout(1000);
  await shot(page, "30_gallery_comparative_scope");

  const scope = page.locator('select[aria-label="Figure gallery scope"]');
  if (await scope.count()) {
    await scope.selectOption({ label: /Gallus/i });
    await page.waitForTimeout(600);
    await shot(page, "31_gallery_gallus_scope");
    await scope.selectOption({ label: /Mus/i });
    await page.waitForTimeout(600);
    await shot(page, "32_gallery_mus_scope");
    await scope.selectOption({ value: "comparative" });
    await page.waitForTimeout(600);
    await shot(page, "33_gallery_comparative_reset");
  }

  // Gene Explorer → Data & Downloads
  await page.getByRole("button", { name: /Gene Explorer/i }).click();
  await page.waitForTimeout(800);
  const dd = page.getByRole("button", { name: /Data & Downloads/i });
  if (await dd.count()) {
    await dd.click();
    await page.waitForTimeout(1000);
    await shot(page, "40_data_downloads_recommended");
    const custom = page.getByRole("button", { name: /^Custom$/i });
    if (await custom.count()) {
      await custom.click();
      await page.waitForTimeout(400);
      await shot(page, "41_data_downloads_custom");
    }
  }

  console.log("OK — gallery/downloads screenshots complete");
} catch (err) {
  console.error(err);
  process.exitCode = 1;
} finally {
  await browser.close();
}
