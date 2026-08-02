import { chromium } from "playwright";
import { mkdirSync } from "node:fs";
import { artifactPath } from "./repo_paths.mjs";

const OUT = artifactPath("fgfr1_repair", "screenshots");
const BASE = "http://localhost:5173";
mkdirSync(OUT, { recursive: true });
const EXEC = process.env.PW_CHROME
  || "/var/folders/bg/yg2g8h9d5m349g98hm90y88h0000gn/T/cursor-sandbox-cache/efee57be99066b7de11bb33275368c15/playwright/chromium_headless_shell-1228/chrome-headless-shell-mac-arm64/chrome-headless-shell";

const browser = await chromium.launch({ headless: true, executablePath: EXEC });
const page = await browser.newPage({ viewport: { width: 1500, height: 950 } });
const errors = [];
page.on("pageerror", (e) => errors.push(`page: ${e.message}`));
page.on("response", (r) => { if (r.status() >= 400 && r.url().includes("/api/runs/")) errors.push(`${r.status()} ${r.url()}`); });

const nav = async (label) => { await page.locator("header .nav button", { hasText: label }).first().click(); await page.waitForTimeout(600); };
const tab = async (label) => {
  const b = page.locator(".workspace .tabs button.tab", { hasText: label }).first();
  if (!await b.count()) throw new Error(`Missing tab: ${label}`);
  await b.click(); await page.waitForTimeout(500);
};
const shot = async (name) => { await page.screenshot({ path: `${OUT}/${name}.png`, fullPage: true }); console.log("saved", name); };

await page.goto(BASE, { waitUntil: "networkidle" });

// FGFR2 example — immutable reference. Capture Gene Explorer, Files and one event tab.
await page.selectOption("select.ds-select", "example");
await page.waitForTimeout(900);
await nav("Gene Explorer");
await shot("ref_fgfr2_gene_explorer");
await tab("Files");
await shot("ref_fgfr2_files");

await browser.close();
if (errors.length) throw new Error("ERRORS:\n" + errors.join("\n"));
console.log("DONE — FGFR2 reference screenshots captured with no page/API errors");
