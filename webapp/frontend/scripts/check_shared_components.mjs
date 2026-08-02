// Static regression guard for FGFR2 ↔ FGFR1 code-level UI parity.
//
// Fails (exit 1) if the FGFR2 pages or the FGFR1 renderer stop importing the
// shared components, or if a parallel Core-only UI is reintroduced as an active
// route. This is a cheap smoke test — no browser required.
//
//   node scripts/check_shared_components.mjs
import { existsSync, readFileSync } from "node:fs";

const ROOT = new URL("../src/", import.meta.url);
const read = (rel) => readFileSync(new URL(rel, ROOT), "utf8");

const checks = [];
const must = (name, ok) => checks.push({ name, ok });

const overview = read("pages/Overview.jsx");
const geneExplorer = read("pages/GeneExplorer.jsx");
const figureGallery = read("pages/FigureGallery.jsx");
const boundary = read("pages/BoundaryPage.jsx");
const app = read("App.jsx");
const api = read("api.js");

const importsShared = (src) => /from ["']\.\.\/components\/shared["']/.test(src)
  || /from ["']\.\.\/\.\.\/components\/shared["']/.test(src);

// The same four controllers are mounted for every dataset.
must("Overview is the only Overview controller",
  importsShared(overview) && /DatasetPageHeader/.test(overview) && /KpiGrid/.test(overview));
must("GeneExplorer is species-centric and uses the one shared shell",
  importsShared(geneExplorer) && /GeneExplorerShell/.test(geneExplorer)
  && /SpeciesPanel/.test(geneExplorer) && /SpeciesCard/.test(geneExplorer)
  && /ExplorerTabs/.test(geneExplorer));
must("FigureGallery is the only gallery controller",
  importsShared(figureGallery) && /FigureCardGrid/.test(figureGallery) && /FigureCard/.test(figureGallery));
must("BoundaryPage is the only boundary controller",
  importsShared(boundary) && /PendingAnalysisCard/.test(boundary) && /BoundaryHeatmap/.test(boundary));

// App reaches the model through the dataset-bound client returned by `forDataset`
// (a proxy over the same `api` surface), so the binding — not the identifier — is
// what this asserts: the canonical method is called and the endpoint is never
// re-spelled by hand next to it.
must("App loads the canonical dataset model",
  /\bforDataset\b/.test(app) && /\.datasetModel\(\)/.test(app)
  && !/dataset-model/.test(app));
must("App mounts exactly one controller per explore page",
  (app.match(/<Overview\b/g) || []).length === 1
  && (app.match(/<GeneExplorer\b/g) || []).length === 1
  && (app.match(/<FigureGallery\b/g) || []).length === 1
  && (app.match(/<BoundaryPage\b/g) || []).length === 1);
must("App has no dataset-specific controller branch",
  !/sharedExploratory|GenericGeneExplorer|CoreExplorer/.test(app));
must("Legacy GenericGeneExplorer controller is deleted",
  !existsSync(new URL("pages/GenericGeneExplorer.jsx", ROOT)));
must("Legacy BoundaryConsistencyPage controller is deleted",
  !existsSync(new URL("pages/BoundaryConsistencyPage.jsx", ROOT)));
must("Canonical endpoint exists in API client",
  /datasetModel:\s*\(\)\s*=>\s*dget\("\/api\/runs\/current\/dataset-model"\)/.test(api));
must("Explore controllers do not call legacy/shared/core APIs",
  !/api\.(summary|species|figures|shared|core|boundaryConsistency|evidenceStack)/.test(
    overview + geneExplorer + figureGallery + boundary));
must("Gene Explorer has one Workspace and no generic workspace fork",
  (geneExplorer.match(/function Workspace\b/g) || []).length === 1
  && !/SharedGeneWorkspace|GenericGeneWorkspace/.test(geneExplorer));
// The generic layer no longer has an Evidence tab: cassette evidence is a property
// of the validated FGFR2 layer, and the exploratory layer carries its own Candidate
// Evidence tab instead. What must stay true is that neither layer gets a private
// fork of a renderer both use, and that the retired generic route still resolves.
must("Both event layers enter the same Summary/Isoforms renderers",
  /<SummaryTab[^>]+eventType=\{eventType\}/.test(geneExplorer)
  && /<IsoformsTab[^>]+eventType=\{eventType\}/.test(geneExplorer)
  && (geneExplorer.match(/function (SummaryTab|IsoformsTab|EvidenceTab)\b/g) || []).length === 3);
must("Evidence renders only for the validated event layer and keeps no dead route",
  /<EvidenceTab[^>]*\/>/.test(geneExplorer)
  && /effTab === "evidence" && isValidated/.test(geneExplorer)
  && /tab === "evidence" \?/.test(geneExplorer));
must("Architecture and synteny use canonical viewer entry points",
  /<CoordinateTrack/.test(geneExplorer) && /<SyntenyViewer/.test(geneExplorer)
  && !/<ProteinArchitectureTrack|<SyntenyTrack/.test(geneExplorer));

let failed = 0;
for (const c of checks) {
  console.log(`${c.ok ? "PASS" : "FAIL"}  ${c.name}`);
  if (!c.ok) failed += 1;
}
if (failed) { console.error(`\n${failed} parity check(s) failed.`); process.exit(1); }
console.log(`\nAll ${checks.length} shared-component parity checks passed.`);
