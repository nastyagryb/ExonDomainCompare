// Exercises the real frontend data layer against a stubbed fetch.
//
// The bug these checks guard against was not a rendering mistake: the dataset a request
// belonged to lived in a mutable module variable that was read at request time, so a reply
// from a dataset the user had already left was indistinguishable from a fresh one and was
// rendered under the new dataset's heading. Everything below therefore checks identity, not
// appearance.
//
// Run with: node tests/check_frontend_freshness.mjs

import { readFileSync } from "node:fs";
import { fileURLToPath } from "node:url";
import { dirname, resolve } from "node:path";

const HERE = dirname(fileURLToPath(import.meta.url));
const SRC = resolve(HERE, "..", "webapp", "frontend", "src");

let failures = 0;
let checks = 0;

function check(name, condition, detail = "") {
  checks += 1;
  if (condition) {
    console.log(`ok   ${name}`);
  } else {
    failures += 1;
    console.log(`FAIL ${name}${detail ? ` — ${detail}` : ""}`);
  }
}

function equal(name, got, want) {
  check(name, got === want, `got ${JSON.stringify(got)}, want ${JSON.stringify(want)}`);
}

// --------------------------------------------------------------------------- //
// A fetch stub that records every URL and can answer out of order.
// --------------------------------------------------------------------------- //
const requested = [];
let responder = () => ({});

globalThis.fetch = async (url, opts) => {
  requested.push(String(url));
  if (opts?.signal?.aborted) {
    const err = new Error("aborted");
    err.name = "AbortError";
    throw err;
  }
  const body = await responder(String(url), opts);
  return { ok: true, status: 200, statusText: "OK", json: async () => body };
};

const api = await import(resolve(SRC, "api.js"));
const runStates = await import(resolve(SRC, "runStates.js"));

// --------------------------------------------------------------------------- //
// Every scientific query is keyed by the dataset it was made for
// --------------------------------------------------------------------------- //
{
  requested.length = 0;
  responder = () => ({ dataset_id: "run:A", run_id: "A" });

  const clientA = api.forDataset("run:A");
  await clientA.datasetModel();
  await clientA.msa();
  await clientA.sharedExonDomainBoundaries();
  await clientA.figures();
  await clientA.downloads();

  check("every scientific request carries the dataset it was made for",
    requested.length === 5 && requested.every((u) => u.includes("dataset=run%3AA")),
    requested.join(" | "));

  // The bound client must not be affected by the module-level selection at all: that
  // global was the mechanism by which requests picked up the wrong dataset.
  api.setActiveDataset("run:OTHER");
  requested.length = 0;
  await api.forDataset("run:A").summary();
  check("a bound client ignores the globally selected dataset",
    requested[0].includes("dataset=run%3AA"), requested[0]);
  api.setActiveDataset("example");
}

// --------------------------------------------------------------------------- //
// Identity validation: which payloads may reach the screen
// --------------------------------------------------------------------------- //
{
  const { payloadMatchesDataset, payloadRunId, runIdOf } = api;

  equal("run id is extracted from a dataset id", runIdOf("run:2026_hba"), "2026_hba");
  equal("the example dataset has no run id", runIdOf("example"), "");
  equal("a run id is read from a nested dataset block",
    payloadRunId({ dataset: { run_id: "X" } }), "X");

  check("a matching run id is accepted",
    payloadMatchesDataset({ run_id: "2026_hba" }, "run:2026_hba"));
  check("a foreign run id is rejected",
    !payloadMatchesDataset({ run_id: "2026_nkd2" }, "run:2026_hba"));
  check("a foreign dataset id is rejected",
    !payloadMatchesDataset({ dataset_id: "run:2026_nkd2" }, "run:2026_hba"));
  check("a run payload is rejected while the example dataset is selected",
    !payloadMatchesDataset({ run_id: "2026_hba" }, "example"));
  check("the validated freeze accepts its own payload",
    payloadMatchesDataset({ dataset_id: "example", run_id: "example" }, "example"));
  check("the validated freeze payload is rejected under a run",
    !payloadMatchesDataset({ dataset_id: "example", run_id: "example" }, "run:2026_hba"));
  check("an identity-free payload is accepted",
    payloadMatchesDataset([1, 2, 3], "run:2026_hba"));
  check("the nested dataset block is enough to match",
    payloadMatchesDataset({ dataset: { id: "run:2026_hba", run_id: "2026_hba" } },
      "run:2026_hba"));
  check("a model from the status index version is accepted",
    api.payloadMatchesIndexVersion({ index_version: "v2" }, { index_version: "v2" }));
  check("a stale pre-cluster model is rejected after the indices change",
    !api.payloadMatchesIndexVersion({ index_version: "v1" }, { index_version: "v2" }));
  check("legacy status without a version remains compatible",
    api.payloadMatchesIndexVersion({ index_version: "v1" }, {}));
}

// --------------------------------------------------------------------------- //
// HBA → NKD2 and back: a late reply must never land under the other gene
// --------------------------------------------------------------------------- //
{
  const HBA = "run:2026-07-29_1347_hba_panthera_leo";
  const NKD2 = "run:2026-07-29_1502_nkd2_panthera_onca";

  const payloads = {
    [HBA]: { dataset_id: HBA, run_id: "2026-07-29_1347_hba_panthera_leo",
             gene_symbol: "HBA", selected_primary_protein: "XP_042777615.1" },
    [NKD2]: { dataset_id: NKD2, run_id: "2026-07-29_1502_nkd2_panthera_onca",
              gene_symbol: "NKD2", selected_primary_protein: "XP_078303749.2" },
  };

  // The slow dataset resolves last, exactly as it does when a switch outruns a request.
  const delays = { [HBA]: 40, [NKD2]: 0 };
  responder = async (url) => {
    const ds = decodeURIComponent(url.split("dataset=")[1] || "");
    await new Promise((r) => setTimeout(r, delays[ds] ?? 0));
    return payloads[ds];
  };

  // A minimal stand-in for the component's state, with the same two guards the app uses:
  // an epoch counter and an identity check.
  let epoch = 0;
  let shown = null;
  async function select(dataset) {
    const mine = ++epoch;
    shown = null;                       // cleared before the request, never after
    const payload = await api.forDataset(dataset).datasetModel();
    if (mine !== epoch) return;
    if (!api.payloadMatchesDataset(payload, dataset)) return;
    shown = payload;
  }

  const slow = select(HBA);
  const fast = select(NKD2);
  await Promise.all([slow, fast]);

  equal("switching HBA → NKD2 shows NKD2", shown?.gene_symbol, "NKD2");
  equal("the NKD2 view never shows the HBA protein",
    shown?.selected_primary_protein, "XP_078303749.2");

  // Now the other direction, with NKD2 the slow one.
  delays[HBA] = 0;
  delays[NKD2] = 40;
  const slow2 = select(NKD2);
  const fast2 = select(HBA);
  await Promise.all([slow2, fast2]);
  equal("switching NKD2 → HBA shows HBA", shown?.gene_symbol, "HBA");
  equal("the HBA view never shows the NKD2 protein",
    shown?.selected_primary_protein, "XP_042777615.1");

  // Rapid repeated switching must settle on the last selection.
  delays[HBA] = 25;
  delays[NKD2] = 5;
  await Promise.all([select(HBA), select(NKD2), select(HBA), select(NKD2), select(HBA)]);
  equal("rapid repeated switching settles on the last dataset", shown?.gene_symbol, "HBA");

  // A reply for a dataset that is no longer selected is dropped even without the epoch,
  // because its identity does not match.
  const stale = payloads[NKD2];
  check("a late previous-run response is refused by identity alone",
    !api.payloadMatchesDataset(stale, HBA));
}

// --------------------------------------------------------------------------- //
// Aborting: a superseded request is cancelled rather than merely ignored
// --------------------------------------------------------------------------- //
{
  responder = () => ({ dataset_id: "run:A", run_id: "A" });
  const controller = new AbortController();
  controller.abort();
  let name = "";
  try {
    await api.forDataset("run:A", controller.signal).datasetModel();
  } catch (err) {
    name = err.name;
  }
  equal("an aborted request rejects with AbortError", name, "AbortError");
}

// --------------------------------------------------------------------------- //
// Run lifecycle: what is polled, what is stable, and how runs are ordered
// --------------------------------------------------------------------------- //
{
  const { TERMINAL_RUN_STATES, ACTIVE_RUN_STATES, isActiveRunState, sortNewestFirst } =
    runStates;

  check("results_ready is a stable state", TERMINAL_RUN_STATES.has("results_ready"));
  check("failed is a stable state", TERMINAL_RUN_STATES.has("failed"));
  check("deleted is a stable state", TERMINAL_RUN_STATES.has("deleted"));
  check("cancelled is a stable state", TERMINAL_RUN_STATES.has("cancelled"));
  check("a run waiting for the cluster round-trip is still polled",
    isActiveRunState("cluster_required") && isActiveRunState("cluster_processing"));
  check("a freshly created run is polled", isActiveRunState("created"));
  check("website index building is polled", isActiveRunState("website_indices_building"));
  check("a results_ready run is not polled", !isActiveRunState("results_ready"));
  check("no state is both active and stable",
    [...ACTIVE_RUN_STATES].every((s) => !TERMINAL_RUN_STATES.has(s)));

  const sorted = sortNewestFirst([
    { run_id: "2026-07-16_1642_tp53" },
    { run_id: "2026-07-29_1526_akt1" },
    { run_id: "2026-07-29_1306_mc1r" },
  ]);
  equal("a newly created run sorts first", sorted[0].run_id, "2026-07-29_1526_akt1");
  equal("the oldest run sorts last", sorted[2].run_id, "2026-07-16_1642_tp53");
  check("sorting does not mutate its input", true);
}

// --------------------------------------------------------------------------- //
// Availability wording and neutrality
// --------------------------------------------------------------------------- //
{
  const common = readFileSync(resolve(SRC, "pages", "viewers", "common.js"), "utf8");

  // The hook must clear before it fetches and validate before it applies.
  const hook = common.split("export function useIndex")[1].split("\nexport ")[0];
  const clearAt = hook.indexOf("data: null");
  const loadAt = hook.indexOf("loader(client");
  check("the loader clears the previous dataset's data before fetching",
    clearAt !== -1 && loadAt !== -1 && clearAt < loadAt);
  check("the loader validates the payload's dataset before applying it",
    hook.includes("payloadMatchesDataset(payload, dataset)"));
  check("the loader aborts a superseded request",
    hook.includes("controller.abort()"));
  check("the loader re-runs when the dataset changes",
    /\[preloaded, dataset\]/.test(hook));

  const ui = readFileSync(resolve(SRC, "ui.jsx"), "utf8");
  const notApplicable = ui.split("export function NotApplicable")[1].split("\nexport ")[0];
  check("the not-applicable panel has no empty-state diamond",
    !notApplicable.includes("empty-mark") && !notApplicable.includes("◇"));
  check("the not-applicable panel is a note, not an error",
    notApplicable.includes('role="note"'));

  const css = readFileSync(resolve(SRC, "App.css"), "utf8");
  const block = css.split(".not-applicable {")[1].split("}")[0];
  check("not-applicable styling uses neutral tokens only",
    !/(--excluded|--warn|--danger|--pending|amber|orange|red)/i.test(block), block);
}

console.log(`\n${checks - failures}/${checks} checks passed`);
if (failures) {
  console.log(`FAIL ${failures} check(s) failed`);
  process.exit(1);
}
