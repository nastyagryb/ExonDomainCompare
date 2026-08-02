import { useEffect, useMemo, useRef, useState } from "react";
import { api } from "../../api";

const SPECIES_ID_RE = /^[a-z][a-z0-9]+_[a-z0-9_]+$/;

function normalizeToken(raw) {
  return String(raw || "")
    .trim()
    .toLowerCase()
    .replace(/\s+/g, "_")
    .replace(/_{2,}/g, "_")
    .replace(/^_+|_+$/g, "");
}

// Fast local pre-check (backend revalidates + resolves scientific name/taxid).
// Used only to enable the button and count species before the async resolve
// call returns; the authoritative resolved preview comes from the backend.
function localSpeciesCount(text) {
  const seen = new Set();
  for (const rawLine of String(text || "").replace(/[;,]/g, "\n").split("\n")) {
    const s = rawLine.trim();
    if (!s || s.startsWith("#")) continue;
    const token = normalizeToken(s);
    if (token && SPECIES_ID_RE.test(token)) seen.add(token);
  }
  return seen.size;
}

// The form asks for the biology and nothing else: a gene, an optional name and
// a species list. Which workflow runs (validated FGFR2 vs. generic exploratory,
// single-species vs. comparative) follows from those answers, so the router is
// not a question the user is asked — it stays in the backend.
export default function CreateRunPanel({ onStart, starting }) {
  const [runName, setRunName] = useState("");
  const [gene, setGene] = useState("FGFR2");
  const [speciesText, setSpeciesText] = useState("");
  const [error, setError] = useState("");
  const [resolution, setResolution] = useState(null);  // backend resolve-inputs result
  const [resolving, setResolving] = useState(false);
  const fileRef = useRef(null);

  const geneSym = gene.trim().toUpperCase();
  const localCount = useMemo(() => localSpeciesCount(speciesText), [speciesText]);

  // Debounced call to the ONE generic backend resolver (single source of truth
  // for gene validity + workflow routing + resolved species panel).
  useEffect(() => {
    let cancelled = false;
    const handle = setTimeout(() => {
      if (cancelled) return;
      if (!geneSym && !speciesText.trim()) { setResolution(null); setResolving(false); return; }
      setResolving(true);
      api.resolveRunInputs({ gene_symbol: geneSym, species_text: speciesText, mode: "auto" })
        .then((r) => { if (!cancelled) setResolution(r); })
        .catch(() => { if (!cancelled) setResolution(null); })
        .finally(() => { if (!cancelled) setResolving(false); });
    }, 300);
    return () => { cancelled = true; clearTimeout(handle); };
  }, [geneSym, speciesText]);

  const speciesPreview = resolution?.species || null;

  function onFile(e) {
    const file = e.target.files?.[0];
    if (!file) return;
    const reader = new FileReader();
    reader.onload = () => setSpeciesText(String(reader.result || ""));
    reader.readAsText(file);
  }

  async function submit() {
    setError("");
    if (!geneSym) { setError("Enter a gene symbol (e.g. FGFR2, TP53 or FOXP1)."); return; }
    if (resolution?.gene && resolution.gene.valid === false) {
      setError(resolution.gene.message || "Enter a valid gene symbol."); return;
    }
    if (speciesPreview?.invalid?.length) { setError(speciesPreview.invalid[0].message); return; }
    const count = speciesPreview ? speciesPreview.count : localCount;
    if (count < 1) { setError("Add at least one species (one scientific name per line)."); return; }
    try {
      // An empty run name stays empty: the run is titled from its gene and
      // species instead of being stamped with a placeholder.
      await onStart?.({
        run_name: runName.trim(),
        gene_symbol: geneSym,
        mode: "auto",
        species_text: speciesText,
        // Kept for backward compatibility; the backend router decides the workflow.
        case_study: geneSym === "FGFR2" ? "FGFR2_IIIb_IIIc" : `${geneSym}_core_only_pilot`,
      });
      setSpeciesText("");
      setRunName("");
    } catch (e) {
      setError(e?.message || "Could not create and start the run.");
    }
  }

  const count = speciesPreview ? speciesPreview.count : localCount;
  const dups = speciesPreview?.duplicates_removed ?? 0;
  const geneInvalid = geneSym && resolution?.gene?.valid === false;

  return (
    <div className="run-setup card">
      <h3>Create new run</h3>
      <p className="muted small">
        Choose a gene and enter one species per line (e.g. <code>Homo sapiens</code>).
      </p>

      <div className="rs-field">
        <label>Gene symbol</label>
        <input className="rs-input" type="text" value={gene} placeholder="e.g. FGFR2, TP53, FOXP1"
               onChange={(e) => setGene(e.target.value)} />
        {geneInvalid && <div className="rs-invalid-row">{resolution.gene.message}</div>}
      </div>

      <div className="rs-field">
        <label>Run name <span className="muted small">(optional)</span></label>
        <input className="rs-input" type="text" value={runName}
               placeholder="e.g. TP53 mammals validation"
               onChange={(e) => setRunName(e.target.value)} />
      </div>

      <div className="rs-field">
        <label>Species</label>
        <textarea className="rs-textarea" rows={6} value={speciesText}
                  placeholder={"Homo sapiens\nMus musculus\n… (one per line or comma-separated)"}
                  onChange={(e) => setSpeciesText(e.target.value)} />
        <div className="rs-textarea-tools">
          <button className="btn ghost small" type="button" onClick={() => fileRef.current?.click()}>Upload .txt</button>
          <input ref={fileRef} type="file" accept=".txt,text/plain" hidden onChange={onFile} />
          <span className="muted small">
            {resolving ? "Resolving species…"
              : count > 0
                ? `${count} species${dups ? ` · ${dups} duplicate(s) merged` : ""}`
                : "blank lines and # comments are ignored"}
          </span>
        </div>

        {speciesPreview?.resolved?.length > 0 && (
          <div className="rs-preview">
            <div className="rs-preview-label muted small">Resolved species (preview)</div>
            <div className="rs-species-list">
              {speciesPreview.resolved.map((sp) => (
                <div key={sp.species_id} className={`rs-species ${sp.known ? "known" : "unknown"}`}>
                  <div className="rs-species-name">
                    {sp.scientific_name}
                    {sp.common_name ? <span className="muted small"> · {sp.common_name}</span> : null}
                  </div>
                  <div className="rs-species-meta muted small">
                    Taxonomy ID: {sp.taxid || (sp.known ? "—" : "unresolved")}
                  </div>
                  {!sp.known && (
                    <div className="rs-species-flag">
                      Not in the offline lookup — will be resolved live; the run reports a clear
                      message if it cannot be retrieved.
                    </div>
                  )}
                </div>
              ))}
            </div>
          </div>
        )}

        {speciesPreview?.invalid?.length > 0 && (
          <div className="rs-invalid">
            {speciesPreview.invalid.map((iv) => <div key={iv.raw} className="rs-invalid-row">{iv.message}</div>)}
          </div>
        )}
      </div>

      {error && <div className="rs-error">{error}</div>}

      <button className="btn primary rs-create" onClick={submit} disabled={starting || geneInvalid}>
        {starting ? "Creating & starting…" : "Create and start run"}
      </button>
    </div>
  );
}
