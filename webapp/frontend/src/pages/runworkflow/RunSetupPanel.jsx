import { useMemo, useRef, useState } from "react";

const PRESETS = [
  { id: "full30", label: "Full 30 species", hint: "Validated FGFR2 reference panel (30 vertebrates)" },
  { id: "pilot", label: "Pilot panel", hint: "Small panel for a quick test run" },
  { id: "custom", label: "Custom species list", hint: "Paste or upload scientific names" },
];

const SPECIES_ID_RE = /^[a-z][a-z0-9]+_[a-z0-9_]+$/;

function normalizeToken(raw) {
  return String(raw || "")
    .trim()
    .toLowerCase()
    .replace(/\s+/g, "_")
    .replace(/_{2,}/g, "_")
    .replace(/^_+|_+$/g, "");
}

// Mirrors scripts/create_new_run.py:species_error_message.
function speciesErrorMessage(raw) {
  const suggestion = normalizeToken(raw);
  if (suggestion && SPECIES_ID_RE.test(suggestion)) {
    return `Invalid species identifier: "${raw}". Did you mean "${suggestion}"?`;
  }
  return `Invalid species identifier: "${raw}". Expected a lowercase genus_species identifier, e.g. "homo_sapiens".`;
}

// Client-side species normalization (mirrors scripts/create_new_run.py). The
// backend re-normalizes and validates on create, so this is only for preview.
function normalizeSpeciesText(text) {
  const seen = new Set();
  const species = [];
  const invalid = [];
  let duplicates = 0;
  for (const rawLine of String(text || "").replace(/[;,]/g, "\n").split("\n")) {
    const s = rawLine.trim();
    if (!s || s.startsWith("#")) continue;
    const token = normalizeToken(s);
    if (!token) continue;
    if (!SPECIES_ID_RE.test(token)) {
      invalid.push({ raw: s, message: speciesErrorMessage(s) });
      continue;
    }
    if (seen.has(token)) { duplicates += 1; continue; }
    seen.add(token);
    species.push(token);
  }
  return { species, count: species.length, duplicates, invalid };
}

export default function RunSetupPanel({ onCreate, creating }) {
  const [preset, setPreset] = useState("full30");
  const [runName, setRunName] = useState("");
  const [speciesText, setSpeciesText] = useState("");
  const [error, setError] = useState("");
  const fileRef = useRef(null);

  const isCustom = preset === "custom";

  const preview = useMemo(
    () => (isCustom && speciesText.trim() ? normalizeSpeciesText(speciesText) : null),
    [speciesText, isCustom],
  );

  function onFile(e) {
    const file = e.target.files?.[0];
    if (!file) return;
    const reader = new FileReader();
    reader.onload = () => setSpeciesText(String(reader.result || ""));
    reader.readAsText(file);
  }

  async function submit() {
    setError("");
    if (isCustom && preview?.invalid?.length) {
      setError(preview.invalid[0].message);
      return;
    }
    if (isCustom && (!preview || preview.count < 1)) {
      setError("Add at least one species (one scientific name per line).");
      return;
    }
    const payload = {
      run_name: runName.trim() || (isCustom ? "custom_run" : `${preset}_run`),
      case_study: "FGFR2_IIIb_IIIc",
    };
    if (isCustom) payload.species_text = speciesText;
    else payload.preset = preset;
    try {
      await onCreate?.(payload);
      if (isCustom) setSpeciesText("");
      setRunName("");
    } catch (e) {
      setError(e?.message || "Could not create the run.");
    }
  }

  const customCount = preview?.count ?? 0;

  return (
    <div className="run-setup card">
      <h3>New run</h3>

      <div className="rs-field">
        <label>Case study</label>
        <div className="rs-fixed">FGFR2 IIIb/IIIc <span className="muted small">(fixed)</span></div>
      </div>

      <div className="rs-field">
        <label>Species panel</label>
        <div className="rs-presets">
          {PRESETS.map((p) => (
            <label key={p.id} className={`rs-preset${preset === p.id ? " active" : ""}`}>
              <input type="radio" name="preset" value={p.id}
                     checked={preset === p.id} onChange={() => setPreset(p.id)} />
              <span className="rs-preset-label">{p.label}</span>
              <span className="rs-preset-hint">{p.hint}</span>
            </label>
          ))}
        </div>
      </div>

      <div className="rs-field">
        <label>Run name</label>
        <input className="rs-input" type="text" value={runName} placeholder="e.g. fgfr2 test"
               onChange={(e) => setRunName(e.target.value)} />
      </div>

      {isCustom && (
        <div className="rs-field">
          <label>Species list</label>
          <textarea className="rs-textarea" rows={7} value={speciesText}
                    placeholder={"homo_sapiens\nmus_musculus\ngallus_gallus\n… (one per line or comma-separated)"}
                    onChange={(e) => setSpeciesText(e.target.value)} />
          <div className="rs-textarea-tools">
            <button className="btn ghost small" type="button" onClick={() => fileRef.current?.click()}>
              Upload .txt
            </button>
            <input ref={fileRef} type="file" accept=".txt,text/plain" hidden onChange={onFile} />
            <span className="muted small">
              {customCount > 0
                ? `${customCount} species${preview?.duplicates ? ` · ${preview.duplicates} duplicate(s) merged` : ""}`
                : "names are normalized to gallus_gallus style; blank lines and # comments are ignored"}
            </span>
          </div>

          {customCount > 0 && (
            <div className="rs-preview">
              <div className="rs-preview-label muted small">Normalized identifiers</div>
              <div className="rs-chips">
                {preview.species.slice(0, 40).map((sp) => (
                  <span key={sp} className="rs-chip">{sp}</span>
                ))}
                {preview.species.length > 40 && (
                  <span className="rs-chip more">+{preview.species.length - 40} more</span>
                )}
              </div>
            </div>
          )}

          {preview?.invalid?.length > 0 && (
            <div className="rs-invalid">
              {preview.invalid.map((iv) => (
                <div key={iv.raw} className="rs-invalid-row">{iv.message}</div>
              ))}
            </div>
          )}
        </div>
      )}

      {error && <div className="rs-error">{error}</div>}

      <button className="btn primary rs-create" onClick={submit} disabled={creating}>
        {creating ? "Creating…" : "Create run"}
      </button>
      <p className="rs-note muted small">
        Creates a local <code>runs/&lt;run_id&gt;/</code> folder only. No pipeline, InterProScan,
        pyTMHMM, SSH or SLURM is executed.
      </p>
    </div>
  );
}
