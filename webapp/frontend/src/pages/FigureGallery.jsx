import { useEffect, useMemo, useState } from "react";
import { API_BASE, fileUrl, runFileUrl } from "../api";
import { Badge, Drawer, Modal } from "../ui";
import { FigureCardGrid, FigureCard } from "../components/shared";
import { useScientificSelection } from "../components/ScientificSelectionContext";
import CassetteExplorer from "./viewers/CassetteExplorer";
import CoordinateTrack from "./viewers/CoordinateTrack";
import MsaExplorer from "./viewers/MsaExplorer";
import SyntenyViewer from "./viewers/SyntenyViewer";
import DomainArchitecture from "./viewers/DomainArchitecture";
import { orderSpeciesIds } from "./viewers/speciesOrder";

// figure number -> interactive viewer
function interactiveKind(f) {
  const n = (f.number || "").toUpperCase();
  if (f.group === "Boundary consistency") return "boundary";
  if (f.group === "Domain & exon-boundary") return "architecture";
  if (["3B", "6B", "6"].includes(n)) return "cassette";
  if (["2", "3C"].includes(n)) return "coordinates";
  if (["5", "7"].includes(n)) return "msa";
  if (f.group === "Synteny") return "synteny";
  return null;
}
const VIEWER_TITLES = {
  cassette: "Cassette Sequence Explorer",
  coordinates: "Exon → Protein Coordinate Track",
  msa: "MSA Explorer",
  synteny: "Synteny Locus Viewer",
  architecture: "Domain Architecture (post-InterPro · pyTMHMM)",
};

function readInitialScope(multiSpecies) {
  if (!multiSpecies) return "comparative";
  try {
    const q = new URLSearchParams(window.location.search).get("figureScope");
    if (q) return q;
    const saved = sessionStorage.getItem("edc.figureScope");
    if (saved) return saved;
  } catch { /* ignore */ }
  return "comparative";
}

export default function FigureGallery({ model, openBoundary }) {
  const index = normalizeFigureIndex(model);
  const runId = model?.dataset?.kind === "run"
    ? (model?.run_id || model?.dataset?.run_id || "") : "";
  const vi = model?.validated_event_indices || {};
  const gene = model?.dataset_info?.gene_symbol || model?.gene_symbol || "";
  const selection = useScientificSelection();
  const setSelectedFigureScope = selection?.setSelectedFigureScope;
  // Scope selector for multi-species datasets: Comparative by default, then one
  // entry per real species. Single-species datasets hide the selector.
  const speciesList = useMemo(() => {
    const seen = new Map();
    for (const f of index.figures) {
      if (f.species_id && !seen.has(f.species_id)) {
        seen.set(f.species_id, f.species || f.species_id);
      }
    }
    // Prefer scientific names from the coordinate / comparative inventory when
    // the figure cards only carried a short label.
    for (const m of model?.coordinate_models || model?.models || []) {
      if (m.species_id) seen.set(m.species_id, m.scientific_name || m.species_id);
    }
    for (const s of model?.comparative_dataset?.species_inventory || []) {
      if (s.species_id) seen.set(s.species_id, s.scientific_name || s.species_id);
    }
    // The canonical taxonomic order, so the scope selector lists species in the
    // same sequence as every figure it opens.
    const order = model?.comparative_dataset?.species_order || null;
    return orderSpeciesIds([...seen.keys()], order)
      .map((id) => ({ id, name: seen.get(id) || id }));
  }, [index.figures, model]);
  const multiSpecies = speciesList.length > 1
    || Boolean(model?.comparative_dataset?.available)
    || (model?.coordinate_models || model?.models || []).length > 1;
  const [scope, setScope] = useState(() => readInitialScope(multiSpecies));
  const visibleScope = scope === "comparative" || speciesList.some((s) => s.id === scope)
    ? scope : "comparative";
  const [group, setGroup] = useState("all");
  const [search, setSearch] = useState("");
  const [showAll, setShowAll] = useState(false);
  const [lightbox, setLightbox] = useState(null);
  const [viewer, setViewer] = useState(null);

  // Persist Scope in the URL and session state, and mirror it into the shared
  // ScientificSelectionContext so linked pages can open the matching species view.
  useEffect(() => {
    if (!multiSpecies) return;
    try {
      sessionStorage.setItem("edc.figureScope", visibleScope);
      const url = new URL(window.location.href);
      url.searchParams.set("figureScope", visibleScope);
      window.history.replaceState({}, "", url.toString());
    } catch { /* ignore */ }
    setSelectedFigureScope?.(visibleScope);
  }, [visibleScope, multiSpecies, setSelectedFigureScope]);

  const figures = useMemo(() => {
    if (!index) return [];
    return index.figures.filter((f) => {
      // Supplements are hidden by default. The validated FGFR2 index marks them
      // per figure with `kind`, the generic indices with the category, so both
      // signals must be honoured or the curated FGFR2 set balloons to every file.
      if (!showAll && (f.kind !== "main" || f.category === SUPPLEMENTS)) return false;
      if (group !== "all" && f.category !== group) return false;
      if (multiSpecies) {
        const isComparative = f.scope === "comparative"
          || COMPARATIVE_CATEGORY_ORDER.includes(f.category)
          || (!f.species_id && String(f.species || "").toLowerCase() === "comparative");
        if (visibleScope === "comparative") {
          if (!isComparative) return false;
        } else {
          // Species-specific Scope shows that species' standalone Gallery structure
          // only — never comparative cards mixed in.
          if (isComparative) return false;
          if (f.species_id && f.species_id !== visibleScope) return false;
        }
      }
      if (search && !`${f.title} ${f.caption} ${f.scientific_question}`
        .toLowerCase().includes(search.toLowerCase())) return false;
      return true;
    });
  }, [index, group, search, showAll, multiSpecies, visibleScope]);

  // Category order depends on Scope: comparative vs species-specific reading order.
  const orderedCategories = useMemo(() => {
    const present = new Set(figures.map((f) => f.category));
    const preferred = (multiSpecies && visibleScope === "comparative")
      ? COMPARATIVE_CATEGORY_ORDER
      : CATEGORY_ORDER;
    const declared = preferred.filter((c) => present.has(c));
    const fromIndex = (index.categories || []).filter(
      (c) => present.has(c) && !declared.includes(c));
    const extra = [...present].filter((c) => !declared.includes(c) && !fromIndex.includes(c));
    return [...declared, ...fromIndex, ...extra];
  }, [figures, index.categories, multiSpecies, visibleScope]);

  const byGroup = useMemo(() => {
    const m = {};
    for (const f of figures) (m[f.category] = m[f.category] || []).push(f);
    return m;
  }, [figures]);

  const openInteractive = (kind) => {
    if (kind === "boundary") { openBoundary?.(); return; }
    setViewer({ kind });
  };

  return (
    <section className="page">
      <div className="page-head">
        <div>
          <p className="eyebrow">Curated figure center</p>
          <h2>Figure Gallery</h2>
        </div>
        <div className="filters">
          <input className="search" placeholder="Search figures…" value={search} onChange={(e) => setSearch(e.target.value)} />
          {multiSpecies && (
            <select value={visibleScope} onChange={(e) => { setScope(e.target.value); setGroup("all"); }}
              title="Comparative figures or one species"
              aria-label="Figure gallery scope">
              <option value="comparative">Comparative</option>
              {speciesList.map((s) => <option key={s.id} value={s.id}>{s.name}</option>)}
            </select>
          )}
          <select value={group} onChange={(e) => setGroup(e.target.value)}>
            <option value="all">All categories</option>
            {orderedCategories.map((g) => <option key={g} value={g}>{g}</option>)}
          </select>
          <label className="check inline">
            <input type="checkbox" checked={showAll} onChange={(e) => setShowAll(e.target.checked)} />
            <span>Show supplements</span>
          </label>
        </div>
      </div>

      {orderedCategories.map((g) => (
        byGroup[g]?.length ? (
          <div key={g} className="fig-group">
            <h3 className="group-title">{g} <span className="muted small">· {byGroup[g].length}</span></h3>
            <FigureCardGrid>
              {byGroup[g].map((f) => <FigCard key={f.id} f={f} gene={gene} runId={runId}
                onOpen={() => setLightbox(f)} onInteractive={openInteractive} />)}
            </FigureCardGrid>
          </div>
        ) : null
      ))}
      {figures.length === 0 && <p className="muted pad">No figures match the current filters.</p>}

      <Drawer
        open={Boolean(lightbox)}
        onClose={() => setLightbox(null)}
        title={lightbox?.title}
        subtitle={lightbox ? `${lightbox.group} · ${lightbox.kind}` : ""}
        footer={lightbox && (
          <div className="fig-downloads">
            {interactiveKind(lightbox) && (
              <button className="btn primary sm" onClick={() => { const k = interactiveKind(lightbox); setLightbox(null); openInteractive(k); }}>Open interactive</button>
            )}
            {["png", "svg", "pdf"].map((fmt) => lightbox.formats[fmt] && (
              <a key={fmt} className="btn ghost sm" href={assetUrl(lightbox.formats[fmt], false, runId)}>{fmt.toUpperCase()}</a>
            ))}
            {lightbox.source_table && <a className="btn ghost sm" href={assetUrl(lightbox.source_table, false, runId)}>Source table</a>}
          </div>
        )}
      >
        {lightbox && (
          <>
            {lightbox.thumbnail && (
              <a href={assetUrl(lightbox.formats.svg || lightbox.thumbnail, true, runId)} target="_blank" rel="noreferrer">
                <img className="lightbox-img" src={assetUrl(lightbox.thumbnail, true, runId)} alt={lightbox.title} />
              </a>
            )}
            {lightbox.scientific_question && (
              <p className="fig-detail-q"><b>Question:</b> {lightbox.scientific_question}</p>)}
            <CaptionEditor f={lightbox} gene={gene} />
            {/* Reproducibility metadata belongs here, not on the card face. */}
            {(lightbox.feature_sources?.length > 0 || lightbox.domain_source
              || lightbox.threshold != null) && (
              <p className="muted small fig-detail-src">
                {lightbox.domain_source && <>Domain source: {lightbox.domain_source}. </>}
                {lightbox.threshold != null && <>Near-edge threshold: {lightbox.threshold} aa. </>}
                {lightbox.feature_sources?.length > 0
                  && <>Feature sources: {lightbox.feature_sources.join(", ")}.</>}
              </p>
            )}
          </>
        )}
      </Drawer>

      <Modal open={Boolean(viewer)} onClose={() => setViewer(null)} title={viewer ? VIEWER_TITLES[viewer.kind] : ""}>
        {viewer?.kind === "cassette" && <CassetteExplorer preloaded={vi.cassette || vi.cassette_residue_index || {}} />}
        {viewer?.kind === "coordinates" && <CoordinateTrack preloaded={vi.coordinates || vi.coordinate_track_index || {}} />}
        {viewer?.kind === "msa" && <MsaExplorer preloaded={vi.msa || vi.msa_index || {}} />}
        {viewer?.kind === "synteny" && <SyntenyViewer preloaded={vi.synteny || model?.synteny || {}} />}
        {viewer?.kind === "architecture" && <DomainArchitecture preloaded={model?.domain_architecture?.species_index
          || model?.legacy_fgfr2_indices?.species_domain_architecture || {}} />}
      </Modal>
    </section>
  );
}

// Figure Gallery information architecture: species-specific and comparative
// reading orders (Part 3).
const SUPPLEMENTS = "Supplements";
const CATEGORY_ORDER = [
  "Exon structure",
  "Isoform analysis",
  // A gene whose analysis is about alternative isoforms reads its cassette
  // evidence right after the exon structure it sits in.
  "Isoform and cassette analysis",
  "Domain architecture",
  "Exon–domain boundaries",
  "Genomic context",
  "Exploratory candidates",
  // An established event is not an exploratory candidate and does not belong in
  // the same section as one.
  "FGFR2 event evidence",
  SUPPLEMENTS,
];
const COMPARATIVE_CATEGORY_ORDER = [
  "Comparative exon structure",
  "FGFR2 cassette evidence",
  "Comparative sequence analysis",
  "Comparative domain architecture",
  // Two boundary sections, deliberately not merged. The first carries the
  // validated cassette-boundary conclusions of the FGFR2 analysis; the second asks
  // where every supported coding-exon boundary of the whole protein falls and makes
  // no validated claim. One heading over both would lend the second the confidence
  // of the first.
  "FGFR2 IIIb/IIIc Boundary Consistency",
  "Comparative exon–domain boundaries",
  "Comparative genomic context",
  "Comparative isoform diversity",
  SUPPLEMENTS,
];

// Map a legacy section/group label onto one of the canonical categories, so runs
// whose indices predate the category field still land in the right place.
const LEGACY_CATEGORY = {
  "Domain architecture": "Domain architecture",
  "Exon–domain boundaries": "Exon–domain boundaries",
  "Exon structure": "Exon structure",
  "Isoform analysis": "Isoform analysis",
  "Genomic context": "Genomic context",
  "Exploratory candidates": "Exploratory candidates",
  "Comparative exon structure": "Comparative exon structure",
  "Comparative sequence analysis": "Comparative sequence analysis",
  "Comparative domain architecture": "Comparative domain architecture",
  "Comparative exon–domain boundaries": "Comparative exon–domain boundaries",
  "Comparative genomic context": "Comparative genomic context",
  "Comparative isoform diversity": "Comparative isoform diversity",
};
function categoryOf(f, fallbackTitle) {
  if (f.category && (CATEGORY_ORDER.includes(f.category)
    || COMPARATIVE_CATEGORY_ORDER.includes(f.category))) return f.category;
  const legacy = LEGACY_CATEGORY[f.section] || LEGACY_CATEGORY[f.group];
  if (legacy) return legacy;
  const id = String(f.figure_id || f.id || "").toLowerCase();
  const t = String(fallbackTitle || f.title || "").toLowerCase();
  if (f.scope === "comparative" || id.startsWith("cmp_")) {
    if (id.includes("boundary") || t.includes("boundar")) {
      return "Comparative exon–domain boundaries";
    }
    if (id.includes("domain")) return "Comparative domain architecture";
    if (id.includes("msa") || id.includes("identity") || id.includes("sequence")) {
      return "Comparative sequence analysis";
    }
    if (id.includes("synteny") || id.includes("genomic")) {
      return "Comparative genomic context";
    }
    if (id.includes("isoform")) return "Comparative isoform diversity";
    if (id.includes("exon")) return "Comparative exon structure";
  }
  if (id.includes("boundary") || t.includes("boundar")) return "Exon–domain boundaries";
  if (id.includes("domain") || t.includes("domain")) return "Domain architecture";
  if (id.includes("alignment") || t.includes("alignment")) return "Isoform analysis";
  if (id.includes("neighbourhood") || id.includes("synteny") || t.includes("neighbourhood")) return "Genomic context";
  if (id.includes("candidate") || t.includes("candidate")) return "Exploratory candidates";
  if (id.includes("exon") || t.includes("exon")) return "Exon structure";
  return SUPPLEMENTS;
}

function normalizeFigureIndex(model) {
  const legacy = model?.figures?.index || model?.figures?.legacy;
  if (legacy?.figures) return withCategories(legacy);
  if (model?.figures?.groups && model?.figures?.figures) return withCategories(model.figures);
  const canonical = model?.figures?.figures || [];
  const available = canonical.length
    ? canonical.filter((f) => f.status === "available")
    : model?.figures?.available || [];
  const pending = canonical.length
    ? canonical.filter((f) => f.status !== "available")
    : model?.figures?.pending || [];
  const group = "Analysis figures";
  const pendingGroup = "Pending analyses";
  const toFigure = (f, isPending = false) => ({
    id: f.figure_id || f.id,
    // The species-independent name of the scientific figure. A species Scope of a
    // multi-species run carries a species-suffixed id but the same figure_type as
    // the corresponding standalone single-species run.
    figure_type: f.figure_type || f.figure_id || f.id || "",
    number: "",
    title: f.title,
    scientific_question: f.scientific_question || "",
    species: f.species || f.scientific_name || "",
    species_id: f.species_id || "",
    scope: f.scope || (f.species_id ? "species" : ""),
    protein_id: f.protein_id || "",
    proteins: f.proteins || "",
    analysis_status: f.analysis_status || f.status || "",
    data_availability: f.data_availability || f.status || "",
    feature_sources: f.feature_sources || [],
    domain_source: f.domain_source || "",
    threshold: f.near_edge_threshold_aa ?? null,
    stage: f.stage || "",
    // The interpretation is the short card text; the caption is the full,
    // editable figure legend shown in the detail view and in the export dialog.
    interpretation: f.interpretation || "",
    caption: f.caption || f.interpretation
      || (isPending ? "Available after the required analysis stage completes." : ""),
    // A pending figure stays in its own scientific category with a pending badge,
    // rather than being demoted into the hidden supplements.
    category: categoryOf(f),
    group: f.section || f.group || (isPending ? pendingGroup : group),
    // A card may declare itself a supplement while still belonging to a scientific
    // category (the member-database signature figure is one). That marking has to
    // survive, or the supplement is shown as a main figure.
    kind: f.kind === "supplement" ? "supplement" : "main",
    thumbnail: f.png_url || f.png_path || f.thumbnail || "",
    // The card's selectable views and the models behind them. A rendered file is a
    // view of a card, so several files collapse onto one card here rather than
    // becoming several cards.
    modes: (f.modes || []).map((m) => ({
      mode_id: m.mode_id, label: m.label, description: m.description || "",
      is_default: Boolean(m.is_default), formats: m.formats || {},
      thumbnail: m.thumbnail || m.formats?.png || "",
    })),
    model_selection: f.model_selection || null,
    unavailable_models: f.availability?.unavailable_models
      || f.model_selection?.unavailable || [],
    formats: {
      png: f.png_url || f.png_path || f.formats?.png || "",
      svg: f.svg_url || f.svg_path || f.formats?.svg || "",
      pdf: f.pdf_url || f.pdf_path || f.formats?.pdf || "",
    },
    source_table: f.source_table || (f.source_files || [])[0] || "",
    pending: isPending,
    status: f.status || (isPending ? "pending_cluster" : "available"),
    error: f.error || "",
  });
  const all = [...available.map((f) => toFigure(f)), ...pending.map((f) => toFigure(f, true))];
  // One card per figure. SVG/PDF/PNG are formats of a single card, so a repeated
  // figure_id must never produce a second card.
  const figures = [];
  const seen = new Set();
  for (const f of all) {
    const key = f.id || `${f.category}:${f.title}`;
    if (seen.has(key)) continue;
    seen.add(key);
    figures.push(f);
  }
  const categories = [];
  for (const c of CATEGORY_ORDER) if (figures.some((f) => f.category === c)) categories.push(c);
  return { categories, groups: categories, figures };
}

// Attach categories (and de-duplicate) to an index that already carries its own
// group structure, e.g. the validated FGFR2 figure index. Its curated group names
// and their order are preserved; only a supplement group is renamed so it sorts
// last like everywhere else.
function withCategories(index) {
  const seen = new Set();
  const figures = [];
  for (const f of index.figures || []) {
    const key = f.id || f.figure_id || `${f.group}:${f.title}`;
    if (seen.has(key)) continue;
    seen.add(key);
    const group = f.group || SUPPLEMENTS;
    figures.push({ ...f, category: f.category || (group === "Supplement" ? SUPPLEMENTS : group) });
  }
  const categories = [];
  for (const f of figures) if (!categories.includes(f.category)) categories.push(f.category);
  const ordered = [...categories.filter((c) => c !== SUPPLEMENTS),
    ...categories.filter((c) => c === SUPPLEMENTS)];
  return { ...index, figures, categories: ordered, groups: ordered };
}

function assetUrl(path, inline = false, runId = "") {
  if (!path) return "";
  if (/^https?:\/\//.test(path)) return path;
  if (path.startsWith("/api/")) return `${API_BASE}${path}`;
  if (runId) return runFileUrl(runId, path, inline);
  return fileUrl(path, inline);
}

// Auto-generated caption, editable before download (Part 18). Captions appear in
// the figure detail, the export list and the downloaded caption file — never as a
// large permanent page element.
function autoCaption(f, gene) {
  if (f.caption) return f.caption;
  const bits = [];
  bits.push(`${f.title}.`);
  const ids = [gene, f.species, f.protein_id].filter(Boolean).join(", ");
  if (ids) bits.push(`${ids}.`);
  if (f.interpretation) bits.push(f.interpretation);
  if (f.domain_source) bits.push(`Domain annotation: ${f.domain_source}.`);
  if (f.threshold != null) bits.push(`Near-edge threshold ${f.threshold} aa.`);
  bits.push(f.stage === "post_cluster"
    ? "Post-cluster analysis stage."
    : f.stage === "pre_cluster" ? "Pre-cluster analysis stage." : "");
  bits.push("Exploratory analysis unless explicitly marked as validated.");
  return bits.filter(Boolean).join(" ");
}

function CaptionEditor({ f, gene }) {
  const [text, setText] = useState(() => autoCaption(f, gene));
  const download = () => {
    const blob = new Blob([`${f.title}\n\n${text}\n`], { type: "text/plain" });
    const url = URL.createObjectURL(blob);
    const a = document.createElement("a");
    a.href = url;
    a.download = `${(f.id || "figure").replace(/[^\w.-]+/g, "_")}.caption.txt`;
    document.body.appendChild(a); a.click(); a.remove();
    setTimeout(() => URL.revokeObjectURL(url), 1000);
  };
  return (
    <div className="fig-caption-editor">
      <div className="fce-head">
        <b className="small">Caption</b>
        <span className="fce-actions">
          <button className="btn ghost sm" onClick={() => setText(autoCaption(f, gene))}>Reset</button>
          <button className="btn ghost sm" onClick={download}>Download .txt</button>
        </span>
      </div>
      <textarea className="fce-text" rows={4} value={text}
        onChange={(e) => setText(e.target.value)}
        aria-label={`Editable caption for ${f.title}`} />
    </div>
  );
}

const STAGE_LABEL = { pre_cluster: "pre-cluster", post_cluster: "post-cluster" };

// A comparative card describes a dataset, not one protein, so a per-protein
// message would be the wrong scientific claim. Availability is reported at the
// level the card actually analyses.
function unavailableText(f) {
  if (f.status === "failed") return "Figure could not be produced";
  if (f.status === "not_applicable") {
    return f.scope === "comparative"
      ? "Not applicable to this dataset"
      : "Not applicable to this protein";
  }
  if (f.pending) {
    return f.scope === "comparative"
      ? "Available after the required analysis stage completes for this dataset"
      : "Available after the required analysis stage completes";
  }
  return f.scope === "comparative"
    ? "Not available for this dataset"
    : "Not available for this protein";
}

// One card per figure: the scientific question, a short cautious interpretation,
// the identity of the analysed molecule, and the download formats. Reproducibility
// provenance (source file names, tool versions) lives in the detail view and in
// the downloadable tables, not on the primary card face (Part 7).
function FigCard({ f, gene, runId, onOpen, onInteractive }) {
  const kind = interactiveKind(f);
  // Which previews failed, by url rather than a single flag: a card can switch
  // between several rendered views, and one broken view must not hide the others.
  const [failedThumbs, setFailedThumbs] = useState(() => new Set());
  // A card's views — an isoform model, a taxon filter, a coordinate system. These
  // are states of one figure, so they are a selector on the card rather than
  // separate cards; a gallery that made a card of each would bury the question the
  // figure answers under its own variants.
  const modes = f.modes || [];
  const [modeId, setModeId] = useState(() =>
    (modes.find((m) => m.is_default) || modes[0])?.mode_id || "");
  const mode = modes.find((m) => m.mode_id === modeId) || null;
  const formats = mode?.formats || f.formats;
  const thumbnail = mode?.thumbnail || f.thumbnail;
  const withheld = f.unavailable_models || [];
  const unavailable = f.pending || f.status === "failed"
    || (modes.length === 0 && withheld.length > 0);
  const idMeta = [gene, f.species, f.protein_id].filter(Boolean).join(" · ");
  const short = f.interpretation || f.caption || "";

  return (
    <FigureCard
      thumbOnClick={onOpen}
      thumb={thumbnail && !failedThumbs.has(thumbnail)
        ? <img key={thumbnail} src={assetUrl(thumbnail, true, runId)} alt={f.title} loading="lazy"
            onError={() => setFailedThumbs((s) => new Set(s).add(thumbnail))} />
        : <div className="fig-noimg">{f.error
          || (unavailable ? unavailableText(f) : "Preview unavailable")}</div>}
      title={f.title}
      badge={<Badge cls={unavailable ? "neutral" : "accepted"} soft>
        {unavailable ? (f.status === "failed" ? "failed" : "unavailable")
          : (STAGE_LABEL[f.stage] || "available")}
      </Badge>}
      caption={<>
        {f.scientific_question && <span className="fig-q">{f.scientific_question}</span>}
        {short ? `${short.slice(0, 170)}${short.length > 170 ? "…" : ""}` : null}
        {modes.length > 1 && (
          <span className="fig-modes">
            {modes.map((m) => (
              <button key={m.mode_id} type="button" title={m.description || m.label}
                className={`chip sm${m.mode_id === modeId ? " on" : ""}`}
                onClick={() => setModeId(m.mode_id)}>{m.label}</button>
            ))}
          </span>
        )}
        {/* An absent model is named with the reason the analysis recorded, so a
            missing view is never mistaken for a broken page. */}
        {withheld.map((u) => (
          <span key={`${u.isoform}${u.reason}`} className="fig-unavailable small muted">
            {u.isoform}: not available — {u.reason}
          </span>
        ))}
        {idMeta && <span className="fig-idmeta">{idMeta}</span>}
      </>}
      actions={<>
        {kind && !unavailable && <button className="btn primary sm" onClick={() => onInteractive(kind)}>Interactive</button>}
        {!unavailable && <button className="btn ghost sm" onClick={onOpen}>Open</button>}
        {formats.svg && <a className="btn ghost sm" href={assetUrl(formats.svg, false, runId)}>SVG</a>}
        {formats.pdf && <a className="btn ghost sm" href={assetUrl(formats.pdf, false, runId)}>PDF</a>}
        {formats.png && <a className="btn ghost sm" href={assetUrl(formats.png, false, runId)}>PNG</a>}
        {f.source_table && <a className="btn ghost sm" href={assetUrl(f.source_table, false, runId)}>TSV</a>}
      </>}
    />
  );
}
