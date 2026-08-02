// One canonical species order for every comparative view.
//
// A reader compares species by scanning down a column, so row order carries
// meaning: neighbouring rows should be related species. Views used to sort
// independently — alphabetically here, by whatever the index returned there —
// and the same two species swapped places between two figures of one dataset.
//
// The order is TAXONOMIC, not phylogenetic. It groups species by clade and
// keeps the approved arrangement of the validated 30-species reference panel.
// No tree is involved, so nothing here may be labelled phylogenetic.
//
// This mirrors scripts/shared_gene_analysis/species_order.py, and a test
// compares the two so they cannot drift. When an index carries its own
// `species_order` document, that document wins: it was built from the dataset.

// The approved reference-panel order, with each species' clade.
const REFERENCE_PANEL = [
  ["homo_sapiens", "mammal"],
  ["pan_troglodytes", "mammal"],
  ["gorilla_gorilla_gorilla", "mammal"],
  ["pongo_abelii", "mammal"],
  ["macaca_mulatta", "mammal"],
  ["callithrix_jacchus", "mammal"],
  ["mus_musculus", "mammal"],
  ["rattus_norvegicus", "mammal"],
  ["oryctolagus_cuniculus", "mammal"],
  ["canis_lupus_familiaris", "mammal"],
  ["felis_catus", "mammal"],
  ["bos_taurus", "mammal"],
  ["sus_scrofa", "mammal"],
  ["equus_caballus", "mammal"],
  ["ovis_aries", "mammal"],
  ["monodelphis_domestica", "mammal"],
  ["ornithorhynchus_anatinus", "mammal"],
  ["gallus_gallus", "bird"],
  ["meleagris_gallopavo", "bird"],
  ["taeniopygia_guttata", "bird"],
  ["anolis_carolinensis", "reptile"],
  ["alligator_mississippiensis", "reptile"],
  ["chrysemys_picta_bellii", "reptile"],
  ["xenopus_tropicalis", "amphibian"],
  ["ambystoma_mexicanum", "amphibian"],
  ["danio_rerio", "fish"],
  ["oryzias_latipes", "fish"],
  ["gasterosteus_aculeatus", "fish"],
  ["takifugu_rubripes", "fish"],
  ["oreochromis_niloticus", "fish"],
];

export const CLADE_ORDER = ["mammal", "bird", "reptile", "amphibian", "fish",
  "invertebrate", "other"];

export const CLADE_LABELS = {
  mammal: "Mammals", bird: "Birds", reptile: "Reptiles",
  amphibian: "Amphibians", fish: "Fishes", invertebrate: "Invertebrates",
  other: "Other",
};

const PANEL_RANK = new Map(REFERENCE_PANEL.map(([id], i) => [id, i]));
const PANEL_CLADE = new Map(REFERENCE_PANEL);

export function normaliseSpeciesId(value) {
  return String(value || "").trim().toLowerCase().replace(/\s+/g, "_");
}

export function cladeOf(speciesId) {
  return PANEL_CLADE.get(normaliseSpeciesId(speciesId)) || "other";
}

// `gallus_gallus` -> `Gallus gallus`; only the genus is capitalised.
export function speciesDisplayName(speciesId) {
  const parts = normaliseSpeciesId(speciesId).split("_").filter(Boolean);
  if (!parts.length) return String(speciesId || "");
  return [parts[0].charAt(0).toUpperCase() + parts[0].slice(1), ...parts.slice(1)].join(" ");
}

// A lookup from species_id to position, built from an index's own species_order
// document when it has one, and from the shared rule otherwise.
export function speciesRank(orderDoc) {
  const rows = orderDoc?.species;
  if (Array.isArray(rows) && rows.length) {
    const rank = new Map();
    rows.forEach((r, i) => {
      const id = normaliseSpeciesId(r.species_id ?? r.species ?? r);
      if (id) rank.set(id, Number.isFinite(r.display_order) ? r.display_order : i);
    });
    return rank;
  }
  return null;
}

function fallbackKey(speciesId) {
  const id = normaliseSpeciesId(speciesId);
  if (PANEL_RANK.has(id)) return [0, PANEL_RANK.get(id), ""];
  const clade = cladeOf(id);
  const rank = CLADE_ORDER.indexOf(clade);
  return [1, rank < 0 ? CLADE_ORDER.length : rank, speciesDisplayName(id)];
}

function compare(a, b, rank) {
  const ida = normaliseSpeciesId(a);
  const idb = normaliseSpeciesId(b);
  if (rank) {
    const ra = rank.has(ida) ? rank.get(ida) : Number.MAX_SAFE_INTEGER;
    const rb = rank.has(idb) ? rank.get(idb) : Number.MAX_SAFE_INTEGER;
    if (ra !== rb) return ra - rb;
    // Species the document does not cover still need a defined order.
  }
  const ka = fallbackKey(ida);
  const kb = fallbackKey(idb);
  if (ka[0] !== kb[0]) return ka[0] - kb[0];
  if (ka[1] !== kb[1]) return ka[1] - kb[1];
  return String(ka[2]).localeCompare(String(kb[2]));
}

// A comparator for the canonical order, for callers that already sort rows.
export function speciesCompare(a, b, orderDoc = null) {
  return compare(a, b, speciesRank(orderDoc));
}

// Sort species identifiers into the canonical order.
export function orderSpeciesIds(ids, orderDoc = null) {
  const rank = speciesRank(orderDoc);
  return [...new Set((ids || []).map(normaliseSpeciesId).filter(Boolean))]
    .sort((a, b) => compare(a, b, rank));
}

// Sort arbitrary rows that carry a species identifier. `getId` defaults to the
// field names the indices use.
export function orderSpeciesRows(rows, getId = defaultId, orderDoc = null) {
  const rank = speciesRank(orderDoc);
  return [...(rows || [])].sort((a, b) => compare(getId(a), getId(b), rank));
}

function defaultId(row) {
  if (typeof row === "string") return row;
  return row?.species_id ?? row?.species ?? row?.speciesId ?? "";
}

// Group ordered rows by clade, keeping the canonical order inside each group.
export function groupByClade(rows, getId = defaultId, orderDoc = null) {
  const ordered = orderSpeciesRows(rows, getId, orderDoc);
  const groups = new Map();
  for (const row of ordered) {
    const clade = cladeOf(getId(row));
    if (!groups.has(clade)) groups.set(clade, []);
    groups.get(clade).push(row);
  }
  return CLADE_ORDER
    .filter((c) => groups.has(c))
    .map((c) => ({ clade: c, label: CLADE_LABELS[c], rows: groups.get(c) }));
}

export const ORDERING_METHOD = "taxonomic";
export const ORDERING_NOTE =
  "Species are ordered taxonomically (clade groups; the validated reference " +
  "panel keeps its approved order). This is not a phylogenetic tree order.";
