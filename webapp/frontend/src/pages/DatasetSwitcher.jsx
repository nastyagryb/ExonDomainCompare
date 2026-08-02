import { Badge } from "../ui";
import { DATASET_STATUS_META, datasetStatusLabel, normalizeDatasetStatus } from "../datasetStatus";

// The switcher reads the same canonical vocabulary as every page, so its badge
// can never disagree with the Summary / Boundary status.
function statusMeta(status) {
  const canonical = normalizeDatasetStatus(status);
  const [cls] = DATASET_STATUS_META[canonical] || DATASET_STATUS_META.unavailable;
  return [cls, datasetStatusLabel(canonical)];
}

function optionLabel(d) {
  if (d.kind === "example") return "Example 30-species dataset";
  if (d.switcher_label) return d.switcher_label;
  return `${d.label} — ${statusMeta(d.status)[1]}`;
}

export default function DatasetSwitcher({ datasets, activeId, onChange }) {
  const active = datasets.find((d) => d.id === activeId);
  // Robustness: if the active dataset is not (yet) in the list, still show a
  // synthetic option for it so the dropdown label can never silently fall back
  // to "Example" while the pages display a different dataset (the state bug).
  const options = active
    ? datasets
    : [{ id: activeId, kind: "run", label: activeId, switcher_label: `${activeId} (loading…)` }, ...datasets];
  const activeMeta = active || options[0];
  return (
    <div className="dataset-switcher" title="Choose which dataset all pages explore">
      <span className="ds-label">Dataset</span>
      <div className="ds-select-wrap">
        <select
          className="ds-select"
          value={activeId}
          onChange={(e) => onChange(e.target.value)}
        >
          {options.map((d) => (
            <option key={d.id} value={d.id}>
              {optionLabel(d)}
            </option>
          ))}
        </select>
      </div>
      {activeMeta && (
        <Badge cls={statusMeta(activeMeta.status)[0]} soft>
          {activeMeta.kind === "example"
            ? "Validated · read-only"
            : statusMeta(activeMeta.status)[1]}
        </Badge>
      )}
    </div>
  );
}
