// Repository-relative path resolution for the developer tooling scripts.
// Absolute paths from a single machine must never be hardcoded here, otherwise
// these scripts only work in one checkout.
import { fileURLToPath } from "node:url";
import { dirname, resolve } from "node:path";

const HERE = dirname(fileURLToPath(import.meta.url));

export const FRONTEND_ROOT = resolve(HERE, "..");
export const REPO_ROOT = resolve(HERE, "../../..");

export const artifactPath = (...parts) => resolve(REPO_ROOT, "artifacts", ...parts);
export const frontendPath = (...parts) => resolve(FRONTEND_ROOT, ...parts);
