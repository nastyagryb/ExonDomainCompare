
from __future__ import annotations

import argparse
import copy
import math
import urllib.request
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml
from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley

from common import ensure_dir, read_tsv, setup_logging, write_json


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Quantitative mapping of isoform-discriminating FGFR2 residues to receptor-ligand interfaces."
    )
    p.add_argument("--positions", required=True, help="structure_mapping_17_discriminating_positions.tsv")
    p.add_argument("--config", required=True, help="YAML structure configuration")
    p.add_argument("--outdir", required=True)
    p.add_argument("--pdb-dir", help="Directory containing/downloading PDB files")
    p.add_argument("--direct-distance", type=float, default=4.5)
    p.add_argument("--near-distance", type=float, default=8.0)
    p.add_argument("--permutations", type=int, default=10000)
    p.add_argument("--position-window", type=int, default=12)
    p.add_argument("--seed", type=int, default=20260804)
    p.add_argument("--no-download", action="store_true")
    p.add_argument("--verbose", action="store_true")
    return p.parse_args()


def download_pdb(pdb_id: str, destination: Path) -> None:
    url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
    destination.parent.mkdir(parents=True, exist_ok=True)
    urllib.request.urlretrieve(url, destination)


def trim_structure(structure, keep_chains: set[str]):
    result = copy.deepcopy(structure)
    for model in list(result):
        for chain in list(model):
            if chain.id not in keep_chains:
                model.detach_child(chain.id)
    return result


def find_residue(chain, number: int):
    matches = [res for res in chain if res.id[1] == number and res.id[0] == " "]
    if not matches:
        return None
    return matches[0]


def heavy_atom_coordinates(residue) -> np.ndarray:
    coords = []
    for atom in residue.get_atoms():
        element = (atom.element or "").upper()
        if element != "H":
            coords.append(atom.coord.astype(float))
    return np.asarray(coords, dtype=float)


def min_heavy_atom_distance(residue, ligand_atoms: np.ndarray) -> float:
    receptor = heavy_atom_coordinates(residue)
    if receptor.size == 0 or ligand_atoms.size == 0:
        return np.nan

    distances = np.linalg.norm(receptor[:, None, :] - ligand_atoms[None, :, :], axis=2)
    return float(distances.min())


def residue_sasa_map(structure, chain_id: str) -> dict[int, float]:
    sr = ShrakeRupley()
    sr.compute(structure, level="R")
    model = next(structure.get_models())
    chain = model[chain_id]
    result = {}
    for residue in chain:
        if residue.id[0] == " ":
            result[int(residue.id[1])] = float(getattr(residue, "sasa", np.nan))
    return result


def matched_permutation(
    universe: pd.DataFrame,
    targets: pd.DataFrame,
    statistic_col: str,
    n_perm: int,
    position_window: int,
    seed: int,
) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    universe = universe[np.isfinite(universe[statistic_col]) & np.isfinite(universe["sasa_free"])].copy()
    targets = targets[np.isfinite(targets[statistic_col]) & np.isfinite(targets["sasa_free"])].copy()
    observed = float(targets[statistic_col].sum())
    controls = universe[~universe["is_discriminating"]].copy()
    if controls.empty or targets.empty:
        return observed, np.nan

    controls["sasa_bin"] = pd.qcut(controls["sasa_free"], q=4, labels=False, duplicates="drop")
    targets = targets.copy()
    bins = np.quantile(universe["sasa_free"].dropna(), [0.25, 0.5, 0.75])
    targets["sasa_bin"] = np.digitize(targets["sasa_free"], bins, right=True)
    null = np.zeros(n_perm, dtype=float)
    for b in range(n_perm):
        selected = []
        for _, target in targets.iterrows():
            pool = controls[
                (controls["sasa_bin"] == target["sasa_bin"])
                & ((controls["resi"] - target["resi"]).abs() <= position_window)
            ]
            if pool.empty:
                pool = controls[controls["sasa_bin"] == target["sasa_bin"]]
            if pool.empty:
                pool = controls
            selected.append(pool.iloc[rng.integers(0, len(pool))][statistic_col])
        null[b] = np.sum(selected)
    p = (1 + np.sum(null >= observed - 1e-12)) / (n_perm + 1)
    return observed, float(p)


def main() -> None:
    args = parse_args()
    setup_logging(args.verbose)
    out = ensure_dir(args.outdir)
    pdb_dir = ensure_dir(args.pdb_dir or out / "pdb")
    positions = read_tsv(args.positions)
    config = yaml.safe_load(Path(args.config).read_text())
    all_rows = []
    summary_rows = []

    for structure_cfg in config.get("structures", []):
        if not structure_cfg.get("enabled", True):
            continue
        pdb_id = str(structure_cfg["pdb_id"]).upper()
        pdb_path = pdb_dir / f"{pdb_id}.pdb"
        if not pdb_path.exists():
            if args.no_download:
                raise FileNotFoundError(f"Missing {pdb_path}; rerun without --no-download or provide the file")
            download_pdb(pdb_id, pdb_path)

        receptor_chain = str(structure_cfg["receptor_chain"])
        ligand_chains = [str(x) for x in structure_cfg["ligand_chains"]]
        mapping_column = str(structure_cfg["mapping_column"])
        cassette_start = int(structure_cfg["cassette_start_resi"])
        cassette_end = int(structure_cfg["cassette_end_resi"])

        parser = PDBParser(QUIET=True)
        original = parser.get_structure(pdb_id, str(pdb_path))
        complex_structure = trim_structure(original, {receptor_chain, *ligand_chains})
        receptor_structure = trim_structure(original, {receptor_chain})
        complex_model = next(complex_structure.get_models())
        receptor_chain_obj = complex_model[receptor_chain]

        ligand_coords = []
        for chain_id in ligand_chains:
            for residue in complex_model[chain_id]:
                if residue.id[0] == " ":
                    coords = heavy_atom_coordinates(residue)
                    if coords.size:
                        ligand_coords.append(coords)
        ligand_atoms = np.vstack(ligand_coords) if ligand_coords else np.empty((0, 3))

        sasa_complex = residue_sasa_map(complex_structure, receptor_chain)
        sasa_free = residue_sasa_map(receptor_structure, receptor_chain)
        mapped_positions = set(
            pd.to_numeric(positions[mapping_column], errors="coerce").dropna().astype(int).tolist()
        )

        rows = []
        for resi in range(cassette_start, cassette_end + 1):
            residue = find_residue(receptor_chain_obj, resi)
            if residue is None:
                continue
            distance = min_heavy_atom_distance(residue, ligand_atoms)
            free = sasa_free.get(resi, np.nan)
            complex_sasa = sasa_complex.get(resi, np.nan)
            delta = free - complex_sasa if np.isfinite(free) and np.isfinite(complex_sasa) else np.nan
            rows.append(
                {
                    "pdb_id": pdb_id,
                    "isoform": structure_cfg.get("isoform", ""),
                    "receptor_chain": receptor_chain,
                    "ligand_chains": ",".join(ligand_chains),
                    "resi": resi,
                    "resname": residue.resname,
                    "is_discriminating": resi in mapped_positions,
                    "min_ligand_distance_A": distance,
                    "direct_contact": bool(distance <= args.direct_distance) if np.isfinite(distance) else False,
                    "near_interface": bool(distance <= args.near_distance) if np.isfinite(distance) else False,
                    "sasa_free": free,
                    "sasa_complex": complex_sasa,
                    "delta_sasa_A2": delta,
                }
            )
        table = pd.DataFrame(rows)
        if table.empty:
            raise ValueError(f"No cassette residues found for {pdb_id} chain {receptor_chain}")
        table["direct_contact_numeric"] = table["direct_contact"].astype(int)
        targets = table[table["is_discriminating"]].copy()
        observed_contacts, contact_p = matched_permutation(
            table,
            targets,
            "direct_contact_numeric",
            args.permutations,
            args.position_window,
            args.seed,
        )
        observed_dsasa, dsasa_p = matched_permutation(
            table,
            targets,
            "delta_sasa_A2",
            args.permutations,
            args.position_window,
            args.seed + 1,
        )
        table.to_csv(out / f"{pdb_id}_residue_interface_metrics.tsv", sep="\t", index=False)
        all_rows.append(table)
        summary_rows.append(
            {
                "pdb_id": pdb_id,
                "isoform": structure_cfg.get("isoform", ""),
                "n_cassette_residues_observed": len(table),
                "n_discriminating_residues_mapped": len(targets),
                "n_discriminating_direct_contacts": int(targets["direct_contact"].sum()),
                "n_discriminating_near_interface": int(targets["near_interface"].sum()),
                "sum_discriminating_delta_sasa_A2": targets["delta_sasa_A2"].sum(),
                "matched_permutation_p_direct_contacts": contact_p,
                "matched_permutation_p_sum_delta_sasa": dsasa_p,
                "direct_distance_A": args.direct_distance,
                "near_distance_A": args.near_distance,
            }
        )

        fig, ax = plt.subplots(figsize=(10, 5))
        ax.scatter(table["resi"], table["min_ligand_distance_A"], label="Cassette residues")
        ax.scatter(
            targets["resi"], targets["min_ligand_distance_A"], marker="o", s=55, label="Discriminating residues"
        )
        ax.axhline(args.direct_distance, linewidth=1, label="Direct-contact threshold")
        ax.axhline(args.near_distance, linewidth=1, linestyle="--", label="Near-interface threshold")
        ax.set_xlabel("PDB residue number")
        ax.set_ylabel("Minimum heavy-atom distance to ligand (Å)")
        ax.set_title(f"{pdb_id}: quantitative ligand-interface mapping")
        ax.legend()
        fig.tight_layout()
        fig.savefig(out / f"figure_{pdb_id}_interface_distances.png", dpi=300)
        fig.savefig(out / f"figure_{pdb_id}_interface_distances.svg")
        plt.close(fig)

    combined = pd.concat(all_rows, ignore_index=True) if all_rows else pd.DataFrame()
    combined.to_csv(out / "structure_interface_metrics_all.tsv", sep="\t", index=False)
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(out / "structure_interface_enrichment_summary.tsv", sep="\t", index=False)
    write_json(
        {
            "n_structures": len(summary),
            "direct_distance_A": args.direct_distance,
            "near_distance_A": args.near_distance,
            "n_permutations": args.permutations,
            "matching": "receptor-alone SASA quartile plus local sequence-position window",
        },
        out / "structure_mapping_summary.json",
    )


if __name__ == "__main__":
    main()
