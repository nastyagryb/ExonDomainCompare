
from __future__ import annotations

import argparse
import re
import shutil
import subprocess
from pathlib import Path

import pandas as pd
from Bio import AlignIO, SeqIO
from Bio.SeqRecord import SeqRecord

from common import ensure_dir, read_tsv, setup_logging, write_json


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Prepare and optionally run unconstrained, isoform-constrained and species-pair-constrained IQ-TREE analyses with AU tests."
    )
    p.add_argument("--alignment", required=True)
    p.add_argument("--truth-table")
    p.add_argument("--outdir", required=True)
    p.add_argument("--iqtree", default="iqtree2")
    p.add_argument("--threads", default="AUTO")
    p.add_argument("--rell-replicates", type=int, default=10000)
    p.add_argument("--bootstrap", type=int, default=1000)
    p.add_argument("--seed-unconstrained", type=int, default=107252)
    p.add_argument("--seed-isoform-constraint", type=int, default=995853)
    p.add_argument("--seed-species-pair-constraint", type=int, default=460166)
    p.add_argument("--seed-au-test", type=int, default=160107)
    p.add_argument("--run", action="store_true")
    p.add_argument("--redo", action="store_true")
    p.add_argument("--verbose", action="store_true")
    return p.parse_args()


def safe_name(species: str, isoform: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]", "_", f"{species}__{isoform}")


def run_command(command: list[str], cwd: Path) -> None:
    subprocess.run(command, cwd=cwd, check=True)


def parse_best_model(iqtree_report: Path) -> str:
    text = iqtree_report.read_text(errors="replace")
    patterns = [
        r"Best-fit model according to BIC:\s*(\S+)",
        r"Best-fit model:\s*(\S+)",
    ]
    for pattern in patterns:
        match = re.search(pattern, text)
        if match:
            return match.group(1)
    raise RuntimeError(f"Could not parse best-fit model from {iqtree_report}")


def main() -> None:
    args = parse_args()
    setup_logging(args.verbose)
    out = ensure_dir(args.outdir)
    aln = AlignIO.read(args.alignment, "fasta")
    entries = []
    for rec in aln:
        parts = rec.id.split("|")
        if len(parts) < 3:
            raise ValueError(f"Expected species|isoform|protein ID: {rec.id}")
        entries.append({"rec": rec, "species": parts[0], "isoform": parts[1], "protein_id": parts[2]})

    if args.truth_table:
        truth = read_tsv(args.truth_table, required=["species", "final_isoform_label", "protein_id"])
        if "recommended_use_post_rescue" in truth.columns:
            truth = truth[truth["recommended_use_post_rescue"] == "main_analysis"]
        keep = set(zip(truth["species"], truth["final_isoform_label"], truth["protein_id"]))
        entries = [e for e in entries if (e["species"], e["isoform"], e["protein_id"]) in keep]

    counts = pd.DataFrame(entries).groupby(["species", "isoform"]).size().unstack(fill_value=0)
    complete = sorted(counts.index[(counts.get("IIIb", 0) == 1) & (counts.get("IIIc", 0) == 1)])
    entries = [e for e in entries if e["species"] in set(complete)]
    entries.sort(key=lambda e: (e["species"], e["isoform"]))

    renamed = []
    mapping_rows = []
    for e in entries:
        name = safe_name(e["species"], e["isoform"])
        renamed.append(SeqRecord(e["rec"].seq, id=name, description=""))
        mapping_rows.append({"safe_id": name, "species": e["species"], "isoform": e["isoform"], "protein_id": e["protein_id"]})
    fasta = out / "cassette_complete_pairs.safe_ids.faa"
    SeqIO.write(renamed, fasta, "fasta")
    pd.DataFrame(mapping_rows).to_csv(out / "sequence_id_map.tsv", sep="\t", index=False)

    b = [safe_name(s, "IIIb") for s in complete]
    c = [safe_name(s, "IIIc") for s in complete]
    iso_constraint = f"(({','.join(b)}),({','.join(c)}));\n"
    pair_constraint = "(" + ",".join(f"({safe_name(s,'IIIb')},{safe_name(s,'IIIc')})" for s in complete) + ");\n"
    (out / "constraint_isoform_monophyly.nwk").write_text(iso_constraint)
    (out / "constraint_species_pairs.nwk").write_text(pair_constraint)

    commands = [
        f"{args.iqtree} -s {fasta.name} -st AA -m MFP -B {args.bootstrap} --alrt {args.bootstrap} -T {args.threads} -seed {args.seed_unconstrained} --prefix unconstrained",
        "BEST_MODEL=$(awk '/Best-fit model according to BIC:/ {print $NF}' unconstrained.iqtree)",
        f'{args.iqtree} -s {fasta.name} -st AA -m "$BEST_MODEL" -g constraint_isoform_monophyly.nwk -T {args.threads} -seed {args.seed_isoform_constraint} --prefix isoform_constraint',
        f'{args.iqtree} -s {fasta.name} -st AA -m "$BEST_MODEL" -g constraint_species_pairs.nwk -T {args.threads} -seed {args.seed_species_pair_constraint} --prefix species_pair_constraint',
        "cat unconstrained.treefile isoform_constraint.treefile species_pair_constraint.treefile > candidate_topologies.trees",
        f'{args.iqtree} -s {fasta.name} -st AA -m "$BEST_MODEL" -n 0 -z candidate_topologies.trees -zb {args.rell_replicates} -au -seed {args.seed_au_test} --prefix topology_AU_test',
    ]
    (out / "run_commands.sh").write_text("#!/usr/bin/env bash\nset -euo pipefail\n" + "\n".join(commands) + "\n")
    (out / "run_commands.sh").chmod(0o755)

    if args.run:
        executable = shutil.which(args.iqtree)
        if not executable:
            raise FileNotFoundError(f"IQ-TREE executable not found: {args.iqtree}")
        redo = ["-redo"] if args.redo else []
        run_command(
            [executable, "-s", fasta.name, "-st", "AA", "-m", "MFP", "-B", str(args.bootstrap), "--alrt", str(args.bootstrap), "-T", args.threads, "-seed", str(args.seed_unconstrained), "--prefix", "unconstrained", *redo],
            out,
        )
        model = parse_best_model(out / "unconstrained.iqtree")
        for prefix, constraint, seed in [
            ("isoform_constraint", "constraint_isoform_monophyly.nwk", args.seed_isoform_constraint),
            ("species_pair_constraint", "constraint_species_pairs.nwk", args.seed_species_pair_constraint),
        ]:
            run_command(
                [executable, "-s", fasta.name, "-st", "AA", "-m", model, "-g", constraint, "-T", args.threads, "-seed", str(seed), "--prefix", prefix, *redo],
                out,
            )
        trees = "".join((out / f).read_text().strip() + "\n" for f in ["unconstrained.treefile", "isoform_constraint.treefile", "species_pair_constraint.treefile"])
        (out / "candidate_topologies.trees").write_text(trees)
        run_command(
            [executable, "-s", fasta.name, "-st", "AA", "-m", model, "-n", "0", "-z", "candidate_topologies.trees", "-zb", str(args.rell_replicates), "-au", "-seed", str(args.seed_au_test), "--prefix", "topology_AU_test", *redo],
            out,
        )
    else:
        existing_report = out / "unconstrained.iqtree"
        model = parse_best_model(existing_report) if existing_report.exists() else None

    write_json(
        {
            "n_complete_species_pairs": len(complete),
            "n_sequences": len(entries),
            "alignment_length": aln.get_alignment_length(),
            "run_executed": args.run,
            "existing_results_detected": (out / "topology_AU_test.iqtree").exists(),
            "best_model": model,
            "random_seeds": {
                "unconstrained": args.seed_unconstrained,
                "isoform_constraint": args.seed_isoform_constraint,
                "species_pair_constraint": args.seed_species_pair_constraint,
                "topology_AU_test": args.seed_au_test,
            },
            "hypotheses": ["unconstrained", "isoform_monophyly", "species_pairs"],
        },
        out / "iqtree_topology_test_manifest.json",
    )


if __name__ == "__main__":
    main()
