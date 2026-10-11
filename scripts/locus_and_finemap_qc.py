#!/usr/bin/env python3
"""Build a simple locus catalog and audit synthetic fine-mapping output."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def make_loci(summary: list[dict[str, str]], threshold: float, merge_distance: int) -> list[dict]:
    required = {
        "variant_id",
        "chr",
        "pos",
        "build",
        "effect_allele",
        "other_allele",
        "beta",
        "se",
        "p",
        "eaf",
        "n",
    }
    if not summary or not required.issubset(summary[0]):
        missing = sorted(required - set(summary[0] if summary else []))
        raise ValueError(f"Missing summary-statistic fields: {missing}")
    builds = {row["build"] for row in summary}
    if len(builds) != 1:
        raise ValueError(f"Mixed genome builds: {sorted(builds)}")
    significant = [row for row in summary if float(row["p"]) < threshold]
    significant.sort(key=lambda row: (int(row["chr"]), int(row["pos"])))
    groups: list[list[dict[str, str]]] = []
    for row in significant:
        if not groups:
            groups.append([row])
            continue
        prior = groups[-1]
        same_chrom = int(row["chr"]) == int(prior[-1]["chr"])
        close = int(row["pos"]) - int(prior[-1]["pos"]) <= merge_distance
        if same_chrom and close:
            prior.append(row)
        else:
            groups.append([row])
    loci = []
    for index, group in enumerate(groups, start=1):
        sentinel = min(group, key=lambda row: float(row["p"]))
        loci.append(
            {
                "locus_id": f"L{index}",
                "chr": sentinel["chr"],
                "start": min(int(row["pos"]) for row in group),
                "end": max(int(row["pos"]) for row in group),
                "sentinel_variant": sentinel["variant_id"],
                "sentinel_p": sentinel["p"],
                "n_significant_variants": len(group),
                "build": sentinel["build"],
            }
        )
    return loci


def audit_finemapping(rows: list[dict[str, str]], valid_loci: set[str]) -> tuple[list[dict], list[dict]]:
    required = {"locus_id", "variant_id", "pip"}
    if not rows or not required.issubset(rows[0]):
        raise ValueError("Fine-mapping input must contain locus_id, variant_id, and pip")
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    seen: set[tuple[str, str]] = set()
    for row in rows:
        key = (row["locus_id"], row["variant_id"])
        if key in seen:
            raise ValueError(f"Duplicate fine-mapping row: {key}")
        seen.add(key)
        pip = float(row["pip"])
        if not 0 <= pip <= 1:
            raise ValueError(f"Invalid PIP for {key}: {pip}")
        grouped[row["locus_id"]].append(row)

    output = []
    qc = []
    for locus_id in sorted(grouped):
        ordered = sorted(grouped[locus_id], key=lambda row: float(row["pip"]), reverse=True)
        total = sum(float(row["pip"]) for row in ordered)
        cumulative = 0.0
        threshold_reached = False
        for row in ordered:
            cumulative += float(row["pip"])
            include = not threshold_reached
            output.append(
                {
                    "locus_id": locus_id,
                    "variant_id": row["variant_id"],
                    "pip": f"{float(row['pip']):.4f}",
                    "cumulative_pip": f"{cumulative:.4f}",
                    "in_95pct_credible_set": int(include),
                }
            )
            if cumulative >= 0.95:
                threshold_reached = True
        qc.append(
            {
                "locus_id": locus_id,
                "locus_in_catalog": int(locus_id in valid_loci),
                "pip_sum": f"{total:.4f}",
                "credible_set_reaches_0.95": int(cumulative >= 0.95),
                "status": "PASS" if locus_id in valid_loci and abs(total - 1.0) < 1e-6 and cumulative >= 0.95 else "REVIEW",
            }
        )
    return output, qc


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--finemapping", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--p-threshold", type=float, default=5e-8)
    parser.add_argument("--merge-distance", type=int, default=1_000_000)
    args = parser.parse_args()

    loci = make_loci(read_tsv(args.summary), args.p_threshold, args.merge_distance)
    credible, qc = audit_finemapping(read_tsv(args.finemapping), {row["locus_id"] for row in loci})
    write_tsv(
        args.outdir / "locus_catalog.tsv",
        ["locus_id", "chr", "start", "end", "sentinel_variant", "sentinel_p", "n_significant_variants", "build"],
        loci,
    )
    write_tsv(
        args.outdir / "credible_sets.tsv",
        ["locus_id", "variant_id", "pip", "cumulative_pip", "in_95pct_credible_set"],
        credible,
    )
    write_tsv(
        args.outdir / "finemapping_qc.tsv",
        ["locus_id", "locus_in_catalog", "pip_sum", "credible_set_reaches_0.95", "status"],
        qc,
    )
    print(f"Defined {len(loci)} synthetic reporting loci; audited {len(qc)} fine-mapped loci")


if __name__ == "__main__":
    main()
