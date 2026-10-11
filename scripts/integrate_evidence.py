#!/usr/bin/env python3
"""Integrate independent-looking evidence channels without declaring causal genes."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--loci", type=Path, required=True)
    parser.add_argument("--credible-sets", type=Path, required=True)
    parser.add_argument("--evidence", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    loci = {row["locus_id"]: row for row in read_tsv(args.loci)}
    cs_by_locus: dict[str, list[str]] = defaultdict(list)
    for row in read_tsv(args.credible_sets):
        if row["in_95pct_credible_set"] == "1":
            cs_by_locus[row["locus_id"]].append(row["variant_id"])

    evidence = read_tsv(args.evidence)
    bonferroni = 0.05 / len(evidence)
    output = []
    for row in evidence:
        locus_id = row["locus_id"]
        if locus_id not in loci:
            raise ValueError(f"Evidence row references unknown locus: {locus_id}")
        fine_mapped = int(bool(cs_by_locus.get(locus_id)))
        twas_support = int(float(row["twas_p"]) < bonferroni)
        coloc_support = int(float(row["coloc_pp4"]) >= 0.80)
        pops_support = int(int(row["pops_rank"]) == 1)
        annotation_support = int(float(row["annotation_score"]) >= 2)
        nearest_gene = int(row["nearest_gene"])
        # Fine-mapping is locus-level context here; do not credit every gene at
        # the locus with the same variant-level evidence. Rank genes only on
        # gene-specific or gene-linked channels.
        evidence_count = twas_support + coloc_support + pops_support + annotation_support
        if evidence_count >= 3:
            tier = "high_priority_follow_up"
        elif evidence_count >= 2:
            tier = "moderate_priority_follow_up"
        else:
            tier = "exploratory"
        output.append(
            {
                "locus_id": locus_id,
                "sentinel_variant": loci[locus_id]["sentinel_variant"],
                "gene": row["gene"],
                "fine_mapped_locus": fine_mapped,
                "twas_bonferroni_support": twas_support,
                "coloc_pp4_ge_0.80": coloc_support,
                "pops_top_rank": pops_support,
                "functional_annotation_support": annotation_support,
                "nearest_gene": nearest_gene,
                "evidence_count": evidence_count,
                "priority_tier": tier,
                "claim_ceiling": "Convergent prioritization for follow-up; not proof that the gene is causal or druggable",
            }
        )

    output.sort(key=lambda row: (-row["evidence_count"], row["locus_id"], row["gene"]))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "locus_id",
        "sentinel_variant",
        "gene",
        "fine_mapped_locus",
        "twas_bonferroni_support",
        "coloc_pp4_ge_0.80",
        "pops_top_rank",
        "functional_annotation_support",
        "nearest_gene",
        "evidence_count",
        "priority_tier",
        "claim_ceiling",
    ]
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(output)
    print(f"Integrated {len(output)} synthetic locus-gene evidence rows")


if __name__ == "__main__":
    main()
