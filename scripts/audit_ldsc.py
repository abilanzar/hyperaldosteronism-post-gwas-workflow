#!/usr/bin/env python3
"""Apply reporting/QC checks to compact synthetic LDSC metrics."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    with args.input.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    output = []
    for row in rows:
        estimate = float(row["estimate"])
        se = float(row["se"])
        z = estimate / se
        p = math.erfc(abs(z) / math.sqrt(2.0))
        status = "PASS"
        notes = []
        if not row["ld_reference"]:
            status = "REVIEW"
            notes.append("missing LD reference")
        if row["analysis"] == "h2_observed":
            if not row["sample_prevalence"] or not row["population_prevalence"]:
                status = "REVIEW"
                notes.append("missing prevalence inputs")
            intercept = float(row["intercept"])
            intercept_se = float(row["intercept_se"])
            if intercept - 1.96 * intercept_se > 1.05:
                status = "REVIEW"
                notes.append("intercept suggests residual inflation")
            ceiling = "Common-variant heritability conditional on QC and scale assumptions"
        else:
            if not -1 <= estimate <= 1:
                status = "FAIL"
                notes.append("genetic correlation outside [-1, 1]")
            ceiling = "Shared polygenic architecture; not causality"
        output.append(
            {
                "analysis": row["analysis"],
                "trait_1": row["trait_1"],
                "trait_2": row["trait_2"],
                "estimate": f"{estimate:.4f}",
                "se": f"{se:.4f}",
                "z": f"{z:.4f}",
                "p": f"{p:.6g}",
                "status": status,
                "interpretation_ceiling": ceiling,
                "notes": "; ".join(notes) if notes else "synthetic demo checks passed",
            }
        )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as handle:
        fields = ["analysis", "trait_1", "trait_2", "estimate", "se", "z", "p", "status", "interpretation_ceiling", "notes"]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(output)
    print(f"Audited {len(output)} synthetic LDSC result rows")


if __name__ == "__main__":
    main()
