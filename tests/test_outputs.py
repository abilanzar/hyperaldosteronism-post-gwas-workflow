#!/usr/bin/env python3
"""Small release-gate checks for the deterministic demo outputs."""

from __future__ import annotations

import csv
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"


def read(name: str) -> list[dict[str, str]]:
    path = RESULTS / name
    assert path.exists() and path.stat().st_size > 0, f"Missing or empty output: {name}"
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def numeric(value: str) -> float:
    result = float(value)
    assert math.isfinite(result)
    return result


def main() -> None:
    loci = read("locus_catalog.tsv")
    assert len(loci) == 3
    assert {row["build"] for row in loci} == {"GRCh37"}

    fine_qc = read("finemapping_qc.tsv")
    assert len(fine_qc) == 3
    assert all(row["status"] == "PASS" for row in fine_qc)

    ldsc = read("ldsc_audit.tsv")
    assert len(ldsc) == 2
    assert all(row["status"] == "PASS" for row in ldsc)

    mr = read("mr_summary.tsv")
    ivw = next(row for row in mr if row["method"] == "IVW")
    assert int(ivw["n_instruments"]) >= 5
    assert numeric(ivw["se"]) > 0

    prs = read("prs_metrics.tsv")
    metric = {row["metric"]: numeric(row["value"]) for row in prs}
    assert metric["n"] == 500
    assert 0 <= metric["auc_covariates"] <= 1
    assert 0 <= metric["auc_covariates_plus_prs"] <= 1

    phewas = read("phewas_results.tsv")
    assert len(phewas) == 6
    assert all(0 <= numeric(row["fdr"]) <= 1 for row in phewas)

    integrated = read("integrated_gene_evidence.tsv")
    assert len(integrated) == 6
    assert any(row["priority_tier"] == "high_priority_follow_up" for row in integrated)
    assert all("not proof" in row["claim_ceiling"] for row in integrated)

    print("All output release-gate checks passed")


if __name__ == "__main__":
    main()
