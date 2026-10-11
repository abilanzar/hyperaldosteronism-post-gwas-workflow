#!/usr/bin/env python3
"""Generate deterministic fictional inputs for the public post-GWAS demo."""

from __future__ import annotations

import csv
import math
import random
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "example"


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def make_summary_stats() -> None:
    rows = []
    specs = [
        (1, 10_000_000, -0.18, 2.0e-10),
        (5, 55_000_000, 0.22, 4.0e-12),
        (12, 90_000_000, -0.15, 8.0e-9),
    ]
    alleles = [("A", "G"), ("C", "T"), ("G", "A"), ("T", "C")]
    for locus_index, (chrom, center, beta0, p0) in enumerate(specs, start=1):
        for j, offset in enumerate((-180_000, -60_000, -10_000, 0, 18_000, 85_000, 210_000)):
            ea, oa = alleles[(locus_index + j) % len(alleles)]
            p = p0 * (1 + j * 2.7)
            rows.append(
                {
                    "variant_id": f"rsSYN{locus_index}{j + 1:02d}",
                    "chr": chrom,
                    "pos": center + offset,
                    "build": "GRCh37",
                    "effect_allele": ea,
                    "other_allele": oa,
                    "beta": f"{beta0 * (1 - j * 0.035):.5f}",
                    "se": f"{0.025 + j * 0.0015:.5f}",
                    "p": f"{p:.8g}",
                    "eaf": f"{0.12 + 0.07 * ((j + locus_index) % 5):.3f}",
                    "n": 120000,
                }
            )
    for j in range(9):
        rows.append(
            {
                "variant_id": f"rsNULL{j + 1:02d}",
                "chr": 2 + j,
                "pos": 20_000_000 + j * 1_500_000,
                "build": "GRCh37",
                "effect_allele": "A",
                "other_allele": "C",
                "beta": f"{(-0.03 + j * 0.007):.5f}",
                "se": "0.03000",
                "p": f"{0.02 + j * 0.07:.5f}",
                "eaf": "0.250",
                "n": 120000,
            }
        )
    write_tsv(
        DATA / "summary_stats.tsv",
        ["variant_id", "chr", "pos", "build", "effect_allele", "other_allele", "beta", "se", "p", "eaf", "n"],
        rows,
    )


def make_finemapping() -> None:
    rows = []
    pips = {
        "L1": [("rsSYN104", 0.62), ("rsSYN105", 0.22), ("rsSYN103", 0.12), ("rsSYN106", 0.04)],
        "L2": [("rsSYN204", 0.48), ("rsSYN205", 0.28), ("rsSYN203", 0.19), ("rsSYN206", 0.05)],
        "L3": [("rsSYN304", 0.79), ("rsSYN303", 0.11), ("rsSYN305", 0.06), ("rsSYN302", 0.04)],
    }
    for locus, values in pips.items():
        for variant, pip in values:
            rows.append({"locus_id": locus, "variant_id": variant, "pip": f"{pip:.4f}"})
    write_tsv(DATA / "fine_mapping.tsv", ["locus_id", "variant_id", "pip"], rows)


def make_ldsc_metrics() -> None:
    rows = [
        {
            "analysis": "h2_observed",
            "trait_1": "synthetic_ha",
            "trait_2": "",
            "estimate": "0.085",
            "se": "0.020",
            "intercept": "1.006",
            "intercept_se": "0.008",
            "sample_prevalence": "0.500",
            "population_prevalence": "0.003",
            "ld_reference": "synthetic_EUR_reference",
        },
        {
            "analysis": "genetic_correlation",
            "trait_1": "synthetic_ha",
            "trait_2": "synthetic_blood_pressure",
            "estimate": "0.310",
            "se": "0.090",
            "intercept": "",
            "intercept_se": "",
            "sample_prevalence": "",
            "population_prevalence": "",
            "ld_reference": "synthetic_EUR_reference",
        },
    ]
    write_tsv(
        DATA / "ldsc_metrics.tsv",
        [
            "analysis",
            "trait_1",
            "trait_2",
            "estimate",
            "se",
            "intercept",
            "intercept_se",
            "sample_prevalence",
            "population_prevalence",
            "ld_reference",
        ],
        rows,
    )


def make_mr() -> None:
    rows = []
    variants = [
        ("rsMR01", "A", "G", 0.12, 0.022, 0.035, 0.014, False),
        ("rsMR02", "C", "T", 0.10, 0.020, 0.029, 0.013, True),
        ("rsMR03", "G", "C", 0.14, 0.026, 0.046, 0.016, False),
        ("rsMR04", "T", "C", 0.09, 0.019, 0.021, 0.012, False),
        ("rsMR05", "A", "C", 0.11, 0.021, 0.031, 0.014, True),
        ("rsMR06", "G", "T", 0.13, 0.023, 0.038, 0.015, False),
        ("rsMR07", "A", "T", 0.10, 0.024, 0.030, 0.016, False),
        ("rsMR08", "C", "A", 0.08, 0.019, 0.018, 0.013, False),
    ]
    for snp, ea, oa, bx, sx, by, sy, flip in variants:
        out_ea, out_oa = (oa, ea) if flip else (ea, oa)
        out_beta = -by if flip else by
        rows.append(
            {
                "variant_id": snp,
                "effect_allele_exposure": ea,
                "other_allele_exposure": oa,
                "beta_exposure": bx,
                "se_exposure": sx,
                "effect_allele_outcome": out_ea,
                "other_allele_outcome": out_oa,
                "beta_outcome": out_beta,
                "se_outcome": sy,
            }
        )
    write_tsv(
        DATA / "mr_instruments.tsv",
        [
            "variant_id",
            "effect_allele_exposure",
            "other_allele_exposure",
            "beta_exposure",
            "se_exposure",
            "effect_allele_outcome",
            "other_allele_outcome",
            "beta_outcome",
            "se_outcome",
        ],
        rows,
    )


def logistic(x: float) -> float:
    return 1.0 / (1.0 + math.exp(-x))


def make_prs_and_phewas() -> None:
    rng = random.Random(20261011)
    cohort = []
    phenotypes = []
    effects = {
        "phe_401": 0.65,
        "phe_250": 0.42,
        "phe_585": -0.15,
        "phe_327": 0.00,
        "phe_710": 0.20,
        "phe_999": 0.00,
    }
    for i in range(1, 501):
        prs = rng.gauss(0, 1)
        age = rng.uniform(40, 75)
        sex = rng.randint(0, 1)
        pc1 = rng.gauss(0, 0.03)
        pc2 = rng.gauss(0, 0.03)
        p_case = logistic(-2.65 + 0.72 * prs + 0.022 * (age - 55) + 0.18 * sex)
        case = int(rng.random() < p_case)
        cohort.append(
            {
                "id": f"SYN{i:04d}",
                "case": case,
                "prs": f"{prs:.6f}",
                "age": f"{age:.3f}",
                "sex": sex,
                "PC1": f"{pc1:.6f}",
                "PC2": f"{pc2:.6f}",
            }
        )
        row = {"id": f"SYN{i:04d}"}
        for j, (phecode, effect) in enumerate(effects.items()):
            baseline = -1.65 - 0.10 * (j % 3)
            probability = logistic(baseline + effect * prs + 0.012 * (age - 55) + 0.10 * sex)
            row[phecode] = int(rng.random() < probability)
        phenotypes.append(row)
    write_tsv(DATA / "target_cohort.tsv", ["id", "case", "prs", "age", "sex", "PC1", "PC2"], cohort)
    write_tsv(DATA / "phecode_matrix.tsv", ["id", *effects.keys()], phenotypes)


def make_gene_evidence() -> None:
    rows = [
        {"locus_id": "L1", "gene": "GENE_A", "nearest_gene": 1, "twas_z": -5.4, "twas_p": "6.7e-8", "coloc_pp4": 0.94, "pops_rank": 1, "annotation_score": 3},
        {"locus_id": "L1", "gene": "GENE_B", "nearest_gene": 0, "twas_z": -1.1, "twas_p": "0.27", "coloc_pp4": 0.18, "pops_rank": 2, "annotation_score": 1},
        {"locus_id": "L2", "gene": "GENE_C", "nearest_gene": 0, "twas_z": 4.8, "twas_p": "1.6e-6", "coloc_pp4": 0.86, "pops_rank": 1, "annotation_score": 2},
        {"locus_id": "L2", "gene": "GENE_D", "nearest_gene": 1, "twas_z": 0.6, "twas_p": "0.55", "coloc_pp4": 0.22, "pops_rank": 2, "annotation_score": 2},
        {"locus_id": "L3", "gene": "GENE_E", "nearest_gene": 1, "twas_z": 2.9, "twas_p": "0.0037", "coloc_pp4": 0.72, "pops_rank": 2, "annotation_score": 3},
        {"locus_id": "L3", "gene": "GENE_F", "nearest_gene": 0, "twas_z": 3.2, "twas_p": "0.0014", "coloc_pp4": 0.81, "pops_rank": 1, "annotation_score": 1},
    ]
    write_tsv(
        DATA / "gene_evidence.tsv",
        ["locus_id", "gene", "nearest_gene", "twas_z", "twas_p", "coloc_pp4", "pops_rank", "annotation_score"],
        rows,
    )


def main() -> None:
    make_summary_stats()
    make_finemapping()
    make_ldsc_metrics()
    make_mr()
    make_prs_and_phewas()
    make_gene_evidence()
    print(f"Synthetic inputs written to {DATA}")


if __name__ == "__main__":
    main()
