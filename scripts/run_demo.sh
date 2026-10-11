#!/usr/bin/env bash
set -euo pipefail

repo_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_dir"

mkdir -p data/example results

python3 scripts/generate_synthetic_data.py

python3 scripts/locus_and_finemap_qc.py \
  --summary data/example/summary_stats.tsv \
  --finemapping data/example/fine_mapping.tsv \
  --outdir results

python3 scripts/audit_ldsc.py \
  --input data/example/ldsc_metrics.tsv \
  --output results/ldsc_audit.tsv

Rscript scripts/mr_ivw.R \
  data/example/mr_instruments.tsv \
  results

Rscript scripts/prs_phewas.R \
  data/example/target_cohort.tsv \
  data/example/phecode_matrix.tsv \
  results

python3 scripts/integrate_evidence.py \
  --loci results/locus_catalog.tsv \
  --credible-sets results/credible_sets.tsv \
  --evidence data/example/gene_evidence.tsv \
  --output results/integrated_gene_evidence.tsv

python3 tests/test_outputs.py

echo "Synthetic post-GWAS demo completed successfully."
