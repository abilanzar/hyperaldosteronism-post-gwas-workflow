# Synthetic post-GWAS evidence workflow

This repository is a public-safe, executable demonstration of how I organize post-GWAS analyses after association discovery. It uses entirely synthetic data and does **not** contain participant-level data, unpublished summary statistics, manuscript results, collaborator files, or controlled-access resources.

The demo connects five questions:

1. Which association signals belong to distinct loci?
2. Which variants enter a 95% credible set, and did the fine-mapping output pass basic QC?
3. What does a harmonized inverse-variance weighted Mendelian-randomization analysis estimate, and what sensitivity diagnostics are needed?
4. Does a polygenic score add information beyond covariates, and which phenotypes associate with that score after multiple-testing correction?
5. How can fine-mapping, transcriptomic, colocalization, and gene-prioritization evidence be combined without calling a gene causal?

## Run the complete demo

Requirements: Python 3.10+ and base R 4.2+.

```bash
bash scripts/run_demo.sh
```

The workflow deterministically generates synthetic inputs, runs each analysis module, writes tabular results, and executes output-level tests.

## Repository map

```text
config/module_manifest.tsv       Question, inputs, assumptions, and claim ceiling
data/example/                    Deterministically generated synthetic inputs
scripts/generate_synthetic_data.py
scripts/locus_and_finemap_qc.py  Locus grouping and credible-set audit
scripts/audit_ldsc.py            LDSC result checks and interpretation limits
scripts/mr_ivw.R                 Allele harmonization, IVW MR, heterogeneity, Egger, LOO
scripts/prs_phewas.R             PRS association/AUC and covariate-adjusted PheWAS
scripts/integrate_evidence.py     Evidence-tier integration across loci and genes
results/                         Reproducible synthetic outputs
tests/test_outputs.py            Release-gate checks
docs/analysis_plan.md            Purpose, inputs, assumptions, failure modes, and go/no-go rules
docs/interview_refresher.md       Concise interview explanations
docs/publication_boundary.md     What was intentionally excluded
```

## Scientific boundaries

- Distance-based grouping in this demo illustrates locus bookkeeping; it is not a replacement for FUMA, conditional analysis, or ancestry-matched LD clumping.
- Credible sets prioritize variants under a specified model and LD structure; they do not prove causality.
- Colocalization supports compatibility with a shared signal under its assumptions; it does not establish mediation.
- Mendelian randomization depends on instrument relevance, independence, and exclusion-restriction assumptions.
- PRS performance is population-, phenotype-, and dataset-specific and is not evidence of clinical utility by itself.
- PheWAS findings are association screens requiring multiplicity control and independent replication.
- LDSC estimates depend on summary-statistic QC, ancestry-matched LD scores, phenotype scale, and prevalence assumptions.

See [docs/analysis_plan.md](docs/analysis_plan.md) for the full interpretation framework.

## Provenance

The code structure reflects experience with a multi-biobank hyperaldosteronism genetics project. All values generated here are fictional and must not be cited as findings from that study.
