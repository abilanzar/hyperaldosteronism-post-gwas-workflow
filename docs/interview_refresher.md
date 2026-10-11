# Post-GWAS interview refresher

## 90-second explanation

After the GWAS meta-analysis, I treated post-GWAS work as a sequence of distinct questions rather than one omnibus analysis. First, I used locus annotation to define reproducible genomic risk loci and separate locus, lead-SNP, and candidate-SNP concepts. I then used fine-mapping to prioritize variants while keeping the LD source, model configuration, and credible-set uncertainty explicit. LDSC addressed common-variant heritability and shared genetic architecture. TWAS, colocalization, functional annotation, and PoPS contributed complementary gene-prioritization evidence, but none alone established a causal gene. Bidirectional MR tested prespecified directional hypotheses with allele harmonization and pleiotropy/sensitivity checks. Finally, PRS and PheWAS evaluated target-cohort association, incremental prediction, broader phenotypic associations, and external EHR-linked validation. The main discipline was provenance: every claim had to trace to the exact input, software/settings, output, population, and interpretation ceiling.

## What each module contributes

| Module | Question answered | What it cannot establish alone |
|---|---|---|
| FUMA/locus annotation | Where are the reportable association loci and what annotations overlap them? | Causal variant or causal gene |
| Fine-mapping | Which variants are most compatible with the association signal under the model? | Proof of causality |
| LDSC | How heritable is the trait on the analyzed scale, and which traits share polygenic architecture? | Direction or causality |
| TWAS | Is genetically predicted expression associated with trait liability in a reference tissue? | Direct effect of changing expression |
| Colocalization | Are two association patterns compatible with a shared causal signal? | Mediation or treatment direction |
| PoPS/gene ranking | Which nearby genes are prioritized by genome-wide features? | Statistical significance or causal-gene status |
| MR | Is genetically proxied exposure liability associated with an outcome under IV assumptions? | Definitive causation when assumptions are uncertain |
| PRS | Does aggregate genetic liability associate with or improve prediction of the target phenotype? | Clinical utility or cross-ancestry portability |
| PheWAS | Which phenotypes associate with the score across the EHR phenome? | Downstream causal effects |
| External validation | Are effect directions or score associations reproducible in another cohort? | Universal generalizability |

## High-yield distinctions

- A genomic risk locus is a reporting region; a lead SNP is a representative associated variant; an independent significant SNP is defined using an LD rule; these are not interchangeable.
- A 95% credible set means the set reaches 0.95 posterior coverage under the fitted model—not that every variant has a 95% probability of being causal.
- PP4 is evidence for a shared signal under colocalization priors and model assumptions. PP3 versus PP4, priors, tissue relevance, LD, and multiple signals must be considered.
- MR instrument harmonization is not clerical: an allele flip changes the sign of the outcome effect, and unresolved palindromic variants can invalidate direction.
- PRS association and PRS prediction are different questions. Incremental AUC should be accompanied by calibration, uncertainty, and independent validation.
- PheWAS requires a defined testing family, minimum case counts, convergence checks, and FDR/Bonferroni control.

## Interview questions to rehearse

### Why did you not call the top gene causal?

Gene prioritization sources are correlated and assumption-dependent. A nearest gene, TWAS hit, high PP4, or top PoPS rank can prioritize follow-up, but causal-gene claims require stronger triangulation and functional validation.

### What were the most important QC points?

Genome build and allele alignment; ancestry-matched LD; exact locus definitions; tool versions and parameters; instrument strength and pleiotropy for MR; discovery-target independence for PRS; minimum case counts and multiplicity for PheWAS; and population/phenotype comparability across validation cohorts.

### What would make you reject an analysis?

Unresolved build/allele errors, ancestry-mismatched LD, failed fine-mapping convergence, weak or pleiotropic MR instruments, PRS leakage or sample overlap, sparse PheWAS outcomes, or a result that cannot be traced to a locked input and exact run.

### How did you avoid overclaiming?

I assigned every module an interpretation ceiling and reconciled the claim across source output, derived table, figure, manuscript text, and resume language. Association, prediction, shared architecture, and causality were kept separate.

## What this repository demonstrates

The synthetic demo shows executable examples of locus bookkeeping, credible-set QC, LDSC output auditing, MR harmonization/IVW diagnostics, PRS model comparison, PheWAS multiple-testing correction, and conservative evidence integration. It intentionally does not publish or reproduce the underlying study results.
