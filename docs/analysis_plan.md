# Recommended post-GWAS analysis plan

This plan starts only after the discovery summary statistics pass build, allele, identifier, frequency, sample-size, imputation-quality, and duplication checks.

## 1. Locus definition and annotation

- **Purpose:** Convert variant-level associations into reproducible reporting units and identify candidate regulatory or coding annotations.
- **Inputs needed:** Harmonized summary statistics, genome build, ancestry-matched LD, significance threshold, and locus-definition rules.
- **Method:** Use FUMA or a documented clumping/conditional strategy. Record genomic risk loci, independent significant variants, lead SNPs, and candidate SNPs as distinct objects.
- **Key assumptions:** LD reference matches the analyzed ancestry; coordinates and alleles are consistent.
- **QC checks:** Confirm variant counts, locus boundaries, build, LD panel, excluded variants, and a summary-statistic-to-annotation crosswalk.
- **Primary outputs:** Locus catalog, lead-variant table, functional annotations, and provenance manifest.
- **Interpretation limit:** Defines reporting loci and candidate annotations; the nearest or mapped gene is not automatically causal.
- **Go/no-go:** Stop if build, allele orientation, or LD reference cannot be established.

## 2. Fine-mapping

- **Purpose:** Prioritize variants compatible with one or more association signals at a locus.
- **Inputs needed:** Dense locus summary statistics, matched LD, sample-size information, and model settings.
- **Method:** Apply an appropriate single- or multi-ancestry fine-mapping method; retain PIPs, 95% credible sets, and diagnostics.
- **Key assumptions:** LD is accurate and the causal-signal model is reasonable.
- **QC checks:** PIP range/sums, allele and variant alignment, credible-set coverage, convergence, seed/version capture, and sensitivity to LD/model choices.
- **Primary outputs:** Variant-level PIPs, credible sets, locus-level status, and no-credible-set table.
- **Interpretation limit:** Credible sets prioritize variants; they do not prove a variant is causal.
- **Go/no-go:** Do not interpret a locus with sparse coverage, severe LD mismatch, failed convergence, or unclear configuration.

## 3. LDSC heritability and genetic correlation

- **Purpose:** Estimate common-variant heritability and shared polygenic architecture with relevant traits.
- **Inputs needed:** Well-calibrated ancestry-specific summary statistics, sample size or effective sample size, prevalence assumptions for binary traits, and ancestry-matched LD scores.
- **Method:** Munge summary statistics, run univariate LDSC, then prespecified cross-trait LDSC with multiplicity correction.
- **Key assumptions:** Polygenic model, appropriate LD scores, valid phenotype scale and prevalence, and controlled sample overlap.
- **QC checks:** Regression SNP count, mean chi-square, lambda GC, intercept, attenuation ratio, h2 Z score, liability-scale sensitivity, and source-release manifest.
- **Primary outputs:** Observed/liability h2, intercept, ratio, genetic correlations, standard errors, and adjusted P values.
- **Interpretation limit:** Genetic correlation is shared architecture, not directionality or causation.
- **Go/no-go:** Do not emphasize estimates with poor calibration, mismatched ancestry, unstable prevalence conversion, or very low h2 Z score.

## 4. TWAS, colocalization, and PoPS/gene prioritization

- **Purpose:** Prioritize genes and tissues whose genetically predicted expression or feature-based evidence aligns with association loci.
- **Inputs needed:** GWAS summary statistics, expression-weight/eQTL resources, matched LD, locus windows, colocalization priors, and gene-feature matrices.
- **Method:** Run tissue-aware TWAS, evaluate colocalization at relevant loci, and compare with feature-based gene prioritization and functional annotations.
- **Key assumptions:** Reference tissues are relevant, alleles/LD are aligned, and the colocalization model is suitable for the number of signals.
- **QC checks:** Weight coverage, multiple testing, locus boundaries, PP3/PP4 and prior sensitivity, conditional signals, and agreement or conflict across tissues.
- **Primary outputs:** Tissue-gene TWAS table, colocalization posterior table, ranked genes, and an evidence matrix.
- **Interpretation limit:** Convergent evidence prioritizes follow-up; it does not establish mediation, causal genes, or treatment direction.
- **Go/no-go:** Treat tissue mismatch, ancestry mismatch, low overlap, or unresolved multiple signals as reasons to downgrade evidence.

Fine-mapping is retained as locus-level context rather than automatically credited to every gene in the locus. A gene requires gene-linked evidence to move up the prioritization tier.

## 5. Mendelian randomization

- **Purpose:** Evaluate whether genetically proxied liability is associated with prespecified outcomes and, when justified, the reverse direction.
- **Inputs needed:** Independent strong instruments, compatible exposure/outcome summary statistics, ancestry/build alignment, and sample-overlap information.
- **Method:** Harmonize alleles; run IVW as primary where appropriate; add weighted-median, MR-Egger, heterogeneity, pleiotropy, leave-one-out, and MR-PRESSO diagnostics.
- **Key assumptions:** Instrument relevance, independence, exclusion restriction, and correct directionality.
- **QC checks:** Instrument strength, palindromic/mismatched variants, harmonization, heterogeneity, Egger intercept, outliers, leave-one-out stability, sample overlap, and multiple testing.
- **Primary outputs:** Main and sensitivity estimates, instrument table, diagnostics, exclusions, and multiplicity-adjusted results.
- **Interpretation limit:** Supports causal inference only under assumptions; isolated nominal estimates are not definitive evidence.
- **Go/no-go:** Downgrade or stop for weak instruments, unresolved allele mismatches, dominant single variants, or strong unaddressed pleiotropy.

## 6. PRS evaluation

- **Purpose:** Test whether polygenic liability associates with the target phenotype and adds prediction beyond covariates.
- **Inputs needed:** Discovery weights, independent target genotypes/phenotypes, covariates, ancestry labels, and a prespecified tuning strategy.
- **Method:** Build scores using a justified approach, tune without test-set leakage, then compare covariate-only and covariate-plus-PRS models.
- **Key assumptions:** Discovery and target samples do not overlap; phenotype and ancestry are sufficiently compatible.
- **QC checks:** Variant matching, allele alignment, score distribution, relatedness, ancestry structure, calibration, AUC with uncertainty, and nested cross-validation where tuning occurs.
- **Primary outputs:** Association per SD, score distributions, calibration, incremental discrimination, and ancestry-stratified results.
- **Interpretation limit:** Dataset-specific performance; not clinical utility without prospective and external validation.
- **Go/no-go:** Do not claim prediction if tuning leaks into testing, overlap is unresolved, or target sample size is inadequate.

## 7. PheWAS and external EHR validation

- **Purpose:** Screen phenotypic consequences of polygenic liability and test reproducibility in an independent EHR-linked biobank.
- **Inputs needed:** PRS, phenotype-code mapping, covariates, minimum case thresholds, and an independent validation design.
- **Method:** Fit a common covariate-adjusted model across eligible phecodes, correct for multiple testing, and compare effect direction/magnitude across cohorts.
- **Key assumptions:** Coding density, ascertainment, model, and phenotype definitions are adequate and reasonably comparable.
- **QC checks:** Case/control counts, separation/convergence, related phenotypes, FDR/Bonferroni families, ancestry stratification, and discovery-validation mapping.
- **Primary outputs:** Full PheWAS table, significance summary, category plot, and cross-cohort replication table.
- **Interpretation limit:** Hypothesis-generating phenotypic associations; not downstream causal effects.
- **Go/no-go:** Suppress unstable or sparse phenotypes and avoid replication claims without a prespecified mapping and directionally consistent external estimate.

## Suggested order

1. Summary-statistic QC and harmonization
2. Locus definition and annotation
3. Fine-mapping
4. LDSC heritability/genetic correlation
5. TWAS, colocalization, PoPS, and functional integration
6. Prespecified bidirectional MR
7. PRS development and independent evaluation
8. PheWAS and external EHR-linked validation
9. Finding-level evidence reconciliation across every output, table, figure, and claim

## Claims allowed

- A locus, variant, gene, tissue, or phenotype can be described as prioritized or associated under the stated method.
- Multiple aligned modules can be described as convergent evidence.
- MR can be described as supporting causal inference only with explicit assumptions and sensitivity evidence.
- PRS can be described using target-specific association, calibration, and discrimination metrics.

## Claims not allowed

- The nearest gene is the causal gene.
- A credible-set variant is proven causal.
- Colocalization proves mediation.
- A TWAS Z-score directly specifies an agonist or antagonist.
- Genetic correlation proves causality.
- A nominal MR result proves a causal effect.
- A statistically associated PRS is clinically useful or portable across ancestries.
- A PheWAS association is a disease consequence of the exposure.
