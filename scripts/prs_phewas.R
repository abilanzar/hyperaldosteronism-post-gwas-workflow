#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) stop("Usage: prs_phewas.R <target_cohort.tsv> <phecode_matrix.tsv> <output_directory>")

cohort <- read.delim(args[[1]], stringsAsFactors = FALSE)
phenotypes <- read.delim(args[[2]], stringsAsFactors = FALSE)
outdir <- args[[3]]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

required <- c("id", "case", "prs", "age", "sex", "PC1", "PC2")
if (!all(required %in% names(cohort))) stop("Target cohort is missing required columns")
if (anyDuplicated(cohort$id) || anyDuplicated(phenotypes$id)) stop("Participant IDs must be unique")

cohort$prs_z <- as.numeric(scale(cohort$prs))
fit_null <- glm(case ~ age + sex + PC1 + PC2, family = binomial(), data = cohort)
fit_full <- glm(case ~ prs_z + age + sex + PC1 + PC2, family = binomial(), data = cohort)

auc_rank <- function(y, score) {
  n1 <- sum(y == 1)
  n0 <- sum(y == 0)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  ranks <- rank(score, ties.method = "average")
  (sum(ranks[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

coef_prs <- summary(fit_full)$coefficients["prs_z", ]
metrics <- data.frame(
  metric = c("n", "cases", "prs_log_or_per_sd", "prs_se", "prs_p", "prs_or_per_sd", "auc_covariates", "auc_covariates_plus_prs", "delta_auc"),
  value = c(
    nrow(cohort), sum(cohort$case), coef_prs["Estimate"], coef_prs["Std. Error"], coef_prs["Pr(>|z|)"],
    exp(coef_prs["Estimate"]), auc_rank(cohort$case, predict(fit_null, type = "response")),
    auc_rank(cohort$case, predict(fit_full, type = "response")),
    auc_rank(cohort$case, predict(fit_full, type = "response")) - auc_rank(cohort$case, predict(fit_null, type = "response"))
  ),
  interpretation = c(
    "synthetic target sample size", "synthetic case count", "association estimate, not causality", "model uncertainty",
    "association test", "odds ratio per SD in this synthetic target", "apparent discrimination",
    "apparent discrimination with PRS", "incremental apparent AUC; external validation required"
  )
)

dat <- merge(cohort, phenotypes, by = "id", all = FALSE)
phecodes <- setdiff(names(phenotypes), "id")
phewas <- data.frame()
for (phecode in phecodes) {
  outcome <- dat[[phecode]]
  cases <- sum(outcome == 1)
  controls <- sum(outcome == 0)
  if (cases < 20 || controls < 20) next
  model_data <- dat
  model_data$outcome <- outcome
  fit <- try(glm(outcome ~ prs_z + age + sex + PC1 + PC2, family = binomial(), data = model_data), silent = TRUE)
  if (inherits(fit, "try-error")) next
  coefficient <- summary(fit)$coefficients["prs_z", ]
  phewas <- rbind(
    phewas,
    data.frame(
      phecode = phecode,
      cases = cases,
      controls = controls,
      beta_prs = coefficient["Estimate"],
      se = coefficient["Std. Error"],
      odds_ratio = exp(coefficient["Estimate"]),
      p = coefficient["Pr(>|z|)"]
    )
  )
}

phewas$fdr <- p.adjust(phewas$p, method = "BH")
phewas$bonferroni <- p.adjust(phewas$p, method = "bonferroni")
phewas$interpretation <- "Covariate-adjusted association screen; replication required"
phewas <- phewas[order(phewas$p), ]

write.table(metrics, file.path(outdir, "prs_metrics.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(phewas, file.path(outdir, "phewas_results.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
cat("PRS/PheWAS complete:", nrow(cohort), "participants and", nrow(phewas), "phenotypes\n")
