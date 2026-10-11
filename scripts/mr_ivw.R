#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("Usage: mr_ivw.R <mr_instruments.tsv> <output_directory>")

input_path <- args[[1]]
outdir <- args[[2]]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

dat <- read.delim(input_path, stringsAsFactors = FALSE, check.names = FALSE)
required <- c(
  "variant_id", "effect_allele_exposure", "other_allele_exposure",
  "beta_exposure", "se_exposure", "effect_allele_outcome",
  "other_allele_outcome", "beta_outcome", "se_outcome"
)
if (!all(required %in% names(dat))) stop("MR input is missing required columns")

is_palindromic <- function(a1, a2) {
  paste(sort(c(a1, a2)), collapse = "") %in% c("AT", "CG")
}

dat$harmonization_status <- "excluded_incompatible"
dat$beta_outcome_harmonized <- NA_real_
for (i in seq_len(nrow(dat))) {
  ex <- c(dat$effect_allele_exposure[i], dat$other_allele_exposure[i])
  out <- c(dat$effect_allele_outcome[i], dat$other_allele_outcome[i])
  if (is_palindromic(ex[1], ex[2])) {
    dat$harmonization_status[i] <- "excluded_palindromic"
  } else if (identical(ex, out)) {
    dat$harmonization_status[i] <- "aligned"
    dat$beta_outcome_harmonized[i] <- dat$beta_outcome[i]
  } else if (identical(ex, rev(out))) {
    dat$harmonization_status[i] <- "flipped"
    dat$beta_outcome_harmonized[i] <- -dat$beta_outcome[i]
  }
}

keep <- dat[dat$harmonization_status %in% c("aligned", "flipped"), ]
if (nrow(keep) < 3) stop("Fewer than three harmonized instruments")

w <- 1 / (keep$se_outcome ^ 2)
denom <- sum(w * keep$beta_exposure ^ 2)
ivw_beta <- sum(w * keep$beta_exposure * keep$beta_outcome_harmonized) / denom
ivw_se <- sqrt(1 / denom)
ivw_z <- ivw_beta / ivw_se
ivw_p <- 2 * pnorm(-abs(ivw_z))

residual <- keep$beta_outcome_harmonized - ivw_beta * keep$beta_exposure
q <- sum(w * residual ^ 2)
q_df <- nrow(keep) - 1
q_p <- pchisq(q, df = q_df, lower.tail = FALSE)

egger <- lm(beta_outcome_harmonized ~ beta_exposure, data = keep, weights = w)
egger_coef <- summary(egger)$coefficients

loo <- data.frame()
for (i in seq_len(nrow(keep))) {
  subset <- keep[-i, ]
  sw <- 1 / (subset$se_outcome ^ 2)
  estimate <- sum(sw * subset$beta_exposure * subset$beta_outcome_harmonized) /
    sum(sw * subset$beta_exposure ^ 2)
  loo <- rbind(loo, data.frame(excluded_variant = keep$variant_id[i], ivw_beta = estimate))
}

summary_table <- data.frame(
  method = c("IVW", "MR_Egger_slope", "MR_Egger_intercept"),
  n_instruments = nrow(keep),
  estimate = c(ivw_beta, egger_coef["beta_exposure", "Estimate"], egger_coef["(Intercept)", "Estimate"]),
  se = c(ivw_se, egger_coef["beta_exposure", "Std. Error"], egger_coef["(Intercept)", "Std. Error"]),
  p = c(ivw_p, egger_coef["beta_exposure", "Pr(>|t|)"], egger_coef["(Intercept)", "Pr(>|t|)"]),
  q = c(q, NA, NA),
  q_df = c(q_df, NA, NA),
  q_p = c(q_p, NA, NA),
  interpretation_ceiling = c(
    "Causal evidence only if MR assumptions hold",
    "Pleiotropy-robust sensitivity estimate with lower precision",
    "Directional-pleiotropy diagnostic, not a causal effect"
  )
)

write.table(dat, file.path(outdir, "mr_harmonized.tsv"), sep = "\t", row.names = FALSE, quote = FALSE, na = "")
write.table(summary_table, file.path(outdir, "mr_summary.tsv"), sep = "\t", row.names = FALSE, quote = FALSE, na = "")
write.table(loo, file.path(outdir, "mr_leave_one_out.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
cat("MR complete:", nrow(keep), "harmonized instruments;", nrow(dat) - nrow(keep), "excluded\n")
