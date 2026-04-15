################################################################################
# Grace Power
# 15 Feb 2026
# Leave-one-study-out MR using stu_out_dat.txt
#
# Rationale:
#   The study-specific outcome file is used to estimate MR separately within each
#   study. I then carry out leave-one-study-out by recomputing a fixed-effect
#   meta-analysis of the remaining study-specific MR estimates.
#
# Notes:
#   - Study-level LOO is implemented for IVW estimates.
#   - Egger/median/mode are still computed per study (where possible) and written
#     out, but the LOO aggregation uses IVW as the primary estimator.
#   - Exposure SE is reconstructed from beta and two-sided p-values where missing/0
#   - Exposure alleles are restricted to A/C/G/T to avoid format_data() excluding rows
#   - Exposure other_allele is set to NA to match the working main MR script
#
# Inputs:
#   /projects/MRC-IEU/research/projects/ieu3/p5/017/working/data/MR-PREG/exposures_pa.txt
#   /projects/MRC-IEU/research/projects/ieu3/p5/017/working/data/MR-PREG/stu_out_dat.txt
#
# Outputs:
#   /projects/MRC-IEU/research/projects/ieu3/p5/017/working/results/mr_study_loo/study_level_mr.tsv
#   /projects/MRC-IEU/research/projects/ieu3/p5/017/working/results/mr_study_loo/study_loo_ivw.tsv
################################################################################

rm(list = ls())

suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(data.table)
  library(dplyr)
  library(metafor)
})

# ----------------------------- Paths ------------------------------------------
base_dir  <- "/projects/MRC-IEU/research/projects/ieu3/p5/017"
base_data <- file.path(base_dir, "working/data/MR-PREG")

exposure_file <- file.path(base_data, "exposures_pa.txt")
stu_out_file  <- file.path(base_data, "stu_out_dat.txt")

outdir <- file.path(base_dir, "working/results/mr_study_loo")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

cat("Exposure file:", exposure_file, "\n")
cat("Exists:", file.exists(exposure_file), "\n")
cat("Study outcome file:", stu_out_file, "\n")
cat("Exists:", file.exists(stu_out_file), "\n")

# ----------------------------- Exposures --------------------------------------
# One row per SNP per PA phenotype.
# SE is reconstructed from beta and pval if missing/0.
# Exposure alleles are restricted to A/C/G/T.
# Exposure other_allele is set to NA to avoid format_data dropping rows.
exp_raw <- fread(exposure_file, data.table = FALSE)

if (!"se" %in% names(exp_raw)) exp_raw$se <- NA_real_

exp_raw <- exp_raw %>%
  mutate(
    SNP = as.character(SNP),
    Phenotype = as.character(Phenotype),
    effect_allele = toupper(as.character(effect_allele)),
    beta = as.numeric(beta),
    eaf  = as.numeric(eaf),
    pval = as.numeric(pval),
    se   = suppressWarnings(as.numeric(se))
  )

# Keep only standard single-base alleles for exposure coding
exp_raw <- exp_raw %>%
  filter(effect_allele %in% c("A", "C", "G", "T"))

# Reconstruct SE where missing or unusable
needs_se <- is.na(exp_raw$se) | exp_raw$se <= 0
if (any(needs_se)) {
  z <- qnorm(1 - exp_raw$pval[needs_se] / 2)
  exp_raw$se[needs_se] <- abs(exp_raw$beta[needs_se] / z)
}

# Drop unusable rows
exp_raw <- exp_raw %>%
  filter(!is.na(SNP), !is.na(Phenotype), !is.na(effect_allele),
         !is.na(beta), !is.na(se), !is.na(pval))

# IMPORTANT: blank other_allele on exposure side, as in your working script
exp_raw$other_allele <- NA_character_

exp_dat <- format_data(
  exp_raw,
  type = "exposure",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  eaf_col = "eaf",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  pval_col = "pval",
  phenotype_col = "Phenotype"
)

exp_list <- split(exp_dat, exp_dat$exposure)

methods <- c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode")

# ----------------------------- Study-specific outcomes ------------------------
stu_raw <- fread(stu_out_file, data.table = FALSE)

stopifnot("study" %in% names(stu_raw))
stopifnot("Phenotype" %in% names(stu_raw))

stu_raw <- stu_raw %>%
  mutate(
    row_id = seq_len(n()),
    SNP = as.character(SNP),
    Phenotype = as.character(Phenotype),
    study = as.character(study),
    effect_allele = toupper(as.character(effect_allele)),
    other_allele  = toupper(as.character(other_allele)),
    beta = as.numeric(beta),
    se   = as.numeric(se),
    eaf  = as.numeric(eaf),
    pval = as.numeric(pval)
  ) %>%
  filter(!is.na(SNP), !is.na(Phenotype), !is.na(study),
         !is.na(effect_allele), !is.na(other_allele),
         !is.na(beta), !is.na(se), !is.na(pval))

# Keep only rows with alleles format_data can use
stu_raw <- stu_raw %>%
  filter(
    effect_allele %in% c("A","C","G","T","D","I") |
      grepl("^[ACGT]+$", effect_allele),
    other_allele %in% c("A","C","G","T","D","I") |
      grepl("^[ACGT]+$", other_allele)
  )

stu_out <- format_data(
  stu_raw,
  type = "outcome",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  eaf_col = "eaf",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  pval_col = "pval",
  phenotype_col = "Phenotype"
)

# Reattach study by matching retained formatted rows back to source rows
# using SNP + phenotype + beta + se + pval
stu_key_raw <- paste(
  stu_raw$SNP,
  stu_raw$Phenotype,
  signif(stu_raw$beta, 12),
  signif(stu_raw$se, 12),
  signif(stu_raw$pval, 12),
  sep = "||"
)

stu_key_fmt <- paste(
  stu_out$SNP,
  stu_out$outcome,
  signif(stu_out$beta.outcome, 12),
  signif(stu_out$se.outcome, 12),
  signif(stu_out$pval.outcome, 12),
  sep = "||"
)

m <- match(stu_key_fmt, stu_key_raw)
stu_out$study <- stu_raw$study[m]

# Drop any rows where study could not be recovered
stu_out <- stu_out %>% filter(!is.na(study))

outcomes <- unique(stu_out$outcome)

# ----------------------------- MR within each study ---------------------------
study_level <- list()
k <- 1

for (e_name in names(exp_list)) {
  for (o_name in outcomes) {

    o_all <- stu_out %>% filter(outcome == o_name)
    o_by_study <- split(o_all, o_all$study)

    for (stu in names(o_by_study)) {

      message("Study MR: ", e_name, " -> ", o_name, " (", stu, ")")

      dat_h <- harmonise_data(exp_list[[e_name]], o_by_study[[stu]], action = 2)

      if (!"mr_keep" %in% names(dat_h)) next
      dat_h <- dat_h %>% filter(mr_keep)

      nsnp <- length(unique(dat_h$SNP))
      if (nsnp == 0) next

      res <- tryCatch(mr(dat_h, method_list = methods), error = function(e) NULL)
      if (is.null(res) || nrow(res) == 0) {
        res <- tryCatch(mr(dat_h), error = function(e) NULL)
      }
      if (is.null(res) || nrow(res) == 0) next

      study_level[[k]] <- res %>%
        transmute(
          exposure = exposure,
          outcome = outcome,
          study = stu,
          method = method,
          nsnp = nsnp,
          b = b,
          se = se,
          pval = pval,
          OR  = exp(b),
          LCI = exp(b - 1.96 * se),
          UCI = exp(b + 1.96 * se)
        )
      k <- k + 1
    }
  }
}

study_level_df <- bind_rows(study_level)

fwrite(
  study_level_df,
  file = file.path(outdir, "study_level_mr.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ----------------------------- Leave-one-study-out (IVW FE meta) --------------
ivw_df <- study_level_df %>%
  filter(method == "Inverse variance weighted") %>%
  filter(!is.na(b), !is.na(se))

loo_rows <- list()
kk <- 1

for (e_name in unique(ivw_df$exposure)) {
  for (o_name in unique(ivw_df$outcome)) {

    dt <- ivw_df %>% filter(exposure == e_name, outcome == o_name)
    if (nrow(dt) < 2) next

    studies <- unique(dt$study)

    for (drop_stu in studies) {

      dt2 <- dt %>% filter(study != drop_stu)
      if (nrow(dt2) < 1) next

      fit <- tryCatch(
        metafor::rma.uni(yi = dt2$b, sei = dt2$se, method = "FE"),
        error = function(e) NULL
      )
      if (is.null(fit)) next

      loo_rows[[kk]] <- data.frame(
        exposure = e_name,
        outcome = o_name,
        removed_study = drop_stu,
        k_remaining = nrow(dt2),
        b = as.numeric(fit$b),
        se = as.numeric(fit$se),
        pval = as.numeric(fit$pval),
        OR  = exp(as.numeric(fit$b)),
        LCI = exp(as.numeric(fit$b) - 1.96 * as.numeric(fit$se)),
        UCI = exp(as.numeric(fit$b) + 1.96 * as.numeric(fit$se)),
        stringsAsFactors = FALSE
      )
      kk <- kk + 1
    }
  }
}

study_loo_ivw <- bind_rows(loo_rows)

fwrite(
  study_loo_ivw,
  file = file.path(outdir, "study_loo_ivw.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("Finished.")
message("Wrote: ", file.path(outdir, "study_level_mr.tsv"))
message("Wrote: ", file.path(outdir, "study_loo_ivw.tsv"))
