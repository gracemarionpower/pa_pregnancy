################################################################################
# Grace Power
# 25 Feb 2026
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
#
# Inputs:
#   data/exposures/exposures_pa.txt
#   data/outcomes/stu_out_dat.txt
#
# Outputs:
#   results/mr_study_loo/study_level_mr.tsv
#   results/mr_study_loo/study_loo_ivw.tsv
################################################################################

rm(list = ls())

suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(data.table)
  library(dplyr)
  library(metafor)
})

# ----------------------------- Paths ------------------------------------------
exposure_file <- "data/exposures/exposures_pa.txt"
stu_out_file  <- "data/outcomes/stu_out_dat.txt"

outdir <- "results/mr_study_loo"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ----------------------------- Exposures --------------------------------------
exp_raw <- fread(exposure_file, data.table = FALSE)

if (!"other_allele" %in% names(exp_raw)) exp_raw$other_allele <- NA_character_

exp_raw <- exp_raw %>%
  mutate(
    SNP = as.character(SNP),
    Phenotype = as.character(Phenotype),
    effect_allele = as.character(effect_allele),
    other_allele = as.character(other_allele)
  ) %>%
  filter(!is.na(SNP), !is.na(Phenotype), !is.na(effect_allele),
         !is.na(beta), !is.na(se), !is.na(pval))

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

# Preserve study column and match rows by SNP + phenotype
stu_out$study <- stu_raw$study[match(paste(stu_out$SNP, stu_out$outcome),
                                    paste(stu_raw$SNP, stu_raw$Phenotype))]

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

      dat_h <- harmonise_data(exp_list[[e_name]], o_by_study[[stu]], action = 2) %>%
        filter(mr_keep)

      nsnp <- length(unique(dat_h$SNP))
      if (nsnp == 0) next

      res <- tryCatch(mr(dat_h, method_list = methods), error = function(e) NULL)
      if (is.null(res) || nrow(res) == 0) res <- tryCatch(mr(dat_h), error = function(e) NULL)
      if (is.null(res) || nrow(res) == 0) next

      study_level[[k]] <- res %>%
        transmute(
          exposure = exposure,
          outcome = outcome,
          study = stu,
          method = method,
          nsnp = nsnp,
          b = b, se = se, pval = pval
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

      fit <- metafor::rma.uni(yi = dt2$b, sei = dt2$se, method = "FE")

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
