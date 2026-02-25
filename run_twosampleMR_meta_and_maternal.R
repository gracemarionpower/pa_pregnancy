################################################################################
# Grace Power
# Plain MR analyses: all PA exposures x all outcomes
#
# What I run
#   (A) MR of each PA exposure phenotype against each meta-analysis outcome (ma_out_dat.txt)
#   (B) MR of each PA exposure phenotype against maternal outcome effects:
#         - maternal_unadjusted (duos+trios): beta_mat
#         - maternal_fetal_adjusted (duos only): beta_mat_donuts
#         - maternal_fetal_paternal_adjusted (trios only): beta_mat_donuts
#
# Methods (where possible)
#   IVW, MR-Egger, weighted median, weighted mode
#
# Notes
#   - Exposure SE is reconstructed from beta and two-sided p-values where missing/0
#   - I only write MR result tables; no plots or SNP-level outputs
#
# Inputs
#   /projects/MRC-IEU/research/projects/ieu3/p5/017/working/data/MR-PREG/exposures_pa.txt
#   /projects/MRC-IEU/research/projects/ieu3/p5/017/working/data/MR-PREG/ma_out_dat.txt
#   /projects/MRC-IEU/research/projects/ieu3/p5/017/working/data/MR-PREG/duos_out_dat.txt
#   /projects/MRC-IEU/research/projects/ieu3/p5/017/working/data/MR-PREG/trios_out_dat.txt
#
# Outputs
#   /projects/.../MR-PREG/results/mr_all/mr_ma_results.tsv
#   /projects/.../MR-PREG/results/mr_all/mr_maternal_results.tsv
################################################################################

rm(list = ls())

suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(data.table)
  library(dplyr)
})

# ----------------------------- Paths ------------------------------------------
base_dir <- "/projects/MRC-IEU/research/projects/ieu3/p5/017/working/data/MR-PREG"

exposure_file   <- file.path(base_dir, "exposures_pa.txt")
ma_outcome_file <- file.path(base_dir, "ma_out_dat.txt")
duos_file       <- file.path(base_dir, "duos_out_dat.txt")
trios_file      <- file.path(base_dir, "trios_out_dat.txt")

outdir <- file.path(base_dir, "results", "mr_all")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ----------------------------- Exposures --------------------------------------
# One row per SNP per PA phenotype. SE is reconstructed from beta and pval if missing/0.
exp_raw <- fread(exposure_file, data.table = FALSE)

if (!"other_allele" %in% names(exp_raw)) exp_raw$other_allele <- NA_character_

exp_raw <- exp_raw %>%
  mutate(
    SNP = as.character(SNP),
    Phenotype = as.character(Phenotype),
    effect_allele = as.character(effect_allele),
    other_allele = as.character(other_allele),
    beta = as.numeric(beta),
    eaf = as.numeric(eaf),
    pval = as.numeric(pval),
    se = suppressWarnings(as.numeric(se))
  )

needs_se <- is.na(exp_raw$se) | exp_raw$se <= 0
if (any(needs_se)) {
  z <- qnorm(1 - exp_raw$pval[needs_se] / 2)
  exp_raw$se[needs_se] <- abs(exp_raw$beta[needs_se] / z)
}

exp_raw <- exp_raw %>%
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

# ==============================================================================
# (A) Meta-analysis outcomes: ma_out_dat.txt
# ==============================================================================

ma_raw <- fread(ma_outcome_file, data.table = FALSE)

ma_out <- format_data(
  ma_raw,
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

ma_list <- split(ma_out, ma_out$outcome)

ma_res <- list()
k <- 1

for (e_name in names(exp_list)) {
  for (o_name in names(ma_list)) {

    message("MA MR: ", e_name, " -> ", o_name)

    dat_h <- harmonise_data(exp_list[[e_name]], ma_list[[o_name]], action = 2)

    if (!"mr_keep" %in% names(dat_h)) next
    dat_h <- dat_h %>% filter(mr_keep)

    nsnp <- length(unique(dat_h$SNP))
    if (nsnp == 0) next

    res <- tryCatch(mr(dat_h, method_list = methods), error = function(e) NULL)
    if (is.null(res) || nrow(res) == 0) res <- tryCatch(mr(dat_h), error = function(e) NULL)
    if (is.null(res) || nrow(res) == 0) next

    ma_res[[k]] <- res %>%
      transmute(
        exposure = exposure,
        outcome = outcome,
        analysis = "meta_analysis",
        method = method,
        nsnp = nsnp,
        b = b, se = se, pval = pval,
        OR  = exp(b),
        LCI = exp(b - 1.96 * se),
        UCI = exp(b + 1.96 * se)
      )
    k <- k + 1
  }
}

ma_final <- bind_rows(ma_res)

fwrite(
  ma_final,
  file = file.path(outdir, "mr_ma_results.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ==============================================================================
# (B) Maternal outcomes: duos + trios
# ==============================================================================

duos_raw  <- fread(duos_file, data.table = FALSE)
trios_raw <- fread(trios_file, data.table = FALSE)

maternal_specs <- data.frame(
  analysis = c(
    "maternal_unadjusted",
    "maternal_fetal_adjusted",
    "maternal_fetal_paternal_adjusted"
  ),
  source = c("both", "duos", "trios"),
  beta_col = c("beta_mat", "beta_mat_donuts", "beta_mat_donuts"),
  se_col   = c("se_mat",   "se_mat_donuts",   "se_mat_donuts"),
  p_col    = c("p_mat",    "p_mat_donuts",    "p_mat_donuts"),
  stringsAsFactors = FALSE
)

maternal_out_list <- list()

for (i in seq_len(nrow(maternal_specs))) {

  spec <- maternal_specs[i, ]

  tmp <- switch(
    spec$source,
    both  = bind_rows(duos_raw, trios_raw),
    duos  = duos_raw,
    trios = trios_raw
  )

  tmp <- tmp %>%
    filter(!is.na(.data[[spec$se_col]])) %>%
    mutate(Phenotype = as.character(Phenotype))

  out_tmp <- format_data(
    tmp,
    type = "outcome",
    snp_col = "SNP",
    beta_col = spec$beta_col,
    se_col   = spec$se_col,
    pval_col = spec$p_col,
    eaf_col  = "eaf_mat",
    effect_allele_col = "effect_allele",
    other_allele_col  = "other_allele",
    phenotype_col = "Phenotype"
  )

  out_tmp$analysis <- spec$analysis
  maternal_out_list[[i]] <- out_tmp
}

maternal_out <- bind_rows(maternal_out_list)

maternal_split <- split(
  maternal_out,
  interaction(maternal_out$outcome, maternal_out$analysis, drop = TRUE)
)

maternal_res <- list()
kk <- 1

for (e_name in names(exp_list)) {
  for (key in names(maternal_split)) {

    o_dat <- maternal_split[[key]]
    a_name <- unique(o_dat$analysis)

    message("Maternal MR: ", e_name, " -> ", unique(o_dat$outcome), " [", a_name, "]")

    dat_h <- harmonise_data(exp_list[[e_name]], o_dat, action = 2)

    if (!"mr_keep" %in% names(dat_h)) next
    dat_h <- dat_h %>% filter(mr_keep)

    nsnp <- length(unique(dat_h$SNP))
    if (nsnp == 0) next

    res <- tryCatch(mr(dat_h, method_list = methods), error = function(e) NULL)
    if (is.null(res) || nrow(res) == 0) res <- tryCatch(mr(dat_h), error = function(e) NULL)
    if (is.null(res) || nrow(res) == 0) next

    maternal_res[[kk]] <- res %>%
      transmute(
        exposure = exposure,
        outcome = outcome,
        analysis = a_name,
        method = method,
        nsnp = nsnp,
        b = b, se = se, pval = pval,
        OR  = exp(b),
        LCI = exp(b - 1.96 * se),
        UCI = exp(b + 1.96 * se)
      )
    kk <- kk + 1
  }
}

maternal_final <- bind_rows(maternal_res)

fwrite(
  maternal_final,
  file = file.path(outdir, "mr_maternal_results.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("Finished.")
message("Wrote: ", file.path(outdir, "mr_ma_results.tsv"))
message("Wrote: ", file.path(outdir, "mr_maternal_results.tsv"))
