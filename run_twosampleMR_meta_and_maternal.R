################################################################################
# Grace Power
# MR: all PA exposures x all outcomes
# Runs:
#   (A) Meta-analysis outcomes (ma_out_dat.txt)
#   (B) Maternal effects (unadjusted) and maternal effects adjusted for fetal (mat_donuts)
#       using duos_out_dat.txt + trios_out_dat.txt
#
# Methods (where possible):
#   IVW, MR-Egger, weighted median, weighted mode
#
# Inputs:
#   data/exposures/exposures_pa.txt
#   data/outcomes/ma_out_dat.txt
#   data/outcomes/duos_out_dat.txt
#   data/outcomes/trios_out_dat.txt
#
# Outputs:
#   results/mr_all/mr_ma_results.tsv
#   results/mr_all/mr_maternal_results.tsv
################################################################################

rm(list = ls())

suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(data.table)
  library(dplyr)
  library(ggplot2)
})

# ----------------------------- Paths ------------------------------------------
base_dir <- "/projects/MRC-IEU/research/projects/ieu3/p5/017/working/data/MR-PREG"

exposure_file   <- file.path(base_dir, "exposures_pa.txt")
ma_outcome_file <- file.path(base_dir, "ma_out_dat.txt")
duos_file       <- file.path(base_dir, "duos_out_dat.txt")
trios_file      <- file.path(base_dir, "trios_out_dat.txt")

outdir <- file.path(base_dir, "results", "mr_all")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(outdir, "plots"), showWarnings = FALSE, recursive = TRUE)

# ----------------------------- Exposures --------------------------------------
# The exposure table is one row per SNP per PA phenotype.
# SE is reconstructed from beta and pval if missing or if rounded to zero.
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

# Reconstruct SE where missing or unusable (e.g., 0.00 from rounding)
needs_se <- is.na(exp_raw$se) | exp_raw$se <= 0
if (any(needs_se)) {
  z <- qnorm(1 - exp_raw$pval[needs_se] / 2)
  exp_raw$se[needs_se] <- abs(exp_raw$beta[needs_se] / z)
}

# Drop any rows that still cannot be used
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
# (A) MR against meta-analysis outcomes (ma_out_dat.txt)
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

    # If harmonisation returns an object without mr_keep (can happen when 0 SNPs survive
    # or allele matching fails), skip this exposure–outcome pair cleanly
    if (!"mr_keep" %in% names(dat_h)) {
      message("Skipping (no mr_keep returned): ", e_name, " -> ", o_name)
      next
    }

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

    if (nsnp >= 2) {
      p <- mr_scatter_plot(res, dat_h)
      ggsave(
        filename = file.path(outdir, "plots",
                             paste0("scatter_ma_", make.names(e_name), "_", make.names(o_name), ".png")),
        plot = p[[1]],
        width = 6.5, height = 5, dpi = 250
      )
    }
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
# (B) MR against maternal outcomes: unadjusted (mat) and fetal-adjusted (mat_donuts)
#      using duos_out_dat.txt + trios_out_dat.txt
# ==============================================================================

duos_raw <- fread(duos_file, data.table = FALSE)
trios_raw <- fread(trios_file, data.table = FALSE)

wlm_pool <- bind_rows(duos_raw, trios_raw)

maternal_specs <- data.frame(
  analysis = c("maternal_unadjusted", "maternal_fetal_adjusted"),
  beta_col = c("beta_mat", "beta_mat_donuts"),
  se_col   = c("se_mat",   "se_mat_donuts"),
  p_col    = c("p_mat",    "p_mat_donuts"),
  stringsAsFactors = FALSE
)

maternal_out_list <- list()

for (i in seq_len(nrow(maternal_specs))) {

  spec <- maternal_specs[i, ]

  tmp <- wlm_pool %>%
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

maternal_split <- split(maternal_out, interaction(maternal_out$outcome, maternal_out$analysis, drop = TRUE))

maternal_res <- list()
kk <- 1

for (e_name in names(exp_list)) {
  for (key in names(maternal_split)) {

    o_dat <- maternal_split[[key]]
    o_name <- unique(o_dat$outcome)
    a_name <- unique(o_dat$analysis)

    message("Maternal MR: ", e_name, " -> ", o_name, " [", a_name, "]")

    dat_h <- harmonise_data(exp_list[[e_name]], o_dat, action = 2) %>%
      filter(mr_keep)

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

    if (nsnp >= 2) {
      p <- mr_scatter_plot(res, dat_h)
      ggsave(
        filename = file.path(outdir, "plots",
                             paste0("scatter_", a_name, "_", make.names(e_name), "_", make.names(o_name), ".png")),
        plot = p[[1]],
        width = 6.5, height = 5, dpi = 250
      )
    }
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

message("Finished.")
message("Wrote: ", file.path(outdir, "mr_ma_results.tsv"))
message("Wrote: ", file.path(outdir, "mr_maternal_results.tsv"))
