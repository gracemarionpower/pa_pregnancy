################################################################################
# Steiger-filtered MR: meta-analysis outcomes only
# Removes rs1160545 and rs7613360 from MVPA_leisure only
################################################################################

rm(list = ls())

suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(data.table)
  library(dplyr)
})

# ----------------------------- Paths ------------------------------------------
base_data <- "/projects/MRC-IEU/research/projects/ieu3/p5/017/working/data/MR-PREG"

exposure_file   <- file.path(base_data, "exposures_pa.txt")
ma_outcome_file <- file.path(base_data, "ma_out_dat.txt")

outdir <- "/projects/MRC-IEU/research/projects/ieu3/p5/017/working/results/mr"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ----------------------------- Exposures --------------------------------------
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

# Remove Steiger-flagged education SNPs from MVPA_leisure only
exp_raw <- exp_raw %>%
  filter(!(Phenotype == "MVPA_leisure" & SNP %in% c("rs1160545", "rs7613360")))

# Keep only standard single-base alleles
exp_raw <- exp_raw %>%
  filter(effect_allele %in% c("A", "C", "G", "T"))

# Reconstruct missing SEs
needs_se <- is.na(exp_raw$se) | exp_raw$se <= 0
if (any(needs_se)) {
  z <- qnorm(1 - exp_raw$pval[needs_se] / 2)
  exp_raw$se[needs_se] <- abs(exp_raw$beta[needs_se] / z)
}

exp_raw <- exp_raw %>%
  filter(
    !is.na(SNP),
    !is.na(Phenotype),
    !is.na(effect_allele),
    !is.na(beta),
    !is.na(se),
    !is.na(pval)
  )

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

methods <- c(
  "mr_ivw",
  "mr_egger_regression",
  "mr_weighted_median",
  "mr_weighted_mode"
)

# ==============================================================================
# Meta-analysis outcomes only
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

    message("Steiger-filtered MA MR: ", e_name, " -> ", o_name)

    dat_h <- harmonise_data(exp_list[[e_name]], ma_list[[o_name]], action = 2)

    if (!"mr_keep" %in% names(dat_h)) next
    dat_h <- dat_h %>% filter(mr_keep)

    nsnp <- length(unique(dat_h$SNP))
    if (nsnp == 0) next

    res <- tryCatch(
      mr(dat_h, method_list = methods),
      error = function(e) NULL
    )

    if (is.null(res) || nrow(res) == 0) {
      res <- tryCatch(mr(dat_h), error = function(e) NULL)
    }

    if (is.null(res) || nrow(res) == 0) next

    ma_res[[k]] <- res %>%
      transmute(
        exposure = exposure,
        outcome = outcome,
        analysis = "meta_analysis_steiger_filtered_education",
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

ma_final <- bind_rows(ma_res)

outfile <- file.path(
  outdir,
  "mr_ma_results_steiger_filtered_education.tsv"
)

fwrite(
  ma_final,
  file = outfile,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("Finished.")
message("Wrote: ", outfile)
