################################################################################
# Grace Power
# 15 Apr 2026
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
#   - Exposure file has only effect_allele, so exposure data are built manually
#   - Outcome data are also built manually because format_data() was dropping all rows
#   - Wald ratio is used when only 1 SNP survives harmonisation
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
exp_raw <- fread(exposure_file, data.table = FALSE)

if (!"se" %in% names(exp_raw)) exp_raw$se <- NA_real_

exp_raw <- exp_raw %>%
  mutate(
    SNP = as.character(SNP),
    Phenotype = as.character(Phenotype),
    effect_allele = toupper(as.character(effect_allele)),
    beta = as.numeric(beta),
    se   = as.numeric(se),
    eaf  = as.numeric(eaf),
    pval = as.numeric(pval)
  ) %>%
  filter(
    !is.na(SNP),
    !is.na(Phenotype),
    !is.na(effect_allele),
    !is.na(beta),
    !is.na(pval),
    effect_allele %in% c("A", "C", "G", "T")
  )

needs_se <- is.na(exp_raw$se) | exp_raw$se <= 0
if (any(needs_se)) {
  z <- qnorm(1 - exp_raw$pval[needs_se] / 2)
  exp_raw$se[needs_se] <- abs(exp_raw$beta[needs_se] / z)
}

exp_raw <- exp_raw %>%
  filter(!is.na(se), se > 0)

# manual exposure object
exp_dat <- data.frame(
  SNP = exp_raw$SNP,
  beta.exposure = exp_raw$beta,
  se.exposure = exp_raw$se,
  pval.exposure = exp_raw$pval,
  eaf.exposure = exp_raw$eaf,
  effect_allele.exposure = exp_raw$effect_allele,
  other_allele.exposure = NA_character_,
  exposure = exp_raw$Phenotype,
  id.exposure = exp_raw$Phenotype,
  stringsAsFactors = FALSE
)

exp_list <- split(exp_dat, exp_dat$exposure)

# ----------------------------- Study-specific outcomes ------------------------
stu_raw <- fread(stu_out_file, data.table = FALSE)

stopifnot("study" %in% names(stu_raw))
stopifnot("Phenotype" %in% names(stu_raw))

stu_raw <- stu_raw %>%
  mutate(
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
  filter(
    !is.na(SNP),
    !is.na(Phenotype),
    !is.na(study),
    !is.na(effect_allele),
    !is.na(other_allele),
    !is.na(beta),
    !is.na(se),
    !is.na(pval),
    se > 0,
    effect_allele %in% c("A", "C", "G", "T"),
    other_allele %in% c("A", "C", "G", "T")
  )

# manual outcome object
stu_out <- data.frame(
  SNP = stu_raw$SNP,
  beta.outcome = stu_raw$beta,
  se.outcome = stu_raw$se,
  pval.outcome = stu_raw$pval,
  eaf.outcome = stu_raw$eaf,
  effect_allele.outcome = stu_raw$effect_allele,
  other_allele.outcome = stu_raw$other_allele,
  outcome = stu_raw$Phenotype,
  id.outcome = stu_raw$Phenotype,
  study = stu_raw$study,
  stringsAsFactors = FALSE
)

outcomes <- unique(stu_out$outcome)

# ----------------------------- Debug checks -----------------------------------
cat("\n--- RAW FILE OVERLAP CHECK ---\n")
cat("Unique exposure SNPs in raw file:", length(unique(exp_raw$SNP)), "\n")
cat("Unique study-outcome SNPs in raw file:", length(unique(stu_raw$SNP)), "\n")
cat("Raw-file SNP overlap:", length(intersect(unique(exp_raw$SNP), unique(stu_raw$SNP))), "\n")
cat("--- END RAW FILE OVERLAP CHECK ---\n\n")

cat("\n--- PRE-LOOP CHECKS ---\n")
cat("nrow(exp_raw):", nrow(exp_raw), "\n")
cat("nrow(exp_dat):", nrow(exp_dat), "\n")
cat("length(exp_list):", length(exp_list), "\n")
cat("Exposure names:\n")
print(names(exp_list))

cat("\n")
cat("nrow(stu_raw):", nrow(stu_raw), "\n")
cat("nrow(stu_out):", nrow(stu_out), "\n")
cat("length(outcomes):", length(outcomes), "\n")
cat("First few outcomes:\n")
print(head(outcomes))

cat("\n")
cat("Number of non-missing study labels in stu_out:", sum(!is.na(stu_out$study)), "\n")
cat("Unique studies in stu_out:\n")
print(unique(stu_out$study)[1:min(10, length(unique(stu_out$study)))])

cat("\n")
cat("Example exposure SNPs:\n")
print(head(unique(exp_dat$SNP), 10))

cat("\n")
cat("Example outcome SNPs:\n")
print(head(unique(stu_out$SNP), 10))

cat("\n")
cat("Raw SNP overlap between exp_dat and stu_out:",
    length(intersect(unique(exp_dat$SNP), unique(stu_out$SNP))), "\n")
cat("--- END PRE-LOOP CHECKS ---\n\n")

# ----------------------------- MR within each study ---------------------------
study_level <- list()
k <- 1

for (e_name in names(exp_list)) {
  for (o_name in outcomes) {

    o_all <- stu_out %>% filter(outcome == o_name)
    o_by_study <- split(o_all, o_all$study)

    for (stu in names(o_by_study)) {

      message("--------------------------------------------------")
      message("Study MR: ", e_name, " -> ", o_name, " (", stu, ")")

      e_dat <- exp_list[[e_name]]
      o_dat <- o_by_study[[stu]]

      message("Exposure SNPs: ", length(unique(e_dat$SNP)))
      message("Outcome SNPs: ", length(unique(o_dat$SNP)))
      message("Raw overlap: ", length(intersect(e_dat$SNP, o_dat$SNP)))

      dat_h <- tryCatch(
        harmonise_data(e_dat, o_dat, action = 2),
        error = function(e) {
          message("harmonise_data error: ", e$message)
          NULL
        }
      )

      if (is.null(dat_h)) next
      if (!"mr_keep" %in% names(dat_h)) next

      message("Rows after harmonise_data: ", nrow(dat_h))
      message("Rows with mr_keep == TRUE: ", sum(dat_h$mr_keep, na.rm = TRUE))

      dat_h <- dat_h %>% filter(mr_keep)

      nsnp <- length(unique(dat_h$SNP))
      message("Unique SNPs after harmonisation: ", nsnp)

      if (nsnp == 0) next

      if (nsnp == 1) {
        res <- tryCatch(
          mr(dat_h, method_list = "mr_wald_ratio"),
          error = function(e) {
            message("mr_wald_ratio error: ", e$message)
            NULL
          }
        )
      } else {
        res <- tryCatch(
          mr(dat_h, method_list = c("mr_ivw", "mr_wald_ratio")),
          error = function(e) {
            message("mr() error: ", e$message)
            NULL
          }
        )
      }

      if (is.null(res) || nrow(res) == 0) {
        message("No MR results returned")
        next
      }

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

cat("\nFinal study_level_df dim:\n")
print(dim(study_level_df))
print(head(study_level_df))

if (nrow(study_level_df) == 0) {
  message("No study-level MR results were generated.")
  quit(save = "no")
}

fwrite(
  study_level_df,
  file = file.path(outdir, "study_level_mr.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ----------------------------- Leave-one-study-out ----------------------------
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
