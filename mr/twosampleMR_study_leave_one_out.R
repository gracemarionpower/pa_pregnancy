################################################################################
# Grace Power
# Leave-one-study-out MR using stu_out_dat.txt
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

# reconstruct se if needed
needs_se <- is.na(exp_raw$se) | exp_raw$se <= 0
if (any(needs_se)) {
  z <- qnorm(1 - exp_raw$pval[needs_se] / 2)
  exp_raw$se[needs_se] <- abs(exp_raw$beta[needs_se] / z)
}

exp_raw <- exp_raw %>%
  filter(!is.na(se), se > 0)

# Build exposure data manually instead of format_data()
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

methods <- c("mr_ivw", "mr_wald_ratio")

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
    !is.na(SNP), !is.na(Phenotype), !is.na(study),
    !is.na(effect_allele), !is.na(other_allele),
    !is.na(beta), !is.na(se), !is.na(pval),
    effect_allele %in% c("A", "C", "G", "T"),
    other_allele %in% c("A", "C", "G", "T")
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

# reattach study by matching on SNP + phenotype + beta + se + pval
raw_key <- paste(
  stu_raw$SNP,
  stu_raw$Phenotype,
  signif(stu_raw$beta, 12),
  signif(stu_raw$se, 12),
  signif(stu_raw$pval, 12),
  sep = "||"
)

fmt_key <- paste(
  stu_out$SNP,
  stu_out$outcome,
  signif(stu_out$beta.outcome, 12),
  signif(stu_out$se.outcome, 12),
  signif(stu_out$pval.outcome, 12),
  sep = "||"
)

stu_out$study <- stu_raw$study[match(fmt_key, raw_key)]
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

      message("Rows after harmonise_data: ", nrow(dat_h))

      if (!"mr_keep" %in% names(dat_h)) {
        message("No mr_keep column returned")
        next
      }

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

# ----------------------------- Leave-one-study-out ----------------------------
if (nrow(study_level_df) == 0) {
  message("No study-level MR results were generated.")
  fwrite(data.frame(), file = file.path(outdir, "study_loo_ivw.tsv"), sep = "\t")
  quit(save = "no")
}

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
