suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

# ----------------------------- Paths ------------------------------------------

base_work <- "/projects/MRC-IEU/research/projects/ieu3/p5/017/working"
ldsc_dir  <- "/user/work/sd20930/project_pa_mrpreg/ldsc"

exposure_file <- file.path(base_work, "data/MR-PREG/exposures_pa.txt")

outdir <- "/user/work/sd20930/project_pa_mrpreg/variant_relevance/results/steiger_confounder_filter"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

bmi_file <- file.path(ldsc_dir, "qc_bmil.sumstats.gz")
edu_file <- file.path(ldsc_dir, "qc_edu_attain.sumstats.gz")
bmd_file <- file.path(ldsc_dir, "qc_bmd.sumstats.gz")

# ----------------------------- Functions --------------------------------------

# Approximate R2 from Z and N:
# R2 = Z^2 / (Z^2 + N)
# adequate for Steiger-style comparison when only Z and N are available
calc_r2_from_z <- function(z, n) {
  z2 <- z^2
  z2 / (z2 + n)
}

read_confounder <- function(file, trait_name) {
  fread(file, data.table = FALSE) %>%
    rename(
      SNP = SNP,
      A1_confounder = A1,
      A2_confounder = A2,
      Z_confounder = Z,
      N_confounder = N
    ) %>%
    mutate(
      confounder = trait_name,
      Z_confounder = as.numeric(Z_confounder),
      N_confounder = as.numeric(N_confounder),
      P_confounder = 2 * pnorm(-abs(Z_confounder)),
      R2_confounder = calc_r2_from_z(Z_confounder, N_confounder)
    ) %>%
    select(
      SNP,
      confounder,
      A1_confounder,
      A2_confounder,
      Z_confounder,
      N_confounder,
      P_confounder,
      R2_confounder
    )
}

# ----------------------------- Exposure data ----------------------------------

exp_raw <- fread(exposure_file, data.table = FALSE)

if (!"se" %in% names(exp_raw)) exp_raw$se <- NA_real_

exp_raw <- exp_raw %>%
  mutate(
    SNP = as.character(SNP),
    Phenotype = as.character(Phenotype),
    effect_allele = toupper(as.character(effect_allele)),
    beta = as.numeric(beta),
    se = suppressWarnings(as.numeric(se)),
    pval = as.numeric(pval),
    eaf = as.numeric(eaf)
  ) %>%
  filter(effect_allele %in% c("A", "C", "G", "T")) %>%
  filter(!is.na(SNP), !is.na(Phenotype), !is.na(beta), !is.na(se), se > 0)

# Exposure R2 approximation from Z and N is unavailable unless N is in exposure file.
# Therefore use F/(F+N) if N exists; otherwise use Z^2 only for ranking.
# For Steiger comparison with confounder R2 we need an exposure N.
# Add N manually below if not in exposures_pa.txt.

if (!"N" %in% names(exp_raw)) {
  exp_raw <- exp_raw %>%
    mutate(
      N = case_when(
        Study == "Wang" ~ NA_real_,
        Study == "Klimentidis" ~ NA_real_,
        Study == "Qi" ~ NA_real_,
        Study == "Doherty" ~ NA_real_,
        Study == "Schoeler" ~ NA_real_,
        TRUE ~ NA_real_
      )
    )
}

# IMPORTANT:
# Fill in exposure GWAS N values here if known.
# If left as NA, script will still output Z and confounder R2, but not final Steiger flags.
# Example:
# Wang leisure MVPA N = XXXXX
# Klimentidis self-report/accelerometer N = XXXXX
# Qi accelerometer N = XXXXX
# Doherty accelerometer N = XXXXX

exp_raw <- exp_raw %>%
  mutate(
    N_exposure = as.numeric(N),
    Z_exposure = beta / se,
    P_exposure_check = 2 * pnorm(-abs(Z_exposure)),
    R2_exposure = if_else(
      !is.na(N_exposure),
      calc_r2_from_z(Z_exposure, N_exposure),
      NA_real_
    )
  )

# Restrict to final MR exposures of interest
exp_keep <- exp_raw %>%
  filter(Phenotype %in% c(
    "MVPA_leisure",
    "Vigorous_PA",
    "Sedentary_time",
    "Total_log_acceleration"
  ))

# ----------------------------- Confounder data --------------------------------

confounders <- bind_rows(
  read_confounder(bmi_file, "BMI"),
  read_confounder(edu_file, "Educational attainment"),
  read_confounder(bmd_file, "Bone mineral density")
)

# ----------------------------- Steiger comparison -----------------------------

steiger_dat <- exp_keep %>%
  select(
    Phenotype,
    SNP,
    Study,
    beta_exposure = beta,
    se_exposure = se,
    pval_exposure = pval,
    eaf_exposure = eaf,
    N_exposure,
    Z_exposure,
    R2_exposure
  ) %>%
  left_join(confounders, by = "SNP") %>%
  mutate(
    steiger_remove = if_else(
      !is.na(R2_exposure) & !is.na(R2_confounder) &
        R2_confounder > R2_exposure,
      TRUE,
      FALSE,
      missing = NA
    ),
    confounder_bonferroni = P_confounder < 0.05 / n(),
    nominal_confounder = P_confounder < 0.05
  )

# ----------------------------- Summaries --------------------------------------

removed_snps <- steiger_dat %>%
  filter(steiger_remove == TRUE) %>%
  distinct(Phenotype, SNP, confounder)

summary_dat <- steiger_dat %>%
  group_by(Phenotype, confounder) %>%
  summarise(
    n_snps_tested = n_distinct(SNP[!is.na(R2_confounder)]),
    n_steiger_removed = n_distinct(SNP[steiger_remove == TRUE], na.rm = TRUE),
    removed_snps = paste(unique(SNP[steiger_remove == TRUE]), collapse = ";"),
    n_nominal_confounder = n_distinct(SNP[nominal_confounder == TRUE], na.rm = TRUE),
    n_bonferroni_confounder = n_distinct(SNP[confounder_bonferroni == TRUE], na.rm = TRUE),
    .groups = "drop"
  )

# Filtered exposure files for rerunning MR
for (conf in unique(steiger_dat$confounder)) {

  snps_to_remove <- steiger_dat %>%
    filter(confounder == conf, steiger_remove == TRUE) %>%
    distinct(Phenotype, SNP)

  filtered <- exp_raw

  if (nrow(snps_to_remove) > 0) {
    filtered <- filtered %>%
      anti_join(snps_to_remove, by = c("Phenotype", "SNP"))
  }

  conf_stub <- conf %>%
    gsub(" ", "_", .) %>%
    gsub("[^A-Za-z0-9_]", "", .)

  fwrite(
    filtered,
    file.path(outdir, paste0("exposures_pa_steiger_filtered_", conf_stub, ".txt")),
    sep = "\t",
    quote = FALSE
  )
}

# Also write combined file removing SNPs flagged by any confounder
snps_remove_any <- steiger_dat %>%
  filter(steiger_remove == TRUE) %>%
  distinct(Phenotype, SNP)

exp_filtered_any <- exp_raw %>%
  anti_join(snps_remove_any, by = c("Phenotype", "SNP"))

fwrite(
  exp_filtered_any,
  file.path(outdir, "exposures_pa_steiger_filtered_any_confounder.txt"),
  sep = "\t",
  quote = FALSE
)

# Outputs
fwrite(
  steiger_dat,
  file.path(outdir, "steiger_confounder_filter_snp_level.tsv"),
  sep = "\t",
  quote = FALSE
)

fwrite(
  summary_dat,
  file.path(outdir, "steiger_confounder_filter_summary.tsv"),
  sep = "\t",
  quote = FALSE
)

fwrite(
  removed_snps,
  file.path(outdir, "steiger_confounder_removed_snps.tsv"),
  sep = "\t",
  quote = FALSE
)

print(summary_dat)
