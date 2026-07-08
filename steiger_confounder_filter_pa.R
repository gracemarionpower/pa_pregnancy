suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

# Paths
base_work <- "/projects/MRC-IEU/research/projects/ieu3/p5/017/working"
ldsc_dir  <- "/user/work/sd20930/project_pa_mrpreg/ldsc"

exposure_file <- file.path(base_work, "data/MR-PREG/exposures_pa.txt")

outdir <- "/user/work/sd20930/project_pa_mrpreg/variant_relevance/results/steiger_confounder_filter"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

bmi_file <- file.path(ldsc_dir, "ieu-b-40_summary_bmi.tsv")
edu_file <- file.path(ldsc_dir, "GCST90029013_edu_attain_summary.tsv")
bmd_file <- file.path(ldsc_dir, "GCST005348.sumstats_bmd.withN.tsv")

# Add exposure GWAS sample sizes here
exposure_N <- tibble::tribble(
  ~Phenotype, ~N_exposure,
  "MVPA_leisure", 606820,
  "Vigorous_PA", 261055,
  "Sedentary_time", 91105,
  "Total_log_acceleration", 88411
)

calc_r2 <- function(beta, se, eaf, n) {
  num <- 2 * (beta^2) * eaf * (1 - eaf)
  den <- num + 2 * n * eaf * (1 - eaf) * (se^2)
  num / den
}

read_confounder <- function(file, trait_name, default_n = NA_real_) {
  x <- fread(file, data.table = FALSE)

  if (!"N" %in% names(x)) {
    x$N <- default_n
  }

  x %>%
    transmute(
      SNP = as.character(SNP),
      confounder = trait_name,
      A1_confounder = A1,
      A2_confounder = A2,
      beta_confounder = as.numeric(BETA),
      se_confounder = as.numeric(SE),
      p_confounder = as.numeric(P),
      eaf_confounder = suppressWarnings(as.numeric(EAF)),
      N_confounder = as.numeric(N),
      R2_confounder = calc_r2(
        beta_confounder,
        se_confounder,
        eaf_confounder,
        N_confounder
      )
    )
}

exp_raw <- fread(exposure_file, data.table = FALSE)

exp_dat <- exp_raw %>%
  mutate(
    SNP = as.character(SNP),
    Phenotype = as.character(Phenotype),
    Study = as.character(Study),
    beta_exposure = as.numeric(beta),
    se_exposure = as.numeric(se),
    pval_exposure = as.numeric(pval),
    eaf_exposure = as.numeric(eaf)
  ) %>%
  left_join(exposure_N, by = "Phenotype") %>%
  filter(
    Phenotype %in% c(
      "MVPA_leisure",
      "Vigorous_PA",
      "Sedentary_time",
      "Total_log_acceleration"
    )
  ) %>%
  mutate(
    R2_exposure = calc_r2(
      beta_exposure,
      se_exposure,
      eaf_exposure,
      N_exposure
    )
  )

confounders <- bind_rows(
  read_confounder(bmi_file, "BMI", default_n = 681275),
  read_confounder(edu_file, "Educational attainment", default_n = 461457),
  read_confounder(bmd_file, "Bone mineral density")
)

steiger_dat <- exp_dat %>%
  select(
    Phenotype,
    SNP,
    Study,
    beta_exposure,
    se_exposure,
    pval_exposure,
    eaf_exposure,
    N_exposure,
    R2_exposure
  ) %>%
  left_join(confounders, by = "SNP") %>%
  mutate(
    steiger_remove = R2_confounder > R2_exposure,
    nominal_confounder = p_confounder < 0.05,
    bonferroni_confounder = p_confounder < 0.05 / sum(!is.na(p_confounder))
  )

summary_dat <- steiger_dat %>%
  group_by(Phenotype, confounder) %>%
  summarise(
    n_snps_tested = n_distinct(SNP[!is.na(R2_confounder)]),
    n_steiger_removed = n_distinct(SNP[steiger_remove == TRUE], na.rm = TRUE),
    removed_snps = paste(unique(SNP[steiger_remove == TRUE]), collapse = ";"),
    n_nominal_confounder = n_distinct(SNP[nominal_confounder == TRUE], na.rm = TRUE),
    n_bonferroni_confounder = n_distinct(SNP[bonferroni_confounder == TRUE], na.rm = TRUE),
    .groups = "drop"
  )

removed_snps <- steiger_dat %>%
  filter(steiger_remove == TRUE) %>%
  distinct(Phenotype, SNP, confounder)

# Write SNP-level and summary files
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

# Create filtered exposure files
for (conf in unique(na.omit(steiger_dat$confounder))) {

  snps_to_remove <- steiger_dat %>%
    filter(confounder == conf, steiger_remove == TRUE) %>%
    distinct(Phenotype, SNP)

  filtered <- exp_raw %>%
    anti_join(snps_to_remove, by = c("Phenotype", "SNP"))

  conf_stub <- gsub("[^A-Za-z0-9]+", "_", conf)

  fwrite(
    filtered,
    file.path(outdir, paste0("exposures_pa_steiger_filtered_", conf_stub, ".txt")),
    sep = "\t",
    quote = FALSE
  )
}

# Any-confounder filtered file
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

print(summary_dat)
