#!/bin/bash
#SBATCH --job-name=PA_condF_MET
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH -p compute
#SBATCH --account=SSCM013902
#SBATCH --mem=32GB
#SBATCH --time=1-00:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=grace.power@bristol.ac.uk
#SBATCH --output=/user/work/sd20930/project_pa_mrpreg/variant_relevance/logs/conditional_f_%j.out
#SBATCH --error=/user/work/sd20930/project_pa_mrpreg/variant_relevance/logs/conditional_f_%j.err

set -euo pipefail

PROJECT_ROOT="/user/work/sd20930/project_pa_mrpreg/variant_relevance"
cd "$PROJECT_ROOT"

PLINK2="/user/work/sd20930/project_pa_mrpreg/plink2"

BGEN_DIR="/group/alspac/gi_1000g_g0m_g1/released/2015-10-30/data/dosage_bgen"
SAMPLE="/group/alspac/gi_1000g_g0m_g1/released/2015-10-30/data/data.sample"

PHENO_DIR="${PROJECT_ROOT}/inputs/traits_dedup"
COVAR="${PROJECT_ROOT}/inputs/covar_age_pcs.txt"

MET_TRAIT="b_mat_exercise_met_sum"
PHENO="${PHENO_DIR}/${MET_TRAIT}.pheno.txt"

OUTDIR="${PROJECT_ROOT}/results/conditional_f_MET"
DOSAGEDIR="${OUTDIR}/dosages"
SNPDIR="${OUTDIR}/snp_lists"
RSCRIPT="${OUTDIR}/conditional_f_MET.R"

mkdir -p "$OUTDIR" "$DOSAGEDIR" "$SNPDIR" "${PROJECT_ROOT}/logs"

[[ -x "$PLINK2" ]] || { echo "ERROR: plink2 not executable: $PLINK2"; exit 1; }
[[ -f "$SAMPLE" ]] || { echo "ERROR: sample file missing: $SAMPLE"; exit 1; }
[[ -f "$PHENO" ]] || { echo "ERROR: phenotype file missing: $PHENO"; exit 1; }
[[ -f "$COVAR" ]] || { echo "ERROR: covariate file missing: $COVAR"; exit 1; }

echo "Phenotype: $PHENO"
echo "Covariates: $COVAR"
echo "Output directory: $OUTDIR"

cat > "${SNPDIR}/MVPA_leisure.txt" <<EOF
rs1160545
rs7613360
rs1691471
rs13201721
rs1625595
rs385301
EOF

cat > "${SNPDIR}/Vigorous_PA.txt" <<EOF
rs1248860
rs2764261
rs13243553
rs3781411
rs328902
EOF

cat > "${SNPDIR}/Total_log_acceleration.txt" <<EOF
rs2532402
rs1268539
EOF

cat > "${SNPDIR}/Sedentary_time.txt" <<EOF
rs26579
rs25981
rs1858242
rs34858520
EOF

for INSTRUMENT in MVPA_leisure Vigorous_PA Total_log_acceleration Sedentary_time; do
  RSIDS="${SNPDIR}/${INSTRUMENT}.txt"

  echo "Exporting dosages for ${INSTRUMENT}"

  for CHR in $(seq -w 1 22); do
    BGEN="${BGEN_DIR}/data_chr${CHR}.bgen"
    OUT="${DOSAGEDIR}/${INSTRUMENT}_chr${CHR}_dosage"

    if "$PLINK2" \
      --threads 2 \
      --bgen "$BGEN" ref-first \
      --sample "$SAMPLE" \
      --extract "$RSIDS" \
      --export A \
      --out "$OUT" >/dev/null 2>&1; then
      echo "  chr${CHR}: exported"
    else
      echo "  chr${CHR}: no matching SNPs"
      rm -f "${OUT}.raw"
    fi
  done
done

cat > "$RSCRIPT" <<'EOF'
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)

dosage_dir <- args[1]
pheno_file <- args[2]
covar_file <- args[3]
outfile <- args[4]

pheno <- fread(pheno_file, data.table = FALSE, header = FALSE)
covar <- fread(covar_file, data.table = FALSE, header = FALSE)

names(pheno) <- c("FID", "IID", "trait")
names(covar) <- c("FID", "IID", "age", paste0("PC", 1:10))

base_dat <- pheno %>%
  mutate(trait = as.numeric(trait)) %>%
  filter(!is.na(trait), trait != -9) %>%
  left_join(covar, by = c("FID", "IID"))

covar_cols <- c("age", paste0("PC", 1:10))

instruments <- c(
  "MVPA_leisure",
  "Vigorous_PA",
  "Total_log_acceleration",
  "Sedentary_time"
)

results <- list()

for (instrument in instruments) {

  message("Processing ", instrument)

  files <- list.files(
    dosage_dir,
    pattern = paste0("^", instrument, "_chr[0-9]+_dosage\\.raw$"),
    full.names = TRUE
  )

  if (length(files) == 0) {
    warning("No dosage files found for ", instrument)
    next
  }

  dosage_list <- lapply(files, function(f) {
    x <- fread(f, data.table = FALSE)
    names(x)[1:2] <- c("FID", "IID")
    x
  })

  dosage <- Reduce(function(x, y) full_join(x, y, by = c("FID", "IID")), dosage_list)

  non_snp_cols <- c("FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE")
  snp_cols <- setdiff(names(dosage), non_snp_cols)

  dat <- base_dat %>%
    left_join(dosage %>% select(FID, IID, all_of(snp_cols)), by = c("FID", "IID"))

  snp_cols <- snp_cols[
    sapply(dat[snp_cols], function(x) {
      x <- suppressWarnings(as.numeric(x))
      sum(!is.na(x)) > 0 && sd(x, na.rm = TRUE) > 0
    })
  ]

  if (length(snp_cols) == 0) {
    warning("No usable SNPs for ", instrument)
    next
  }

  model_dat <- dat %>%
    select(trait, all_of(covar_cols), all_of(snp_cols))

  names(model_dat) <- make.names(names(model_dat), unique = TRUE)

  trait_name <- "trait"
  covar_names <- make.names(covar_cols, unique = TRUE)
  snp_names <- setdiff(names(model_dat), c(trait_name, covar_names))

  model_dat <- model_dat %>%
    mutate(across(everything(), as.numeric)) %>%
    filter(complete.cases(.))

  if (nrow(model_dat) == 0) {
    warning("No complete cases for ", instrument)
    next
  }

  reduced_formula <- as.formula(
    paste(trait_name, "~", paste(covar_names, collapse = " + "))
  )

  full_formula <- as.formula(
    paste(trait_name, "~", paste(c(covar_names, snp_names), collapse = " + "))
  )

  fit_reduced <- lm(reduced_formula, data = model_dat)
  fit_full <- lm(full_formula, data = model_dat)

  a <- anova(fit_reduced, fit_full)

  r2_reduced <- summary(fit_reduced)$r.squared
  r2_full <- summary(fit_full)$r.squared
  partial_r2 <- (r2_full - r2_reduced) / (1 - r2_reduced)

  results[[instrument]] <- tibble(
    instrument = instrument,
    n = nrow(model_dat),
    n_snps_available = length(snp_names),
    conditional_F = a$F[2],
    p_value = a$`Pr(>F)`[2],
    r2_reduced = r2_reduced,
    r2_full = r2_full,
    partial_r2 = partial_r2,
    snps_included = paste(snp_cols, collapse = ";")
  )
}

final <- bind_rows(results)

fwrite(final, outfile, sep = "\t", quote = FALSE)

print(final)
EOF

Rscript "$RSCRIPT" \
  "$DOSAGEDIR" \
  "$PHENO" \
  "$COVAR" \
  "${OUTDIR}/conditional_F_b_mat_exercise_met_sum_results.tsv"

echo "Finished."
echo "Results written to:"
echo "${OUTDIR}/conditional_F_b_mat_exercise_met_sum_results.tsv"
