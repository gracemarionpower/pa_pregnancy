#!/bin/bash
#SBATCH --job-name=PA_PLINK2_alltraits
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH -p compute
#SBATCH --account=SSCM013902
#SBATCH --mem=32GB
#SBATCH --time=10-00:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=grace.power@bristol.ac.uk
#SBATCH --output=logs/alltraits_allchr/plink2_%j.out
#SBATCH --error=logs/alltraits_allchr/plink2_%j.err

set -euo pipefail

cd /user/work/sd20930/project_pa_mrpreg/variant_relevance

mkdir -p logs/alltraits_allchr results/per_chr

PLINK2="/user/work/sd20930/project_pa_mrpreg/plink2"
BGEN_DIR="/group/alspac/gi_1000g_g0m_g1/released/2015-10-30/data/dosage_bgen"
SAMPLE="/group/alspac/gi_1000g_g0m_g1/released/2015-10-30/data/data.sample"

TRAIT_LIST="inputs/list_traits_all.txt"
PHENO_DIR="inputs/traits_dedup"
COVAR="inputs/covar_age_pcs.txt"

OUTDIR="results/per_chr"

# chr01..chr22
for CHR in $(seq -w 1 22); do
  BGEN="${BGEN_DIR}/data_chr${CHR}.bgen"

  while read -r TRAIT; do
    [[ -z "${TRAIT}" ]] && continue
    PHENO="${PHENO_DIR}/${TRAIT}.pheno.txt"

    echo "Running trait=${TRAIT} chr=${CHR}"

    "${PLINK2}" \
      --bgen "${BGEN}" ref-first \
      --sample "${SAMPLE}" \
      --pheno "${PHENO}" \
      --pheno-col-nums 3 \
      --covar "${COVAR}" \
      --covar-col-nums 3-13 \
      --covar-variance-standardize \
      --glm linear hide-covar \
      --out "${OUTDIR}/gwas_${TRAIT}_chr${CHR}"

  done < "${TRAIT_LIST}"
done
