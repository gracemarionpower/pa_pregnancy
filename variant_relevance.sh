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
#SBATCH --output=/user/work/sd20930/project_pa_mrpreg/variant_relevance/logs/alltraits_allchr/plink2_%j.out
#SBATCH --error=/user/work/sd20930/project_pa_mrpreg/variant_relevance/logs/alltraits_allchr/plink2_%j.err

set -euo pipefail

# ---- Set project root explicitly (most robust) ----
PROJECT_ROOT="/user/work/sd20930/project_pa_mrpreg/variant_relevance"
cd "$PROJECT_ROOT"

mkdir -p results/per_chr logs/alltraits_allchr

# ---- Paths ----
PLINK2="/user/work/sd20930/project_pa_mrpreg/plink2"

BGEN_DIR="/group/alspac/gi_1000g_g0m_g1/released/2015-10-30/data/dosage_bgen"
SAMPLE="/group/alspac/gi_1000g_g0m_g1/released/2015-10-30/data/data.sample"

TRAIT_LIST="inputs/list_traits_all.txt"
PHENO_DIR="inputs/traits_dedup"
COVAR="inputs/covar_age_pcs.txt"
OUTDIR="results/per_chr"

# ---- Fast fail if key files missing (helps debug in logs) ----
[[ -x "$PLINK2" ]] || { echo "ERROR: plink2 not executable: $PLINK2"; exit 1; }
[[ -f "$SAMPLE" ]] || { echo "ERROR: sample file missing: $SAMPLE"; exit 1; }
[[ -f "$TRAIT_LIST" ]] || { echo "ERROR: trait list missing: $TRAIT_LIST"; exit 1; }
[[ -f "$COVAR" ]] || { echo "ERROR: covar file missing: $COVAR"; exit 1; }

echo "Running in: $(pwd)"
echo "Using plink2: $PLINK2"
echo "Trait list: $TRAIT_LIST"
echo "Phenos: $PHENO_DIR"
echo "Covar: $COVAR"

# ---- Run GWAS: all traits × chr01–22 ----
for CHR in $(seq -w 1 22); do
  BGEN="${BGEN_DIR}/data_chr${CHR}.bgen"
  [[ -f "$BGEN" ]] || { echo "ERROR: bgen missing: $BGEN"; exit 1; }

  while read -r TRAIT; do
    [[ -z "$TRAIT" ]] && continue

    PHENO="${PHENO_DIR}/${TRAIT}.pheno.txt"
    if [[ ! -f "$PHENO" ]]; then
      echo "WARNING: phenotype missing, skipping: $PHENO"
      continue
    fi

    echo "Running trait=${TRAIT} chr=${CHR}"

    "$PLINK2" \
      --bgen "$BGEN" ref-first \
      --sample "$SAMPLE" \
      --pheno "$PHENO" \
      --pheno-col-nums 3 \
      --covar "$COVAR" \
      --covar-col-nums 3-13 \
      --covar-variance-standardize \
      --glm hide-covar firth-fallback \
      --out "${OUTDIR}/gwas_${TRAIT}_chr${CHR}"

  done < "$TRAIT_LIST"
done
