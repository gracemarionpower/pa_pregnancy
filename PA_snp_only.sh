#!/bin/bash
#SBATCH --job-name=PA_snp_only
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH -p compute
#SBATCH --account=SSCM013902
#SBATCH --mem=32GB
#SBATCH --time=2-00:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=grace.power@bristol.ac.uk
#SBATCH --output=/user/work/sd20930/project_pa_mrpreg/variant_relevance/logs/snp_only/plink2_%j.out
#SBATCH --error=/user/work/sd20930/project_pa_mrpreg/variant_relevance/logs/snp_only/plink2_%j.err

set -euo pipefail

# Submit with: sbatch PA_snp_only.sh

PROJECT_ROOT="/user/work/sd20930/project_pa_mrpreg/variant_relevance"
cd "$PROJECT_ROOT"

# ---- Paths ----
PLINK2="/user/work/sd20930/project_pa_mrpreg/plink2"

BGEN_DIR="/group/alspac/gi_1000g_g0m_g1/released/2015-10-30/data/dosage_bgen"
SAMPLE="/group/alspac/gi_1000g_g0m_g1/released/2015-10-30/data/data.sample"

TRAIT_LIST="${PROJECT_ROOT}/inputs/list_traits_all.txt"
PHENO_DIR="${PROJECT_ROOT}/inputs/traits_dedup"
COVAR="${PROJECT_ROOT}/inputs/covar_age_pcs.txt"
RSIDS="${PROJECT_ROOT}/inputs/rsids_pa.txt"

OUTDIR="${PROJECT_ROOT}/results/snp_only"
LOGDIR="${PROJECT_ROOT}/logs/snp_only"
PRESDIR="${OUTDIR}/presence"

mkdir -p "$OUTDIR" "$LOGDIR" "$PRESDIR"

# ---- Fast fail ----
[[ -x "$PLINK2" ]] || { echo "ERROR: plink2 not executable: $PLINK2"; exit 1; }
[[ -f "$SAMPLE" ]] || { echo "ERROR: sample file missing: $SAMPLE"; exit 1; }
[[ -f "$TRAIT_LIST" ]] || { echo "ERROR: trait list missing: $TRAIT_LIST"; exit 1; }
[[ -f "$COVAR" ]] || { echo "ERROR: covar file missing: $COVAR"; exit 1; }
[[ -f "$RSIDS" ]] || { echo "ERROR: RSID list missing: $RSIDS"; exit 1; }

echo "Running in: $(pwd)"
echo "Using plink2: $PLINK2"
echo "RSIDs file: $RSIDS"
echo "Trait list: $TRAIT_LIST"
echo "Pheno dir:  $PHENO_DIR"
echo "Covar:      $COVAR"
echo "Outdir:     $OUTDIR"

# ---- 1) Presence check per chromosome (robust to 0 hits) ----
echo "Building per-chromosome SNP presence lists..."
for CHR in $(seq -w 1 22); do
  BGEN="${BGEN_DIR}/data_chr${CHR}.bgen"
  [[ -f "$BGEN" ]] || { echo "ERROR: bgen missing: $BGEN"; exit 1; }

  out="${PRESDIR}/chr${CHR}"

  # If extract finds 0 variants, plink2 exits nonzero; we catch that and write an empty snplist.
  if "$PLINK2" \
      --bgen "$BGEN" ref-first \
      --sample "$SAMPLE" \
      --extract "$RSIDS" \
      --write-snplist \
      --out "$out" >/dev/null 2>&1; then
    :
  else
    : > "${out}.snplist"
  fi

  n=$(wc -l < "${out}.snplist" 2>/dev/null || echo 0)
  echo "chr${CHR}: ${n} SNPs found"
done

# ---- 2) Run SNP-only tests: traits × chr01–22 (skip empty chr) ----
while read -r TRAIT; do
  [[ -z "$TRAIT" ]] && continue

  PHENO="${PHENO_DIR}/${TRAIT}.pheno.txt"
  if [[ ! -f "$PHENO" ]]; then
    echo "WARNING: phenotype missing, skipping: $PHENO"
    continue
  fi

  for CHR in $(seq -w 1 22); do
    snplist="${PRESDIR}/chr${CHR}.snplist"
    if [[ ! -s "$snplist" ]]; then
      echo "Trait=${TRAIT} chr=${CHR}: no SNPs on this chr, skipping"
      continue
    fi

    BGEN="${BGEN_DIR}/data_chr${CHR}.bgen"
    nsnps=$(wc -l < "$snplist")

    echo "Trait=${TRAIT} chr=${CHR}: running on ${nsnps} SNPs"

    "$PLINK2" \
      --threads 2 \
      --bgen "$BGEN" ref-first \
      --sample "$SAMPLE" \
      --extract "$snplist" \
      --pheno "$PHENO" --pheno-col-nums 3 \
      --covar "$COVAR" --covar-col-nums 3-13 \
      --covar-variance-standardize \
      --glm hide-covar firth-fallback \
      --out "${OUTDIR}/${TRAIT}_chr${CHR}_snpsOnly"
  done
done < "$TRAIT_LIST"


