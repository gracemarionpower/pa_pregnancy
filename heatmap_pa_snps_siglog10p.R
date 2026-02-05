#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(pheatmap)
})

# -------- Set project root (robust) --------
args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg) else ""
script_dir <- if (nzchar(script_path)) dirname(normalizePath(script_path)) else getwd()
project_root <- normalizePath(file.path(script_dir, "..", ".."))
setwd(project_root)
cat("Project root:", project_root, "\n")

MERGED_DIR <- "results/snp_only/merged_by_trait"
SNP_FILE   <- "inputs/rsids_pa.txt"
OUT_PDF    <- "results/snp_only/plots/heatmap_Z_beta_over_se.pdf"
dir.create(dirname(OUT_PDF), recursive = TRUE, showWarnings = FALSE)

# Load SNP list
snps <- fread(SNP_FILE, header = FALSE)$V1

# Read all merged trait files
files <- list.files(MERGED_DIR, pattern = "_ALLCHR\\.tsv$", full.names = TRUE)
if (length(files) == 0) stop("No merged *_ALLCHR.tsv files found in: ", MERGED_DIR)

dat_list <- lapply(files, function(f) {
  dt <- fread(f)
  trait <- sub("_ALLCHR\\.tsv$", "", basename(f))
  dt[, trait := trait]
  dt
})
dat <- rbindlist(dat_list, fill = TRUE)

# Filter and compute Z
dat <- dat[TEST == "ADD"]
dat[, ID := as.character(ID)]
dat[, BETA := suppressWarnings(as.numeric(BETA))]
dat[, SE   := suppressWarnings(as.numeric(SE))]
dat <- dat[ID %in% snps]

if (nrow(dat) == 0) stop("No rows left after filtering to SNP list. Check ID matching.")

dat[, Z := BETA / SE]
dat[!is.finite(Z), Z := NA_real_]

# Trait × SNP matrix of Z
wide <- dcast(dat, trait ~ ID, value.var = "Z")
mat <- as.matrix(wide[, -1, with = FALSE])
rownames(mat) <- wide$trait

# Order SNP columns as in SNP list
common_snps <- intersect(snps, colnames(mat))
mat <- mat[, common_snps, drop = FALSE]

# Drop rows/cols that are completely missing (keep partial missingness)
mat <- mat[rowSums(!is.na(mat)) > 0, colSums(!is.na(mat)) > 0, drop = FALSE]
if (nrow(mat) == 0 || ncol(mat) == 0) stop("Matrix empty after filtering.")

# Robust symmetric colour limits using 99th percentile of |Z|
vals <- as.numeric(mat)
vals <- vals[is.finite(vals)]
lim <- as.numeric(quantile(abs(vals), 0.99, na.rm = TRUE))
if (!is.finite(lim) || lim <= 0) lim <- 1

breaks <- seq(-lim, lim, length.out = 101)
cols <- colorRampPalette(c("blue", "white", "red"))(100)

# Plot (legend is included by pheatmap)
pdf(OUT_PDF, width = 14, height = 10)
pheatmap(
  mat,
  color = cols,
  breaks = breaks,
  na_col = "grey90",
  border_color = NA,
  main = "SNP × Trait associations (Z = BETA/SE)",
  fontsize = 9,
  angle_col = 45,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete"
)
dev.off()

cat("Heatmap written to:", OUT_PDF, "\n")
cat("Traits:", nrow(mat), " SNPs:", ncol(mat), "\n")
cat("Z colour limit (approx 99th %ile |Z|):", lim, "\n")

