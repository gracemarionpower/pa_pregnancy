#!/usr/bin/env Rscript
suppressPackageStartupMessages({ library(data.table) })

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
OUT_PDF    <- "results/snp_only/plots/heatmap_signedlog10p.pdf"
dir.create(dirname(OUT_PDF), recursive = TRUE, showWarnings = FALSE)

# -------- Load SNP list --------
snps <- fread(SNP_FILE, header = FALSE)$V1

# -------- Load merged trait files --------
files <- list.files(MERGED_DIR, pattern = "_ALLCHR\\.tsv$", full.names = TRUE)
if (length(files) == 0) stop("No merged *_ALLCHR.tsv files found in: ", MERGED_DIR)

dat_list <- lapply(files, function(f) {
  dt <- fread(f)
  trait <- sub("_ALLCHR\\.tsv$", "", basename(f))
  dt[, trait := trait]
  dt
})
dat <- rbindlist(dat_list, fill = TRUE)

# -------- Filter + numeric conversion --------
dat <- dat[TEST == "ADD"]
dat[, ID := as.character(ID)]
dat[, BETA := suppressWarnings(as.numeric(BETA))]
dat[, P := suppressWarnings(as.numeric(P))]

# Keep SNPs of interest
dat <- dat[ID %in% snps]
if (nrow(dat) == 0) stop("No rows left after filtering to SNP list. Check ID matching.")

# -------- Evidence metric: signed -log10(P) --------
# Direction from BETA sign, evidence from P-value magnitude
dat[, signed_log10P := sign(BETA) * (-log10(P))]

# -------- Create trait × SNP matrix --------
wide <- dcast(dat, trait ~ ID, value.var = "signed_log10P")
mat <- as.matrix(wide[, -1, with = FALSE])
rownames(mat) <- wide$trait

# Order SNP columns as in SNP list
common_snps <- intersect(snps, colnames(mat))
mat <- mat[, common_snps, drop = FALSE]

# -------- Clean matrix for clustering --------
# Drop rows/cols completely missing
keep_rows <- rowSums(!is.na(mat)) > 0
keep_cols <- colSums(!is.na(mat)) > 0
mat <- mat[keep_rows, keep_cols, drop = FALSE]

# Replace NA/NaN/Inf with 0 so dist()/hclust() works
mat[!is.finite(mat)] <- NA
mat[is.na(mat)] <- 0

# Symmetric color limits
lim <- max(abs(mat), na.rm = TRUE)
if (!is.finite(lim) || lim == 0) lim <- 1

cols <- colorRampPalette(c("blue", "white", "red"))(100)

# -------- Plot heatmap + legend --------
pdf(OUT_PDF, width = 12, height = 9)

# Make room on the right for a legend
op <- par(mar = c(5, 8, 4, 8))

heatmap(
  mat,
  scale = "none",
  col = cols,
  margins = c(10, 12),
  main = "SNP × Trait signed -log10(P)",
  zlim = c(-lim, lim)
)

# Add legend (color key) on the right
par(new = TRUE, mar = c(5, 1, 4, 6))
image(
  z = matrix(seq(-lim, lim, length.out = 100), ncol = 1),
  col = cols,
  axes = FALSE
)
axis(
  side = 4,
  at = seq(0, 1, length.out = 5),
  labels = round(seq(-lim, lim, length.out = 5), 1)
)
mtext("signed -log10(P)", side = 4, line = 3)

par(op)
dev.off()

cat("Heatmap written to:", OUT_PDF, "\n")
cat("Traits:", nrow(mat), " SNPs:", ncol(mat), "\n")
