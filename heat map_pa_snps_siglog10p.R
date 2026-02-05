#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
})

# ============================================================
# heatmap_pa_snps_signedlog10p.R
# SNP × Trait heatmap using signed -log10(P)
# No clustering; explicit legend; BluePebble-safe
# ============================================================

# -------- Set project root robustly --------
args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg) else ""
script_dir <- if (nzchar(script_path)) dirname(normalizePath(script_path)) else getwd()

# Script lives in: <project_root>/scripts/pa_pregnancy
project_root <- normalizePath(file.path(script_dir, "..", ".."))
setwd(project_root)

cat("Project root:", project_root, "\n")

# -------- Paths --------
MERGED_DIR <- "results/snp_only/merged_by_trait"
SNP_FILE   <- "inputs/rsids_pa.txt"
OUT_PDF    <- "results/snp_only/plots/heatmap_signedlog10p.pdf"

dir.create(dirname(OUT_PDF), recursive = TRUE, showWarnings = FALSE)

# -------- Load SNP list --------
snps <- fread(SNP_FILE, header = FALSE)$V1

# -------- Load merged trait files --------
files <- list.files(MERGED_DIR, pattern = "_ALLCHR\\.tsv$", full.names = TRUE)
if (length(files) == 0) stop("No merged *_ALLCHR.tsv files found")

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
if (nrow(dat) == 0) stop("No rows left after filtering to SNP list")

# -------- Evidence metric --------
dat[, signed_log10P := sign(BETA) * (-log10(P))]

# -------- Trait × SNP matrix --------
wide <- dcast(dat, trait ~ ID, value.var = "signed_log10P")
mat <- as.matrix(wide[, -1, with = FALSE])
rownames(mat) <- wide$trait

# Preserve SNP order
common_snps <- intersect(snps, colnames(mat))
mat <- mat[, common_snps, drop = FALSE]

# -------- Clean matrix --------
# Drop rows/cols completely missing
mat <- mat[rowSums(!is.na(mat)) > 0, colSums(!is.na(mat)) > 0, drop = FALSE]

# Replace NA/NaN/Inf with 0 (means "no evidence / not estimated")
mat[!is.finite(mat)] <- NA
mat[is.na(mat)] <- 0

# -------- Cap extremes for visibility --------
cap <- quantile(abs(mat), 0.98, na.rm = TRUE)
if (!is.finite(cap) || cap == 0) cap <- max(abs(mat), na.rm = TRUE)
if (!is.finite(cap) || cap == 0) cap <- 1

mat[mat >  cap] <-  cap
mat[mat < -cap] <- -cap

# -------- Colour scale --------
cols <- colorRampPalette(c("blue", "white", "red"))(101)
brks <- seq(-cap, cap, length.out = length(cols) + 1)

# -------- Plot (image + explicit legend) --------
pdf(OUT_PDF, width = 14, height = 10)

layout(matrix(c(1,2), nrow = 1), widths = c(5,1))

# Heatmap panel
par(mar = c(10, 10, 4, 2))
image(
  x = 1:ncol(mat),
  y = 1:nrow(mat),
  z = t(mat[nrow(mat):1, , drop = FALSE]),  # flip so first trait at top
  col = cols,
  breaks = brks,
  xaxt = "n",
  yaxt = "n",
  xlab = "",
  ylab = "",
  main = "SNP × Trait signed -log10(P)"
)

# X-axis (SNPs)
x_labs <- colnames(mat)
x_step <- if (length(x_labs) <= 60) 1 else if (length(x_labs) <= 120) 2 else 3
axis(
  1,
  at = seq(1, length(x_labs), by = x_step),
  labels = x_labs[seq(1, length(x_labs), by = x_step)],
  las = 2,
  cex.axis = 0.6
)

# Y-axis (Traits)
axis(
  2,
  at = 1:nrow(mat),
  labels = rev(rownames(mat)),
  las = 2,
  cex.axis = 0.7
)

# Legend panel
par(mar = c(10, 2, 4, 6))
image(
  x = 1,
  y = seq(-cap, cap, length.out = 200),
  z = matrix(seq(-cap, cap, length.out = 200), ncol = 1),
  col = cols,
  breaks = brks,
  xaxt = "n",
  yaxt = "n",
  xlab = "",
  ylab = ""
)
axis(4, at = pretty(c(-cap, cap), 5), las = 2, cex.axis = 0.8)
mtext("signed -log10(P)", side = 4, line = 3)

dev.off()

cat("Heatmap written to:", OUT_PDF, "\n")
cat("Traits:", nrow(mat), " SNPs:", ncol(mat), "\n")
