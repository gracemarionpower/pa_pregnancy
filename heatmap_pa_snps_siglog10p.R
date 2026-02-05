#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
})

# ==============================
# Set project root (robust)
# ==============================
args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg) else ""
script_dir <- if (nzchar(script_path)) dirname(normalizePath(script_path)) else getwd()
project_root <- normalizePath(file.path(script_dir, "..", ".."))
setwd(project_root)

cat("Project root:", project_root, "\n")

# ==============================
# Paths
# ==============================
MERGED_DIR <- "results/snp_only/merged_by_trait"
SNP_FILE   <- "inputs/rsids_pa.txt"
OUT_PDF    <- "results/snp_only/plots/heatmap_pa_snps_signedlog10p.pdf"
dir.create(dirname(OUT_PDF), recursive = TRUE, showWarnings = FALSE)

# ==============================
# Load SNP list
# ==============================
snps <- fread(SNP_FILE, header = FALSE)$V1

# ==============================
# Load merged PLINK outputs
# ==============================
files <- list.files(MERGED_DIR, pattern = "_ALLCHR\\.tsv$", full.names = TRUE)
if (length(files) == 0) stop("No merged *_ALLCHR.tsv files found in ", MERGED_DIR)

dat <- rbindlist(lapply(files, function(f) {
  dt <- fread(f)
  dt[, trait := sub("_ALLCHR\\.tsv$", "", basename(f))]
  dt
}), fill = TRUE)

# ==============================
# Filter + compute signed -log10(P)
# ==============================
dat <- dat[TEST == "ADD"]
dat[, ID := as.character(ID)]
dat[, BETA := suppressWarnings(as.numeric(BETA))]
dat[, P    := suppressWarnings(as.numeric(P))]
dat <- dat[ID %in% snps]

if (nrow(dat) == 0) stop("No rows left after filtering to SNP list. Check ID matching.")

# guard P=0/NA/Inf
dat[P > 0, Pcap := P]
dat[!is.finite(Pcap) | Pcap <= 0, Pcap := .Machine$double.xmin]

dat[, value := sign(BETA) * (-log10(Pcap))]
dat[!is.finite(value), value := NA_real_]

# ==============================
# Trait × SNP matrix
# ==============================
wide <- dcast(dat, trait ~ ID, value.var = "value")
mat <- as.matrix(wide[, -1, with = FALSE])
rownames(mat) <- wide$trait

# Order SNP columns to match input list
common_snps <- intersect(snps, colnames(mat))
mat <- mat[, common_snps, drop = FALSE]

# Drop completely empty rows/cols (keep partial missingness as NA)
mat <- mat[rowSums(!is.na(mat)) > 0, colSums(!is.na(mat)) > 0, drop = FALSE]
if (nrow(mat) == 0 || ncol(mat) == 0) stop("Matrix empty after filtering.")

# ==============================
# Colour scale (robust)
# ==============================
vals <- abs(as.numeric(mat))
vals <- vals[is.finite(vals)]
lim <- as.numeric(quantile(vals, 0.99, na.rm = TRUE))
if (!is.finite(lim) || lim <= 0) lim <- 1

cols   <- colorRampPalette(c("blue", "white", "red"))(100)
breaks <- seq(-lim, lim, length.out = 101)

# ==============================
# Plot with legend (base R)
# ==============================
pdf(OUT_PDF, width = 14, height = 10)

# 2-panel layout: heatmap + legend strip
layout(matrix(c(1, 2), ncol = 2), widths = c(4, 1))

# ---- Heatmap panel ----
par(mar = c(6, 12, 4, 1))

# Base image needs a numeric matrix; we also want NA to be grey:
# Strategy: plot NA as 0 but then overlay grey rectangles for NA cells.
mat_plot <- mat
mat_plot[is.na(mat_plot)] <- 0

image(
  x = 1:ncol(mat_plot),
  y = 1:nrow(mat_plot),
  z = t(mat_plot[nrow(mat_plot):1, ]),
  col = cols,
  breaks = breaks,
  axes = FALSE,
  xlab = "SNP",
  ylab = "Trait",
  main = "PA SNP × Pregnancy Trait associations\n(signed −log10(P))"
)

axis(1, at = 1:ncol(mat_plot), labels = colnames(mat_plot), las = 2, cex.axis = 0.7)
axis(2, at = 1:nrow(mat_plot), labels = rev(rownames(mat_plot)), las = 2, cex.axis = 0.7)

# Overlay NA cells as grey rectangles
na_idx <- which(is.na(mat), arr.ind = TRUE)
if (nrow(na_idx) > 0) {
  for (k in seq_len(nrow(na_idx))) {
    i <- na_idx[k, 1]  # row in mat (trait)
    j <- na_idx[k, 2]  # col in mat (snp)
    # convert to plotted coords (rows are reversed in image call)
    y_plot <- nrow(mat) - i + 1
    rect(j - 0.5, y_plot - 0.5, j + 0.5, y_plot + 0.5, col = "grey90", border = NA)
  }
}

# ---- Legend panel ----
par(mar = c(6, 2, 4, 4))

# Make y a sequence of 101 boundaries so z is 1x100 (valid for image)
y <- seq(-lim, lim, length.out = 101)
z <- matrix(seq(-lim, lim, length.out = 100), nrow = 1)

image(
  x = c(0, 1),
  y = y,
  z = z,
  col = cols,
  breaks = breaks,
  axes = FALSE,
  xlab = "",
  ylab = ""
)
axis(4)
mtext("signed −log10(P)", side = 4, line = 2)

dev.off()

cat("Heatmap written to:", OUT_PDF, "\n")
cat("Traits:", nrow(mat), " SNPs:", ncol(mat), "\n")
cat("Colour limit (99th %ile |value|):", lim, "\n")
