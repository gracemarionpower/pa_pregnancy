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
PHENO_DIR  <- "inputs/traits_dedup"
SNP_FILE   <- "inputs/rsids_pa.txt"

OUT_PDF    <- "results/snp_only/plots/heatmap_betaSD_clustered.pdf"
OUT_SD_TSV <- "results/snp_only/plots/trait_sd_table.tsv"

dir.create(dirname(OUT_PDF), recursive = TRUE, showWarnings = FALSE)

# ---- user-facing settings ----
P_HILITE <- 0.05      # outline cells with P < this
FIX_LIM  <- 0.05      # show colour scale as ±0.05 SD (set NA to auto)

snps <- fread(SNP_FILE, header = FALSE)$V1

files <- list.files(MERGED_DIR, pattern = "_ALLCHR\\.tsv$", full.names = TRUE)
if (length(files) == 0) stop("No merged *_ALLCHR.tsv files found in: ", MERGED_DIR)

dat <- rbindlist(lapply(files, function(f) {
  dt <- fread(f)
  dt[, trait := sub("_ALLCHR\\.tsv$", "", basename(f))]
  dt
}), fill = TRUE)

# Filter + numeric conversion
dat <- dat[TEST == "ADD"]
dat[, ID := as.character(ID)]
dat[, `:=`(
  BETA = suppressWarnings(as.numeric(BETA)),
  SE   = suppressWarnings(as.numeric(SE)),
  P    = suppressWarnings(as.numeric(P))
)]
dat <- dat[ID %in% snps]
if (nrow(dat) == 0) stop("No rows left after filtering to SNP list.")

# Compute SD per trait from phenotype files (column 3)
traits <- sort(unique(dat$trait))

sd_dt <- rbindlist(lapply(traits, function(tr) {
  pheno_path <- file.path(PHENO_DIR, paste0(tr, ".pheno.txt"))
  if (!file.exists(pheno_path)) {
    return(data.table(trait = tr, pheno_sd = NA_real_, n = NA_integer_, note="pheno_missing"))
  }
  ph <- fread(pheno_path)
  if (ncol(ph) < 3) {
    return(data.table(trait = tr, pheno_sd = NA_real_, n = NA_integer_, note="pheno_<3cols"))
  }
  y <- suppressWarnings(as.numeric(ph[[3]]))
  y <- y[is.finite(y)]
  n <- length(y)
  s <- if (n >= 2) sd(y) else NA_real_
  note <- if (is.finite(s) && s > 0) "ok" else "sd_bad"
  data.table(trait = tr, pheno_sd = s, n = n, note = note)
}), fill = TRUE)

fwrite(sd_dt, OUT_SD_TSV, sep = "\t")
cat("Trait SD table written to:", OUT_SD_TSV, "\n")

dat <- merge(dat, sd_dt[, .(trait, pheno_sd)], by = "trait", all.x = TRUE)

# SD-scaled effect + SD-scaled SE (useful for CI later)
dat[, BETA_SD := BETA / pheno_sd]
dat[, SE_SD   := SE   / pheno_sd]
dat[!is.finite(BETA_SD), BETA_SD := NA_real_]
dat[!is.finite(SE_SD),   SE_SD   := NA_real_]

# Matrices: effect + p for highlighting
wide_eff <- dcast(dat, trait ~ ID, value.var = "BETA_SD")
wide_p   <- dcast(dat, trait ~ ID, value.var = "P")

mat <- as.matrix(wide_eff[, -1, with = FALSE])
pmat <- as.matrix(wide_p[, -1, with = FALSE])
rownames(mat) <- wide_eff$trait
rownames(pmat) <- wide_p$trait

# Order SNPs
common_snps <- intersect(snps, colnames(mat))
mat <- mat[, common_snps, drop = FALSE]
pmat <- pmat[, common_snps, drop = FALSE]

# Drop empty rows/cols
keep_r <- rowSums(!is.na(mat)) > 0
keep_c <- colSums(!is.na(mat)) > 0
mat <- mat[keep_r, keep_c, drop = FALSE]
pmat <- pmat[keep_r, keep_c, drop = FALSE]
if (nrow(mat) == 0 || ncol(mat) == 0) stop("Matrix empty after filtering.")

# Clustering using NA->0 only for ordering
mat_cl <- mat
mat_cl[is.na(mat_cl)] <- 0
row_hc <- hclust(dist(mat_cl), method = "complete")
col_hc <- hclust(dist(t(mat_cl)), method = "complete")

mat_ord <- mat[row_hc$order, col_hc$order, drop = FALSE]
pmat_ord <- pmat[row_hc$order, col_hc$order, drop = FALSE]

# Set colour limits
vals <- abs(as.numeric(mat_ord))
vals <- vals[is.finite(vals)]
auto_lim <- as.numeric(quantile(vals, 0.99, na.rm = TRUE))
if (!is.finite(auto_lim) || auto_lim <= 0) auto_lim <- 0.05

lim <- if (is.finite(FIX_LIM) && FIX_LIM > 0) FIX_LIM else auto_lim

# If fixed limit is too small (would saturate everything), bump it
if (lim < auto_lim * 0.5) lim <- auto_lim

cols   <- colorRampPalette(c("blue", "white", "red"))(100)
breaks <- seq(-lim, lim, length.out = 101)

# Plot
pdf(OUT_PDF, width = 18, height = 10)
layout(matrix(c(1, 2), ncol = 2), widths = c(6, 1))

# Heatmap panel
par(mar = c(7, 14, 4, 1))
mat_plot <- mat_ord
mat_plot[is.na(mat_plot)] <- 0

image(
  x = 1:ncol(mat_plot),
  y = 1:nrow(mat_plot),
  z = t(mat_plot[nrow(mat_plot):1, ]),
  col = cols,
  breaks = breaks,
  axes = FALSE,
  xlab = "SNP (clustered)",
  ylab = "Trait (clustered)",
  main = sprintf("SNP × Trait effects (SD units = BETA/SD)\nColour scale fixed at ±%.3f SD; outlined if P<%.2g", lim, P_HILITE)
)

axis(1, at = 1:ncol(mat_plot), labels = colnames(mat_plot), las = 2, cex.axis = 0.75)
axis(2, at = 1:nrow(mat_plot), labels = rev(rownames(mat_plot)), las = 2, cex.axis = 0.75)

# Grey out true NA cells
na_idx <- which(is.na(mat_ord), arr.ind = TRUE)
if (nrow(na_idx) > 0) {
  for (k in seq_len(nrow(na_idx))) {
    i <- na_idx[k, 1]; j <- na_idx[k, 2]
    y_plot <- nrow(mat_ord) - i + 1
    rect(j - 0.5, y_plot - 0.5, j + 0.5, y_plot + 0.5, col = "grey90", border = NA)
  }
}

# Outline significant cells
sig_idx <- which(is.finite(pmat_ord) & pmat_ord < P_HILITE & is.finite(mat_ord), arr.ind = TRUE)
if (nrow(sig_idx) > 0) {
  for (k in seq_len(nrow(sig_idx))) {
    i <- sig_idx[k, 1]; j <- sig_idx[k, 2]
    y_plot <- nrow(mat_ord) - i + 1
    rect(j - 0.5, y_plot - 0.5, j + 0.5, y_plot + 0.5, border = "black", lwd = 1.2)
  }
}

# Legend panel (big)
par(mar = c(7, 3, 4, 6))
y <- seq(-lim, lim, length.out = 101)
z <- matrix(seq(-lim, lim, length.out = 100), nrow = 1)
image(x = c(0, 1), y = y, z = z, col = cols, breaks = breaks, axes = FALSE)
axis(4, at = pretty(c(-lim, lim), 5))
mtext("Effect (SD units)", side = 4, line = 3, cex = 1.1)
mtext("blue = negative, red = positive", side = 4, line = 2, cex = 0.9)
mtext("grey = missing", side = 4, line = 1, cex = 0.9)

dev.off()

cat("Heatmap written to:", OUT_PDF, "\n")
cat("Traits:", nrow(mat_ord), " SNPs:", ncol(mat_ord), "\n")
cat("Auto lim (99th %ile |BETA_SD|):", auto_lim, "  Used lim:", lim, "\n")
