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

# ==============================
# Inputs / outputs
# ==============================
MERGED_DIR <- "results/snp_only/merged_by_trait"
PHENO_DIR  <- "inputs/traits_dedup"
SNP_FILE   <- "inputs/rsids_pa.txt"
OUT_PDF    <- "results/snp_only/plots/heatmap_betaSD_clustered.pdf"
dir.create(dirname(OUT_PDF), recursive = TRUE, showWarnings = FALSE)

# ==============================
# Load SNP list
# ==============================
snps <- fread(SNP_FILE, header = FALSE)$V1

# ==============================
# Load merged results
# ==============================
files <- list.files(MERGED_DIR, pattern = "_ALLCHR\\.tsv$", full.names = TRUE)
if (length(files) == 0) stop("No merged *_ALLCHR.tsv files found in: ", MERGED_DIR)

dat <- rbindlist(lapply(files, function(f) {
  dt <- fread(f)
  trait <- sub("_ALLCHR\\.tsv$", "", basename(f))
  dt[, trait := trait]
  dt
}), fill = TRUE)

# Filter to additive model + SNP list
dat <- dat[TEST == "ADD"]
dat[, ID := as.character(ID)]
dat[, BETA := suppressWarnings(as.numeric(BETA))]
dat <- dat[ID %in% snps]

if (nrow(dat) == 0) stop("No rows left after filtering to SNP list. Check ID matching.")

# ==============================
# Compute SD for each trait from phenotype files (column 3)
# ==============================
traits_in_results <- sort(unique(dat$trait))

sd_dt <- rbindlist(lapply(traits_in_results, function(tr) {
  pheno_path <- file.path(PHENO_DIR, paste0(tr, ".pheno.txt"))
  if (!file.exists(pheno_path)) {
    return(data.table(trait = tr, pheno_sd = NA_real_, n_pheno = NA_integer_, note = "pheno_file_missing"))
  }
  ph <- fread(pheno_path)
  if (ncol(ph) < 3) {
    return(data.table(trait = tr, pheno_sd = NA_real_, n_pheno = NA_integer_, note = "pheno_has_<3_cols"))
  }
  y <- suppressWarnings(as.numeric(ph[[3]]))
  y <- y[is.finite(y)]
  n <- length(y)
  s <- if (n >= 2) sd(y) else NA_real_
  note <- if (is.finite(s) && s > 0) "ok" else "sd_not_finite_or_zero"
  data.table(trait = tr, pheno_sd = s, n_pheno = n, note = note)
}), fill = TRUE)

# Merge SD into association results
dat <- merge(dat, sd_dt[, .(trait, pheno_sd)], by = "trait", all.x = TRUE)

# Convert to SD units
dat[, BETA_SD := BETA / pheno_sd]
dat[!is.finite(BETA_SD), BETA_SD := NA_real_]

# If SD missing for some traits, they will drop later if all NA
n_missing_sd <- sum(!is.finite(sd_dt$pheno_sd) | sd_dt$pheno_sd <= 0, na.rm = TRUE)
if (n_missing_sd > 0) {
  cat("Warning:", n_missing_sd, "traits had missing/zero SD and will likely drop.\n")
}

# ==============================
# Build Trait × SNP matrix (BETA in SD units)
# ==============================
wide <- dcast(dat, trait ~ ID, value.var = "BETA_SD")
mat <- as.matrix(wide[, -1, with = FALSE])
rownames(mat) <- wide$trait

# Keep SNP column order as SNP list (but only those present)
common_snps <- intersect(snps, colnames(mat))
mat <- mat[, common_snps, drop = FALSE]

# Drop empty rows/cols (all NA)
mat <- mat[rowSums(!is.na(mat)) > 0, colSums(!is.na(mat)) > 0, drop = FALSE]
if (nrow(mat) == 0 || ncol(mat) == 0) stop("Matrix empty after filtering. Check SD calculation + SNP matching.")

# ==============================
# Clustering (need numeric matrix with no NA)
# We impute NA as 0 ONLY for clustering / ordering.
# We still show true NA as grey in the plot.
# ==============================
mat_for_clust <- mat
mat_for_clust[!is.finite(mat_for_clust)] <- NA
mat_for_clust[is.na(mat_for_clust)] <- 0

row_hc <- hclust(dist(mat_for_clust), method = "complete")
col_hc <- hclust(dist(t(mat_for_clust)), method = "complete")

mat_ord <- mat[row_hc$order, col_hc$order, drop = FALSE]

# ==============================
# Colour scale & legend (OBVIOUS)
# Use robust symmetric limits (99th percentile of |BETA_SD|)
# ==============================
vals <- abs(as.numeric(mat_ord))
vals <- vals[is.finite(vals)]
lim <- as.numeric(quantile(vals, 0.99, na.rm = TRUE))
if (!is.finite(lim) || lim <= 0) lim <- max(abs(vals), na.rm = TRUE)
if (!is.finite(lim) || lim <= 0) lim <- 1

# Make the legend ticks easy to read
tick_vals <- pretty(c(-lim, lim), n = 5)

cols   <- colorRampPalette(c("blue", "white", "red"))(100)
breaks <- seq(-lim, lim, length.out = 101)

# ==============================
# Plot (base R) with big legend
# ==============================
pdf(OUT_PDF, width = 16, height = 10)
layout(matrix(c(1, 2), ncol = 2), widths = c(5, 1))

# Heatmap panel
par(mar = c(7, 14, 4, 1))

# Plot matrix for image: NAs replaced with 0 temporarily
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
  main = "SNP × Trait effects in SD units (BETA / SD(trait))"
)

axis(1, at = 1:ncol(mat_plot), labels = colnames(mat_plot), las = 2, cex.axis = 0.75)
axis(2, at = 1:nrow(mat_plot), labels = rev(rownames(mat_plot)), las = 2, cex.axis = 0.75)

# Overlay true NA cells as grey blocks (so missingness is obvious)
na_idx <- which(is.na(mat_ord), arr.ind = TRUE)
if (nrow(na_idx) > 0) {
  for (k in seq_len(nrow(na_idx))) {
    i <- na_idx[k, 1]
    j <- na_idx[k, 2]
    y_plot <- nrow(mat_ord) - i + 1
    rect(j - 0.5, y_plot - 0.5, j + 0.5, y_plot + 0.5, col = "grey90", border = NA)
  }
}

# Legend panel (big + labelled clearly)
par(mar = c(7, 3, 4, 6))
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
axis(4, at = tick_vals, labels = format(tick_vals, digits = 2))
mtext("Effect (SD units)", side = 4, line = 3, cex = 1.1)
mtext("blue = negative, red = positive", side = 4, line = 2, cex = 0.9)
mtext("grey = missing SNP–trait result", side = 4, line = 1, cex = 0.9)

dev.off()

cat("Heatmap written to:", OUT_PDF, "\n")
cat("Traits:", nrow(mat_ord), " SNPs:", ncol(mat_ord), "\n")
cat("Colour limit (approx 99th %ile |BETA_SD|):", lim, "\n")
