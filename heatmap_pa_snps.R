#!/usr/bin/env Rscript

# heatmap_pa_snps.R
# Create SNP × trait heatmaps from merged PLINK2 outputs:
#   results/snp_only/merged_by_trait/*_ALLCHR.tsv
#
# Run from anywhere, e.g.:
#   Rscript heatmap_pa_snps.R
#   Rscript heatmap_pa_snps.R --metric signedlogp --out results/snp_only/plots/heatmap_signedlogp.pdf
#
# Metrics:
#   beta       = BETA
#   z          = BETA/SE
#   signedlogp = sign(BETA) * -log10(P)

suppressPackageStartupMessages({
  library(data.table)
  library(tools)
  library(pheatmap)
})

# -------- helpers --------
get_script_dir <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  if (length(file_arg) == 0) return(getwd())
  normalizePath(sub("^--file=", "", file_arg[1]))
  dirname(normalizePath(sub("^--file=", "", file_arg[1])))
}

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  out <- list(metric = "beta", out = "", cluster_rows = TRUE, cluster_cols = TRUE, width = 12, height = 9)
  if (length(args) == 0) return(out)

  i <- 1
  while (i <= length(args)) {
    key <- args[i]
    val <- if (i < length(args)) args[i + 1] else ""
    if (key == "--metric") { out$metric <- val; i <- i + 2; next }
    if (key == "--out")    { out$out <- val;    i <- i + 2; next }
    if (key == "--no-cluster-rows") { out$cluster_rows <- FALSE; i <- i + 1; next }
    if (key == "--no-cluster-cols") { out$cluster_cols <- FALSE; i <- i + 1; next }
    if (key == "--width")  { out$width <- as.numeric(val);  i <- i + 2; next }
    if (key == "--height") { out$height <- as.numeric(val); i <- i + 2; next }
    i <- i + 1
  }
  out
}

ensure_dir <- function(path) {
  d <- dirname(path)
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# -------- locate project root robustly --------
script_dir <- get_script_dir()

# Assumes this script lives in: <PROJECT_ROOT>/scripts/pa_pregnancy/
project_root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = FALSE)
setwd(project_root)

# -------- inputs --------
args <- parse_args()
metric <- tolower(args$metric)

merged_dir <- "results/snp_only/merged_by_trait"
snplist_path <- "inputs/rsids_pa.txt"
plots_dir <- "results/snp_only/plots"

if (!dir.exists(merged_dir)) stop("Missing directory: ", merged_dir)
if (!file.exists(snplist_path)) stop("Missing SNP list: ", snplist_path)

files <- list.files(merged_dir, pattern = "_ALLCHR\\.tsv$", full.names = TRUE)
if (length(files) == 0) stop("No *_ALLCHR.tsv files found in: ", merged_dir)

snps <- fread(snplist_path, header = FALSE)$V1

# Default output name
if (args$out == "") {
  if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
  args$out <- file.path(plots_dir, paste0("heatmap_", metric, ".pdf"))
}
ensure_dir(args$out)

# -------- read & combine --------
dat_list <- lapply(files, function(f) {
  dt <- fread(f)
  trait <- sub("_ALLCHR\\.tsv$", "", basename(f))
  dt[, trait := trait]
  dt
})
dat <- rbindlist(dat_list, fill = TRUE)

# Keep only additive model rows (recommended)
if ("TEST" %in% names(dat)) dat <- dat[TEST == "ADD"]

# Coerce numeric safely (PLINK can sometimes write "." etc.)
if ("BETA" %in% names(dat)) dat[, BETA := suppressWarnings(as.numeric(BETA))]
if ("SE" %in% names(dat))   dat[, SE   := suppressWarnings(as.numeric(SE))]
if ("P" %in% names(dat))    dat[, P    := suppressWarnings(as.numeric(P))]

# Filter to SNP list
if (!("ID" %in% names(dat))) stop("Expected column 'ID' not found in merged files.")
dat <- dat[ID %in% snps]

if (nrow(dat) == 0) stop("After filtering to rsids_pa.txt, there are 0 rows to plot. (Check that IDs match.)")

# Choose value for heatmap
value_col <- NULL
if (metric == "beta") {
  value_col <- "BETA"
} else if (metric == "z") {
  if (!all(c("BETA", "SE") %in% names(dat))) stop("Need BETA and SE columns for metric=z")
  dat[, Z := BETA / SE]
  value_col <- "Z"
} else if (metric == "signedlogp") {
  if (!all(c("BETA", "P") %in% names(dat))) stop("Need BETA and P columns for metric=signedlogp")
  dat[, signedlogp := sign(BETA) * (-log10(P))]
  value_col <- "signedlogp"
} else {
  stop("Unknown --metric: ", metric, " (use beta | z | signedlogp)")
}

# Make trait × SNP matrix
wide <- dcast(dat, trait ~ ID, value.var = value_col)
mat <- as.matrix(wide[, -1, with = FALSE])
rownames(mat) <- wide$trait

# Optional: order SNP columns as in snplist (nice for readability)
common_snps <- intersect(snps, colnames(mat))
mat <- mat[, common_snps, drop = FALSE]

# -------- plot --------
main_title <- switch(
  metric,
  beta = "SNP × Trait effects (BETA)",
  z = "SNP × Trait effects (Z = BETA/SE)",
  signedlogp = "SNP × Trait signed -log10(P)",
  "SNP × Trait heatmap"
)

# Choose output device based on extension
ext <- tolower(file_ext(args$out))
if (ext %in% c("pdf")) {
  pdf(args$out, width = args$width, height = args$height)
} else if (ext %in% c("png")) {
  png(args$out, width = args$width, height = args$height, units = "in", res = 300)
} else {
  stop("Unsupported output extension: .", ext, " (use .pdf or .png)")
}

pheatmap(
  mat,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  na_col = "grey90",
  cluster_rows = args$cluster_rows,
  cluster_cols = args$cluster_cols,
  fontsize = 9,
  main = main_title
)

dev.off()

cat("Wrote heatmap to:", args$out, "\n")
cat("Project root:", project_root, "\n")
cat("Files read:", length(files), "\n")
cat("Traits:", nrow(mat), " SNPs:", ncol(mat), " Metric:", metric, "\n")
