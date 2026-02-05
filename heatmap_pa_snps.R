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
OUT_PDF    <- "results/snp_only/plots/heatmap_beta.pdf"
dir.create(dirname(OUT_PDF), recursive = TRUE, showWarnings = FALSE)

snps <- fread(SNP_FILE, header = FALSE)$V1

files <- list.files(MERGED_DIR, pattern = "_ALLCHR\\.tsv$", full.names = TRUE)
if (length(files) == 0) stop("No merged *_ALLCHR.tsv files found in: ", MERGED_DIR)

dat_list <- lapply(files, function(f) {
  dt <- fread(f)
  trait <- sub("_ALLCHR\\.tsv$", "", basename(f))
  dt[, trait := trait]
  dt
})
dat <- rbindlist(dat_list, fill = TRUE)

# Filter + numeric conversion
dat <- dat[TEST == "ADD"]
dat[, ID := as.character(ID)]
dat[, BETA := suppressWarnings(as.numeric(BETA))]
dat <- dat[ID %in% snps]

if (nrow(dat) == 0) stop("No rows left after filtering to SNP list. Check ID matching.")

# Trait × SNP matrix
wide <- dcast(dat, trait ~ ID, value.var = "BETA")
mat <- as.matrix(wide[, -1, with = FALSE])
rownames(mat) <- wide$trait

# Order SNP columns as in SNP list
common_snps <- intersect(snps, colnames(mat))
mat <- mat[, common_snps, drop = FALSE]

# ---- Clean matrix for clustering ----
# Drop rows/cols completely missing
keep_rows <- rowSums(!is.na(mat)) > 0
keep_cols <- colSums(!is.na(mat)) > 0
mat <- mat[keep_rows, keep_cols, drop = FALSE]

# Replace NA/NaN/Inf with 0 so dist()/hclust() works
mat[!is.finite(mat)] <- NA
mat[is.na(mat)] <- 0

# Plot
pdf(OUT_PDF, width = 12, height = 9)
heatmap(
  mat,
  scale = "none",
  col = colorRampPalette(c("blue", "white", "red"))(100),
  margins = c(10, 12),
  main = "SNP × Trait effects (BETA)"
)
dev.off()

cat("Heatmap written to:", OUT_PDF, "\n")
cat("Traits:", nrow(mat), " SNPs:", ncol(mat), "\n")

