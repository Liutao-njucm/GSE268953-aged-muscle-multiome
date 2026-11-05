## 29a_prep_regulon_input.R
## Prepare donor-level log2 CPM pseudobulk matrices for TF regulon activity analysis
## Input : results/GSE268953/deg/pseudobulk_<celltype>.rds  (list: res / pb counts / coldata)
## Output: results/GSE268953/regulon/pb_logcpm_<celltype>.csv  + donor_group.csv

.libPaths(c("C:/Users/wang/R/win-library/4.5", .libPaths()))
suppressPackageStartupMessages(library(Matrix))

root <- "D:/workbuddy/outputs/2026-10-25-20-23-31"
degdir <- file.path(root, "results/GSE268953/deg")
outdir <- file.path(root, "results/GSE268953/regulon")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

cts <- c("FAP", "Myonuclei", "MuSC", "EC", "T_NK_cells", "Macrophage", "SMC_Pericyte")

grp_ref <- NULL
for (ct in cts) {
  f <- file.path(degdir, paste0("pseudobulk_", ct, ".rds"))
  if (!file.exists(f)) { message("missing: ", f); next }
  x <- readRDS(f)
  pb <- x$pb
  cd <- x$coldata
  ## group vector: coldata may be a 1-col df or a vector
  grp <- if (is.data.frame(cd)) cd[[1]] else cd
  names(grp) <- rownames(cd)
  if (is.null(grp_ref)) grp_ref <- data.frame(donor = colnames(pb), celltype = ct, group = grp,
                                              stringsAsFactors = FALSE)
  ## CPM then log2
  lib <- colSums(pb)
  cpm <- sweep(pb, 2, lib, "/") * 1e6
  lc <- log2(cpm + 1)
  ## drop genes never expressed in this cell type
  keep <- rowSums(pb >= 5) >= 2
  lc <- lc[keep, , drop = FALSE]
  write.csv(data.frame(gene = rownames(lc), lc, check.names = FALSE),
            file.path(outdir, paste0("pb_logcpm_", ct, ".csv")), row.names = FALSE)
  cat(sprintf("%-14s genes=%5d donors=%2d  libsize=%.2f-%.2fM\n", ct, nrow(lc), ncol(lc),
              min(lib)/1e6, max(lib)/1e6))
}

## consolidated donor table (donor id is shared across cell types)
donors <- unique(grp_ref[, c("donor", "group")])
write.csv(donors, file.path(outdir, "donor_group.csv"), row.names = FALSE)
cat("\ndonor table:\n"); print(donors)
cat("\nDONE ->", outdir, "\n")
