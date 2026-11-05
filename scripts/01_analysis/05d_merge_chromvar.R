# ============================================================
# Phase 2A - Step 1d: merge per-sample chromVAR results
# Input : atac_tmp/dev_<sample>.rds
# Output: chromvar_deviations.rds (motifs x cells)
#         atac_pseudobulk_per_sample.rds (pb list + peak GRanges list)
# ============================================================
.libPaths(c("C:/Users/wang/R/win-library/4.5", .libPaths()))
setwd("D:/workbuddy/outputs/2026-10-25-20-23-31/results/GSE268953")
files <- list.files("atac_tmp", pattern = "^dev_.*\\.rds$", full.names = TRUE)
samples <- gsub("^dev_|\\.rds$", "", basename(files))
cat("samples found:", paste(samples, collapse = ","), "\n")
objs <- setNames(lapply(files, readRDS), samples)

common <- Reduce(intersect, lapply(objs, function(x) rownames(x$dev)))
cat("motifs shared by all samples:", length(common), "\n")
dev.mat <- do.call(cbind, lapply(objs, function(x) x$dev[common, , drop = FALSE]))
cat("deviations:", nrow(dev.mat), "x", ncol(dev.mat), "\n")
saveRDS(dev.mat, "chromvar_deviations.rds")

saveRDS(list(pb = lapply(objs, `[[`, "pb"),
             gr = lapply(objs, `[[`, "gr")),
        "atac_pseudobulk_per_sample.rds")
cat("MERGE_DONE\n")
