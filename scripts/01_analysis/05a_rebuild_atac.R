# ============================================================
# Phase 2A - Step 0 (v2): rebuild multiome object, fast path
# - RNA: cbind (same 36601 genes across samples)
# - ATAC: map per-sample peaks onto union peak set, then cbind
# Log: scripts/05a_log.txt
# ============================================================
.libPaths(c("C:/Users/wang/R/win-library/4.5", .libPaths()))
suppressMessages({library(Seurat); library(Matrix)})
setwd("D:/workbuddy/outputs/2026-10-25-20-23-31/results/GSE268953")
logf <- file("../../scripts/05a_log.txt", open = "wt")
sink(logf, split = TRUE); sink(logf, type = "message")
t0 <- Sys.time(); cat("start:", format(t0), "\n")

cat("loading merged_annotated for metadata...\n"); flush.console()
m <- readRDS("merged_annotated.rds")
md <- data.frame(barcode = colnames(m),
                 celltype = as.character(m$celltype),
                 group = as.character(m$group),
                 sample = as.character(m@meta.data$sample),
                 stringsAsFactors = FALSE)
rm(m); invisible(gc())
write.csv(md, "cell_metadata.csv", row.names = FALSE)
cat("metadata ok:", nrow(md), "cells\n"); flush.console()

cat("loading samples_qc_list...\n"); flush.console()
qcl <- readRDS("samples_qc_list.rds")
cat("samples:", paste(names(qcl), collapse = ","), "\n"); flush.console()

# ---- per-sample: subset cells, collect matrices + peak ranges ----
rna.list <- list(); atac.list <- list(); peak.gr <- list()
for (s in names(qcl)) {
  o <- qcl[[s]]
  sel <- md$barcode[md$sample == s]
  keep <- gsub(paste0("^", s, "_"), "", sel)
  keep <- intersect(keep, colnames(o))
  cat(s, "singlets found:", length(keep), "/", length(sel), "\n"); flush.console()
  o <- subset(o, cells = keep)
  rna.list[[s]] <- GetAssayData(o[["RNA"]], layer = "counts")
  colnames(rna.list[[s]]) <- paste(s, colnames(rna.list[[s]]), sep = "_")
  aa <- o[["ATAC"]]
  atac.list[[s]] <- GetAssayData(aa, layer = "counts")
  colnames(atac.list[[s]]) <- paste(s, colnames(atac.list[[s]]), sep = "_")
  peak.gr[[s]] <- Signac::granges(aa)
  rm(o, aa); invisible(gc())
  cat(s, "done\n"); flush.console()
}
rm(qcl); invisible(gc())

# ---- union peak set ----
cat("building union peaks...\n"); flush.console()
allpk <- unique(unlist(lapply(peak.gr, function(g) Signac::GRangesToString(g))))
upk <- do.call(rbind, strsplit(allpk, "-"))
upk.gr <- GenomicRanges::GRanges(upk[, 1],
        IRanges::IRanges(as.numeric(upk[, 2]), as.numeric(upk[, 3])))
names(upk.gr) <- allpk
cat("union peaks:", length(upk.gr), "\n"); flush.console()

# ---- map each sample's counts onto union peaks (exact coord match) ----
atac.u <- lapply(names(atac.list), function(s) {
  pk <- Signac::GRangesToString(peak.gr[[s]])
  idx <- match(pk, allpk)                       # union row for each original peak
  stopifnot(!any(is.na(idx)))
  M <- Matrix::sparseMatrix(i = idx, j = seq_along(idx), x = rep(1, length(idx)),
                            dims = c(length(allpk), length(idx)))
  M %*% atac.list[[s]]
})
names(atac.u) <- names(atac.list)
atac.mat <- do.call(cbind, atac.u)
rm(atac.list, atac.u); invisible(gc())
rna.mat <- do.call(cbind, rna.list)
rm(rna.list); invisible(gc())
cat("matrices: RNA", nrow(rna.mat), "x", ncol(rna.mat),
    "| ATAC", nrow(atac.mat), "x", ncol(atac.mat), "\n"); flush.console()

# ---- build one Seurat object ----
cat("creating Seurat object...\n"); flush.console()
obj <- CreateSeuratObject(counts = rna.mat, project = "GSE268953")
chrom <- CreateChromatinAssay(counts = atac.mat, genome = "hg38", ranges = upk.gr)
obj[["ATAC"]] <- chrom
md2 <- md[match(colnames(obj), md$barcode), ]
obj$celltype <- md2$celltype; obj$group <- md2$group; obj$sample <- md2$sample
cat("object:", ncol(obj), "nuclei | assays:", paste(Assays(obj), collapse = ","), "\n")
flush.console()
saveRDS(obj, "merged_multiome.rds")
cat("REBUILD_DONE in", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")
sink(type = "message"); sink(); close(logf)
