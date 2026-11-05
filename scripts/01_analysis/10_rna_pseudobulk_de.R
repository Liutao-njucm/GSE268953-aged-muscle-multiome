# ============================================================
# 10: RNA donor-level pseudobulk differential expression (5 vs 5)
#   - per cell type: sum raw counts across nuclei of each donor
#   - DESeq2 ~ group (Aged vs Young), Wald test
#   - outputs DEG tables, summary, and a compact object for
#     downstream chromatin-transcriptome coupling analysis
# ============================================================
.libPaths(c("C:/Users/wang/R/win-library/4.5", .libPaths()))
suppressMessages({library(Seurat); library(Matrix); library(DESeq2)})

proj <- "D:/workbuddy/outputs/2026-10-25-20-23-31"
out.dir <- file.path(proj, "results/GSE268953")
deg.dir <- file.path(out.dir, "deg"); dir.create(deg.dir, showWarnings = FALSE)

obj <- readRDS(file.path(out.dir, "merged_annotated.rds"))
DefaultAssay(obj) <- "RNA"
cnt <- LayerData(obj, assay = "RNA", layer = "counts")
meta <- obj@meta.data
meta$sample <- factor(as.character(meta$sample),
                      levels = c(paste0("Young_", 1:5), paste0("Aged_", 1:5)))
meta$group <- factor(ifelse(grepl("^Aged", as.character(meta$sample)), "Aged", "Young"),
                     levels = c("Young", "Aged"))
rm(obj); gc()

MIN_CELLS <- 30          # donor must contribute >=30 nuclei to a cell type
MIN_COUNT <- 10          # gene must reach >=10 counts ...
MIN_DONOR <- 5           # ... in at least 5 donors
celltypes <- sort(unique(as.character(meta$celltype)))
celltypes <- setdiff(celltypes, character(0))

summary_list <- list(); deg_all <- list()

for (ct in celltypes) {
  idx <- which(meta$celltype == ct)
  tab <- table(meta$sample[idx])
  keep.samples <- names(tab)[tab >= MIN_CELLS]
  if (length(keep.samples) < 6 ||
      sum(grepl("^Aged", keep.samples)) < 3 ||
      sum(grepl("^Young", keep.samples)) < 3) {
    cat(sprintf("skip %-14s : donors with >=%d nuclei = %d\n", ct, MIN_CELLS, length(keep.samples)))
    next
  }
  # donor x gene pseudobulk (sum of UMI)
  grp <- factor(ifelse(grepl("^Aged", keep.samples), "Aged", "Young"), levels = c("Young", "Aged"))
  pb <- sapply(keep.samples, function(s) {
    Matrix::rowSums(cnt[, idx[as.character(meta$sample[idx]) == s], drop = FALSE])
  })
  pb <- as.matrix(pb)
  # keep expressed genes
  keep.g <- rowSums(pb >= MIN_COUNT) >= MIN_DONOR
  pb <- pb[keep.g, , drop = FALSE]
  cat(sprintf("%-14s : %6d nuclei | %2d donors (%dY/%dA) | %5d genes tested\n",
              ct, length(idx), ncol(pb), sum(grp == "Young"), sum(grp == "Aged"), nrow(pb)))

  coldata <- data.frame(group = grp, row.names = colnames(pb))
  dds <- DESeqDataSetFromMatrix(round(pb), coldata, design = ~ group)
  dds <- DESeq(dds, quiet = TRUE)
  res <- results(dds, contrast = c("group", "Aged", "Young"))
  res <- as.data.frame(res)
  res$gene <- rownames(res); res$celltype <- ct
  res$baseMean_donor <- res$baseMean
  write.csv(res, file.path(deg.dir, paste0("DEG_", gsub("[/ ]", "_", ct), ".csv")), row.names = FALSE)
  deg_all[[ct]] <- res

  sig <- !is.na(res$padj) & res$padj < 0.05 & abs(res$log2FoldChange) > 1
  summary_list[[ct]] <- data.frame(
    celltype = ct, n_nuclei = length(idx), n_donors = ncol(pb),
    n_Young = sum(grp == "Young"), n_Aged = sum(grp == "Aged"),
    n_genes_tested = nrow(res),
    n_sig = sum(sig), n_up_Aged = sum(sig & res$log2FoldChange > 0),
    n_down_Aged = sum(sig & res$log2FoldChange < 0),
    n_sig_loose = sum(!is.na(res$padj) & res$padj < 0.1 & abs(res$log2FoldChange) > 0.58)
  )
  saveRDS(list(res = res, pb = pb, coldata = coldata),
          file.path(deg.dir, paste0("pseudobulk_", gsub("[/ ]", "_", ct), ".rds")))
}

summ <- do.call(rbind, summary_list)
rownames(summ) <- NULL
write.csv(summ, file.path(deg.dir, "DEG_summary_by_celltype.csv"), row.names = FALSE)
cat("\n=== summary ===\n"); print(summ)

all.res <- do.call(rbind, deg_all)
saveRDS(all.res, file.path(deg.dir, "DEG_all_celltypes.rds"))
write.csv(all.res, file.path(deg.dir, "DEG_all_celltypes.csv"), row.names = FALSE)
cat("\nsaved:", deg.dir, "\n")
