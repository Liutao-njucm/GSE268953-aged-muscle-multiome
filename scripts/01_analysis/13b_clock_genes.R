.libPaths(c("C:/Users/wang/R/win-library/4.5", .libPaths()))
suppressMessages(library(data.table))
out.dir <- "D:/workbuddy/outputs/2026-10-25-20-23-31/results/GSE268953"
deg.dir <- file.path(out.dir, "deg")

clock <- c("ARNTL","ARNTL2","CLOCK","NPAS2","PER1","PER2","PER3","CRY1","CRY2","NR1D1","NR1D2",
           "DBP","TEF","HLF","NFIL3","CIART","RORA","NR2F6","BHLHE40","BHLHE41","SIRT1","NAMPT")
cts <- c("FAP","EC","Myonuclei","MuSC","SMC_Pericyte")
res <- list()
for (ct in cts) {
  f <- file.path(deg.dir, paste0("DEG_", ct, ".csv"))
  if (!file.exists(f)) next
  d <- fread(f)
  d <- d[gene %in% clock]
  res[[ct]] <- d[, .(celltype = ct, gene, log2FoldChange, padj)]
}
r <- rbindlist(res)
w <- dcast(r, gene ~ celltype, value.var = "log2FoldChange")
cat("=== clock gene log2FC (Aged vs Young) ===\n"); print(as.data.frame(w), digits = 2, row.names = FALSE)

p <- dcast(r, gene ~ celltype, value.var = "padj")
cat("\n=== padj ===\n"); print(format(as.data.frame(p), digits = 2), row.names = FALSE)
fwrite(r, file.path(out.dir, "clock_genes_by_celltype.csv"))
