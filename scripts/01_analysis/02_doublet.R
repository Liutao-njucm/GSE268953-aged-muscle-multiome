# 02_doublet.R — per-sample doublet removal via scDblFinder (Seurat5-compatible)
lib <- "C:/Users/wang/R/win-library/4.5"
.libPaths(c(lib, .libPaths()))
suppressMessages({
  library(Seurat); library(Signac); library(scDblFinder)
})
setwd("D:/workbuddy/outputs/2026-10-25-20-23-31")
out.dir <- "results/GSE268953"

obj.list <- readRDS(file.path(out.dir, "samples_qc_list.rds"))
clean.list <- list()

for (s in names(obj.list)) {
  cat("==== scDblFinder:", s, "====\n"); flush.console()
  obj <- DietSeurat(obj.list[[s]], assays = "RNA")
  # scDblFinder works on raw counts SCE-like; run on Seurat object directly
  sce <- as.SingleCellExperiment(obj, assay = "RNA")
  sce <- scDblFinder(sce, dbr = NULL, verbose = FALSE)
  cls <- sce$scDblFinder.class
  names(cls) <- colnames(sce)
  tab <- table(cls)
  cat(s, ": Singlet =", tab[["singlet"]], " Doublet =", tab[["doublet"]], "\n")
  obj <- subset(obj, cells = names(cls)[cls == "singlet"])
  clean.list[[s]] <- obj
  flush.console()
}

saveRDS(clean.list, file.path(out.dir, "samples_doublet_free.rds"))
tot <- sapply(clean.list, ncol)
print(tot)
cat("TOTAL SINGLET NUCLEI:", sum(tot), "\n")
cat("PHASE1B DONE\n")
