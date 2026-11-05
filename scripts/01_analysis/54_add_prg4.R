# 54_add_prg4.R — compute per-donor FAP CPM10k for PRG4 (+ validation genes)
#   CPM10k = donor FAP pseudobulk counts / total FAP counts * 1e4
#   Validation: IL6, CCL2, NR1D1 must reproduce existing FAP_clock_per_donor.csv values.
.libPaths(c("C:/Users/wang/R/win-library/4.5", .libPaths()))
suppressMessages({library(Seurat); library(data.table); library(Matrix)})
setwd("D:/workbuddy/outputs/2026-10-25-20-23-31")
obj <- readRDS("results/GSE268953/merged_annotated.rds")
cat("object loaded:", ncol(obj), "cells\n")
sub <- subset(obj, celltype == "FAP")
cat("FAP cells:", ncol(sub), "\n")
m <- tryCatch(GetAssayData(sub, assay = "RNA", slot = "counts"),
              error = function(e) GetAssayData(sub, layer = "counts"))
samp <- as.character(sub$sample)
genes <- c("PRG4", "IL6", "CCL2", "NR1D1", "MUSK", "PDGFRA")
genes <- intersect(genes, rownames(m))
cat("genes found:", genes, "\n")
tot <- Matrix::colSums(m)
res <- data.table(gene = character(), CPM10k = numeric(), donor = character(), group = character())
for (dn in sort(unique(samp))) {
  idx <- which(samp == dn)
  denom <- sum(tot[idx])
  for (g in genes) {
    cnt <- sum(m[g, idx])
    res <- rbind(res, data.table(gene = g, CPM10k = cnt / denom * 1e4, donor = dn,
                                 group = ifelse(grepl("^Young", dn), "Young", "Aged")))
  }
}
fwrite(res, "results/GSE268953/_prg4_check.csv")
cat("DONE\n")
