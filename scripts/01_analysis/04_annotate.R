# 04_annotate.R — rule-based scoring annotation + group proportions + plots
lib <- "C:/Users/wang/R/win-library/4.5"
.libPaths(c(lib, .libPaths()))
suppressMessages({ library(Seurat); library(ggplot2) })
setwd("D:/workbuddy/outputs/2026-10-25-20-23-31")
out.dir <- "results/GSE268953"

merged <- readRDS(file.path(out.dir, "merged_harmony.rds"))
Idents(merged) <- "RNA_snn_res.0.4"
cat("cluster sizes:\n"); print(table(Idents(merged)))
cat("cluster x group:\n"); print(table(Idents(merged), merged$group))

sets <- list(
  Myonuclei = c("TTN","TNNT3","CKM","MYBPC1","ACTN3","TNNT1"),
  MuSC      = c("PAX7","MYF5","MYOD1","CALCR","GREM1"),
  FAP       = c("DCN","PDGFRA","COL1A1","COL6A3","PDK4"),
  EC        = c("PECAM1","VWF","KDR","CDH5"),
  SMC_Pericyte = c("ACTA2","TAGLN","RGS5","MYH11","MYOCD"),
  Macrophage   = c("C1QA","C1QB","CD68","LYZ","FCGR3A"),
  T_NK         = c("CD3D","CD3E","NKG7","GNLY","PRF1"),
  Schwann_Glia = c("PLP1","MPZ","MBP","NRXN1"),
  Mast         = c("HDC","TPSAB1","CPA3","CTSG"),
  Proliferating= c("MKI67","TOP2A","AURKB")
)
sets <- lapply(sets, function(g) intersect(g, rownames(merged)))
for (nm in names(sets)) {
  merged <- AddModuleScore(merged, features = list(sets[[nm]]), name = paste0("Score_", nm))
}
score.cols <- grep("^Score_", colnames(merged@meta.data), value = TRUE)
names(sets) <- sub("^Score_", "", score.cols)   # keep mapping aligned
sc <- FetchData(merged, c(score.cols, "RNA_snn_res.0.4"))
sc$cluster <- sc$RNA_snn_res.0.4
mean.sc <- aggregate(sc[, score.cols], by = list(cluster = sc$cluster), FUN = mean)
rownames(mean.sc) <- as.character(mean.sc$cluster)
mean.sc$cluster <- NULL
assigned <- colnames(mean.sc)[max.col(mean.sc, ties.method = "first")]
names(assigned) <- rownames(mean.sc)
print(mean.sc); print(assigned)

merged$celltype <- unname(assigned[as.character(merged$RNA_snn_res.0.4)])
ct <- gsub("^Score_|[0-9]+$", "", merged$celltype)
ct[ct == "SMC_Pericyte"] <- "SMC/Pericyte"
ct[ct == "T_NK"] <- "T/NK cells"
ct[ct == "Schwann_Glia"] <- "Schwann/Glia"
merged$celltype <- factor(ct)
write.csv(mean.sc, file.path(out.dir, "cluster_celltype_scores.csv"))

# group proportions
prop <- prop.table(table(merged$celltype, merged$group), margin = 2)
prop.df <- as.data.frame(prop); colnames(prop.df) <- c("celltype","group","prop")
write.csv(prop.df, file.path(out.dir, "celltype_proportions.csv"), row.names = FALSE)
print(round(prop * 100, 1))
# group comparison per celltype
ct <- table(merged$celltype, merged$group)
write.csv(ct, file.path(out.dir, "celltype_counts.csv"))

# UMAP plots
p1 <- DimPlot(merged, group.by = "celltype", label = TRUE, repel = TRUE) +
  ggtitle("Cell types (GSE268953, n=42172)")
ggsave(file.path(out.dir, "umap_celltype.png"), p1, width = 10, height = 7, dpi = 150)
p2 <- DimPlot(merged, group.by = "group") + ggtitle("Young vs Aged")
ggsave(file.path(out.dir, "umap_group.png"), p2, width = 8, height = 6, dpi = 150)

saveRDS(merged, file.path(out.dir, "merged_annotated.rds"))
cat("PHASE1 ANNOTATION DONE\n")
