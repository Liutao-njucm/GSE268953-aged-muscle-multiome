# 03_integrate.R — merge, Harmony integration, clustering, marker genes
lib <- "C:/Users/wang/R/win-library/4.5"
.libPaths(c(lib, .libPaths()))
suppressMessages({
  library(Seurat); library(harmony); library(presto)
})
setwd("D:/workbuddy/outputs/2026-10-25-20-23-31")
out.dir <- "results/GSE268953"

clean.list <- readRDS(file.path(out.dir, "samples_doublet_free.rds"))
clean.list <- lapply(clean.list, function(o) DietSeurat(o, assays = "RNA"))
merged <- merge(clean.list[[1]], y = clean.list[-1], add.cell.ids = names(clean.list))
merged <- JoinLayers(merged)

merged <- NormalizeData(merged, verbose = FALSE)
merged <- FindVariableFeatures(merged, nfeatures = 3000, verbose = FALSE)
merged <- ScaleData(merged, vars.to.regress = c("percent.mt", "nCount_RNA"),
                    verbose = FALSE)
merged <- RunPCA(merged, npcs = 50, verbose = FALSE)
merged <- RunHarmony(merged, group.by.vars = "sample", verbose = FALSE)
merged <- RunUMAP(merged, reduction = "harmony", dims = 1:30, verbose = FALSE)
merged <- FindNeighbors(merged, reduction = "harmony", dims = 1:30, verbose = FALSE)
for (r in c(0.2, 0.4, 0.6)) {
  merged <- FindClusters(merged, resolution = r, verbose = FALSE)
}
cat("clusters (res 0.4):", length(unique(merged$RNA_snn_res.0.4)), "\n")

cat("finding markers (presto)...\n"); flush.console()
Idents(merged) <- "RNA_snn_res.0.4"
markers <- FindAllMarkers(merged, only.pos = TRUE, min.pct = 0.25,
                          logfc.threshold = 0.25, verbose = FALSE)
write.csv(markers, file.path(out.dir, "markers_res04.csv"), row.names = FALSE)

top5 <- markers[markers$avg_log2FC > 1, ]
top5 <- top5[order(top5$cluster, -top5$avg_log2FC), ]
agg <- aggregate(avg_log2FC ~ cluster, data = top5,
                 FUN = function(x) paste(head(x[order(-x)], 5), collapse = ", "))
genes <- aggregate(gene ~ cluster, data = top5,
                   FUN = function(x) paste(head(x, 8), collapse = ", "))
print(merge(agg, genes, by = "cluster"))

# canonical marker mean expression per cluster for manual annotation
mk <- c("MYH7","MYH2","MYH1","MYH4","TNNT3","CKM","TTN",
        "PAX7","MYOD1","MYF5","MKI67",
        "DCN","PDGFRA","COL1A1","COL6A3","ADIPOQ",
        "PECAM1","VWF","KDR","ACTA2","TAGLN","RGS5",
        "PTPRC","LYZ","C1QA","CD68","CD3D","NKG7",
        "PLP1","MPZ","SCX","TNMD")
mk <- intersect(mk, rownames(merged))
Idents(merged) <- "RNA_snn_res.0.4"
cm <- AverageExpression(merged, features = mk, layer = "data")$RNA
write.csv(cm, file.path(out.dir, "cluster_marker_means.csv"))

saveRDS(merged, file.path(out.dir, "merged_harmony.rds"))
cat("TOTAL NUCLEI:", ncol(merged), "\n")
cat("PHASE1C DONE\n")
