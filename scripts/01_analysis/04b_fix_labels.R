# 04b_fix_labels.R
.libPaths(c("C:/Users/wang/R/win-library/4.5", .libPaths()))
suppressMessages({ library(Seurat); library(ggplot2) })
setwd("D:/workbuddy/outputs/2026-10-25-20-23-31")
out.dir <- "results/GSE268953"
merged <- readRDS(file.path(out.dir, "merged_annotated.rds"))
ct <- gsub("^Score_|[0-9]+$", "", as.character(merged$celltype))
ct[ct == "SMC_Pericyte"] <- "SMC/Pericyte"
ct[ct == "T_NK"] <- "T/NK cells"
ct[ct == "Schwann_Glia"] <- "Schwann/Glia"
merged$celltype <- factor(ct)
prop <- prop.table(table(merged$celltype, merged$group), margin = 2)
print(round(prop * 100, 1))
write.csv(as.data.frame(prop), file.path(out.dir, "celltype_proportions.csv"))
p1 <- DimPlot(merged, group.by = "celltype", label = TRUE, repel = TRUE) +
  ggtitle("Cell types (GSE268953, n=42172)")
ggsave(file.path(out.dir, "umap_celltype.png"), p1, width = 10, height = 7, dpi = 150)
p2 <- DimPlot(merged, group.by = "group") + ggtitle("Young vs Aged")
ggsave(file.path(out.dir, "umap_group.png"), p2, width = 8, height = 6, dpi = 150)
saveRDS(merged, file.path(out.dir, "merged_annotated.rds"))
cat("FIXED\n")
