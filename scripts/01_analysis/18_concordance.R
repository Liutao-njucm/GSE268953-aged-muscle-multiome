# ============================================================
# 18: cross-dataset concordance — snRNA FAP log2FC vs GTEx bulk
#     skeletal muscle ageing log2FC
# ============================================================
.libPaths(c("C:/Users/wang/R/win-library/4.5", .libPaths()))
suppressMessages({library(data.table); library(ggplot2); library(ggrepel)})

proj <- "D:/workbuddy/outputs/2026-10-25-20-23-31"
out.dir <- file.path(proj, "results/GSE268953")

gx <- fread(file.path(out.dir, "gtex/gtex_muscle_validation.csv"))
rna <- fread(file.path(out.dir, "deg/DEG_FAP.csv"))
clock <- c("NR1D1","NR1D2","PER1","PER2","PER3","CRY1","CRY2","DBP","TEF","HLF","NFIL3","RORA","BHLHE40","NR2F1","NR2F2")

m <- merge(gx[lfc_aged_vs_young != ""], rna[, .(gene, sn_lfc = log2FoldChange, sn_padj = padj)],
           by.x = "gene", by.y = "gene")
m[, lfc_aged_vs_young := as.numeric(lfc_aged_vs_young)]
m <- m[!is.na(sn_lfc)]
cat("shared genes:", nrow(m), "\n")
rho <- suppressWarnings(cor(m$sn_lfc, m$lfc_aged_vs_young, method = "spearman"))
ct <- suppressWarnings(cor.test(m$sn_lfc, m$lfc_aged_vs_young, method = "spearman"))
cat(sprintf("Spearman rho = %.3f, p = %.3g\n", rho, ct$p.value))

m[, set := fifelse(gene %in% clock, "circadian clock",
            fifelse(grepl("IL6|CXCL|CCL[0-9]|MMP|SERPINE1|IGFBP7|TGFB1|CSF1|CDKN|GLB1", gene), "SASP/senescence",
            "other"))]
fwrite(m, file.path(out.dir, "gtex_snRNA_concordance.csv"))

pal <- c("circadian clock" = "#8E44AD", "SASP/senescence" = "#E67E22", other = "grey60")
p <- ggplot(m, aes(sn_lfc, lfc_aged_vs_young, colour = set)) +
  geom_point(size = 2.6, alpha = 0.85) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey50") +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
  ggrepel::geom_text_repel(aes(label = gene), size = 2.9, max.overlaps = 18, show.legend = FALSE) +
  scale_colour_manual(values = pal, name = NULL) +
  labs(x = "log2FC aged vs young (FAP, snMultiome, 5v5 donors)",
       y = "log2FC aged vs young (GTEx skeletal muscle, bulk RNA-seq)",
       title = "External validation in GTEx (n = 818 donors, 20-39 vs 60-79 y)",
       subtitle = sprintf("Spearman rho = %.2f, p = %.2g across %d shared genes", rho, ct$p.value, nrow(m))) +
  theme_bw(base_size = 11)
ggsave(file.path(out.dir, "figures", "gtex_concordance_FAP.png"), p, width = 7.4, height = 5.6, dpi = 300)
cat("saved gtex_concordance_FAP.png\n")
