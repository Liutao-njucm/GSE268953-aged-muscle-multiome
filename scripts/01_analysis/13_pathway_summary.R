# ============================================================
# 13: pathway-level summary of the ageing transcriptome
#   curated gene sets (no external DB needed) + hypergeometric test
# ============================================================
.libPaths(c("C:/Users/wang/R/win-library/4.5", .libPaths()))
suppressMessages({library(data.table); library(ggplot2)})

out.dir <- "D:/workbuddy/outputs/2026-10-25-20-23-31/results/GSE268953"
deg.dir <- file.path(out.dir, "deg")
fig.dir <- file.path(out.dir, "figures")

gene_sets <- list(
  Circadian_core = c("ARNTL","ARNTL2","CLOCK","NPAS2","PER1","PER2","PER3","CRY1","CRY2",
                     "NR1D1","NR1D2","DBP","TEF","HLF","NFIL3","CIART","DEC1","DEC2","RORA","NR2F6"),
  Circadian_output = c("NAMPT","SIRT1","PPARGC1A","GYS1","SLC2A4","FABP3","PDK4","UCP3","MYOD1"),
  SASP_inflammatory = c("IL6","CXCL1","CXCL2","CXCL8","CXCL9","CXCL10","CCL2","CCL3","CCL4","CCL5",
                        "IL1B","IL1A","TNF","MMP1","MMP3","MMP9","MMP13","IGFBP3","IGFBP7","SERPINE1",
                        "TGFB1","CSF1","CSF2","ICAM1","VCAM1","NFKB1","RELA"),
  ECM_fibrosis = c("COL1A1","COL1A2","COL3A1","COL5A1","COL6A1","COL6A2","COL6A3","COL8A1","COL12A1",
                   "FN1","POSTN","LUM","DCN","FBN1","ELN","LOX","LOXL1","LOXL2","TIMP1","TIMP3",
                   "CTGF","CCN2","SPARC","THBS1","ADAMTS2","ADAMTS16","SERPINH1"),
  Adipogenic = c("PPARG","CEBPA","CEBPD","FABP4","ADIPOQ","PLIN1","PLIN4","LPL","DGAT2","CIDEC","LEP"),
  OXPHOS_mito = c("NDUFA1","NDUFA2","NDUFB1","NDUFB2","NDUFS1","NDUFS2","NDUFV1","NDUFV2",
                  "SDHA","SDHB","SDHC","SDHD","UQCRC1","UQCRC2","UQCRFS1","COX4I1","COX5A","COX6A1",
                  "COX7A1","ATP5F1A","ATP5F1B","ATP5MC1","ATP5MC2","ATP5PO","CYCS","TOMM20","TFAM"),
  Glycolysis = c("HK1","HK2","GPI","PFKL","PFKM","PFKP","ALDOA","GAPDH","PGK1","PGAM1","ENO1",
                 "ENO3","PKM","LDHA","LDHB","SLC2A1","SLC2A4","PDK1","PDK4"),
  Proteostasis = c("HSPA1A","HSPA1B","HSPA5","HSPA8","HSPB1","HSPB5","CRYAB","HSP90AA1","HSP90AB1",
                   "DNAJA1","DNAJB1","BAG3","STIP1","HSF1","VCP"),
  Autophagy_lysosome = c("MAP1LC3B","SQSTM1","BECN1","ATG5","ATG7","ATG12","ULK1","LAMP1","LAMP2",
                         "CTSB","CTSD","TFEB","TFE3","BNIP3","BNIP3L","FUNDC1"),
  Senescence_markers = c("CDKN1A","CDKN2A","CDKN2B","TP53","RB1","GLB1","SERPINE1","IL6","HMGB1","LMNB1"),
  Myogenesis = c("PAX7","PAX3","MYF5","MYOD1","MYOG","MYF6","MEF2A","MEF2C","DES","SIX1","SIX4","CHRDL2"),
  NMJ_synaptic = c("CHRNA1","CHRNB1","CHRND","CHRNE","MUSK","LRP4","DOK7","RAPSN","AGRN","NCAM1","SYT2","SCN4A"),
  Sarcomere = c("ACTA1","ACTN2","MYH1","MYH2","MYH7","TNNT1","TNNT3","TNNI1","TNNI2","TNNC1","TNNC2",
                "MYBPH","MYL1","MYL2","LMOD2","LMOD3","NEB","TTN","RYR1","CASQ1","ATP2A1"),
  Fibro_adipo_precursor = c("PDGFRA","PDGFRB","DCN","LUM","GSN","WISP2","CCN5","DLK1","PIEZO1","TCF21"),
  Immune_MHC = c("HLA-A","HLA-B","HLA-C","HLA-DRA","HLA-DRB1","HLA-E","B2M","TAP1","TAP2","PSMB8","PSMB9","CIITA")
)

summ.list <- list(); detail.list <- list()
fs <- list.files(deg.dir, pattern = "^DEG_.*\\.csv$")
fs <- fs[!grepl("summary|all_celltypes", fs)]
for (f in fs) {
  ct <- sub("^DEG_", "", sub("\\.csv$", "", f))
  rna <- fread(file.path(deg.dir, f))
  rna <- rna[!is.na(padj)]
  bg <- unique(rna$gene)
  sig_up <- rna[padj < 0.05 & log2FoldChange > 1, gene]
  sig_dn <- rna[padj < 0.05 & log2FoldChange < -1, gene]
  for (gs in names(gene_sets)) {
    g <- intersect(gene_sets[[gs]], bg)
    if (length(g) < 3) next
    hit_up <- sum(g %in% sig_up); hit_dn <- sum(g %in% sig_dn)
    p_up <- phyper(hit_up - 1, length(sig_up), length(bg) - length(sig_up), length(g), lower.tail = FALSE)
    p_dn <- phyper(hit_dn - 1, length(sig_dn), length(bg) - length(sig_dn), length(g), lower.tail = FALSE)
    lfc <- rna[gene %in% g, log2FoldChange]
    summ.list[[paste(ct, gs)]] <- data.table(
      celltype = ct, gene_set = gs, n_genes_in_set = length(g),
      n_up = hit_up, n_down = hit_dn,
      mean_lfc = mean(lfc, na.rm = TRUE), median_lfc = median(lfc, na.rm = TRUE),
      frac_up = mean(lfc > 0, na.rm = TRUE),
      p_up = p_up, p_down = p_dn,
      direction = ifelse(mean(lfc, na.rm = TRUE) > 0, "up", "down"),
      p_dir = ifelse(mean(lfc, na.rm = TRUE) > 0, p_up, p_dn)
    )
    detail.list[[paste(ct, gs)]] <- data.table(celltype = ct, gene_set = gs, gene = g, lfc = lfc)
  }
}
summ <- rbindlist(summ.list)
summ[, padj_dir := p.adjust(p_dir, "BH"), by = celltype]
fwrite(summ, file.path(out.dir, "pathway_summary_by_celltype.csv"))
det <- rbindlist(detail.list, fill = TRUE)
fwrite(det, file.path(out.dir, "pathway_gene_lfc_long.csv"))

## heatmap of mean LFC per celltype x gene set
mat <- dcast(summ, gene_set ~ celltype, value.var = "mean_lfc")
m2 <- as.matrix(mat[, -1]); rownames(m2) <- mat$gene_set
ord <- hclust(dist(m2))$order
long <- data.table(gene_set = rep(rownames(m2), times = ncol(m2)),
                   celltype = rep(colnames(m2), each = nrow(m2)),
                   mean_lfc = as.vector(m2))
long[, gene_set := factor(gene_set, levels = rownames(m2)[ord])]
sig <- summ[, .(gene_set, celltype, padj_dir)]
long <- merge(long, sig, by = c("gene_set", "celltype"), all.x = TRUE)
long[, star := ifelse(!is.na(padj_dir) & padj_dir < 0.05, "*", "")]

p <- ggplot(long, aes(celltype, gene_set, fill = mean_lfc)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = star), size = 4, vjust = 0.78, colour = "black") +
  scale_fill_gradient2(low = "#1F618D", mid = "white", high = "#C0392B",
                       midpoint = 0, name = "mean\nlog2FC") +
  labs(x = NULL, y = NULL, title = "Pathway-level shifts in aged muscle (donor pseudobulk)",
       subtitle = "* BH-adjusted directional hypergeometric p < 0.05") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1),
        panel.grid = element_blank())
ggsave(file.path(fig.dir, "pathway_heatmap.png"), p, width = 8.2, height = 6.4, dpi = 300)

cat("=== top directional gene-set shifts (padj_dir<0.05) ===\n")
print(as.data.frame(summ[padj_dir < 0.05][order(padj_dir)][1:min(30, .N)]))
cat("\nsaved pathway_summary_by_celltype.csv + pathway_heatmap.png\n")
