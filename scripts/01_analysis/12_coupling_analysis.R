# ============================================================
# 12: chromatin-transcriptome coupling (per cell type)
#   peak-level DA  ->  gene-level chromatin score
#   gene-level RNA  ->  pseudobulk DESeq2 LFC
#   classify genes into coupled / decoupled states
# ============================================================
.libPaths(c("C:/Users/wang/R/win-library/4.5", .libPaths()))
suppressMessages({library(data.table); library(ggplot2); library(ggrepel)})

out.dir <- "D:/workbuddy/outputs/2026-10-25-20-23-31/results/GSE268953"
deg.dir <- file.path(out.dir, "deg")
fig.dir <- file.path(out.dir, "figures"); dir.create(fig.dir, showWarnings = FALSE)

ann <- readRDS(file.path(out.dir, "peak_annotation.rds"))
da  <- fread(file.path(out.dir, "DA_peaks_Aged_vs_Young_by_celltype.csv"))
setDT(ann)

## link peaks to genes (drop far intergenic)
ann2 <- ann[!is.na(linked_gene) & region_class != "Intergenic(>100kb)",
            .(peak_id, linked_gene, region_class, dist_to_tss)]
da <- merge(da, ann2, by = "peak_id")
cat("DA rows with gene link:", nrow(da), "\n")

## significant DA peaks
da_sig <- da[!is.na(padj) & padj < 0.05 & abs(log2FoldChange) > 1]
cat("significant DA peaks:", nrow(da_sig), "\n")

## gene-level chromatin score from significant peaks
chrom <- da_sig[, .(n_peak_up = sum(log2FoldChange > 0),
                    n_peak_down = sum(log2FoldChange < 0),
                    chrom_score = mean(log2FoldChange),
                    max_abs_lfc = max(abs(log2FoldChange)),
                    n_promoter = sum(grepl("Promoter", region_class))),
                by = .(celltype, gene = linked_gene)]

cpl.list <- list(); summ.list <- list()
classes_all <- list()

for (ct in unique(chrom$celltype)) {
  f <- file.path(deg.dir, paste0("DEG_", gsub("[/ ]", "_", ct), ".csv"))
  if (!file.exists(f)) next
  rna <- fread(f)
  ch <- chrom[celltype == ct]
  m <- merge(ch, rna[, .(gene, rna_lfc = log2FoldChange, rna_padj = padj, rna_baseMean = baseMean)],
             by = "gene", all.x = TRUE)
  ## also bring in RNA-significant genes without any linked DA peak
  only_rna <- rna[!gene %in% m$gene & !is.na(padj) & padj < 0.05 & abs(log2FoldChange) > 1,
                  .(gene, rna_lfc = log2FoldChange, rna_padj = padj, rna_baseMean = baseMean)]
  if (nrow(only_rna)) {
    m <- rbind(m, data.table(celltype = ct, gene = only_rna$gene,
                             n_peak_up = 0L, n_peak_down = 0L,
                             chrom_score = 0, max_abs_lfc = 0, n_promoter = 0L,
                             rna_lfc = only_rna$rna_lfc, rna_padj = only_rna$rna_padj,
                             rna_baseMean = only_rna$rna_baseMean), fill = TRUE)
  }
  m$chrom_sig <- m$n_peak_up + m$n_peak_down > 0
  m$rna_sig <- !is.na(m$rna_padj) & m$rna_padj < 0.05 & abs(m$rna_lfc) > 0.5

  m$class <- "Unchanged"
  m[chrom_sig & chrom_score > 0 & rna_sig & rna_lfc > 0, class := "Coupled_up"]
  m[chrom_sig & chrom_score < 0 & rna_sig & rna_lfc < 0, class := "Coupled_down"]
  m[chrom_sig & chrom_score > 0 & !(rna_sig & rna_lfc > 0), class := "Chromatin_primed"]   # open but not transcribed
  m[chrom_sig & chrom_score < 0 & !(rna_sig & rna_lfc < 0), class := "Chromatin_poised_down"]
  m[!chrom_sig & rna_sig & rna_lfc > 0, class := "RNA_only_up"]
  m[!chrom_sig & rna_sig & rna_lfc < 0, class := "RNA_only_down"]

  ## concordance statistics (genes with both signals measurable)
  sub <- m[!is.na(rna_lfc)]
  tab <- table(factor(sub$chrom_score > 0, c(FALSE, TRUE)),
               factor(sub$rna_lfc > 0, c(FALSE, TRUE)))
  ft <- tryCatch(fisher.test(tab), error = function(e) NULL)
  rho <- suppressWarnings(cor(sub$chrom_score, sub$rna_lfc, method = "spearman"))

  ## among chromatin-significant genes only: how often does RNA follow?
  cs <- m[chrom_sig == TRUE]
  conc_up <- mean(cs$rna_lfc > 0 & cs$rna_sig, na.rm = TRUE)
  conc_dn <- mean(cs$rna_lfc < 0 & cs$rna_sig, na.rm = TRUE)

  summ.list[[ct]] <- data.table(
    celltype = ct,
    n_genes_linked = nrow(m),
    n_chrom_sig_genes = sum(m$chrom_sig),
    n_rna_sig_genes = sum(m$rna_sig),
    n_Coupled_up = sum(m$class == "Coupled_up"),
    n_Coupled_down = sum(m$class == "Coupled_down"),
    n_Chromatin_primed = sum(m$class == "Chromatin_primed"),
    n_Chromatin_poised_down = sum(m$class == "Chromatin_poised_down"),
    n_RNA_only_up = sum(m$class == "RNA_only_up"),
    n_RNA_only_down = sum(m$class == "RNA_only_down"),
    fisher_OR = ifelse(is.null(ft), NA, unname(ft$estimate)),
    fisher_p = ifelse(is.null(ft), NA, ft$p.value),
    spearman_rho = rho,
    frac_chrom_up_transcribed = conc_up,
    frac_chrom_down_transcribed = conc_dn
  )
  cpl.list[[ct]] <- m
  classes_all[[ct]] <- as.data.table(table(m$class))
}

cpl <- rbindlist(cpl.list, fill = TRUE)
summ <- rbindlist(summ.list, fill = TRUE)
setorder(summ, -n_Coupled_up)

fwrite(cpl, file.path(out.dir, "coupling_gene_level_by_celltype.csv"))
fwrite(summ, file.path(out.dir, "coupling_summary_by_celltype.csv"))
cat("\n=== coupling summary ===\n"); print(as.data.frame(summ))

## ---------------- figures ----------------
pal <- c(Coupled_up = "#C0392B", Coupled_down = "#1F618D",
         Chromatin_primed = "#E67E22", Chromatin_poised_down = "#7FB3D5",
         RNA_only_up = "#F5B7B1", RNA_only_down = "#AED6F1", Unchanged = "grey85")

for (ct in c("FAP", "Myonuclei", "MuSC", "EC")) {
  m <- cpl[celltype == ct]
  if (!nrow(m)) next
  m[, label := ifelse(class %in% c("Coupled_up", "Coupled_down", "Chromatin_primed") &
                        max_abs_lfc >= 1.5, gene, "")]
  p <- ggplot(m, aes(chrom_score, rna_lfc)) +
    geom_point(aes(colour = class), size = 1.6, alpha = 0.75) +
    geom_hline(yintercept = 0, linetype = 2, colour = "grey40") +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey40") +
    scale_colour_manual(values = pal, name = NULL) +
    labs(x = "chromatin accessibility change (mean log2FC of linked DA peaks)",
         y = "RNA expression change (donor pseudobulk log2FC)",
         title = paste0(ct, ": chromatin vs transcriptome in ageing"),
         subtitle = sprintf("Spearman rho = %.3f | Fisher OR = %.2f (p = %.3g)",
                            summ[celltype == ct, spearman_rho],
                            summ[celltype == ct, fisher_OR],
                            summ[celltype == ct, fisher_p])) +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom", legend.text = element_text(size = 8))
  if (any(m$label != "")) p <- p + geom_text_repel(aes(label = label), size = 2.6, max.overlaps = 25)
  ggsave(file.path(fig.dir, paste0("coupling_", gsub("[/ ]", "_", ct), ".png")),
         p, width = 7.2, height = 6, dpi = 300)
}

## stacked bar of classes
bar <- cpl[, .N, by = .(celltype, class)]
bar[, frac := N / sum(N), by = celltype]
p2 <- ggplot(bar[class != "Unchanged"], aes(celltype, frac, fill = class)) +
  geom_col(width = 0.72) +
  scale_fill_manual(values = pal, name = NULL) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = NULL, y = "fraction of genes with a linked DA peak",
       title = "Chromatin-transcriptome coupling classes") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "right")
ggsave(file.path(fig.dir, "coupling_class_fraction.png"), p2, width = 7.5, height = 4.4, dpi = 300)

cat("\nsaved coupling tables + figures\n")
