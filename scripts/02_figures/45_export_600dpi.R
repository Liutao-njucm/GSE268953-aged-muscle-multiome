## ============================================================
## 45_export_600dpi.R  —  re-render all Fig.1-5 panels at 600 dpi
##   and re-compose the 5 main figures as TIFF (LZW) + PNG.
##   Originals in submission/Figures/ are NOT touched.
##   Usage: Rscript 45_export_600dpi.R <stage>
##   stages: seurat | motif | csv | locus | compose
## ============================================================
.libPaths(c("C:/Users/wang/R/win-library/4.5", .libPaths()))
args <- commandArgs(trailingOnly = TRUE)
stage <- if (length(args)) args[1] else stop("give a stage: seurat|motif|csv|locus|compose")

proj    <- "D:/workbuddy/outputs/2026-10-25-20-23-31"
res     <- file.path(proj, "results/GSE268953")
pan.dir <- file.path(res, "panels_600dpi");   dir.create(pan.dir, showWarnings = FALSE)
out.dir <- file.path(res, "submission/Figures_600dpi"); dir.create(out.dir, showWarnings = FALSE)
DPI  <- 600
PG   <- function(pkg) suppressMessages(library(pkg, character.only = TRUE))

## ============================================================
## Stage 1: UMAP panels (needs the 1.4 GB Seurat object)
## ============================================================
if (stage == "seurat") {
  PG("Seurat"); PG("ggplot2")
  merged <- readRDS(file.path(res, "merged_annotated.rds"))
  ct <- gsub("^Score_|[0-9]+$", "", as.character(merged$celltype))
  ct[ct == "SMC_Pericyte"] <- "SMC/Pericyte"
  ct[ct == "T_NK"] <- "T/NK cells"
  ct[ct == "Schwann_Glia"] <- "Schwann/Glia"
  merged$celltype <- factor(ct)
  Idents(merged) <- "celltype"
  p1 <- DimPlot(merged, group.by = "celltype", label = TRUE, repel = TRUE) +
    ggtitle("Cell types (GSE268953, n=42172)")
  ggsave(file.path(pan.dir, "umap_celltype_600.png"), p1, width = 10, height = 7, dpi = DPI)
  p2 <- DimPlot(merged, group.by = "group") + ggtitle("Young vs Aged")
  ggsave(file.path(pan.dir, "umap_group_600.png"), p2, width = 8, height = 6, dpi = DPI)
  ## save coordinates once so future re-exports never need the big object again
  emb <- as.data.frame(Embeddings(merged, "umap"))
  emb$celltype <- as.character(merged$celltype)
  emb$group <- as.character(merged$group)
  emb$sample <- as.character(merged$orig.ident)
  write.csv(emb, file.path(res, "umap_coords_per_cell.csv"), row.names = TRUE)
  rm(merged); gc()
  cat("SEURAT STAGE DONE\n")
}

## ============================================================
## Stage 2: motif-activity panels (needs 423 MB motif matrix)
##   chromatin_aging_score_by_celltype + aging_TF_heatmap
## ============================================================
if (stage == "motif") {
  PG("Matrix"); PG("ggplot2"); PG("matrixStats")
  dat <- readRDS(file.path(res, "motif_activity_merged.rds"))
  z <- dat$z; meta <- dat$meta
  rm(dat); gc()
  sg <- read.csv(file.path(res, "aging_TF_signature.csv"), stringsAsFactors = FALSE)
  up <- intersect(sg$motif[sg$direction == "up_in_Aged"], rownames(z))
  dn <- intersect(sg$motif[sg$direction == "down_in_Aged"], rownames(z))

  ## Fig 2C: chromatin aging score by cell type (recompute from z, as in script 07)
  score <- colMeans(z[up, , drop = FALSE]) - colMeans(z[dn, , drop = FALSE])
  score <- (score - mean(score)) / sd(score)
  meta$aging_score <- score
  cs <- meta[, c("celltype", "sample", "aging_score")]
  agg <- aggregate(aging_score ~ celltype + sample, data = cs, FUN = mean)
  agg$group <- ifelse(grepl("^Aged", agg$sample), "Aged", "Young")
  ## keep cell types with >=3 donors per group (same rule as script 07)
  keep <- names(which(sapply(split(agg, agg$celltype), function(d)
    sum(d$group == "Aged") >= 3 && sum(d$group == "Young") >= 3)))
  agg <- agg[agg$celltype %in% keep, ]

  p1 <- ggplot(agg, aes(x = reorder(celltype, aging_score, FUN = mean),
                        y = aging_score, fill = group)) +
    geom_boxplot(outlier.size = 0.6, alpha = 0.85) +
    scale_fill_manual(values = c(Young = "#4C9BE8", Aged = "#E8604C")) +
    labs(x = NULL, y = "chromatin aging score (donor mean)",
         title = "ATAC-based chromatin aging score") +
    theme_bw(base_size = 12) + coord_flip()
  ggsave(file.path(pan.dir, "chromatin_aging_score_by_celltype_600.png"),
         p1, width = 8, height = 5.5, dpi = DPI)

  ## Fig 2B: aging-signature TF motif heatmap (donor-averaged z)
  cts <- sort(unique(meta$celltype))
  hm <- do.call(rbind, lapply(cts, function(ct) {
    cells <- which(meta$celltype == ct)
    if (length(cells) < 50) return(NULL)
    data.frame(celltype = ct, motif = c(up, dn),
               z = rowMeans(z[c(up, dn), cells, drop = FALSE]))
  }))
  top <- c(up, dn)
  p3 <- ggplot(hm, aes(x = celltype, y = factor(motif, levels = rev(top)), fill = z)) +
    geom_tile() + scale_fill_gradient2(low = "#2C6FA6", mid = "white", high = "#B22222") +
    labs(x = NULL, y = NULL, title = "Aging signature TF motif activity across cell types") +
    theme_bw(base_size = 10) + theme(axis.text.x = element_text(angle = 45, hjust = 1),
                                     axis.text.y = element_text(size = 5.2))
  ggsave(file.path(pan.dir, "aging_TF_heatmap_600.png"), p3, width = 8.5, height = 8, dpi = DPI)
  rm(z, meta); gc()
  cat("MOTIF STAGE DONE\n")
}

## ============================================================
## Stage 3: all CSV-driven panels
## ============================================================
if (stage == "csv") {
  PG("data.table"); PG("ggplot2"); PG("ggrepel"); PG("scales")
  setDTthreads(4)
  COL_UP <- "#C0392B"; COL_DOWN <- "#2E86C1"
  ct.order <- c("MuSC","Myonuclei","SMC/Pericyte","EC","Macrophage","T/NK cells","FAP")
  Rf <- function(...) fread(file.path(res, ...))

  ## ---- Fig 2A: DA peak counts (script 19A) ----
  da <- Rf("DA_peaks_Aged_vs_Young_by_celltype.csv")
  da.sig <- da[!is.na(padj) & padj < 0.05 & abs(log2FoldChange) > 1]
  da.sig[, dir := ifelse(log2FoldChange > 0, "Up in Aged", "Down in Aged")]
  da.sum <- da.sig[, .N, by = .(celltype, dir)]
  pDA <- ggplot(da.sum, aes(factor(celltype, levels = ct.order), N, fill = dir)) +
    geom_col(width = 0.72) +
    scale_fill_manual(values = c("Up in Aged" = COL_UP, "Down in Aged" = COL_DOWN), name = NULL) +
  labs(x = NULL, y = "Differentially accessible peaks\n(padj<0.05, |LFC|>1)",
       title = "Chromatin accessibility changes with age") +
    theme_classic(11) + theme(axis.text.x = element_text(angle = 30, hjust = 1),
                              legend.position = "top", plot.title = element_text(face = "bold", size = 11))
  ggsave(file.path(pan.dir, "p_DA_counts_600.png"), pDA, width = 6.2, height = 3.4, dpi = DPI)

  ## ---- Fig 3A: DEG counts (script 19B) ----
  ds <- Rf("deg", "DEG_summary_by_celltype.csv")
  ds.long <- melt(ds[, .(celltype, `n_up_Aged`, `n_down_Aged`)],
                  id.vars = "celltype", variable.name = "dir", value.name = "N")
  ds.long[, dir := ifelse(dir == "n_up_Aged", "Up in Aged", "Down in Aged")]
  pDEG <- ggplot(ds.long, aes(factor(celltype, levels = ct.order[ct.order %in% celltype]), N, fill = dir)) +
    geom_col(width = 0.72) +
    geom_text(aes(label = N), position = position_stack(vjust = 0.5), size = 2.9, colour = "white") +
    scale_fill_manual(values = c("Up in Aged" = COL_UP, "Down in Aged" = COL_DOWN), name = NULL) +
    labs(x = NULL, y = "Differentially expressed genes (padj<0.05)",
         title = "Transcriptional changes with age (pseudobulk, 5v5)") +
    theme_classic(11) + theme(axis.text.x = element_text(angle = 30, hjust = 1),
                              legend.position = "top", plot.title = element_text(face = "bold", size = 11))
  ggsave(file.path(pan.dir, "p_DEG_counts_600.png"), pDEG, width = 6.2, height = 3.4, dpi = DPI)

  ## ---- Fig 4D: regulatory chain (script 19C) ----
  rc <- Rf("regulatory_chain_FAP.csv")
  rc.sig <- rc[p_binom < 0.05 & DA_direction == "up"][order(median_rna_lfc)]
  rc.sig[, TF := factor(TF, levels = TF)]
  pCH <- ggplot(rc.sig, aes(TF, frac_rna_up)) +
    geom_col(aes(fill = frac_rna_up), width = 0.72) +
    geom_hline(yintercept = 0.5, linetype = 2, colour = "grey40") +
    scale_fill_gradient2(low = "#2E86C1", mid = "#F4D03F", high = "#C0392B",
                         midpoint = 0.65, name = NULL) +
    coord_flip() +
    labs(x = NULL, y = "Fraction of linked target genes up-regulated in Aged FAP",
         title = "TF -> open enhancer -> target-gene chain (FAP)") +
    theme_classic(11) + theme(legend.position = "none",
                              plot.title = element_text(face = "bold", size = 11))
  ggsave(file.path(pan.dir, "p_regchain_600.png"), pCH, width = 4.6, height = 3.4, dpi = DPI)

  ## ---- Fig 5C: model schematic (script 19D) ----
  sch <- data.frame(
    x = c(0.02, 0.36, 0.36, 0.70, 0.36, 0.70, 0.02),
    y = c(0.60, 0.86, 0.34, 0.86, 0.60, 0.34, 0.08),
    w = c(0.30, 0.28, 0.28, 0.28, 0.28, 0.28, 0.96),
    h = c(0.16, 0.14, 0.14, 0.14, 0.14, 0.14, 0.13),
    lab = c("Ageing + disuse\n(motor-neuron decline)",
            "Clock repressor-arm\ncollapse  NR1D1/PER/DBP\n(FAP, padj 3e-9)",
            "Loss of clock licence\n-> AP-1 / CEBP / NFIA\ntake over enhancers",
            "SASP / fibrosis enhancers\nopen  (1,122 DA peaks)",
            "Chromatin primed\n>90% not transcribed\n(OR 1.85-2.55)",
            "FAP senescence &\nfibrotic niche\n(IL6, LUM, POSTN...)",
            "Myonuclear niche: glycolytic / sarcomeric decline ; MUSK up-regulation"),
    fill = c("#E8E8E8", "#D6EAF8", "#FCF3CF", "#FADBD8", "#FADBD8", "#F5B7B1", "#E8F8F5"),
    stringsAsFactors = FALSE)
  pMO <- ggplot(sch, aes(x, y)) +
    geom_rect(aes(xmin = x, xmax = x + w, ymin = y - h/2, ymax = y + h/2, fill = fill),
              colour = "grey30", linewidth = 0.35) +
    geom_text(aes(label = lab), size = 2.55, lineheight = 0.85) +
    annotate("segment", x = 0.32, xend = 0.36, y = 0.66, yend = 0.84, linewidth = 0.5,
             arrow = arrow(length = unit(2.2, "mm"))) +
    annotate("segment", x = 0.50, xend = 0.50, y = 0.79, yend = 0.70, linewidth = 0.5,
             arrow = arrow(length = unit(2.2, "mm"))) +
    annotate("segment", x = 0.64, xend = 0.70, y = 0.84, yend = 0.66, linewidth = 0.5,
             arrow = arrow(length = unit(2.2, "mm"))) +
    annotate("segment", x = 0.64, xend = 0.70, y = 0.36, yend = 0.54, linewidth = 0.5,
             arrow = arrow(length = unit(2.2, "mm"))) +
    annotate("segment", x = 0.50, xend = 0.50, y = 0.27, yend = 0.36, linewidth = 0.5,
             arrow = arrow(length = unit(2.2, "mm"))) +
    scale_fill_identity() + xlim(0, 1.05) + ylim(0, 1) +
    labs(title = "Working model: clock-licensed senescence enhancers drive the aged FAP niche") +
    theme_void(11) + theme(plot.title = element_text(face = "bold", size = 11, hjust = 0.5))
  ggsave(file.path(pan.dir, "p_model_600.png"), pMO, width = 7.4, height = 3.2, dpi = DPI)

  ## ---- Fig 3B: pathway heatmap (script 13, plot only) ----
  summ <- Rf("pathway_summary_by_celltype.csv")
  mat <- dcast(summ, gene_set ~ celltype, value.var = "mean_lfc")
  m2 <- as.matrix(mat[, -1]); rownames(m2) <- mat$gene_set
  ord <- hclust(dist(m2))$order
  long <- data.table(gene_set = rep(rownames(m2), times = ncol(m2)),
                     celltype = rep(colnames(m2), each = nrow(m2)),
                     mean_lfc = as.vector(m2))
  long[, gene_set := factor(gene_set, levels = rownames(m2)[ord])]
  sigp <- summ[, .(gene_set, celltype, padj_dir)]
  long <- merge(long, sigp, by = c("gene_set", "celltype"), all.x = TRUE)
  long[, star := ifelse(!is.na(padj_dir) & padj_dir < 0.05, "*", "")]
  pPW <- ggplot(long, aes(celltype, gene_set, fill = mean_lfc)) +
    geom_tile(colour = "white", linewidth = 0.4) +
    geom_text(aes(label = star), size = 4, vjust = 0.78, colour = "black") +
    scale_fill_gradient2(low = "#1F618D", mid = "white", high = "#C0392B",
                         midpoint = 0, name = "mean\nlog2FC") +
    labs(x = NULL, y = NULL, title = "Pathway-level shifts in aged muscle (donor pseudobulk)",
         subtitle = "* BH-adjusted directional hypergeometric p < 0.05") +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1), panel.grid = element_blank())
  ggsave(file.path(pan.dir, "pathway_heatmap_600.png"), pPW, width = 8.2, height = 6.4, dpi = DPI)

  ## ---- Fig 4A/4B: coupling FAP scatter + class fractions (script 12, plot only) ----
  cpl <- Rf("coupling_gene_level_by_celltype.csv")
  summC <- Rf("coupling_summary_by_celltype.csv")
  pal <- c(Coupled_up = "#C0392B", Coupled_down = "#1F618D",
           Chromatin_primed = "#E67E22", Chromatin_poised_down = "#7FB3D5",
           RNA_only_up = "#F5B7B1", RNA_only_down = "#AED6F1", Unchanged = "grey85")
  m <- cpl[celltype == "FAP"]
  m[, label := ifelse(class %in% c("Coupled_up", "Coupled_down", "Chromatin_primed") &
                        max_abs_lfc >= 1.5, gene, "")]
  pCP <- ggplot(m, aes(chrom_score, rna_lfc)) +
    geom_point(aes(colour = class), size = 1.6, alpha = 0.75) +
    geom_hline(yintercept = 0, linetype = 2, colour = "grey40") +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey40") +
    scale_colour_manual(values = pal, name = NULL) +
    labs(x = "chromatin accessibility change (mean log2FC of linked DA peaks)",
         y = "RNA expression change (donor pseudobulk log2FC)",
         title = "FAP: chromatin vs transcriptome in ageing",
         subtitle = sprintf("Spearman rho = %.3f | Fisher OR = %.2f (p = %.3g)",
                            summC[celltype == "FAP", spearman_rho],
                            summC[celltype == "FAP", fisher_OR],
                            summC[celltype == "FAP", fisher_p])) +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom", legend.text = element_text(size = 8))
  if (any(m$label != "")) pCP <- pCP + geom_text_repel(aes(label = label), size = 2.6, max.overlaps = 25)
  ggsave(file.path(pan.dir, "coupling_FAP_600.png"), pCP, width = 7.2, height = 6, dpi = DPI)

  bar <- cpl[, .N, by = .(celltype, class)]
  bar[, frac := N / sum(N), by = celltype]
  p2 <- ggplot(bar[class != "Unchanged"], aes(celltype, frac, fill = class)) +
    geom_col(width = 0.72) +
    scale_fill_manual(values = pal, name = NULL) +
    scale_y_continuous(labels = percent) +
    labs(x = NULL, y = "fraction of genes with a linked DA peak",
         title = "Chromatin-transcriptome coupling classes") +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "right")
  ggsave(file.path(pan.dir, "coupling_class_fraction_600.png"), p2, width = 7.5, height = 4.4, dpi = DPI)

  ## ---- Fig 4C: clock TF decoupling (script 14, plot only) ----
  res14 <- Rf("regulatory_chain_FAP.csv")
  clock_tfs <- c("NR1D1","NR1D2","NR2F1","NR2F2","DBP","TEF","HLF","NFIL3")
  fg <- res14[TF %in% clock_tfs]
  mo <- Rf("motif_enrichment_DA_peaks.csv")
  fg2 <- merge(fg[DA_direction == "up", .(TF, n_DA_peaks_with_motif, median_rna_lfc)],
               mo[celltype == "FAP" & direction == "up" & tf_name %in% clock_tfs,
                  .(TF = tf_name, odds_ratio, fdr)], by = "TF")
  rna <- Rf("deg", "DEG_FAP.csv")
  tf_lfc <- rna[gene %in% clock_tfs, .(TF = gene, tf_rna_lfc = log2FoldChange, tf_padj = padj)]
  fg2 <- merge(fg2, tf_lfc, by = "TF", all.x = TRUE)
  pCD <- ggplot(fg2, aes(odds_ratio, tf_rna_lfc)) +
    geom_point(aes(size = n_DA_peaks_with_motif), colour = "#C0392B", alpha = 0.85) +
    geom_hline(yintercept = 0, linetype = 2, colour = "grey40") +
    geom_vline(xintercept = 1, linetype = 2, colour = "grey40") +
    geom_text_repel(aes(label = TF), size = 3.2, max.overlaps = 20) +
    labs(x = "motif enrichment in FAP ageing-open enhancers (odds ratio)",
         y = "TF mRNA change in aged FAP (log2FC)",
         title = "Clock-factor motifs open while their transcripts fall",
         subtitle = "Aged vs Young FAP, donor pseudobulk (5v5)") +
    theme_bw(base_size = 11)
  ggsave(file.path(pan.dir, "clock_TF_decoupling_FAP_600.png"), pCD, width = 7, height = 5.2, dpi = DPI)

  ## ---- Fig 1C: proportion stats (script 17, plot only) ----
  meta <- Rf("cell_metadata.csv")
  meta[, sample := as.character(sample)]
  tab <- meta[, .N, by = .(sample, celltype)]
  tot <- meta[, .N, by = sample]; setnames(tot, "N", "n_total")
  tab <- merge(tab, tot, by = "sample")
  tab[, prop := N / n_total]
  tab[, group := ifelse(grepl("^Aged", sample), "Aged", "Young")]
  dsum <- tab[, .(mean = mean(prop), sd = sd(prop)), by = .(celltype, group)]
  r17 <- Rf("celltype_proportion_stats.csv")
  lv <- r17[order(-log2FC), celltype]
  dsum[, celltype := factor(celltype, levels = lv)]
  pPS <- ggplot(dsum, aes(celltype, mean * 100, fill = group)) +
    geom_col(position = position_dodge(0.75), width = 0.7) +
    geom_errorbar(aes(ymin = (mean - sd) * 100, ymax = (mean + sd) * 100),
                  position = position_dodge(0.75), width = 0.25) +
    scale_fill_manual(values = c(Young = "#2E86C1", Aged = "#C0392B"), name = NULL) +
    labs(x = NULL, y = "proportion (%)",
         title = "Cell-type composition shift (donor mean \u00b1 SD)",
         subtitle = "Wilcoxon 5v5: FAP p = 0.008, EC p = 0.008") +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "top")
  ggsave(file.path(pan.dir, "celltype_proportion_stats_600.png"), pPS, width = 7.4, height = 4.6, dpi = DPI)

  ## ---- Fig 3C: FAP clock per donor (script 18b) ----
  d <- Rf("FAP_clock_per_donor.csv")
  d <- d[!is.na(donor)]
  d[, gene := factor(gene, levels = c("NR1D1","NR1D2","PER1","PER3","DBP","TEF","HLF","CRY2",
                                      "ARNTL","NPAS2","IL6","CCL2","LUM","COL1A1","MUSK","PDGFRA"))]
  d[, lab := ifelse(grepl("^Young", donor), "Young (20-27 y)", "Aged (60-77 y)")]
  pFP <- ggplot(d, aes(lab, log2(CPM10k + 0.05), colour = lab)) +
    geom_point(position = position_jitter(width = 0.12), size = 2, alpha = 0.8) +
    stat_summary(fun = mean, geom = "crossbar", width = 0.5, colour = "black") +
    facet_wrap(~ gene, scales = "free_y", ncol = 8) +
    scale_colour_manual(values = c(`Young (20-27 y)` = "#2E86C1", `Aged (60-77 y)` = "#C0392B"), name = NULL) +
    labs(x = NULL, y = "log2 expression (CPM/10k, FAP pseudobulk)",
         title = "Per-donor expression in FAP: clock repression is donor-consistent") +
    theme_bw(base_size = 9) +
    theme(legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1, size = 7))
  ggsave(file.path(pan.dir, "FAP_clock_per_donor_600.png"), pFP, width = 10.5, height = 4.4, dpi = DPI)

  ## ---- Fig 5A: GTEx module concordance (script 18c, plot only) ----
  m2c <- Rf("gtex_snRNA_concordance.csv")
  mods <- list(
    Clock_repressive_arm = c("NR1D1","NR1D2","PER1","PER2","PER3","CRY1","CRY2","DBP","TEF","HLF","NR2F1","NR2F2","BHLHE40","BHLHE41","NFIL3","RORA","CIART"),
    SASP_inflammatory = c("IL6","CXCL2","CXCL8","CCL2","CSF1","TGFB1","SERPINE1","IGFBP7","MMP2","TIMP1","CDKN1A","CDKN2A","GLB1","LMNB1"),
    ECM_fibrosis = c("COL1A1","COL3A1","DCN","LUM","POSTN","PDGFRA","PRG4","ADAMTS16","ITGB8","COL19A1","SAA2","ABCA8","MMP3","WISP2","CFD"),
    Muscle_fiber = c("PFKM","PGK1","LMOD2","STAC3","ATP5PD","CRYAB","MUSK"))
  mod_of <- unlist(lapply(names(mods), function(k) setNames(rep(k, length(mods[[k]])), mods[[k]])))
  m2c[, module := ifelse(gene %in% names(mod_of), mod_of[gene], "other")]
  pal2 <- c(Clock_repressive_arm = "#8E44AD", SASP_inflammatory = "#E67E22",
            ECM_fibrosis = "#27AE60", Muscle_fiber = "#2E86C1", other = "grey70")
  pGT <- ggplot(m2c, aes(sn_lfc, lfc_aged_vs_young, colour = module)) +
    geom_point(size = 2.7, alpha = 0.9) +
    geom_hline(yintercept = 0, linetype = 2, colour = "grey50") +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
    geom_text_repel(data = m2c[module != "other"], aes(label = gene), size = 2.7,
                    max.overlaps = 20, show.legend = FALSE) +
    scale_colour_manual(values = pal2, name = NULL) +
    labs(x = "log2FC Aged vs Young \u2014 FAP, snMultiome (5v5 donors)",
         y = "log2FC Aged vs Young \u2014 GTEx skeletal muscle (bulk, n=818)",
         title = "External check in GTEx: module-wise concordance",
         subtitle = "FAP SASP/ECM programs concordant; clock arm discordant (bulk = mixed cell types + time-of-day sampling)") +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom")
  ggsave(file.path(pan.dir, "gtex_concordance_by_module_600.png"), pGT, width = 7.8, height = 6, dpi = DPI)

  cat("CSV STAGE DONE\n")
}

## ============================================================
## Stage 4: clock locus tracks (script 16, plot only)
## ============================================================
if (stage == "locus") {
  PG("Matrix"); PG("data.table"); PG("ggplot2"); PG("GenomicRanges")
  PG("jsonlite"); PG("patchwork")
  peaks <- readRDS(file.path(res, "consensus_peaks.rds"))
  pk <- data.table(chr = as.character(seqnames(peaks)), start = start(peaks), end = end(peaks))
  pk[, idx := .I]
  samples <- c(paste0("Young_", 1:5), paste0("Aged_", 1:5))
  get_fap <- function(s) {
    pb <- readRDS(file.path(res, "atac_tmp", paste0("pb_", s, ".rds")))$pb
    v <- pb[, "FAP"]; as.numeric(v) / sum(v) * 1e6
  }
  cpm <- sapply(samples, get_fap); colnames(cpm) <- samples
  grp <- ifelse(grepl("^Young", samples), "Young", "Aged")
  fap_young <- rowMeans(cpm[, grp == "Young", drop = FALSE])
  fap_aged  <- rowMeans(cpm[, grp == "Aged",  drop = FALSE])
  loci <- data.table(
    gene = c("NR1D1", "PER1", "DBP", "HLF"),
    chr  = c("chr17", "chr17", "chr19", "chr17"),
    start = c(40065000, 8110000, 48615000, 55240000),
    end   = c(40135000, 8195000, 48655000, 55395000))
  gene_model <- function(chr, genes) {
    j <- fromJSON(file.path(proj, "data/annotation", paste0("refGene_", chr, ".json")),
                  simplifyDataFrame = TRUE)
    rg <- as.data.table(j$refGene)[name2 %in% genes]
    rbindlist(lapply(seq_len(nrow(rg)), function(i) {
      s <- as.integer(strsplit(gsub(",$", "", rg$exonStarts[i]), ",")[[1]])
      e <- as.integer(strsplit(gsub(",$", "", rg$exonEnds[i]), ",")[[1]])
      data.table(symbol = rg$name2[i], txStart = rg$txStart[i], txEnd = rg$txEnd[i],
                 start = s, end = e)
    }))
  }
  plots <- list()
  for (i in seq_len(nrow(loci))) {
    L <- loci[i]
    sel <- pk[chr == L$chr & end > L$start & start < L$end]
    if (!nrow(sel)) next
    d <- data.table(mid = (sel$start + sel$end) / 2, w = sel$end - sel$start,
                    Young = fap_young[sel$idx], Aged = fap_aged[sel$idx],
                    start = sel$start, end = sel$end)
    win_counts <- colSums(cpm[sel$idx, , drop = FALSE])
    y <- log2(win_counts[grp == "Young"] + 1); a <- log2(win_counts[grp == "Aged"] + 1)
    p <- tryCatch(wilcox.test(a, y)$p.value, error = function(e) NA)
    long <- rbind(data.table(mid = d$mid, cpm = d$Young, grp = "Young"),
                  data.table(mid = d$mid, cpm = d$Aged, grp = "Aged"))
    p1 <- ggplot(long, aes(mid, cpm, fill = grp)) +
      geom_area(data = subset(long, grp == "Aged"), fill = "#C0392B", alpha = 0.45, colour = NA) +
      geom_line(data = subset(long, grp == "Young"), aes(colour = grp), linewidth = 0.6) +
      scale_colour_manual(values = c(Young = "#2E86C1", Aged = "#C0392B"), name = NULL) +
      scale_x_continuous(labels = function(x) sprintf("%.2f Mb", x / 1e6), expand = c(0.005, 0)) +
      labs(y = "accessibility (CPM)", x = NULL,
           title = sprintf("%s : FAP chromatin accessibility", L$gene),
           subtitle = sprintf("locus log2FC = %.2f, Wilcoxon P = %.2g",
                              log2((mean(d$Aged) + .1) / (mean(d$Young) + .1)), p)) +
      theme_bw(base_size = 9) +
      theme(panel.grid.minor = element_blank(), legend.position = "top",
            legend.key.size = unit(0.4, "cm"), plot.title = element_text(size = 9, face = "bold"))
    gm <- gene_model(L$chr, L$gene)
    p2 <- ggplot(gm) +
      geom_segment(aes(x = min(txStart), xend = max(txEnd), y = 0, yend = 0), linewidth = 0.4) +
      geom_rect(aes(xmin = start, xmax = end, ymin = -0.2, ymax = 0.2), fill = "#34495E", colour = NA) +
      scale_x_continuous(labels = function(x) sprintf("%.2f Mb", x / 1e6),
                         limits = c(L$start, L$end), expand = c(0.005, 0)) +
      scale_y_continuous(limits = c(-0.6, 0.6)) +
      labs(x = NULL, y = NULL) +
      theme_bw(base_size = 9) +
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
            panel.grid = element_blank(), axis.line.y = element_blank())
    plots[[i]] <- (p1 / p2) + plot_layout(heights = c(3, 1))
  }
  p <- wrap_plots(plots, ncol = 2)
  ggsave(file.path(pan.dir, "clock_loci_tracks_600.png"), p, width = 9.6, height = 6.4, dpi = DPI)
  cat("LOCUS STAGE DONE\n")
}

## ============================================================
## Stage 5: compose the 5 main figures -> TIFF(LZW) + PNG at 600 dpi
## ============================================================
if (stage == "compose") {
  PG("png"); PG("grid"); PG("gridExtra")
  F  <- pan.dir          # 600-dpi panels
  W  <- round(3000 * DPI / 210)   # 8571 px, same physical size as the 210-dpi originals
  compose <- function(cells, files, out, widths, heights, labels) {
    for (dev in c("tiff", "png")) {
      if (dev == "tiff") {
        tiff(out, width = W, height = round(W * sum(heights) / sum(widths)),
             res = DPI, compression = "lzw")
      } else {
        png(sub("\\.tiff$", ".png", out), width = W,
            height = round(W * sum(heights) / sum(widths)), res = DPI)
      }
      lay <- grid.layout(nrow = length(heights), ncol = length(widths),
                         widths = unit(widths, "null"), heights = unit(heights, "null"))
      pushViewport(viewport(layout = lay))
      for (i in seq_along(files)) {
        if (is.na(files[i]) || !file.exists(files[i])) { cat("WARN missing:", files[i], "\n"); next }
        r <- cells[[i]][1]; cl <- cells[[i]][2]
        pushViewport(viewport(layout.pos.row = r, layout.pos.col = cl))
        img <- readPNG(files[i])
        imgA <- nrow(img) / ncol(img)
        cellA <- heights[r] / widths[cl]
        if (imgA > cellA) {
          grid.raster(img, x = 0.5, y = 0.5, height = 1, width = unit(cellA / imgA, "npc"))
        } else {
          grid.raster(img, x = 0.5, y = 0.5, width = 1, height = unit(imgA / cellA, "npc"))
        }
        grid.text(labels[i], x = 0.015, y = 0.97, just = c("left","top"),
                  gp = gpar(fontsize = 20, fontface = "bold"))
        popViewport()
      }
      dev.off()
      cat("wrote", out, if (dev == "png") "(png copy)\n" else "\n")
    }
  }
  compose(cells = list(c(1,1), c(1,2), c(2,1)),
          files = c(file.path(F, "umap_celltype_600.png"),
                    file.path(F, "umap_group_600.png"),
                    file.path(F, "celltype_proportion_stats_600.png")),
          out = file.path(out.dir, "Figure1.tiff"),
          widths = c(1, 1, 2), heights = c(1, 1.05), labels = c("A","B","C"))
  compose(cells = list(c(1,1), c(1,2), c(2,1)),
          files = c(file.path(F, "p_DA_counts_600.png"),
                    file.path(F, "aging_TF_heatmap_600.png"),
                    file.path(F, "chromatin_aging_score_by_celltype_600.png")),
          out = file.path(out.dir, "Figure2.tiff"),
          widths = c(1, 1, 2), heights = c(1, 1.05), labels = c("A","B","C"))
  compose(cells = list(c(1,1), c(2,1), c(2,2)),
          files = c(file.path(F, "p_DEG_counts_600.png"),
                    file.path(F, "pathway_heatmap_600.png"),
                    file.path(F, "FAP_clock_per_donor_600.png")),
          out = file.path(out.dir, "Figure3.tiff"),
          widths = c(1.15, 1), heights = c(1, 1.25), labels = c("A","B","C"))
  compose(cells = list(c(1,1), c(1,2), c(2,1), c(2,2)),
          files = c(file.path(F, "coupling_FAP_600.png"),
                    file.path(F, "coupling_class_fraction_600.png"),
                    file.path(F, "clock_TF_decoupling_FAP_600.png"),
                    file.path(F, "p_regchain_600.png")),
          out = file.path(out.dir, "Figure4.tiff"),
          widths = c(1.25, 1), heights = c(1.1, 1), labels = c("A","B","C","D"))
  compose(cells = list(c(1,1), c(2,1), c(3,1)),
          files = c(file.path(F, "gtex_concordance_by_module_600.png"),
                    file.path(F, "clock_loci_tracks_600.png"),
                    file.path(F, "p_model_600.png")),
          out = file.path(out.dir, "Figure5.tiff"),
          widths = c(1), heights = c(1, 1, 0.9), labels = c("A","B","C"))
  cat("COMPOSE STAGE DONE\n")
}
