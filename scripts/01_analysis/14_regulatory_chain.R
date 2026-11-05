# ============================================================
# 14: TF -> DA peak -> target gene regulatory chain (FAP focus)
#   which TFs' motifs sit in ageing-open/closed enhancers,
#   and what happens to the linked genes' expression
# ============================================================
.libPaths(c("C:/Users/wang/R/win-library/4.5", .libPaths()))
suppressMessages({library(Matrix); library(data.table); library(ggplot2)})

out.dir <- "D:/workbuddy/outputs/2026-10-25-20-23-31/results/GSE268953"
mm <- readRDS(file.path(out.dir, "consensus_motif_matches.rds"))
map <- fread(file.path(out.dir, "motif_id_to_tf_name.csv"))
ann <- fread(file.path(out.dir, "peak_annotation.csv"))
da  <- fread(file.path(out.dir, "DA_peaks_Aged_vs_Young_by_celltype.csv"))
rna <- fread(file.path(out.dir, "deg/DEG_FAP.csv"))

cat("motif matrix:", dim(mm), "| rownames:", !is.null(rownames(mm)), "\n")

## motif matrix covers all consensus peaks; match DA peaks by rowname (peak_id)
stopifnot(!is.null(rownames(mm)))
peak_id_sub <- rownames(mm)

## FAP significant DA peaks
fap <- da[celltype == "FAP" & !is.na(padj) & padj < 0.05 & abs(log2FoldChange) > 1]
fap[, dir := ifelse(log2FoldChange > 0, "up", "down")]
cat("FAP sig DA peaks:", nrow(fap), " (up", sum(fap$dir == "up"), "/ down", sum(fap$dir == "down"), ")\n")

ji <- match(fap$peak_id, peak_id_sub)
stopifnot(!any(is.na(ji)))
mot <- as.matrix(mm[ji, ])          # peaks x motifs logical

ann_map <- setNames(ann$linked_gene, ann$peak_id)

colnames(mot) <- colnames(mm)
tf_of <- setNames(map$tf_name, map$motif_id)
colnames(mot) <- ifelse(colnames(mot) %in% names(tf_of), tf_of[colnames(mot)], colnames(mot))

rna_map <- setNames(rna$log2FoldChange, rna$gene)
padj_map <- setNames(rna$padj, rna$gene)

tfs <- c("NR1D1","NR1D2","NR2F1","NR2F2","DBP","TEF","HLF","NFIL3",   # clock
         "FOS","JUN","JUNB","JUND","FOSL1","FOSL2","BATF","MAFK","MAFB","MAF",  # AP-1
         "CEBPA","CEBPD","CEBPE","CEBPB","RUNX2","RUNX3","RUNX1","NFIA","TEAD1","MEF2C")

rows <- list()
for (tf in intersect(tfs, colnames(mot))) {
  for (d in c("up", "down")) {
    sel <- fap$dir == d & mot[, tf]
    n <- sum(sel)
    if (n < 5) next
    genes <- unique(na.omit(ann_map[fap$peak_id[sel]]))
    lfc <- rna_map[genes]; lfc <- lfc[!is.na(lfc)]
    rows[[length(rows) + 1]] <- data.table(
      TF = tf, DA_direction = d, n_DA_peaks_with_motif = n,
      n_linked_genes = length(genes),
      n_genes_measured = length(lfc),
      median_rna_lfc = if (length(lfc)) median(lfc) else NA,
      frac_rna_up = if (length(lfc)) mean(lfc > 0) else NA,
      top_genes = paste(head(genes[order(-abs(rna_map[genes]))], 12), collapse = ", ")
    )
  }
}
res <- rbindlist(rows)
res[, p_binom := mapply(function(k, n) binom.test(k, n, 0.5)$p.value, round(frac_rna_up * n_genes_measured), n_genes_measured)]
fwrite(res, file.path(out.dir, "regulatory_chain_FAP.csv"))
cat("\n=== regulatory chain (FAP) ===\n")
print(res[order(DA_direction, -n_DA_peaks_with_motif)], digits = 2)

## ---- figure: clock TFs - motif opens, mRNA falls ----
clock_tfs <- c("NR1D1","NR1D2","NR2F1","NR2F2","DBP","TEF","HLF","NFIL3")
fg <- res[TF %in% clock_tfs]
motif_or <- fread(file.path(out.dir, "motif_enrichment_DA_peaks.csv"))
fg2 <- merge(fg[DA_direction == "up", .(TF, n_DA_peaks_with_motif, median_rna_lfc)],
             motif_or[celltype == "FAP" & direction == "up" & tf_name %in% clock_tfs,
                      .(TF = tf_name, odds_ratio, fdr)], by = "TF")
tf_lfc <- rna[gene %in% clock_tfs, .(TF = gene, tf_rna_lfc = log2FoldChange, tf_padj = padj)]
fg2 <- merge(fg2, tf_lfc, by = "TF", all.x = TRUE)

p <- ggplot(fg2, aes(odds_ratio, tf_rna_lfc)) +
  geom_point(aes(size = n_DA_peaks_with_motif), colour = "#C0392B", alpha = 0.85) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey40") +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey40") +
  ggrepel::geom_text_repel(aes(label = TF), size = 3.2, max.overlaps = 20) +
  labs(x = "motif enrichment in FAP ageing-open enhancers (odds ratio)",
       y = "TF mRNA change in aged FAP (log2FC)",
       title = "Clock-factor motifs open while their transcripts fall",
       subtitle = "Aged vs Young FAP, donor pseudobulk (5v5)") +
  theme_bw(base_size = 11)
ggsave(file.path(out.dir, "figures/clock_TF_decoupling_FAP.png"), p, width = 7, height = 5.2, dpi = 300)
cat("\nsaved regulatory_chain_FAP.csv + clock_TF_decoupling_FAP.png\n")
