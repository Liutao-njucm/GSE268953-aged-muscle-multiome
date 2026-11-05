# ============================================================
# 16: locus accessibility tracks for clock genes
#     Young vs Aged FAP pseudobulk ATAC (CPM), + refGene models
# ============================================================
.libPaths(c("C:/Users/wang/R/win-library/4.5", .libPaths()))
suppressMessages({library(Matrix); library(data.table); library(ggplot2); library(GenomicRanges)
  library(jsonlite); library(patchwork)})

proj <- "D:/workbuddy/outputs/2026-10-25-20-23-31"
res <- file.path(proj, "results/GSE268953")
ann.dir <- file.path(proj, "data/annotation")
fig.dir <- file.path(res, "figures")

peaks <- readRDS(file.path(res, "consensus_peaks.rds"))
pk <- data.table(chr = as.character(seqnames(peaks)), start = start(peaks), end = end(peaks))
pk[, idx := .I]

samples <- c(paste0("Young_", 1:5), paste0("Aged_", 1:5))
get_fap <- function(s) {
  pb <- readRDS(file.path(res, "atac_tmp", paste0("pb_", s, ".rds")))$pb
  v <- pb[, "FAP"]
  tot <- sum(v)
  as.numeric(v) / tot * 1e6
}
cpm <- sapply(samples, get_fap)   # peaks x samples
colnames(cpm) <- samples
grp <- ifelse(grepl("^Young", samples), "Young", "Aged")
fap_young <- rowMeans(cpm[, grp == "Young", drop = FALSE])
fap_aged  <- rowMeans(cpm[, grp == "Aged",  drop = FALSE])
cat("tracks built:", length(fap_young), "peaks\n")

loci <- data.table(
  gene = c("NR1D1", "PER1", "DBP", "HLF"),
  chr  = c("chr17", "chr17", "chr19", "chr17"),
  start = c(40065000, 8110000, 48615000, 55240000),
  end   = c(40135000, 8195000, 48655000, 55395000)
)

gene_model <- function(chr, genes) {
  j <- fromJSON(file.path(ann.dir, paste0("refGene_", chr, ".json")), simplifyDataFrame = TRUE)
  rg <- as.data.table(j$refGene)[name2 %in% genes]
  rbindlist(lapply(seq_len(nrow(rg)), function(i) {
    s <- as.integer(strsplit(gsub(",$", "", rg$exonStarts[i]), ",")[[1]])
    e <- as.integer(strsplit(gsub(",$", "", rg$exonEnds[i]), ",")[[1]])
    data.table(symbol = rg$name2[i], txStart = rg$txStart[i], txEnd = rg$txEnd[i],
               start = s, end = e)
  }))
}

plots <- list(); stats_rows <- list()
for (i in seq_len(nrow(loci))) {
  L <- loci[i]
  sel <- pk[chr == L$chr & end > L$start & start < L$end]
  if (!nrow(sel)) next
  d <- data.table(mid = (sel$start + sel$end) / 2, w = sel$end - sel$start,
                  Young = fap_young[sel$idx], Aged = fap_aged[sel$idx],
                  start = sel$start, end = sel$end)
  d[, lfc := log2((Aged + 0.1) / (Young + 0.1))]
  ## locus-level test across donors
  win_counts <- colSums(cpm[sel$idx, , drop = FALSE])
  y <- log2(win_counts[grp == "Young"] + 1); a <- log2(win_counts[grp == "Aged"] + 1)
  p <- tryCatch(wilcox.test(a, y)$p.value, error = function(e) NA)
  stats_rows[[L$gene]] <- data.table(gene = L$gene, n_peaks = nrow(d),
                                     mean_cpm_young = mean(d$Young), mean_cpm_aged = mean(d$Aged),
                                     log2FC_locus = log2((mean(d$Aged) + .1) / (mean(d$Young) + .1)),
                                     p_wilcox = p)
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
ggsave(file.path(fig.dir, "clock_loci_tracks.png"), p, width = 9.6, height = 6.4, dpi = 300)
st <- rbindlist(stats_rows)
fwrite(st, file.path(res, "clock_locus_accessibility.csv"))
cat("\n=== locus accessibility (FAP, Young vs Aged) ===\n"); print(as.data.frame(st), digits = 3)
cat("saved clock_loci_tracks.png\n")
