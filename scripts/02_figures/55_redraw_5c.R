# ============================================================
# 55 (v2): redraw Fig 5C clock-locus tracks, genome-browser style
#
#  v1 problems fixed here:
#   a) window was gene +/- ~30-48 kb -> gene body occupied only 11-19% of the
#      window, so the gene model looked "short" and detatched from the tracks.
#      Now: window = gene body +/- max(12 kb, 0.30 * gene span)
#   b) y axis had ticks but no axis title, x axis unlabelled -> not acceptable
#      in a journal figure. Now: y = "ATAC-seq (CPM)", x = "chrN position (Mb)"
#   c) only the longest transcript is drawn (HLF had 2 overlapping transcripts
#      painted on top of each other -> the "irregular thick/thin bars" artefact)
#   d) gene model now uses UTR (thin) vs CDS (thick) boxes + 3'-end arrow,
#      i.e. the standard Ensembl/IGV convention
#   e) a single legend is printed once under the whole panel
# ============================================================
.libPaths(c("C:/Users/wang/R/win-library/4.5", .libPaths()))
suppressMessages({library(Matrix); library(data.table); library(ggplot2)
  library(GenomicRanges); library(jsonlite); library(patchwork)})

proj <- "D:/workbuddy/outputs/2026-10-25-20-23-31"
res  <- file.path(proj, "results/GSE268953")
pan.dir <- file.path(res, "panels_600dpi")
DPI <- 600

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

loci <- data.table(gene = c("NR1D1", "PER1", "DBP", "HLF"),
                   chr  = c("chr17", "chr17", "chr19", "chr17"))

gene_model <- function(chr, gene) {
  j <- fromJSON(file.path(proj, "data/annotation", paste0("refGene_", chr, ".json")),
                simplifyDataFrame = TRUE)
  rg <- as.data.table(j$refGene)[name2 == gene]
  rg[, span := txEnd - txStart]
  rg <- rg[order(-span)][1]            # longest transcript only
  s <- as.integer(strsplit(gsub(",$", "", rg$exonStarts), ",")[[1]])
  e <- as.integer(strsplit(gsub(",$", "", rg$exonEnds),   ",")[[1]])
  list(tx = rg$name, txStart = rg$txStart, txEnd = rg$txEnd,
       cdsStart = rg$cdsStart, cdsEnd = rg$cdsEnd, strand = rg$strand,
       exon = data.table(start = s, end = e))
}

COL_Y <- "#2E86C1"; COL_A <- "#C0392B"
plots <- list(); stat_rows <- list()

for (i in seq_len(nrow(loci))) {
  L  <- loci[i]
  gm <- gene_model(L$chr, L$gene)
  span <- gm$txEnd - gm$txStart
  pad  <- max(12000, round(0.30 * span))
  win  <- c(gm$txStart - pad, gm$txEnd + pad)
  strand <- gm$strand
  tss <- if (strand == "+") gm$txStart else gm$txEnd
  dir <- if (strand == "+") 1 else -1

  sel <- pk[chr == L$chr & end > win[1] & start < win[2]]
  sel[, mid := (start + end) / 2]
  sel[, Young := fap_young[idx]]; sel[, Aged := fap_aged[idx]]
  sel[, d_tss := abs(mid - tss)]
  pro <- sel[which.min(d_tss)]
  yv <- cpm[pro$idx, grp == "Young"]; av <- cpm[pro$idx, grp == "Aged"]
  pw <- tryCatch(wilcox.test(av, yv)$p.value, error = function(e) NA)
  pro_lfc <- log2((mean(av) + .1) / (mean(yv) + .1))
  stat_rows[[i]] <- data.table(gene = L$gene, chr = L$chr, strand = strand, tss = tss,
                               win_start = win[1], win_end = win[2],
                               win_kb = round(diff(win) / 1e3),
                               gene_kb = round(span / 1e3, 1),
                               prom_start = pro$start, prom_end = pro$end,
                               prom_d_tss = round(pro$d_tss),
                               prom_young = round(mean(yv), 1), prom_aged = round(mean(av), 1),
                               prom_log2FC = round(pro_lfc, 3), prom_p = signif(pw, 2),
                               n_peaks = nrow(sel), n_prom = sum(sel$d_tss <= 1000))

  ## x axis in kb (comma-separated), browser style; exact TSS printed in subtitle
  brk <- pretty(win, n = 4); brk <- brk[brk >= win[1] & brk <= win[2]]
  xlab_fun <- function(x) formatC(round(x / 1000), format = "d", big.mark = ",")

  ymax <- max(sel$Young, sel$Aged) * 1.20
  track <- function(dat, col, lab, show_y) {
    ggplot(dat) +
      annotate("rect", xmin = tss - 1000, xmax = tss + 1000, ymin = -Inf, ymax = Inf,
               fill = "grey55", alpha = 0.13) +
      geom_rect(aes(xmin = start, xmax = end, ymin = 0, ymax = .data[[lab]]),
                fill = col, alpha = 0.85, colour = NA) +
      geom_vline(xintercept = tss, linetype = "22", linewidth = 0.32, colour = "grey15") +
      annotate("text", x = -Inf, y = Inf, label = lab, hjust = -0.25, vjust = 1.5,
               size = 2.4, colour = col, fontface = "bold") +
      scale_x_continuous(limits = win, breaks = brk, expand = c(0.004, 0)) +
      scale_y_continuous(limits = c(0, ymax), expand = expansion(mult = c(0, 0)),
                         breaks = pretty(c(0, ymax), 3)) +
      labs(y = if (show_y) "ATAC-seq\n(CPM)" else NULL, x = NULL) +
      theme_bw(base_size = 7.5) +
      theme(panel.grid = element_blank(),
            axis.text.x  = element_blank(), axis.ticks.x = element_blank(),
            axis.text.y  = element_text(size = 5.6, colour = "grey25"),
            axis.title.y = element_text(size = 6.2, colour = "grey20"),
            plot.margin = margin(1, 3, 0, 1))
  }

  pA <- track(sel, COL_A, "Aged", TRUE) +
    labs(title = sprintf("%s  (%s strand, %s, %d exons)", L$gene, strand, gm$tx, nrow(gm$exon)),
         subtitle = sprintf(paste0("TSS (hg38): %s:%s   |   promoter peak %d bp: ",
                                   "log2FC = %+.2f, P = %.2g   |   %d peaks / %d kb"),
                            L$chr, formatC(tss, format = "d", big.mark = ","),
                            round(pro$d_tss), pro_lfc, pw, nrow(sel),
                            round(diff(win) / 1e3))) +
    theme(plot.title    = element_text(size = 8.6, face = "bold", hjust = 0),
          plot.subtitle = element_text(size = 6.0, colour = "grey25", hjust = 0),
          plot.margin   = margin(3, 3, 0, 1))
  pY <- track(sel, COL_Y, "Young", FALSE)

  ## ---- gene model: grey intron line, thin UTR boxes, thick CDS boxes ----
  ## UTR clipping must be strictly per-exon: cdsStart/cdsEnd may fall OUTSIDE
  ## the terminal exon (e.g. PER1 exon 1 is pure 68 bp UTR), so clamp both ends
  ex <- copy(gm$exon)
  ex[, `:=`(cs = pmax(start, gm$cdsStart), ce = pmin(end, gm$cdsEnd))]
  thick <- ex[cs < ce, .(xmin = cs, xmax = ce)]
  thinA <- ex[start < pmin(end, gm$cdsStart),
              .(xmin = start, xmax = pmin(end, gm$cdsStart))]
  thinB <- ex[end > pmax(start, gm$cdsEnd),
              .(xmin = pmax(start, gm$cdsEnd), xmax = end)]
  thin  <- rbind(thinA, thinB)

  ## direction arrow placed INSIDE the largest intron (IGV convention) so it
  ## never overlays exon boxes -- the old txEnd-anchored arrow fused the two
  ## leftmost CDS boxes of NR1D1 into one bar
  exs <- sort(ex$start); exe <- sort(ex$end)
  gaps <- data.table(gfrom = exe[-length(exe)], gto = exs[-1])
  gaps[, g := gto - gfrom]
  big <- gaps[which.max(g)]
  gmid <- (big$gfrom + big$gto) / 2
  half <- 0.35 * big$g / 2
  ## arrow points toward the 3' end of the transcript
  a3_from <- gmid - dir * half; a3_to <- gmid + dir * half

  pG <- ggplot() +
    annotate("rect", xmin = tss - 1000, xmax = tss + 1000, ymin = -Inf, ymax = Inf,
             fill = "grey55", alpha = 0.13) +
    annotate("rect", xmin = gm$txStart, xmax = gm$txEnd, ymin = -0.055, ymax = 0.055,
             fill = "grey72", colour = NA) +
    annotate("rect", xmin = thin$xmin,  xmax = thin$xmax,  ymin = -0.13, ymax = 0.13,
             fill = "#5D6D7E", colour = NA) +
    annotate("rect", xmin = thick$xmin, xmax = thick$xmax, ymin = -0.26, ymax = 0.26,
             fill = "#1B2631", colour = NA) +
    annotate("segment", x = a3_from, xend = a3_to, y = 0, yend = 0, linewidth = 0.45,
             colour = "grey10", arrow = arrow(length = unit(1.5, "mm"), type = "closed")) +
    annotate("segment", x = tss, xend = tss + dir * diff(win) * 0.045, y = 0.62, yend = 0.62,
             linewidth = 0.4, colour = "grey10",
             arrow = arrow(length = unit(1.5, "mm"), type = "closed")) +
    annotate("text", x = tss, y = 0.78, label = "TSS", hjust = 0.5, vjust = 0,
             size = 2.5, colour = "grey10", fontface = "bold") +
    ## 5' sits on the UPSTREAM side of the TSS, 3' beyond the transcript end;
    ## both are placed outside the gene body so they never sit "inside" the model
    annotate("text", x = tss - dir * diff(win) * 0.030, y = 0.62, label = "5'",
             size = 2.4, colour = "grey30", hjust = 0.5) +
    annotate("text", x = (if (dir == 1) gm$txEnd else gm$txStart) + dir * diff(win) * 0.030,
             y = 0.62, label = "3'", size = 2.4, colour = "grey30", hjust = 0.5) +
    geom_vline(xintercept = tss, linetype = "22", linewidth = 0.32, colour = "grey15") +
    scale_x_continuous(limits = win, breaks = brk, expand = c(0.004, 0),
                       labels = xlab_fun) +
    scale_y_continuous(limits = c(-0.55, 1.15)) +
    labs(x = sprintf("%s position (kb)", L$chr), y = NULL) +
    theme_bw(base_size = 7.5) +
    theme(panel.grid = element_blank(),
          axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          axis.line.y = element_blank(),
          axis.text.x  = element_text(size = 5.8, colour = "grey25"),
          axis.title.x = element_text(size = 6.2, colour = "grey20"),
          plot.margin = margin(0, 3, 2, 1))

  plots[[i]] <- (pA / pY / pG) + plot_layout(heights = c(1, 1, 0.92))
}

## ---------------- one shared legend, printed once ----------------
lg <- ggplot() + xlim(0, 1) + ylim(0, 1) + theme_void() +
  annotate("rect", xmin = 0.012, xmax = 0.030, ymin = 0.62, ymax = 0.88,
           fill = COL_A, alpha = 0.85) +
  annotate("text", x = 0.038, y = 0.75, hjust = 0, size = 2.5, colour = "grey15",
           label = "Aged FAP (n = 5)") +
  annotate("rect", xmin = 0.012, xmax = 0.030, ymin = 0.16, ymax = 0.42,
           fill = COL_Y, alpha = 0.85) +
  annotate("text", x = 0.038, y = 0.29, hjust = 0, size = 2.5, colour = "grey15",
           label = "Young FAP (n = 5)") +
  annotate("rect", xmin = 0.250, xmax = 0.282, ymin = 0.62, ymax = 0.88,
           fill = "grey45", alpha = 0.85) +
  annotate("text", x = 0.290, y = 0.75, hjust = 0, size = 2.5, colour = "grey15",
           label = "one consensus peak (true width)") +
  annotate("rect", xmin = 0.250, xmax = 0.282, ymin = 0.16, ymax = 0.42,
           fill = "grey55", alpha = 0.13) +
  annotate("segment", x = 0.266, xend = 0.266, y = 0.16, yend = 0.42,
           linetype = "22", linewidth = 0.32, colour = "grey15") +
  annotate("text", x = 0.290, y = 0.29, hjust = 0, size = 2.5, colour = "grey15",
           label = "TSS  and  TSS \u00b1 1 kb promoter band") +
  annotate("rect", xmin = 0.660, xmax = 0.678, ymin = 0.68, ymax = 0.84,
           fill = "#1B2631") +
  annotate("rect", xmin = 0.660, xmax = 0.678, ymin = 0.20, ymax = 0.36,
           fill = "#5D6D7E") +
  annotate("segment", x = 0.640, xend = 0.700, y = 0.28, yend = 0.28,
           linewidth = 0.5, colour = "grey72") +
  annotate("text", x = 0.690, y = 0.78, hjust = 0, size = 2.5, colour = "grey15",
           label = "CDS") +
  annotate("text", x = 0.690, y = 0.28, hjust = 0, size = 2.5, colour = "grey15",
           label = "UTR / intron") +
  annotate("segment", x = 0.640, xend = 0.700, y = 0.58, yend = 0.58,
           linewidth = 0.45, colour = "grey10",
           arrow = arrow(length = unit(1.2, "mm"), type = "closed")) +
  annotate("text", x = 0.710, y = 0.58, hjust = 0, size = 2.5, colour = "grey15",
           label = "5'\u21923' direction") +
  annotate("text", x = 0.012, y = 0.02, hjust = 0, size = 2.3, colour = "grey30",
           label = paste0("Gene models are drawn at true genomic scale, so CDS blocks narrower than ~2 pixels appear as ",
                          "narrow bars rather than full-width boxes (NR1D1 exon 1: 31 bp CDS + 496 bp 5'UTR)."))

p5C <- (wrap_plots(plots, ncol = 2) / lg) + plot_layout(heights = c(1, 0.105))
ggsave(file.path(pan.dir, "clock_loci_tracks_600.png"), p5C,
       width = 9.6, height = 9.2, dpi = DPI, bg = "white")

st <- rbindlist(stat_rows)
fwrite(st, file.path(res, "clock_locus_track_stats.csv"))
cat("5C rebuilt (v2: tight window, axis titles, UTR/CDS model, legend)\n")
print(st)
