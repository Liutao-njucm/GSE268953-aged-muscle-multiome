# ============================================================
# 57 (v2): teaching figure for Fig 5C — NR1D1 locus, numbered callouts
#   All annotations are placed in DATA (genomic) coordinates.
#   Strand-aware: NR1D1 is on the minus strand, so upstream (5') is at
#   HIGHER coordinates, i.e. to the RIGHT of the TSS.
#     (1) promoter peak, straddles the TSS
#     (2) strongest peak inside the gene body (intron 1 enhancer) 1.3 kb
#         downstream of the TSS — the tall bar, NOT the promoter
#     (3) proximal 5' flanking peak (7.8 kb upstream = to the RIGHT)
#     (4) distal peak near the window edge (11 kb upstream)
#     (5) exon boxes (not an ATAC feature)
# ============================================================
.libPaths(c("C:/Users/wang/R/win-library/4.5", .libPaths()))
suppressMessages({library(Matrix); library(data.table); library(ggplot2)
  library(GenomicRanges); library(jsonlite); library(patchwork)})

proj <- "D:/workbuddy/outputs/2026-10-25-20-23-31"
res  <- file.path(proj, "results/GSE268953")
DPI  <- 600

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

## ---- NR1D1 locus ----
j  <- fromJSON(file.path(proj, "data/annotation", "refGene_chr17.json"), simplifyDataFrame = TRUE)
rg <- as.data.table(j$refGene)[name2 == "NR1D1"][order(-(txEnd - txStart))][1]
txStart <- rg$txStart; txEnd <- rg$txEnd; strand <- rg$strand
tss <- if (strand == "+") txStart else txEnd
dir <- if (strand == "+") 1 else -1
pad <- 12000
win <- c(txStart - pad, txEnd + pad)

s_ex <- as.integer(strsplit(gsub(",$", "", rg$exonStarts), ",")[[1]])
e_ex <- as.integer(strsplit(gsub(",$", "", rg$exonEnds),   ",")[[1]])
ex <- data.table(start = s_ex, end = e_ex)
ex[, `:=`(cs = pmax(start, rg$cdsStart), ce = pmin(end, rg$cdsEnd))]

sel <- pk[chr == "chr17" & end > win[1] & start < win[2]]
sel[, mid := (start + end) / 2]
sel[, Young := fap_young[idx]]; sel[, Aged := fap_aged[idx]]
## strand-aware offset: positive = upstream (5') of the TSS in transcription space
sel[, up := if (strand == "-") mid - tss else tss - mid]
sel[, in_body := mid > txStart & mid < txEnd]

i_pro <- which.min(abs(sel$mid - tss))
prom  <- sel[i_pro]

## (2) strongest peak inside the gene body, excluding the promoter window
c2 <- sel[in_body == TRUE & abs(up) > 1000]
p2 <- c2[which.max(c2$Aged)]
## (3) proximal 5' flanking peak: 1-10 kb upstream
c3 <- sel[up > 1000 & up <= 10000]
p3 <- c3[which.max(c3$Aged)]
## (4) distal peak: >10 kb upstream (near the window edge)
c4 <- sel[up > 10000]
p4 <- c4[which.max(c4$Aged)]

ymax <- max(sel$Young, sel$Aged) * 1.32

ann <- data.table(
  n    = 1:4,
  x    = c(prom$mid, p2$mid, p3$mid, p4$mid),
  ytop = c(prom$Aged, p2$Aged, p3$Aged, p4$Aged))
ann[, ytxt := ytop + 9]
ann[2, ytxt := ymax * 0.975]
ann[3, ytxt := ymax * 0.80]
ann[4, ytxt := ymax * 0.935]

brk <- pretty(win, n = 5); brk <- brk[brk >= win[1] & brk <= win[2]]

track <- function(dat, col, lab, show_y) {
  ggplot(dat) +
    annotate("rect", xmin = tss - 1000, xmax = tss + 1000, ymin = -Inf, ymax = Inf,
             fill = "grey55", alpha = 0.13) +
    geom_rect(aes(xmin = start, xmax = end, ymin = 0, ymax = .data[[lab]]),
              fill = col, alpha = 0.85, colour = NA) +
    geom_vline(xintercept = tss, linetype = "22", linewidth = 0.35, colour = "grey15") +
    annotate("text", x = -Inf, y = Inf, label = lab, hjust = -0.2, vjust = 1.4,
             size = 3.1, colour = col, fontface = "bold") +
    geom_point(data = ann, aes(x = x, y = ytxt), shape = 21, size = 4.4,
               fill = "white", stroke = 0.6, colour = "grey10") +
    geom_text(data = ann, aes(x = x, y = ytxt, label = n), size = 2.9,
              fontface = "bold", colour = "grey10") +
    geom_segment(data = ann[n != 2], inherit.aes = FALSE,
                 aes(x = x, y = ytxt - 0.085 * ymax, xend = x, yend = ytop + 1.5),
                 linewidth = 0.3, colour = "grey35") +
    scale_x_continuous(limits = win, breaks = brk, expand = c(0.004, 0)) +
    scale_y_continuous(limits = c(0, ymax), expand = expansion(mult = c(0, 0))) +
    labs(y = if (show_y) "ATAC-seq\n(CPM)" else NULL, x = NULL) +
    theme_bw(base_size = 8.5) +
    theme(panel.grid = element_blank(),
          axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          plot.margin = margin(1, 3, 0, 1))
}
pA <- track(sel, "#C0392B", "Aged", TRUE) +
  labs(title = sprintf("NR1D1 locus (%s strand) \u2014 how to read the tracks", strand)) +
  theme(plot.title = element_text(size = 10.5, face = "bold", hjust = 0),
        plot.margin = margin(4, 3, 0, 1))
pY <- track(sel, "#2E86C1", "Young", FALSE)

## ---- gene model ----
ex_mid <- ex[which.max(end - start)]
pG <- ggplot() +
  annotate("rect", xmin = tss - 1000, xmax = tss + 1000, ymin = -Inf, ymax = Inf,
           fill = "grey55", alpha = 0.13) +
  annotate("rect", xmin = txStart, xmax = txEnd, ymin = -0.055, ymax = 0.055,
           fill = "grey72", colour = NA) +
  annotate("rect", xmin = ex[cs < ce, cs], xmax = ex[cs < ce, ce],
           ymin = -0.26, ymax = 0.26, fill = "#1B2631", colour = NA) +
  annotate("segment", x = txStart + (txEnd - txStart) * 0.16, xend = txStart, y = 0, yend = 0,
           linewidth = 0.5, colour = "grey10",
           arrow = arrow(length = unit(1.7, "mm"), type = "closed")) +
  annotate("point", x = ex_mid$cs + 100, y = 0.26, shape = 21, size = 4.4,
           fill = "white", stroke = 0.6, colour = "grey10") +
  annotate("text", x = ex_mid$cs + 100, y = 0.26, label = "5", size = 2.9,
           fontface = "bold", colour = "grey10") +
  ## TSS marker + 5'/3' labels
  annotate("segment", x = tss, xend = tss + dir * diff(win) * 0.05, y = 0.62, yend = 0.62,
           linewidth = 0.4, colour = "grey10",
           arrow = arrow(length = unit(1.6, "mm"), type = "closed")) +
  annotate("text", x = tss, y = 0.80, label = "TSS", hjust = 0.5, size = 3.0,
           colour = "grey10", fontface = "bold") +
  annotate("text", x = tss + dir * diff(win) * 0.085, y = 0.62,
           label = "5'", size = 2.8, colour = "grey25", hjust = 0.5) +
  annotate("text", x = txStart - dir * diff(win) * 0.015, y = 0.62,
           label = "3'", size = 2.8, colour = "grey25", hjust = 0.5) +
  geom_vline(xintercept = tss, linetype = "22", linewidth = 0.35, colour = "grey15") +
  scale_x_continuous(limits = win, breaks = brk, expand = c(0.004, 0),
                     labels = function(x) sprintf("%.3f", x / 1e6)) +
  scale_y_continuous(limits = c(-0.55, 1.15)) +
  labs(x = "chr17 position (Mb)", y = NULL) +
  theme_bw(base_size = 8.5) +
  theme(panel.grid = element_blank(), axis.text.y = element_blank(),
        axis.ticks.y = element_blank(), axis.line.y = element_blank(),
        plot.margin = margin(0, 3, 2, 1))

## ---- key ----
key <- data.table(n = 1:5, txt = c(
  "promoter peak straddling the TSS (Aged 33.6 vs Young 38.0, P = 0.15)",
  "strongest peak of the locus, but 1.3 kb INSIDE the gene body (intron-1 enhancer) = 3' of the TSS, not the promoter",
  "proximal 5' flanking peak, 7.8 kb upstream - on the RIGHT of the TSS because NR1D1 is on the minus strand",
  "distal peak 11 kb upstream near the window edge - still not this gene's promoter",
  "exon boxes: ATAC peaks are not expected on exons, so peaks never have to 'sit on' an exon"))
key[, ypos := 1 - (seq_len(.N) - 0.5) / .N]
klab <- ggplot() + xlim(0, 1) + ylim(0, 1) + theme_void() +
  geom_text(data = key, aes(x = 0.012, y = ypos, label = paste0(n, "   ", txt)),
            hjust = 0, size = 2.6, colour = "grey10")

out <- (pA / pY / pG / klab) + plot_layout(heights = c(1, 1, 0.85, 0.30))
ggsave(file.path(res, "teaching_NR1D1_howto.png"), out,
       width = 7.6, height = 6.9, dpi = DPI, bg = "white")
cat("teaching figure written (v2, strand-aware)\n")
cat(sprintf("TSS=%d strand=%s\n (1) promoter  mid=%d  up=%+.0f\n (2) intronic  mid=%d  up=%+.0f  in_body=%s\n (3) prox 5'   mid=%d  up=%+.0f\n (4) distal    mid=%d  up=%+.0f\n",
            tss, strand, prom$mid, prom$up, p2$mid, p2$up, p2$in_body, p3$mid, p3$up, p4$mid, p4$up))
