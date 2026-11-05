## ============================================================
## 48_fix_panels.R  —  redraw the five problem panels
##   1. Fig 2C  chromatin aging score: donor dots + box (n>=30 filter, no Mast)
##   2. Fig 4B  coupling: two-panel (class fractions | % silent among open)
##   3. Fig 5C  clock locus tracks: both groups in legend, aligned x axes
##   4. Fig 5D  working model: text centred in boxes, genes corrected
##   5. Fig 4D  regulatory chain: colour by TF family, n annotations
##   6. Fig 3C  per-donor FAP: subtitle clarifying IL6 up / structural ECM not
## All panels rendered at 600 dpi into panels_600dpi/ (same file names the
## compose scripts already expect, except the new coupling panel).
## ============================================================
.libPaths(c("C:/Users/wang/R/win-library/4.5", .libPaths()))
suppressMessages({library(data.table); library(ggplot2); library(ggrepel); library(patchwork)})
setDTthreads(4)

res     <- "D:/workbuddy/outputs/2026-10-25-20-23-31/results/GSE268953"
pan.dir <- file.path(res, "panels_600dpi"); dir.create(pan.dir, showWarnings = FALSE)
DPI <- 600
Rf  <- function(...) fread(file.path(res, ...))

## ------------------------------------------------------------
## 1. Fig 2C — chromatin aging score, donor dots + box
## ------------------------------------------------------------
sc <- Rf("cell_aging_score.csv")
sc[, n := .N, by = .(celltype, sample)]
d30 <- sc[n >= 30]
## same rule as script 07: keep cell types with >=3 donors (n>=30 cells) per group
keep_ct <- d30[, .(k = uniqueN(sample)), by = .(celltype, group)][, .(ok = all(k >= 3)), by = celltype][ok == TRUE, celltype]
d30 <- d30[celltype %in% keep_ct]
dsum <- d30[, .(m = mean(aging_score)), by = .(celltype, sample, group)]
dsum[, ctmean := mean(m), by = celltype]
lv <- dsum[, .(ctmean = mean(m)), by = celltype][order(ctmean), celltype]
dsum[, celltype := factor(celltype, levels = lv)]

p2C <- ggplot(dsum, aes(celltype, m, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.62, linewidth = 0.3,
               position = position_dodge(0.72)) +
  geom_point(aes(colour = group), position = position_dodge2(0.72, padding = 0.35),
             size = 1.7, alpha = 0.9) +
  scale_fill_manual(values = c(Young = "#A9CCE3", Aged = "#F5B7B1"), name = NULL) +
  scale_colour_manual(values = c(Young = "#1A5276", Aged = "#922B21"), name = NULL) +
  labs(x = NULL, y = "chromatin aging score (donor mean)",
       title = "ATAC-based chromatin aging score",
       subtitle = "Each dot = one donor (5 v 5); box = median + IQR.\nSMC/Pericyte: young donor means nearly tied, so the box collapses to a line.") +
  theme_bw(base_size = 12) + coord_flip() +
  theme(legend.position = "top")
ggsave(file.path(pan.dir, "chromatin_aging_score_by_celltype_600.png"),
       p2C, width = 8, height = 5.8, dpi = DPI)
cat("2C done\n")

## ------------------------------------------------------------
## 2. Fig 4B — coupling classes (left) + % silent among open (right)
## ------------------------------------------------------------
cpl <- Rf("coupling_gene_level_by_celltype.csv")
pal <- c(Coupled_up = "#B03A2E", Coupled_down = "#1A5276",
         Chromatin_primed = "#E67E22", Chromatin_poised_down = "#5499C7",
         RNA_only_up = "#F5B7B1", RNA_only_down = "#D6EAF8")
bar <- cpl[, .N, by = .(celltype, class)]
bar[, frac := N / sum(N), by = celltype]
pL <- ggplot(bar[class != "Unchanged"], aes(celltype, frac, fill = class)) +
  geom_col(width = 0.72, colour = "grey25", linewidth = 0.15) +
  scale_fill_manual(values = pal, name = "coupling class") +
  scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
  labs(x = NULL, y = "fraction of genes with a linked DA peak",
       title = "Coupling classes among all linked genes") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "bottom")

op <- bar[class %in% c("Coupled_up", "Chromatin_primed"), .(n_open = sum(N)), by = celltype]
tr <- bar[class == "Coupled_up", .(n_out = sum(N)), by = celltype]
dR <- merge(op, tr, by = "celltype")[n_open >= 20]
dR[, silent := 100 * (1 - n_out / n_open)]
dR[, lab := sprintf("%.1f%%\n(n=%d)", silent, n_open)]
cat("silent fractions:\n"); print(dR)
pR <- ggplot(dR, aes(reorder(celltype, -silent), silent)) +
  geom_col(fill = "#E67E22", width = 0.62, colour = "grey25", linewidth = 0.15) +
  geom_hline(yintercept = 90, linetype = 2, colour = "grey40") +
  geom_text(aes(label = lab), vjust = -0.25, size = 3.1, lineheight = 0.85) +
  scale_y_continuous(limits = c(0, 108), expand = c(0, 0)) +
  labs(x = NULL, y = "% of age-open-enhancer genes",
       title = "Open-enhancer genes WITHOUT concordant\nRNA up-regulation (>90%, dashed line)") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
p4B <- pL + pR + plot_layout(widths = c(1.3, 1))
ggsave(file.path(pan.dir, "coupling_class_fraction_600.png"),
       p4B, width = 9.8, height = 4.7, dpi = DPI)
cat("4B done\n")

## ------------------------------------------------------------
## 3. Fig 5C — clock locus tracks, both groups in legend, aligned axes
## ------------------------------------------------------------
suppressMessages({library(Matrix); library(GenomicRanges); library(jsonlite)})
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
  j <- fromJSON(file.path("D:/workbuddy/outputs/2026-10-25-20-23-31/data/annotation",
                          paste0("refGene_", chr, ".json")), simplifyDataFrame = TRUE)
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
  d <- data.table(mid = (sel$start + sel$end) / 2,
                  Young = fap_young[sel$idx], Aged = fap_aged[sel$idx])
  win_counts <- colSums(cpm[sel$idx, , drop = FALSE])
  y <- log2(win_counts[grp == "Young"] + 1); a <- log2(win_counts[grp == "Aged"] + 1)
  p.w <- tryCatch(wilcox.test(a, y)$p.value, error = function(e) NA)
  long <- rbind(data.table(mid = d$mid, cpm = d$Young, grp = "Young"),
                data.table(mid = d$mid, cpm = d$Aged, grp = "Aged"))
  brk <- pretty(c(L$start, L$end), n = 3)
  brk <- brk[brk >= L$start & brk <= L$end]
  p1 <- ggplot(long, aes(mid, cpm)) +
    geom_ribbon(aes(ymin = 0, ymax = cpm, fill = grp), alpha = 0.4, colour = NA) +
    geom_line(aes(colour = grp), linewidth = 0.55) +
    scale_fill_manual(values = c(Young = "#2E86C1", Aged = "#C0392B"), name = NULL) +
    scale_colour_manual(values = c(Young = "#2E86C1", Aged = "#C0392B"), name = NULL) +
    scale_x_continuous(limits = c(L$start, L$end), breaks = brk, expand = c(0.01, 0),
                       labels = function(x) sprintf("%.2f Mb", x / 1e6)) +
    labs(y = "accessibility (CPM)", x = NULL,
         title = sprintf("%s : FAP chromatin accessibility", L$gene),
         subtitle = sprintf("locus log2FC = %.2f, Wilcoxon P = %.2g",
                            log2((mean(d$Aged) + .1) / (mean(d$Young) + .1)), p.w)) +
    theme_bw(base_size = 9) +
    theme(panel.grid.minor = element_blank(), legend.position = "top",
          legend.key.size = unit(0.35, "cm"), plot.title = element_text(size = 9, face = "bold"),
          axis.text.x = element_blank(), axis.ticks.x = element_blank())
  gm <- gene_model(L$chr, L$gene)
  p2 <- ggplot(gm) +
    geom_segment(aes(x = min(txStart), xend = max(txEnd), y = 0, yend = 0), linewidth = 0.4) +
    geom_rect(aes(xmin = start, xmax = end, ymin = -0.2, ymax = 0.2), fill = "#34495E", colour = NA) +
    scale_x_continuous(limits = c(L$start, L$end), breaks = brk, expand = c(0.01, 0),
                       labels = function(x) sprintf("%.2f Mb", x / 1e6)) +
    scale_y_continuous(limits = c(-0.6, 0.6)) +
    labs(x = NULL, y = NULL) +
    theme_bw(base_size = 9) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          panel.grid = element_blank(), axis.line.y = element_blank(),
          plot.margin = margin(1, 5, 1, 1))
  plots[[i]] <- (p1 / p2) + plot_layout(heights = c(3, 1))
}
p5C <- wrap_plots(plots, ncol = 2)
ggsave(file.path(pan.dir, "clock_loci_tracks_600.png"), p5C, width = 9.6, height = 6.6, dpi = DPI)
cat("5C done\n")

## ------------------------------------------------------------
## 4. Fig 5D — working model, text centred, genes corrected
## ------------------------------------------------------------
sch <- data.frame(
  x = c(0.02, 0.36, 0.36, 0.70, 0.36, 0.70, 0.02),
  y = c(0.60, 0.86, 0.34, 0.86, 0.60, 0.34, 0.08),
  w = c(0.30, 0.28, 0.28, 0.28, 0.28, 0.28, 0.96),
  h = c(0.16, 0.14, 0.14, 0.14, 0.14, 0.14, 0.13),
  lab = c("Ageing + disuse\n(motor-neuron decline)",
          "Clock repressor-arm\ncollapse NR1D1/PER/DBP\n(FAP, padj 3e-9)",
          "Loss of clock licence\n-> AP-1 / CEBP / NFIA\ntake over enhancers",
          "SASP / fibrosis enhancers\nopen (1,122 DA peaks)",
          "Chromatin primed\n>90% not transcribed\n(OR 1.85-2.55)",
          "FAP senescence &\nfibrotic niche\n(IL6, PRG4, ADAMTS16)",
          "Myonuclear niche: glycolytic / sarcomeric decline ; MUSK up-regulation"),
  fill = c("#E8E8E8", "#D6EAF8", "#FCF3CF", "#FADBD8", "#FADBD8", "#F5B7B1", "#E8F8F5"),
  stringsAsFactors = FALSE)
pMO <- ggplot(sch, aes(x, y)) +
  geom_rect(aes(xmin = x, xmax = x + w, ymin = y - h/2, ymax = y + h/2, fill = fill),
            colour = "grey30", linewidth = 0.35) +
  geom_text(aes(x = x + w/2, y = y, label = lab), size = 2.45, lineheight = 0.85) +
  annotate("segment", x = 0.32, xend = 0.36, y = 0.66, yend = 0.84, linewidth = 0.5,
           arrow = arrow(length = unit(2.2, "mm"))) +
  annotate("segment", x = 0.50, xend = 0.50, y = 0.79, yend = 0.70, linewidth = 0.5,
           arrow = arrow(length = unit(2.2, "mm"))) +
  annotate("segment", x = 0.64, xend = 0.70, y = 0.84, yend = 0.66, linewidth = 0.5,
           arrow = arrow(length = unit(2.2, "mm"))) +
  annotate("segment", x = 0.64, xend = 0.70, y = 0.36, yend = 0.54, linewidth = 0.5,
           arrow = arrow(length = unit(2.2, "mm"))) +
  annotate("segment", x = 0.50, xend = 0.50, y = 0.53, yend = 0.41, linewidth = 0.5,
           arrow = arrow(length = unit(2.2, "mm"))) +
  scale_fill_identity() + xlim(0, 1.02) + ylim(0, 1) +
  labs(title = "Working model: clock-licensed senescence enhancers drive the aged FAP niche") +
  theme_void(11) + theme(plot.title = element_text(face = "bold", size = 11, hjust = 0.5))
ggsave(file.path(pan.dir, "p_model_600.png"), pMO, width = 7.4, height = 3.2, dpi = DPI)
cat("model done\n")

## ------------------------------------------------------------
## 5. Fig 4D — regulatory chain by TF family
## ------------------------------------------------------------
rc <- Rf("regulatory_chain_FAP.csv")
rc.sig <- rc[p_binom < 0.05 & DA_direction == "up"][order(median_rna_lfc)]
fam <- function(t) {
  ifelse(t %in% c("NR1D1","NR1D2","NR2F1","NR2F2","DBP","TEF","HLF","NFIL3"), "Clock",
  ifelse(t %in% c("FOS","JUN","JUNB","JUND","FOSL1","FOSL2","BATF","MAFK","MAFB","MAF"), "AP-1",
  ifelse(t %in% c("CEBPA","CEBPB","CEBPD","CEBPE"), "CEBP", "Other")))
}
rc.sig[, family := fam(TF)]
rc.sig[, TF := factor(TF, levels = TF)]
rc.sig[, nlab := paste0("n=", n_linked_genes)]
fpal <- c(Clock = "#8E44AD", "AP-1" = "#C0392B", CEBP = "#E67E22", Other = "#5D6D7E")
pCH <- ggplot(rc.sig, aes(TF, frac_rna_up, fill = family)) +
  geom_col(width = 0.72, colour = "grey25", linewidth = 0.15) +
  geom_hline(yintercept = 0.5, linetype = 2, colour = "grey40") +
  geom_text(aes(label = nlab), hjust = -0.12, size = 2.3, colour = "grey30") +
  scale_fill_manual(values = fpal, name = "TF family") +
  coord_flip(clip = "off") +
  labs(x = NULL, y = "Fraction of linked target genes up-regulated in Aged FAP",
       title = "TF -> open enhancer -> target-gene chain (FAP)",
       subtitle = "Significant chains (exact binomial p < 0.05;\nall also pass Fisher's exact test);\nn = linked genes; dashed line = 49.3% background") +
  theme_classic(11) +
  theme(legend.position = "top", plot.title = element_text(face = "bold", size = 11))
ggsave(file.path(pan.dir, "p_regchain_600.png"), pCH, width = 5.4, height = 4.6, dpi = DPI)
cat("4D done\n")

## ------------------------------------------------------------
## 6. Fig 3C — per-donor FAP, clarifying subtitle
## ------------------------------------------------------------
d <- Rf("FAP_clock_per_donor.csv")
d <- d[!is.na(donor)]
d[, gene := factor(gene, levels = c("NR1D1","NR1D2","PER1","PER3","DBP","TEF","HLF","CRY2",
                                    "ARNTL","NPAS2","IL6","CCL2","LUM","COL1A1","PRG4","PDGFRA"))]
d <- d[!is.na(gene)]
d[, lab := ifelse(grepl("^Young", donor), "Young (20-27 y)", "Aged (60-77 y)")]
pFP <- ggplot(d, aes(lab, log2(CPM10k + 0.05), colour = lab)) +
  geom_point(position = position_jitter(width = 0.12), size = 2, alpha = 0.8) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.5, colour = "black") +
  facet_wrap(~ gene, scales = "free_y", ncol = 8) +
  scale_colour_manual(values = c(`Young (20-27 y)` = "#2E86C1", `Aged (60-77 y)` = "#C0392B"), name = NULL) +
  labs(x = NULL, y = "log2 expression (CPM/10k, FAP pseudobulk)",
       title = "Per-donor expression in FAP: clock repression is donor-consistent",
       subtitle = "Repressor arm falls 5/5 donors; IL6 and PRG4 rise (padj < 0.05); CCL2 trends up (ns); structural ECM genes (LUM, COL1A1) are NOT induced within FAP - bulk ECM signal comes from other sources (see Fig. 5a, b)") +
  theme_bw(base_size = 9) +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        plot.subtitle = element_text(size = 7.5))
ggsave(file.path(pan.dir, "FAP_clock_per_donor_600.png"), pFP, width = 10.5, height = 4.8, dpi = DPI)
cat("3C done\n")
cat("ALL PANELS DONE\n")
