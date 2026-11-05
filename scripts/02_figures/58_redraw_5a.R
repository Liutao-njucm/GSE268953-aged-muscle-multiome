## 58_redraw_5a.R - redraw only the Fig5A GTEx concordance panel
## fix: subtitle was clipped at right edge (one long line) -> wrap to 2 lines, left-align
proj  <- "D:/workbuddy/outputs/2026-10-25-20-23-31"
res   <- file.path(proj, "results/GSE268953")
.libPaths("C:/Users/wang/R/win-library/4.5")
suppressPackageStartupMessages({library(ggplot2); library(ggrepel); library(data.table)})
DPI <- 600
pan.dir <- file.path(res, "panels_600dpi")

m2c <- fread(file.path(res, "gtex_snRNA_concordance.csv"))
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
       subtitle = paste0("FAP SASP/ECM programs concordant; clock arm discordant\n",
                         "(bulk = mixed cell types + time-of-day sampling)")) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom",
        plot.subtitle = element_text(size = 8.5, hjust = 0, colour = "grey25"))

ggsave(file.path(pan.dir, "gtex_concordance_by_module_600.png"), pGT,
       width = 7.8, height = 6, dpi = DPI)
cat("5A REDRAW DONE\n")
