# 18c: module-wise GTEx concordance (honest split: SASP/ECM vs clock vs other)
.libPaths(c("C:/Users/wang/R/win-library/4.5", .libPaths()))
suppressMessages({library(data.table); library(ggplot2); library(ggrepel)})
proj <- "D:/workbuddy/outputs/2026-10-25-20-23-31"
out.dir <- file.path(proj, "results/GSE268953")
m <- fread(file.path(out.dir, "gtex_snRNA_concordance.csv"))

mods <- list(
  Clock_repressive_arm = c("NR1D1","NR1D2","PER1","PER2","PER3","CRY1","CRY2","DBP","TEF","HLF","NR2F1","NR2F2","BHLHE40","BHLHE41","NFIL3","RORA","CIART"),
  SASP_inflammatory = c("IL6","CXCL2","CXCL8","CCL2","CSF1","TGFB1","SERPINE1","IGFBP7","MMP2","TIMP1","CDKN1A","CDKN2A","GLB1","LMNB1"),
  ECM_fibrosis = c("COL1A1","COL3A1","DCN","LUM","POSTN","PDGFRA","PRG4","ADAMTS16","ITGB8","COL19A1","SAA2","ABCA8","MMP3","WISP2","CFD"),
  Muscle_fiber = c("PFKM","PGK1","LMOD2","STAC3","ATP5PD","CRYAB","MUSK")
)
mod_of <- unlist(lapply(names(mods), function(k) setNames(rep(k, length(mods[[k]])), mods[[k]])))
m[, module := ifelse(gene %in% names(mod_of), mod_of[gene], "other")]

res <- list()
for (k in c(names(mods), "other")) {
  d <- m[module == k]
  if (nrow(d) < 5) next
  ct <- suppressWarnings(cor.test(d$sn_lfc, d$lfc_aged_vs_young, method = "spearman"))
  res[[k]] <- data.table(module = k, n = nrow(d),
                         mean_sn_lfc = mean(d$sn_lfc), mean_gtex_lfc = mean(d$lfc_aged_vs_young),
                         spearman_rho = unname(ct$estimate), p = ct$p.value,
                         frac_same_dir = mean(sign(d$sn_lfc) == sign(d$lfc_aged_vs_young)))
}
r <- rbindlist(res)
fwrite(r, file.path(out.dir, "gtex_concordance_by_module.csv"))
cat("=== module-wise concordance ===\n"); print(as.data.frame(r), digits = 3)

pal <- c(Clock_repressive_arm = "#8E44AD", SASP_inflammatory = "#E67E22",
         ECM_fibrosis = "#27AE60", Muscle_fiber = "#2E86C1", other = "grey70")
p <- ggplot(m, aes(sn_lfc, lfc_aged_vs_young, colour = module)) +
  geom_point(size = 2.7, alpha = 0.9) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey50") +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
  ggrepel::geom_text_repel(data = m[module != "other"], aes(label = gene), size = 2.7,
                           max.overlaps = 20, show.legend = FALSE) +
  scale_colour_manual(values = pal, name = NULL) +
  labs(x = "log2FC Aged vs Young — FAP, snMultiome (5v5 donors)",
       y = "log2FC Aged vs Young — GTEx skeletal muscle (bulk, n=818)",
       title = "External check in GTEx: module-wise concordance",
       subtitle = "FAP SASP/ECM programs concordant; clock arm discordant (bulk = mixed cell types + time-of-day sampling)") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")
ggsave(file.path(out.dir, "figures", "gtex_concordance_by_module.png"), p, width = 7.8, height = 6, dpi = 300)
cat("saved gtex_concordance_by_module.png\n")
