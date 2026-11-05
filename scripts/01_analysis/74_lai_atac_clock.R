# 74_lai_atac_clock.R --- Lai OMIX004305 snATAC 时钟位点可及性复制分析
# 主检验 = FAP（CD55+/GPC3+/MME+，排除 Tenocyte）供体级
#   每基因两口径：gene window（body±12kb）与 promoter（TSS±2kb）
#   log2(CPM10k+1) 供体伪批量；MWU old vs young（双尾+单尾方向）+ Spearman vs age
#   模块分 = 7 抑制臂基因 z 均值；敏感性：>=50 核 / 男-only / China-only
# 次要 = MuSC（三态合并）
# 输出 results/GSE268953/lai_atac/
.libPaths("C:/Users/wang/R/win-library/4.5")
suppressMessages({
  library(SummarizedExperiment); library(GenomicRanges); library(jsonlite); library(Matrix)
})
options(width = 200)
set.seed(1)

DIR <- "data/Lai2024_ATAC"
OUT <- "results/GSE268953/lai_atac"; dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

REP <- c("NR1D1","NR1D2","PER1","PER3","DBP","HLF","TEF")
ACT <- c("ARNTL","NPAS2","CLOCK","PER2","CRY1")
GENES <- c(REP, ACT)

## ---- 基因窗口（hg38 refGene）----
files <- Sys.glob("data/annotation/refGene_chr*.json")
genes <- do.call(rbind, lapply(files, function(f) {
  d <- fromJSON(f)$refGene
  d <- d[d$name2 %in% GENES, c("name2","chrom","txStart","txEnd","strand")]
  d[!duplicated(d$name2), ]
}))
rownames(genes) <- genes$name2
genes <- genes[GENES, ]
genes$tss <- ifelse(genes$strand == "+", genes$txStart, genes$txEnd)
cat("gene windows:\n")
print(genes[, c("chrom","txStart","txEnd","strand","tss")])

gr_body <- GRanges(genes$chrom, IRanges(pmax(0, genes$txStart - 12000), genes$txEnd + 12000))
gr_prom <- GRanges(genes$chrom, IRanges(pmax(0, genes$tss - 2000), genes$tss + 2000))
names(gr_body) <- names(gr_prom) <- GENES

## ---- 供体级位点可及性 ----
locus_donor <- function(se, keep_anno) {
  cd <- as.data.frame(colData(se))
  cd$sample <- as.character(cd$sample); cd$group <- as.character(cd$group)
  cd$Sex <- as.character(cd$Sex); cd$country <- as.character(cd$country)
  cd$age <- suppressWarnings(as.numeric(as.character(cd$age)))
  keep <- cd$Annotation %in% keep_anno
  se2 <- se[, keep]; cd2 <- cd[keep, ]
  peaks <- rowRanges(se)
  M <- assay(se2, "PeakMatrix")
  tot <- Matrix::colSums(M)
  okc <- !is.na(tot) & tot > 0
  M <- M[, okc]; tot <- tot[okc]; cd2 <- cd2[okc, ]
  donors <- unique(cd2[, c("sample","age","Sex","country","group")])
  donors <- donors[order(donors$sample), ]
  ncell <- table(cd2$sample)
  loci_body <- list(); loci_prom <- list()
  for (g in GENES) {
    hb <- unique(queryHits(findOverlaps(peaks, gr_body[g])))
    hp <- unique(queryHits(findOverlaps(peaks, gr_prom[g])))
    dcnt_b <- tapply(Matrix::colSums(M[hb, , drop = FALSE]), cd2$sample, sum)
    dcnt_p <- tapply(Matrix::colSums(M[hp, , drop = FALSE]), cd2$sample, sum)
    dtot   <- tapply(tot, cd2$sample, sum)
    v_b <- dcnt_b[names(dtot)] / dtot          # 按供体名对齐，保留名字
    v_p <- dcnt_p[names(dtot)] / dtot
    loci_body[[g]] <- log2(v_b * 1e4 + 1)
    loci_prom[[g]] <- log2(v_p * 1e4 + 1)
  }
  names(loci_body) <- GENES; names(loci_prom) <- GENES
  list(donors = donors, ncell = ncell, body = loci_body, prom = loci_prom)
}

mk_mat <- function(L, which) {
  ll <- L[[which]]
  M <- sapply(ll, function(x) x[L$donors$sample])
  rownames(M) <- L$donors$sample
  colnames(M) <- names(ll)
  M
}

run_test <- function(M, donors, tag) {
  old <- donors$group == "old"; yng <- !old
  out <- list()
  for (g in colnames(M)) {
    v <- M[, g]; o <- v[old]; y <- v[yng]
    o <- o[!is.na(o)]; y <- y[!is.na(y)]
    if (length(o) < 3 || length(y) < 2) { cat("  [skip]", tag, g, "n_old =", length(o), "n_young =", length(y), "\n"); next }
    w  <- suppressWarnings(wilcox.test(o, y, exact = FALSE))
    wl <- suppressWarnings(wilcox.test(o, y, alternative = "less", exact = FALSE))
    wg <- suppressWarnings(wilcox.test(o, y, alternative = "greater", exact = FALSE))
    sp <- suppressWarnings(cor.test(donors$age, v, method = "spearman"))
    out[[g]] <- data.frame(gene = g,
      arm = ifelse(g %in% REP, "repressor", "activator"), analysis = tag,
      n_old = sum(old), n_young = sum(yng),
      med_old = median(o), med_young = median(y), delta = median(o) - median(y),
      p_two = w$p.value, p_less = wl$p.value, p_greater = wg$p.value,
      rho = unname(sp$estimate), p_age = sp$p.value)
  }
  do.call(rbind, out)
}

module_stats <- function(M, donors, gene_set, tag) {
  z <- scale(M[, gene_set])
  mod <- rowMeans(z, na.rm = TRUE)
  d <- data.frame(donor = rownames(M), age = donors$age, sex = donors$Sex,
                  country = donors$country, group = donors$group,
                  module = as.numeric(mod), stringsAsFactors = FALSE)
  old <- d$group == "old"; yng <- !old
  wt <- suppressWarnings(wilcox.test(d$module[old], d$module[yng], exact = FALSE))
  wl <- suppressWarnings(wilcox.test(d$module[old], d$module[yng],
                                     alternative = "less", exact = FALSE))
  wg <- suppressWarnings(wilcox.test(d$module[old], d$module[yng],
                                     alternative = "greater", exact = FALSE))
  sp <- suppressWarnings(cor.test(d$age, d$module, method = "spearman"))
  stat <- data.frame(analysis = tag, n_old = sum(old), n_young = sum(yng),
                     med_old = median(d$module[old]), med_young = median(d$module[yng]),
                     p_two = wt$p.value, p_less = wl$p.value, p_greater = wg$p.value,
                     rho = unname(sp$estimate), p_age = sp$p.value)
  d$analysis <- tag
  list(donors = d, stat = stat)
}

analyse <- function(rds, keep_anno, label) {
  cat("\n=====", label, "=====\n")
  se <- readRDS(file.path(DIR, rds))
  L  <- locus_donor(se, keep_anno)
  L$donors$n_fap_nuclei <- as.integer(L$ncell[L$donors$sample])
  Mbody <- mk_mat(L, "body"); Mprom <- mk_mat(L, "prom")
  cat("donors:", nrow(L$donors), " cells kept:", sum(L$ncell), "\n")
  print(L$donors)
  out_dir <- file.path(OUT, label); dir.create(out_dir, showWarnings = FALSE)
  write.csv(L$donors, file.path(out_dir, "donors.csv"), row.names = FALSE)

  allrows <- rbind(run_test(Mbody, L$donors, "gene_window"),
                   run_test(Mprom, L$donors, "promoter"))
  sel50 <- L$donors$n_fap_nuclei >= 50
  selM  <- L$donors$Sex == "Male"
  selC  <- L$donors$country == "China"
  sens <- list(ge50_nuclei = sel50, male_only = selM, china_only = selC)
  mod_sens <- list()
  for (nm in names(sens)) {
    sel <- sens[[nm]]
    if (sum(sel) >= 7 && sum(!sel) >= 2) {
      allrows  <- rbind(allrows, run_test(Mbody[sel, ], L$donors[sel, ], nm))
      mod_sens[[nm]] <- module_stats(Mbody[sel, ], L$donors[sel, ],
                                     REP, paste0("repressor_module_", nm))$stat
    } else cat("sensitivity", nm, "skipped: n_old =", sum(sel & L$donors$group=="old"),
               "n_young =", sum(sel & !L$donors$group=="old"), "\n")
  }
  mods <- list(
    rep_window = module_stats(Mbody, L$donors, REP, "repressor_module_gene_window"),
    rep_prom   = module_stats(Mprom, L$donors, REP, "repressor_module_promoter"),
    act_window = module_stats(Mbody, L$donors, ACT, "activator_module_gene_window")
  )
  mod_tab <- rbind(mods$rep_window$stat, mods$rep_prom$stat,
                   mods$act_window$stat, do.call(rbind, mod_sens))
  write.csv(allrows, file.path(out_dir, "per_gene_tests.csv"), row.names = FALSE)
  write.csv(mod_tab, file.path(out_dir, "module_tests.csv"), row.names = FALSE)
  dn <- rbind(mods$rep_window$donors, mods$act_window$donors)
  write.csv(dn, file.path(out_dir, "donor_modules.csv"), row.names = FALSE)
  cat("\n-- per-gene (gene window, all donors) --\n")
  print(allrows[allrows$analysis == "gene_window",
                c("gene","arm","n_old","n_young","delta","p_two","p_less","rho","p_age")])
  cat("\n-- module tests --\n"); print(mod_tab)
  invisible(list(per_gene = allrows, module = mod_tab, donors = L$donors))
}

res_fap  <- analyse("FAPs_snATAC_11486.RDS",
                    c("CD55+ FAP","GPC3+ FAP","MME+ FAP"), "FAP")
res_musc <- analyse("MuSC_snATAC_2275.RDS",
                    c("Diff.MuSC","Early Primed MuSC","Quiescent MuSC"), "MuSC")

write.csv(rbind(res_fap$per_gene, res_musc$per_gene),
          file.path(OUT, "lai_atac_per_gene_all.csv"), row.names = FALSE)
write.csv(rbind(res_fap$module, res_musc$module),
          file.path(OUT, "lai_atac_module_all.csv"), row.names = FALSE)
cat("\nDONE -> ", OUT, "\n")
