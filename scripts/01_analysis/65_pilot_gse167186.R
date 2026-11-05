## 65_pilot_gse167186.R --- GSE167186 pilot: detection-rate gate for clock repressor arm in FAP
## 1 Young (HM5, age 20) + 1 Old (HM1, age 75), vastus lateralis, 10x 5' v1 snRNA
## Gate: >=5 of 7 repressor genes detected in >=10% of FAP nuclei  -> proceed to full 17-donor run
.libPaths(c("C:/Users/wang/R/win-library/4.5", .libPaths()))
suppressPackageStartupMessages({library(Seurat); library(data.table)})
set.seed(1)
d.save <- function(x, f) fwrite(x, file.path("results/GSE167186_pilot", f))

dir.create("results/GSE167186_pilot", showWarnings = FALSE)

files <- c(Young = "HM5_Young20", Old = "HM1_Old75")
obj <- lapply(names(files), function(g) {
  d <- file.path("data/GSE167186", files[[g]])
  m <- ReadMtx(mtx = file.path(d, "matrix.mtx"), cells = file.path(d, "barcodes.tsv"),
               features = file.path(d, "features.tsv.gz"), feature.column = 2)
  o <- CreateSeuratObject(m, project = g, min.cells = 3, min.features = 200)
  o$group <- g
  o
})
sc <- merge(obj[[1]], obj[[2]], add.cell.ids = names(files))
sc <- JoinLayers(sc)
cat("cells per sample:\n"); print(table(sc$group))

## QC (nuclei: mitochondrial fraction normally low)
sc[["percent.mt"]] <- PercentageFeatureSet(sc, pattern = "^MT-")
print(summary(sc$nFeature_RNA))
print(summary(sc$percent.mt))
sc <- subset(sc, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & percent.mt < 15)
cat("after QC:", ncol(sc), "cells\n")

sc <- NormalizeData(sc, verbose = FALSE) |>
  FindVariableFeatures(nfeatures = 2500, verbose = FALSE) |>
  ScaleData(verbose = FALSE) |>
  RunPCA(npcs = 30, verbose = FALSE) |>
  RunUMAP(dims = 1:20, verbose = FALSE) |>
  FindNeighbors(dims = 1:20, verbose = FALSE) |>
  FindClusters(resolution = 0.5, verbose = FALSE)
cat("clusters:", length(unique(Idents(sc))), "\n")

## ---- marker-based annotation ----------------------------------------------------
mark <- list(
  FAP   = c("PDGFRA", "DCN", "LUM", "COL1A1", "PI16", "MFAP5"),
  Myo   = c("MYH2", "MYH7", "ACTA1", "TNNT3", "TNNI2"),
  EC    = c("PECAM1", "VWF", "CLDN5"),
  SMC   = c("ACTA2", "TAGLN", "MYH11"),
  MuSC  = c("MYOD1", "MYF5", "MYCN"),
  Mac   = c("CD68", "LYZ", "C1QA"),
  TNK   = c("CD3D", "NKG7", "GNLY")
)
cm <- AggregateExpression(sc, features = intersect(unique(unlist(mark)), rownames(sc)), group.by = "seurat_clusters")$RNA
zs <- t(scale(t(cm)))
anno <- sapply(colnames(cm), function(cl) {
  s <- sapply(mark, function(g) {
    gg <- intersect(g, rownames(zs)); if (!length(gg)) return(-Inf)
    mean(zs[gg, cl])
  })
  names(which.max(s))
})
cat("\ncluster marker z-scores (context only, too coarse at this depth):\n")
print(round(t(zs), 2))

## marker-gated FAP: robust at shallow depth (clustering too coarse here)
cnt <- GetAssayData(sc, layer = "counts")
has <- function(g) Matrix::colSums(cnt[intersect(g, rownames(cnt)), , drop = FALSE]) > 0
fap_like  <- has("PDGFRA") | (has("DCN") & has("LUM"))
non_fap   <- has(c("MYH2", "MYH7", "ACTA1", "TNNT3")) | has(c("PECAM1", "VWF")) | has(c("ACTA2", "TAGLN"))
sc$is_fap <- fap_like & !non_fap
cat("\nFAP-gated nuclei:\n"); print(table(sc$is_fap, sc$group))

## ---- the gate: detection at TWO levels ------------------------------------------
repress <- c("NR1D1", "NR1D2", "PER1", "PER3", "DBP", "HLF", "TEF")
activ   <- c("ARNTL", "NPAS2", "CLOCK", "PER2", "CRY1")

getdet <- function(genes, o) {
  g2 <- intersect(genes, rownames(o))
  mm <- GetAssayData(o, layer = "counts")[g2, , drop = FALSE]
  data.table(gene = g2,
             det_pct = round(Matrix::rowMeans(mm > 0) * 100, 1),
             mean_cnt = round(Matrix::rowMeans(mm), 3))
}
fap <- subset(sc, subset = is_fap)
cat("\nFAP-gated nuclei:", ncol(fap), "\n")
dr <- getdet(c(repress, activ), fap)
dr$arm <- ifelse(dr$gene %in% repress, "repressor", "activator")
dr <- dr[order(arm, -det_pct)]
print(dr)
d.save(dr, "pilot_clock_detection_FAP.csv")
n_pass <- sum(dr$arm == "repressor" & dr$det_pct >= 10)
cat("\nLEVEL-1 gate (per-nucleus det>=10%):", n_pass, "/ 7\n")

## pseudobulk per group (this is the actual full-run primary endpoint)
fap <- AddMetaData(fap, fap$group, col.name = "grp")
pbm <- AggregateExpression(fap, features = c(repress, activ), group.by = "grp")$RNA
pb <- round(log1p(1e4 * pbm / rep(colSums(pbm), each = nrow(pbm))), 2)
out <- data.table(gene = rownames(pb), Young = pb[, 1], Old = pb[, 2])
out$direction <- ifelse(out$Old < out$Young, "down_in_old", "up_in_old")
cat("\npseudobulk log1p(CPM10k), FAP-gated:\n"); print(out)
d.save(out, "pilot_pseudobulk_direction.csv")
n_pass2 <- sum(out$direction == "down_in_old" & out$gene %in% repress)
cat("\nLEVEL-2 gate (pseudobulk signal, >=5/7 repressor nonzero & directional):",
    sum(out$gene %in% repress & (out$Young > 0 | out$Old > 0)), "/ 7 nonzero;",
    n_pass2, "/ 7 down_in_old\n")

## per-gene detection across ALL nuclei (context: maybe signal sits elsewhere)
drall <- getdet(c(repress, activ), sc)
drall$arm <- ifelse(drall$gene %in% repress, "repressor", "activator")
d.save(drall[order(arm, -det_pct)], "pilot_clock_detection_allcells.csv")

saveRDS(sc, "results/GSE167186_pilot/pilot_seurat.rds")
cat("\ndone\n")
