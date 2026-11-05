# ============================================================
# Phase 2A - Step 1 (v3): per-sample chromVAR, peak-capped, incremental
# - fixed top 60,000 peaks by accessibility per sample (equal peak count
#   across samples -> comparable deviations, and ~4x faster)
# - each sample saved immediately to atac_tmp/dev_<sample>.rds
# Log: ../../scripts/05b_log.txt
# ============================================================
.libPaths(c("C:/Users/wang/R/win-library/4.5", .libPaths()))
suppressMessages({
  library(Seurat); library(Signac); library(chromVAR); library(motifmatchr)
  library(JASPAR2024); library(TFBSTools); library(BSgenome.Hsapiens.UCSC.hg38)
  library(DBI); library(RSQLite)
})
setwd("D:/workbuddy/outputs/2026-10-25-20-23-31/results/GSE268953")
dir.create("atac_tmp", showWarnings = FALSE)
logf <- file("../../scripts/05b_log.txt", open = "wt")
sink(logf, split = TRUE); sink(logf, type = "message")
t0 <- Sys.time(); cat("start:", format(t0), "\n"); flush.console()

PEAK_CAP <- 60000
meta <- read.csv("cell_metadata.csv", stringsAsFactors = FALSE)
jaspar.con <- DBI::dbConnect(RSQLite::SQLite(), JASPAR2024::JASPAR2024()@db)
pfm <- getMatrixSet(jaspar.con,
                    opts = list(species = 9606, collection = "CORE", all_versions = FALSE))
bgs <- BSgenome.Hsapiens.UCSC.hg38
keep.chr <- paste0("chr", c(1:22, "X", "Y"))
cat("motifs:", length(pfm), "| peak cap:", PEAK_CAP, "\n"); flush.console()

qcl <- readRDS("samples_qc_list.rds")
for (s in names(qcl)) {
  out.f <- file.path("atac_tmp", paste0("dev_", s, ".rds"))
  if (file.exists(out.f)) { cat(s, "already done, skip\n"); next }
  o0 <- qcl[[s]]
  sel <- meta$barcode[meta$sample == s]
  keep <- intersect(gsub(paste0("^", s, "_"), "", sel), colnames(o0))
  cnt <- GetAssayData(o0[["ATAC"]], layer = "counts")[, keep]
  pk <- Signac::granges(o0[["ATAC"]])
  rs <- Matrix::rowSums(cnt)
  ok <- which((as.character(seqnames(pk)) %in% keep.chr) & rs >= 10)
  ord <- ok[order(rs[ok], decreasing = TRUE)]
  sel.f <- sort(head(ord, PEAK_CAP))
  cnt <- cnt[sel.f, ]; pk <- pk[sel.f]
  cat(s, "cells:", ncol(cnt), "| peaks:", nrow(cnt), "\n"); flush.console()

  bc.g <- setNames(meta$group[meta$sample == s],
                   gsub(paste0("^", s, "_"), "", meta$barcode[meta$sample == s]))
  grp <- factor(bc.g[keep], levels = c("Young", "Aged"))
  pb <- sapply(levels(grp), function(g) Matrix::rowSums(cnt[, grp == g, drop = FALSE]))
  colnames(pb) <- paste(s, colnames(pb), sep = "_")

  motif.ix <- motifmatchr::matchMotifs(pfm, pk, genome = bgs, out = "matches")
  cat(s, "motif matching done, computing deviations...\n"); flush.console()
  dvo <- chromVAR::computeDeviations(object = cnt, motifs = motif.ix)
  dv <- as.matrix(SummarizedExperiment::assay(dvo, "deviations"))
  colnames(dv) <- paste(s, keep, sep = "_")

  saveRDS(list(dev = dv, pb = pb, gr = pk), out.f)
  cat(s, "DONE:", nrow(dv), "motifs x", ncol(dv), "cells | elapsed",
      round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n"); flush.console()
  rm(o0, cnt, pk, pb, dv, dvo, motif.ix); invisible(gc())
}
cat("ALL_SAMPLES_DONE in", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")
sink(type = "message"); sink(); close(logf)
