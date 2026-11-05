# ============================================================
# 05f: ATAC per-sample processing (chromVAR + consensus pseudobulk)
# Direct raw-mtx reading -> avoids loading 1.8GB Seurat objects
# ============================================================
.libPaths(c("C:/Users/wang/R/win-library/4.5", .libPaths()))
suppressMessages({
  library(Matrix); library(GenomicRanges); library(IRanges)
  library(SummarizedExperiment); library(chromVAR); library(motifmatchr)
  library(JASPAR2024); library(TFBSTools); library(BSgenome.Hsapiens.UCSC.hg38)
  library(DBI); library(RSQLite); library(GenomeInfoDb)
})

proj   <- "D:/workbuddy/outputs/2026-10-25-20-23-31"
data.dir <- file.path(proj, "data/GSE268953")
out.dir  <- file.path(proj, "results/GSE268953")
tmp.dir  <- file.path(out.dir, "atac_tmp")
dir.create(tmp.dir, showWarnings = FALSE, recursive = TRUE)
logf <- file(file.path(proj, "scripts/05f_log.txt"), open = "wt")
say <- function(...) { cat(format(Sys.time(), "%H:%M:%S"), "|", ..., "\n");
                       cat(format(Sys.time(), "%H:%M:%S"), "|", ..., "\n", file = logf); flush(logf) }

say("=== START 05f ===")
gc()

# Windows: BiocParallel bplapply deadlocks with default backend -> force serial
suppressMessages(library(BiocParallel))
register(SerialParam(progressbar = FALSE), default = TRUE)
say("BiocParallel backend:", class(bpparam())[1])

# ---- args: which samples to run ----
args <- commandArgs(trailingOnly = TRUE)
all.samples <- c(paste0("Aged_", 1:5), paste0("Young_", 1:5))
todo <- if (length(args) > 0) args else all.samples
say("samples to run:", paste(todo, collapse = ","))

# ---- cell metadata (from Phase 1) ----
meta <- read.csv(file.path(out.dir, "cell_metadata.csv"), stringsAsFactors = FALSE)
meta$bc <- sub("^[A-Za-z]+_[0-9]+_", "", meta$barcode)
say("metadata cells:", nrow(meta), "| celltypes:", length(unique(meta$celltype)))

# ---- motifs (JASPAR2024 human CORE) ----
bgs <- BSgenome.Hsapiens.UCSC.hg38
db  <- JASPAR2024::JASPAR2024()
con <- DBI::dbConnect(RSQLite::SQLite(), db@db)
pfm <- getMatrixSet(con, opts = list(species = 9606, collection = "CORE", all_versions = FALSE))
DBI::dbDisconnect(con)
say("motifs loaded:", length(pfm))

std.chr <- paste0("chr", c(1:22, "X", "Y"))

# ---- read features helper ----
read_feats <- function(s) {
  read.delim(gzfile(file.path(data.dir, paste0("GSE268953_", s, "_features.tsv.gz"))),
             header = FALSE, stringsAsFactors = FALSE)
}

# ---- PASS 1: consensus peak set (reduce of union, supported by >=2 samples) ----
cons.file <- file.path(out.dir, "consensus_peaks.rds")
if (file.exists(cons.file)) {
  cons <- readRDS(cons.file)
  say("consensus loaded from cache:", length(cons), "regions")
} else {
  say("building consensus peaks from all 10 samples ...")
  all.gr <- GRangesList()
  for (s in all.samples) {
    f <- read_feats(s); pr <- f[f$V3 == "Peaks", ]
    g <- GRanges(pr$V4, IRanges(pr$V5, pr$V6))
    g <- keepSeqlevels(g, intersect(std.chr, unique(as.character(seqnames(g)))),
                       pruning.mode = "coarse")
    all.gr[[s]] <- g
    say("  ", s, "peaks:", length(g))
  }
  cons <- reduce(unlist(all.gr))
  cons <- keepSeqlevels(cons, intersect(std.chr, unique(as.character(seqnames(cons)))),
                        pruning.mode = "coarse")
  seqinfo(cons) <- seqinfo(bgs)[seqlevels(cons)]
  saveRDS(cons, cons.file)
  say("consensus regions:", length(cons))
  rm(all.gr); gc()
}

# ---- PASS 2: per sample ----
for (s in todo) {
  dev.file <- file.path(tmp.dir, paste0("dev_", s, ".rds"))
  pb.file  <- file.path(tmp.dir, paste0("pb_", s, ".rds"))
  if (file.exists(dev.file) && file.exists(pb.file)) { say(s, "already done, skip"); next }

  t0 <- Sys.time()
  say("---- ", s, " ----")
  f <- read_feats(s)
  peak.rows <- which(f$V3 == "Peaks")
  pk.all <- GRanges(f$V4[peak.rows], IRanges(f$V5[peak.rows], f$V6[peak.rows]))
  say(s, "features: genes", sum(f$V3 == "Gene Expression"), "| peaks", length(peak.rows))

  bcs <- readLines(gzfile(file.path(data.dir, paste0("GSE268953_", s, "_barcodes.tsv.gz"))))
  say(s, "raw barcodes:", length(bcs))

  m <- Matrix::readMM(gzfile(file.path(data.dir, paste0("GSE268953_", s, "_matrix.mtx.gz"))))
  say(s, "matrix read:", nrow(m), "x", ncol(m),
      "| elapsed", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min")

  # subset to retained barcodes
  keep.bc <- meta$bc[meta$sample == s]
  j <- match(keep.bc, bcs)
  ok <- !is.na(j)
  j <- j[ok]; keep.bc <- keep.bc[ok]
  m <- m[, j, drop = FALSE]
  cnt <- as(m[peak.rows, , drop = FALSE], "CsparseMatrix")
  rm(m); gc()
  colnames(cnt) <- paste(s, keep.bc, sep = "_")
  say(s, "cells used:", ncol(cnt))

  # ---- standard chromosomes ----
  kin <- as.character(seqnames(pk.all)) %in% std.chr
  cnt <- cnt[kin, , drop = FALSE]; pk <- pk.all[kin]
  pk <- keepSeqlevels(pk, intersect(std.chr, unique(as.character(seqnames(pk)))),
                      pruning.mode = "coarse")
  seqinfo(pk) <- seqinfo(bgs)[seqlevels(pk)]

  # ---- pseudobulk on consensus peaks (for differential accessibility) ----
  if (!file.exists(pb.file)) {
    ov <- findOverlaps(pk, cons)
    A <- sparseMatrix(i = subjectHits(ov), j = queryHits(ov), x = 1,
                      dims = c(length(cons), length(pk)))
    agg <- as(A %*% cnt, "CsparseMatrix")   # consensus x cells
    rm(A, ov); gc()
    cts <- unique(meta$celltype[meta$sample == s])
    pbs <- sapply(cts, function(ct) {
      cells <- intersect(colnames(agg),
                         meta$barcode[meta$sample == s & meta$celltype == ct])
      if (length(cells) < 30) return(rep(0, nrow(agg)))
      Matrix::rowSums(agg[, cells, drop = FALSE])
    })
    pbs <- as(pbs, "CsparseMatrix")
    colnames(pbs) <- cts
    saveRDS(list(pb = pbs, celltype_n = sapply(cts, function(ct)
      sum(meta$sample == s & meta$celltype == ct))), pb.file)
    say(s, "pseudobulk saved:", nrow(pbs), "peaks x", ncol(pbs), "celltypes")
    rm(agg, pbs); gc()
  }

  # ---- cap peaks for chromVAR: top 60k by accessibility ----
  rs <- Matrix::rowSums(cnt)
  keep1 <- rs >= 10
  cnt <- cnt[keep1, , drop = FALSE]; pk <- pk[keep1]
  say(s, "peaks with >=10 fragments:", nrow(cnt))
  if (nrow(cnt) > 60000) {
    top <- order(Matrix::rowSums(cnt), decreasing = TRUE)[1:60000]
    top <- sort(top)
    cnt <- cnt[top, , drop = FALSE]; pk <- pk[top]
  }
  cnt <- as(cnt, "CsparseMatrix")
  say(s, "chromVAR peaks:", nrow(cnt), "cells:", ncol(cnt))

  # ---- chromVAR ----
  t1 <- Sys.time()
  mi <- matchMotifs(pfm, pk, genome = bgs, out = "matches")
  say(s, "motif matching done:", round(as.numeric(difftime(Sys.time(), t1, units = "mins")), 1), "min")

  se <- SummarizedExperiment(assays = list(counts = cnt), rowRanges = pk)
  se <- addGCBias(se, genome = bgs)
  bg <- getBackgroundPeaks(se, niterations = 25)
  say(s, "background peaks done:", round(as.numeric(difftime(Sys.time(), t1, units = "mins")), 1), "min")

  dv <- tryCatch(
    computeDeviations(object = cnt, annotations = mi, background_peaks = bg),
    error = function(e) { say(s, "chromVAR error:", conditionMessage(e), "-> sparse fallback"); NULL })

  if (!is.null(dv)) {
    dev <- assay(dv, "deviations"); zz <- assay(dv, "z")
    motif.ids <- as.character(rowData(dv)$name)
  } else {
    M <- as(assay(mi, "matches"), "CsparseMatrix")      # peaks x motifs
    expct <- computeExpectations(cnt)                    # per-peak expectation
    frag <- Matrix::colSums(cnt)
    obs <- as.matrix(Matrix::crossprod(M, cnt))          # motifs x cells
    exp.vec <- as.vector(Matrix::crossprod(M, expct))
    expected <- outer(exp.vec, as.vector(frag))
    dev <- (obs - expected) / expected
    zz  <- dev / sqrt(rowMeans(dev^2) + 1e-8)
    motif.ids <- colnames(M)
  }
  colnames(dev) <- colnames(zz) <- colnames(cnt)
  say(s, "chromVAR DONE: ", nrow(dev), "motifs x", ncol(dev), "cells | stage",
      round(as.numeric(difftime(Sys.time(), t1, units = "mins")), 1), "min | total",
      round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min")

  saveRDS(list(dev = dev, z = zz, motif = motif.ids, peaks = pk), dev.file)
  rm(dev, zz, dv, se, bg, mi, cnt, pk, pk.all); invisible(gc())
  say(s, "saved:", basename(dev.file))
}

say("=== ALL DONE ===")
close(logf)
