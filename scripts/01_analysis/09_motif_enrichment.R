# ============================================================
# 09: motif enrichment on differential accessibility peaks
#   foreground: DA peaks (Aged vs Young, per cell type, up/down)
#   background: all tested consensus peaks
#   -> defines the aging TF signature with far more power than
#      per-cell motif activity comparison (5v5 donor design)
# ============================================================
.libPaths(c("C:/Users/wang/R/win-library/4.5", .libPaths()))
suppressMessages({
  library(Matrix); library(GenomicRanges); library(SummarizedExperiment)
  library(motifmatchr); library(JASPAR2024); library(TFBSTools)
  library(BSgenome.Hsapiens.UCSC.hg38); library(DBI); library(RSQLite)
})

proj <- "D:/workbuddy/outputs/2026-10-25-20-23-31"
out.dir <- file.path(proj, "results/GSE268953")

cons <- readRDS(file.path(out.dir, "consensus_peaks.rds"))
cat("consensus peaks:", length(cons), "\n")
cons.id <- paste0(as.character(seqnames(cons)), "-", start(cons), "-", end(cons))

# ---- motif x peak incidence on consensus peaks (cached) ----
mm.file <- file.path(out.dir, "consensus_motif_matches.rds")
if (file.exists(mm.file)) {
  MI <- readRDS(mm.file); cat("motif matches loaded from cache:", dim(MI), "\n")
} else {
  bgs <- BSgenome.Hsapiens.UCSC.hg38
  db <- JASPAR2024::JASPAR2024(); con <- DBI::dbConnect(RSQLite::SQLite(), db@db)
  pfm <- getMatrixSet(con, opts = list(species = 9606, collection = "CORE", all_versions = FALSE))
  DBI::dbDisconnect(con)
  cat("motifs:", length(pfm), "| scanning", length(cons), "peaks ...\n")
  t0 <- Sys.time()
  mi <- matchMotifs(pfm, cons, genome = bgs, out = "matches")
  MI <- as(assay(mi, "motifMatches"), "CsparseMatrix")   # peaks x motifs
  rownames(MI) <- cons.id; colnames(MI) <- colnames(mi)
  cat("scan done in", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")
  saveRDS(MI, mm.file)
}

da <- read.csv(file.path(out.dir, "DA_peaks_Aged_vs_Young_by_celltype.csv"),
               stringsAsFactors = FALSE)
da <- da[!is.na(da$padj), ]
cat("DA table rows:", nrow(da), "\n")

tested <- unique(da$peak)
cat("tested peaks:", length(tested), "\n")
mi.t <- MI[as.integer(tested), , drop = FALSE]
motif.names <- colnames(mi.t)
motif.size <- Matrix::colSums(mi.t)
cat("motifs with >=10 peaks in tested set:", sum(motif.size >= 10), "\n")

fisher_motifs <- function(fg, bg.pool = tested, mi.mat = mi.t) {
  fg <- intersect(fg, bg.pool)
  if (length(fg) < 30) return(NULL)
  fg.idx <- match(fg, bg.pool)
  in.fg <- rep(FALSE, nrow(mi.mat)); in.fg[fg.idx] <- TRUE
  fg.hits <- Matrix::colSums(mi.mat[in.fg, , drop = FALSE])
  bg.hits <- Matrix::colSums(mi.mat[!in.fg, , drop = FALSE])
  n.fg <- sum(in.fg); n.bg <- sum(!in.fg)
  or <- (fg.hits / pmax(n.fg - fg.hits, 1)) / (bg.hits / pmax(n.bg - bg.hits, 1))
  p <- pbinom(fg.hits - 1, n.fg, (fg.hits + bg.hits) / (n.fg + n.bg), lower.tail = FALSE)
  data.frame(motif = motif.names, n_motif = motif.size,
             hits_fg = fg.hits, fg_size = n.fg,
             expected = (fg.hits + bg.hits) / (n.fg + n.bg) * n.fg,
             odds_ratio = as.vector(or), p = as.vector(p),
             fdr = p.adjust(as.vector(p), "BH"), stringsAsFactors = FALSE)
}

# ---- per cell type, up and down ----
res <- list()
for (ct in unique(da$celltype)) {
  sub <- da[da$celltype == ct, ]
  for (dir in c("up", "down")) {
    fg <- sub$peak[sub$padj < 0.05 & sub$log2FoldChange > 0.58]
    if (dir == "down") fg <- sub$peak[sub$padj < 0.05 & sub$log2FoldChange < -0.58]
    if (length(fg) < 30) { cat("  ", ct, dir, ": only", length(fg), "peaks, skip\n"); next }
    r <- fisher_motifs(fg)
    if (is.null(r)) next
    r$celltype <- ct; r$direction <- dir
    res[[paste(ct, dir)]] <- r
    cat("  ", ct, dir, ": ", length(fg), "peaks | motifs FDR<0.05:",
        sum(r$fdr < 0.05), "\n")
  }
}
all.res <- do.call(rbind, res)
all.res <- all.res[order(all.res$celltype, all.res$p), ]
write.csv(all.res, file.path(out.dir, "motif_enrichment_DA_peaks.csv"), row.names = FALSE)

# ---- global signature: motifs enriched consistently across cell types ----
agg <- aggregate(cbind(odds_ratio, p) ~ motif + direction, data = all.res, FUN = mean)
cnt <- aggregate(celltype ~ motif + direction, data = all.res[all.res$fdr < 0.05, ],
                 FUN = function(x) length(unique(x)))
agg <- merge(agg, cnt, by = c("motif", "direction"), all.x = TRUE)
names(agg)[names(agg) == "celltype"] <- "n_celltypes_sig"
agg$n_celltypes_sig[is.na(agg$n_celltypes_sig)] <- 0
agg$fdr <- p.adjust(agg$p, "BH")
agg <- agg[order(-agg$odds_ratio), ]

sig <- agg[agg$n_celltypes_sig >= 2, ]
if (nrow(sig) < 10) sig <- head(agg, 40)
up <- head(sig$motif[sig$direction == "up"], 40)
dn <- head(sig$motif[sig$direction == "down"], 40)
write.csv(data.frame(motif = c(up, dn),
                     direction = c(rep("up_in_Aged", length(up)), rep("down_in_Aged", length(dn)))),
          file.path(out.dir, "aging_TF_signature.csv"), row.names = FALSE)
cat("\nsignature: up =", length(up), "| down =", length(dn), "\n")
print(head(agg, 12), row.names = FALSE)
