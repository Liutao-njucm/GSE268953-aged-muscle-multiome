# ============================================================
# 11: annotate consensus peaks to genes (hg38, UCSC refGene)
#   - nearest TSS, distance, region class
#   - gene-body overlap
# ============================================================
.libPaths(c("C:/Users/wang/R/win-library/4.5", .libPaths()))
suppressMessages({library(GenomicRanges); library(jsonlite); library(data.table)})

proj <- "D:/workbuddy/outputs/2026-10-25-20-23-31"
out.dir <- file.path(proj, "results/GSE268953")
ann.dir <- file.path(proj, "data/annotation")

peaks <- readRDS(file.path(out.dir, "consensus_peaks.rds"))
support <- readRDS(file.path(out.dir, "consensus_support.rds"))
cat("peaks:", length(peaks), "\n")

## ---- load refGene (per chromosome json) ----
fs <- list.files(ann.dir, pattern = "^refGene_chr.*\\.json$", full.names = TRUE)
cat("refGene files:", length(fs), "\n")
rg <- rbindlist(lapply(fs, function(f) {
  j <- fromJSON(f, simplifyDataFrame = TRUE)
  if (is.null(j$refGene) || length(j$refGene) == 0) return(NULL)
  as.data.table(j$refGene)[, .(chrom, txStart, txEnd, cdsStart, strand, name, name2)]
}), fill = TRUE)
setnames(rg, "name2", "symbol", skip_absent = TRUE)
rg <- rg[!is.na(symbol) & symbol != "" & !is.na(txStart)]
cat("transcripts:", nrow(rg), "| unique symbols:", uniqueN(rg$symbol), "\n")

## TSS (1 bp)
tss <- rg[, .(chrom, pos = ifelse(strand == "+", txStart + 1L, txEnd), symbol, strand)]
tss[, pos := as.integer(pos)]
tss <- unique(tss[, .(chrom, pos, symbol, strand)])
tss.gr <- GRanges(tss$chrom, IRanges(tss$pos, tss$pos), symbol = tss$symbol)
cat("unique TSS:", length(tss.gr), "\n")

## gene bodies (union per symbol)
gb <- rg[, .(start = min(txStart), end = max(txEnd)), by = .(symbol, chrom)]
gb.gr <- GRanges(gb$chrom, IRanges(gb$start, gb$end), symbol = gb$symbol)

## ---- nearest TSS ----
d <- distanceToNearest(peaks, tss.gr)
ann <- data.table(
  peak_id = paste0(seqnames(peaks), "-", start(peaks), "-", end(peaks)),
  chr = as.character(seqnames(peaks)),
  start = start(peaks), end = end(peaks), width = width(peaks),
  support = as.integer(support),
  symbol = NA_character_, dist_to_tss = NA_integer_,
  tss_pos = NA_integer_, tss_strand = NA_character_
)
qi <- queryHits(d); si <- subjectHits(d)
ann$symbol[qi] <- mcols(tss.gr)$symbol[si]
ann$dist_to_tss[qi] <- mcols(d)$distance
ann$tss_pos[qi] <- start(tss.gr)[si]
ann$tss_strand[qi] <- as.character(strand(tss.gr))[si]

## ---- gene body overlap ----
ov <- findOverlaps(peaks, gb.gr)
ann$in_gene_body <- FALSE
ann$body_symbol <- NA_character_
if (length(ov)) {
  byp <- split(subjectHits(ov), queryHits(ov))
  qidx <- as.integer(names(byp))
  ann$in_gene_body[qidx] <- TRUE
  ann$body_symbol[qidx] <- vapply(byp, function(i) mcols(gb.gr)$symbol[i[1]], character(1))
}

ann[, region_class := cut(abs(dist_to_tss),
                          breaks = c(-1, 2000, 10000, 50000, 100000, Inf),
                          labels = c("Promoter(<=2kb)", "Proximal(2-10kb)",
                                     "Distal(10-50kb)", "Far(50-100kb)", "Intergenic(>100kb)"))]
ann[, linked_gene := ifelse(!is.na(symbol), symbol, body_symbol)]

cat("\n=== region class ===\n"); print(table(ann$region_class, useNA = "ifany"))
cat("\n=== support x region ===\n"); print(table(ann$support, ann$region_class))
cat("\npeaks with a gene within 100kb:", sum(!is.na(ann$linked_gene)), "/", nrow(ann), "\n")
cat("peaks with support==10 and a gene:", ann[support == 10 & !is.na(linked_gene), .N], "\n")

fwrite(ann, file.path(out.dir, "peak_annotation.csv"))
saveRDS(ann, file.path(out.dir, "peak_annotation.rds"))
cat("\nsaved peak_annotation.csv\n")
