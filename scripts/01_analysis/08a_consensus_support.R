# ============================================================
# 08a: consensus peak sample-support
#   how many of the 10 samples called a peak overlapping each
#   consensus region -> used to remove regions with missingness bias
# ============================================================
.libPaths(c("C:/Users/wang/R/win-library/4.5", .libPaths()))
suppressMessages({library(GenomicRanges); library(IRanges)
  library(BSgenome.Hsapiens.UCSC.hg38); library(GenomeInfoDb)})

proj <- "D:/workbuddy/outputs/2026-10-25-20-23-31"
data.dir <- file.path(proj, "data/GSE268953")
out.dir <- file.path(proj, "results/GSE268953")

cons <- readRDS(file.path(out.dir, "consensus_peaks.rds"))
cat("consensus regions:", length(cons), "\n")

samples <- c(paste0("Aged_", 1:5), paste0("Young_", 1:5))
std.chr <- paste0("chr", c(1:22, "X", "Y"))
support <- integer(length(cons))
for (s in samples) {
  f <- read.delim(gzfile(file.path(data.dir, paste0("GSE268953_", s, "_features.tsv.gz"))),
                  header = FALSE, stringsAsFactors = FALSE)
  pr <- f[f$V3 == "Peaks", ]
  g <- GRanges(pr$V4, IRanges(pr$V5, pr$V6))
  g <- g[as.character(seqnames(g)) %in% std.chr]
  support <- support + as.integer(countOverlaps(cons, g) > 0)
  cat(" ", s, "peaks:", length(g), "\n")
}

saveRDS(support, file.path(out.dir, "consensus_support.rds"))
cat("\nsupport distribution:\n"); print(table(support))
cat("regions called in all 10 samples:", sum(support == 10), "\n")
cat("regions called in >=9 samples:", sum(support >= 9), "\n")
