# 01_load_qc.R — GSE268953 snMultiome: load 10 samples, split GEX/ATAC, per-sample QC
lib <- "C:/Users/wang/R/win-library/4.5"
.libPaths(c(lib, .libPaths()))
suppressMessages({
  library(Seurat); library(Signac); library(GenomicRanges); library(Matrix)
})
setwd("D:/workbuddy/outputs/2026-10-25-20-23-31")
data.dir <- "data/GSE268953"; out.dir <- "results/GSE268953"
dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)

samples <- c(paste0("Young_", 1:5), paste0("Aged_", 1:5))
group.map <- setNames(rep(c("Young", "Aged"), each = 5), samples)

read_multiome <- function(fp) {
  f <- read.delim(gzfile(paste0(fp, "_features.tsv.gz")), header = FALSE,
                  stringsAsFactors = FALSE)
  bc <- readLines(gzfile(paste0(fp, "_barcodes.tsv.gz")))
  m <- Matrix::readMM(gzfile(paste0(fp, "_matrix.mtx.gz")))
  m <- as(m, "dgCMatrix")   # 10x mtx: rows=features, cols=cells (no transpose)
  stopifnot(nrow(m) == nrow(f), ncol(m) == length(bc))
  rownames(m) <- make.unique(f$V2); colnames(m) <- bc
  idx.ge <- f$V3 == "Gene Expression"; idx.pk <- f$V3 == "Peaks"
  list(rna = m[idx.ge, ], peaks = m[idx.pk, ],
       ge = f[idx.ge, ], pk = f[idx.pk, ])
}

qc.list <- list(); obj.list <- list()
for (s in samples) {
  cat("==== loading", s, "====\n"); flush.console()
  fp <- file.path(data.dir, paste0("GSE268953_", s))
  parts <- read_multiome(fp)
  rna <- parts$rna; atac <- parts$peaks
  # peak ranges from features file columns (chr, start, end)
  peaks.gr <- GRanges(parts$pk$V4, IRanges(as.numeric(parts$pk$V5),
                                           as.numeric(parts$pk$V6)))
  pk.names <- paste(parts$pk$V4, parts$pk$V5, parts$pk$V6, sep = "-")
  dimnames(atac) <- list(pk.names, colnames(atac))

  obj <- CreateSeuratObject(counts = rna, project = s,
                            min.cells = 3, min.features = 200)
  obj$sample <- s; obj$group <- group.map[[s]]
  atac <- atac[, colnames(obj)]   # sync cells with RNA object
  chrom <- CreateChromatinAssay(counts = atac, genome = "hg38",
                                min.cells = 5, ranges = peaks.gr)
  obj[["ATAC"]] <- chrom
  DefaultAssay(obj) <- "RNA"

  # gene metadata (biotypes from features file)
  gi <- match(rownames(obj), parts$ge$V2)
  obj$type <- parts$ge$V3[gi]

  # RNA QC
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  obj[["percent.ribo"]] <- PercentageFeatureSet(obj, pattern = "^RP[SL]")
  # ATAC QC
  obj$nCount_ATAC <- colSums(GetAssayData(obj, assay = "ATAC", layer = "counts")[, colnames(obj)])
  obj <- subset(obj,
    subset = nFeature_RNA > 200 & nFeature_RNA < 8000 &
             percent.mt < 10 & nCount_RNA < 50000 &
             nCount_ATAC > 500 & nCount_ATAC < 100000)

  qc.list[[s]] <- data.frame(
    sample = s, group = group.map[[s]],
    n_raw = ncol(obj) + 0, n_keep = ncol(obj),
    median_genes = median(obj$nFeature_RNA), median_counts = median(obj$nCount_RNA),
    median_mt = median(obj$percent.mt), median_atac = median(obj$nCount_ATAC))
  obj.list[[s]] <- obj
  cat(s, "kept", ncol(obj), "nuclei\n"); flush.console()
  rm(rna, atac, chrom); invisible(gc())
}

qc.df <- do.call(rbind, qc.list)
write.csv(qc.df, file.path(out.dir, "per_sample_qc.csv"), row.names = FALSE)
saveRDS(obj.list, file.path(out.dir, "samples_qc_list.rds"))
print(qc.df)
cat("PHASE1A DONE\n")
