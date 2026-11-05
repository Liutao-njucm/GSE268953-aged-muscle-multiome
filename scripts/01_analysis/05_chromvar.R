# ============================================================
# Phase 2A - Step 1: chromVAR TF motif activity on ATAC
# GSE268953 skeletal muscle snMultiome (Li et al. Nat Commun 2025 re-analysis)
# Input : results/GSE268953/merged_annotated.rds (from Phase 1)
# Output: motif deviation matrix + LSI embedding saved to results/GSE268953/
# ============================================================
.libPaths(c("C:/Users/wang/R/win-library/4.5", .libPaths()))
suppressMessages({
  library(Seurat); library(Signac); library(chromVAR); library(motifmatchr)
  library(JASPAR2024); library(TFBSTools)
  library(BSgenome.Hsapiens.UCSC.hg38)
})
setwd("D:/workbuddy/outputs/2026-10-25-20-23-31")
out.dir <- "results/GSE268953"

merged <- readRDS(file.path(out.dir, "merged_multiome.rds"))
cat("loaded:", ncol(merged), "nuclei | assays:", paste(Assays(merged), collapse=","), "\n")
DefaultAssay(merged) <- "ATAC"

# ---- 1. keep standard chromosomes only ----
peaks <- granges(merged[["ATAC"]])
keep <- as.character(seqnames(peaks)) %in% paste0("chr", c(1:22, "X", "Y"))
cat("peaks total:", length(peaks), "| keep:", sum(keep), "\n")
merged[["ATAC"]] <- subset(merged[["ATAC"]], features = which(keep))

# ---- 2. motif matching (JASPAR 2024 CORE vertebrates) ----
pfm <- getMatrixSet(JASPAR2024,
                    opts = list(species = 9606, collection = "CORE",
                                all_versions = FALSE))
cat("motifs:", length(pfm), "\n"); flush.console()
merged <- AddMotifs(object = merged, genome = BSgenome.Hsapiens.UCSC.hg38, pfm = pfm)
cat("motifs added\n"); flush.console()

# ---- 3. LSI (TFIDF -> SVD) for ATAC latent space ----
merged <- RunTFIDF(merged)
merged <- RunSVD(merged, n = 30, reduction.name = "lsi")
cat("LSI done\n"); flush.console()

# ---- 4. chromVAR deviations ----
merged <- RunChromVAR(merged, genome = BSgenome.Hsapiens.UCSC.hg38)
cat("chromVAR done\n"); flush.console()

dev <- GetAssayData(merged[["chromvar"]], layer = "chromvar_motion")
saveRDS(as.matrix(dev), file.path(out.dir, "chromvar_deviations.rds"))
cat("deviation matrix:", nrow(dev), "motifs x", ncol(dev), "nuclei\n")

# ---- 5. save metadata (celltype/group/sample for downstream) ----
meta <- data.frame(
  barcode   = colnames(merged),
  celltype  = as.character(merged$celltype),
  group     = as.character(merged$group),
  sample    = as.character(merged$sample),
  stringsAsFactors = FALSE)
write.csv(meta, file.path(out.dir, "cell_metadata.csv"), row.names = FALSE)

saveRDS(merged, file.path(out.dir, "merged_chromvar.rds"))
cat("PHASE2A_STEP1_DONE\n")
