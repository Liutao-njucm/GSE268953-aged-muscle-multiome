# ============================================================
# 06: merge per-sample chromVAR deviations into one matrix
# ============================================================
.libPaths(c("C:/Users/wang/R/win-library/4.5", .libPaths()))
suppressMessages(library(Matrix))

proj <- "D:/workbuddy/outputs/2026-10-25-20-23-31"
out.dir <- file.path(proj, "results/GSE268953")
tmp.dir <- file.path(out.dir, "atac_tmp")

samples <- c(paste0("Aged_", 1:5), paste0("Young_", 1:5))
files <- file.path(tmp.dir, paste0("dev_", samples, ".rds"))
missing <- samples[!file.exists(files)]
if (length(missing)) stop("missing: ", paste(missing, collapse = ","))

cat("merging", length(samples), "samples ...\n")
z.list <- list(); dev.list <- list(); motif.ids <- NULL
for (i in seq_along(samples)) {
  o <- readRDS(files[i])
  if (is.null(motif.ids)) motif.ids <- o$motif
  stopifnot(identical(o$motif, motif.ids))
  z.list[[samples[i]]] <- o$z
  dev.list[[samples[i]]] <- o$dev
  cat(" ", samples[i], dim(o$z), "\n")
  rm(o); invisible(gc())
}

z.all <- do.call(cbind, z.list)
dev.all <- do.call(cbind, dev.list)
rownames(z.all) <- rownames(dev.all) <- motif.ids
cat("merged:", nrow(z.all), "motifs x", ncol(z.all), "cells\n")

meta <- read.csv(file.path(out.dir, "cell_metadata.csv"), stringsAsFactors = FALSE)
idx <- match(colnames(z.all), meta$barcode)
cat("cells matched to metadata:", sum(!is.na(idx)), "/", ncol(z.all), "\n")
stopifnot(!any(is.na(idx)))
meta <- meta[idx, ]

saveRDS(list(z = z.all, dev = dev.all, motif = motif.ids, meta = meta),
        file.path(out.dir, "motif_activity_merged.rds"))
cat("saved motif_activity_merged.rds\n")
