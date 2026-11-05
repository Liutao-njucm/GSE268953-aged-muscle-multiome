# ============================================================
# 09b: map JASPAR matrix IDs -> TF names, rebuild aging signature
# ============================================================
.libPaths(c("C:/Users/wang/R/win-library/4.5", .libPaths()))
suppressMessages({library(JASPAR2024); library(TFBSTools); library(motifmatchr)
  library(DBI); library(RSQLite)})

out.dir <- "D:/workbuddy/outputs/2026-10-25-20-23-31/results/GSE268953"

db <- JASPAR2024::JASPAR2024()
con <- DBI::dbConnect(RSQLite::SQLite(), db@db)
pfm <- getMatrixSet(con, opts = list(species = 9606, collection = "CORE", all_versions = FALSE))
DBI::dbDisconnect(con)

map <- data.frame(
  motif_id = names(pfm),
  tf_name  = vapply(pfm, function(x) x@name, character(1)),
  stringsAsFactors = FALSE
)
map$tf_name[is.na(map$tf_name) | map$tf_name == ""] <- map$motif_id[is.na(map$tf_name) | map$tf_name == ""]
cat("motif->TF map:", nrow(map), "matrices |", length(unique(map$tf_name)), "unique TF names\n")
write.csv(map, file.path(out.dir, "motif_id_to_tf_name.csv"), row.names = FALSE)

e <- read.csv(file.path(out.dir, "motif_enrichment_DA_peaks.csv"), stringsAsFactors = FALSE)
e$tf_name <- map$tf_name[match(e$motif, map$motif_id)]
e <- e[!is.na(e$tf_name), ]
write.csv(e, file.path(out.dir, "motif_enrichment_DA_peaks.csv"), row.names = FALSE)
cat("enrichment rows mapped:", nrow(e), "\n")

# ---- rebuild signature at TF level ----
sig.cells <- aggregate(celltype ~ tf_name + direction,
                       data = e[e$fdr < 0.05, ],
                       FUN = function(x) length(unique(x)))
if (nrow(sig.cells) == 0) {
  cat("no motif passes FDR<0.05; using top odds ratio\n")
  agg <- aggregate(odds_ratio ~ tf_name + direction, data = e, FUN = mean)
  agg <- agg[order(-agg$odds_ratio), ]
  agg$n_celltypes_sig <- 0
} else {
  agg <- aggregate(cbind(odds_ratio, p) ~ tf_name + direction, data = e, FUN = mean)
  agg <- merge(agg, sig.cells, by = c("tf_name", "direction"), all.x = TRUE)
  names(agg)[names(agg) == "celltype"] <- "n_celltypes_sig"
  agg$n_celltypes_sig[is.na(agg$n_celltypes_sig)] <- 0
  agg <- agg[order(-agg$odds_ratio), ]
}
cat("TF-level enrichment summary:\n")
print(head(agg, 15), row.names = FALSE)

agg2 <- agg[order(-agg$n_celltypes_sig, -agg$odds_ratio), ]
up <- head(unique(agg2$tf_name[agg2$direction == "up"]), 40)
dn <- head(unique(agg2$tf_name[agg2$direction == "down"]), 40)
sig.df <- data.frame(motif = c(up, dn),
                     direction = c(rep("up_in_Aged", length(up)),
                                   rep("down_in_Aged", length(dn))),
                     stringsAsFactors = FALSE)
write.csv(sig.df, file.path(out.dir, "aging_TF_signature.csv"), row.names = FALSE)
cat("\nsignature saved: up =", length(up), "| down =", length(dn), "\n")
cat("up TFs:", paste(up, collapse = ", "), "\n")
cat("down TFs:", paste(dn, collapse = ", "), "\n")
