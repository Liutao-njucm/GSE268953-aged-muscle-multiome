# 73_inspect_lai_atac.R --- 解包 OMIX004305-05 并检查 per-celltype peak matrix RDS 结构
.libPaths("C:/Users/wang/R/win-library/4.5")
suppressMessages({library(Seurat)})
options(width = 200)

TAR  <- "data/Lai2024_ATAC/OMIX004305-05.tar.gz"
DIR  <- "data/Lai2024_ATAC"

cat("== tar listing ==\n")
tl <- untar(TAR, list = TRUE)
print(tl)

# 只解压小文件（FAPs/MuSC/Immune/Vascular peak matrix；跳过 Myofiber 与 All_CellType）
wanted <- tl[grepl("^(FAPs|MuSC|Immune|Vascular)_snATAC_[0-9]+\\.RDS$", tl)]
cat("wanted:", wanted, "\n")
for (w in wanted) untar(TAR, files = w, exdir = DIR)

f <- file.path(DIR, wanted[1])
cat("\n== inspect:", f, "==\n")
obj <- readRDS(f)
cat("class:", class(obj), "\n")
if (inherits(obj, "Seurat")) {
  print(obj)
  cat("\nassays:", Assays(obj), "\n")
  cat("\nmeta.data columns:\n"); print(colnames(obj@meta.data))
  for (cn in intersect(c("age_group","age","sample","sample_id","donor","sex","orig.ident","celltype","annotation","tech"), colnames(obj@meta.data))) {
    cat("\n--", cn, "--\n"); print(table(obj@meta.data[[cn]], useNA = "ifany"))
  }
  a <- obj[[Assays(obj)[1]]]
  cat("\nassay:", class(a), "features:", nrow(a), "cells:", ncol(a), "\n")
  cat("feature head:\n"); print(head(rownames(a), 20))
  cat("slot names:", slotNames(a), "\n")
} else {
  str(obj, max.level = 2, list.len = 30)
}
