.libPaths(c("C:/Users/wang/R/win-library/4.5", .libPaths()))
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
options(BioC_mirror = "https://mirrors.westlake.edu.cn/bioconductor")
options(timeout = 600)
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
need <- c("org.Hs.eg.db", "TxDb.Hsapiens.UCSC.hg38.knownGene")
for (p in need) {
  if (!requireNamespace(p, quietly = TRUE)) {
    cat("installing", p, "\n"); flush.console()
    try(BiocManager::install(p, ask = FALSE, update = FALSE))
  }
  cat(p, "->", requireNamespace(p, quietly = TRUE), "\n"); flush.console()
}
