lib <- "C:/Users/wang/R/win-library/4.5"
repos <- c(
  CRAN     = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/",
  BioCsoft = "https://mirrors.westlake.edu.cn/bioconductor/packages/3.22/bioc",
  BioCann  = "https://mirrors.westlake.edu.cn/bioconductor/packages/3.22/data/annotation"
)
options(repos = repos, Ncpus = 4, timeout = 1800)
.libPaths(c(lib, .libPaths()))

ap <- available.packages()
cat("available packages:", nrow(ap), "\n")
cat("DESeq2 available:", "DESeq2" %in% rownames(ap), "\n")

install.packages("DESeq2", lib = lib, dependencies = TRUE)
cat("DESeq2 loadable:", requireNamespace("DESeq2", quietly = TRUE), "\n")
if (requireNamespace("DESeq2", quietly = TRUE))
  cat("version:", as.character(packageVersion("DESeq2")), "\n")
