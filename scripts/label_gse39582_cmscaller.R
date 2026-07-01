#!/usr/bin/env Rscript
# Derive CMS labels for the GSE39582 tumor cohort using CMScaller, since the
# original CRCSC network-consensus label source for this cohort
# (Sage-Bionetworks/crcsc) is no longer reachable (as of this revision).
# Microarray mode (RNAseq=FALSE) since values are already log2 RMA-normalized.

suppressMessages(library(CMScaller))

args <- commandArgs(trailingOnly = TRUE)
expr_path <- if (length(args) >= 1) args[[1]] else file.path("generated", "gse39582", "gse39582_tumor_gene_expr.tsv.gz")
out_dir <- if (length(args) >= 2) args[[2]] else file.path("generated", "gse39582")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("Loading GSE39582 tumor expression:", expr_path, "\n")
expr <- read.delim(gzfile(expr_path), row.names = 1, check.names = FALSE)
emat <- as.matrix(expr)
mode(emat) <- "numeric"
cat("emat dims:", nrow(emat), "genes x", ncol(emat), "samples\n")

cat("\nRunning CMScaller (RNAseq=FALSE, rowNames='symbol', seed=42)...\n")
res <- CMScaller(emat, rowNames = "symbol", RNAseq = FALSE, seed = 42, doPlot = FALSE, verbose = TRUE)

pred <- as.character(res$prediction)
names(pred) <- colnames(emat)
n_na <- sum(is.na(pred))
cat("\nCMScaller produced", sum(!is.na(pred)), "confident calls;", n_na, "NA.\n")

out <- data.frame(
  sample = colnames(emat),
  cms_label = pred,
  p_value = res$p.value,
  FDR = res$FDR,
  stringsAsFactors = FALSE
)
out_path <- file.path(out_dir, "gse39582_cmscaller_labels.csv")
write.csv(out, out_path, row.names = FALSE)
cat("Wrote", out_path, "\n")
print(table(pred, useNA = "ifany"))
