#!/usr/bin/env Rscript
# Real head-to-head benchmark: run CMScaller (Eide et al. 2017) on the exact
# same 511-sample TCGA CRC-CMS cohort matrix used by MLHeatmap, and compare
# CMScaller's calls against the CRCSC gold labels (dataset ground truth).
# Cross-comparison against MLHeatmap's own out-of-fold predictions is done
# afterwards in Python (scripts/compare_cmscaller_vs_mlheatmap.py) once both
# result files exist.
#
# Usage: Rscript scripts/run_cmscaller_benchmark.R [cohort_dir] [out_dir]

suppressMessages(library(CMScaller))

args <- commandArgs(trailingOnly = TRUE)
cohort_dir <- if (length(args) >= 1) args[[1]] else file.path("generated", "crc_cms_public")
out_dir <- if (length(args) >= 2) args[[2]] else file.path("generated", "analysis")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

counts_path <- file.path(cohort_dir, "tcga_crc_cms_gold_counts.tsv.gz")
stopifnot(file.exists(counts_path))

cat("Loading cohort:", counts_path, "\n")
counts <- read.delim(gzfile(counts_path), row.names = 1, check.names = FALSE)
emat <- as.matrix(counts)
mode(emat) <- "numeric"
cat("emat dims:", nrow(emat), "genes x", ncol(emat), "samples\n")

sample_names <- colnames(emat)
gold_labels <- sub("__.*$", "", sample_names)
cat("Gold label distribution:\n")
print(table(gold_labels))

cat("\nRunning CMScaller (RNAseq=TRUE, rowNames='symbol', seed=42)...\n")
res <- CMScaller(emat, rowNames = "symbol", RNAseq = TRUE, seed = 42, doPlot = FALSE, verbose = TRUE)

pred <- as.character(res$prediction)
names(pred) <- sample_names
n_na <- sum(is.na(pred))
cat("\nCMScaller produced", sum(!is.na(pred)), "confident calls;", n_na, "NA (below FDR threshold / ambiguous).\n")

comparison <- data.frame(
  sample = sample_names,
  gold_label = gold_labels,
  cmscaller_prediction = pred,
  cmscaller_p_value = res$p.value,
  cmscaller_FDR = res$FDR,
  stringsAsFactors = FALSE
)

out_csv <- file.path(out_dir, "cmscaller_predictions.csv")
write.csv(comparison, out_csv, row.names = FALSE)
cat("Wrote", out_csv, "\n")

confident <- comparison[!is.na(comparison$cmscaller_prediction), ]
agreement <- confident$gold_label == confident$cmscaller_prediction
overall_agreement <- mean(agreement)
cat(sprintf(
  "\nCMScaller vs CRCSC gold labels: %d/%d confident calls agree (%.1f%%); %d/%d of all samples called confidently (%.1f%%)\n",
  sum(agreement), nrow(confident), 100 * overall_agreement,
  nrow(confident), nrow(comparison), 100 * nrow(confident) / nrow(comparison)
))

confusion <- table(gold = confident$gold_label, cmscaller = confident$cmscaller_prediction)
cat("\nConfusion matrix (gold rows x CMScaller columns, confident calls only):\n")
print(confusion)

summary_list <- list(
  n_samples = nrow(comparison),
  n_confident_calls = nrow(confident),
  n_na_calls = n_na,
  overall_agreement_confident = overall_agreement,
  overall_agreement_all = sum(agreement) / nrow(comparison),
  gold_distribution = as.list(table(gold_labels)),
  cmscaller_distribution = as.list(table(pred, useNA = "ifany")),
  confusion_matrix = as.data.frame.matrix(confusion)
)

library(jsonlite)
writeLines(toJSON(summary_list, auto_unbox = TRUE, pretty = TRUE, na = "null"),
           file.path(out_dir, "cmscaller_summary.json"))
cat("Wrote", file.path(out_dir, "cmscaller_summary.json"), "\n")
