#!/usr/bin/env Rscript
# MaAsLin2 analysis script for KO genes
# Usage:
#   Rscript masslin/maaslin_ko.R --features dataset/maaslin_ko_features.tsv --meta dataset/maaslin_metadata.tsv --outdir results/maaslin2_ko

outdir <- "results/maaslin2_ko"
features_path <- "dataset/maaslin/maaslin_ko_features.tsv"
meta_path <- "dataset/maaslin/maaslin_metadata.tsv"

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  for (i in seq(1, length(args), by = 2)) {
    if (i + 1 <= length(args)) {
      key <- args[i]
      val <- args[i + 1]
      if (key == "--outdir") outdir <- val
      if (key == "--features") features_path <- val
      if (key == "--meta") meta_path <- val
    }
  }
}

if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

if (!requireNamespace("Maaslin2", quietly = TRUE)) {
  stop("Please install Maaslin2 in R: install.packages('BiocManager'); BiocManager::install('Maaslin2')")
}
if (!requireNamespace("data.table", quietly = TRUE)) {
  stop("Please install data.table in R: install.packages('data.table')")
}

message("Loading packages...")
library(Maaslin2)
library(data.table)

read_with_sample_id <- function(path, dataset_name) {
  df <- fread(path, sep = "\t", header = TRUE, data.table = FALSE, check.names = FALSE)

  if (ncol(df) < 2) {
    stop(paste(dataset_name, "has fewer than 2 columns:"), path)
  }

  first_col_name <- colnames(df)[1]
  sample_col_candidates <- c("SAMPLE_ID", "sample_id", "sample", "SampleID", "sampleID")
  sample_col <- if (first_col_name %in% sample_col_candidates) first_col_name else colnames(df)[1]

  sample_ids <- trimws(as.character(df[[sample_col]]))
  keep <- !is.na(sample_ids) & sample_ids != ""
  if (sum(!keep) > 0) {
    message(sprintf("%s: dropped %d rows with empty sample ID", dataset_name, sum(!keep)))
  }

  df <- df[keep, , drop = FALSE]
  sample_ids <- sample_ids[keep]

  dup <- duplicated(sample_ids)
  if (sum(dup) > 0) {
    message(sprintf("%s: dropped %d duplicated sample IDs", dataset_name, sum(dup)))
    df <- df[!dup, , drop = FALSE]
    sample_ids <- sample_ids[!dup]
  }

  rownames(df) <- sample_ids
  df[[sample_col]] <- NULL

  return(df)
}

message("Reading metadata...")
meta <- read_with_sample_id(meta_path, "metadata")

message("Reading KO features...")
features <- read_with_sample_id(features_path, "features")

common <- intersect(rownames(features), rownames(meta))
if (length(common) == 0) stop("No overlapping samples between KO features and metadata")

features <- features[common, , drop = FALSE]
meta <- meta[common, , drop = FALSE]

features[] <- lapply(features, function(x) suppressWarnings(as.numeric(x)))
features[is.na(features)] <- 0

message("Preparing covariates...")
covariates <- c("age", "gender_label", "diff_stage")
missing_covariates <- covariates[!covariates %in% colnames(meta)]
if (length(missing_covariates) > 0) {
  stop(paste("Missing covariates in metadata:", paste(missing_covariates, collapse = ", ")))
}

if ("diff_stage" %in% colnames(meta)) {
  meta$diff_stage <- as.factor(meta$diff_stage)
}

message(sprintf("Running KO MaAsLin2 analysis with %d samples and %d KO features...", nrow(features), ncol(features)))
fit <- Maaslin2(
  input_data = as.data.frame(features),
  input_metadata = as.data.frame(meta),
  output = outdir,
  fixed_effects = covariates,
  reference = c("diff_stage,0"),
  transform = "LOG",
  normalization = "TSS",
  correction = "BH",
  min_prevalence = 0.1,
  plot_heatmap = FALSE,
  plot_scatter = FALSE
)

message("KO MaAsLin2 analysis completed.")
