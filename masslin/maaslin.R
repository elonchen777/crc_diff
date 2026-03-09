#!/usr/bin/env Rscript
# MaAsLin2 analysis script
# Usage:
#   Rscript maaslin.r --features dataset/maaslin_features.tsv --meta dataset/maaslin_metadata.tsv --outdir results

outdir <- "results/maaslin2"

if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

# load MaAsLin2
if (!requireNamespace("Maaslin2", quietly = TRUE)) {
	stop("Please install Maaslin2 in R: install.packages('BiocManager'); BiocManager::install('Maaslin2') or devtools::install_github('biobakery/biobakery_workflows') etc.")
}
message("Loading Maaslin2 package...")
library(Maaslin2)
library("data.table")

# read data
message("Reading metadata")
# meta <- read.table("dataset/maaslin_metadata.tsv", sep="\t", header=TRUE, row.names=1, check.names = FALSE, stringsAsFactors = FALSE)
meta <- fread("dataset/maaslin_metadata.tsv", sep="\t", header=TRUE, data.table=FALSE)

message("Reading features")
# features <- read.table("dataset/maaslin_features.tsv", sep="\t", header=TRUE, row.names=1, check.names = FALSE)
features <- fread("dataset/maaslin_features.tsv", sep="\t", header=TRUE, data.table=FALSE)

# Ensure samples match
common <- intersect(rownames(features), rownames(meta))
if (length(common) == 0) stop("No overlapping samples between features and metadata")
features <- features[common, , drop = FALSE]
meta <- meta[common, , drop = FALSE]

# Example: choose covariates from metadata (use all columns except those you want to exclude)
message("Preparing covariates...")
covariates <- c("age", "gender_label","diff_stage")

message("Running MaAsLin2 analysis...")
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

