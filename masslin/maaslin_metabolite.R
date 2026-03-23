#!/usr/bin/env Rscript
# MaAsLin2 analysis script for metabolites
# Usage:
#   Rscript maaslin_metabolite.R --features dataset/maaslin_metabolite_features.tsv --meta dataset/maaslin_metadata.tsv --outdir results

outdir <- "results/maaslin2_metabolite"

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
meta <- fread("dataset/maaslin_metadata.tsv", sep="\t", header=TRUE, data.table=FALSE)

message("Reading metabolite features")
features <- fread("dataset/maaslin_metabolite_features.tsv", sep="\t", header=TRUE, data.table=FALSE)

# Ensure row names are set
rownames(meta) <- meta$SAMPLE_ID
rownames(features) <- features$SAMPLE_ID
meta$SAMPLE_ID <- NULL
features$SAMPLE_ID <- NULL

# Ensure samples match
common <- intersect(rownames(features), rownames(meta))
if (length(common) == 0) stop("No overlapping samples between features and metadata")
features <- features[common, , drop = FALSE]
meta <- meta[common, , drop = FALSE]

# choose covariates from metadata
message("Preparing covariates...")
covariates <- c("age", "gender_label", "diff_stage")

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
