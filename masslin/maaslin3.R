#!/usr/bin/env Rscript
# MaasLin3 analysis script
# Usage:
#   Rscript maaslin3.R

outdir <- "results/maaslin3"

if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

# load MaasLin3
if (!requireNamespace("maaslin3", quietly = TRUE)) {
	stop("Please install maaslin3 in R: install.packages('BiocManager'); BiocManager::install('biobakery/maaslin3')")
}
message("Loading maaslin3 package...")
library(maaslin3)
library("data.table")

# read data
message("Reading metadata")
meta <- fread("dataset/maaslin_metadata.tsv", sep="\t", header=TRUE, data.table=FALSE)
# 将第一列设为行名
rownames(meta) <- meta[[1]]
meta <- meta[, -1, drop = FALSE]

message("Reading features")
features <- fread("dataset/maaslin_features.tsv", sep="\t", header=TRUE, data.table=FALSE)
# 将第一列设为行名
rownames(features) <- features[[1]]
features <- features[, -1, drop = FALSE]

# Ensure samples match
common <- intersect(rownames(features), rownames(meta))
if (length(common) == 0) stop("No overlapping samples between features and metadata")
message(paste("共同样本数:", length(common)))

features <- features[common, , drop = FALSE]
meta <- meta[common, , drop = FALSE]

# Convert categorical variables to factors
message("Preparing covariates...")
# 使用正确的列名（从 merged_dataset 转换而来）
meta$differentiation <- factor(meta$differentiation, labels = c("well", "poor"))
meta$gender_label <- factor(meta$gender_label, labels = c("female", "male"))
meta$tnm_stage <- factor(meta$tnm_stage, levels = c(0, 1, 2, 3, 4))

# Run MaasLin3 analysis
message("Running MaasLin3 analysis...")
fit_out <- maaslin3(
	input_data = as.data.frame(features),
	input_metadata = as.data.frame(meta),
	output = outdir,
	formula = '~ differentiation + age + gender_label + tnm_stage',
	normalization = 'TSS',
	transform = 'LOG',
	cores = 4
)

message("MaasLin3 analysis completed!")
message(paste("Results saved to:", outdir))
