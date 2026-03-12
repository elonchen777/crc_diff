# ============================================================
# 顶刊级：三组菌-菌 Spearman 相关性热图（ComplexHeatmap 风格）
# ============================================================

library(ComplexHeatmap); library(circlize); library(data.table)
library(Hmisc); library(grid)

merged_data <- fread("dataset/merged_dataset_processed.csv",
                     stringsAsFactors = FALSE, data.table = FALSE)
output_dir <- "results/R_plots/heatmap_correlation_plot"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

merged_data$group <- with(merged_data, ifelse(
  crc_label == 0, "control",
  ifelse(differentiation == 1, "CRC_poor_diff", "CRC_well_diff")
))

target_features <- c(
  "Peptostreptococcus_stomatis", "Porphyromonas_gingivalis",
  "Prevotella_intermedia",       "Fusobacterium_periodonticum",
  "Campylobacter_rectus",        "Faecalibacterium_prausnitzii",
  "Roseburia_intestinalis",      "Eubacterium_rectale",
  "Coprococcus_comes",           "Ruminococcus_lactaris"
)

match_species_to_cols <- function(target_list, available_cols) {
  matched <- character(0)
  missing <- character(0)
  
  for (species_name in target_list) {
    pattern <- tolower(species_name)
    matches <- available_cols[grepl(pattern, tolower(available_cols))]
    
    if (length(matches) > 0) {
      matched <- c(matched, matches[1])
    } else {
      missing <- c(missing, species_name)
    }
  }
  
  if (length(missing) > 0) {
    stop("Species not found in dataset: ", paste(missing, collapse = ", "))
  }
  
  message("Matched species: ", paste(gsub("^tax_s__", "", matched), collapse = ", "))
  return(matched)
}

all_tax_cols <- colnames(merged_data)[grepl("^tax_s__", colnames(merged_data))]
available_cols <- match_species_to_cols(target_features, all_tax_cols)

simplify_name <- function(name) {
  x <- gsub("^tax_s__", "", name)
  parts <- strsplit(x, "_")[[1]]
  if (length(parts) >= 2) paste0(substr(parts[1], 1, 1), ". ", parts[2])
  else x
}
sp_labels <- sapply(available_cols, simplify_name)

make_corr_heatmap <- function(group_name, title_label, col_title_color) {
  df_g <- merged_data[merged_data$group == group_name, available_cols]
  colnames(df_g) <- sp_labels
  rc <- rcorr(as.matrix(df_g), type = "spearman")
  r <- rc$r; p <- rc$P
  
  # 星号矩阵
  sig_mat <- matrix("", nrow = nrow(r), ncol = ncol(r))
  sig_mat[p < 0.001] <- "***"
  sig_mat[p >= 0.001 & p < 0.01]  <- "**"
  sig_mat[p >= 0.01  & p < 0.05]  <- "*"
  
  col_fun <- colorRamp2(c(-1, -0.5, 0, 0.5, 1),
                        c("#4575B4", "#91BFDB", "white", "#FC8D59", "#D73027"))
  
  Heatmap(r,
    name = "Spearman r",
    col  = col_fun,
    cell_fun = function(j, i, x, y, width, height, fill) {
      grid.text(sig_mat[i, j], x, y, gp = gpar(fontsize = 9, col = "black"))
    },
    row_names_gp    = gpar(fontsize = 9, fontface = "italic"),
    column_names_gp = gpar(fontsize = 9, fontface = "italic"),
    column_names_rot = 45,
    show_row_dend    = TRUE,
    show_column_dend = TRUE,
    column_title     = title_label,
    column_title_gp  = gpar(fontsize = 12, fontface = "bold",
                             col = col_title_color),
    width  = unit(6, "cm"),
    height = unit(6, "cm")
  )
}

ht_ctrl      <- make_corr_heatmap("control",       "Control",       "#2E86AB")
ht_well      <- make_corr_heatmap("CRC_well_diff", "CRC well-diff", "#A23B72")
ht_poor      <- make_corr_heatmap("CRC_poor_diff", "CRC poor-diff", "#F18F01")


png(file.path(output_dir, "heatmap_3groups.png"),
    width = 18, height = 7, units = "in", res = 300)
draw(ht_ctrl + ht_well + ht_poor,
     column_title = "Microbial Correlation Heatmaps across CRC Subtypes",
     column_title_gp = gpar(fontsize = 15, fontface = "bold"))
dev.off()

pdf(file.path(output_dir, "heatmap_3groups.pdf"),
    width = 18, height = 7)
draw(ht_ctrl + ht_well + ht_poor,
     column_title = "Microbial Correlation Heatmaps across CRC Subtypes",
     column_title_gp = gpar(fontsize = 15, fontface = "bold"))
dev.off()
message("Done: ", output_dir)