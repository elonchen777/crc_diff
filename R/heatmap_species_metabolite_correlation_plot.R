# ============================================================
# Three-group species-metabolite Spearman correlation heatmaps
# ============================================================

library(ComplexHeatmap)
library(circlize)
library(data.table)
library(Hmisc)
library(grid)

merged_data <- fread(
  "dataset/merged_dataset_processed.csv",
  stringsAsFactors = FALSE,
  data.table = FALSE
)

output_dir <- "results/R_plots/heatmap_species_metabolite_correlation_plot"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

merged_data$group <- with(merged_data, ifelse(
  crc_label == 0, "control",
  ifelse(differentiation == 1, "CRC_poor_diff", "CRC_well_diff")
))

FIXED_SPECIES_LIST <- c(
  "Peptostreptococcus_stomatis",
  "Porphyromonas_gingivalis",
  "Prevotella_intermedia",
  "Fusobacterium_periodonticum",
  "Campylobacter_rectus",
  "Faecalibacterium_prausnitzii",
  "Roseburia_intestinalis",
  "Eubacterium_rectale",
  "Coprococcus_comes",
  "Ruminococcus_lactaris"
)

FIXED_METABOLITES_LIST <- c(
  "SQDG 26:2; SQDG(13:1/13:1)",
  "Cytosine",
  "Perfluorooctanesulfonic acid",
  "Methyl dihydrojasmonate",
  "Pyrogallol-2-O-sulphate",
  "5'-(3',4'-Dihydroxyphenyl)-gamma-valerolactone sulfate",
  "2-Hydroxy-4,7-dimethoxy-2H-1,4-benzoxazin-3(4H)-one",
  "trans-3,5-Dimethoxy-4-hydroxycinnamaldehyde",
  "(R)-3-Hydroxy-5-phenylpentanoic acid",
  "N-Methyl-D-glucamine",
  "Chenodeoxycholic acid sulfate",
  "Creatinine",
  "Lucidenic acid F",
  "Demissidine",
  "Alpha-Hydroxyisobutyric acid",
  "Pyrocatechol",
  "Gentisic acid",
  "D-Galacturonic acid",
  "1,3-Dimethyluric acid",
  "4-Hydroxy-5-(phenyl)-valeric acid-O-sulphate"
)

normalize_text <- function(x) {
  tolower(gsub("[^a-z0-9]", "", x))
}

match_targets_to_cols <- function(target_list, available_cols, remove_prefix = "") {
  display_names <- gsub(remove_prefix, "", available_cols)
  norm_display <- normalize_text(display_names)

  matched_cols <- character(0)
  missing_targets <- character(0)

  for (target in target_list) {
    norm_target <- normalize_text(target)

    exact_idx <- which(norm_display == norm_target)
    if (length(exact_idx) > 0) {
      matched_cols <- c(matched_cols, available_cols[exact_idx[1]])
      next
    }

    partial_idx <- which(grepl(norm_target, norm_display, fixed = TRUE))
    if (length(partial_idx) > 0) {
      matched_cols <- c(matched_cols, available_cols[partial_idx[1]])
      next
    }

    missing_targets <- c(missing_targets, target)
  }

  if (length(missing_targets) > 0) {
    stop("Target features not found: ", paste(missing_targets, collapse = ", "))
  }

  matched_cols
}

species_cols_all <- colnames(merged_data)[grepl("^tax_s__", colnames(merged_data))]
met_cols_all <- colnames(merged_data)[grepl("^met_", colnames(merged_data))]

species_cols <- match_targets_to_cols(
  target_list = FIXED_SPECIES_LIST,
  available_cols = species_cols_all,
  remove_prefix = "^tax_s__"
)

met_cols <- match_targets_to_cols(
  target_list = FIXED_METABOLITES_LIST,
  available_cols = met_cols_all,
  remove_prefix = "^met_"
)

simplify_species_name <- function(x) {
  raw <- gsub("^tax_s__", "", x)
  parts <- strsplit(raw, "_")[[1]]
  if (length(parts) >= 2) {
    paste0(substr(parts[1], 1, 1), ". ", parts[2])
  } else {
    raw
  }
}

species_labels <- vapply(species_cols, simplify_species_name, character(1))
met_labels <- gsub("^met_", "", met_cols)

get_control_met_order <- function() {
  df_species <- merged_data[merged_data$group == "control", species_cols, drop = FALSE]
  df_met <- merged_data[merged_data$group == "control", met_cols, drop = FALSE]

  colnames(df_species) <- species_labels
  colnames(df_met) <- met_labels

  all_mat <- cbind(df_species, df_met)
  rc <- rcorr(as.matrix(all_mat), type = "spearman")

  n_sp <- ncol(df_species)
  n_met <- ncol(df_met)
  r_cross <- rc$r[1:n_sp, (n_sp + 1):(n_sp + n_met), drop = FALSE]

  if (ncol(r_cross) <= 1) {
    return(colnames(r_cross))
  }

  colnames(r_cross)[hclust(dist(t(r_cross)))$order]
}

ctrl_met_order <- get_control_met_order()

message("Matched species columns: ", paste(species_cols, collapse = "; "))
message("Matched metabolite columns: ", paste(met_cols, collapse = "; "))

make_cross_corr_heatmap <- function(
  group_name,
  title_label,
  col_title_color,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  column_order = NULL
) {
  df_species <- merged_data[merged_data$group == group_name, species_cols, drop = FALSE]
  df_met <- merged_data[merged_data$group == group_name, met_cols, drop = FALSE]

  colnames(df_species) <- species_labels
  colnames(df_met) <- met_labels

  all_mat <- cbind(df_species, df_met)
  rc <- rcorr(as.matrix(all_mat), type = "spearman")

  n_sp <- ncol(df_species)
  n_met <- ncol(df_met)

  r_cross <- rc$r[1:n_sp, (n_sp + 1):(n_sp + n_met), drop = FALSE]
  p_cross <- rc$P[1:n_sp, (n_sp + 1):(n_sp + n_met), drop = FALSE]

  sig_mat <- matrix("", nrow = nrow(p_cross), ncol = ncol(p_cross))
  sig_mat[p_cross < 0.001] <- "***"
  sig_mat[p_cross >= 0.001 & p_cross < 0.01] <- "**"
  sig_mat[p_cross >= 0.01 & p_cross < 0.05] <- "*"

  write.csv(
    cbind(Species = rownames(r_cross), as.data.frame(r_cross, check.names = FALSE)),
    file.path(output_dir, paste0("cross_corr_r_", group_name, ".csv")),
    row.names = FALSE
  )
  write.csv(
    cbind(Species = rownames(p_cross), as.data.frame(p_cross, check.names = FALSE)),
    file.path(output_dir, paste0("cross_corr_p_", group_name, ".csv")),
    row.names = FALSE
  )

  col_fun <- colorRamp2(
    c(-1, -0.5, 0, 0.5, 1),
    c("#4575B4", "#91BFDB", "white", "#FC8D59", "#D73027")
  )

  Heatmap(
    r_cross,
    name = "Spearman r",
    col = col_fun,
    cell_fun = function(j, i, x, y, width, height, fill) {
      grid.text(sig_mat[i, j], x, y, gp = gpar(fontsize = 8, col = "black"))
    },
    row_names_gp = gpar(fontsize = 9, fontface = "italic"),
    column_names_gp = gpar(fontsize = 8),
    column_names_rot = 45,
    cluster_rows = cluster_rows,
    cluster_columns = cluster_columns,
    column_order = column_order,
    show_row_dend = FALSE,
    show_column_dend = FALSE,
    column_title = title_label,
    column_title_gp = gpar(fontsize = 12, fontface = "bold", col = col_title_color),
    width = unit(10, "cm"),
    height = unit(6.5, "cm"),
    show_heatmap_legend = (group_name == "control")
  )
}

ht_ctrl <- make_cross_corr_heatmap("control", "Ctrl", "#2E86AB")
ht_well <- make_cross_corr_heatmap(
  "CRC_well_diff",
  "CRC-Well",
  "#F18F01",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  column_order = ctrl_met_order
)
ht_poor <- make_cross_corr_heatmap(
  "CRC_poor_diff",
  "CRC-Poor",
  "#D7263D",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  column_order = ctrl_met_order
)

lgd_sig <- Legend(
  labels = c("p <= 0.05", "p <= 0.01", "p <= 0.001"),
  title = "Significance:",
  direction = "horizontal",
  nrow = 1,
  grid_height = unit(6, "mm"),
  grid_width = unit(6, "mm"),
  labels_gp = gpar(fontsize = 10),
  title_gp = gpar(fontsize = 11, fontface = "bold"),
  type = "points",
  pch = NA,
  graphics = list(
    function(x, y, w, h) grid.text("*", x, y, gp = gpar(fontsize = 12)),
    function(x, y, w, h) grid.text("**", x, y, gp = gpar(fontsize = 12)),
    function(x, y, w, h) grid.text("***", x, y, gp = gpar(fontsize = 12))
  )
)

png(
  file.path(output_dir, "heatmap_species_metabolite_3groups.png"),
  width = 30,
  height = 8,
  units = "in",
  res = 300
)
draw(
  ht_ctrl + ht_well + ht_poor,
  column_title = "Species-Metabolite Correlation Heatmaps across CRC Subtypes",
  column_title_gp = gpar(fontsize = 15, fontface = "bold"),
  heatmap_legend_side = "right",
  annotation_legend_side = "bottom",
  annotation_legend_list = list(lgd_sig)
)
dev.off()

pdf(
  file.path(output_dir, "heatmap_species_metabolite_3groups.pdf"),
  width = 30,
  height = 8
)
draw(
  ht_ctrl + ht_well + ht_poor,
  column_title = "Species-Metabolite Correlation Heatmaps across CRC Subtypes",
  column_title_gp = gpar(fontsize = 15, fontface = "bold"),
  heatmap_legend_side = "right",
  annotation_legend_side = "bottom",
  annotation_legend_list = list(lgd_sig)
)
dev.off()

message("Done: ", output_dir)