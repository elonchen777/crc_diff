# ============================================================
# Three-group species-KO Spearman correlation heatmaps
# KO columns are ordered by KEGG KOEntry order.
# ============================================================

library(ComplexHeatmap)
library(circlize)
library(data.table)
library(Hmisc)
library(grid)
library(KEGGREST)

merged_data <- fread(
  "dataset/merged_dataset_processed.csv",
  stringsAsFactors = FALSE,
  data.table = FALSE
)

output_dir <- "results/R_plots/heatmap_species_ko_correlation_plot"
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

KO_SIGNIFICANT_PATH <- "results/maaslin2_ko/significant_results.tsv"
KEGG_CACHE_DIR <- "dataset/KEGG/cache"
N_KO <- 50

normalize_text <- function(x) {
  tolower(gsub("[^a-z0-9]", "", x))
}

match_targets_to_cols <- function(target_list, available_cols, remove_prefix = "") {
  display_names <- gsub(remove_prefix, "", available_cols)
  norm_display <- normalize_text(display_names)

  matched_cols <- character(0)
  missing_targets <- character(0)

  for (target in target_list) {
    target_core <- gsub(remove_prefix, "", target)
    norm_target <- normalize_text(target_core)

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

select_ko_targets <- function(significant_path, n_ko = 40) {
  if (!file.exists(significant_path)) {
    stop("KO significant file not found: ", significant_path)
  }

  sig <- fread(significant_path, sep = "\t", data.table = FALSE)
  feature_col <- if ("feature" %in% colnames(sig)) "feature" else colnames(sig)[1]

  ko_features <- unique(sig[[feature_col]])
  ko_features <- ko_features[grepl("^kegg_K[0-9]{5}$", ko_features)]

  if (length(ko_features) == 0) {
    stop("No KO features like 'kegg_Kxxxxx' found in: ", significant_path)
  }

  gsub("^kegg_", "", ko_features[seq_len(min(n_ko, length(ko_features)))])
}

fetch_kegg_pathway_classification <- function(ko_ids, cache_file, cache_dir = KEGG_CACHE_DIR) {
  if (!dir.exists(cache_dir)) dir.create(cache_dir, recursive = TRUE)

  link_cache_file <- file.path(cache_dir, "kegg_ko_pathway_links.csv")
  pathway_cache_file <- file.path(cache_dir, "kegg_pathway_metadata.csv")

  if (file.exists(cache_file)) {
    ann_local <- fread(cache_file, data.table = FALSE)
    if (all(c("ko_id", "pathway_id", "pathway_name", "level1", "level2") %in% colnames(ann_local))) {
      miss_ko <- setdiff(ko_ids, ann_local$ko_id)
      if (length(miss_ko) == 0) {
        message("Using per-run KEGG annotation cache: ", cache_file)
        return(ann_local)
      }
    }
  }

  if (file.exists(link_cache_file)) {
    link_df_all <- fread(link_cache_file, data.table = FALSE)
    if (!all(c("ko_id", "pathway_id") %in% colnames(link_df_all))) {
      link_df_all <- NULL
    }
  } else {
    link_df_all <- NULL
  }

  if (is.null(link_df_all)) {
    message("Building KEGG KO-pathway cache from KEGGREST (first run may take longer)...")
    link_vec <- tryCatch(
      keggLink("pathway", "ko"),
      error = function(e) {
        stop(
            "Details: ", conditionMessage(e),
          call. = FALSE
        )
      }
    )
    link_df_all <- unique(data.frame(
      ko_id = sub("^ko:", "", names(link_vec)),
      pathway_id = sub("^path:", "", unname(link_vec)),
      stringsAsFactors = FALSE
    ))
    write.csv(link_df_all, link_cache_file, row.names = FALSE)
  }

  link_df <- unique(link_df_all[link_df_all$ko_id %in% ko_ids, , drop = FALSE])
  if (nrow(link_df) == 0) {
    stop("No KO-pathway links found for selected KOs (from cache/KEGGREST)")
  }

  if (file.exists(pathway_cache_file)) {
    pathway_meta <- fread(pathway_cache_file, data.table = FALSE)
    if (!all(c("pathway_id", "pathway_name", "level1", "level2") %in% colnames(pathway_meta))) {
      pathway_meta <- NULL
    }
  } else {
    pathway_meta <- NULL
  }

  if (is.null(pathway_meta)) {
    message("Building KEGG pathway metadata cache from KEGGREST (first run may take longer)...")
    pathway_name_vec <- tryCatch(
      keggList("pathway"),
      error = function(e) {
        stop(
          "Unable to fetch pathway metadata list from KEGGREST (rest.kegg.jp). ",
          "If running offline, prepare cache file first: ", pathway_cache_file,
          ". Details: ", conditionMessage(e),
          call. = FALSE
        )
      }
    )
    pathway_meta <- data.frame(
      pathway_id = sub("^path:", "", names(pathway_name_vec)),
      pathway_name = unname(pathway_name_vec),
      level1 = NA_character_,
      level2 = NA_character_,
      stringsAsFactors = FALSE
    )

    for (i in seq_len(nrow(pathway_meta))) {
      pid <- pathway_meta$pathway_id[i]
      entry <- tryCatch(
        keggGet(paste0("path:", pid))[[1]],
        error = function(e) {
          stop(
            "Unable to fetch pathway class details from KEGGREST for ", pid,
            ". If running offline, prepare cache file first: ", pathway_cache_file,
            ". Details: ", conditionMessage(e),
            call. = FALSE
          )
        }
      )
      cls <- if (!is.null(entry) && !is.null(entry$CLASS)) entry$CLASS[1] else NA_character_
      cls_parts <- trimws(strsplit(cls, ";", fixed = TRUE)[[1]])
      pathway_meta$level1[i] <- if (length(cls_parts) >= 1) cls_parts[1] else NA_character_
      pathway_meta$level2[i] <- if (length(cls_parts) >= 2) cls_parts[2] else NA_character_
      if (i %% 25 == 0) Sys.sleep(0.2)
    }

    write.csv(pathway_meta, pathway_cache_file, row.names = FALSE)
  }

  ann <- merge(link_df, pathway_meta, by = "pathway_id", all.x = TRUE)
  ann$level1[is.na(ann$level1) | ann$level1 == ""] <- "Unclassified"
  ann$level2[is.na(ann$level2) | ann$level2 == ""] <- "Unclassified"
  ann$pathway_name[is.na(ann$pathway_name) | ann$pathway_name == ""] <- ann$pathway_id

  ann <- ann[order(ann$ko_id, ann$level1, ann$level2, ann$pathway_name, ann$pathway_id), ]
  ann_primary <- ann[!duplicated(ann$ko_id), c("ko_id", "pathway_id", "pathway_name", "level1", "level2"), drop = FALSE]

  write.csv(ann_primary, cache_file, row.names = FALSE)
  ann_primary
}

order_kos_by_keggrest <- function(ko_cols, output_dir) {
  ko_ids <- gsub("^kegg_", "", ko_cols)
  cache_file <- file.path(output_dir, "ko_keggrest_pathway_annotation.csv")

  ann <- fetch_kegg_pathway_classification(ko_ids, cache_file, cache_dir = KEGG_CACHE_DIR)

  ann <- merge(
    data.frame(ko_id = ko_ids, ko_col = ko_cols, stringsAsFactors = FALSE),
    ann,
    by = "ko_id",
    all.x = TRUE
  )

  unmapped_ko <- ann$ko_id[is.na(ann$pathway_id) | ann$pathway_id == ""]
  if (length(unmapped_ko) > 0) {
    unmapped_ko <- unique(unmapped_ko)
    message("KOs without KEGG pathway mapping: ", paste(unmapped_ko, collapse = ", "))
    write.csv(
      data.frame(ko_id = unmapped_ko, stringsAsFactors = FALSE),
      file.path(output_dir, "ko_without_pathway_mapping.csv"),
      row.names = FALSE
    )
  }

  ann$level1[is.na(ann$level1)] <- "Unclassified"
  ann$level2[is.na(ann$level2)] <- "Unclassified"
  ann$pathway_name[is.na(ann$pathway_name)] <- "Unclassified"

  ann <- ann[order(ann$level1, ann$level2, ann$pathway_name, ann$ko_id), ]
  ordered_cols <- ann$ko_col

  ann_out <- ann[, c("ko_col", "ko_id", "level1", "level2", "pathway_name", "pathway_id"), drop = FALSE]
  write.csv(ann_out, file.path(output_dir, "ko_column_order_with_kegg_class.csv"), row.names = FALSE)

  list(ordered_cols = ordered_cols, annotation = ann_out)
}

species_cols_all <- colnames(merged_data)[grepl("^tax_s__", colnames(merged_data))]
ko_cols_all <- colnames(merged_data)[grepl("^kegg_", colnames(merged_data))]

species_cols <- match_targets_to_cols(
  target_list = FIXED_SPECIES_LIST,
  available_cols = species_cols_all,
  remove_prefix = "^tax_s__"
)

ko_targets <- select_ko_targets(KO_SIGNIFICANT_PATH, n_ko = N_KO)
ko_cols <- match_targets_to_cols(
  target_list = ko_targets,
  available_cols = ko_cols_all,
  remove_prefix = "^kegg_"
)
ko_order_info <- order_kos_by_keggrest(ko_cols, output_dir)
ko_cols <- ko_order_info$ordered_cols

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
ko_labels <- gsub("^kegg_", "", ko_cols)
ko_level2_labels <- setNames(
  ko_order_info$annotation$level2,
  ko_order_info$annotation$ko_col
)[ko_cols]

ko_row_split <- factor(ko_level2_labels, levels = unique(ko_level2_labels))
ko_row_order <- seq_along(ko_cols)

ko_relative_mat <- as.matrix(merged_data[, ko_cols, drop = FALSE])
mode(ko_relative_mat) <- "numeric"
# KO columns in merged_dataset_processed.csv are log1p-transformed; convert back to relative abundance.
ko_relative_mat <- pmax(expm1(ko_relative_mat), 0)

ko_abundance_mat <- sapply(c("control", "CRC_well_diff", "CRC_poor_diff"), function(g) {
  idx <- merged_data$group == g
  colMeans(ko_relative_mat[idx, , drop = FALSE], na.rm = TRUE)
})
colnames(ko_abundance_mat) <- c("Ctrl", "Well", "Poor")
rownames(ko_abundance_mat) <- ko_labels

# Row-wise normalization: each KO sums to 1 across three groups.
ko_abundance_prop <- ko_abundance_mat
ko_row_sum <- rowSums(ko_abundance_mat, na.rm = TRUE)
nonzero_idx <- !is.na(ko_row_sum) & ko_row_sum > 0
ko_abundance_prop[nonzero_idx, ] <- ko_abundance_mat[nonzero_idx, , drop = FALSE] / ko_row_sum[nonzero_idx]
ko_abundance_prop[!nonzero_idx, ] <- 0

level2_levels <- unique(ko_level2_labels)
level2_colors <- structure(
  grDevices::hcl(
    h = seq(15, 375, length.out = length(level2_levels) + 1)[-1],
    c = 85,
    l = 65
  ),
  names = level2_levels
)

ko_level2_factor <- factor(ko_level2_labels, levels = level2_levels)

ko_left_annotation <- rowAnnotation(
  Level2 = ko_level2_factor,
  Abundance = anno_barplot(
    ko_abundance_prop,
    gp = gpar(fill = c("#2E86AB", "#F18F01", "#D7263D"), col = NA),
    bar_width = 0.85,
    ylim = c(0, 1),
    axis_param = list(side = "top", at = c(0, 0.5, 1), labels = c("0", "0.5", "1"), gp = gpar(fontsize = 7)),
    border = FALSE,
    width = unit(24, "mm")
  ),
  col = list(Level2 = level2_colors),
  annotation_name_gp = gpar(fontsize = 8, fontface = "bold"),
  annotation_name_rot = 0,
  annotation_name_side = "top",
  gap = unit(2, "mm"),
  show_annotation_name = c(FALSE, TRUE),
  show_legend = FALSE
)

message("Matched species columns: ", paste(species_cols, collapse = "; "))
message("Matched KO columns: ", paste(ko_cols, collapse = "; "))

make_cross_corr_heatmap <- function(group_name, title_label, col_title_color, show_left_anno = FALSE, show_row_names = TRUE) {
  df_species <- merged_data[merged_data$group == group_name, species_cols, drop = FALSE]
  df_ko <- merged_data[merged_data$group == group_name, ko_cols, drop = FALSE]

  colnames(df_species) <- species_labels
  colnames(df_ko) <- ko_labels

  all_mat <- cbind(df_species, df_ko)
  rc <- rcorr(as.matrix(all_mat), type = "spearman")

  n_sp <- ncol(df_species)
  n_ko <- ncol(df_ko)

  r_cross <- rc$r[1:n_sp, (n_sp + 1):(n_sp + n_ko), drop = FALSE]
  p_cross <- rc$P[1:n_sp, (n_sp + 1):(n_sp + n_ko), drop = FALSE]

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

  r_plot <- t(r_cross)
  p_plot <- t(p_cross)
  cell_size_mm <- 5

  sig_mat <- matrix("", nrow = nrow(p_plot), ncol = ncol(p_plot))
  sig_mat[p_plot < 0.001] <- "***"
  sig_mat[p_plot >= 0.001 & p_plot < 0.01] <- "**"
  sig_mat[p_plot >= 0.01 & p_plot < 0.05] <- "*"

  Heatmap(
    r_plot,
    name = "Spearman r",
    col = col_fun,
    cell_fun = function(j, i, x, y, width, height, fill) {
      grid.text(sig_mat[i, j], x, y, gp = gpar(fontsize = 5, col = "black"))
    },
    row_names_gp = gpar(fontsize = 6),
    show_row_names = show_row_names,
    row_names_side = "left",
    column_names_gp = gpar(fontsize = 9, fontface = "italic"),
    column_names_rot = 45,
    show_row_dend = FALSE,
    show_column_dend = FALSE,
    row_split = ko_row_split,
    row_order = ko_row_order,
    row_gap = unit(2.5, "mm"),
    row_title_rot = 0,
    row_title_gp = gpar(fontsize = 8, fontface = "bold"),
    border = FALSE,
    rect_gp = gpar(col = NA),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    column_order = seq_along(species_labels),
    left_annotation = if (show_left_anno) ko_left_annotation else NULL,
    column_title = title_label,
    column_title_gp = gpar(fontsize = 12, fontface = "bold", col = col_title_color),
    width = unit(ncol(r_plot) * cell_size_mm, "mm"),
    height = unit(nrow(r_plot) * cell_size_mm, "mm"),
    show_heatmap_legend = (group_name == "control")
  )
}

ht_ctrl <- make_cross_corr_heatmap("control", "Ctrl", "#2E86AB", show_left_anno = TRUE, show_row_names = TRUE)
ht_well <- make_cross_corr_heatmap("CRC_well_diff", "CRC-Well", "#F18F01", show_left_anno = FALSE, show_row_names = FALSE)
ht_poor <- make_cross_corr_heatmap("CRC_poor_diff", "CRC-Poor", "#D7263D", show_left_anno = FALSE, show_row_names = FALSE)

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
    function(x, y, w, h) grid.text("*", x, y, gp = gpar(fontsize = 9)),
    function(x, y, w, h) grid.text("**", x, y, gp = gpar(fontsize = 9)),
    function(x, y, w, h) grid.text("***", x, y, gp = gpar(fontsize = 9))
  )
)

png(
  file.path(output_dir, "heatmap_species_ko_3groups.png"),
  width = 26,
  height = 12,
  units = "in",
  res = 300
)
draw(
  ht_ctrl + ht_well + ht_poor,
  column_title = "Species-KO (KEGG-ordered) Correlation Heatmaps across CRC Subtypes",
  column_title_gp = gpar(fontsize = 15, fontface = "bold"),
  heatmap_legend_side = "right",
  annotation_legend_side = "bottom",
  annotation_legend_list = list(lgd_sig),
  merge_legend = TRUE
)
dev.off()

pdf(
  file.path(output_dir, "heatmap_species_ko_3groups.pdf"),
  width = 26,
  height = 12
)
draw(
  ht_ctrl + ht_well + ht_poor,
  column_title = "Species-KO (KEGG-ordered) Correlation Heatmaps across CRC Subtypes",
  column_title_gp = gpar(fontsize = 15, fontface = "bold"),
  heatmap_legend_side = "right",
  annotation_legend_side = "bottom",
  annotation_legend_list = list(lgd_sig),
  merge_legend = TRUE
)
dev.off()

message("Done: ", output_dir)

