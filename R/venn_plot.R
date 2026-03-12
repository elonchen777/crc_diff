# ============================================================
# Publication-style Venn plot for CRC differentiation subtypes
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggVennDiagram)
  library(gridExtra)
  library(grid)
})

lda_cutoff <- 2.5
padj_cutoff <- 0.05
max_rows_per_column <- 12
well_color <- "#A23B72"
poor_color <- "#F18F01"
shared_color <- "#C9653A"

well_diff_file <- "results/R_plots/lefse_analysis/LEFSE_control_vs_CRC_well_diff_results.csv"
poor_diff_file <- "results/R_plots/lefse_analysis/LEFSE_control_vs_CRC_poor_diff_results.csv"

if (!file.exists(well_diff_file) || !file.exists(poor_diff_file)) {
  stop("Input files not found. Please ensure LEfSe analysis has been run.")
}

df_well <- read.csv(well_diff_file, stringsAsFactors = FALSE)
df_poor <- read.csv(poor_diff_file, stringsAsFactors = FALSE)

extract_species <- function(taxa_str) {
  if (is.na(taxa_str) || taxa_str == "") {
    return(NA_character_)
  }

  parts <- strsplit(taxa_str, "\\|")[[1]]
  species <- trimws(parts[length(parts)])
  species <- sub("^s__", "", species)

  if (species == "" || species == "NA") {
    return(NA_character_)
  }

  species
}

get_sig_species <- function(df, target_group, lda_threshold, padj_threshold) {
  required_cols <- c("Taxa", "Group", "LDA", "P.adj")
  missing_cols <- setdiff(required_cols, names(df))

  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  selected <- df[
    !is.na(df$P.adj) &
      !is.na(df$LDA) &
      df$P.adj < padj_threshold &
      df$LDA > lda_threshold &
      df$Group == target_group,
    , drop = FALSE
  ]

  species <- unique(vapply(selected$Taxa, extract_species, character(1)))
  sort(stats::na.omit(species))
}

format_species <- function(species) {
  gsub("_", " ", species)
}

split_into_columns <- function(species, max_rows = 10) {
  if (length(species) == 0) {
    return(list("None"))
  }

  split(species, ceiling(seq_along(species) / max_rows))
}

make_metric_box <- function(title, value, fill, subtitle) {
  grobTree(
    roundrectGrob(
      x = 0.5,
      y = 0.5,
      width = 0.98,
      height = 0.94,
      r = unit(0.08, "snpc"),
      gp = gpar(fill = fill, col = NA)
    ),
    textGrob(
      label = title,
      x = 0.08,
      y = 0.78,
      just = "left",
      gp = gpar(col = "white", fontsize = 11, fontface = "bold", fontfamily = "serif")
    ),
    textGrob(
      label = value,
      x = 0.08,
      y = 0.44,
      just = "left",
      gp = gpar(col = "white", fontsize = 25, fontface = "bold", fontfamily = "serif")
    ),
    textGrob(
      label = paste(strwrap(subtitle, width = 34), collapse = "\n"),
      x = 0.08,
      y = 0.15,
      just = c("left", "bottom"),
      gp = gpar(col = "white", fontsize = 8.8, lineheight = 1.1, fontfamily = "sans")
    )
  )
}

make_note_box <- function(lines) {
  grobTree(
    roundrectGrob(
      x = 0.5,
      y = 0.5,
      width = 0.98,
      height = 0.94,
      r = unit(0.06, "snpc"),
      gp = gpar(fill = "#F7F3EC", col = "#D9CFBE", lwd = 1)
    ),
    textGrob(
      label = "Figure notes",
      x = 0.08,
      y = 0.82,
      just = "left",
      gp = gpar(col = "#3A342F", fontsize = 11, fontface = "bold", fontfamily = "serif")
    ),
    textGrob(
      label = paste(lines, collapse = "\n"),
      x = 0.08,
      y = 0.63,
      just = c("left", "top"),
      gp = gpar(col = "#4F4943", fontsize = 8.8, lineheight = 1.25, fontfamily = "sans")
    )
  )
}

make_species_block <- function(species, title, accent, description, max_rows = 10) {
  display_species <- format_species(sort(species))
  columns <- split_into_columns(display_species, max_rows = max_rows)

  text_columns <- lapply(columns, function(entries) {
    textGrob(
      paste(entries, collapse = "\n"),
      x = 0,
      y = 1,
      just = c("left", "top"),
      gp = gpar(
        col = "#2F2A24",
        fontsize = 7.2,
        fontface = if (length(species) == 0) "plain" else "italic",
        lineheight = 1.08,
        fontfamily = "serif"
      )
    )
  })

  body_table <- arrangeGrob(grobs = text_columns, ncol = length(text_columns))

  header <- grobTree(
    roundrectGrob(
      x = 0.5,
      y = 0.5,
      width = 0.98,
      height = 0.94,
      r = unit(0.06, "snpc"),
      gp = gpar(fill = accent, col = NA)
    ),
    textGrob(
      sprintf("%s (n = %d)", title, length(species)),
      x = 0.05,
      y = 0.68,
      just = "left",
      gp = gpar(col = "white", fontsize = 10.5, fontface = "bold", fontfamily = "serif")
    ),
    textGrob(
      paste(strwrap(description, width = 38), collapse = "\n"),
      x = 0.05,
      y = 0.2,
      just = c("left", "bottom"),
      gp = gpar(col = "white", fontsize = 7.4, lineheight = 1.05, fontfamily = "sans")
    )
  )

  body <- grobTree(
    roundrectGrob(
      x = 0.5,
      y = 0.5,
      width = 0.98,
      height = 0.96,
      r = unit(0.05, "snpc"),
      gp = gpar(fill = "#FBF9F5", col = "#DDD5C9", lwd = 1)
    ),
    editGrob(body_table, vp = viewport(x = 0.5, y = 0.5, width = 0.9, height = 0.86))
  )

  arrangeGrob(header, body, ncol = 1, heights = c(0.6, 1.7))
}

well_set <- get_sig_species(df_well, "CRC_well_diff", lda_cutoff, padj_cutoff)
poor_set <- get_sig_species(df_poor, "CRC_poor_diff", lda_cutoff, padj_cutoff)

if (length(well_set) == 0 && length(poor_set) == 0) {
  stop("No species passed the selected thresholds for either CRC subtype.")
}

common_set <- intersect(well_set, poor_set)
well_only <- setdiff(well_set, poor_set)
poor_only <- setdiff(poor_set, well_set)

venn_labels <- c("CRC well-diff", "CRC poor-diff")
venn_colors <- c(well_color, poor_color)

venn_list <- list(
  "CRC well-diff" = well_set,
  "CRC poor-diff" = poor_set
)

p_venn <- ggVennDiagram(
  venn_list,
  label = "count",
  label_alpha = 0,
  label_geom = "label",
  label_color = "#2E2924",
  label_size = 5.2,
  category.names = venn_labels,
  set_color = venn_colors,
  set_size = 5.3,
  edge_size = 1.5
) +
  scale_fill_gradient(low = "#FFFDF9", high = "#D8C6A4") +
  theme_void() +
  theme(
    legend.position = "none",
    plot.title = element_text(
      hjust = 0.5,
      size = 18,
      face = "bold",
      family = "serif",
      color = "#2F2A24",
      margin = margin(b = 6)
    ),
    plot.subtitle = element_text(
      hjust = 0.5,
      size = 10.5,
      family = "sans",
      color = "#58514A",
      lineheight = 1.15,
      margin = margin(b = 10)
    ),
    plot.caption = element_text(
      hjust = 0.5,
      size = 8.6,
      family = "sans",
      color = "#6B645D",
      margin = margin(t = 8)
    ),
    plot.margin = margin(12, 12, 8, 12),
    plot.background = element_rect(fill = "#FCFBF8", color = NA),
    panel.background = element_rect(fill = "#FCFBF8", color = NA)
  ) +
  labs(
    title = "Shared and Subtype-Specific Differential Species",
    subtitle = sprintf(
      "CRC-enriched species identified by LEfSe under LDA > %.1f and adjusted P < %.2f",
      lda_cutoff,
      padj_cutoff
    ),
    caption = sprintf(
      "CRC well-diff: n = %d | CRC poor-diff: n = %d | Shared species: n = %d",
      length(well_set),
      length(poor_set),
      length(common_set)
    )
  )

summary_panel <- arrangeGrob(
  make_species_block(
    species = well_only,
    title = "CRC well-diff specific",
    accent = well_color,
    description = "Subtype-specific taxa listed with smaller typography for compact manuscript display.",
    max_rows = max_rows_per_column
  ),
  make_species_block(
    species = common_set,
    title = "Intersection",
    accent = shared_color,
    description = "Shared taxa retained in both CRC well-diff and CRC poor-diff.",
    max_rows = max_rows_per_column
  ),
  make_species_block(
    species = poor_only,
    title = "CRC poor-diff specific",
    accent = poor_color,
    description = "Subtype-specific taxa listed with smaller typography for compact manuscript display.",
    max_rows = max_rows_per_column
  ),
  make_note_box(c(
    "Input: LEfSe differential species tables versus control.",
    sprintf("Selection rule: LDA > %.1f and adjusted P < %.2f.", lda_cutoff, padj_cutoff),
    "Species names are shown in italic below for direct manuscript interpretation."
  )),
  ncol = 1,
  heights = c(0.72, 0.72, 0.72, 1.45, 1.45, 1.45, 1.02)
)

top_panel <- arrangeGrob(p_venn, summary_panel, ncol = 2, widths = c(1.25, 0.9))
final_plot <- top_panel

output_dir <- "results/R_plots/venn_analysis"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

output_file_base <- file.path(output_dir, "venn_publication_style_lda2.5")

write.csv(
  data.frame(Species = format_species(sort(well_only))),
  file.path(output_dir, "CRC_well_diff_specific_species.csv"),
  row.names = FALSE
)
write.csv(
  data.frame(Species = format_species(sort(common_set))),
  file.path(output_dir, "CRC_shared_species.csv"),
  row.names = FALSE
)
write.csv(
  data.frame(Species = format_species(sort(poor_only))),
  file.path(output_dir, "CRC_poor_diff_specific_species.csv"),
  row.names = FALSE
)

max_columns_needed <- max(
  length(split_into_columns(well_only, max_rows_per_column)),
  length(split_into_columns(common_set, max_rows_per_column)),
  length(split_into_columns(poor_only, max_rows_per_column))
)
figure_height <- 12.2 + max(0, max_columns_needed - 1) * 0.5

png(
  paste0(output_file_base, ".png"),
  width = 13,
  height = figure_height,
  units = "in",
  res = 400,
  bg = "white"
)
grid.newpage()
grid.draw(final_plot)
dev.off()

pdf(
  paste0(output_file_base, ".pdf"),
  width = 13,
  height = figure_height,
  bg = "white"
)
grid.newpage()
grid.draw(final_plot)
dev.off()

message("Success! Plots saved to: ", output_dir)
message("Thresholds used: LDA > ", lda_cutoff, "; adjusted P < ", padj_cutoff)
message(
  "Counts - Well differentiated: ", length(well_set),
  ", Poorly differentiated: ", length(poor_set),
  ", Shared: ", length(common_set)
)