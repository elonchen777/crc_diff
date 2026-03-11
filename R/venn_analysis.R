well_diff_file <- "d:/cyl/biotech/crc_diff/results/R_plots/lefse_analysis/LEFSE_control_vs_CRC_well_diff_results.csv"
poor_diff_file <- "d:/cyl/biotech/crc_diff/results/R_plots/lefse_analysis/LEFSE_control_vs_CRC_poor_diff_results.csv"

df_well <- read.csv(well_diff_file, stringsAsFactors = FALSE)
df_poor <- read.csv(poor_diff_file, stringsAsFactors = FALSE)

df_well_filtered <- df_well[df_well$P.adj < 0.05 & df_well$LDA > 2, ]
df_poor_filtered <- df_poor[df_poor$P.adj < 0.05 & df_poor$LDA > 2, ]

extract_species <- function(taxa_str) {
  parts <- strsplit(taxa_str, "\\|")[[1]]
  return(parts[length(parts)])
}

well_species <- sapply(df_well_filtered$Taxa, extract_species)
poor_species <- sapply(df_poor_filtered$Taxa, extract_species)

well_set <- unique(well_species)
poor_set <- unique(poor_species)

common <- intersect(well_set, poor_set)
well_only <- setdiff(well_set, poor_set)
poor_only <- setdiff(poor_set, well_set)

cat(strrep("=", 60), "\n")
cat("Venn图分析结果: CRC_well_diff vs CRC_poor_diff\n")
cat("筛选条件: P.adj < 0.05 且 LDA > 2\n")
cat(strrep("=", 60), "\n\n")
cat("CRC_well_diff 特有物种数:", length(well_only), "\n")
cat("CRC_poor_diff 特有物种数:", length(poor_only), "\n")
cat("共有物种数:", length(common), "\n")
cat("CRC_well_diff 总物种数:", length(well_set), "\n")
cat("CRC_poor_diff 总物种数:", length(poor_set), "\n")

cat("\n", strrep("=", 60), "\n")
cat("共有物种 (两个数据集都包含):\n")
cat(strrep("=", 60), "\n")
for (i in seq_along(sort(common))) {
  cat(i, ". ", sort(common)[i], "\n", sep = "")
}

cat("\n", strrep("=", 60), "\n")
cat("CRC_well_diff 特有物种:\n")
cat(strrep("=", 60), "\n")
for (i in seq_along(sort(well_only))) {
  cat(i, ". ", sort(well_only)[i], "\n", sep = "")
}

cat("\n", strrep("=", 60), "\n")
cat("CRC_poor_diff 特有物种:\n")
cat(strrep("=", 60), "\n")
for (i in seq_along(sort(poor_only))) {
  cat(i, ". ", sort(poor_only)[i], "\n", sep = "")
}

library(ggplot2)
library(ggVennDiagram)

venn_data <- list(
  CRC_well_diff = well_set,
  CRC_poor_diff = poor_set
)

p <- ggVennDiagram(venn_data, label_alpha = 0.5, category.names = c("CRC_well_diff", "CRC_poor_diff")) +
  scale_fill_gradient(low = "white", high = "#B3E2CD") +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.margin = margin(1, 1, 1, 1, "cm")
  ) +
  labs(title = "Venn Diagram: CRC_well_diff vs CRC_poor_diff\n(P.adj < 0.05, LDA > 2)")

format_italic_list <- function(species_vec) {
  if (length(species_vec) == 0) return("None")
  species_vec <- gsub("_", " ", species_vec)
  return(paste(species_vec, collapse = "\n"))
}

common_formatted <- format_italic_list(sort(common))
well_only_formatted <- format_italic_list(sort(well_only))
poor_only_formatted <- format_italic_list(sort(poor_only))

p_venn <- p +
  annotate("text", x = 0, y = -0.6, label = common_formatted, size = 2.5, fontface = "italic", hjust = 0.5) +
  annotate("text", x = -0.6, y = -0.6, label = well_only_formatted, size = 2.5, fontface = "italic", hjust = 1) +
  annotate("text", x = 0.6, y = -0.6, label = poor_only_formatted, size = 2.5, fontface = "italic", hjust = 0)

output_file <- "d:/cyl/biotech/crc_diff/results/R_plots/lefse_analysis/Venn_CRC_well_vs_poor_diff.pdf"
ggsave(output_file, p_venn, width = 16, height = 12)
cat("\nVenn图已保存至:", output_file, "\n")

output_png <- "d:/cyl/biotech/crc_diff/results/R_plots/lefse_analysis/Venn_CRC_well_vs_poor_diff.png"
ggsave(output_png, p_venn, width = 16, height = 12, dpi = 300)
cat("Venn图(PNG)已保存至:", output_png, "\n")

library(gridExtra)

make_species_table <- function(species_vec, title) {
  if (length(species_vec) == 0) {
    return(data.frame(Species = "None"))
  }
  species_vec <- gsub("_", " ", sort(species_vec))
  return(data.frame(Species = species_vec))
}

common_df <- make_species_table(common, "Common Species")
well_only_df <- make_species_table(well_only, "CRC_well_diff Only")
poor_only_df <- make_species_table(poor_only, "CRC_poor_diff Only")

p_common <- ggplot(common_df, aes(x = 1, y = seq_along(Species))) +
  geom_text(aes(label = Species), hjust = 0, fontface = "italic", size = 3) +
  labs(title = "Common Species", x = NULL, y = NULL) +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12))

p_well <- ggplot(well_only_df, aes(x = 1, y = seq_along(Species))) +
  geom_text(aes(label = Species), hjust = 0, fontface = "italic", size = 3) +
  labs(title = "CRC_well_diff Only", x = NULL, y = NULL) +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12))

p_poor <- ggplot(poor_only_df, aes(x = 1, y = seq_along(Species))) +
  geom_text(aes(label = Species), hjust = 0, fontface = "italic", size = 3) +
  labs(title = "CRC_poor_diff Only", x = NULL, y = NULL) +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12))

combined <- grid.arrange(p, p_well, p_poor, p_common, ncol = 2, nrow = 2,
                         widths = c(1, 1), heights = c(1, 1))

output_combined <- "d:/cyl/biotech/crc_diff/results/R_plots/lefse_analysis/Venn_CRC_well_vs_poor_diff_combined.pdf"
ggsave(output_combined, combined, width = 16, height = 14)
cat("组合Venn图已保存至:", output_combined, "\n")

output_combined_png <- "d:/cyl/biotech/crc_diff/results/R_plots/lefse_analysis/Venn_CRC_well_vs_poor_diff_combined.png"
ggsave(output_combined_png, combined, width = 16, height = 14, dpi = 300)
cat("组合Venn图(PNG)已保存至:", output_combined_png, "\n")
