#!/usr/bin/env Rscript

#' Read Maaslin2 output, filter for diff_stage, and draw a lollipop chart of coefficients.
library(tidyverse)
library(RColorBrewer)

theme_microbiome <- function(){ 
  theme_bw(base_size = 14) + 
  theme( 
    panel.grid = element_blank(), 
    panel.border = element_blank(), 
    
    axis.text = element_text(color="black"), 
    axis.title = element_text(size=14), 
    
    legend.background = element_blank(), 
    legend.key = element_blank(), 
    
    strip.background = element_blank(), 
    strip.text = element_text(size=13, face="bold") 
  ) 
} 

group_colors <- c("control" = "#2E86AB", "CRC_well_diff" = "#F18F01", 
                  "CRC_poor_diff" = "#D7263D")

input_path <- "results/maaslin2/significant_results.tsv"
target_metadata <- "diff_stage"
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

maaslin_data <- readr::read_tsv(input_path)
available_features <- unique(maaslin_data$feature)
target_features <- match_species_to_cols(target_features, available_features)

output_path <- "results/R_plots/lollipop_plot/diff_lollipop_plot.png"

# helpers ------------------------------------------------------------------
star_from_qval <- function(qval) {
  case_when(
    qval <= 1e-3 ~ "***",
    qval <= 1e-2 ~ "**",
    qval <= 5e-2 ~ "*",
    TRUE ~ ""
  )
}

readr::read_tsv(input_path) %>%
  filter(metadata == target_metadata, feature %in% target_features) %>%
  mutate(
    feature = factor(feature, levels = target_features),
    feature = fct_relabel(feature, ~gsub("^tax_s__", "", .x)),
    qval_star = star_from_qval(qval)
  ) -> plot_data

if (nrow(plot_data) == 0) {
  stop("No features found for the specified metadata and feature list.")
}

plot_data <- plot_data %>%
  mutate(feature = fct_reorder(feature, coef),
         coef_label = sprintf("%+.2f", coef))

max_coef <- max(abs(plot_data$coef))
label_offset <- max_coef * 0.06
plot_data <- plot_data %>%
  mutate(
    direction = ifelse(coef >= 0, 1, -1),
    star_x = coef + direction * label_offset
  )

plot <- ggplot(plot_data, aes(x = coef, y = feature)) +
  geom_vline(xintercept = 0, linetype = "solid", color = "black", size = 0.5) +
  geom_segment(aes(x = 0, xend = coef, yend = feature, color = coef > 0), size = 1.2, alpha = 0.7) +
  geom_point(aes(color = coef > 0, size = -log10(qval)), alpha = 0.9) +
  geom_text(aes(x = ifelse(coef < 0, 0.02, -0.02), label = feature),
            hjust = ifelse(plot_data$coef < 0, 0, 1), size = 3, color = "black", fontface = "italic") +
  geom_text(data = filter(plot_data, qval_star != ""),
            aes(x = star_x, label = qval_star),
            color = "black", size = 5, vjust = 0.7) +
  scale_color_manual(values = c("TRUE" = "#D55E00", "FALSE" = "#0072B2"), 
                     labels = c("TRUE" = "Enriched", "FALSE" = "Depleted"),
                     name = "Association") +
  scale_size_continuous(range = c(3, 7), name = "-log10(q-value)") +
  labs(
    title = "Differential Taxa by Disease Stage",
    subtitle = "Maaslin2 Linear Model Coefficients",
    x = "Coefficient (Effect Size)",
    y = NULL,
    caption = "* q ≤ 0.05, ** q ≤ 0.01, *** q ≤ 0.001"
  ) +
  theme_microbiome() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.x = element_line(color = "black", linewidth = 0.5),
    axis.text.y = element_blank(),
    axis.text.x = element_text(color = "black"),
    legend.position = "right",
    legend.box = "vertical",
    legend.spacing.y = unit(0.2, "cm"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, color = "gray30"),
    plot.caption = element_text(hjust = 1, size = 8, face = "italic")
  )

plot_height <- max(5, nrow(plot_data) * 0.5)

ggsave(output_path, plot, width = 10, height = plot_height, dpi = 300)

message("Saved lollipop plot to ", output_path)
