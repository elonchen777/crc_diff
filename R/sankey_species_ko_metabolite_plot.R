# ============================================================
# Species -> KO -> Metabolite Sankey plot
# Fixed KO set provided by user
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggalluvial)
  library(ggnewscale)
  library(scales)
})

input_file <- "dataset/merged_dataset_processed.csv"
output_dir <- "results/R_plots/sankey_species_ko_metabolite"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# User-specified KOs and functional labels.
selected_ko_desc <- c(
  K12688 = "Lipid metabolism",
  K02548 = "Transport system",
  K07091 = "Lipid / stress response",
  K10200 = "Oxidative phosphorylation"
#   K10201 = "Oxidative phosphorylation",
#   K02114 = "Electron transport chain",
#   K01512 = "Carbohydrate metabolism",
#   K01647 = "SCFA metabolism",
#   K02804 = "Ribosome",
#   K03785 = "Transcription"
)

FIXED_SPECIES_LIST <- c(
  "Peptostreptococcus_stomatis",
  "Porphyromonas_gingivalis",
  "Prevotella_intermedia",
  "Fusobacterium_periodonticum",
  "Campylobacter_rectus"
#   "Faecalibacterium_prausnitzii",
#   "Roseburia_intestinalis",
#   "Eubacterium_rectale",
#   "Coprococcus_comes",
#   "Ruminococcus_lactaris"
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

N_SPECIES_PER_KO <- 8
N_METABOLITES_PER_KO <- 8
COR_THRESHOLD <- 0.25

if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file)
}

dat <- fread(input_file, data.table = FALSE)

normalize_feature_name <- function(x) {
  x <- trimws(as.character(x))
  x <- tolower(x)
  x <- gsub("^tax_s__", "", x)
  x <- gsub("^met_", "", x)
  x <- gsub("^s__", "", x)
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

match_targets_to_cols <- function(target_list, available_cols, remove_prefix_pattern) {
  display_names <- gsub(remove_prefix_pattern, "", available_cols)
  norm_display <- normalize_feature_name(display_names)

  out <- character(0)
  missing <- character(0)
  for (target in target_list) {
    norm_target <- normalize_feature_name(target)
    idx <- which(norm_display == norm_target)
    if (length(idx) > 0) {
      out <- c(out, available_cols[idx[1]])
    } else {
      missing <- c(missing, target)
    }
  }
  list(cols = unique(out), missing = unique(missing))
}

species_cols_all <- colnames(dat)[grepl("^tax_s__", colnames(dat))]
met_cols_all <- colnames(dat)[grepl("^met_", colnames(dat))]
ko_cols_needed <- paste0("kegg_", names(selected_ko_desc))

species_match <- match_targets_to_cols(FIXED_SPECIES_LIST, species_cols_all, "^tax_s__")
met_match <- match_targets_to_cols(FIXED_METABOLITES_LIST, met_cols_all, "^met_")
species_cols <- species_match$cols
met_cols <- met_match$cols

missing_ko <- setdiff(ko_cols_needed, colnames(dat))
if (length(missing_ko) > 0) {
  stop("These selected KO columns are missing in dataset: ", paste(missing_ko, collapse = ", "))
}

if (length(species_cols) == 0 || length(met_cols) == 0) {
  stop("No matched fixed species or metabolites in dataset.")
}

if (length(species_match$missing) > 0) {
  warning("Fixed species not matched (ignored): ", paste(species_match$missing, collapse = "; "))
}
if (length(met_match$missing) > 0) {
  warning("Fixed metabolites not matched (ignored): ", paste(met_match$missing, collapse = "; "))
}

N_SPECIES_PER_KO <- min(N_SPECIES_PER_KO, length(species_cols))
N_METABOLITES_PER_KO <- min(N_METABOLITES_PER_KO, length(met_cols))

species_mat <- as.matrix(dat[, species_cols, drop = FALSE])
ko_mat <- as.matrix(dat[, ko_cols_needed, drop = FALSE])
met_mat <- as.matrix(dat[, met_cols, drop = FALSE])

mode(species_mat) <- "numeric"
mode(ko_mat) <- "numeric"
mode(met_mat) <- "numeric"

# Correlation matrices (pairwise complete).
cor_sp_ko <- cor(species_mat, ko_mat, use = "pairwise.complete.obs", method = "spearman")
cor_ko_met <- cor(ko_mat, met_mat, use = "pairwise.complete.obs", method = "spearman")

pick_top_edges <- function(cor_vec, candidates, top_n, threshold) {
  keep <- which(!is.na(cor_vec) & abs(cor_vec) >= threshold)
  if (length(keep) == 0) {
    keep <- which(!is.na(cor_vec))
  }
  if (length(keep) == 0) {
    return(integer(0))
  }
  ord <- keep[order(abs(cor_vec[keep]), decreasing = TRUE)]
  ord[seq_len(min(top_n, length(ord)))]
}

simplify_species <- function(x) {
  raw <- sub("^tax_s__", "", x)
  parts <- strsplit(raw, "_", fixed = TRUE)[[1]]
  if (length(parts) >= 2) {
    paste0(substr(parts[1], 1, 1), ". ", parts[2])
  } else {
    raw
  }
}

simplify_met <- function(x) {
  sub("^met_", "", x)
}

ko_labels <- paste0(
  sub("^kegg_", "", ko_cols_needed),
  " | ",
  unname(selected_ko_desc[sub("^kegg_", "", ko_cols_needed)])
)
names(ko_labels) <- ko_cols_needed

links_sp_ko <- data.frame(
  source_name = character(0),
  target_name = character(0),
  value = numeric(0),
  sign = character(0),
  stringsAsFactors = FALSE
)

for (ko_col in ko_cols_needed) {
  v <- cor_sp_ko[, ko_col]
  pick_idx <- pick_top_edges(v, rownames(cor_sp_ko), N_SPECIES_PER_KO, COR_THRESHOLD)
  if (length(pick_idx) == 0) next

  block <- data.frame(
    source_name = vapply(rownames(cor_sp_ko)[pick_idx], simplify_species, character(1)),
    target_name = ko_labels[[ko_col]],
    value = abs(v[pick_idx]),
    sign = ifelse(v[pick_idx] >= 0, "positive", "negative"),
    stringsAsFactors = FALSE
  )
  links_sp_ko <- rbind(links_sp_ko, block)
}

links_ko_met <- data.frame(
  source_name = character(0),
  target_name = character(0),
  value = numeric(0),
  sign = character(0),
  stringsAsFactors = FALSE
)

for (ko_col in ko_cols_needed) {
  v <- cor_ko_met[ko_col, ]
  pick_idx <- pick_top_edges(v, colnames(cor_ko_met), N_METABOLITES_PER_KO, COR_THRESHOLD)
  if (length(pick_idx) == 0) next

  met_names <- colnames(cor_ko_met)[pick_idx]
  block <- data.frame(
    source_name = ko_labels[[ko_col]],
    target_name = vapply(met_names, simplify_met, character(1)),
    value = abs(v[pick_idx]),
    sign = ifelse(v[pick_idx] >= 0, "positive", "negative"),
    stringsAsFactors = FALSE
  )
  links_ko_met <- rbind(links_ko_met, block)
}

links_all <- rbind(links_sp_ko, links_ko_met)
if (nrow(links_all) == 0) {
  stop("No links were selected. Try lowering COR_THRESHOLD.")
}

# Remove duplicate edges by keeping the strongest one.
links_all$key <- paste(links_all$source_name, links_all$target_name, sep = "__")
links_all <- links_all[order(links_all$key, -links_all$value), ]
links_all <- links_all[!duplicated(links_all$key), ]
links_all$key <- NULL

if (nrow(links_sp_ko) == 0 || nrow(links_ko_met) == 0) {
  stop("Insufficient links for three-layer Sankey. Please lower COR_THRESHOLD or increase top-N.")
}

species_nodes <- unique(links_sp_ko$source_name)
ko_nodes <- unique(c(links_sp_ko$target_name, links_ko_met$source_name))
met_nodes <- unique(links_ko_met$target_name)

nodes <- data.frame(
  name = c(species_nodes, ko_nodes, met_nodes),
  group = c(
    rep("Species", length(species_nodes)),
    rep("KO", length(ko_nodes)),
    rep("Metabolite", length(met_nodes))
  ),
  stringsAsFactors = FALSE
)

# Build Species -> KO -> Metabolite triplets through each KO.
triplets <- data.frame(
  species = character(0),
  ko = character(0),
  metabolite = character(0),
  value = numeric(0),
  ko_id = character(0),
  stringsAsFactors = FALSE
)

for (ko_node in ko_nodes) {
  sp_block <- links_sp_ko[links_sp_ko$target_name == ko_node, , drop = FALSE]
  met_block <- links_ko_met[links_ko_met$source_name == ko_node, , drop = FALSE]
  if (nrow(sp_block) == 0 || nrow(met_block) == 0) next

  for (i in seq_len(nrow(sp_block))) {
    for (j in seq_len(nrow(met_block))) {
      triplets <- rbind(
        triplets,
        data.frame(
          species = sp_block$source_name[i],
          ko = ko_node,
          metabolite = met_block$target_name[j],
          value = sqrt(sp_block$value[i] * met_block$value[j]),
          ko_id = sub(" \\|.*$", "", ko_node),
          stringsAsFactors = FALSE
        )
      )
    }
  }
}

if (nrow(triplets) == 0) {
  stop("No Sankey triplets generated. Please relax edge filtering parameters.")
}

triplets <- aggregate(value ~ species + ko + metabolite + ko_id, data = triplets, FUN = sum)

ko_palette <- c(
  K12688 = "#A61C3C",
  K02548 = "#D96C06",
  K07091 = "#E38F2D",
  K10200 = "#2E6F95",
  K10201 = "#4F9DA6",
  K02114 = "#6CA6C1",
  K01512 = "#4C956C",
  K01647 = "#8AAE5D",
  K02804 = "#7A5EA6",
  K03785 = "#A85D8E"
)

node_names <- c(species_nodes, ko_nodes, met_nodes)
species_palette <- setNames(alpha(hue_pal(l = 80, c = 40)(length(species_nodes)), 0.42), species_nodes)
ko_palette_nodes <- setNames(hue_pal(l = 52, c = 95)(length(ko_nodes)), ko_nodes)
met_palette <- setNames(alpha(hue_pal(l = 82, c = 35, h.start = 200)(length(met_nodes)), 0.42), met_nodes)
node_palette <- c(species_palette, ko_palette_nodes, met_palette)

p <- ggplot(
  triplets,
  aes(
    axis1 = species,
    axis2 = ko,
    axis3 = metabolite,
    y = value
  )
) +
  geom_alluvium(aes(fill = ko_id), width = 0.12, alpha = 0.78, knot.pos = 0.35, decreasing = FALSE) +
  scale_fill_manual(values = ko_palette, guide = "none") +
  ggnewscale::new_scale_fill() +
  geom_stratum(aes(fill = after_stat(stratum)), width = 0.12, color = "white", linewidth = 0.6) +
  scale_fill_manual(values = node_palette, guide = "none") +
  geom_text(
    stat = "stratum",
    aes(label = after_stat(stratum), fontface = after_stat(ifelse(stratum %in% ko_nodes, "bold", "plain"))),
    size = 2.4,
    family = "sans",
    lineheight = 0.95,
    color = "#1F1F1F"
  ) +
  scale_x_discrete(limits = c("Species", "KO", "Metabolite"), expand = c(0.08, 0.08)) +
  labs(
    title = "Species → KO → Metabolite Sankey",
    subtitle = "Fixed species, fixed KO panel, and fixed metabolite set",
    x = NULL,
    y = "Association Weight"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    plot.title = element_text(size = 15, face = "bold", color = "#111111", hjust = 0.5),
    plot.subtitle = element_text(size = 8.5, color = "#4D4D4D", hjust = 0.5, margin = margin(b = 5)),
    axis.text.x = element_text(size = 9.5, face = "bold", color = "#222222"),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "#EDEDED", linewidth = 0.25),
    plot.margin = margin(12, 18, 12, 18),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )

ggsave(
  filename = file.path(output_dir, "sankey_species_ko_metabolite.png"),
  plot = p,
  width = 28,
  height = 8.5,
  dpi = 800,
  bg = "white"
)

ggsave(
  filename = file.path(output_dir, "sankey_species_ko_metabolite.pdf"),
  plot = p,
  width = 28,
  height = 8.5,
  device = cairo_pdf,
  bg = "white"
)

write.csv(nodes, file.path(output_dir, "sankey_nodes.csv"), row.names = FALSE)
write.csv(links_all,
          file.path(output_dir, "sankey_links.csv"), row.names = FALSE)
write.csv(triplets, file.path(output_dir, "sankey_triplets.csv"), row.names = FALSE)

message("Done: ", output_dir)
