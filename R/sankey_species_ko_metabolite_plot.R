suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(ggalluvial)
  library(ggnewscale)
  library(scales)
})

input_file <- "dataset/merged_dataset_processed.csv"
output_dir <- "results/R_plots/sankey_species_ko_metabolite"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

ko_link_file <- "dataset/KEGG/cache/kegg_ko_pathway_links.csv"
pathway_meta_file <- "dataset/KEGG/cache/kegg_pathway_metadata.csv"

FIXED_KO_LIST <- c(
  "K12688",
  "K02548",
  "K07091",
  "K10200",
  "K10201",
  "K02114",
  "K01512",
  "K01647",
  "K02804",
  "K03785"
)

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

N_SPECIES_PER_KO <- 8
N_METABOLITES_PER_KO <- 8
COR_THRESHOLD <- 0.2
P_THRESHOLD <- 0.05
FLOW_WIDTH_MODE <- "effect_size"
EPS <- 1e-6
NODE_UP <- "#B2182B"
NODE_MID <- "#F7F7F2"
NODE_DOWN <- "#2166AC"
GROUP_ORDER <- c("control", "CRC_well_diff", "CRC_poor_diff")

if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file)
}

dat <- fread(input_file, data.table = FALSE)

nice_kegg_pathway <- function(x) {
  x <- as.character(x)
  x[is.na(x) | x == ""] <- "Unknown pathway"
  x
}

load_ko_pathway_label <- function(ko_ids) {
  if (!file.exists(ko_link_file) || !file.exists(pathway_meta_file)) {
    warning("KEGG cache 文件缺失，KO 标签将仅显示 KO 编号")
    out <- rep("Unknown pathway", length(ko_ids))
    names(out) <- ko_ids
    return(out)
  }

  ko_links <- fread(ko_link_file, data.table = FALSE)
  pathway_meta <- fread(pathway_meta_file, data.table = FALSE)

  pathway_meta$pathway_id_norm <- gsub("^ko", "map", pathway_meta$pathway_id)
  ko_links$pathway_id_norm <- gsub("^ko", "map", ko_links$pathway_id)

  merged <- ko_links %>%
    inner_join(pathway_meta, by = "pathway_id_norm") %>%
    filter(ko_id %in% ko_ids)

  ko_pathway <- merged %>%
    mutate(
      pathway_label = case_when(
        !is.na(level2) & level2 != "" ~ level2,
        !is.na(pathway_name) & pathway_name != "" ~ pathway_name,
        !is.na(level1) & level1 != "" ~ level1,
        TRUE ~ "Unknown pathway"
      )
    ) %>%
    group_by(ko_id, pathway_label) %>%
    summarise(n = n(), .groups = "drop") %>%
    arrange(ko_id, desc(n), pathway_label) %>%
    group_by(ko_id) %>%
    summarise(pathway_label = first(pathway_label), .groups = "drop")

  out <- rep("Unknown pathway", length(ko_ids))
  names(out) <- ko_ids
  if (nrow(ko_pathway) > 0) {
    out[ko_pathway$ko_id] <- ko_pathway$pathway_label
  }
  out
}

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

species_cols_all <- colnames(dat)[grepl("^tax_s__", colnames(dat))]
met_cols_all <- colnames(dat)[grepl("^met_", colnames(dat))]
ko_cols_needed <- paste0("kegg_", FIXED_KO_LIST)

species_cols <- match_targets_to_cols(FIXED_SPECIES_LIST, species_cols_all, "^tax_s__")
met_cols <- match_targets_to_cols(FIXED_METABOLITES_LIST, met_cols_all, "^met_")

missing_ko <- setdiff(ko_cols_needed, colnames(dat))
if (length(missing_ko) > 0) {
  stop("These selected KO columns are missing in dataset: ", paste(missing_ko, collapse = ", "))
}

if (length(species_cols) == 0 || length(met_cols) == 0) {
  stop("No matched fixed species or metabolites in dataset.")
}

N_SPECIES_PER_KO <- min(N_SPECIES_PER_KO, length(species_cols))
N_METABOLITES_PER_KO <- min(N_METABOLITES_PER_KO, length(met_cols))

dat$group <- with(dat, ifelse(
  crc_label == 0, "control",
  ifelse(differentiation == 1, "CRC_poor_diff", "CRC_well_diff")
))

present_groups <- unique(dat$group)
present_groups <- GROUP_ORDER[GROUP_ORDER %in% present_groups]
if (length(present_groups) == 0) {
  stop("No valid groups found. Expected: control / CRC_well_diff / CRC_poor_diff")
}

message("Group sample sizes:")
print(table(dat$group))

pairwise_spearman_stats <- function(mat_x, mat_y, x_type, y_type) {
  x_names <- colnames(mat_x)
  y_names <- colnames(mat_y)
  grid <- expand.grid(x = x_names, y = y_names, stringsAsFactors = FALSE)

  test_one <- function(xx, yy) {
    vx <- mat_x[, xx]
    vy <- mat_y[, yy]
    ok <- is.finite(vx) & is.finite(vy)
    if (sum(ok) < 5 || sd(vx[ok]) == 0 || sd(vy[ok]) == 0) {
      return(c(r = NA_real_, p = NA_real_))
    }
    ct <- suppressWarnings(cor.test(vx[ok], vy[ok], method = "spearman", exact = FALSE))
    c(r = unname(ct$estimate), p = ct$p.value)
  }

  stats <- purrr::map2_dfr(grid$x, grid$y, function(a, b) {
    vv <- test_one(a, b)
    tibble(from = a, to = b, r = vv[["r"]], p = vv[["p"]])
  })

  stats %>%
    mutate(
      fdr = p.adjust(p, method = "fdr"),
      edge_sign = ifelse(r >= 0, "positive", "negative"),
      edge_weight = abs(r),
      neglog10_p = -log10(pmax(p, EPS)),
      neglog10_fdr = -log10(pmax(fdr, EPS)),
      from_type = x_type,
      to_type = y_type
    )
}

select_top_edges <- function(stats_df, top_n, r_threshold, p_threshold) {
  stats_df %>%
    filter(!is.na(r), !is.na(p), abs(r) > r_threshold, p < p_threshold) %>%
    arrange(desc(abs(r)), p) %>%
    slice_head(n = top_n)
}

choose_flow_width <- function(r, p, fdr, mode = c("effect_size", "neglog10_p", "neglog10_fdr")) {
  mode <- match.arg(mode)
  if (mode == "effect_size") {
    return(abs(r))
  }
  if (mode == "neglog10_p") {
    return(-log10(pmax(p, EPS)))
  }
  -log10(pmax(fdr, EPS))
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

wrap_node_label <- function(x, type) {
  if (type == "KO") {
    # 针对 KO，将 "Kxxxxx | PathwayName" 拆分为两行
    if (grepl(" \\| ", x)) {
      parts <- strsplit(x, " \\| ")[[1]]
      ko_id <- parts[1]
      pw_name <- parts[2]
      # 对 pathway 名称进行换行处理，然后合并
      wrapped_pw <- str_wrap(pw_name, width = 20)
      return(paste0(ko_id, "\n", wrapped_pw))
    }
    return(str_wrap(x, width = 20))
  }
  if (type == "species") {
    return(str_wrap(x, width = 16))
  }
  # 代谢物名称通常较长，缩窄宽度强制换行以留在节点内
  str_wrap(x, width = 15, whitespace_only = FALSE)
}

build_node_fc_table <- function(dat, present_groups, species_cols, ko_cols_needed, met_cols, ko_labels) {
  ctrl_df <- dat[dat$group == "control", , drop = FALSE]
  ctrl_mean <- colMeans(ctrl_df[, c(species_cols, ko_cols_needed, met_cols), drop = FALSE], na.rm = TRUE)

  out <- data.frame(
    group = character(0),
    stratum = character(0),
    type = character(0),
    log2_fc = numeric(0),
    stringsAsFactors = FALSE
  )

  for (g in present_groups) {
    grp_df <- dat[dat$group == g, , drop = FALSE]
    grp_mean <- colMeans(grp_df[, c(species_cols, ko_cols_needed, met_cols), drop = FALSE], na.rm = TRUE)

    sp_df <- data.frame(
      group = g,
      stratum = vapply(species_cols, simplify_species, character(1)),
      type = "species",
      log2_fc = log2((grp_mean[species_cols] + EPS) / (ctrl_mean[species_cols] + EPS)),
      stringsAsFactors = FALSE
    )
    ko_df <- data.frame(
      group = g,
      stratum = unname(ko_labels[ko_cols_needed]),
      type = "KO",
      log2_fc = log2((grp_mean[ko_cols_needed] + EPS) / (ctrl_mean[ko_cols_needed] + EPS)),
      stringsAsFactors = FALSE
    )
    met_df <- data.frame(
      group = g,
      stratum = vapply(met_cols, simplify_met, character(1)),
      type = "metabolite",
      log2_fc = log2((grp_mean[met_cols] + EPS) / (ctrl_mean[met_cols] + EPS)),
      stringsAsFactors = FALSE
    )

    out <- rbind(out, sp_df, ko_df, met_df)
  }

  out$log2_fc[!is.finite(out$log2_fc)] <- 0
  out
}

ko_ids_needed <- sub("^kegg_", "", ko_cols_needed)
ko_pathway_map <- load_ko_pathway_label(ko_ids_needed)
ko_pathway_only <- nice_kegg_pathway(unname(ko_pathway_map[ko_ids_needed]))
ko_labels <- paste0(
  ko_ids_needed,
  " | ",
  ko_pathway_only
)
names(ko_labels) <- ko_cols_needed

ko_pathway_by_id <- setNames(ko_pathway_only, ko_ids_needed)
pathway_levels <- unique(ko_pathway_only)
pathway_legend_levels <- unique(c(pathway_levels, "Unknown pathway"))

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

pathway_palette <- setNames(unname(ko_palette[seq_along(pathway_levels)]), pathway_levels)

species_levels <- vapply(FIXED_SPECIES_LIST, simplify_species, character(1))
ko_levels <- unname(ko_labels[ko_cols_needed])
met_levels <- vapply(FIXED_METABOLITES_LIST, simplify_met, character(1))
node_order_levels <- c(species_levels, ko_levels, met_levels)

node_fc_tbl <- build_node_fc_table(
  dat = dat,
  present_groups = present_groups,
  species_cols = species_cols,
  ko_cols_needed = ko_cols_needed,
  met_cols = met_cols,
  ko_labels = ko_labels
)

links_sp_ko <- data.frame(
  source_name = character(0),
  target_name = character(0),
  value = numeric(0),
  r = numeric(0),
  p = numeric(0),
  fdr = numeric(0),
  sign = character(0),
  group = character(0),
  stringsAsFactors = FALSE
)

links_ko_met <- data.frame(
  source_name = character(0),
  target_name = character(0),
  value = numeric(0),
  r = numeric(0),
  p = numeric(0),
  fdr = numeric(0),
  sign = character(0),
  group = character(0),
  stringsAsFactors = FALSE
)

for (g in present_groups) {
  sub <- dat[dat$group == g, , drop = FALSE]

  species_mat <- as.matrix(sub[, species_cols, drop = FALSE])
  ko_mat <- as.matrix(sub[, ko_cols_needed, drop = FALSE])
  met_mat <- as.matrix(sub[, met_cols, drop = FALSE])

  mode(species_mat) <- "numeric"
  mode(ko_mat) <- "numeric"
  mode(met_mat) <- "numeric"

  stats_sp_ko <- pairwise_spearman_stats(species_mat, ko_mat, x_type = "species", y_type = "KO")
  stats_ko_met <- pairwise_spearman_stats(ko_mat, met_mat, x_type = "KO", y_type = "metabolite")

  for (ko_col in ko_cols_needed) {
    ko_stats <- stats_sp_ko %>%
      filter(to == ko_col) %>%
      select_top_edges(top_n = N_SPECIES_PER_KO, r_threshold = COR_THRESHOLD, p_threshold = P_THRESHOLD)
    if (nrow(ko_stats) == 0) next

    block <- ko_stats %>%
      transmute(
        source_name = vapply(from, simplify_species, character(1)),
        target_name = ko_labels[[ko_col]],
        value = choose_flow_width(r, p, fdr, FLOW_WIDTH_MODE),
        r = r,
        p = p,
        fdr = fdr,
        sign = edge_sign,
        group = g
      )
    links_sp_ko <- rbind(links_sp_ko, block)
  }

  for (ko_col in ko_cols_needed) {
    ko_stats <- stats_ko_met %>%
      filter(from == ko_col) %>%
      select_top_edges(top_n = N_METABOLITES_PER_KO, r_threshold = COR_THRESHOLD, p_threshold = P_THRESHOLD)
    if (nrow(ko_stats) == 0) next

    block <- ko_stats %>%
      transmute(
        source_name = ko_labels[[ko_col]],
        target_name = vapply(to, simplify_met, character(1)),
        value = choose_flow_width(r, p, fdr, FLOW_WIDTH_MODE),
        r = r,
        p = p,
        fdr = fdr,
        sign = edge_sign,
        group = g
      )
    links_ko_met <- rbind(links_ko_met, block)
  }
}

links_all <- rbind(links_sp_ko, links_ko_met)
if (nrow(links_all) == 0) {
  stop("No links were selected. Try lowering COR_THRESHOLD or relaxing P_THRESHOLD.")
}

# Keep only complete Species -> KO -> Metabolite paths within each group.
valid_ko_group <- merge(
  unique(links_sp_ko[, c("group", "target_name"), drop = FALSE]),
  unique(links_ko_met[, c("group", "source_name"), drop = FALSE]),
  by.x = c("group", "target_name"),
  by.y = c("group", "source_name")
)

if (nrow(valid_ko_group) == 0) {
  stop("No complete Species->KO->Metabolite paths found under current thresholds.")
}

links_sp_ko <- links_sp_ko %>%
  inner_join(valid_ko_group, by = c("group", "target_name"))

links_ko_met <- links_ko_met %>%
  inner_join(valid_ko_group, by = c("group", "source_name" = "target_name"))

links_all <- rbind(links_sp_ko, links_ko_met)
if (nrow(links_all) == 0) {
  stop("No complete-path links remain after filtering.")
}

message("Edge counts after complete-path filtering:")
for (g in present_groups) {
  n_sp_ko <- nrow(links_sp_ko[links_sp_ko$group == g, , drop = FALSE])
  n_ko_met <- nrow(links_ko_met[links_ko_met$group == g, , drop = FALSE])
  message(sprintf("  %s: species->KO=%d, KO->metabolite=%d", g, n_sp_ko, n_ko_met))
}

# Remove duplicate edges by keeping the strongest one.
links_all$key <- paste(links_all$group, links_all$source_name, links_all$target_name, sep = "__")
links_all <- links_all[order(links_all$key, -links_all$value), ]
links_all <- links_all[!duplicated(links_all$key), ]
links_all$key <- NULL

if (nrow(links_sp_ko) == 0 || nrow(links_ko_met) == 0) {
  stop("Insufficient links for three-layer Sankey. Please lower COR_THRESHOLD / P_THRESHOLD or increase top-N.")
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

links_sp_ko_balanced <- links_sp_ko[0, , drop = FALSE]
links_ko_met_balanced <- links_ko_met[0, , drop = FALSE]

for (ko_node in ko_nodes) {
  for (g in present_groups) {
    sp_block <- links_sp_ko[links_sp_ko$target_name == ko_node & links_sp_ko$group == g, , drop = FALSE]
    met_block <- links_ko_met[links_ko_met$source_name == ko_node & links_ko_met$group == g, , drop = FALSE]
    if (nrow(sp_block) == 0 || nrow(met_block) == 0) next

    # Balance both sides to a shared KO throughput so flow widths match node heights.
    sp_vals <- pmax(sp_block$value, EPS)
    met_vals <- pmax(met_block$value, EPS)
    sp_total <- sum(sp_vals)
    met_total <- sum(met_vals)
    shared_total <- min(sp_total, met_total)
    if (!is.finite(shared_total) || shared_total <= 0) next

    sp_scaled <- sp_vals / sp_total * shared_total
    met_scaled <- met_vals / met_total * shared_total

    sp_block$value <- sp_scaled
    met_block$value <- met_scaled
    links_sp_ko_balanced <- rbind(links_sp_ko_balanced, sp_block)
    links_ko_met_balanced <- rbind(links_ko_met_balanced, met_block)

    flow_mat <- outer(sp_scaled, met_scaled, FUN = "*") / shared_total

    for (i in seq_len(nrow(sp_block))) {
      for (j in seq_len(nrow(met_block))) {
        flow_val <- flow_mat[i, j]
        if (!is.finite(flow_val) || flow_val <= 0) next
        triplets <- rbind(
          triplets,
          data.frame(
            species = sp_block$source_name[i],
            ko = ko_node,
            metabolite = met_block$target_name[j],
            value = flow_val,
            ko_id = sub(" \\|.*$", "", ko_node),
            group = g,
            stringsAsFactors = FALSE
          )
        )
      }
    }
  }
}

links_sp_ko <- links_sp_ko_balanced
links_ko_met <- links_ko_met_balanced
links_all <- rbind(links_sp_ko, links_ko_met)

if (nrow(triplets) == 0) {
  stop("No Sankey triplets generated. Please relax edge filtering parameters.")
}

triplets <- aggregate(value ~ group + species + ko + metabolite + ko_id, data = triplets, FUN = sum)
triplets$flow_id <- seq_len(nrow(triplets))

triplets_lodes <- to_lodes_form(
  triplets,
  axes = c("species", "ko", "metabolite"),
  id = "flow_id",
  key = "axis",
  value = "stratum"
)

triplets_lodes$ko_id <- triplets$ko_id[match(triplets_lodes$flow_id, triplets$flow_id)]

if (any(is.na(triplets_lodes$ko_id))) {
  stop("Failed to map ko_id onto lodes data.")
}

triplets_lodes$axis <- factor(triplets_lodes$axis, levels = c("species", "ko", "metabolite"), labels = c("Species", "KO", "Metabolite"))
triplets_lodes$type <- dplyr::case_when(
  triplets_lodes$axis == "Species" ~ "species",
  triplets_lodes$axis == "KO" ~ "KO",
  TRUE ~ "metabolite"
)

triplets_lodes$stratum <- factor(triplets_lodes$stratum, levels = node_order_levels)
triplets_lodes$pathway <- unname(ko_pathway_by_id[triplets_lodes$ko_id])
triplets_lodes$pathway[is.na(triplets_lodes$pathway) | triplets_lodes$pathway == ""] <- "Unknown pathway"
triplets_lodes$pathway <- factor(triplets_lodes$pathway, levels = pathway_legend_levels)
present_pathways <- unique(as.character(triplets_lodes$pathway))
present_pathways <- present_pathways[!is.na(present_pathways) & present_pathways != ""]
present_pathways <- intersect(pathway_legend_levels, present_pathways)

triplets_lodes$label <- mapply(function(s, t) wrap_node_label(as.character(s), as.character(t)), 
                               triplets_lodes$stratum, triplets_lodes$type, USE.NAMES = FALSE)
triplets_lodes$fontface <- ifelse(triplets_lodes$type == "species", "italic", "plain")

triplets_lodes <- triplets_lodes %>%
  left_join(node_fc_tbl, by = c("group", "stratum", "type"))
triplets_lodes$log2_fc[is.na(triplets_lodes$log2_fc)] <- 0

fc_lim <- max(abs(triplets_lodes$log2_fc), na.rm = TRUE)
if (!is.finite(fc_lim) || fc_lim <= 0) fc_lim <- 1

p <- ggplot(
  triplets_lodes,
  aes(
    x = axis,
    stratum = stratum,
    alluvium = flow_id,
    y = value
  )
) +
  geom_alluvium(aes(fill = pathway), width = 0.12, alpha = 0.46, knot.pos = 0.35, decreasing = FALSE) +
  scale_fill_manual(
    values = c(pathway_palette, "Unknown pathway" = "#BDBDBD"),
    breaks = present_pathways,
    drop = TRUE,
    name = "Pathway",
    guide = guide_legend(override.aes = list(alpha = 0.9))
  ) +
  ggnewscale::new_scale_fill() +
  geom_stratum(aes(fill = log2_fc), width = 0.12, color = "white", linewidth = 0.6, decreasing = FALSE) +
  scale_fill_gradient2(
    low = NODE_DOWN,
    mid = NODE_MID,
    high = NODE_UP,
    midpoint = 0,
    limits = c(-fc_lim, fc_lim),
    oob = squish,
    name = "Node log2FC",
    guide = guide_colourbar(barheight = unit(45, "mm"), barwidth = unit(4, "mm"))
  ) +
  geom_text(
    stat = "stratum",
    aes(label = label, fontface = fontface),
    decreasing = FALSE,
    size = 2.2,
    family = "sans",
    lineheight = 0.9,
    color = "#1F1F1F"
  ) +
  scale_x_discrete(limits = c("Species", "KO", "Metabolite"), expand = c(0.08, 0.08)) +
  labs(
    title = "Species -> KO -> Metabolite Sankey by Group",
    subtitle = "Node color: red enriched vs CTRL, blue depleted vs CTRL; flow width follows selected effect metric",
    x = NULL,
    y = ifelse(FLOW_WIDTH_MODE == "effect_size", "Association strength", "-log10(p-value)")
  ) +
  facet_wrap(~group, ncol = 1, scales = "free_y") +
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
    plot.background = element_rect(fill = "white", color = NA),
    strip.text = element_text(size = 10, face = "bold", color = "#222222"),
    strip.background = element_rect(fill = "#F3F5F7", color = NA)
  )

ggsave(
  filename = file.path(output_dir, "sankey_species_ko_metabolite.png"),
  plot = p,
  width = 28,
  height = 20,
  dpi = 800,
  bg = "white"
)

ggsave(
  filename = file.path(output_dir, "sankey_species_ko_metabolite.pdf"),
  plot = p,
  width = 28,
  height = 20,
  device = cairo_pdf,
  bg = "white"
)

for (g in present_groups) {
  g_lodes <- triplets_lodes[triplets_lodes$group == g, , drop = FALSE]
  if (nrow(g_lodes) == 0) {
    warning(sprintf("Skip plotting %s: no valid triplets after filtering (|r| > %.2f, p < %.2f)", g, COR_THRESHOLD, P_THRESHOLD))
    next
  }

  g_fc_lim <- max(abs(g_lodes$log2_fc), na.rm = TRUE)
  if (!is.finite(g_fc_lim) || g_fc_lim <= 0) g_fc_lim <- 1

  p_g <- ggplot(
    g_lodes,
    aes(
      x = axis,
      stratum = stratum,
      alluvium = flow_id,
      y = value
    )
  ) +
    geom_alluvium(aes(fill = pathway), width = 0.12, alpha = 0.46, knot.pos = 0.35, decreasing = FALSE) +
    scale_fill_manual(
      values = c(pathway_palette, "Unknown pathway" = "#BDBDBD"),
      breaks = present_pathways,
      drop = TRUE,
      name = "Pathway",
      guide = guide_legend(override.aes = list(alpha = 0.9))
    ) +
    ggnewscale::new_scale_fill() +
    geom_stratum(aes(fill = log2_fc), width = 0.12, color = "white", linewidth = 0.6, decreasing = FALSE) +
    scale_fill_gradient2(
      low = NODE_DOWN,
      mid = NODE_MID,
      high = NODE_UP,
      midpoint = 0,
      limits = c(-g_fc_lim, g_fc_lim),
      oob = squish,
      name = "Node log2FC",
      guide = guide_colourbar(barheight = unit(45, "mm"), barwidth = unit(4, "mm"))
    ) +
    geom_text(
      stat = "stratum",
      aes(label = label, fontface = fontface),
      decreasing = FALSE,
      size = 2.2,
      family = "sans",
      lineheight = 0.9,
      color = "#1F1F1F"
    ) +
    scale_x_discrete(limits = c("Species", "KO", "Metabolite"), expand = c(0.08, 0.08)) +
    labs(
      title = paste0("Species -> KO -> Metabolite Sankey (", g, ")"),
      subtitle = paste0("Group-specific flows, |r| > ", COR_THRESHOLD, ", p < ", P_THRESHOLD),
      x = NULL,
      y = ifelse(FLOW_WIDTH_MODE == "effect_size", "Association strength", "-log10(p-value)")
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
    filename = file.path(output_dir, paste0("sankey_species_ko_metabolite_", g, ".png")),
    plot = p_g,
    width = 28,
    height = 8.5,
    dpi = 800,
    bg = "white"
  )

  ggsave(
    filename = file.path(output_dir, paste0("sankey_species_ko_metabolite_", g, ".pdf")),
    plot = p_g,
    width = 28,
    height = 8.5,
    device = cairo_pdf,
    bg = "white"
  )
}

write.csv(nodes, file.path(output_dir, "sankey_nodes.csv"), row.names = FALSE)
write.csv(links_all,
          file.path(output_dir, "sankey_links.csv"), row.names = FALSE)
write.csv(triplets, file.path(output_dir, "sankey_triplets.csv"), row.names = FALSE)

for (g in present_groups) {
  write.csv(
    links_all[links_all$group == g, , drop = FALSE],
    file.path(output_dir, paste0("sankey_links_", g, ".csv")),
    row.names = FALSE
  )
  write.csv(
    triplets[triplets$group == g, , drop = FALSE],
    file.path(output_dir, paste0("sankey_triplets_", g, ".csv")),
    row.names = FALSE
  )
}

message("Done: ", output_dir)
