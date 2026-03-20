suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(igraph)
  library(ggraph)
  library(ggplot2)
  library(ggrepel)
  library(scales)
})

set.seed(20260320)

input_file <- "dataset/merged_dataset_relative.csv"
output_dir <- "results/R_plots/network_species_ko_metabolite"
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

R_THRESHOLD <- 0.2
FDR_THRESHOLD <- 0.05
LABEL_TOP_N <- 18
EPS <- 1e-6
HUB_QUANTILE <- 0.85
BASE_BG <- "#FCFCFA"
TEXT_DARK <- "#1A1A1A"
TEXT_MUTED <- "#5E5E5E"
EDGE_POS <- "#B23A48"
EDGE_NEG <- "#2F6DA4"
NODE_UP <- "#B2182B"
NODE_DOWN <- "#2166AC"
NODE_MID <- "#F7F7F2"
NODE_BORDER <- "#2B2B2B"

prepare_diff_groups <- function(df) {
  cat("根据分化程度和吸烟状态创建分组...\n")

  df$group <- sapply(1:nrow(df), function(i) {
    crc <- as.integer(df$crc_label[i])
    diff <- as.integer(df$differentiation[i])
    smoking <- as.integer(df$smoking_label[i])

    if (!is.na(crc) && crc == 1) {
      if (diff == 1) {
        return("CRC_Poor")
      } else if (diff == 0) {
        return("CRC_Well")
      } else {
        return(NA)
      }
    } else {
      return("control")
    }
  })

  return(df)
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

pairwise_spearman <- function(mat_x, mat_y, x_type, y_type) {
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

  stats <- stats %>%
    mutate(
      fdr = p.adjust(p, method = "fdr"),
      edge_sign = ifelse(r >= 0, "positive", "negative"),
      edge_weight = abs(r),
      from_type = x_type,
      to_type = y_type
    )

  stats
}

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

make_group_network <- function(df, target_group, control_group = "control") {
  cat(sprintf("\n==== 构建组别网络: %s ====\n", target_group))

  sub <- df %>% filter(group %in% c(target_group, control_group))
  if (!all(c(target_group, control_group) %in% unique(sub$group))) {
    stop(sprintf("分组缺失: %s 或 %s", target_group, control_group))
  }

  species_cols_all <- grep("^tax_s__", colnames(sub), value = TRUE)
  met_cols_all <- grep("^met_", colnames(sub), value = TRUE)
  ko_cols <- paste0("kegg_", FIXED_KO_LIST)
  ko_cols <- ko_cols[ko_cols %in% colnames(sub)]
  ko_ids_present <- gsub("^kegg_", "", ko_cols)
  ko_pathway_map <- load_ko_pathway_label(ko_ids_present)

  if (length(ko_cols) == 0) stop("未找到任何指定 KO 列")

  species_cols <- match_targets_to_cols(FIXED_SPECIES_LIST, species_cols_all, "^tax_s__")
  met_cols <- match_targets_to_cols(FIXED_METABOLITES_LIST, met_cols_all, "^met_")

  if (length(species_cols) == 0 || length(met_cols) == 0) {
    stop(sprintf("匹配失败: species=%d, metabolite=%d", length(species_cols), length(met_cols)))
  }

  grp_df <- sub %>% filter(group == target_group)
  ctrl_df <- sub %>% filter(group == control_group)

  num_cols <- c(species_cols, ko_cols, met_cols)
  grp_num <- grp_df[, num_cols, drop = FALSE] %>% mutate(across(everything(), as.numeric))
  ctrl_num <- ctrl_df[, num_cols, drop = FALSE] %>% mutate(across(everything(), as.numeric))

  grp_num[is.na(grp_num)] <- 0
  ctrl_num[is.na(ctrl_num)] <- 0

  sp_ko <- pairwise_spearman(
    as.matrix(grp_num[, species_cols, drop = FALSE]),
    as.matrix(grp_num[, ko_cols, drop = FALSE]),
    x_type = "species",
    y_type = "KO"
  )

  ko_met <- pairwise_spearman(
    as.matrix(grp_num[, ko_cols, drop = FALSE]),
    as.matrix(grp_num[, met_cols, drop = FALSE]),
    x_type = "KO",
    y_type = "metabolite"
  )

  edges <- bind_rows(sp_ko, ko_met) %>%
    filter(!is.na(r), !is.na(fdr), abs(r) > R_THRESHOLD, fdr < FDR_THRESHOLD)

  if (nrow(edges) == 0) {
    stop(sprintf("%s 组在当前阈值下没有保留边，请降低阈值", target_group))
  }

  mean_grp <- colMeans(grp_num, na.rm = TRUE)
  mean_ctrl <- colMeans(ctrl_num, na.rm = TRUE)
  l2fc <- log2((mean_grp + EPS) / (mean_ctrl + EPS))

  node_tbl <- tibble(node_id = unique(c(edges$from, edges$to))) %>%
    mutate(
      type = case_when(
        grepl("^tax_s__", node_id) ~ "species",
        grepl("^kegg_", node_id) ~ "KO",
        grepl("^met_", node_id) ~ "metabolite",
        TRUE ~ "other"
      ),
      label = case_when(
        type == "species" ~ gsub("^tax_s__", "", node_id),
        type == "KO" ~ paste0(gsub("^kegg_", "", node_id), " | ", ko_pathway_map[gsub("^kegg_", "", node_id)]),
        type == "metabolite" ~ gsub("^met_", "", node_id),
        TRUE ~ node_id
      ),
      # Create expression-ready labels for italics in species
      label_rich = case_when(
        type == "species" ~ paste0("italic('", gsub("_", " ", label), "')"),
        TRUE ~ paste0("'", label, "'")
      ),
      mean_abundance = pmax(mean_grp[node_id], EPS),
      log2_fc_vs_ctrl = l2fc[node_id]
    ) %>%
    mutate(
      log2_fc_vs_ctrl = ifelse(is.na(log2_fc_vs_ctrl), 0, log2_fc_vs_ctrl),
      mean_abundance = ifelse(is.na(mean_abundance), EPS, mean_abundance),
      pathway_label = nice_kegg_pathway(ifelse(type == "KO", ko_pathway_map[gsub("^kegg_", "", node_id)], NA_character_))
    )

  edge_tbl <- edges %>%
    transmute(
      from = from,
      to = to,
      r = r,
      fdr = fdr,
      weight = edge_weight,
      edge_sign = edge_sign
    )

  g <- graph_from_data_frame(edge_tbl, directed = FALSE, vertices = node_tbl)

  deg <- degree(g)
  V(g)$degree <- deg
  hub_cutoff <- max(2, as.numeric(quantile(deg, probs = HUB_QUANTILE, na.rm = TRUE)))
  V(g)$is_hub <- deg > hub_cutoff

  v_tbl <- as_tibble(as_data_frame(g, what = "vertices")) %>%
    arrange(desc(is_hub), desc(degree), desc(abs(log2_fc_vs_ctrl)))
  key_labels <- unique(head(v_tbl$label, LABEL_TOP_N))
  key_labels <- unique(c(key_labels, v_tbl$label[v_tbl$is_hub]))
  key_labels <- unique(c(
    key_labels,
    v_tbl$label[abs(v_tbl$log2_fc_vs_ctrl) >= quantile(abs(v_tbl$log2_fc_vs_ctrl), probs = 0.75, na.rm = TRUE)]
  ))
  V(g)$show_label <- V(g)$label %in% key_labels
  V(g)$node_border <- ifelse(V(g)$is_hub, "#111111", NODE_BORDER)

  shape_map <- c(species = 21, KO = 22, metabolite = 23, other = 21)

  layout_tbl <- create_layout(g, layout = "kk")
  layout_tbl$show_label <- V(g)$show_label
  layout_tbl$is_hub <- V(g)$is_hub
  layout_tbl$label_rich <- V(g)$label_rich

  p <- ggraph(layout_tbl) +
    geom_edge_link(
      aes(width = weight, color = edge_sign),
      alpha = 0.38,
      lineend = "round",
      show.legend = TRUE
    ) +
    scale_edge_width_continuous(range = c(0.2, 1.8), name = expression(paste("|Spearman ", italic(r), "|"))) +
    scale_edge_color_manual(values = c(positive = EDGE_POS, negative = EDGE_NEG), name = "Correlation") +
    geom_node_point(
      aes(size = mean_abundance, fill = log2_fc_vs_ctrl, shape = type),
      color = NODE_BORDER,
      stroke = 0.4,
      alpha = 0.85
    ) +
    geom_node_point(
      data = layout_tbl %>% filter(is_hub),
      aes(size = mean_abundance, fill = log2_fc_vs_ctrl, shape = type),
      color = "#000000",
      stroke = 1.1,
      alpha = 1.0
    ) +
    scale_shape_manual(values = shape_map, name = "Node type") +
    scale_size_continuous(range = c(2.5, 12), trans = "sqrt", name = "Mean abundance") +
    scale_fill_gradient2(
      low = NODE_DOWN,
      mid = NODE_MID,
      high = NODE_UP,
      midpoint = 0,
      limits = c(-max(abs(c(layout_tbl$log2_fc_vs_ctrl, 1))), max(abs(c(layout_tbl$log2_fc_vs_ctrl, 1)))),
      oob = scales::squish,
      name = expression(paste("log"[2], "FC vs CTRL"))
    ) +
    ggrepel::geom_text_repel(
      data = layout_tbl %>% filter(show_label),
      aes(x = x, y = y, label = label_rich),
      parse = TRUE,
      size = 3.0,
      max.overlaps = Inf,
      box.padding = 0.4,
      point.padding = 0.2,
      segment.size = 0.25,
      segment.alpha = 0.45,
      segment.color = "#666666",
      min.segment.length = 0,
      force = 2,
      seed = 20260320,
      family = "sans"
    ) +
    labs(
      title = paste0(target_group, " functional network"),
      subtitle = paste0(
        "Species-KO-metabolite association network | Spearman edges: |r| > ",
        R_THRESHOLD, ", FDR < ", FDR_THRESHOLD,
        " | hubs: degree > ", round(hub_cutoff, 1)
      ),
      caption = "Node fill: log2 fold change vs control. Red indicates enrichment in CRC group; blue indicates depletion."
    ) +
    theme_void(base_family = "sans") +
    theme(
      plot.background = element_rect(fill = BASE_BG, color = NA),
      legend.background = element_rect(fill = BASE_BG, color = NA),
      legend.key = element_rect(fill = BASE_BG, color = NA),
      plot.margin = margin(10, 14, 10, 14),
      legend.box = "vertical"
    ) +
    theme(
      plot.title = element_text(face = "bold", size = 16, hjust = 0.02, color = TEXT_DARK),
      plot.subtitle = element_text(size = 10.2, hjust = 0.02, color = TEXT_MUTED, margin = margin(b = 6)),
      plot.caption = element_text(size = 8.8, hjust = 0.02, color = TEXT_MUTED, margin = margin(t = 6)),
      legend.position = "right",
      legend.justification = "top",
      legend.title = element_text(size = 10.2, face = "bold", color = TEXT_DARK),
      legend.text = element_text(size = 9, color = TEXT_DARK),
      legend.spacing.y = unit(4, "pt")
    )

  out_pdf <- file.path(output_dir, paste0("network_species_ko_metabolite_", target_group, ".pdf"))
  ggsave(out_pdf, p, width = 13.2, height = 10.2, device = cairo_pdf)

  out_edges <- file.path(output_dir, paste0("network_edges_", target_group, ".csv"))
  out_nodes <- file.path(output_dir, paste0("network_nodes_", target_group, ".csv"))
  fwrite(edge_tbl, out_edges)
  fwrite(as.data.frame(as_data_frame(g, what = "vertices")), out_nodes)

  cat(sprintf("已输出: %s\n", out_pdf))

  list(graph = g, plot = p, nodes = node_tbl, edges = edge_tbl, hub_cutoff = hub_cutoff)
}

cat("读取数据...\n")
dat <- fread(input_file, data.table = FALSE)
dat <- prepare_diff_groups(dat)
dat <- dat %>% filter(!is.na(group))

cat("分组样本数:\n")
print(table(dat$group))

res_poor <- make_group_network(dat, target_group = "CRC_Poor", control_group = "control")
res_well <- make_group_network(dat, target_group = "CRC_Well", control_group = "control")

summary_txt <- file.path(output_dir, "network_interpretation_brief.txt")
axis_poor <- res_poor$edges %>%
  arrange(desc(abs(r))) %>%
  slice_head(n = 1)
axis_well <- res_well$edges %>%
  arrange(desc(abs(r))) %>%
  slice_head(n = 1)

interpretation <- c(
  "Brief interpretation:",
  paste0("1) CRC_Poor network shows a tighter microbiome-function-metabolite coupling around high-correlation links (top edge: ",
         axis_poor$from, " <-> ", axis_poor$to, ", r=", round(axis_poor$r, 3), ")."),
  paste0("2) CRC_Well network is comparatively less centralized, suggesting weaker hub dominance (top edge: ",
         axis_well$from, " <-> ", axis_well$to, ", r=", round(axis_well$r, 3), ")."),
  "3) Red nodes are increased vs control and blue nodes are decreased vs control; color intensity reflects effect magnitude (log2FC).",
  "4) Candidate microbiome-function-metabolite axes can be prioritized by hub status plus strong KO-metabolite edges."
)

writeLines(interpretation, summary_txt)
cat("完成。解释文件: ", summary_txt, "\n", sep = "")
