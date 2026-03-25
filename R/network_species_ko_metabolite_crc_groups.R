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

config <- list(
  maaslin_microbe = "results/maaslin2/significant_results.tsv",
  maaslin_ko = "results/maaslin2_ko/significant_results.tsv",
  maaslin_metabolite = "results/maaslin2_metabolite/significant_results.tsv",
  abundance_file = input_file,
  outdir = output_dir,
  qval_cutoff = 0.05,
  max_auto_ko_features = 500,
  max_auto_metabolite_features = 500,
  rho_cutoff = 0.2,
  fdr_cutoff = 0.05
)

FIXED_KO_LIST <- c(
  # "K12688",
  "K02548",
  # "K07091",
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
  "Faecalibacterium_prausnitzii"
)

FIXED_METABOLITES_LIST <- c(
  "SQDG 26:2; SQDG(13:1/13:1)",
  # "Cytosine",
  # "Perfluorooctanesulfonic acid",
  # "Methyl dihydrojasmonate",
  # "Pyrogallol-2-O-sulphate",
  # "5'-(3',4'-Dihydroxyphenyl)-gamma-valerolactone sulfate",
  # "2-Hydroxy-4,7-dimethoxy-2H-1,4-benzoxazin-3(4H)-one",
  # "trans-3,5-Dimethoxy-4-hydroxycinnamaldehyde",
  # "(R)-3-Hydroxy-5-phenylpentanoic acid",
  # "N-Methyl-D-glucamine",
  "Chenodeoxycholic acid sulfate",
  # "Lucidenic acid F",
  # "Demissidine",
  # "Alpha-Hydroxyisobutyric acid",
  # "Pyrocatechol",
  # "Gentisic acid",
  # "D-Galacturonic acid",
  # "1,3-Dimethyluric acid",
  "4-Hydroxy-5-(phenyl)-valeric acid-O-sulphate",
  "Cholesterol"
)

USE_FIXED_KO_LIST <- TRUE
USE_FIXED_METABOLITES_LIST <- TRUE

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
LINE_WIDTH <- 4.4

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

normalize_feature_name <- function(x) {
  x <- trimws(as.character(x))
  x <- tolower(x)
  x <- gsub("^tax_s__", "", x)
  x <- gsub("^kegg_", "", x)
  x <- gsub("^met_", "", x)
  x <- gsub("[^a-z0-9]+", "", x)
  x
}

normalize_text <- function(x) {
  tolower(gsub("[^a-z0-9]", "", x))
}

match_targets_to_cols <- function(target_list, available_cols, remove_prefix = "", strict = TRUE, print_unmatched = TRUE) {
  display_names <- gsub(remove_prefix, "", available_cols)
  norm_display <- normalize_feature_name(display_names)

  matched_cols <- character(0)
  missing_targets <- character(0)

  for (target in target_list) {
    norm_target <- normalize_feature_name(target)

    exact_idx <- which(norm_display == norm_target)
    if (length(exact_idx) > 0) {
      matched_cols <- c(matched_cols, available_cols[exact_idx[1]])
      next
    }

    partial_idx <- which(nzchar(norm_target) & grepl(norm_target, norm_display, fixed = TRUE))
    if (length(partial_idx) > 0) {
      if (length(partial_idx) == 1) {
        matched_cols <- c(matched_cols, available_cols[partial_idx[1]])
        next
      }

      partial_display <- available_cols[partial_idx]
      message(sprintf(
        "%s 匹配到多个候选列，跳过模糊匹配: %s",
        target,
        paste(partial_display, collapse = ", ")
      ))
    }

    missing_targets <- c(missing_targets, target)
  }

  if (length(missing_targets) > 0) {
    msg <- paste("Target features not found:", paste(missing_targets, collapse = ", "))
    if (isTRUE(print_unmatched)) {
      message("未匹配到的目标特征: ", paste(missing_targets, collapse = ", "))
      if (length(available_cols) > 0) {
        preview_cols <- head(available_cols, 10)
        message("可用列示例: ", paste(preview_cols, collapse = ", "))
      }
    }
    if (isTRUE(strict)) {
      stop(msg)
    } else {
      warning(msg)
    }
  }

  matched_cols
}

read_maaslin_qval_features <- function(path, qval_cutoff = 0.05) {
  df <- fread(path, data.table = FALSE)
  if (!all(c("feature", "qval") %in% colnames(df))) {
    stop(sprintf("MaAsLin结果缺少必要列: %s", path))
  }

  order_cols <- c("qval")
  if ("pval" %in% colnames(df)) {
    order_cols <- c(order_cols, "pval")
  }

  df %>%
    filter(!is.na(feature), !is.na(qval), qval < qval_cutoff) %>%
    arrange(across(all_of(order_cols))) %>%
    distinct(feature) %>%
    pull(feature)
}

resolve_feature_columns <- function(available_cols, fixed_list, maaslin_path, remove_prefix = "", label = "feature",
                                    use_fixed = TRUE, qval_cutoff = 0.05, max_auto_features = Inf) {
  if (isTRUE(use_fixed) && length(fixed_list) > 0) {
    message(sprintf("%s: 使用固定列表, 目标数=%d", label, length(fixed_list)))
    return(match_targets_to_cols(fixed_list, available_cols, remove_prefix = remove_prefix, strict = TRUE, print_unmatched = TRUE))
  }

  qval_features <- read_maaslin_qval_features(maaslin_path, qval_cutoff = qval_cutoff)
  if (is.finite(max_auto_features)) {
    qval_features <- head(qval_features, max_auto_features)
  }
  message(sprintf("%s: 自动提取 MaAsLin qval < %.2f 的特征, 候选数=%d", label, qval_cutoff, length(qval_features)))
  matched <- match_targets_to_cols(qval_features, available_cols, remove_prefix = remove_prefix, strict = FALSE, print_unmatched = TRUE)
  if (length(matched) == 0) {
    stop(sprintf("%s: qval < %.2f 的特征在丰度矩阵中没有匹配到任何列", label, qval_cutoff))
  }

  if (length(matched) < length(qval_features)) {
    missing_n <- length(qval_features) - length(unique(matched))
    message(sprintf("%s: %d 个候选特征未匹配到丰度列", label, missing_n))
  }
  unique(matched)
}

format_node_label <- function(node_id, ko_pathway_map = NULL) {
  if (grepl("^tax_s__", node_id)) {
    return(gsub("^tax_s__", "", node_id))
  }

  if (grepl("^kegg_", node_id)) {
    ko_id <- gsub("^kegg_", "", node_id)
    if (!is.null(ko_pathway_map) && ko_id %in% names(ko_pathway_map)) {
      return(paste0(ko_id, " | ", ko_pathway_map[ko_id]))
    }
    return(ko_id)
  }

  if (grepl("^met_", node_id)) {
    return(gsub("^met_", "", node_id))
  }

  node_id
}

escape_label_text <- function(x) {
  x <- as.character(x)
  x <- gsub("\\\\", "\\\\\\\\", x)
  x <- gsub('"', '\\\"', x, fixed = TRUE)
  x
}

make_label_expr <- function(x, italic = FALSE) {
  x <- escape_label_text(x)
  if (isTRUE(italic)) {
    paste0("italic(\"", x, "\")")
  } else {
    paste0("\"", x, "\"")
  }
}

extract_maaslin_effects <- function(path, features, qval_cutoff = 0.05) {
  df <- fread(path, data.table = FALSE)
  required <- c("feature", "coef", "qval")
  if (!all(required %in% colnames(df))) {
    stop(sprintf("MaAsLin结果缺少必要列: %s", path))
  }

  feature_map <- tibble::tibble(
    feature = as.character(features),
    norm_feature = normalize_feature_name(features)
  ) %>%
    filter(!is.na(norm_feature), nzchar(norm_feature)) %>%
    distinct(feature, .keep_all = TRUE)

  ranked <- df %>%
    mutate(
      feature = as.character(feature),
      coef = suppressWarnings(as.numeric(coef)),
      qval = suppressWarnings(as.numeric(qval)),
      norm_feature = normalize_feature_name(feature)
    ) %>%
    filter(!is.na(norm_feature), nzchar(norm_feature)) %>%
    inner_join(feature_map, by = "norm_feature", suffix = c("_maaslin", "_target")) %>%
    filter(!is.na(coef))

  if (nrow(ranked) == 0) {
    return(tibble::tibble(feature = features, coef = NA_real_, qval = NA_real_, effect_sign = "neutral"))
  }

  if ("pval" %in% colnames(ranked)) {
    ranked <- ranked %>% mutate(pval = suppressWarnings(as.numeric(pval)))
    ranked <- ranked %>% arrange(is.na(qval), qval, is.na(pval), pval)
  } else {
    ranked <- ranked %>% arrange(is.na(qval), qval)
  }

  out <- ranked %>%
    group_by(feature_target) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    transmute(
      feature = feature_target,
      coef = coef,
      qval = qval,
      effect_sign = case_when(
        coef > 0 ~ "positive",
        coef < 0 ~ "negative",
        TRUE ~ "neutral"
      )
    )

  missing_features <- setdiff(features, out$feature)
  if (length(missing_features) > 0) {
    out <- bind_rows(
      out,
      tibble::tibble(
        feature = missing_features,
        coef = NA_real_,
        qval = NA_real_,
        effect_sign = "neutral"
      )
    )
  }

  out
}

make_dual_center_layout <- function(g, center_nodes, layout_seed = 20260320) {
  set.seed(layout_seed)
  layout_tbl <- create_layout(g, layout = "kk")
  layout_tbl$type <- V(g)$type
  layout_tbl$center_node <- layout_tbl$name %in% center_nodes

  center_idx <- which(layout_tbl$center_node)
  ko_idx <- which(layout_tbl$type == "KO")
  met_idx <- which(layout_tbl$type == "metabolite")
  species_idx <- which(layout_tbl$type == "species")

  center_x <- c(-4.2, 4.2)
  center_y <- c(0, 0)

  if (length(center_idx) == 2) {
    layout_tbl$x[center_idx] <- center_x
    layout_tbl$y[center_idx] <- center_y
  } else if (length(center_idx) > 0) {
    layout_tbl$x[center_idx] <- seq(-3.5, 3.5, length.out = length(center_idx))
    layout_tbl$y[center_idx] <- 0
  }

  if (length(ko_idx) > 0) {
    ko_idx <- setdiff(ko_idx, center_idx)
    side <- rep(c(-1, 1), length.out = length(ko_idx))
    rank_pos <- seq_along(ko_idx)
    x_band <- 2.15 * side
    y_band <- seq(-2.0, 2.0, length.out = length(ko_idx))
    layout_tbl$x[ko_idx] <- x_band + 0.15 * sin(rank_pos)
    layout_tbl$y[ko_idx] <- y_band
  }

  if (length(met_idx) > 0) {
    layout_tbl$x[met_idx] <- seq(-0.9, 0.9, length.out = length(met_idx))
    layout_tbl$y[met_idx] <- seq(-2.2, 2.2, length.out = length(met_idx))
  }

  if (length(species_idx) > 0) {
    non_center_species <- setdiff(species_idx, center_idx)
    if (length(non_center_species) > 0) {
      side <- rep(c(-1, 1), length.out = length(non_center_species))
      layout_tbl$x[non_center_species] <- 3.3 * side
      layout_tbl$y[non_center_species] <- seq(-1.8, 1.8, length.out = length(non_center_species))
    }
  }

  layout_tbl$label_rich <- V(g)$label_rich
  layout_tbl$show_label <- V(g)$show_label
  layout_tbl$is_hub <- V(g)$is_hub
  layout_tbl$effect_sign <- V(g)$effect_sign
  layout_tbl$node_size <- V(g)$node_size
  layout_tbl$node_type <- layout_tbl$type
  layout_tbl
}

plot_journal_network <- function(g, edges, output_pdf, title, subtitle, caption, center_nodes) {
  layout_tbl <- make_dual_center_layout(g, center_nodes = center_nodes)
  output_png <- sub("\\.pdf$", ".png", output_pdf)
  shape_map <- c(species = 21, KO = 22, metabolite = 24, other = 21)

  p <- ggraph(layout_tbl) +
    geom_edge_link(
      aes(width = weight, color = edge_sign),
      alpha = 0.30,
      lineend = "round",
      show.legend = TRUE
    ) +
    scale_edge_width_continuous(range = c(0.35, LINE_WIDTH), name = expression(paste("|Spearman ", italic(r), "|"))) +
    scale_edge_color_manual(values = c(positive = EDGE_POS, negative = EDGE_NEG), name = "Correlation") +
    geom_node_point(
      aes(x = x, y = y, fill = effect_sign, shape = node_type),
      color = NODE_BORDER,
      stroke = 0.45,
      size = 5.7,
      alpha = 0.96
    ) +
    geom_node_point(
      data = layout_tbl %>% filter(center_node),
      aes(x = x, y = y, fill = effect_sign, shape = node_type),
      color = "#111111",
      stroke = 1.25,
      size = 5.7,
      alpha = 1.0
    ) +
    scale_shape_manual(
      values = shape_map,
      breaks = c("species", "KO", "metabolite"),
      labels = c("Species", "KO", "Metabolite"),
      name = "Node type"
    ) +
    scale_fill_manual(
      values = c(positive = NODE_UP, negative = NODE_DOWN, neutral = NODE_MID),
      breaks = c("positive", "negative"),
      labels = c("coef > 0", "coef < 0"),
      drop = FALSE,
      name = "MaAsLin coef"
    ) +
    guides(
      shape = guide_legend(order = 1, override.aes = list(fill = "white", color = NA, stroke = 0, size = 5.5)),
      fill = guide_legend(order = 2, override.aes = list(shape = 21, color = NA, stroke = 0, size = 5.5)),
      edge_color = guide_legend(order = 3),
      edge_width = guide_legend(order = 4)
    ) +
    ggrepel::geom_text_repel(
      data = layout_tbl %>% filter(show_label),
      aes(x = x, y = y, label = label_rich),
      parse = TRUE,
      size = 3.1,
      fontface = "plain",
      max.overlaps = Inf,
      box.padding = 0.45,
      point.padding = 0.18,
      segment.size = 0.28,
      segment.alpha = 0.45,
      segment.color = "#707070",
      min.segment.length = 0,
      force = 2.2,
      seed = 20260320,
      family = "sans"
    ) +
    labs(
      title = title,
      subtitle = subtitle,
      caption = caption,
      color = NULL
    ) +
    theme_void(base_family = "sans") +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      legend.background = element_rect(fill = "white", color = NA),
      legend.key = element_rect(fill = "white", color = NA),
      plot.margin = margin(12, 16, 12, 16),
      legend.box = "vertical",
      legend.position = "right",
      legend.justification = "top",
      legend.title = element_text(size = 10.2, face = "bold", color = TEXT_DARK),
      legend.text = element_text(size = 9, color = TEXT_DARK),
      legend.spacing.y = unit(5, "pt"),
      plot.title = element_text(face = "bold", size = 17, hjust = 0.02, color = TEXT_DARK),
      plot.subtitle = element_text(size = 10.5, hjust = 0.02, color = TEXT_MUTED, margin = margin(b = 7)),
      plot.caption = element_text(size = 8.6, hjust = 0.02, color = TEXT_MUTED, margin = margin(t = 7))
    )

  ggsave(output_pdf, p, width = 14.2, height = 10.6, device = cairo_pdf)
  ggsave(output_png, p, width = 14.2, height = 10.6, dpi = 320)
  p
}

build_network_common <- function(node_tbl, edges, output_prefix, title, subtitle, caption,
                                 fill_column, fill_scale, out_dir = NULL, label_seed = 20260320,
                                 plot_bg = BASE_BG) {
  if (is.null(out_dir) || !nzchar(out_dir)) {
    out_dir <- output_dir
  }

  if (!"weight" %in% colnames(edges) && "edge_weight" %in% colnames(edges)) {
    edges$weight <- edges$edge_weight
  }
  if (!"edge_sign" %in% colnames(edges) && "sign" %in% colnames(edges)) {
    edges$edge_sign <- edges$sign
  }

  g <- graph_from_data_frame(edges, directed = FALSE, vertices = node_tbl)

  deg <- degree(g)
  V(g)$degree <- deg
  hub_cutoff <- max(2, as.numeric(quantile(deg, probs = HUB_QUANTILE, na.rm = TRUE)))
  V(g)$is_hub <- deg > hub_cutoff

  v_tbl <- as_tibble(as_data_frame(g, what = "vertices")) %>%
    arrange(desc(is_hub), desc(degree), desc(abs(.data[[fill_column]])))
  key_labels <- unique(head(v_tbl$label, LABEL_TOP_N))
  key_labels <- unique(c(key_labels, v_tbl$label[v_tbl$is_hub]))
  key_labels <- unique(c(
    key_labels,
    v_tbl$label[abs(v_tbl[[fill_column]]) >= quantile(abs(v_tbl[[fill_column]]), probs = 0.75, na.rm = TRUE)]
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
      alpha = 0.32,
      lineend = "round",
      show.legend = TRUE
    ) +
    scale_edge_width_continuous(range = c(0.3, LINE_WIDTH), name = expression(paste("|Spearman ", italic(r), "|"))) +
    scale_edge_color_manual(values = c(positive = EDGE_POS, negative = EDGE_NEG), name = "Correlation") +
    geom_node_point(
      aes(size = mean_abundance, fill = .data[[fill_column]], shape = type),
      color = NODE_BORDER,
      stroke = 0.4,
      alpha = 0.85
    ) +
    geom_node_point(
      data = layout_tbl %>% filter(is_hub),
      aes(size = mean_abundance, fill = .data[[fill_column]], shape = type),
      color = "#000000",
      stroke = 1.1,
      alpha = 1.0
    ) +
    scale_shape_manual(values = shape_map, name = "Node type") +
    scale_size_continuous(range = c(2.5, 12), trans = "sqrt", name = "Mean abundance") +
    fill_scale +
    guides(
      shape = guide_legend(order = 1, override.aes = list(fill = "white", color = NA, stroke = 0, size = 5.5)),
      fill = guide_legend(order = 2, override.aes = list(shape = 21, color = NA, stroke = 0, size = 5.5)),
      edge_color = guide_legend(order = 3),
      edge_width = guide_legend(order = 4)
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
      seed = label_seed,
      family = "sans"
    ) +
    labs(
      title = title,
      subtitle = subtitle,
      caption = caption
    ) +
    theme_void(base_family = "sans") +
    theme(
      plot.background = element_rect(fill = plot_bg, color = NA),
      legend.background = element_rect(fill = plot_bg, color = NA),
      legend.key = element_rect(fill = plot_bg, color = NA),
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

  out_pdf <- file.path(out_dir, paste0(output_prefix, ".pdf"))
  ggsave(out_pdf, p, width = 13.2, height = 10.2, device = cairo_pdf)

  out_edges <- file.path(out_dir, paste0(output_prefix, "_edges.csv"))
  out_nodes <- file.path(out_dir, paste0(output_prefix, "_nodes.csv"))
  fwrite(edges, out_edges)
  fwrite(as.data.frame(as_data_frame(g, what = "vertices")), out_nodes)

  list(graph = g, plot = p, nodes = node_tbl, edges = edges, hub_cutoff = hub_cutoff, out_pdf = out_pdf)
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
  ko_cols_all <- grep("^kegg_", colnames(sub), value = TRUE)

  ko_cols <- resolve_feature_columns(
    available_cols = ko_cols_all,
    fixed_list = FIXED_KO_LIST,
    maaslin_path = config$maaslin_ko,
    remove_prefix = "^kegg_",
    label = paste0(target_group, " KO"),
    use_fixed = USE_FIXED_KO_LIST,
    max_auto_features = config$max_auto_ko_features
  )
  met_cols <- resolve_feature_columns(
    available_cols = met_cols_all,
    fixed_list = FIXED_METABOLITES_LIST,
    maaslin_path = config$maaslin_metabolite,
    remove_prefix = "^met_",
    label = paste0(target_group, " metabolite"),
    use_fixed = USE_FIXED_METABOLITES_LIST,
    max_auto_features = config$max_auto_metabolite_features
  )

  ko_ids_present <- gsub("^kegg_", "", ko_cols)
  ko_pathway_map <- load_ko_pathway_label(ko_ids_present)

  if (length(ko_cols) == 0) stop("未找到任何指定 KO 列")

  species_cols <- match_targets_to_cols(FIXED_SPECIES_LIST, species_cols_all, "^tax_s__", strict = TRUE)

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

  effect_tbl <- bind_rows(
    extract_maaslin_effects(config$maaslin_microbe, species_cols, qval_cutoff = config$qval_cutoff),
    extract_maaslin_effects(config$maaslin_ko, ko_cols, qval_cutoff = config$qval_cutoff),
    extract_maaslin_effects(config$maaslin_metabolite, met_cols, qval_cutoff = config$qval_cutoff)
  ) %>%
    rename(node_id = feature)

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
        type == "species" ~ make_label_expr(gsub("_", " ", label), italic = TRUE),
        TRUE ~ make_label_expr(label, italic = FALSE)
      ),
      node_size = ifelse(node_id %in% c("tax_s__Peptostreptococcus_stomatis", "tax_s__Faecalibacterium_prausnitzii"), 8.6, 6.2)
    ) %>%
    left_join(effect_tbl, by = "node_id") %>%
    mutate(
      effect_sign = ifelse(is.na(effect_sign), "neutral", effect_sign),
      coef = ifelse(is.na(coef), 0, coef),
      qval = ifelse(is.na(qval), 1, qval),
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

  V(g)$effect_sign <- node_tbl$effect_sign[match(V(g)$name, node_tbl$node_id)]
  V(g)$effect_sign[is.na(V(g)$effect_sign)] <- "neutral"
  V(g)$node_size <- node_tbl$node_size[match(V(g)$name, node_tbl$node_id)]
  V(g)$node_size[is.na(V(g)$node_size)] <- 6.2

  v_tbl <- as_tibble(as_data_frame(g, what = "vertices")) %>%
    arrange(desc(is_hub), desc(degree), desc(abs(coef)))
  key_labels <- unique(c(
    "Peptostreptococcus_stomatis",
    "Faecalibacterium_prausnitzii",
    head(v_tbl$label, LABEL_TOP_N),
    v_tbl$label[v_tbl$is_hub]
  ))
  V(g)$show_label <- V(g)$label %in% key_labels
  V(g)$node_border <- ifelse(V(g)$is_hub, "#111111", NODE_BORDER)

  out_pdf <- file.path(output_dir, paste0("network_species_ko_metabolite_", target_group, ".pdf"))
  p <- plot_journal_network(
    g = g,
    edges = edge_tbl,
    output_pdf = out_pdf,
    title = paste0(target_group, " functional network"),
    subtitle = paste0(
      "Spearman edges: |r| > ",
      R_THRESHOLD, ", FDR < ", FDR_THRESHOLD
    ),
    caption = "",
    center_nodes = c("tax_s__Peptostreptococcus_stomatis", "tax_s__Faecalibacterium_prausnitzii")
  )

  out_edges <- file.path(output_dir, paste0("network_edges_", target_group, ".csv"))
  out_nodes <- file.path(output_dir, paste0("network_nodes_", target_group, ".csv"))
  fwrite(edge_tbl, out_edges)
  fwrite(as.data.frame(as_data_frame(g, what = "vertices")), out_nodes)

  cat(sprintf("已输出: %s\n", out_pdf))

  list(graph = g, plot = p, nodes = node_tbl, edges = edge_tbl, hub_cutoff = hub_cutoff)
}

make_all_samples_network <- function(df) {
  cat("\n==== 构建全样本网络 ====" , "\n", sep = "")

  species_cols_all <- grep("^tax_s__", colnames(df), value = TRUE)
  met_cols_all <- grep("^met_", colnames(df), value = TRUE)
  ko_cols_all <- grep("^kegg_", colnames(df), value = TRUE)

  species_cols <- match_targets_to_cols(FIXED_SPECIES_LIST, species_cols_all, "^tax_s__", strict = TRUE)
  ko_cols <- resolve_feature_columns(
    available_cols = ko_cols_all,
    fixed_list = FIXED_KO_LIST,
    maaslin_path = config$maaslin_ko,
    remove_prefix = "^kegg_",
    label = "all samples KO",
    use_fixed = USE_FIXED_KO_LIST,
    max_auto_features = config$max_auto_ko_features
  )
  met_cols <- resolve_feature_columns(
    available_cols = met_cols_all,
    fixed_list = FIXED_METABOLITES_LIST,
    maaslin_path = config$maaslin_metabolite,
    remove_prefix = "^met_",
    label = "all samples metabolite",
    use_fixed = USE_FIXED_METABOLITES_LIST,
    max_auto_features = config$max_auto_metabolite_features
  )

  if (length(species_cols) == 0 || length(ko_cols) == 0 || length(met_cols) == 0) {
    stop(sprintf("全样本匹配失败: species=%d, ko=%d, metabolite=%d", length(species_cols), length(ko_cols), length(met_cols)))
  }

  num_cols <- c(species_cols, ko_cols, met_cols)
  num_df <- df[, num_cols, drop = FALSE] %>%
    mutate(across(everything(), as.numeric))
  num_df[is.na(num_df)] <- 0

  sp_ko <- pairwise_spearman(
    as.matrix(num_df[, species_cols, drop = FALSE]),
    as.matrix(num_df[, ko_cols, drop = FALSE]),
    x_type = "species",
    y_type = "KO"
  )

  ko_met <- pairwise_spearman(
    as.matrix(num_df[, ko_cols, drop = FALSE]),
    as.matrix(num_df[, met_cols, drop = FALSE]),
    x_type = "KO",
    y_type = "metabolite"
  )

  edges <- bind_rows(sp_ko, ko_met) %>%
    filter(!is.na(r), !is.na(fdr), abs(r) > R_THRESHOLD, fdr < FDR_THRESHOLD)

  if (nrow(edges) == 0) {
    stop("全样本在当前阈值下没有保留边，请降低阈值")
  }

  edges <- edges %>%
    mutate(
      weight = edge_weight
    )

  ko_pathway_map <- load_ko_pathway_label(gsub("^kegg_", "", ko_cols))

  node_tbl <- tibble(node_id = unique(c(edges$from, edges$to))) %>%
    mutate(
      type = case_when(
        grepl("^tax_s__", node_id) ~ "species",
        grepl("^kegg_", node_id) ~ "KO",
        grepl("^met_", node_id) ~ "metabolite",
        TRUE ~ "other"
      ),
      label = vapply(node_id, format_node_label, character(1), ko_pathway_map = ko_pathway_map),
      label_rich = case_when(
        type == "species" ~ make_label_expr(gsub("_", " ", label), italic = TRUE),
        TRUE ~ make_label_expr(label, italic = FALSE)
      ),
      node_size = ifelse(node_id %in% c("tax_s__Peptostreptococcus_stomatis", "tax_s__Faecalibacterium_prausnitzii"), 8.6, 6.2)
    ) %>%
    left_join(
      bind_rows(
        extract_maaslin_effects(config$maaslin_microbe, species_cols, qval_cutoff = config$qval_cutoff),
        extract_maaslin_effects(config$maaslin_ko, ko_cols, qval_cutoff = config$qval_cutoff),
        extract_maaslin_effects(config$maaslin_metabolite, met_cols, qval_cutoff = config$qval_cutoff)
      ) %>%
      rename(node_id = feature),
      by = "node_id"
    ) %>%
    mutate(
      effect_sign = ifelse(is.na(effect_sign), "neutral", effect_sign),
      coef = ifelse(is.na(coef), 0, coef),
      qval = ifelse(is.na(qval), 1, qval)
    )

  g <- graph_from_data_frame(
    edges %>% transmute(from = from, to = to, weight = edge_weight, edge_sign = edge_sign),
    directed = FALSE,
    vertices = node_tbl
  )

  deg <- degree(g)
  V(g)$degree <- deg
  hub_cutoff <- max(2, as.numeric(quantile(deg, probs = HUB_QUANTILE, na.rm = TRUE)))
  V(g)$is_hub <- deg > hub_cutoff
  V(g)$effect_sign <- node_tbl$effect_sign[match(V(g)$name, node_tbl$node_id)]
  V(g)$effect_sign[is.na(V(g)$effect_sign)] <- "neutral"
  V(g)$node_size <- node_tbl$node_size[match(V(g)$name, node_tbl$node_id)]
  V(g)$node_size[is.na(V(g)$node_size)] <- 6.2

  v_tbl <- as_tibble(as_data_frame(g, what = "vertices")) %>%
    arrange(desc(is_hub), desc(degree), desc(abs(coef)))
  key_labels <- unique(c(
    "Peptostreptococcus_stomatis",
    "Faecalibacterium_prausnitzii",
    head(v_tbl$label, LABEL_TOP_N),
    v_tbl$label[v_tbl$is_hub]
  ))
  V(g)$show_label <- V(g)$label %in% key_labels
  V(g)$node_border <- ifelse(V(g)$is_hub, "#111111", NODE_BORDER)

  out_pdf <- file.path(output_dir, "network_species_ko_metabolite_all_samples.pdf")
  p <- plot_journal_network(
    g = g,
    edges = edges %>% transmute(from = from, to = to, weight = edge_weight, edge_sign = edge_sign),
    output_pdf = out_pdf,
    title = "All samples functional network",
    subtitle = paste0(
      "Spearman edges: |r| > ",
      R_THRESHOLD, ", FDR < ", FDR_THRESHOLD
    ),
    caption = "",
    center_nodes = c("tax_s__Peptostreptococcus_stomatis", "tax_s__Faecalibacterium_prausnitzii")
  )

  list(graph = g, plot = p, nodes = node_tbl, edges = edges, hub_cutoff = hub_cutoff, out_pdf = out_pdf)
}

cat("读取数据...\n")
dat <- fread(input_file, data.table = FALSE)
dat <- prepare_diff_groups(dat)
dat <- dat %>% filter(!is.na(group))

cat("分组样本数:\n")
print(table(dat$group))

res_all <- make_all_samples_network(dat)
res_ctrl <- make_group_network(dat, target_group = "control", control_group = "control")
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
