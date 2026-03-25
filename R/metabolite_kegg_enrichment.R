suppressPackageStartupMessages({
  library(KEGGREST)
  library(ggplot2)
})

read_tsv_flexible <- function(path) {
  if (requireNamespace("data.table", quietly = TRUE)) {
    return(data.table::fread(path, sep = "\t", data.table = FALSE, check.names = FALSE))
  }
  read.delim(
    path,
    sep = "\t",
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    quote = "",
    comment.char = "",
    fill = TRUE
  )
}

extract_kegg_ids <- function(x) {
  if (is.na(x) || x == "") {
    return(character(0))
  }
  ids <- regmatches(x, gregexpr("C[0-9]{5}", x, perl = TRUE))[[1]]
  unique(ids)
}

kegg_link_compound_pathway <- function(compounds, chunk_size = 100) {
  out <- list()
  if (length(compounds) == 0) {
    return(data.frame(kegg_id = character(0), pathway_id = character(0), stringsAsFactors = FALSE))
  }

  idx <- seq(1, length(compounds), by = chunk_size)
  for (i in idx) {
    batch <- compounds[i:min(i + chunk_size - 1, length(compounds))]
    query <- paste0("cpd:", batch)
    links <- tryCatch(KEGGREST::keggLink("pathway", query), error = function(e) NULL)
    if (is.null(links) || length(links) == 0) {
      next
    }

    out[[length(out) + 1]] <- data.frame(
      kegg_id = sub("^cpd:", "", names(links)),
      pathway_id = sub("^path:", "", unname(links)),
      stringsAsFactors = FALSE
    )
  }

  if (length(out) == 0) {
    return(data.frame(kegg_id = character(0), pathway_id = character(0), stringsAsFactors = FALSE))
  }

  unique(do.call(rbind, out))
}

args <- commandArgs(trailingOnly = TRUE)
q_cutoff <- ifelse(length(args) >= 1, as.numeric(args[1]), 0.05)

sig_file <- "results/maaslin2_metabolite/significant_results.tsv"
combine_file <- "dataset/rawdata/combine.intensity.xls"
out_dir <- "results/pathway_metabolite"
KEGG_CACHE_DIR <- "dataset/KEGG/cache"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(KEGG_CACHE_DIR, recursive = TRUE, showWarnings = FALSE)

get_cached_compound_pathway_links <- function(compounds, cache_dir = KEGG_CACHE_DIR) {
  link_cache_file <- file.path(cache_dir, "kegg_compound_pathway_links.csv")

  link_df_all <- NULL
  if (file.exists(link_cache_file)) {
    link_df_all <- read.csv(link_cache_file, stringsAsFactors = FALSE)
    if (!all(c("kegg_id", "pathway_id") %in% colnames(link_df_all))) {
      link_df_all <- NULL
    }
  }

  if (is.null(link_df_all)) {
    link_df_all <- data.frame(kegg_id = character(0), pathway_id = character(0), stringsAsFactors = FALSE)
  }

  known_compounds <- unique(link_df_all$kegg_id)
  miss_compounds <- setdiff(compounds, known_compounds)

  if (length(miss_compounds) > 0) {
    cat(sprintf("  KEGG cache 缺失 compound 数: %d，正在增量拉取...\n", length(miss_compounds)))
    miss_links <- kegg_link_compound_pathway(miss_compounds)
    if (nrow(miss_links) > 0) {
      link_df_all <- unique(rbind(link_df_all, miss_links))
      write.csv(link_df_all, link_cache_file, row.names = FALSE)
    }
  }

  unique(link_df_all[link_df_all$kegg_id %in% compounds, c("kegg_id", "pathway_id"), drop = FALSE])
}

get_cached_pathway_names <- function(pathway_ids, cache_dir = KEGG_CACHE_DIR) {
  pathway_cache_file <- file.path(cache_dir, "kegg_pathway_metadata.csv")

  pathway_meta <- NULL
  if (file.exists(pathway_cache_file)) {
    pathway_meta <- read.csv(pathway_cache_file, stringsAsFactors = FALSE)
    if (!all(c("pathway_id", "pathway_name") %in% colnames(pathway_meta))) {
      pathway_meta <- NULL
    }
  }

  if (is.null(pathway_meta)) {
    pathway_name_vec <- tryCatch(KEGGREST::keggList("pathway"), error = function(e) NULL)
    if (!is.null(pathway_name_vec)) {
      pathway_meta <- data.frame(
        pathway_id = sub("^path:", "", names(pathway_name_vec)),
        pathway_name = unname(pathway_name_vec),
        stringsAsFactors = FALSE
      )
      write.csv(pathway_meta, pathway_cache_file, row.names = FALSE)
    } else {
      pathway_meta <- data.frame(pathway_id = character(0), pathway_name = character(0), stringsAsFactors = FALSE)
    }
  }

  miss_pathways <- setdiff(pathway_ids, pathway_meta$pathway_id)
  if (length(miss_pathways) > 0) {
    pathway_meta <- rbind(
      pathway_meta,
      data.frame(pathway_id = miss_pathways, pathway_name = miss_pathways, stringsAsFactors = FALSE)
    )
  }

  unique(pathway_meta[pathway_meta$pathway_id %in% pathway_ids, c("pathway_id", "pathway_name"), drop = FALSE])
}

cat("[1/7] 读取差异结果...\n")
sig <- read_tsv_flexible(sig_file)
if (!all(c("feature", "qval") %in% colnames(sig))) {
  stop("差异结果缺少必要列: feature/qval")
}

sig$qval <- suppressWarnings(as.numeric(sig$qval))
sig$coef <- suppressWarnings(as.numeric(sig$coef))
sig_filt <- sig[!is.na(sig$qval) & sig$qval < q_cutoff, , drop = FALSE]
write.table(
  sig_filt,
  file = file.path(out_dir, "metabolite_qval_lt_0.05.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

sig_features <- unique(as.character(sig_filt$feature))
cat(sprintf("  qval < %.3f 的特征数: %d\n", q_cutoff, length(sig_features)))

cat("[2/7] 读取并构建 combine 注释映射...\n")
comb <- read_tsv_flexible(combine_file)
for (nm in c("ID", "MS2Metabolite", "MS2kegg")) {
  if (!nm %in% colnames(comb)) {
    comb[[nm]] <- ""
  }
}

ms2 <- trimws(as.character(comb$MS2Metabolite))
ms2kegg <- trimws(as.character(comb$MS2kegg))

keep <- ms2 != "" & ms2kegg != ""
map_tbl <- data.frame(
  feature = make.names(paste0("met_", ms2), unique = FALSE),
  annotation = ms2,
  ms2kegg = ms2kegg,
  keep = keep,
  stringsAsFactors = FALSE
)
map_tbl <- map_tbl[map_tbl$keep, c("feature", "annotation", "ms2kegg"), drop = FALSE]

kegg_list <- lapply(map_tbl$ms2kegg, extract_kegg_ids)
map_tbl$kegg_ids <- vapply(kegg_list, function(x) paste(x, collapse = ";"), character(1))

expand_idx <- lengths(kegg_list)
map_expand <- map_tbl[rep(seq_len(nrow(map_tbl)), expand_idx), c("feature", "annotation"), drop = FALSE]
if (nrow(map_expand) > 0) {
  map_expand$kegg_id <- unlist(kegg_list)
  map_expand <- unique(map_expand)
} else {
  map_expand <- data.frame(feature = character(0), annotation = character(0), kegg_id = character(0), stringsAsFactors = FALSE)
}

cat(sprintf("  可映射到 KEGG 的背景代谢物注释条目: %d\n", nrow(map_expand)))

cat("[3/7] 匹配显著特征与 KEGG...\n")
sig_key <- data.frame(feature = sig_features, stringsAsFactors = FALSE)
sig_map <- merge(sig_key, map_expand, by = "feature", all.x = TRUE)
sig_coef <- unique(sig_filt[, c("feature", "coef"), drop = FALSE])
sig_map <- merge(sig_map, sig_coef, by = "feature", all.x = TRUE)

metabolite_kegg_map <- sig_map[, c("feature", "annotation", "kegg_id", "coef"), drop = FALSE]
colnames(metabolite_kegg_map) <- c("feature", "metabolite_name", "kegg_id", "coef")
metabolite_kegg_map$qval <- sig_filt$qval[match(metabolite_kegg_map$feature, sig_filt$feature)]
metabolite_kegg_map <- metabolite_kegg_map[order(metabolite_kegg_map$metabolite_name, metabolite_kegg_map$kegg_id), , drop = FALSE]

mapped_feature_n <- length(unique(sig_map$feature[!is.na(sig_map$kegg_id)]))
unmapped_features <- setdiff(sig_features, unique(sig_map$feature[!is.na(sig_map$kegg_id)]))

write.table(
  sig_map,
  file = file.path(out_dir, "metabolite_feature_kegg_mapping.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  metabolite_kegg_map,
  file = file.path(out_dir, "metabolite_name_keggid_map.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  data.frame(feature = unmapped_features, stringsAsFactors = FALSE),
  file = file.path(out_dir, "metabolite_unmapped_features.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat(sprintf("  显著特征中成功映射 KEGG 的特征数: %d/%d\n", mapped_feature_n, length(sig_features)))

cat("[4/7] KEGGREST 获取化合物-通路关联...\n")
bg_compounds <- unique(map_expand$kegg_id)
sig_compounds <- unique(sig_map$kegg_id[!is.na(sig_map$kegg_id)])

cpd2path <- get_cached_compound_pathway_links(bg_compounds, cache_dir = KEGG_CACHE_DIR)
write.table(
  cpd2path,
  file = file.path(out_dir, "kegg_compound_pathway_links.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

if (nrow(cpd2path) == 0) {
  stop("未获取到 KEGG compound-pathway 关联，无法进行富集分析。")
}

cat("[5/7] 计算通路富集（超几何检验）...\n")
bg_mapped <- unique(cpd2path$kegg_id)
sig_in_bg <- intersect(sig_compounds, bg_mapped)
M <- length(bg_mapped)
n <- length(sig_in_bg)

if (M == 0 || n == 0) {
  stop("可用于富集计算的化合物为空（M 或 n 为 0）。")
}

pathways <- sort(unique(cpd2path$pathway_id))
enrich_rows <- lapply(pathways, function(pid) {
  path_cpd <- unique(cpd2path$kegg_id[cpd2path$pathway_id == pid])
  K <- length(path_cpd)
  k <- length(intersect(sig_in_bg, path_cpd))
  if (k == 0) {
    return(NULL)
  }
  pval <- phyper(k - 1, K, M - K, n, lower.tail = FALSE)
  # pathway impact 使用富集命中占通路化合物比例作为拓扑影响近似值
  impact <- k / K

  path_hits <- intersect(sig_in_bg, path_cpd)
  hit_rows <- sig_map[!is.na(sig_map$coef) & sig_map$kegg_id %in% path_hits, c("feature", "kegg_id", "coef"), drop = FALSE]
  dir_score <- if (nrow(hit_rows) == 0) 0 else mean(hit_rows$coef, na.rm = TRUE)
  dir_label <- ifelse(dir_score >= 0, "Up", "Down")

  data.frame(
    pathway_id = pid,
    hit_count = k,
    pathway_size = K,
    pathway_impact = impact,
    sig_in_bg = n,
    bg_size = M,
    pval = pval,
    direction_score = dir_score,
    direction_label = dir_label,
    hit_compounds = paste(sort(intersect(sig_in_bg, path_cpd)), collapse = ";"),
    stringsAsFactors = FALSE
  )
})

enrich <- do.call(rbind, enrich_rows)
if (is.null(enrich) || nrow(enrich) == 0) {
  enrich <- data.frame(
    pathway_id = character(0),
    hit_count = integer(0),
    pathway_size = integer(0),
    sig_in_bg = integer(0),
    bg_size = integer(0),
    pval = numeric(0),
    hit_compounds = character(0),
    stringsAsFactors = FALSE
  )
} else {
  enrich$qval <- p.adjust(enrich$pval, method = "BH")
  enrich <- enrich[order(enrich$pval, -enrich$hit_count), , drop = FALSE]
}

if (nrow(enrich) > 0) {
  path_meta <- get_cached_pathway_names(unique(enrich$pathway_id), cache_dir = KEGG_CACHE_DIR)
  enrich <- merge(enrich, path_meta, by = "pathway_id", all.x = TRUE)
  enrich$pathway_name[is.na(enrich$pathway_name) | enrich$pathway_name == ""] <- enrich$pathway_id
}

if (nrow(enrich) > 0) {
  enrich <- enrich[order(enrich$pval, -enrich$hit_count), , drop = FALSE]
  enrich <- enrich[, c("pathway_id", "pathway_name", "hit_count", "pathway_size", "pathway_impact", "sig_in_bg", "bg_size", "pval", "qval", "direction_score", "direction_label", "hit_compounds")]
}

write.table(
  enrich,
  file = file.path(out_dir, "kegg_pathway_enrichment.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("[6/7] 输出汇总...\n")
summary_lines <- c(
  sprintf("qval cutoff: %.4f", q_cutoff),
  sprintf("Significant features (qval < cutoff): %d", length(sig_features)),
  sprintf("Significant features mapped to KEGG: %d", mapped_feature_n),
  sprintf("Unmapped significant features: %d", length(unmapped_features)),
  sprintf("Background compounds with pathway mapping (M): %d", M),
  sprintf("Significant compounds in background (n): %d", n),
  sprintf("Enriched pathways tested: %d", nrow(enrich)),
  sprintf("Pathways with qval < 0.05: %d", ifelse(nrow(enrich) == 0, 0, sum(enrich$qval < 0.05, na.rm = TRUE)))
)
writeLines(summary_lines, con = file.path(out_dir, "kegg_enrichment_summary.txt"))

cat("[7/7] 绘制可视化图（气泡图与富集圈图）...\n")

if (nrow(enrich) > 0) {
  enrich_plot <- enrich
  enrich_plot$neg_log10_q <- -log10(pmax(enrich_plot$qval, .Machine$double.xmin))
  enrich_plot$gene_ratio <- enrich_plot$hit_count / pmax(enrich_plot$sig_in_bg, 1)

  top_n <- min(50, nrow(enrich_plot))
  top_idx <- order(enrich_plot$qval, -enrich_plot$hit_count)
  top_plot <- enrich_plot[top_idx[seq_len(top_n)], , drop = FALSE]
  top_plot <- top_plot[order(top_plot$neg_log10_q), , drop = FALSE]

  top_plot$pathway_label <- top_plot$pathway_name
  top_plot$pathway_label <- ifelse(
    nchar(top_plot$pathway_label) > 60,
    paste0(substr(top_plot$pathway_label, 1, 57), "..."),
    top_plot$pathway_label
  )
  top_plot$pathway_label <- factor(top_plot$pathway_label, levels = top_plot$pathway_label)

  direction_colors <- c("Up" = "#C23B22", "Down" = "#2F5D8A")

  bubble_plot <- ggplot(
    top_plot,
    aes(
      x = gene_ratio,
      y = pathway_label,
      size = hit_count,
      color = direction_label
    )
  ) +
    geom_point(alpha = 0.88) +
    scale_color_manual(values = direction_colors, drop = FALSE) +
    scale_size_continuous(range = c(2.8, 10), breaks = pretty(top_plot$hit_count, n = 4)) +
    labs(
      title = "KEGG Pathway Enrichment Bubble Plot",
      subtitle = paste0("Top ", top_n, " pathways ranked by adjusted P-value"),
      x = "Gene Ratio (Hit / Significant Compounds)",
      y = NULL,
      size = "Hit Count",
      color = "Direction"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14, color = "#222222"),
      plot.subtitle = element_text(size = 10, color = "#555555"),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.y = element_text(size = 9, color = "#222222"),
      axis.text.x = element_text(color = "#222222"),
      legend.position = "right",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9)
    )

  ggsave(
    filename = file.path(out_dir, "kegg_enrichment_bubble_plot.pdf"),
    plot = bubble_plot,
    width = 10,
    height = 7.5,
    dpi = 320,
    bg = "white"
  )
  ggsave(
    filename = file.path(out_dir, "kegg_enrichment_bubble_plot.png"),
    plot = bubble_plot,
    width = 10,
    height = 7.5,
    dpi = 320,
    bg = "white"
  )

  circle_plot_df <- top_plot[order(top_plot$neg_log10_q, decreasing = TRUE), , drop = FALSE]
  circle_plot_df$pathway_label <- factor(circle_plot_df$pathway_label, levels = circle_plot_df$pathway_label)

  circle_plot <- ggplot(
    circle_plot_df,
    aes(x = pathway_label, y = neg_log10_q, fill = direction_label)
  ) +
    geom_col(width = 0.92, alpha = 0.94, color = "white", linewidth = 0.2) +
    geom_point(
      aes(y = neg_log10_q + 0.2, size = hit_count),
      shape = 21,
      stroke = 0.3,
      color = "#333333",
      fill = "white"
    ) +
    scale_fill_manual(values = direction_colors, drop = FALSE) +
    scale_size_continuous(range = c(2.2, 7.8), breaks = pretty(circle_plot_df$hit_count, n = 4)) +
    coord_polar(start = 0) +
    labs(
      title = "KEGG Enrichment Circle Plot",
      subtitle = "Bar height: -log10(q-value); Point size: hit count",
      x = NULL,
      y = "-log10(q-value)",
      fill = "Direction",
      size = "Hit Count"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14, color = "#222222"),
      plot.subtitle = element_text(size = 10, color = "#555555"),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      panel.grid = element_blank(),
      axis.title = element_blank(),
      legend.position = "right",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9)
    )

  ggsave(
    filename = file.path(out_dir, "kegg_enrichment_circle_plot.pdf"),
    plot = circle_plot,
    width = 9,
    height = 9,
    dpi = 320,
    bg = "white"
  )
  ggsave(
    filename = file.path(out_dir, "kegg_enrichment_circle_plot.png"),
    plot = circle_plot,
    width = 9,
    height = 9,
    dpi = 320,
    bg = "white"
  )

  write.table(
    top_plot,
    file = file.path(out_dir, "kegg_enrichment_top20_for_plot.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  cat("  已输出: 气泡图与富集圈图（PDF + PNG）\n")
} else {
  cat("  enrich 结果为空，跳过绘图。\n")
}

cat("完成。输出目录: ", out_dir, "\n", sep = "")
