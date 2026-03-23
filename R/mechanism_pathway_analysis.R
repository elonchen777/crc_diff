#!/usr/bin/env Rscript

# 机制通路分析脚本
# 目标: Microbe -> KO -> Metabolite, 可选扩展为 Microbe -> KO -> Metabolite -> Differentiation
#
# 输出:
# 1) results/pathway_mechanism/pathway_table.csv
# 2) results/pathway_mechanism/top_pathways.csv
# 3) results/pathway_mechanism/module_summary.csv
# 4) results/pathway_mechanism/hierarchical_network.png
# 5) results/pathway_mechanism/sankey_plot.png

suppressPackageStartupMessages({
  library(tidyverse)
  library(Hmisc)
  library(igraph)
  library(ggraph)
  library(ggalluvial)
})

# ------------------------------
# 参数配置
# ------------------------------
config <- list(
  maaslin_microbe = "results/maaslin2/significant_results.tsv",
  maaslin_ko = "results/maaslin2_ko/significant_results.tsv",
  maaslin_metabolite = "results/maaslin2_metabolite/significant_results.tsv",
  abundance_file = "dataset/merged_dataset_processed.csv",
  outdir = "results/pathway_mechanism",

  # Step 1: 从3个MaAsLin2结果中各取前n个特征
  top_n_microbe = 1000,
  top_n_ko = 1000,
  top_n_metabolite = 1000,

  enable_met_diff = TRUE,
  enable_sign_consistency = FALSE,

  rho_cutoff = 0.2,
  fdr_cutoff = 0.1,

  top_n = 100,
  bootstrap_n = 100,
  bootstrap_frac = 0.8,
  bootstrap_keep_freq = 0.7,
  module_top_k = 5,

  seed = 20260323
)

set.seed(config$seed)
dir.create(config$outdir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------
# 工具函数
# ------------------------------
stop_if_missing <- function(paths) {
  missing <- paths[!file.exists(paths)]
  if (length(missing) > 0) {
    stop(paste("以下文件不存在:", paste(missing, collapse = ", ")))
  }
}

read_maaslin_features <- function(path, top_n = 200) {
  df <- suppressMessages(readr::read_tsv(path, show_col_types = FALSE))
  required <- c("feature", "qval")
  if (!all(required %in% colnames(df))) {
    stop(sprintf("MaAsLin结果缺少必要列: %s", path))
  }

  ranked <- df %>% dplyr::filter(!is.na(feature), !is.na(qval))
  if ("pval" %in% colnames(ranked)) {
    ranked <- ranked %>% dplyr::arrange(qval, pval)
  } else {
    ranked <- ranked %>% dplyr::arrange(qval)
  }

  ranked %>%
    dplyr::slice_head(n = top_n) %>%
    dplyr::distinct(feature) %>%
    dplyr::pull(feature)
}

as_numeric_vector <- function(x) {
  y <- suppressWarnings(as.numeric(x))
  y
}

get_differentiation_column <- function(df) {
  candidates <- c("differentiation_stage", "diff_stage", "differentiation")
  hit <- candidates[candidates %in% colnames(df)]
  if (length(hit) == 0) {
    stop("丰度矩阵中未找到 differentiation_stage/diff_stage/differentiation 列")
  }
  hit[1]
}

build_diff_stage_from_crc_diff <- function(df) {
  if (!all(c("crc_label", "differentiation") %in% colnames(df))) {
    return(NULL)
  }

  crc <- suppressWarnings(as.integer(df$crc_label))
  diff <- suppressWarnings(as.integer(df$differentiation))

  # 参考分组函数:
  # control = 0, well-diff = 1 (differentiation == 0), poor-diff = 2 (differentiation == 1)
  out <- rep(NA_real_, nrow(df))
  out[!is.na(crc) & crc != 1] <- 0
  out[!is.na(crc) & crc == 1 & !is.na(diff) & diff == 0] <- 1
  out[!is.na(crc) & crc == 1 & !is.na(diff) & diff == 1] <- 2
  out
}

standardize_columns <- function(df) {
  scaled <- as.data.frame(scale(df))
  scaled[] <- lapply(scaled, function(x) {
    x[is.na(x) | is.infinite(x)] <- 0
    x
  })
  scaled
}

prepare_abundance_data <- function(path, microbe_features, ko_features, met_features, enable_met_diff = FALSE) {
  df <- suppressMessages(readr::read_csv(path, show_col_types = FALSE))

  if (!("SAMPLE_ID" %in% colnames(df))) {
    stop("丰度矩阵缺少 SAMPLE_ID 列")
  }

  diff_col <- NULL
  if (isTRUE(enable_met_diff)) {
    stage <- build_diff_stage_from_crc_diff(df)
    if (!is.null(stage)) {
      df$diff_stage <- stage
      diff_col <- "diff_stage"
      message("met-diff使用diff_stage编码: control=0, well-diff=1, poor-diff=2")
    } else {
      diff_col <- get_differentiation_column(df)
      message(sprintf("未找到crc_label+differentiation，回退使用列: %s", diff_col))
    }
  }

  keep_microbe <- intersect(microbe_features, colnames(df))
  keep_ko <- intersect(ko_features, colnames(df))
  keep_met <- intersect(met_features, colnames(df))

  message(sprintf("筛选后特征数: microbe=%d, ko=%d, metabolite=%d", length(keep_microbe), length(keep_ko), length(keep_met)))

  if (length(keep_microbe) == 0 || length(keep_ko) == 0 || length(keep_met) == 0) {
    stop("筛选后至少有一类特征为空，请检查MaAsLin结果和丰度矩阵前缀是否一致")
  }

  select_cols <- c(keep_microbe, keep_ko, keep_met)
  if (!is.null(diff_col)) {
    select_cols <- c(select_cols, diff_col)
  }

  sub <- df %>%
    dplyr::select(SAMPLE_ID, dplyr::all_of(select_cols)) %>%
    dplyr::distinct(SAMPLE_ID, .keep_all = TRUE)

  # 统一样本顺序与可比性: 仅保留在所有分析中都有效的样本
  # 这里按完整案例筛选，避免不同相关分析使用不同样本集合
  numeric_part <- sub %>% dplyr::select(-SAMPLE_ID)
  numeric_part[] <- lapply(numeric_part, as_numeric_vector)

  complete_idx <- complete.cases(numeric_part)
  cleaned <- dplyr::bind_cols(sub %>% dplyr::select(SAMPLE_ID), numeric_part) %>%
    dplyr::filter(complete_idx)

  feature_cols <- setdiff(colnames(cleaned), c("SAMPLE_ID", diff_col))
  cleaned[, feature_cols] <- standardize_columns(cleaned[, feature_cols, drop = FALSE])

  if (nrow(cleaned) < 20) {
    stop("完整案例样本数过少(<20)，无法稳定进行相关分析")
  }

  list(
    data = cleaned,
    microbe = keep_microbe,
    ko = keep_ko,
    metabolite = keep_met,
    diff_col = diff_col
  )
}

calc_pairwise_spearman <- function(df, x_cols, y_cols, x_name, y_name) {
  mat_x <- as.matrix(df[, x_cols, drop = FALSE])
  mat_y <- as.matrix(df[, y_cols, drop = FALSE])

  combined <- cbind(mat_x, mat_y)
  rc <- Hmisc::rcorr(combined, type = "spearman")

  nx <- ncol(mat_x)
  ny <- ncol(mat_y)

  rho <- rc$r[seq_len(nx), nx + seq_len(ny), drop = FALSE]
  p <- rc$P[seq_len(nx), nx + seq_len(ny), drop = FALSE]

  out <- expand.grid(
    x = colnames(mat_x),
    y = colnames(mat_y),
    stringsAsFactors = FALSE
  ) %>%
    dplyr::mutate(
      rho = as.vector(rho),
      pval = as.vector(p)
    ) %>%
    dplyr::mutate(fdr = p.adjust(pval, method = "fdr")) %>%
    dplyr::rename(!!x_name := x, !!y_name := y)

  out
}

filter_significant_edges <- function(df, rho_cutoff, fdr_cutoff) {
  df %>%
    dplyr::filter(!is.na(rho), !is.na(fdr)) %>%
    dplyr::filter(abs(rho) > rho_cutoff, fdr < fdr_cutoff)
}

build_pathways <- function(microbe_ko_sig, ko_met_sig, met_diff_sig = NULL, enable_sign_consistency = TRUE) {
  pathways <- microbe_ko_sig %>%
    dplyr::rename(r1 = rho, p1 = pval, fdr1 = fdr) %>%
    dplyr::inner_join(
      ko_met_sig %>% dplyr::rename(r2 = rho, p2 = pval, fdr2 = fdr),
      by = "ko"
    ) %>%
    dplyr::mutate(
      score = abs(r1) * abs(r2)
    )

  if (!is.null(met_diff_sig)) {
    pathways <- pathways %>%
      dplyr::inner_join(
        met_diff_sig %>% dplyr::rename(r3 = rho, p3 = pval, fdr3 = fdr),
        by = "metabolite"
      ) %>%
      dplyr::mutate(
        sign_consistent = if (isTRUE(enable_sign_consistency)) {
          sign(r1) == sign(r2) & sign(r2) == sign(r3)
        } else {
          TRUE
        },
        score = abs(r1) * abs(r2) * abs(r3)
      )
  } else {
    pathways <- pathways %>%
      dplyr::mutate(sign_consistent = if (isTRUE(enable_sign_consistency)) sign(r1) == sign(r2) else TRUE)
  }

  pathways
}

compute_bootstrap_frequency <- function(df, pathways, diff_col = NULL, n_boot = 100, frac = 0.8,
                                        rho_cutoff = 0.3, fdr_cutoff = 0.05,
                                        enable_sign_consistency = TRUE) {
  if (nrow(pathways) == 0) {
    return(tibble::tibble(path_id = character(), bootstrap_freq = numeric()))
  }

  candidate <- pathways %>%
    dplyr::mutate(path_id = paste(microbe, ko, metabolite, sep = "||")) %>%
    dplyr::select(path_id, microbe, ko, metabolite) %>%
    dplyr::distinct()

  mk_pairs <- candidate %>% dplyr::distinct(microbe, ko)
  km_pairs <- candidate %>% dplyr::distinct(ko, metabolite)
  met_list <- candidate %>% dplyr::distinct(metabolite)

  n <- nrow(df)
  b_size <- max(10, floor(n * frac))

  hit_count <- stats::setNames(rep(0L, nrow(candidate)), candidate$path_id)

  pair_corr <- function(dat, pairs, x_col, y_col) {
    if (nrow(pairs) == 0) {
      return(tibble::tibble())
    }

    res <- purrr::pmap_dfr(pairs, function(...) {
      row <- list(...)
      x <- as.numeric(dat[[row[[x_col]]]])
      y <- as.numeric(dat[[row[[y_col]]]])

      tst <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
      tibble::tibble(
        !!x_col := row[[x_col]],
        !!y_col := row[[y_col]],
        rho = unname(tst$estimate),
        pval = tst$p.value
      )
    })

    res %>% dplyr::mutate(fdr = p.adjust(pval, method = "fdr"))
  }

  pair_corr_single <- function(dat, x_cols, y_col) {
    if (length(x_cols) == 0) {
      return(tibble::tibble())
    }

    res <- purrr::map_dfr(x_cols, function(x_col) {
      x <- as.numeric(dat[[x_col]])
      y <- as.numeric(dat[[y_col]])

      tst <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
      tibble::tibble(
        metabolite = x_col,
        rho = unname(tst$estimate),
        pval = tst$p.value
      )
    })

    res %>% dplyr::mutate(fdr = p.adjust(pval, method = "fdr"))
  }

  for (b in seq_len(n_boot)) {
    idx <- sample(seq_len(n), size = b_size, replace = TRUE)
    dat_b <- df[idx, , drop = FALSE]

    mk_b <- pair_corr(dat_b, mk_pairs, "microbe", "ko") %>%
      dplyr::filter(abs(rho) > rho_cutoff, fdr < fdr_cutoff) %>%
      dplyr::rename(r1 = rho)

    km_b <- pair_corr(dat_b, km_pairs, "ko", "metabolite") %>%
      dplyr::filter(abs(rho) > rho_cutoff, fdr < fdr_cutoff) %>%
      dplyr::rename(r2 = rho)

    if (!is.null(diff_col)) {
      md_b <- pair_corr_single(dat_b, unique(met_list$metabolite), diff_col) %>%
        dplyr::filter(abs(rho) > rho_cutoff, fdr < fdr_cutoff) %>%
        dplyr::rename(r3 = rho)
    }

    hit_b <- candidate %>%
      dplyr::inner_join(mk_b %>% dplyr::select(microbe, ko, r1), by = c("microbe", "ko")) %>%
      dplyr::inner_join(km_b %>% dplyr::select(ko, metabolite, r2), by = c("ko", "metabolite")) %>%
      {
        if (!is.null(diff_col)) {
          tmp <- dplyr::inner_join(., md_b %>% dplyr::select(metabolite, r3), by = "metabolite")
          if (isTRUE(enable_sign_consistency)) {
            dplyr::filter(tmp, sign(r1) == sign(r2), sign(r2) == sign(r3))
          } else {
            tmp
          }
        } else {
          if (isTRUE(enable_sign_consistency)) dplyr::filter(., sign(r1) == sign(r2)) else .
        }
      } %>%
      dplyr::pull(path_id)

    if (length(hit_b) > 0) {
      hit_count[hit_b] <- hit_count[hit_b] + 1L
    }
  }

  tibble::tibble(
    path_id = names(hit_count),
    bootstrap_freq = as.numeric(hit_count) / n_boot
  )
}

infer_module_label <- function(metabolites, kos) {
  txt <- paste(c(metabolites, kos), collapse = " ") %>% tolower()

  if (str_detect(txt, "bile|cholic|chenodeoxy|tauro|glyco")) return("Bile acid related")
  if (str_detect(txt, "lipid|fatty|carnitine|cer|sphingo|phospho|chol|linole|arachid|prostagland")) return("Lipid metabolism")
  if (str_detect(txt, "amino|gly|ala|val|leu|ile|phe|trp|tyr|his|orn|glu|gln")) return("Amino acid related")
  if (str_detect(txt, "nucleotide|purine|pyrimid|inosine|uridine|cytidine")) return("Nucleotide related")
  return("Other")
}

build_modules <- function(pathways, top_k = 5) {
  if (nrow(pathways) == 0) {
    return(list(pathways = pathways, summary = tibble::tibble()))
  }

  pw <- pathways %>%
    dplyr::mutate(path_id = dplyr::row_number()) %>%
    dplyr::select(path_id, microbe, ko, metabolite, score)

  # 若两条路径共享KO或代谢物，则连接到同一模块
  edges_ko <- pw %>%
    dplyr::select(path_id, ko) %>%
    dplyr::inner_join(., ., by = "ko", suffix = c("_a", "_b")) %>%
    dplyr::filter(path_id_a < path_id_b) %>%
    dplyr::transmute(from = as.character(path_id_a), to = as.character(path_id_b))

  edges_met <- pw %>%
    dplyr::select(path_id, metabolite) %>%
    dplyr::inner_join(., ., by = "metabolite", suffix = c("_a", "_b")) %>%
    dplyr::filter(path_id_a < path_id_b) %>%
    dplyr::transmute(from = as.character(path_id_a), to = as.character(path_id_b))

  comp_membership <- if (nrow(edges_ko) + nrow(edges_met) == 0) {
    stats::setNames(seq_len(nrow(pw)), as.character(pw$path_id))
  } else {
    g <- igraph::graph_from_data_frame(dplyr::bind_rows(edges_ko, edges_met), directed = FALSE,
                                       vertices = tibble::tibble(name = as.character(pw$path_id)))
    comps <- igraph::components(g)
    comps$membership
  }

  pw_mod <- pw %>%
    dplyr::mutate(module_id = as.integer(comp_membership[as.character(path_id)]))

  module_summary <- pw_mod %>%
    dplyr::group_by(module_id) %>%
    dplyr::summarise(
      n_pathways = dplyr::n(),
      mean_score = mean(score, na.rm = TRUE),
      top_kos = paste(head(sort(table(ko), decreasing = TRUE), 5) %>% names(), collapse = "; "),
      top_metabolites = paste(head(sort(table(metabolite), decreasing = TRUE), 5) %>% names(), collapse = "; "),
      .groups = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(n_pathways), dplyr::desc(mean_score)) %>%
    dplyr::slice_head(n = top_k) %>%
    dplyr::mutate(
      module_label = purrr::map2_chr(top_metabolites, top_kos, ~ infer_module_label(.x, .y))
    )

  pw_mod <- pw_mod %>%
    dplyr::inner_join(module_summary %>% dplyr::select(module_id), by = "module_id")

  list(pathways = pw_mod, summary = module_summary)
}

plot_hierarchical_network <- function(pathways, outfile) {
  if (nrow(pathways) == 0) {
    message("无可绘制路径，跳过分层网络图")
    return(invisible(NULL))
  }

  edges <- dplyr::bind_rows(
    pathways %>% dplyr::transmute(from = microbe, to = ko, weight = abs(r1), sign = sign(r1), edge_type = "microbe_ko"),
    pathways %>% dplyr::transmute(from = ko, to = metabolite, weight = abs(r2), sign = sign(r2), edge_type = "ko_metabolite")
  ) %>%
    dplyr::distinct()

  nodes <- tibble::tibble(name = unique(c(edges$from, edges$to))) %>%
    dplyr::mutate(
      node_type = dplyr::case_when(
        str_starts(name, "tax_") ~ "Microbe",
        str_starts(name, "kegg_") ~ "KO",
        str_starts(name, "met_") ~ "Metabolite",
        TRUE ~ "Other"
      )
    )

  g <- igraph::graph_from_data_frame(edges, directed = TRUE, vertices = nodes)

  p <- ggraph::ggraph(g, layout = "sugiyama") +
    ggraph::geom_edge_link(aes(width = weight, color = sign), alpha = 0.6) +
    ggraph::geom_node_point(aes(color = node_type), size = 3) +
    ggraph::geom_node_text(aes(label = name), size = 2.2, repel = TRUE) +
    scale_edge_color_gradient2(low = "#2c7fb8", mid = "#666666", high = "#d7191c", midpoint = 0) +
    scale_edge_width(range = c(0.2, 1.8)) +
    theme_graph(base_family = "sans") +
    labs(
      title = "Hierarchical Network: Microbe -> KO -> Metabolite",
      edge_color = "Correlation Sign",
      edge_width = "|rho|",
      color = "Node Type"
    )

  ggplot2::ggsave(outfile, p, width = 16, height = 9, dpi = 300)
}

plot_sankey <- function(pathways, outfile) {
  if (nrow(pathways) == 0) {
    message("无可绘制路径，跳过Sankey图")
    return(invisible(NULL))
  }

  sankey_df <- pathways %>%
    dplyr::mutate(path_id = dplyr::row_number()) %>%
    dplyr::select(path_id, microbe, ko, metabolite, score) %>%
    dplyr::mutate(
      microbe = as.factor(microbe),
      ko = as.factor(ko),
      metabolite = as.factor(metabolite)
    )

  p <- ggplot2::ggplot(
    sankey_df,
    aes(axis1 = microbe, axis2 = ko, axis3 = metabolite, y = score)
  ) +
    ggalluvial::geom_alluvium(aes(fill = score), alpha = 0.7, width = 0.15) +
    ggalluvial::geom_stratum(width = 0.15, color = "grey30", fill = "grey90") +
    ggplot2::geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 2.2) +
    scale_x_discrete(limits = c("Microbe", "KO", "Metabolite"), expand = c(0.08, 0.08)) +
    scale_fill_viridis_c(option = "C") +
    theme_minimal(base_family = "sans") +
    labs(title = "Sankey: Microbe -> KO -> Metabolite", y = "Pathway Score", x = NULL, fill = "Score")

  ggplot2::ggsave(outfile, p, width = 16, height = 9, dpi = 300)
}

# ------------------------------
# 主流程
# ------------------------------
message("Step 0: 检查输入文件...")
stop_if_missing(c(
  config$maaslin_microbe,
  config$maaslin_ko,
  config$maaslin_metabolite,
  config$abundance_file
))

message("Step 1: 从MaAsLin2结果按显著性排序提取前n个特征...")
microbe_features <- read_maaslin_features(config$maaslin_microbe, config$top_n_microbe)
ko_features <- read_maaslin_features(config$maaslin_ko, config$top_n_ko)
met_features <- read_maaslin_features(config$maaslin_metabolite, config$top_n_metabolite)

message(sprintf("Step 1完成: microbe_top=%d, ko_top=%d, metabolite_top=%d",
                length(microbe_features), length(ko_features), length(met_features)))

message("Step 2: 读取丰度矩阵并统一样本顺序...")
prep <- prepare_abundance_data(
  config$abundance_file,
  microbe_features,
  ko_features,
  met_features,
  enable_met_diff = config$enable_met_diff
)
abd <- prep$data

message(sprintf("有效样本数: %d", nrow(abd)))

message("Step 3: Spearman相关分析与FDR校正（Microbe-KO / KO-Metabolite）...")
microbe_ko <- calc_pairwise_spearman(abd, prep$microbe, prep$ko, "microbe", "ko")
ko_met <- calc_pairwise_spearman(abd, prep$ko, prep$metabolite, "ko", "metabolite")
met_diff <- NULL

if (isTRUE(config$enable_met_diff)) {
  message("Step 3补充: 启用met-diff相关性分析（diff编码: ctrl=0, well=1, poor=2）...")
  met_diff <- calc_pairwise_spearman(abd, prep$metabolite, prep$diff_col, "metabolite", "differentiation")
}

microbe_ko_sig <- filter_significant_edges(microbe_ko, config$rho_cutoff, config$fdr_cutoff)
ko_met_sig <- filter_significant_edges(ko_met, config$rho_cutoff, config$fdr_cutoff)
met_diff_sig <- NULL

if (!is.null(met_diff)) {
  met_diff_sig <- filter_significant_edges(met_diff, config$rho_cutoff, config$fdr_cutoff)
}

if (is.null(met_diff_sig)) {
  message(sprintf("显著相关边数: microbe-ko=%d, ko-met=%d",
                  nrow(microbe_ko_sig), nrow(ko_met_sig)))
} else {
  message(sprintf("显著相关边数: microbe-ko=%d, ko-met=%d, met-diff=%d",
                  nrow(microbe_ko_sig), nrow(ko_met_sig), nrow(met_diff_sig)))
}

message("Step 4-6: 路径构建 + 方向一致性筛选 + 打分...")
pathway_table <- build_pathways(
  microbe_ko_sig,
  ko_met_sig,
  met_diff_sig,
  enable_sign_consistency = config$enable_sign_consistency
)

pathway_sig <- pathway_table %>% dplyr::arrange(dplyr::desc(score))
if (isTRUE(config$enable_sign_consistency)) {
  pathway_sig <- pathway_sig %>% dplyr::filter(sign_consistent)
}

top_pathways <- pathway_sig %>% dplyr::slice_head(n = config$top_n)

message(sprintf("方向一致路径数: %d", nrow(pathway_sig)))
message(sprintf("Top路径数: %d", nrow(top_pathways)))

message("Step 7: Bootstrap稳定性分析...")
boot_freq <- compute_bootstrap_frequency(
  abd,
  top_pathways,
  diff_col = prep$diff_col,
  n_boot = config$bootstrap_n,
  frac = config$bootstrap_frac,
  rho_cutoff = config$rho_cutoff,
  fdr_cutoff = config$fdr_cutoff,
  enable_sign_consistency = config$enable_sign_consistency
)

if (nrow(top_pathways) > 0) {
  top_pathways <- top_pathways %>%
    dplyr::mutate(path_id = paste(microbe, ko, metabolite, sep = "||")) %>%
    dplyr::left_join(boot_freq, by = "path_id") %>%
    dplyr::mutate(bootstrap_freq = tidyr::replace_na(bootstrap_freq, 0)) %>%
    dplyr::select(-path_id)
}

stable_top <- top_pathways %>%
  dplyr::filter(bootstrap_freq >= config$bootstrap_keep_freq) %>%
  dplyr::arrange(dplyr::desc(score), dplyr::desc(bootstrap_freq))

message(sprintf("稳定路径数(频率>=%.2f): %d", config$bootstrap_keep_freq, nrow(stable_top)))

message("Step 8: 模块构建...")
module_input <- if (nrow(stable_top) > 0) stable_top else top_pathways
modules <- build_modules(module_input, top_k = config$module_top_k)
module_pathways <- modules$pathways
module_summary <- modules$summary

message(sprintf("模块数(输出前%d): %d", config$module_top_k, nrow(module_summary)))

message("Step 9: 导出结果...")
readr::write_csv(pathway_table, file.path(config$outdir, "pathway_table.csv"))
readr::write_csv(
  if (nrow(stable_top) > 0) stable_top else top_pathways,
  file.path(config$outdir, "top_pathways.csv")
)
readr::write_csv(module_summary, file.path(config$outdir, "module_summary.csv"))

message("Step 10: 生成可视化...")
plot_input <- module_input

plot_hierarchical_network(
  plot_input,
  file.path(config$outdir, "hierarchical_network.png")
)
plot_sankey(
  plot_input,
  file.path(config$outdir, "sankey_plot.png")
)

message("分析完成。输出目录:")
message(config$outdir)
