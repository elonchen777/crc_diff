library(ggplot2)
library(dplyr)
library(data.table)

if (!requireNamespace("ropls", quietly = TRUE)) {
  stop("缺少 ropls 包，请先安装: BiocManager::install('ropls')")
}

theme_microbiome <- function() {
  theme_bw(base_size = 14) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(linewidth = 0.8),
      axis.text = element_text(color = "black"),
      axis.title = element_text(size = 14),
      legend.background = element_blank(),
      legend.key = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(size = 13, face = "bold")
    )
}

cat("读取合并后的数据...\n")
merged_data <- fread("dataset/merged_dataset_relative.csv", stringsAsFactors = FALSE, data.table = FALSE)

output_dir <- file.path("results", "R_plots", "oplsda_plot")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

prepare_diff_groups <- function(df) {
  cat("根据分化程度和吸烟状态创建分组...\n")

  df$group <- sapply(1:nrow(df), function(i) {
    crc <- as.integer(df$crc_label[i])
    diff <- as.integer(df$differentiation[i])

    if (!is.na(crc) && crc == 1) {
      if (diff == 1) {
        return("CRC_poor_diff")
      } else if (diff == 0) {
        return("CRC_well_diff")
      } else {
        return(NA)
      }
    } else {
      return("control")
    }
  })

  df
}

grouped_data <- prepare_diff_groups(merged_data)

group_colors <- c(
  "control" = "#2E86AB",
  "CRC_well_diff" = "#F18F01",
  "CRC_poor_diff" = "#D7263D"
)

group_labels <- c(
  "control" = "CTRL",
  "CRC_well_diff" = "CRC-Well",
  "CRC_poor_diff" = "CRC-Poor"
)

group_levels <- c("control", "CRC_well_diff", "CRC_poor_diff")

cat("========== 开始代谢物两两 OPLS-DA 分析 ==========\n")

metabolite_cols <- grep("^met_", colnames(grouped_data), value = TRUE)
cat(sprintf("找到 %d 个代谢物特征\n", length(metabolite_cols)))

if (length(metabolite_cols) < 2) {
  stop("代谢物特征数量不足，无法进行 OPLS-DA")
}

metab_data <- grouped_data[, c("SAMPLE_ID", "group", metabolite_cols)]
metab_data <- metab_data[metab_data$group %in% group_levels, ]

if (nrow(metab_data) < 6) {
  stop("样本数过少，无法进行两两 OPLS-DA")
}

comparison_groups <- list(
  c("control", "CRC_well_diff"),
  c("control", "CRC_poor_diff"),
  c("CRC_well_diff", "CRC_poor_diff")
)

pairwise_summary <- list()

for (i in seq_along(comparison_groups)) {
  group1 <- comparison_groups[[i]][1]
  group2 <- comparison_groups[[i]][2]
  comparison_name <- paste(group1, "vs", group2, sep = "_")

  cat(sprintf("\n=== 代谢物 OPLS-DA 两两比较: %s ===\n", comparison_name))

  current_data <- metab_data[metab_data$group %in% c(group1, group2), ]
  current_data$group <- factor(current_data$group, levels = c(group1, group2))

  group_counts <- table(current_data$group)
  print(group_counts)

  if (any(group_counts < 3)) {
    cat("至少有一组样本量小于3，跳过该比较\n")
    next
  }

  metab_matrix <- as.matrix(current_data[, metabolite_cols])
  rownames(metab_matrix) <- current_data$SAMPLE_ID

  row_sums <- rowSums(metab_matrix, na.rm = TRUE)
  empty_rows <- row_sums == 0
  if (any(empty_rows, na.rm = TRUE)) {
    cat(sprintf("移除 %d 个全零样本\n", sum(empty_rows, na.rm = TRUE)))
    metab_matrix <- metab_matrix[!empty_rows, , drop = FALSE]
    current_data <- current_data[!empty_rows, ]
  }

  if (nrow(metab_matrix) < 6) {
    cat("过滤后样本数不足，跳过该比较\n")
    next
  }

  col_sums <- colSums(metab_matrix, na.rm = TRUE)
  zero_cols <- col_sums == 0
  if (any(zero_cols, na.rm = TRUE)) {
    cat(sprintf("移除 %d 个全零代谢物\n", sum(zero_cols, na.rm = TRUE)))
    metab_matrix <- metab_matrix[, !zero_cols, drop = FALSE]
  }

  if (ncol(metab_matrix) < 2) {
    cat("可用代谢物特征不足2个，跳过该比较\n")
    next
  }

  if (any(is.na(metab_matrix))) {
    cat("警告: 代谢物矩阵中存在 NA，将按列中位数填充\n")
    for (j in seq_len(ncol(metab_matrix))) {
      if (any(is.na(metab_matrix[, j]))) {
        metab_matrix[is.na(metab_matrix[, j]), j] <- median(metab_matrix[, j], na.rm = TRUE)
      }
    }
  }

  cat("构建 OPLS-DA 模型...\n")
  set.seed(123)

  model <- tryCatch(
    {
      ropls::opls(
        x = metab_matrix,
        y = current_data$group,
        predI = 1,
        orthoI = 1,
        permI = 200,
        scaleC = "standard",
        fig.pdfC = "none",
        info.txtC = "none"
      )
    },
    error = function(e) {
      cat("OPLS-DA(predI=1, orthoI=1)失败，回退到 PLS-DA(predI=2, orthoI=0) 用于二维作图\n")
      ropls::opls(
        x = metab_matrix,
        y = current_data$group,
        predI = 2,
        orthoI = 0,
        permI = 200,
        scaleC = "standard",
        fig.pdfC = "none",
        info.txtC = "none"
      )
    }
  )

  score_mn <- as.data.frame(ropls::getScoreMN(model))
  
  # 检查得分矩阵是否为空
  if (nrow(score_mn) == 0) {
    cat("模型得分矩阵为空，跳过该比较\n")
    next
  }

  x_comp <- colnames(score_mn)[1]
  y_comp <- if ("o1" %in% colnames(score_mn)) {
    "o1"
  } else if (ncol(score_mn) >= 2) {
    colnames(score_mn)[2]
  } else {
    # 如果只有一维得分（可能是由于 ropls 的某些模型特性），伪造一个 Y 轴用于绘图
    cat("模型只有一维得分，将添加随机抖动作为 Y 轴以进行二维可视化\n")
    "jitter_y"
  }

  score_df <- data.frame(
    SAMPLE_ID = rownames(score_mn),
    X = score_mn[[x_comp]],
    group = current_data$group[match(rownames(score_mn), current_data$SAMPLE_ID)]
  )
  
  if (y_comp == "jitter_y") {
    score_df$Y <- rnorm(nrow(score_df), 0, 0.05)
  } else {
    score_df$Y <- score_mn[[y_comp]]
  }

  summary_df <- as.data.frame(ropls::getSummaryDF(model))
  r2y <- if ("R2Y(cum)" %in% colnames(summary_df)) summary_df[1, "R2Y(cum)"] else NA
  q2 <- if ("Q2(cum)" %in% colnames(summary_df)) summary_df[1, "Q2(cum)"] else NA
  subtitle_text <- sprintf("R2Y(cum) = %.3f, Q2(cum) = %.3f", as.numeric(r2y), as.numeric(q2))

  active_levels <- c(group1, group2)

  p <- ggplot(score_df, aes(x = X, y = Y, color = group)) +
    geom_point(size = 3, alpha = 0.8) +
    stat_ellipse(level = 0.95, linetype = 1, linewidth = 1) +
    scale_color_manual(
      values = group_colors[active_levels],
      labels = group_labels[active_levels],
      drop = FALSE
    ) +
    labs(
      title = "OPLS-DA (Metabolites)",
      subtitle = subtitle_text,
      x = paste0(x_comp, " score"),
      y = paste0(y_comp, " score"),
      color = "Group"
    ) +
    theme_microbiome() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 10),
      legend.position = "right"
    )

  plot_path <- file.path(output_dir, paste0("OPLSDA_metabolites_", comparison_name, ".png"))
  ggsave(plot_path, p, width = 8, height = 6, dpi = 300)
  cat(sprintf("OPLS-DA图已保存: %s\n", plot_path))

  score_out <- file.path(output_dir, paste0("oplsda_metabolites_", comparison_name, "_scores.csv"))
  write.csv(score_df, score_out, row.names = FALSE)
  cat(sprintf("得分数据已保存: %s\n", score_out))

  vip_vec <- tryCatch(ropls::getVipVn(model), error = function(e) NULL)
  if (!is.null(vip_vec)) {
    vip_df <- data.frame(
      metabolite_name = names(vip_vec),
      VIP = as.numeric(vip_vec),
      stringsAsFactors = FALSE
    ) %>% arrange(desc(VIP))

    vip_out <- file.path(output_dir, paste0("oplsda_metabolites_", comparison_name, "_vip.csv"))
    write.csv(vip_df, vip_out, row.names = FALSE)
    cat(sprintf("VIP结果已保存: %s\n", vip_out))
  }

  summary_out <- file.path(output_dir, paste0("oplsda_metabolites_", comparison_name, "_model_summary.csv"))
  write.csv(summary_df, summary_out, row.names = TRUE)
  cat(sprintf("模型摘要已保存: %s\n", summary_out))

  pairwise_summary[[length(pairwise_summary) + 1]] <- data.frame(
    comparison = comparison_name,
    n_samples = nrow(score_df),
    n_features = ncol(metab_matrix),
    R2Y_cum = as.numeric(r2y),
    Q2_cum = as.numeric(q2),
    stringsAsFactors = FALSE
  )
}

if (length(pairwise_summary) > 0) {
  pairwise_summary_df <- do.call(rbind, pairwise_summary)
  summary_all_out <- file.path(output_dir, "oplsda_metabolites_pairwise_summary.csv")
  write.csv(pairwise_summary_df, summary_all_out, row.names = FALSE)
  cat(sprintf("两两比较汇总已保存: %s\n", summary_all_out))
} else {
  cat("未生成可用的两两比较结果，请检查样本量与特征质量\n")
}

cat("\n========== 代谢物两两 OPLS-DA 分析完成 ==========\n")
