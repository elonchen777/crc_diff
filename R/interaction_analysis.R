# 加载必要的库
library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(data.table)
library(patchwork)

# 1. 读取数据
cat("读取合并后的数据...\n")
merged_data <- fread("dataset/merged_dataset.csv", stringsAsFactors = FALSE, data.table = FALSE)
output_dir <- file.path("results", "interaction_analysis")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# 2. 准备数据 - 只保留CRC患者
cat("筛选CRC患者数据...\n")
crc_data <- merged_data[!is.na(merged_data$crc_label) & merged_data$crc_label == 1, ]
cat(sprintf("CRC患者数量: %d\n", nrow(crc_data)))

# 3. 处理TNM分期和分化程度
crc_data$tnm_stage <- as.factor(crc_data$tnm_stage)
crc_data$differentiation <- factor(crc_data$differentiation, levels = c(0, 1), labels = c("well_diff", "poor_diff"))
crc_data$age <- as.numeric(crc_data$age)
crc_data$gender <- factor(crc_data$gender_label, levels = c(0, 1), labels = c("Female", "Male"))

cat("TNM分期分布:\n")
print(table(crc_data$tnm_stage))
cat("\n分化程度分布:\n")
print(table(crc_data$differentiation))

# 4. 提取宏基因组特征
taxonomy_cols <- grep("^tax_", colnames(crc_data), value = TRUE)
cat(sprintf("\n找到 %d 个宏基因组特征\n", length(taxonomy_cols)))

# 5. 提取代谢组特征
metabolite_cols <- grep("^met_", colnames(crc_data), value = TRUE)
cat(sprintf("找到 %d 个代谢组特征\n", length(metabolite_cols)))

# 6. 交互效应分析函数
run_interaction_analysis <- function(data_df, feature_cols, data_type, output_prefix) {
  cat(sprintf("\n========================================\n"))
  cat(sprintf("开始 %s 交互效应分析\n", data_type))
  cat(sprintf("========================================\n"))
  
  results_list <- list()
  
  cat("数据预处理...\n")
  
  feature_matrix <- as.matrix(data_df[, feature_cols])
  rownames(feature_matrix) <- data_df$SAMPLE_ID
  
  filtered_features <- colnames(feature_matrix)
  
  cat(sprintf("\n开始对 %d 个特征进行交互效应分析...\n", length(filtered_features)))
  
  results <- data.frame(
    feature = filtered_features,
    tnm_stage_coef = NA,
    tnm_stage_pvalue = NA,
    differentiation_coef = NA,
    differentiation_pvalue = NA,
    interaction_coef = NA,
    interaction_pvalue = NA,
    age_coef = NA,
    age_pvalue = NA,
    gender_coef = NA,
    gender_pvalue = NA,
    r_squared = NA,
    adj_r_squared = NA,
    stringsAsFactors = FALSE
  )
  
  pb <- txtProgressBar(min = 0, max = length(filtered_features), style = 3)
  
  for (i in seq_along(filtered_features)) {
    feat <- filtered_features[i]
    y <- feature_matrix[, feat]
    
    model_data <- data.frame(
      y = y,
      tnm_stage = data_df$tnm_stage,
      differentiation = data_df$differentiation,
      age = data_df$age,
      gender = data_df$gender
    )

    model_data <- model_data[complete.cases(model_data), ]

    if (nrow(model_data) < 10) {
      setTxtProgressBar(pb, i)
      next
    }

    fit <- try(lm(log2(y + 1e-5) ~ tnm_stage + differentiation + tnm_stage:differentiation + age + gender, 
                  data = model_data), silent = TRUE)
    
    if (!inherits(fit, "try-error")) {
      coef_sum <- try(summary(fit)$coefficients, silent = TRUE)
      
      if (!inherits(coef_sum, "try-error") && nrow(coef_sum) > 0) {
        r2 <- summary(fit)$r.squared
        adj_r2 <- summary(fit)$adj.r.squared
        results$r_squared[i] <- r2
        results$adj_r_squared[i] <- adj_r2
        
        for (term in rownames(coef_sum)) {
          if (grepl("tnm_stage", term) && !grepl("differentiation", term)) {
            results$tnm_stage_coef[i] <- coef_sum[term, "Estimate"]
            results$tnm_stage_pvalue[i] <- coef_sum[term, "Pr(>|t|)"]
          } else if (grepl("differentiation", term) && !grepl("tnm_stage", term)) {
            results$differentiation_coef[i] <- coef_sum[term, "Estimate"]
            results$differentiation_pvalue[i] <- coef_sum[term, "Pr(>|t|)"]
          } else if (grepl("tnm_stage.*differentiation|differentiation.*tnm_stage", term)) {
            results$interaction_coef[i] <- coef_sum[term, "Estimate"]
            results$interaction_pvalue[i] <- coef_sum[term, "Pr(>|t|)"]
          } else if (term == "age") {
            results$age_coef[i] <- coef_sum[term, "Estimate"]
            results$age_pvalue[i] <- coef_sum[term, "Pr(>|t|)"]
          } else if (grepl("gender", term)) {
            results$gender_coef[i] <- coef_sum[term, "Estimate"]
            results$gender_pvalue[i] <- coef_sum[term, "Pr(>|t|)"]
          }
        }
      }
    }
    
    setTxtProgressBar(pb, i)
  }
  
  close(pb)
  
  cat("\n应用多重检验校正...\n")
  results$tnm_stage_padj <- p.adjust(results$tnm_stage_pvalue, method = "BH")
  results$differentiation_padj <- p.adjust(results$differentiation_pvalue, method = "BH")
  results$interaction_padj <- p.adjust(results$interaction_pvalue, method = "BH")
  
  results$interaction_significant <- results$interaction_padj < 0.05 & !is.na(results$interaction_padj)
  
  sig_interaction <- results[results$interaction_significant, ]
  sig_interaction <- sig_interaction[order(sig_interaction$interaction_padj), ]
  
  cat(sprintf("\n找到 %d 个具有显著交互效应的特征\n", nrow(sig_interaction)))
  
  result_file <- file.path(output_dir, paste0(output_prefix, "_interaction_results.csv"))
  write.csv(results, result_file, row.names = FALSE)
  cat(sprintf("完整结果已保存: %s\n", result_file))
  
  if (nrow(sig_interaction) > 0) {
    sig_file <- file.path(output_dir, paste0(output_prefix, "_significant_interaction.csv"))
    write.csv(sig_interaction, sig_file, row.names = FALSE)
    cat(sprintf("显著交互效应结果已保存: %s\n", sig_file))
    
    cat("\n前10个显著交互效应特征:\n")
    print(head(sig_interaction[, c("feature", "interaction_coef", "interaction_pvalue", "interaction_padj")], 10))
  }
  
  return(results)
}

# 7. 对物种数据进行交互效应分析
species_results <- run_interaction_analysis(
  data_df = crc_data,
  feature_cols = taxonomy_cols,
  data_type = "宏基因组",
  output_prefix = "species"
)

# 8. 对代谢物数据进行交互效应分析
metabolite_results <- run_interaction_analysis(
  data_df = crc_data,
  feature_cols = metabolite_cols,
  data_type = "代谢组",
  output_prefix = "metabolite"
)

# 9. 交互效应可视化函数（斜率图）
plot_interaction_effect <- function(data_df, feature_name, data_type, output_dir) {
  
  if (!feature_name %in% colnames(data_df)) {
    cat(sprintf("特征 %s 不存在\n", feature_name))
    return(NULL)
  }
  
  plot_data <- data.frame(
    abundance = as.numeric(data_df[[feature_name]]),
    tnm_stage = data_df$tnm_stage,
    differentiation = data_df$differentiation,
    SAMPLE_ID = data_df$SAMPLE_ID
  )
  
  plot_data <- plot_data[complete.cases(plot_data), ]
  
  plot_data$stage_group <- ifelse(plot_data$tnm_stage %in% c("0", "1", 0, 1), "Early", "Late")
  plot_data$stage_group <- factor(plot_data$stage_group, levels = c("Early", "Late"))
  
  summary_data <- plot_data %>%
    group_by(stage_group, differentiation) %>%
    summarise(
      mean_abundance = mean(abundance, na.rm = TRUE),
      se = sd(abundance, na.rm = TRUE) / sqrt(n()),
      n = n(),
      .groups = "drop"
    )
  
  diff_colors <- c("well_diff" = "#2E86AB", "poor_diff" = "#A23B72")
  
  p <- ggplot(summary_data, aes(x = stage_group, y = mean_abundance, 
                                 color = differentiation, group = differentiation)) +
    geom_line(size = 1.2) +
    geom_point(size = 4) +
    geom_errorbar(aes(ymin = mean_abundance - se, ymax = mean_abundance + se),
                  width = 0.1, size = 0.8) +
    scale_color_manual(values = diff_colors, name = "分化程度",
                       labels = c("well_diff" = "中-高分化", "poor_diff" = "低分化")) +
    labs(title = paste0("交互效应图: ", feature_name),
         subtitle = paste0("数据类型: ", data_type),
         x = "CRC分期",
         y = "平均丰度") +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 11),
      axis.text = element_text(size = 11),
      axis.title = element_text(size = 12),
      legend.position = "top",
      legend.title = element_text(size = 11),
      legend.text = element_text(size = 10),
      panel.grid.minor = element_blank()
    ) +
    annotate("text", x = 1.5, y = max(summary_data$mean_abundance) * 1.1,
             label = "斜率差异 = 交互效应", size = 3.5, color = "gray30")
  
  return(p)
}

# 10. 批量生成交互效应可视化
generate_interaction_plots <- function(data_df, results, data_type, output_prefix, top_n = 10) {
  cat(sprintf("\n生成 %s 交互效应可视化...\n", data_type))
  
  sig_results <- results[results$interaction_significant & !is.na(results$interaction_coef), ]
  
  if (nrow(sig_results) == 0) {
    cat("没有显著交互效应特征，跳过可视化\n")
    return(NULL)
  }
  
  top_features <- sig_results[order(abs(sig_results$interaction_coef), decreasing = TRUE), ]
  top_features <- top_features[1:min(top_n, nrow(top_features)), ]
  top_features$feature <- as.character(top_features$feature)
  
  plot_dir <- file.path(output_dir, paste0(output_prefix, "_interaction_plots"))
  if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)
  
  for (i in 1:nrow(top_features)) {
    feat <- top_features$feature[i]
    
    p <- plot_interaction_effect(data_df, feat, data_type, plot_dir)
    
    if (!is.null(p)) {
      safe_name <- gsub("[^a-zA-Z0-9_]", "_", feat)
      plot_file <- file.path(plot_dir, paste0(output_prefix, "_interaction_", safe_name, ".png"))
      ggsave(plot_file, p, width = 8, height = 6, dpi = 300)
    }
  }
  
  cat(sprintf("交互效应图已保存到: %s\n", plot_dir))
  
  if (nrow(top_features) >= 4) {
    plot_list <- list()
    for (i in 1:min(4, nrow(top_features))) {
      feat <- top_features$feature[i]
      p <- plot_interaction_effect(data_df, feat, data_type, output_dir)
      if (!is.null(p)) {
        plot_list[[i]] <- p + 
          labs(title = feat, subtitle = NULL) +
          theme(plot.title = element_text(size = 10))
      }
    }
    
    if (length(plot_list) > 0) {
      library(patchwork)
      combined_plot <- wrap_plots(plot_list, ncol = 2) +
        plot_annotation(
          title = paste0("Top 4 显著交互效应特征 (", data_type, ")"),
          theme = theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
        )
      
      combined_file <- file.path(output_dir, paste0(output_prefix, "_interaction_combined.png"))
      ggsave(combined_file, combined_plot, width = 14, height = 12, dpi = 300)
      cat(sprintf("组合交互效应图已保存: %s\n", combined_file))
    }
  }
}

# 11. 其他可视化函数
plot_interaction_results <- function(results, data_type, output_prefix) {
  cat(sprintf("\n生成 %s 可视化结果...\n", data_type))
  
  sig_results <- results[results$interaction_significant & !is.na(results$interaction_coef), ]
  
  if (nrow(sig_results) == 0) {
    cat("没有显著交互效应特征，跳过可视化\n")
    return(NULL)
  }
  
  top_n <- min(20, nrow(sig_results))
  top_features <- sig_results[order(abs(sig_results$interaction_coef), decreasing = TRUE), ]
  top_features <- top_features[1:top_n, ]
  top_features$feature <- as.character(top_features$feature)
  top_features$interaction_coef <- as.numeric(top_features$interaction_coef)
  top_features <- top_features[!is.na(top_features$feature) & !is.na(top_features$interaction_coef), ]
  
  if (nrow(top_features) == 0) {
    cat("没有有效数据用于可视化\n")
    return(NULL)
  }
  
  top_features$direction <- ifelse(top_features$interaction_coef > 0, "Positive", "Negative")
  
  fill_colors <- c("Positive" = "#E74C3C", "Negative" = "#3498DB")
  
  p <- ggplot(top_features, aes(x = reorder(feature, interaction_coef), y = interaction_coef, fill = direction)) +
    geom_bar(stat = "identity", width = 0.7) +
    coord_flip() +
    scale_fill_manual(values = fill_colors, name = "交互效应方向") +
    labs(title = paste0("TNM分期 x 分化 交互效应 (", data_type, ")"),
         subtitle = paste0("Top ", top_n, " 显著交互效应特征"),
         x = "特征",
         y = "交互效应系数") +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      axis.text.y = element_text(size = 9),
      axis.text.x = element_text(size = 10),
      legend.position = "top",
      panel.grid.minor = element_blank()
    )
  
  plot_file <- file.path(output_dir, paste0(output_prefix, "_interaction_barplot.png"))
  ggsave(plot_file, p, width = 12, height = 10, dpi = 300)
  cat(sprintf("交互效应柱状图已保存: %s\n", plot_file))
  
  volcano_data <- results[!is.na(results$interaction_coef) & !is.na(results$interaction_pvalue), ]
  volcano_data$log10p <- -log10(volcano_data$interaction_pvalue)
  volcano_data$significance <- "Not Significant"
  volcano_data$significance[volcano_data$interaction_padj < 0.05 & volcano_data$interaction_coef > 0] <- "Positive"
  volcano_data$significance[volcano_data$interaction_padj < 0.05 & volcano_data$interaction_coef < 0] <- "Negative"
  
  volcano_colors <- c("Positive" = "#E74C3C", "Negative" = "#3498DB", "Not Significant" = "#95A5A6")
  
  p_volcano <- ggplot(volcano_data, aes(x = interaction_coef, y = log10p, color = significance)) +
    geom_point(size = 2, alpha = 0.6) +
    scale_color_manual(values = volcano_colors) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50") +
    geom_vline(xintercept = 0, linetype = "solid", color = "gray30") +
    labs(title = paste0("交互效应火山图 (", data_type, ")"),
         subtitle = sprintf("显著交互效应: Positive = %d, Negative = %d",
                           sum(volcano_data$significance == "Positive"),
                           sum(volcano_data$significance == "Negative")),
         x = "交互效应系数",
         y = "-log10(p-value)",
         color = "显著性") +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 10),
      legend.position = "right"
    )
  
  volcano_file <- file.path(output_dir, paste0(output_prefix, "_interaction_volcano.png"))
  ggsave(volcano_file, p_volcano, width = 10, height = 8, dpi = 300)
  cat(sprintf("交互效应火山图已保存: %s\n", volcano_file))
  
  effect_summary <- data.frame(
    effect_type = c("TNM分期主效应", "分化主效应", "TNM x 分化交互"),
    n_significant = c(
      sum(results$tnm_stage_padj < 0.05, na.rm = TRUE),
      sum(results$differentiation_padj < 0.05, na.rm = TRUE),
      sum(results$interaction_padj < 0.05, na.rm = TRUE)
    )
  )
  
  p_summary <- ggplot(effect_summary, aes(x = effect_type, y = n_significant, fill = effect_type)) +
    geom_bar(stat = "identity", width = 0.6) +
    scale_fill_manual(values = c("#2E86AB", "#A23B72", "#F18F01")) +
    labs(title = paste0("显著效应数量统计 (", data_type, ")"),
         x = "效应类型",
         y = "显著特征数量") +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      legend.position = "none"
    ) +
    geom_text(aes(label = n_significant), vjust = -0.5, size = 4)
  
  summary_file <- file.path(output_dir, paste0(output_prefix, "_effect_summary.png"))
  ggsave(summary_file, p_summary, width = 8, height = 6, dpi = 300)
  cat(sprintf("效应统计图已保存: %s\n", summary_file))
}

# 12. 生成可视化
plot_interaction_results(species_results, "宏基因组", "species")
plot_interaction_results(metabolite_results, "代谢组", "metabolite")

# 13. 生成交互效应斜率图
generate_interaction_plots(crc_data, species_results, "宏基因组", "species", top_n = 10)
generate_interaction_plots(crc_data, metabolite_results, "代谢组", "metabolite", top_n = 10)

# 14. 筛选分化相关的显著特征
cat("\n========================================\n")
cat("筛选分化相关的显著特征\n")
cat("========================================\n")

filter_diff_related <- function(results, data_type, output_prefix) {
  cat(sprintf("\n筛选 %s 分化相关特征...\n", data_type))
  
  diff_main <- results[!is.na(results$differentiation_padj) & results$differentiation_padj < 0.05, ]
  diff_main <- diff_main[order(diff_main$differentiation_padj), ]
  cat(sprintf("分化主效应显著特征: %d 个\n", nrow(diff_main)))
  
  interaction_sig <- results[!is.na(results$interaction_padj) & results$interaction_padj < 0.05, ]
  interaction_sig <- interaction_sig[order(interaction_sig$interaction_padj), ]
  cat(sprintf("TNM x 分化交互效应显著特征: %d 个\n", nrow(interaction_sig)))
  
  diff_related <- results[
    (!is.na(results$differentiation_padj) & results$differentiation_padj < 0.05) |
    (!is.na(results$interaction_padj) & results$interaction_padj < 0.05), 
  ]
  diff_related <- diff_related[order(diff_related$differentiation_padj, na.last = FALSE), ]
  cat(sprintf("分化相关特征（主效应或交互效应）: %d 个\n", nrow(diff_related)))
  
  if (nrow(diff_main) > 0) {
    file_main <- file.path(output_dir, paste0(output_prefix, "_diff_main_effect.csv"))
    write.csv(diff_main, file_main, row.names = FALSE)
    cat(sprintf("分化主效应结果已保存: %s\n", file_main))
  }
  
  if (nrow(interaction_sig) > 0) {
    file_inter <- file.path(output_dir, paste0(output_prefix, "_diff_interaction.csv"))
    write.csv(interaction_sig, file_inter, row.names = FALSE)
    cat(sprintf("分化交互效应结果已保存: %s\n", file_inter))
  }
  
  if (nrow(diff_related) > 0) {
    file_all <- file.path(output_dir, paste0(output_prefix, "_diff_related_all.csv"))
    write.csv(diff_related, file_all, row.names = FALSE)
    cat(sprintf("分化相关所有特征已保存: %s\n", file_all))
    
    cat(sprintf("\nTop 10 分化相关 %s 特征:\n", data_type))
    top_diff <- head(diff_related[, c("feature", "differentiation_coef", "differentiation_padj", 
                                           "interaction_coef", "interaction_padj")], 10)
    print(top_diff)
  }
  
  return(list(
    diff_main = diff_main,
    interaction_sig = interaction_sig,
    diff_related = diff_related
  ))
}

species_diff_results <- filter_diff_related(species_results, "宏基因组", "species")
metabolite_diff_results <- filter_diff_related(metabolite_results, "代谢组", "metabolite")

# 15. 生成分化相关特征汇总表
cat("\n生成分化相关特征汇总...\n")

diff_summary <- data.frame(
  data_type = c("宏基因组", "宏基因组", "宏基因组", "代谢组", "代谢组", "代谢组"),
  effect_type = c("分化主效应", "TNM x 分化交互", "分化相关总计", 
                  "分化主效应", "TNM x 分化交互", "分化相关总计"),
  n_features = c(
    nrow(species_diff_results$diff_main),
    nrow(species_diff_results$interaction_sig),
    nrow(species_diff_results$diff_related),
    nrow(metabolite_diff_results$diff_main),
    nrow(metabolite_diff_results$interaction_sig),
    nrow(metabolite_diff_results$diff_related)
  )
)

write.csv(diff_summary, file.path(output_dir, "diff_related_summary.csv"), row.names = FALSE)
cat("分化相关特征汇总已保存\n")
print(diff_summary)

# 16. 保存分析参数
analysis_params <- data.frame(
  parameter = c("分析类型", "样本数", "物种特征数", "代谢物特征数", "显著性阈值"),
  value = c(
    "TNM分期 x 分化 交互效应分析",
    nrow(crc_data),
    length(taxonomy_cols),
    length(metabolite_cols),
    "FDR < 0.05"
  )
)
write.csv(analysis_params, file.path(output_dir, "analysis_parameters.csv"), row.names = FALSE)

cat("\n========================================\n")
cat("交互效应分析完成！\n")
cat(sprintf("结果保存在: %s\n", output_dir))
cat("========================================\n")
