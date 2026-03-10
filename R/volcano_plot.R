library(ape)
library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(data.table)
library(vegan)

theme_microbiome <- function(){ 
  theme_bw(base_size = 14) + 
  theme( 
    panel.grid = element_blank(), 
    panel.border = element_rect(linewidth=0.8), 
    
    axis.text = element_text(color="black"), 
    axis.title = element_text(size=14), 
    
    legend.background = element_blank(), 
    legend.key = element_blank(), 
    
    strip.background = element_blank(), 
    strip.text = element_text(size=13, face="bold") 
  ) 
} 

cat("读取合并后的数据...\n")
merged_data <- fread("dataset/merged_dataset_processed.csv", stringsAsFactors = FALSE, data.table = FALSE)
output_dir <- file.path("results", "R_plots/volcano_plot")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

prepare_diff_groups <- function(df) {
  cat("根据分化程度和吸烟状态创建分组...\n")
  
  df$group <- sapply(1:nrow(df), function(i) {
    crc <- as.integer(df$crc_label[i])
    diff <- as.integer(df$differentiation[i])
    smoking <- as.integer(df$smoking_label[i])
    
    if (!is.na(crc) && crc == 1) {
      if (diff == 1) {
        return('CRC_poor_diff')
      } else if (diff == 0) {
        return('CRC_well_diff')
      } else {
        return(NA)
      }
    } else {
      return('control')
    }
  })
  
  return(df)
}

grouped_data <- prepare_diff_groups(merged_data)

taxonomy_cols <- grep("^tax_", colnames(grouped_data), value = TRUE)
cat(sprintf("找到 %d 个宏基因组特征\n", length(taxonomy_cols)))

taxonomy_data <- grouped_data[, c("SAMPLE_ID", "group", "age", "gender_label", taxonomy_cols)]

species_matrix <- as.matrix(taxonomy_data[, taxonomy_cols])
rownames(species_matrix) <- taxonomy_data$SAMPLE_ID

species_by_sample <- t(species_matrix)
species_names <- rownames(species_by_sample)

group_colors <- c("control" = "#2E86AB", "CRC_well_diff" = "#A23B72", 
                  "CRC_poor_diff" = "#F18F01")

cat("准备分析数据...\n")

abundance_data <- as.data.frame(t(species_by_sample))
colnames(abundance_data) <- species_names
abundance_data$SAMPLE_ID <- rownames(abundance_data)

analysis_data <- merge(taxonomy_data[, c("SAMPLE_ID", "group", "age", "gender_label")], 
                       abundance_data, by = "SAMPLE_ID")

analysis_data$age <- as.numeric(analysis_data$age)
analysis_data$gender <- as.factor(analysis_data$gender_label)

comparison_groups <- list(
  c("control", "CRC_well_diff"),
  c("control", "CRC_poor_diff"),
  c("CRC_well_diff", "CRC_poor_diff")
)

for (i in 1:length(comparison_groups)) {
  group1 <- comparison_groups[[i]][1]
  group2 <- comparison_groups[[i]][2]
  comparison_name <- paste(group1, "vs", group2, sep = "_")
  
  cat(sprintf("\n=== 物种火山图分析组: %s ===\n", comparison_name))
  
  current_data <- analysis_data[analysis_data$group %in% c(group1, group2), ]
  
  if (nrow(current_data) < 5) {
    cat(sprintf("样本数不足，跳过分析: %s\n", comparison_name))
    next
  }
  
  current_data$group_factor <- factor(current_data$group, levels = c(group1, group2))
  
  species_cols <- species_names
  species_matrix <- as.matrix(current_data[, species_cols])
  rownames(species_matrix) <- current_data$SAMPLE_ID
  
  row_sums <- rowSums(species_matrix, na.rm = TRUE)
  empty_rows <- row_sums == 0
  if (any(empty_rows, na.rm = TRUE)) {
    cat(sprintf("移除 %d 个全零样本\n", sum(empty_rows, na.rm = TRUE)))
    species_matrix <- species_matrix[!empty_rows, , drop = FALSE]
    current_data <- current_data[!empty_rows, ]
  }
  
  if (nrow(species_matrix) < 5) {
    cat(sprintf("过滤后样本数不足，跳过分析: %s\n", comparison_name))
    next
  }
  
  col_sums <- colSums(species_matrix, na.rm = TRUE)
  zero_cols <- col_sums == 0
  if (any(zero_cols, na.rm = TRUE)) {
    cat(sprintf("移除 %d 个全零物种\n", sum(zero_cols, na.rm = TRUE)))
    species_matrix <- species_matrix[, !zero_cols, drop = FALSE]
  }
  
  if (any(is.na(species_matrix))) {
    cat("警告: 物种矩阵中存在NA值，将使用0替换\n")
    species_matrix[is.na(species_matrix)] <- 0
  }
  
  cat("对每个物种运行带协变量的线性模型...\n")
  filtered_species_cols <- colnames(species_matrix)
  species_results <- data.frame(
    species = filtered_species_cols,
    log2FC = NA,
    pvalue = NA,
    stringsAsFactors = FALSE
  )
  
  for (j in 1:length(filtered_species_cols)) {
    species_col <- filtered_species_cols[j]
    y <- current_data[[species_col]]
    
    if (sum(y > 0) < 3) next
    
    fit <- try(lm(log2(y + 1e-5) ~ group_factor + age + gender, data = current_data), silent = TRUE)
    
    if (!inherits(fit, "try-error")) {
      coef_sum <- summary(fit)$coefficients
      if (paste0("group_factor", group2) %in% rownames(coef_sum)) {
        species_results$log2FC[j] <- coef_sum[paste0("group_factor", group2), "Estimate"]
        species_results$pvalue[j] <- coef_sum[paste0("group_factor", group2), "Pr(>|t|)"]
      }
    }
  }
  
  species_results$p_adj <- p.adjust(species_results$pvalue, method = "BH")
  species_results$significant <- species_results$p_adj < 0.05
  
  significant_species <- species_results[species_results$significant, ]
  significant_species <- significant_species[order(abs(significant_species$log2FC), decreasing = TRUE), ]
  
  cat(sprintf("找到 %d 个差异显著的物种\n", nrow(significant_species)))
  
  cat("生成火山图...\n")
  
  volcano_data <- species_results
  volcano_data$log10p <- -log10(volcano_data$p_adj)
  volcano_data$significance <- "Not Significant"
  volcano_data$significance[volcano_data$p_adj < 0.05 & volcano_data$log2FC > 0] <- "Up"
  volcano_data$significance[volcano_data$p_adj < 0.05 & volcano_data$log2FC < 0] <- "Down"
  
  volcano_colors <- c("Up" = "#E74C3C", "Down" = "#3498DB", "Not Significant" = "#95A5A6")
  
  volcano_plot <- ggplot(volcano_data, aes(x = log2FC, y = log10p, color = significance)) +
    geom_point(size = 2, alpha = 0.6) +
    scale_color_manual(values = volcano_colors) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50") +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray50") +
    labs(title = paste("Volcano Plot -", comparison_name),
         subtitle = sprintf("Significant species: Up = %d, Down = %d", 
                           sum(volcano_data$significance == "Up"),
                           sum(volcano_data$significance == "Down")),
         x = "log2 Fold Change",
         y = "-log10(q-value)",
         color = "Significance") +
    theme_microbiome() +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 10),
          legend.position = "right")
  
  if (nrow(significant_species) > 0) {
    top_species <- head(significant_species, 10)
    volcano_plot <- volcano_plot +
      geom_text(data = volcano_data[volcano_data$species %in% top_species$species, ],
                aes(label = species), size = 3, vjust = -0.5, hjust = 0.5, color = "black")
  }
  
  ggsave(file.path(output_dir, paste0("Volcano_species_", comparison_name, ".png")), 
         volcano_plot, width = 8, height = 6, dpi = 300)
  cat(sprintf("火山图已保存: Volcano_species_%s.png\n", comparison_name))
  
  if (nrow(significant_species) > 0) {
    cat("前10个差异最大的物种:\n")
    print(head(significant_species, 10))
    
    write.csv(significant_species, 
              file.path(output_dir, paste0("significant_species_", comparison_name, ".csv")))
  } else {
    cat("没有找到差异显著的物种\n")
  }
}

write.csv(data.frame(species_name = species_names), 
          file.path(output_dir, "species_names.csv"), row.names = FALSE)

cat("\n========== 开始代谢物火山图分析 ==========\n")

cat("读取代谢物数据...\n")
metabolome_data <- fread("dataset/metabolome_data.csv", stringsAsFactors = FALSE, data.table = FALSE)

metabolite_names <- metabolome_data[, 1]
metabolome_matrix <- as.matrix(metabolome_data[, -1])
rownames(metabolome_matrix) <- metabolite_names

cat(sprintf("代谢物数据: %d 个代谢物, %d 个样本\n", nrow(metabolome_matrix), ncol(metabolome_matrix)))

cat("处理代谢物数据...\n")
metabolome_by_sample <- t(metabolome_matrix)
metabolome_by_sample <- metabolome_by_sample[, colSums(metabolome_by_sample, na.rm = TRUE) > 0, drop = FALSE]
cat(sprintf("移除全零代谢物后: %d个代谢物\n", ncol(metabolome_by_sample)))

mean_metab_abundance <- colMeans(metabolome_by_sample, na.rm = TRUE)
metab_abundance_threshold <- 100
metabolome_by_sample <- metabolome_by_sample[, mean_metab_abundance >= metab_abundance_threshold, drop = FALSE]
cat(sprintf("过滤掉平均丰度 < %.2f 的代谢物后: %d个代谢物\n", metab_abundance_threshold, ncol(metabolome_by_sample)))

final_metabolite_names <- colnames(metabolome_by_sample)

metabolome_df <- as.data.frame(metabolome_by_sample)
metabolome_df$SAMPLE_ID <- rownames(metabolome_by_sample)

metab_analysis_data <- merge(
  taxonomy_data[, c("SAMPLE_ID", "group", "age", "gender_label")],
  metabolome_df,
  by = "SAMPLE_ID"
)

metab_analysis_data$age <- as.numeric(metab_analysis_data$age)
metab_analysis_data$gender <- as.factor(metab_analysis_data$gender_label)

cat(sprintf("合并后代谢物分析数据: %d 个样本\n", nrow(metab_analysis_data)))

for (i in 1:length(comparison_groups)) {
  group1 <- comparison_groups[[i]][1]
  group2 <- comparison_groups[[i]][2]
  comparison_name <- paste(group1, "vs", group2, sep = "_")
  
  cat(sprintf("\n=== 代谢物火山图分析组: %s ===\n", comparison_name))
  
  current_metab_data <- metab_analysis_data[metab_analysis_data$group %in% c(group1, group2), ]
  
  if (nrow(current_metab_data) < 5) {
    cat(sprintf("样本数不足，跳过代谢物分析: %s\n", comparison_name))
    next
  }
  
  current_metab_data$group_factor <- factor(current_metab_data$group, levels = c(group1, group2))
  
  metab_cols <- final_metabolite_names
  metab_matrix <- as.matrix(current_metab_data[, metab_cols])
  rownames(metab_matrix) <- current_metab_data$SAMPLE_ID
  
  row_sums_metab <- rowSums(metab_matrix, na.rm = TRUE)
  empty_rows_metab <- row_sums_metab == 0
  if (any(empty_rows_metab, na.rm = TRUE)) {
    cat(sprintf("移除 %d 个全零样本\n", sum(empty_rows_metab, na.rm = TRUE)))
    metab_matrix <- metab_matrix[!empty_rows_metab, , drop = FALSE]
    current_metab_data <- current_metab_data[!empty_rows_metab, ]
  }
  
  if (nrow(metab_matrix) < 5) {
    cat(sprintf("过滤后样本数不足，跳过代谢物分析: %s\n", comparison_name))
    next
  }
  
  col_sums_metab <- colSums(metab_matrix, na.rm = TRUE)
  zero_cols_metab <- col_sums_metab == 0
  if (any(zero_cols_metab, na.rm = TRUE)) {
    cat(sprintf("移除 %d 个全零代谢物\n", sum(zero_cols_metab, na.rm = TRUE)))
    metab_matrix <- metab_matrix[, !zero_cols_metab, drop = FALSE]
  }
  
  if (any(is.na(metab_matrix))) {
    cat("警告: 代谢物矩阵中存在NA值，将使用中位数填充\n")
    for (col in colnames(metab_matrix)) {
      if (any(is.na(metab_matrix[, col]))) {
        metab_matrix[is.na(metab_matrix[, col]), col] <- median(metab_matrix[, col], na.rm = TRUE)
      }
    }
  }
  
  cat("对每个代谢物运行带协变量的线性模型...\n")
  filtered_metab_cols <- colnames(metab_matrix)
  metab_results <- data.frame(
    metabolite = filtered_metab_cols,
    log2FC = NA,
    pvalue = NA,
    stringsAsFactors = FALSE
  )
  
  for (j in 1:length(filtered_metab_cols)) {
    metab_col <- filtered_metab_cols[j]
    y <- current_metab_data[[metab_col]]
    
    if (sum(y > 0, na.rm = TRUE) < 3) next
    
    fit <- try(lm(log2(y + 1) ~ group_factor + age + gender, data = current_metab_data), silent = TRUE)
    
    if (!inherits(fit, "try-error")) {
      coef_sum <- summary(fit)$coefficients
      if (paste0("group_factor", group2) %in% rownames(coef_sum)) {
        metab_results$log2FC[j] <- coef_sum[paste0("group_factor", group2), "Estimate"]
        metab_results$pvalue[j] <- coef_sum[paste0("group_factor", group2), "Pr(>|t|)"]
      }
    }
  }
  
  metab_results$p_adj <- p.adjust(metab_results$pvalue, method = "BH")
  metab_results$significant <- metab_results$p_adj < 0.05
  
  significant_metabolites <- metab_results[metab_results$significant, ]
  significant_metabolites <- significant_metabolites[order(abs(significant_metabolites$log2FC), decreasing = TRUE), ]
  
  cat(sprintf("找到 %d 个差异显著的代谢物\n", nrow(significant_metabolites)))
  
  cat("生成代谢物火山图...\n")
  
  metab_volcano_data <- metab_results
  metab_volcano_data$log10p <- -log10(metab_volcano_data$pvalue)
  metab_volcano_data$significance <- "Not Significant"
  metab_volcano_data$significance[metab_volcano_data$p_adj < 0.05 & metab_volcano_data$log2FC > 0] <- "Up"
  metab_volcano_data$significance[metab_volcano_data$p_adj < 0.05 & metab_volcano_data$log2FC < 0] <- "Down"
  
  metab_volcano_colors <- c("Up" = "#E74C3C", "Down" = "#3498DB", "Not Significant" = "#95A5A6")
  
  metab_volcano_plot <- ggplot(metab_volcano_data, aes(x = log2FC, y = log10p, color = significance)) +
    geom_point(size = 2, alpha = 0.6) +
    scale_color_manual(values = metab_volcano_colors) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50") +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray50") +
    labs(title = paste("Volcano Plot (Metabolites) -", comparison_name),
         subtitle = sprintf("Significant metabolites: Up = %d, Down = %d", 
                           sum(metab_volcano_data$significance == "Up"),
                           sum(metab_volcano_data$significance == "Down")),
         x = "log2 Fold Change",
         y = "-log10(p-value)",
         color = "Significance") +
    theme_microbiome() +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 10),
          legend.position = "right")
  
  if (nrow(significant_metabolites) > 0) {
    top_metabolites <- head(significant_metabolites, 10)
    metab_volcano_plot <- metab_volcano_plot +
      geom_text(data = metab_volcano_data[metab_volcano_data$metabolite %in% top_metabolites$metabolite, ],
                aes(label = metabolite), size = 3, vjust = -0.5, hjust = 0.5, color = "black")
  }
  
  ggsave(file.path(output_dir, paste0("Volcano_metabolites_", comparison_name, ".png")), 
         metab_volcano_plot, width = 8, height = 6, dpi = 300)
  cat(sprintf("代谢物火山图已保存: Volcano_metabolites_%s.png\n", comparison_name))
  
  if (nrow(significant_metabolites) > 0) {
    cat("前10个差异最大的代谢物:\n")
    print(head(significant_metabolites, 10))
    
    write.csv(significant_metabolites, 
              file.path(output_dir, paste0("significant_metabolites_", comparison_name, ".csv")))
  } else {
    cat("没有找到差异显著的代谢物\n")
  }
  
  write.csv(metab_results, 
            file.path(output_dir, paste0("all_metabolites_results_", comparison_name, ".csv")))
}

write.csv(data.frame(metabolite_name = final_metabolite_names), 
          file.path(output_dir, "metabolite_names.csv"), row.names = FALSE)

cat("\n========== 火山图分析完成 ==========\n")
cat("\n分析完成！结果保存在", output_dir, "目录中\n")
