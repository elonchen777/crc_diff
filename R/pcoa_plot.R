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
output_dir <- file.path("results", "R_plots/pcoa_plot")
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
  
  cat(sprintf("\n=== 物种PCoA分析组: %s ===\n", comparison_name))
  
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
  
  cat("运行adonis多变量回归分析...\n")
  adonis_result <- adonis2(
    species_matrix ~ group_factor + age + gender,
    data = current_data,
    permutations = 999,
    method = "bray"
  )
  
  cat("Adonis分析结果:\n")
  print(adonis_result)
  
  write.csv(adonis_result, file.path(output_dir, paste0("adonis_", comparison_name, ".csv")))
  
  cat("生成PCoA图...\n")
  
  dist_matrix <- vegdist(species_matrix, method = "bray")
  
  pcoa_result <- cmdscale(dist_matrix, k = 2, eig = TRUE)
  pcoa_points <- as.data.frame(pcoa_result$points)
  colnames(pcoa_points) <- c("PCo1", "PCo2")
  
  eig_values <- pcoa_result$eig
  var_explained <- eig_values / sum(eig_values) * 100
  pc1_var <- round(var_explained[1], 1)
  pc2_var <- round(var_explained[2], 1)
  
  pcoa_points$group_factor <- current_data$group_factor[match(rownames(pcoa_points), current_data$SAMPLE_ID)]
  
  adonis_pvalue <- adonis_result$`Pr(>F)`[1]
  adonis_r2 <- adonis_result$R2[1]
  
  active_colors <- group_colors[c(group1, group2)]
  
  cat(sprintf("PCoA数据点数: %d, 分组数: %s\n", nrow(pcoa_points), 
              paste(unique(pcoa_points$group_factor), collapse = ", ")))
  
  pcoa_plot <- ggplot(pcoa_points, aes(x = PCo1, y = PCo2, color = group_factor)) +
    geom_point(size = 3, alpha = 0.7) +
    stat_ellipse(level = 0.95, linetype = 1, size = 1) +
    scale_color_manual(values = active_colors, drop = FALSE) +
    labs(title = paste("PCoA -", comparison_name),
         subtitle = sprintf("PERMANOVA: R² = %.3f, p = %.4f", adonis_r2, adonis_pvalue),
         x = paste0("PCo1 (", pc1_var, "%)"),
         y = paste0("PCo2 (", pc2_var, "%)"),
         color = "Group") +
    theme_microbiome() +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 10),
          legend.position = "right")
  
  ggsave(file.path(output_dir, paste0("PCoA_species_", comparison_name, ".png")), 
         pcoa_plot, width = 8, height = 6, dpi = 300)
  cat(sprintf("PCoA图已保存: PCoA_species_%s.png\n", comparison_name))
}

write.csv(data.frame(species_name = species_names), 
          file.path(output_dir, "species_names.csv"), row.names = FALSE)

cat("\n========== 开始代谢物PCoA分析 ==========\n")

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
  
  cat(sprintf("\n=== 代谢物PCoA分析组: %s ===\n", comparison_name))
  
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
  
  cat("运行代谢物adonis多变量回归分析...\n")
  adonis_metab_result <- adonis2(
    metab_matrix ~ group_factor + age + gender,
    data = current_metab_data,
    permutations = 999,
    method = "euclidean"
  )
  
  cat("代谢物Adonis分析结果:\n")
  print(adonis_metab_result)
  
  write.csv(adonis_metab_result, file.path(output_dir, paste0("adonis_metabolites_", comparison_name, ".csv")))
  
  cat("生成代谢物PCoA图...\n")
  
  metab_dist_matrix <- dist(metab_matrix, method = "euclidean")
  
  metab_pcoa_result <- cmdscale(metab_dist_matrix, k = 2, eig = TRUE)
  metab_pcoa_points <- as.data.frame(metab_pcoa_result$points)
  colnames(metab_pcoa_points) <- c("PCo1", "PCo2")
  
  metab_eig_values <- metab_pcoa_result$eig
  metab_var_explained <- metab_eig_values / sum(metab_eig_values) * 100
  metab_pc1_var <- round(metab_var_explained[1], 1)
  metab_pc2_var <- round(metab_var_explained[2], 1)
  
  cat(sprintf("代谢物PCoA行名示例: %s\n", paste(head(rownames(metab_pcoa_points)), collapse = ", ")))
  cat(sprintf("current_metab_data SAMPLE_ID示例: %s\n", paste(head(current_metab_data$SAMPLE_ID), collapse = ", ")))
  
  metab_pcoa_points$group_factor <- current_metab_data$group_factor[match(rownames(metab_pcoa_points), current_metab_data$SAMPLE_ID)]
  
  adonis_metab_pvalue <- adonis_metab_result$`Pr(>F)`[1]
  adonis_metab_r2 <- adonis_metab_result$R2[1]
  
  metab_active_colors <- group_colors[c(group1, group2)]
  
  cat(sprintf("代谢物PCoA数据点数: %d, 分组数: %s\n", nrow(metab_pcoa_points), 
              paste(unique(metab_pcoa_points$group_factor), collapse = ", ")))
  
  metab_pcoa_plot <- ggplot(metab_pcoa_points, aes(x = PCo1, y = PCo2, color = group_factor)) +
    geom_point(size = 3, alpha = 0.7) +
    stat_ellipse(level = 0.95, linetype = 1, size = 1) +
    scale_color_manual(values = metab_active_colors, drop = FALSE) +
    labs(title = paste("PCoA (Metabolites) -", comparison_name),
         subtitle = sprintf("PERMANOVA: R² = %.3f, p = %.4f", adonis_metab_r2, adonis_metab_pvalue),
         x = paste0("PCo1 (", metab_pc1_var, "%)"),
         y = paste0("PCo2 (", metab_pc2_var, "%)"),
         color = "Group") +
    theme_microbiome() +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 10),
          legend.position = "right")
  
  ggsave(file.path(output_dir, paste0("PCoA_metabolites_", comparison_name, ".png")), 
         metab_pcoa_plot, width = 8, height = 6, dpi = 300)
  cat(sprintf("代谢物PCoA图已保存: PCoA_metabolites_%s.png\n", comparison_name))
}

write.csv(data.frame(metabolite_name = final_metabolite_names), 
          file.path(output_dir, "metabolite_names.csv"), row.names = FALSE)

cat("\n========== PCoA分析完成 ==========\n")
cat("\n分析完成！结果保存在", output_dir, "目录中\n")
