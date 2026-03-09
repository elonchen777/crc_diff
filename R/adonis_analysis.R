# 加载必要的库
library(ape)
library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(data.table)
library(vegan)  # 用于adonis多变量回归分析
library(ropls)  # 用于OPLS-DA分析

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

# 1. 读取数据
cat("读取合并后的数据...\n")
merged_data <- fread("dataset/merged_dataset_processed.csv", stringsAsFactors = FALSE, data.table = FALSE)
output_dir <- file.path("results", "R_plots")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# 2. 准备分组函数
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

# 3. 应用分组函数
grouped_data <- prepare_diff_groups(merged_data)

# 4. 提取宏基因组数据
taxonomy_cols <- grep("^tax_", colnames(grouped_data), value = TRUE)
cat(sprintf("找到 %d 个宏基因组特征\n", length(taxonomy_cols)))

taxonomy_data <- grouped_data[, c("SAMPLE_ID", "group", "age", "gender_label", taxonomy_cols)]

# 5. 提取物种丰度矩阵
species_matrix <- as.matrix(taxonomy_data[, taxonomy_cols])
rownames(species_matrix) <- taxonomy_data$SAMPLE_ID

species_by_sample <- t(species_matrix)
species_names <- rownames(species_by_sample)

# 定义分组颜色映射（全局变量，供后续使用）
group_colors <- c("control" = "#2E86AB", "CRC_well_diff" = "#A23B72", 
                  "CRC_poor_diff" = "#F18F01")

# 6. 准备分析数据
cat("准备分析数据...\n")

# 转置丰度数据为样本x物种
abundance_data <- as.data.frame(t(species_by_sample))
colnames(abundance_data) <- species_names
abundance_data$SAMPLE_ID <- rownames(abundance_data)

# 合并样本信息
analysis_data <- merge(taxonomy_data[, c("SAMPLE_ID", "group", "age", "gender_label")], 
                       abundance_data, by = "SAMPLE_ID")

# 确保年龄和性别是数值型
analysis_data$age <- as.numeric(analysis_data$age)
analysis_data$gender <- as.factor(analysis_data$gender_label)

# 7. 定义比较组
comparison_groups <- list(
  c("control", "CRC_well_diff"),
  c("control", "CRC_poor_diff"),
  c("CRC_well_diff", "CRC_poor_diff")
)

# 8. 对每个比较组进行分析
for (i in 1:length(comparison_groups)) {
  group1 <- comparison_groups[[i]][1]
  group2 <- comparison_groups[[i]][2]
  comparison_name <- paste(group1, "vs", group2, sep = "_")
  
  cat(sprintf("\n=== 分析组: %s ===\n", comparison_name))
  
  # 筛选当前比较的两组样本
  current_data <- analysis_data[analysis_data$group %in% c(group1, group2), ]
  
  if (nrow(current_data) < 5) {
    cat(sprintf("样本数不足，跳过分析: %s\n", comparison_name))
    next
  }
  
  # 创建分组因子
  current_data$group_factor <- factor(current_data$group, levels = c(group1, group2))
  
  # 提取物种丰度数据
  species_cols <- species_names
  species_matrix <- as.matrix(current_data[, species_cols])
  # 确保行名是SAMPLE_ID
  rownames(species_matrix) <- current_data$SAMPLE_ID
  
  # 检查并处理缺失值和空行
  # 1. 移除全零行（空样本）
  row_sums <- rowSums(species_matrix, na.rm = TRUE)
  empty_rows <- row_sums == 0
  if (any(empty_rows, na.rm = TRUE)) {
    cat(sprintf("移除 %d 个全零样本\n", sum(empty_rows, na.rm = TRUE)))
    species_matrix <- species_matrix[!empty_rows, , drop = FALSE]
    current_data <- current_data[!empty_rows, ]
  }
  
  # 2. 检查是否有足够的样本继续分析
  if (nrow(species_matrix) < 5) {
    cat(sprintf("过滤后样本数不足，跳过分析: %s\n", comparison_name))
    next
  }
  
  # 3. 将剩余的全零列（物种）移除
  col_sums <- colSums(species_matrix, na.rm = TRUE)
  zero_cols <- col_sums == 0
  if (any(zero_cols, na.rm = TRUE)) {
    cat(sprintf("移除 %d 个全零物种\n", sum(zero_cols, na.rm = TRUE)))
    species_matrix <- species_matrix[, !zero_cols, drop = FALSE]
  }
  
  # 4. 检查是否有NA值
  if (any(is.na(species_matrix))) {
    cat("警告: 物种矩阵中存在NA值，将使用0替换\n")
    species_matrix[is.na(species_matrix)] <- 0
  }
  
  # 运行adonis多变量回归分析
  cat("运行adonis多变量回归分析...\n")
  adonis_result <- adonis2(
    species_matrix ~ group_factor + age + gender,
    data = current_data,
    permutations = 999,
    method = "bray"
  )
  
  cat("Adonis分析结果:\n")
  print(adonis_result)
  
  # 保存adonis结果
  write.csv(adonis_result, file.path(output_dir, paste0("adonis_", comparison_name, ".csv")))
  
  # ==================== PCoA可视化 ====================
  cat("生成PCoA图...\n")
  
  # 计算Bray-Curtis距离矩阵
  dist_matrix <- vegdist(species_matrix, method = "bray")
  
  # 进行PCoA分析
  pcoa_result <- cmdscale(dist_matrix, k = 2, eig = TRUE)
  pcoa_points <- as.data.frame(pcoa_result$points)
  colnames(pcoa_points) <- c("PCo1", "PCo2")
  
  # 计算解释方差比例
  eig_values <- pcoa_result$eig
  var_explained <- eig_values / sum(eig_values) * 100
  pc1_var <- round(var_explained[1], 1)
  pc2_var <- round(var_explained[2], 1)
  
  # 直接添加分组信息（使用索引匹配）
  pcoa_points$group_factor <- current_data$group_factor[match(rownames(pcoa_points), current_data$SAMPLE_ID)]
  
  # 获取adonis的p值
  adonis_pvalue <- adonis_result$`Pr(>F)`[1]
  adonis_r2 <- adonis_result$R2[1]
  
  # 绘制PCoA图
  # 直接使用当前比较的分组设置颜色
  active_colors <- group_colors[c(group1, group2)]
  
  cat(sprintf("PCoA数据点数: %d, 分组数: %s\n", nrow(pcoa_points), 
              paste(unique(pcoa_points$group_factor), collapse = ", ")))
  
  pcoa_plot <- ggplot(pcoa_points, aes(x = PCo1, y = PCo2, color = group_factor)) +
    geom_point(size = 3, alpha = 0.7) +
    stat_ellipse(level = 0.95, linetype = "dashed", size = 1) +
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
  
  # 保存PCoA图
  ggsave(file.path(output_dir, paste0("PCoA_species_", comparison_name, ".png")), 
         pcoa_plot, width = 8, height = 6, dpi = 300)
  cat(sprintf("PCoA图已保存: PCoA_species_%s.png\n", comparison_name))
  
  # 对每个物种运行带协变量的线性模型
  cat("对每个物种运行带协变量的线性模型...\n")
  # 使用过滤后的物种列名
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
    
    # 跳过全零或几乎全零的物种
    if (sum(y > 0) < 3) next
    
    # 拟合带协变量的线性模型
    fit <- try(lm(log2(y + 1e-5) ~ group_factor + age + gender, data = current_data), silent = TRUE)
    
    if (!inherits(fit, "try-error")) {
      coef_sum <- summary(fit)$coefficients
      if (paste0("group_factor", group2) %in% rownames(coef_sum)) {
        species_results$log2FC[j] <- coef_sum[paste0("group_factor", group2), "Estimate"]
        species_results$pvalue[j] <- coef_sum[paste0("group_factor", group2), "Pr(>|t|)"]
      }
    }
  }
  
  # 应用Benjamini-Hochberg校正
  species_results$p_adj <- p.adjust(species_results$pvalue, method = "BH")
  species_results$significant <- species_results$p_adj < 0.05
  
  # 筛选差异显著的物种并按log2FC绝对值排序
  significant_species <- species_results[species_results$significant, ]
  significant_species <- significant_species[order(abs(significant_species$log2FC), decreasing = TRUE), ]
  
  cat(sprintf("找到 %d 个差异显著的物种\n", nrow(significant_species)))
  
  # ==================== 火山图可视化 ====================
  cat("生成火山图...\n")
  
  # 准备火山图数据
  volcano_data <- species_results
  volcano_data$log10p <- -log10(volcano_data$p_adj)
  volcano_data$significance <- "Not Significant"
  volcano_data$significance[volcano_data$p_adj < 0.05 & volcano_data$log2FC > 0] <- "Up"
  volcano_data$significance[volcano_data$p_adj < 0.05 & volcano_data$log2FC < 0] <- "Down"
  
  # 设置颜色
  volcano_colors <- c("Up" = "#E74C3C", "Down" = "#3498DB", "Not Significant" = "#95A5A6")
  
  # 绘制火山图
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
  
  # 添加显著物种的标签（前10个）
  if (nrow(significant_species) > 0) {
    top_species <- head(significant_species, 10)
    volcano_plot <- volcano_plot +
      geom_text(data = volcano_data[volcano_data$species %in% top_species$species, ],
                aes(label = species), size = 3, vjust = -0.5, hjust = 0.5, color = "black")
  }
  
  # 保存火山图
  ggsave(file.path(output_dir, paste0("Volcano_species_", comparison_name, ".png")), 
         volcano_plot, width = 8, height = 6, dpi = 300)
  cat(sprintf("火山图已保存: Volcano_species_%s.png\n", comparison_name))
  
  if (nrow(significant_species) > 0) {
    cat("前10个差异最大的物种:\n")
    print(head(significant_species, 10))
    
    # 保存结果
    write.csv(significant_species, 
              file.path(output_dir, paste0("significant_species_", comparison_name, ".csv")))
  } else {
    cat("没有找到差异显著的物种\n")
  }
}

# 8. 保存总体结果
cat("\n保存总体分析结果...\n")

# 保存物种名称列表
write.csv(data.frame(species_name = species_names), 
          file.path(output_dir, "species_names.csv"), row.names = FALSE)

# ==================== 代谢物差异分析 ====================
cat("\n========== 开始代谢物差异分析 ==========\n")

# 1. 读取代谢物数据
cat("读取代谢物数据...\n")
metabolome_data <- fread("dataset/metabolome_data.csv", stringsAsFactors = FALSE, data.table = FALSE)

# 转置代谢物数据（行为代谢物，列为样本）
metabolite_names <- metabolome_data[, 1]
metabolome_matrix <- as.matrix(metabolome_data[, -1])
rownames(metabolome_matrix) <- metabolite_names


cat(sprintf("代谢物数据: %d 个代谢物, %d 个样本\n", nrow(metabolome_matrix), ncol(metabolome_matrix)))

# 2. 处理代谢物数据
cat("处理代谢物数据...\n")
# 转置为样本x代谢物
metabolome_by_sample <- t(metabolome_matrix)
# 移除全零代谢物
metabolome_by_sample <- metabolome_by_sample[, colSums(metabolome_by_sample, na.rm = TRUE) > 0, drop = FALSE]
cat(sprintf("移除全零代谢物后: %d个代谢物\n", ncol(metabolome_by_sample)))

# 过滤低丰度代谢物（平均丰度 < 100）
mean_metab_abundance <- colMeans(metabolome_by_sample, na.rm = TRUE)
metab_abundance_threshold <- 100
metabolome_by_sample <- metabolome_by_sample[, mean_metab_abundance >= metab_abundance_threshold, drop = FALSE]
cat(sprintf("过滤掉平均丰度 < %.2f 的代谢物后: %d个代谢物\n", metab_abundance_threshold, ncol(metabolome_by_sample)))

final_metabolite_names <- colnames(metabolome_by_sample)

# 3. 准备代谢物分析数据
metabolome_df <- as.data.frame(metabolome_by_sample)
metabolome_df$SAMPLE_ID <- rownames(metabolome_by_sample)

# 合并样本信息（使用与物种分析相同的分组数据）
metab_analysis_data <- merge(
  taxonomy_data[, c("SAMPLE_ID", "group", "age", "gender_label")],
  metabolome_df,
  by = "SAMPLE_ID"
)

# 确保年龄和性别是数值型
metab_analysis_data$age <- as.numeric(metab_analysis_data$age)
metab_analysis_data$gender <- as.factor(metab_analysis_data$gender_label)

cat(sprintf("合并后代谢物分析数据: %d 个样本\n", nrow(metab_analysis_data)))

# 4. 对代谢物进行相同的比较组分析
for (i in 1:length(comparison_groups)) {
  group1 <- comparison_groups[[i]][1]
  group2 <- comparison_groups[[i]][2]
  comparison_name <- paste(group1, "vs", group2, sep = "_")
  
  cat(sprintf("\n=== 代谢物分析组: %s ===\n", comparison_name))
  
  # 筛选当前比较的两组样本
  current_metab_data <- metab_analysis_data[metab_analysis_data$group %in% c(group1, group2), ]
  
  if (nrow(current_metab_data) < 5) {
    cat(sprintf("样本数不足，跳过代谢物分析: %s\n", comparison_name))
    next
  }
  
  # 创建分组因子
  current_metab_data$group_factor <- factor(current_metab_data$group, levels = c(group1, group2))
  
  # 提取代谢物丰度数据
  metab_cols <- final_metabolite_names
  metab_matrix <- as.matrix(current_metab_data[, metab_cols])
  # 确保行名是SAMPLE_ID
  rownames(metab_matrix) <- current_metab_data$SAMPLE_ID
  
  # 检查并处理缺失值和空行
  # 1. 移除全零行（空样本）
  row_sums_metab <- rowSums(metab_matrix, na.rm = TRUE)
  empty_rows_metab <- row_sums_metab == 0
  if (any(empty_rows_metab, na.rm = TRUE)) {
    cat(sprintf("移除 %d 个全零样本\n", sum(empty_rows_metab, na.rm = TRUE)))
    metab_matrix <- metab_matrix[!empty_rows_metab, , drop = FALSE]
    current_metab_data <- current_metab_data[!empty_rows_metab, ]
  }
  
  # 2. 检查是否有足够的样本继续分析
  if (nrow(metab_matrix) < 5) {
    cat(sprintf("过滤后样本数不足，跳过代谢物分析: %s\n", comparison_name))
    next
  }
  
  # 3. 将剩余的全零列（代谢物）移除
  col_sums_metab <- colSums(metab_matrix, na.rm = TRUE)
  zero_cols_metab <- col_sums_metab == 0
  if (any(zero_cols_metab, na.rm = TRUE)) {
    cat(sprintf("移除 %d 个全零代谢物\n", sum(zero_cols_metab, na.rm = TRUE)))
    metab_matrix <- metab_matrix[, !zero_cols_metab, drop = FALSE]
  }
  
  # 4. 检查是否有NA值
  if (any(is.na(metab_matrix))) {
    cat("警告: 代谢物矩阵中存在NA值，将使用中位数填充\n")
    for (col in colnames(metab_matrix)) {
      if (any(is.na(metab_matrix[, col]))) {
        metab_matrix[is.na(metab_matrix[, col]), col] <- median(metab_matrix[, col], na.rm = TRUE)
      }
    }
  }
  
  # 运行adonis多变量回归分析
  cat("运行代谢物adonis多变量回归分析...\n")
  adonis_metab_result <- adonis2(
    metab_matrix ~ group_factor + age + gender,
    data = current_metab_data,
    permutations = 999,
    method = "euclidean"
  )
  
  cat("代谢物Adonis分析结果:\n")
  print(adonis_metab_result)
  
  # 保存adonis结果
  write.csv(adonis_metab_result, file.path(output_dir, paste0("adonis_metabolites_", comparison_name, ".csv")))
  
  # ==================== 代谢物PCoA可视化 ====================
  cat("生成代谢物PCoA图...\n")
  
  # 计算欧氏距离矩阵
  metab_dist_matrix <- dist(metab_matrix, method = "euclidean")
  
  # 进行PCoA分析
  metab_pcoa_result <- cmdscale(metab_dist_matrix, k = 2, eig = TRUE)
  metab_pcoa_points <- as.data.frame(metab_pcoa_result$points)
  colnames(metab_pcoa_points) <- c("PCo1", "PCo2")
  
  # 计算解释方差比例
  metab_eig_values <- metab_pcoa_result$eig
  metab_var_explained <- metab_eig_values / sum(metab_eig_values) * 100
  metab_pc1_var <- round(metab_var_explained[1], 1)
  metab_pc2_var <- round(metab_var_explained[2], 1)
  
  # 直接添加分组信息（使用索引匹配）
  # 调试：检查行名和SAMPLE_ID
  cat(sprintf("代谢物PCoA行名示例: %s\n", paste(head(rownames(metab_pcoa_points)), collapse = ", ")))
  cat(sprintf("current_metab_data SAMPLE_ID示例: %s\n", paste(head(current_metab_data$SAMPLE_ID), collapse = ", ")))
  
  metab_pcoa_points$group_factor <- current_metab_data$group_factor[match(rownames(metab_pcoa_points), current_metab_data$SAMPLE_ID)]
  
  # 获取adonis的p值
  adonis_metab_pvalue <- adonis_metab_result$`Pr(>F)`[1]
  adonis_metab_r2 <- adonis_metab_result$R2[1]
  
  # 绘制代谢物PCoA图
  # 直接使用当前比较的分组设置颜色
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
  
  # 保存代谢物PCoA图
  ggsave(file.path(output_dir, paste0("PCoA_metabolites_", comparison_name, ".png")), 
         metab_pcoa_plot, width = 8, height = 6, dpi = 300)
  cat(sprintf("代谢物PCoA图已保存: PCoA_metabolites_%s.png\n", comparison_name))
  
  # 对每个代谢物运行带协变量的线性模型
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
    
    # 跳过全零或几乎全零的代谢物
    if (sum(y > 0, na.rm = TRUE) < 3) next
    
    # 拟合带协变量的线性模型（使用log2转换，添加1避免log(0)）
    fit <- try(lm(log2(y + 1) ~ group_factor + age + gender, data = current_metab_data), silent = TRUE)
    
    if (!inherits(fit, "try-error")) {
      coef_sum <- summary(fit)$coefficients
      if (paste0("group_factor", group2) %in% rownames(coef_sum)) {
        metab_results$log2FC[j] <- coef_sum[paste0("group_factor", group2), "Estimate"]
        metab_results$pvalue[j] <- coef_sum[paste0("group_factor", group2), "Pr(>|t|)"]
      }
    }
  }
  
  # 应用Benjamini-Hochberg校正
  metab_results$p_adj <- p.adjust(metab_results$pvalue, method = "BH")
  metab_results$significant <- metab_results$p_adj < 0.05
  
  # 筛选差异显著的代谢物并按log2FC绝对值排序
  significant_metabolites <- metab_results[metab_results$significant, ]
  significant_metabolites <- significant_metabolites[order(abs(significant_metabolites$log2FC), decreasing = TRUE), ]
  
  cat(sprintf("找到 %d 个差异显著的代谢物\n", nrow(significant_metabolites)))
  
  # ==================== 代谢物火山图可视化 ====================
  cat("生成代谢物火山图...\n")
  
  # 准备火山图数据
  metab_volcano_data <- metab_results
  metab_volcano_data$log10p <- -log10(metab_volcano_data$pvalue)
  metab_volcano_data$significance <- "Not Significant"
  metab_volcano_data$significance[metab_volcano_data$p_adj < 0.05 & metab_volcano_data$log2FC > 0] <- "Up"
  metab_volcano_data$significance[metab_volcano_data$p_adj < 0.05 & metab_volcano_data$log2FC < 0] <- "Down"
  
  # 设置颜色
  metab_volcano_colors <- c("Up" = "#E74C3C", "Down" = "#3498DB", "Not Significant" = "#95A5A6")
  
  # 绘制代谢物火山图
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
  
  # 添加显著代谢物的标签（前10个）
  if (nrow(significant_metabolites) > 0) {
    top_metabolites <- head(significant_metabolites, 10)
    metab_volcano_plot <- metab_volcano_plot +
      geom_text(data = metab_volcano_data[metab_volcano_data$metabolite %in% top_metabolites$metabolite, ],
                aes(label = metabolite), size = 3, vjust = -0.5, hjust = 0.5, color = "black")
  }
  
  # 保存代谢物火山图
  ggsave(file.path(output_dir, paste0("Volcano_metabolites_", comparison_name, ".png")), 
         metab_volcano_plot, width = 8, height = 6, dpi = 300)
  cat(sprintf("代谢物火山图已保存: Volcano_metabolites_%s.png\n", comparison_name))
  
  if (nrow(significant_metabolites) > 0) {
    cat("前10个差异最大的代谢物:\n")
    print(head(significant_metabolites, 10))
    
    # 保存结果
    write.csv(significant_metabolites, 
              file.path(output_dir, paste0("significant_metabolites_", comparison_name, ".csv")))
  } else {
    cat("没有找到差异显著的代谢物\n")
  }
  
  # 保存所有代谢物结果（包括不显著的）
  write.csv(metab_results, 
            file.path(output_dir, paste0("all_metabolites_results_", comparison_name, ".csv")))
}

# 保存代谢物名称列表
write.csv(data.frame(metabolite_name = final_metabolite_names), 
          file.path(output_dir, "metabolite_names.csv"), row.names = FALSE)

cat("\n========== 代谢物差异分析完成 ==========\n")

cat("\n分析完成！结果保存在", output_dir, "目录中\n")
