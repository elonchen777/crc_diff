# 加载必要的库
library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(data.table)
library(microeco)
library(magrittr)

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
output_dir <- file.path("results", "lefse_analysis")
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

# 5. 准备microtable所需的三个数据表
# 5.1 特征表 (otu_table): 行为物种，列为样本 (注意：microeco要求列名是样本名)
cat("准备microtable数据...\n")

# 提取物种丰度数据
taxonomy_data <- grouped_data[, c("SAMPLE_ID", taxonomy_cols), drop = FALSE]
rownames(taxonomy_data) <- taxonomy_data$SAMPLE_ID
taxonomy_data <- taxonomy_data[, -1, drop = FALSE]

# 转置：行为物种，列为样本
otu_table <- t(taxonomy_data)

# 过滤低丰度物种：移除在所有样本中丰度都为0的物种
row_sums <- rowSums(otu_table, na.rm = TRUE)
valid_species <- row_sums > 0
otu_table <- otu_table[valid_species, , drop = FALSE]
cat(sprintf("过滤后保留 %d 个物种\n", nrow(otu_table)))

# 5.2 样本表 (sample_table): 包含分组信息
sample_table <- grouped_data[, c("SAMPLE_ID", "group", "age", "gender_label"), drop = FALSE]
rownames(sample_table) <- sample_table$SAMPLE_ID
sample_table <- sample_table[, -1, drop = FALSE]
colnames(sample_table) <- c("Group", "Age", "Gender")

# 确保Group是因子
sample_table$Group <- as.factor(sample_table$Group)

# 5.3 分类表 (tax_table): 物种分类信息
# 由于数据中只有species level (tax_s__)，我们构建简化的分类表
species_names <- rownames(otu_table)
tax_table <- data.frame(
  Species = species_names,
  stringsAsFactors = FALSE
)

# 解析分类信息 (格式: tax_s__SpeciesName)
tax_table$Kingdom <- "Bacteria"
tax_table$Phylum <- NA
tax_table$Class <- NA
tax_table$Order <- NA
tax_table$Family <- NA
tax_table$Genus <- NA

for (i in 1:nrow(tax_table)) {
  sp_name <- gsub("^tax_s__", "", tax_table$Species[i])
  tax_table$Genus[i] <- sp_name
}

# 设置tax_table的row names
rownames(tax_table) <- tax_table$Species
tax_table <- tax_table[, -1, drop = FALSE]

# 6. 创建microtable对象
cat("创建microtable对象...\n")
# 转换为data.frame格式
otu_table_df <- as.data.frame(otu_table)
dataset <- microtable$new(
  sample_table = sample_table,
  otu_table = otu_table_df,
  tax_table = tax_table
)

cat(sprintf("microtable对象创建成功: %d 样本, %d 物种\n", 
            nrow(dataset$sample_table), nrow(dataset$otu_table)))

# 7. 定义比较组
comparison_groups <- list(
  c("control", "CRC_well_diff"),
  c("control", "CRC_poor_diff"),
  c("CRC_well_diff", "CRC_poor_diff")
)

# 定义分组颜色映射
group_colors <- c("control" = "#2E86AB", "CRC_well_diff" = "#A23B72", 
                  "CRC_poor_diff" = "#F18F01")

# 8. 对每个比较组进行LEFSE分析
cat("\n========== 开始LEFSE差异分析 (microeco) ==========\n")

for (i in 1:length(comparison_groups)) {
  group1 <- comparison_groups[[i]][1]
  group2 <- comparison_groups[[i]][2]
  comparison_name <- paste(group1, "vs", group2, sep = "_")
  
  cat(sprintf("\n=== LEFSE分析组: %s ===\n", comparison_name))
  
  # 筛选当前比较的两组样本
  current_sample_names <- rownames(dataset$sample_table[
    dataset$sample_table$Group %in% c(group1, group2), 
  ])
  
  if (length(current_sample_names) < 5) {
    cat(sprintf("样本数不足，跳过分析: %s\n", comparison_name))
    next
  }
  
  # 创建子数据集
  dataset_sub <- dataset$clone()
  dataset_sub$sample_table <- dataset$sample_table[current_sample_names, , drop = FALSE]
  # otu_table: 行为物种，列为样本，需要按列（样本名）筛选，转换为data.frame
  dataset_sub$otu_table <- as.data.frame(dataset$otu_table[, current_sample_names, drop = FALSE])
  
  # 确保两组都有足够的样本
  group1_n <- sum(dataset_sub$sample_table$Group == group1)
  group2_n <- sum(dataset_sub$sample_table$Group == group2)
  
  cat(sprintf("样本数: %s = %d, %s = %d\n", 
              group1, group1_n, group2, group2_n))
  
  if (group1_n < 3 || group2_n < 3) {
    cat(sprintf("每组至少需要3个样本，当前: %s = %d, %s = %d\n", 
                group1, group1_n, group2, group2_n))
    next
  }
  
  # 运行LEFSE分析
  cat("运行LEFSE分析...\n")
  tryCatch({
    lefse_result <- trans_diff$new(
      dataset = dataset_sub,
      method = "lefse",
      group = "Group",
      alpha = 0.01,
      lefse_subgroup = NULL
    )
    
    # 检查是否有显著结果
    if (is.null(lefse_result$res_diff) || nrow(lefse_result$res_diff) == 0) {
      cat("没有显著差异物种\n")
      next
    }
    
    cat(sprintf("显著差异物种数: %d\n", nrow(lefse_result$res_diff)))
    
    # 保存LEFSE结果
    write.csv(lefse_result$res_diff, 
              file.path(output_dir, paste0("LEFSE_", comparison_name, "_results.csv")),
              row.names = FALSE)
    
    # 打印前10个
    cat("差异最大的物种:\n")
    print(head(lefse_result$res_diff, 10))
    
    # 绘制LEFSE条形图 - 各组前20个LDA物种，左右分布
    cat("绘制LEFSE条形图...\n")
    
    # 获取LEFSE结果
    lefse_df <- lefse_result$res_diff
    
    # 根据分组筛选各组前20个
    group1_species <- lefse_df[lefse_df$Group == group1, ]
    group2_species <- lefse_df[lefse_df$Group == group2, ]
    
    # 按LDA得分排序，取前20
    if (nrow(group1_species) > 0) {
      group1_species <- group1_species[order(group1_species$LDA, decreasing = TRUE), ]
      group1_top20 <- head(group1_species, 20)
    } else {
      group1_top20 <- group1_species
    }
    
    if (nrow(group2_species) > 0) {
      group2_species <- group2_species[order(group2_species$LDA, decreasing = TRUE), ]
      group2_top20 <- head(group2_species, 20)
    } else {
      group2_top20 <- group2_species
    }
    
    # 合并数据
    plot_data <- rbind(group1_top20, group2_top20)
    
    if (nrow(plot_data) > 0) {
      # 提取物种名称（只保留最后一部分，去掉分类前缀）
      plot_data$Species_short <- sapply(plot_data$Taxa, function(x) {
        parts <- strsplit(as.character(x), "\\|")[[1]]
        return(tail(parts, 1))
      })
      
      # 设置颜色
      plot_data$color <- ifelse(plot_data$Group == group1, 
                                group_colors[group1], 
                                group_colors[group2])
      
      # 为group1的LDA值取负，使其显示在左侧
      plot_data$LDA_plot <- ifelse(plot_data$Group == group1, 
                                   -plot_data$LDA, 
                                   plot_data$LDA)
      
      # 按LDA_plot排序，确保图形美观
      plot_data <- plot_data[order(plot_data$LDA_plot), ]
      plot_data$Species_short <- factor(plot_data$Species_short, levels = plot_data$Species_short)
      
      # 绘制条形图
      p_bar <- ggplot(plot_data, aes(x = LDA_plot, y = Species_short, fill = color)) +
        geom_bar(stat = "identity", width = 0.7) +
        scale_fill_identity() +
        geom_vline(xintercept = 0, linetype = "solid", color = "black", size = 0.5) +
        labs(title = paste("LEFSe Analysis -", comparison_name),
             x = "LDA Score (log10)",
             y = "Species") +
        theme_microbiome() +
        theme(
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          axis.text.y = element_text(size = 8),
          plot.margin = margin(t = 20, r = 60, b = 20, l = 10)
        ) +
        # 添加分组标签在侧边
        annotate("text", x = -max(abs(plot_data$LDA_plot)) * 1.15, y = nrow(plot_data)/2, 
                 label = group1, color = group_colors[group1], 
                 fontface = "bold", size = 5, angle = 90) +
        annotate("text", x = max(abs(plot_data$LDA_plot)) * 1.15, y = nrow(plot_data)/2, 
                 label = group2, color = group_colors[group2], 
                 fontface = "bold", size = 5, angle = 90)
      
      ggsave(file.path(output_dir, paste0("LEFSE_barplot_", comparison_name, ".png")),
             p_bar, width = 14, height = 10, dpi = 300)
      cat(sprintf("LEFSE条形图已保存: LEFSE_barplot_%s.png\n", comparison_name))
    }
    
    # 绘制LEFSE Cladogram
    cat("绘制LEFSE Cladogram...\n")
    p_clado <- lefse_result$plot_diff_cladogram(
      use_taxa_num = 100,
      use_feature_num = 30,
      clade_label_level = 6,
      group_order = c(group1, group2)
    )
    
    if (!is.null(p_clado)) {
      ggsave(file.path(output_dir, paste0("LEFSE_cladogram_", comparison_name, ".png")),
             p_clado, width = 12, height = 10, dpi = 300)
      cat(sprintf("LEFSE Cladogram已保存: LEFSE_cladogram_%s.png\n", comparison_name))
    }
    
  }, error = function(e) {
    cat(sprintf("LEFSE分析出错: %s\n", e$message))
  })
}

cat("\n========== LEFSE分析完成 ==========\n")
