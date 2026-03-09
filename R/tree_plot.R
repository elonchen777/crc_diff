# 加载必要的库
library(ape)
library(ggplot2)
library(ggtree)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(ggnewscale)
library(data.table)
library(ggtreeExtra)

# 1. 读取数据
cat("读取合并后的数据...\n")
merged_data <- fread("dataset/merged_dataset.csv", stringsAsFactors = FALSE, data.table = FALSE)
output_dir <- file.path("results", "tree_plots")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# 2. 准备分组函数（基于 prepare_crc_smoke_groups 逻辑）
prepare_early_late_groups <- function(df) {
  cat("根据CRC状态和吸烟状态创建分组...\n")
  
  df$group <- sapply(1:nrow(df), function(i) {
    crc <- as.integer(df$crc_label[i])
    smoking <- as.integer(df$smoking_label[i])
    
    if (!is.na(crc) && crc == 1) {
      if (!is.na(smoking) && smoking == 1) {
        return('CRC_smoking')
      } else {
        return('CRC_nonsmoking')
      }
    } else {
      return('CTRL')
    }
  })
  
  df <- df[!is.na(df$group), ]
  
  # 添加 CRC_total 组
  df_total <- df[df$crc_label == 1, ]
  if (nrow(df_total) > 0) {
    df_total$group <- "CRC_total"
    df <- rbind(df, df_total)
  }
  
  cat("分组统计:\n")
  print(table(df$group))
  
  return(df)
}

# 3. 应用分组函数
grouped_data <- prepare_early_late_groups(merged_data)

# 4. 提取宏基因组数据
taxonomy_cols <- grep("^tax_", colnames(grouped_data), value = TRUE)
cat(sprintf("找到 %d 个宏基因组特征\n", length(taxonomy_cols)))

taxonomy_data <- grouped_data[, c("SAMPLE_ID", "group", taxonomy_cols)]

# 5. 提取物种丰度矩阵
species_matrix <- as.matrix(taxonomy_data[, taxonomy_cols])
rownames(species_matrix) <- taxonomy_data$SAMPLE_ID

# 5.1 处理物种数据
cat("处理物种数据...\n")
species_by_sample <- t(species_matrix)
species_by_sample <- species_by_sample[rowSums(species_by_sample) > 0, ]
cat(sprintf("移除全零物种后: %d个物种\n", nrow(species_by_sample)))

# 5.2 预处理物种
cat("预处理物种数据...\n")
# 过滤掉丰度过低的物种（例如，平均丰度 < 0.001）
mean_abundance <- rowMeans(species_by_sample)
abundance_threshold <- 0.001
species_by_sample <- species_by_sample[mean_abundance >= abundance_threshold, ]
cat(sprintf("过滤掉平均丰度 < %.4f 的物种后: %d个物种\n", abundance_threshold, nrow(species_by_sample)))

species_names <- rownames(species_by_sample)

# 6. 解析物种名称
cat("解析物种名称...\n")

simplify_species_name <- function(tax_name) {
  name_clean <- gsub("^tax_", "", tax_name)
  parts <- strsplit(name_clean, "__")[[1]]
  
  if (length(parts) >= 2) {
    species_part <- parts[2]
    species_std <- gsub("_", " ", species_part)
    genus_species <- strsplit(species_std, " ")[[1]]
    
    if (length(genus_species) >= 2) {
      genus <- genus_species[1]
      species <- paste(genus_species[-1], collapse = " ")
      return(list(
        original_name = tax_name,
        display_name = paste(genus, species),
        genus = genus,
        species = species,
        is_valid = TRUE
      ))
    }
  }
  
  return(list(
    original_name = tax_name,
    display_name = tax_name,
    genus = "Unknown",
    species = "sp.",
    is_valid = FALSE
  ))
}

parsed_names <- lapply(species_names, simplify_species_name)
valid_indices <- which(sapply(parsed_names, function(x) x$is_valid))

if (length(valid_indices) == 0) {
  stop("没有有效的物种名称")
}

species_by_sample_valid <- species_by_sample[valid_indices, ]

# 7. 创建分类信息数据框
taxonomy_info <- do.call(rbind, lapply(parsed_names[valid_indices], function(x) {
  data.frame(
    original_name = x$original_name,
    display_name = x$display_name,
    genus = x$genus,
    species = x$species,
    stringsAsFactors = FALSE
  )
}))

# 8. 用户选择要显示的多个属（在这里修改您感兴趣的属）
selected_genera <- c("Fusobacteria", "Actinobacteria", "Proteobacteria",  "Firmicutes")

cat("\n=== 用户选择的属 ===\n")
cat(paste(selected_genera, collapse = ", "), "\n")

# 检查这些属是否存在于数据中
available_genera <- unique(taxonomy_info$genus)
missing_genera <- setdiff(selected_genera, available_genera)

if (length(missing_genera) > 0) {
  cat("警告: 以下属在数据中不存在:", paste(missing_genera, collapse = ", "), "\n")
  selected_genera <- intersect(selected_genera, available_genera)
  cat("实际显示的属:", paste(selected_genera, collapse = ", "), "\n")
}

if (length(selected_genera) == 0) {
  stop("没有找到任何选定的属，请检查属名拼写")
}

# 9. 筛选出属于选定属的物种
cat("\n筛选属于选定属的物种...\n")
selected_species_info <- taxonomy_info[taxonomy_info$genus %in% selected_genera, ]

if (nrow(selected_species_info) == 0) {
  stop("没有找到属于选定属的物种")
}

cat(sprintf("找到 %d 个物种，来自 %d 个属:\n", 
            nrow(selected_species_info), 
            length(unique(selected_species_info$genus))))

# 显示每个属的物种数量
for (genus_name in selected_genera) {
  n_species <- sum(selected_species_info$genus == genus_name)
  if (n_species > 0) {
    species_list <- selected_species_info$display_name[selected_species_info$genus == genus_name]
    cat(sprintf("  %s: %d 个物种 (%s)\n", genus_name, n_species, paste(species_list, collapse = ", ")))
  }
}

# 10. 创建综合树
cat("\n创建综合系统发育树...\n")

# 获取所有物种的显示名称
all_species_names <- selected_species_info$display_name

# 创建基于属-种层次的简易分类/进化树
# 接受物种显示名（可能含空格）和对应的属向量，按属将物种分组构建 Newick 字符串
create_taxonomic_tree <- function(species_names, genera) {
  if (length(species_names) == 0) stop("没有物种用于构建树")

  # 按属分组
  genus_groups <- split(species_names, genera)

  # 为 Newick 构建安全标签（用下划线代替空格）并为每个属生成分组文本
  genus_texts <- vapply(names(genus_groups), FUN.VALUE = "",
    FUN = function(g) {
      spp <- genus_groups[[g]]
      safe_spp <- gsub(" ", "_", spp)
      if (length(safe_spp) == 1) {
        return(safe_spp)
      } else {
        return(paste0("(", paste(safe_spp, collapse = ","), ")"))
      }
    }
  )

  # 把所有属组合并成一棵树
  tree_text <- paste0("(", paste(genus_texts, collapse = ","), ");")
  tree <- read.tree(text = tree_text)

  # 将安全标签（下划线）还原为空格，确保 tip.label 使用原始显示名
  tree$tip.label <- gsub("_", " ", tree$tip.label)

  # 统一设置分支长度为 1（可根据需要调整为更复杂的尺度）
  tree$edge.length <- rep(1, nrow(tree$edge))

  return(tree)
}

# 创建树（基于属-种层次的分类/进化树）
phylo_tree <- create_taxonomic_tree(all_species_names, selected_species_info$genus)
cat(sprintf("创建了包含 %d 个物种的树\n", length(phylo_tree$tip.label)))

# 11. 计算物种总丰度
cat("计算物种总丰度...\n")

# 获取物种原始名称对应的丰度数据
selected_species_original_names <- selected_species_info$original_name
species_abundance_data <- species_by_sample_valid[selected_species_original_names, , drop = FALSE]

# 计算每个物种的总丰度
species_total_abundance <- rowSums(species_abundance_data)
species_abundance_df <- data.frame(
  species = selected_species_info$display_name,
  total_abundance = log10(species_total_abundance + 1),
  genus = selected_species_info$genus,
  stringsAsFactors = FALSE
)

# 确保物种顺序与树一致
species_abundance_df <- species_abundance_df[match(phylo_tree$tip.label, species_abundance_df$species), ]

# 12. 创建属颜色映射
cat("创建属颜色映射...\n")
unique_genera <- unique(selected_species_info$genus)
genus_colors <- setNames(
  colorRampPalette(brewer.pal(8, "Set2"))(length(unique_genera)),
  unique_genera
)

# 13. 准备差异分析数据（用于外圈）
cat("进行差异分析（吸烟 vs 非吸烟，分别对早期和晚期组）...\n")

# 提取样本信息
sample_info <- taxonomy_data[, c("SAMPLE_ID", "group")]

# 将丰度数据转置为样本x物种
abundance_for_analysis <- as.data.frame(t(species_abundance_data))
colnames(abundance_for_analysis) <- selected_species_info$display_name
abundance_for_analysis$SAMPLE_ID <- rownames(abundance_for_analysis)

# 合并样本信息和丰度
analysis_df <- merge(sample_info, abundance_for_analysis, by = "SAMPLE_ID")

# 定义差异分析函数
perform_diff_analysis <- function(df, case_group, control_group = "CTRL") {
  # case_group: 病例组（如 "CRC_smoking", "CRC_nonsmoking", "CRC_total"）
  # control_group: 对照组（默认 "CTRL"）
  
  # 筛选病例组和对照组样本
  sub_df <- df[df$group %in% c(case_group, control_group), ]
  if (nrow(sub_df) == 0) return(NULL)
  
  # 创建分组指示变量（1: 病例组, 0: 对照组）
  sub_df$case <- ifelse(sub_df$group == case_group, 1, 0)
  
  species_cols <- selected_species_info$display_name
  results <- data.frame(species = species_cols, log2FC = NA, pvalue = NA, stringsAsFactors = FALSE)
  
  for (sp in species_cols) {
    y <- sub_df[[sp]]
    if (sum(y > 0) < 3) next  # 太少非零样本跳过
    # 拟合单变量线性模型：丰度 ~ case （不校正年龄性别）
    y_log2 <- log2(y + 1e-5)
    fit <- try(lm(y_log2 ~ case, data = sub_df), silent = TRUE)
    if (!inherits(fit, "try-error")) {
      coef_sum <- summary(fit)$coefficients
      if ("case" %in% rownames(coef_sum)) {
        results[results$species == sp, "log2FC"] <- coef_sum["case", "Estimate"]
        results[results$species == sp, "pvalue"] <- coef_sum["case", "Pr(>|t|)"]
      }
    }
  }
  # BH校正
  results$p_adj <- p.adjust(results$pvalue, method = "BH")
  results$significant <- results$p_adj < 0.05
  results$direction <- ifelse(results$significant, ifelse(results$log2FC > 0, "up", "down"), "ns")
  return(results)
}

# 分别对各组进行分析
diff_smoking <- perform_diff_analysis(analysis_df, "CRC_smoking")
diff_nonsmoking <- perform_diff_analysis(analysis_df, "CRC_nonsmoking")
diff_total <- perform_diff_analysis(analysis_df, "CRC_total")

# 14. 准备环形数据
cat("准备环形数据...\n")

# 最内圈：所有样本平均相对丰度（log10）
mean_abundance <- rowMeans(species_abundance_data)
log10_mean_abundance <- log10(mean_abundance + 1e-5)  # 加伪计数避免负无穷
inner_data <- data.frame(
  species = selected_species_info$display_name,
  value = log10_mean_abundance
)

# 中间圈：CRC_smoking vs CTRL的log2FC及显著性
middle_data <- diff_smoking[, c("species", "log2FC", "direction")]
colnames(middle_data) <- c("species", "value", "direction")

# 外圈：CRC_nonsmoking vs CTRL的log2FC及显著性
outer_data <- diff_nonsmoking[, c("species", "log2FC", "direction")]
colnames(outer_data) <- c("species", "value", "direction")

# 最外圈：CRC_total vs CTRL的log2FC及显著性
extra_outer_data <- diff_total[, c("species", "log2FC", "direction")]
colnames(extra_outer_data) <- c("species", "value", "direction")

# 确保物种顺序与树一致
inner_data <- inner_data[match(phylo_tree$tip.label, inner_data$species), ]
middle_data <- middle_data[match(phylo_tree$tip.label, middle_data$species), ]
outer_data <- outer_data[match(phylo_tree$tip.label, outer_data$species), ]
extra_outer_data <- extra_outer_data[match(phylo_tree$tip.label, extra_outer_data$species), ]

# 15. 绘制带环形条带的综合树
cat("绘制综合环形树...\n")

# 创建基础树
p <- ggtree(phylo_tree, layout = "circular", size = 0.6, open.angle = 5) +
  theme(plot.margin = margin(1, 1, 1, 1, "cm"))

# 添加最内圈：平均丰度（蓝色渐变）
p <- p + 
  new_scale_fill() +
  geom_fruit(
    data = inner_data,
    geom = geom_tile,
    mapping = aes(y = species, x = 0, fill = value),
    width = 0.2,
    offset = 0.1,
    color = NA
  ) +
  scale_fill_gradient(
    name = "log10(Mean Abundance)",
    low = "lightblue",
    high = "steelblue",
    na.value = "grey90"
  ) +
  new_scale_fill()

# 添加中间圈：CRC_smoking vs CTRL的log2FC，用方向着色
p <- p + 
  geom_fruit(
    data = middle_data,
    geom = geom_tile,
    mapping = aes(y = species, x = 0, fill = direction),
    width = 0.2,
    offset = 0.1,
    color = NA
  ) +
  scale_fill_manual(
    name = "CRC_smoking\nvs CTRL",
    values = c("up" = "red", "down" = "blue", "ns" = "grey80"),
    na.value = "grey90",
    labels = c("up" = "Higher in CRC smokers", "down" = "Lower in CRC smokers", "ns" = "Not significant")
  ) +
  new_scale_fill()

# 添加外圈：CRC_nonsmoking vs CTRL的log2FC
p <- p + 
  geom_fruit(
    data = outer_data,
    geom = geom_tile,
    mapping = aes(y = species, x = 0, fill = direction),
    width = 0.2,
    offset = 0.1,
    color = NA
  ) +
  scale_fill_manual(
    name = "CRC_nonsmoking\nvs CTRL",
    values = c("up" = "red", "down" = "blue", "ns" = "grey80"),
    na.value = "grey90",
    labels = c("up" = "Higher in CRC non-smokers", "down" = "Lower in CRC non-smokers", "ns" = "Not significant")
  ) +
  new_scale_fill()

# 添加最外圈：CRC_total vs CTRL的log2FC
p <- p + 
  geom_fruit(
    data = extra_outer_data,
    geom = geom_tile,
    mapping = aes(y = species, x = 0, fill = direction),
    width = 0.2,
    offset = 0.1,
    color = NA
  ) +
  scale_fill_manual(
    name = "CRC_total\nvs CTRL",
    values = c("up" = "red", "down" = "blue", "ns" = "grey80"),
    na.value = "grey90",
    labels = c("up" = "Higher in CRC", "down" = "Lower in CRC", "ns" = "Not significant")
  )

# 添加标题和主题
p <- p +
  labs(
    title = "Combined Phylogenetic Tree with CRC vs CTRL Rings",
    subtitle = paste("Inner: mean abundance; Middle: CRC_smoking vs CTRL; Outer: CRC_nonsmoking vs CTRL; Extra outer: CRC_total vs CTRL")
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18, margin = margin(b = 10)),
    plot.subtitle = element_text(hjust = 0.5, size = 12, margin = margin(b = 20)),
    legend.position = "right",
    legend.box = "vertical",
    legend.spacing.y = unit(0.2, "cm"),
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 10, face = "bold"),
    legend.key.size = unit(0.6, "cm")
  )

# 14. 保存图形
cat("保存图形...\n")

ggsave(
  filename = file.path(output_dir, "combined_phylogenetic_tree.png"),
  plot = p,
  width = 16,
  height = 14,
  dpi = 300,
  bg = "white"
)

cat("综合树保存完成！\n")

# 17. 保存数据
cat("\n保存数据文件...\n")

# 保存树文件
write.tree(phylo_tree, file.path(output_dir, "combined_tree.nwk"))

# 保存分类信息
write.csv(selected_species_info, file.path(output_dir, "selected_species_info.csv"), row.names = FALSE)

# 保存物种丰度数据
if (nrow(species_abundance_df) > 0) {
  write.csv(species_abundance_df, file.path(output_dir, "species_abundance.csv"), row.names = FALSE)
}

# 保存差异分析结果
if (!is.null(diff_smoking)) {
  write.csv(diff_smoking, file.path(output_dir, "diff_crc_smoking_vs_ctrl.csv"), row.names = FALSE)
}
if (!is.null(diff_nonsmoking)) {
  write.csv(diff_nonsmoking, file.path(output_dir, "diff_crc_nonsmoking_vs_ctrl.csv"), row.names = FALSE)
}
if (!is.null(diff_total)) {
  write.csv(diff_total, file.path(output_dir, "diff_crc_total_vs_ctrl.csv"), row.names = FALSE)
}

# 保存属统计信息
genus_summary <- selected_species_info %>%
  group_by(genus) %>%
  summarise(n_species = n(), species_list = paste(display_name, collapse = "; "))
write.csv(genus_summary, file.path(output_dir, "genus_summary.csv"), row.names = FALSE)

cat("所有处理完成！\n")