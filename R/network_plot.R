library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(data.table)
library(Hmisc)
library(igraph)
library(ggraph)
library(ggforce)

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
output_dir <- file.path("results", "R_plots/correlation_network")
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

target_species <- c(
  "Peptostreptococcus_stomatis",
  "Porphyromonas_gingivalis", 
  "Prevotella_intermedia",
  "Fusobacterium_periodonticum",
  "Campylobacter_rectus",
  "Faecalibacterium_prausnitzii",
  "Roseburia_intestinalis",
  "Eubacterium_rectale",
  "Coprococcus_comes",
  "Ruminococcus_lactaris"
)

target_cols <- paste0("tax_s__", target_species)
available_cols <- target_cols[target_cols %in% colnames(grouped_data)]
missing_cols <- target_cols[!target_cols %in% colnames(grouped_data)]

cat(sprintf("目标物种数量: %d\n", length(target_species)))
cat(sprintf("数据中可用的物种数量: %d\n", length(available_cols)))
if (length(missing_cols) > 0) {
  cat(sprintf("数据中缺失的物种 (%d个):\n", length(missing_cols)))
  for (col in missing_cols) {
    cat(sprintf("  - %s\n", gsub("^tax_s__", "", col)))
  }
}

if (length(available_cols) < 2) {
  stop("可用物种数量不足，无法进行相关性分析")
}

species_data <- grouped_data[, c("SAMPLE_ID", "group", available_cols)]
cat(sprintf("提取了 %d 个物种的数据用于相关性分析\n", length(available_cols)))

simplify_name <- function(name) {
  name_clean <- gsub("^tax_s__", "", name)
  name_clean <- gsub("_", " ", name_clean)
  return(name_clean)
}

all_comparisons <- list(
  c("control", "CRC_well_diff"),
  c("control", "CRC_poor_diff"),
  c("CRC_well_diff", "CRC_poor_diff")
)

for (comp_idx in 1:length(all_comparisons)) {
  group1 <- all_comparisons[[comp_idx]][1]
  group2 <- all_comparisons[[comp_idx]][2]
  comparison_name <- paste(group1, "vs", group2, sep = "_")
  
  cat(sprintf("\n=== 相关性网络分析: %s ===\n", comparison_name))
  
  current_data <- species_data[species_data$group %in% c(group1, group2), ]
  
  if (nrow(current_data) < 10) {
    cat(sprintf("样本数不足，跳过分析: %s (样本数: %d)\n", comparison_name, nrow(current_data)))
    next
  }
  
  cat(sprintf("当前分组样本数: %d\n", nrow(current_data)))
  
  species_matrix <- as.matrix(current_data[, available_cols])
  rownames(species_matrix) <- current_data$SAMPLE_ID
  
  col_sums <- colSums(species_matrix, na.rm = TRUE)
  zero_cols <- col_sums == 0
  if (any(zero_cols, na.rm = TRUE)) {
    cat(sprintf("移除 %d 个全零物种\n", sum(zero_cols, na.rm = TRUE)))
    species_matrix <- species_matrix[, !zero_cols, drop = FALSE]
    available_cols_filtered <- available_cols[!zero_cols]
  } else {
    available_cols_filtered <- available_cols
  }
  
  if (ncol(species_matrix) < 2) {
    cat("过滤后物种数量不足，跳过分析\n")
    next
  }
  
  species_matrix[is.na(species_matrix)] <- 0
  
  cat("计算Spearman相关系数矩阵...\n")
  cor_matrix <- rcorr(species_matrix, type = "spearman")
  cor_r <- cor_matrix$r
  cor_p <- cor_matrix$P
  
  cat("保存相关系数矩阵...\n")
  write.csv(cor_r, file.path(output_dir, paste0("correlation_matrix_r_", comparison_name, ".csv")))
  write.csv(cor_p, file.path(output_dir, paste0("correlation_matrix_p_", comparison_name, ".csv")))
  
  cat("构建网络图...\n")
  
  n_species <- ncol(cor_r)
  threshold_r <- 0.3
  threshold_p <- 0.05
  
  cor_r[is.na(cor_r)] <- 0
  cor_p[is.na(cor_p)] <- 1
  
  adj_matrix <- ifelse(abs(cor_r) > threshold_r & cor_p < threshold_p, 
                       ifelse(cor_r > 0, cor_r, -cor_r), 
                       0)
  adj_matrix[is.na(adj_matrix)] <- 0
  
  G <- graph_from_adjacency_matrix(
    adj_matrix,
    mode = "undirected",
    weighted = TRUE,
    diag = FALSE
  )
  
  E(G)$correlation <- E(G)$weight
  E(G)$abs_correlation <- abs(E(G)$weight)
  
  pos_edges <- sum(E(G)$weight > 0)
  neg_edges <- sum(E(G)$weight < 0)
  cat(sprintf("网络边数: 正相关=%d, 负相关=%d\n", pos_edges, neg_edges))
  
  V(G)$name_display <- sapply(V(G)$name, simplify_name)
  
  vertex_colors <- rep("#808080", nrow(cor_r))
  species_short <- sapply(available_cols_filtered, function(x) {
    name <- gsub("^tax_s__", "", x)
    if (grepl("Peptostreptococcus|Porphyromonas|Prevotella_intermedia|Fusobacterium|Campylobacter", name)) {
      return("#E74C3C")
    } else if (grepl("Faecalibacterium|Roseburia|Eubacterium|Coprococcus|Ruminococcus", name)) {
      return("#27AE60")
    } else {
      return("#3498DB")
    }
  })
  vertex_colors <- species_short
  
  group_labels <- c(
    "control" = "Control",
    "CRC_well_diff" = "CRC\nWell-diff", 
    "CRC_poor_diff" = "CRC\nPoor-diff"
  )
  
  pdf(file.path(output_dir, paste0("correlation_network_", comparison_name, ".pdf")), 
      width = 14, height = 12)
  
  layout_use <- layout_with_kk(G, dim = 2)
  
  edge_colors <- ifelse(E(G)$weight > 0, 
                        rgb(0.2, 0.6, 0.2, 0.6), 
                        rgb(0.8, 0.2, 0.2, 0.6))
  edge_widths <- pmax(0.5, abs(E(G)$weight) * 4)
  
  par(mar = c(1, 1, 3, 1))
  plot(G, 
       layout = layout_use,
       vertex.color = vertex_colors,
       vertex.size = 20,
       vertex.label = V(G)$name_display,
       vertex.label.cex = 0.7,
       vertex.label.dist = 0.5,
       vertex.frame.color = "white",
       edge.color = edge_colors,
       edge.width = edge_widths,
       main = paste("Species Correlation Network -", comparison_name),
       sub = sprintf("Spearman |r| > %.1f, p < %.2f", threshold_r, threshold_p))
  
  legend("topright", 
         legend = c("Periodontal pathogens (CRC-assoc.)", "Butyrate producers (Beneficial)"),
         col = c("#E74C3C", "#27AE60"), 
         pch = 19, 
         bty = "n",
         pt.cex = 1.5)
  
  dev.off()
  
  cat(sprintf("网络图(PDF)已保存: correlation_network_%s.pdf\n", comparison_name))
  
  p_plot <- ggplot() +
    geom_net(
      aes(x = x, y = y, from = from, to = to, 
          colour = correlation, linewidth = abs_correlation),
      data = tidygraph::as_tidygraph(G) %>% tidygraph::activate(edges) %>% 
        tibble::as_tibble() %>% mutate(x = rep(0, nrow(.)), y = rep(0, nrow(.))),
      layout = layout_in_circle(G),
      curved = TRUE,
      alpha = 0.7
    ) +
    geom_node_point(aes(size = degree(G)), color = vertex_colors, alpha = 0.8) +
    geom_node_text(aes(label = V(G)$name_display), repel = TRUE, size = 3) +
    scale_color_gradient2(low = "#C0392B", mid = "#FDFEFE", high = "#27AE60", 
                          midpoint = 0, name = "Correlation") +
    labs(title = paste("Correlation Network -", comparison_name),
         subtitle = sprintf("Spearman |r| > %.1f, p < %.2f | Edges: %d (pos:%d, neg:%d)", 
                           threshold_r, threshold_p, 
                           pos_edges + neg_edges, pos_edges, neg_edges)) +
    theme_void() +
    theme(legend.position = "right")
  
  ggsave(file.path(output_dir, paste0("correlation_network_ggplot_", comparison_name, ".png")),
         p_plot, width = 14, height = 12, dpi = 300)
  
  cat(sprintf("网络图(PNG)已保存: correlation_network_ggplot_%s.png\n", comparison_name))
  
  corrplot_data <- cor_r
  col_func <- colorRampPalette(c("#C0392B", "#FDFEFE", "#27AE60"))
  
  png(file.path(output_dir, paste0("correlation_heatmap_", comparison_name, ".png")),
      width = 12, height = 10, res = 300)
  
  par(mar = c(8, 8, 3, 3))
  corrplot::corrplot(corrplot_data,
                     method = "color",
                     type = "full",
                     order = "hclust",
                     addrect = 2,
                     col = col_func(100),
                     tl.col = "black",
                     tl.srt = 45,
                     is.corr = FALSE,
                     main = paste("Spearman Correlation -", comparison_name))
  
  dev.off()
  
  cat(sprintf("热图已保存: correlation_heatmap_%s.png\n", comparison_name))
}

cat("\n========== 相关性网络分析完成 ==========\n")
cat("\n分析完成！结果保存在", output_dir, "目录中\n")
