# 导入所需库
library(dplyr)
library(tidyr)
library(tibble)
library(Hmisc)
library(igraph)
library(ggraph)
library(data.table)

plot_multi_omics_network <- function(r_threshold = 0.2, p_threshold = 0.05, 
                                     species_filter = NULL, ko_filter = NULL, metabolite_filter = NULL) {
  
  cat("加载数据集...\n")
  merged_data <- fread(
  "dataset/merged_dataset_processed.csv",
  stringsAsFactors = FALSE,
  data.table = FALSE
)
  
  # 如果我们有独立的特征列表也可以用来筛选：
  species_cols <- grep("^s__", colnames(merged_data), value = TRUE)
  ko_cols <- grep("^K\\d{5}", colnames(merged_data), value = TRUE)
  # 其他视为代谢物（排除分组或元数据列）
  metadata_cols <- c("Group", "SampleID", "TNM_Stage", "Smoking_Status", "Smoking_Status_New")
  metabolite_cols <- setdiff(colnames(merged_data), c(species_cols, ko_cols, metadata_cols))
  
  if(!is.null(species_filter)) species_cols <- intersect(species_cols, species_filter)
  if(!is.null(ko_filter)) ko_cols <- intersect(ko_cols, ko_filter)
  if(!is.null(metabolite_filter)) metabolite_cols <- intersect(metabolite_cols, metabolite_filter)
  
  cat("计算相关性中...\n")
  all_features <- c(species_cols, ko_cols, metabolite_cols)
  mat <- as.matrix(merged_data[, all_features])
  
  # 使用 rcorr 计算 Spearman 相关系数。可以根据需要改为 Pearson。
  res <- rcorr(mat, type = "spearman")
  
  r_mat <- res$r
  p_mat <- res$P
  
  # 转换为长格式
  r_long <- as.data.frame(as.table(r_mat)) %>% 
    rename(Var1 = 1, Var2 = 2, r = Freq)
  p_long <- as.data.frame(as.table(p_mat)) %>% 
    rename(p = Freq)
  
  edges <- r_long %>%
    left_join(p_long, by = c("Var1", "Var2")) %>%
    filter(Var1 != Var2) %>% 
    # 只保留下三角或上三角以防重复
    rowwise() %>%
    mutate(edge_id = paste(sort(c(as.character(Var1), as.character(Var2))), collapse = "_")) %>%
    distinct(edge_id, .keep_all = TRUE) %>%
    ungroup()
  
  # 筛选边
  edges_filtered <- edges %>%
    filter(!is.na(p), p < p_threshold, abs(r) > r_threshold) %>%
    mutate(
      Correlation = ifelse(r > 0, "Positive", "Negative"),
      weight = abs(r)
    )
  
  # 提取存在的节点
  active_nodes <- unique(c(as.character(edges_filtered$Var1), as.character(edges_filtered$Var2)))
  
  nodes <- data.frame(Id = active_nodes, stringsAsFactors = FALSE) %>%
    mutate(
      Type = case_when(
        Id %in% species_cols ~ "Species",
        Id %in% ko_cols ~ "KO Gene",
        Id %in% metabolite_cols ~ "Metabolite",
        TRUE ~ "Other"
      )
    )
    
  # 构建网络图对象
  cat("构建绘图对象...\n")
  net <- graph_from_data_frame(d = edges_filtered, vertices = nodes, directed = FALSE)
  
  # 绘制
  p <- ggraph(net, layout = 'fr') + 
    geom_edge_link(aes(color = Correlation, edge_width = weight), alpha = 0.6) +
    scale_edge_color_manual(values = c("Positive" = "red", "Negative" = "blue")) +
    scale_edge_width_continuous(range = c(0.5, 2), guide = "none") +
    geom_node_point(aes(shape = Type, color = Type), size = 3, alpha = 0.9) +
    scale_shape_manual(values = c("Species" = 16,     # 圆形
                                  "KO Gene" = 17,     # 三角形
                                  "Metabolite" = 15)) + # 方形
    theme_void() +
    theme(legend.position = "right") +
    ggtitle("Integrative Network: Species, KO Genes, and Metabolites")
  
  print(p)
  
  dir.create("results/network_plots", showWarnings = FALSE, recursive = TRUE)
  ggsave("results/network_plots/multi_omics_network.pdf", plot = p, width = 10, height = 8)
  ggsave("results/network_plots/multi_omics_network.png", plot = p, width = 10, height = 8, dpi = 300)
  cat("网络图已保存至 results/network_plots 目录。\n")
  return(net)
}

# 预定义的 KO、Species 和 Metabolite 列表
selected_ko_desc <- c(
  K12688 = "Lipid metabolism",
  K02548 = "Transport system",
  K07091 = "Lipid / stress response",
  K10200 = "Oxidative phosphorylation"
)

FIXED_SPECIES_LIST <- paste0("s__", c(
  "Peptostreptococcus_stomatis",
  "Porphyromonas_gingivalis",
  "Prevotella_intermedia",
  "Fusobacterium_periodonticum",
  "Campylobacter_rectus"
))

FIXED_METABOLITES_LIST <- c(
  "SQDG 26:2; SQDG(13:1/13:1)",
  "Cytosine",
  "Perfluorooctanesulfonic acid",
  "Methyl dihydrojasmonate",
  "Pyrogallol-2-O-sulphate",
  "5'-(3',4'-Dihydroxyphenyl)-gamma-valerolactone sulfate",
  "2-Hydroxy-4,7-dimethoxy-2H-1,4-benzoxazin-3(4H)-one",
  "trans-3,5-Dimethoxy-4-hydroxycinnamaldehyde",
  "(R)-3-Hydroxy-5-phenylpentanoic acid",
  "N-Methyl-D-glucamine",
  "Chenodeoxycholic acid sulfate",
  "Creatinine",
  "Lucidenic acid F",
  "Demissidine",
  "Alpha-Hydroxyisobutyric acid",
  "Pyrocatechol",
  "Gentisic acid",
  "D-Galacturonic acid",
  "1,3-Dimethyluric acid",
  "4-Hydroxy-5-(phenyl)-valeric acid-O-sulphate"
)

# 使用预定义列表绘制网络图
plot_multi_omics_network(
  r_threshold = 0.2, 
  p_threshold = 0.05,
  species_filter = FIXED_SPECIES_LIST,
  ko_filter = names(selected_ko_desc),
  metabolite_filter = FIXED_METABOLITES_LIST
)
