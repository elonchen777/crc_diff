# 加载必要的库
library(curatedMetagenomicData)
library(dplyr)
library(tidyr)

# 设置输出目录
output_dir <- "results/curatedMetagenomicData/control_smokers"
if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
}

# 函数：获取所有可用的数据集信息
get_dataset_info <- function() {
    cat("正在获取数据集信息...\n")
    
    # 获取所有样本的元数据
    all_metadata <- curatedMetagenomicData::sampleMetadata
    
    cat(sprintf("找到 %d 个样本\n", nrow(all_metadata)))
    
    # 获取唯一的数据集名称
    dataset_names <- unique(all_metadata$study_name)
    
    cat(sprintf("找到 %d 个数据集\n", length(dataset_names)))
    
    # 检查每个数据集的元数据
    datasets_with_smoking <- c()
    datasets_with_condition <- c()
    
    for (dataset in dataset_names) {
        # cat(sprintf("数据集 %s 元数据列: %s\n", dataset, colnames(all_metadata)))
        
        # 过滤该数据集的样本
        dataset_metadata <- all_metadata %>% filter(study_name == dataset)
        
        # 检查是否包含smoker和study_condition字段
        if ("smoker" %in% names(dataset_metadata) && any(!is.na(dataset_metadata$smoker))) {
            datasets_with_smoking <- c(datasets_with_smoking, dataset)
        }
        
        if ("study_condition" %in% names(dataset_metadata) && any(!is.na(dataset_metadata$study_condition))) {
            datasets_with_condition <- c(datasets_with_condition, dataset)
        }
    }
    
    cat(sprintf("包含smoking字段的数据集: %d\n", length(datasets_with_smoking)))
    cat(sprintf("包含study_condition字段的数据集: %d\n", length(datasets_with_condition)))
    
    # 找到同时包含两个字段的数据集
    common_datasets <- intersect(datasets_with_smoking, datasets_with_condition)
    cat(sprintf("同时包含smoker和study_condition字段的数据集: %d\n", length(common_datasets)))
    
    if (length(common_datasets) > 0) {
        cat("可用数据集:\n")
        for (ds in common_datasets) {
            cat(paste0("- ", ds, "\n"))
        }
    }
    
    return(list(
        all_datasets = dataset_names,
        with_smoking = datasets_with_smoking,
        with_condition = datasets_with_condition,
        common = common_datasets
    ))
}

# 函数：提取符合条件的数据
extract_control_smokers <- function(dataset_name) {
    cat(sprintf("正在处理数据集: %s\n", dataset_name))
    
    tryCatch({
        # 获取完整数据集 - 使用正确的pattern格式
        dataset <- curatedMetagenomicData(paste0(dataset_name, ".relative_abundance"), dryrun = FALSE) %>% mergeData()
        
        # 检查是否成功获取数据
        if (length(dataset) == 0) {
            cat(sprintf("警告: 数据集 %s 没有返回数据，跳过\n", dataset_name))
            return(NULL)
        }
        
        # 提取物种丰度数据和元数据
        abundance_data <- dataset[[1]]
        metadata <- colData(abundance_data)
        
        # 转换为数据框以便于处理
        metadata_df <- as.data.frame(metadata)
        
        # 检查必要的列是否存在
        if (!"study_condition" %in% colnames(metadata_df)) {
            cat("警告: 该数据集不包含study_condition字段\n")
            return(NULL)
        }
        
        if (!"smoker" %in% colnames(metadata_df)) {
            cat("警告: 该数据集不包含smoker字段\n")
            return(NULL)
        }
        
        filtered_samples <- metadata_df %>%
            filter(study_condition == "control" & smoker == "yes" & body_site == "stool")
        
        if (nrow(filtered_samples) == 0) {
            cat("警告: 该数据集中没有符合条件的样本\n")
            return(NULL)
        }
        
        cat(sprintf("找到 %d 个符合条件的样本\n", nrow(filtered_samples)))
        
        # 提取对应的丰度数据
        sample_names <- rownames(filtered_samples)
        abundance_filtered <- abundance_data[, sample_names]
        
        # 转换为数据框
        abundance_df <- as.data.frame(assay(abundance_filtered))
        
        # 添加物种名称
        abundance_df$species <- rownames(abundance_df)
        
        # 重新排列列
        abundance_df <- abundance_df %>%
            select(species, everything())
        
        # 保存结果
        output_file_abundance <- file.path(output_dir, paste0(dataset_name, "_control_smokers_abundance.csv"))
        output_file_metadata <- file.path(output_dir, paste0(dataset_name, "_control_smokers_metadata.csv"))
        
        write.csv(abundance_df, output_file_abundance, row.names = FALSE)
        write.csv(filtered_samples, output_file_metadata, row.names = TRUE)
        
        cat(sprintf("数据已保存到: %s 和 %s\n", output_file_abundance, output_file_metadata))
        
        return(list(
            abundance = abundance_df,
            metadata = filtered_samples,
            dataset_name = dataset_name
        ))
        
    }, error = function(e) {
        cat(sprintf("处理数据集 %s 时出错: %s\n", dataset_name, e$message))
        return(NULL)
    })
}

# 主函数
main <- function() {
    cat("=== 开始提取控制组吸烟者数据 ===\n")
    
    # 获取数据集信息
    dataset_info <- get_dataset_info()
    
    if (length(dataset_info$common) == 0) {
        cat("警告: 没有找到同时包含smoking和study_condition字段的数据集\n")
        cat("将尝试在所有数据集中搜索符合条件的样本...\n")
        
        # 如果没有共同的数据集，尝试在所有数据集中搜索
        all_results <- list()
        for (dataset in dataset_info$all_datasets) {
            result <- extract_control_smokers(dataset)
            if (!is.null(result)) {
                all_results[[dataset]] <- result
            }
        }
        
        if (length(all_results) == 0) {
            cat("在所有数据集中都没有找到符合条件的样本\n")
        } else {
            cat(sprintf("成功从 %d 个数据集中提取了数据\n", length(all_results)))
        }
        
    } else {
        # 在共同的数据集中提取数据
        all_results <- list()
        for (dataset in dataset_info$common) {
            result <- extract_control_smokers(dataset)
            if (!is.null(result)) {
                all_results[[dataset]] <- result
            }
        }
        
        if (length(all_results) == 0) {
            cat("在共同数据集中没有找到符合条件的样本\n")
        } else {
            cat(sprintf("成功从 %d 个数据集中提取了数据\n", length(all_results)))
        }
    }
    
    # 生成汇总报告
    if (length(all_results) > 0) {
        summary_file <- file.path(output_dir, "extraction_summary.txt")
        
        summary_content <- paste(
            "=== 数据提取汇总报告 ===",
            paste("提取时间:", Sys.time()),
            paste("成功提取的数据集数量:", length(all_results)),
            "",
            "提取的数据集详情:",
            sep = "\n"
        )
        
        for (dataset_name in names(all_results)) {
            result <- all_results[[dataset_name]]
            summary_content <- paste(
                summary_content,
                paste0("数据集: ", dataset_name),
                paste0("样本数量: ", nrow(result$metadata)),
                paste0("物种数量: ", nrow(result$abundance)),
                "",
                sep = "\n"
            )
        }
        
        writeLines(summary_content, summary_file)
        cat(sprintf("汇总报告已保存到: %s\n", summary_file))
    }
    
    cat("=== 数据提取完成 ===\n")
}

# 执行主函数
if (interactive()) {
    main()
} else {
    # 如果作为脚本运行
    main()
}