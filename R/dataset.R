# 安装必要的包（如果还没有安装）
# install.packages(c("readxl", "dplyr", "tidyr", "stringr"))

# 加载必要的库
library(readxl)
library(dplyr)
library(tidyr)
library(stringr)

# 设置文件路径
sample_file <- 'dataset/id_sample.xlsx'
taxonomy_file <- 'dataset/taxonomy_Species_abund.txt'
metabolomics_file <- 'dataset/metabolome_data.csv'

# 1. 从Excel文件加载样本ID
load_sample_ids <- function() {
  cat("加载样本信息...\n")
  
  tryCatch({
    # 读取两个sheet
    df_patients <- read_excel(sample_file, sheet = 1)  # 患者
    df_controls <- read_excel(sample_file, sheet = 2)  # 健康对照组
    
    # 提取SAMPLE_ID列，并去除空格
    patient_ids <- trimws(as.character(na.omit(df_patients$SAMPLE_ID)))
    control_ids <- trimws(as.character(na.omit(df_controls$SAMPLE_ID)))
    
    # 合并所有样本ID
    sample_ids <- c(patient_ids, control_ids)
    
    cat(sprintf("加载了 %d 个患者样本和 %d 个健康对照样本\n", 
                length(patient_ids), length(control_ids)))
    
    return(list(
      df_patients = df_patients,
      df_controls = df_controls,
      sample_ids = sample_ids
    ))
    
  }, error = function(e) {
    stop(paste("加载样本ID失败:", e$message))
  })
}

# 2. 加载宏基因组数据
load_taxonomy_data <- function() {
  cat("加载宏基因组数据...\n")
  
  tryCatch({
    # 读取宏基因组数据
    taxonomy_data <- read.delim(taxonomy_file, sep = "\t", stringsAsFactors = FALSE)
    
    # 去除列名中的空格
    colnames(taxonomy_data) <- trimws(colnames(taxonomy_data))
    
    # 去除unclassified
    species_col <- colnames(taxonomy_data)[1]
    taxonomy_data <- taxonomy_data[!grepl('unclassified', taxonomy_data[[species_col]]), ]
    taxonomy_data <- taxonomy_data[!grepl('_sp', taxonomy_data[[species_col]]), ]
    taxonomy_data <- taxonomy_data[!grepl('\\[',  taxonomy_data[[species_col]]), ]
    
    cat(sprintf("宏基因组数据形状: %d行 x %d列\n", 
                nrow(taxonomy_data), ncol(taxonomy_data)))
    
    return(taxonomy_data)
    
  }, error = function(e) {
    stop(paste("加载宏基因组数据失败:", e$message))
  })
}

# 3. 加载代谢组数据
load_metabolomics_data <- function() {
  cat("加载代谢组数据...\n")
  
  tryCatch({
    # 读取代谢组数据
    metabolomics_data <- read.csv(metabolomics_file, 
                                  sep = ",", 
                                  stringsAsFactors = FALSE,
                                  check.names = FALSE)
    
    # 去除列名中的空格
    colnames(metabolomics_data) <- trimws(colnames(metabolomics_data))
    
    cat(sprintf("代谢组数据形状: %d行 x %d列\n", 
                nrow(metabolomics_data), ncol(metabolomics_data)))
    
    # 检查第一行是否为标签行
    if (nrow(metabolomics_data) >= 1 && 
        tolower(trimws(metabolomics_data[1, 1])) %in% c("lable", "label")) {
      cat("检测到标签行，将在处理时跳过\n")
    }
    
    return(metabolomics_data)
    
  }, error = function(e) {
    stop(paste("加载代谢组数据失败:", e$message))
  })
}

# 4. 提取吸烟标签
extract_smoking_labels <- function(df_patients, df_controls) {
  cat("提取吸烟标签...\n")
  
  smoking_labels <- list()
  
  tryCatch({
    # 对于患者样本，从sheet1获取吸烟信息
    if (!is.null(df_patients)) {
      # 查找包含'吸烟'的列
      smoking_cols <- grep('吸烟', colnames(df_patients), value = TRUE)
      
      if (length(smoking_cols) > 0) {
        # 使用第一个包含'吸烟'的列
        smoking_col <- smoking_cols[1]
        cat(sprintf("使用列 '%s' 获取吸烟信息\n", smoking_col))
        
        # 提取患者样本的吸烟标签
        for (i in 1:nrow(df_patients)) {
          sample_id <- trimws(as.character(df_patients$SAMPLE_ID[i]))
          smoking_value <- df_patients[[smoking_col]][i]
          
          # 处理缺失值
          if (is.na(smoking_value)) {
            smoking_labels[[sample_id]] <- 0  # 默认不吸烟
          } else {
            # 转换为整数
            tryCatch({
              smoking_labels[[sample_id]] <- as.integer(as.numeric(smoking_value))
            }, error = function(e) {
              smoking_labels[[sample_id]] <- 0
            })
          }
        }
      } else {
        cat("警告: 未找到吸烟信息列，所有患者样本标记为0\n")
        patient_ids <- trimws(as.character(na.omit(df_patients$SAMPLE_ID)))
        for (sample_id in patient_ids) {
          smoking_labels[[sample_id]] <- 0
        }
      }
    }
    
    if (!is.null(df_controls)) {
      control_ids <- trimws(as.character(na.omit(df_controls$SAMPLE_ID)))
      for (sample_id in control_ids) {
        smoking_labels[[sample_id]] <- 0
      }
    }
    
    cat(sprintf("提取了 %d 个样本的吸烟标签\n", length(smoking_labels)))
    
  }, error = function(e) {
    cat(paste("提取吸烟标签时出错:", e$message, "\n"))
  })
  
  return(smoking_labels)
}

# 5. 提取性别标签
extract_gender_labels <- function(df_patients, df_controls) {
  cat("提取性别标签...\n")
  
  gender_labels <- list()
  
  tryCatch({
    # 对于患者样本，从sheet1获取性别信息
    if (!is.null(df_patients)) {
      # 查找包含'性别'的列
      gender_cols <- grep('性别', colnames(df_patients), value = TRUE)
      
      if (length(gender_cols) > 0) {
        # 使用第一个包含'性别'的列
        gender_col <- gender_cols[1]
        cat(sprintf("使用列 '%s' 获取性别信息\n", gender_col))
        
        # 提取患者样本的性别标签
        for (i in 1:nrow(df_patients)) {
          sample_id <- trimws(as.character(df_patients$SAMPLE_ID[i]))
          gender_value <- df_patients[[gender_col]][i]
          
          # 处理缺失值
          if (is.na(gender_value)) {
            gender_labels[[sample_id]] <- 0  # 默认女性
          } else {
            # 转换为整数：0=女性, 1=男性
            tryCatch({
              if (is.character(gender_value)) {
                if (grepl('男', gender_value) || 
                    grepl('male', tolower(gender_value)) ||
                    gender_value == '1') {
                  gender_labels[[sample_id]] <- 1
                } else {
                  gender_labels[[sample_id]] <- 0
                }
              } else {
                gender_labels[[sample_id]] <- as.integer(as.numeric(gender_value))
              }
            }, error = function(e) {
              gender_labels[[sample_id]] <- 0
            })
          }
        }
      } else {
        cat("警告: 未找到性别信息列，所有患者样本标记为0\n")
        patient_ids <- trimws(as.character(na.omit(df_patients$SAMPLE_ID)))
        for (sample_id in patient_ids) {
          gender_labels[[sample_id]] <- 0
        }
      }
    }
    
    # 对于健康对照组，从sheet2获取性别信息
    if (!is.null(df_controls)) {
      # 查找包含'性别'的列
      gender_cols <- grep('性别', colnames(df_controls), value = TRUE)
      
      if (length(gender_cols) > 0) {
        # 使用第一个包含'性别'的列
        gender_col <- gender_cols[1]
        cat(sprintf("使用列 '%s' 获取健康对照组性别信息\n", gender_col))
        
        # 提取健康对照样本的性别标签
        for (i in 1:nrow(df_controls)) {
          sample_id <- trimws(as.character(df_controls$SAMPLE_ID[i]))
          gender_value <- df_controls[[gender_col]][i]
          
          # 处理缺失值
          if (is.na(gender_value)) {
            gender_labels[[sample_id]] <- 0  # 默认女性
          } else {
            tryCatch({
              if (is.character(gender_value)) {
                if (grepl('男', gender_value) || 
                    grepl('male', tolower(gender_value)) ||
                    gender_value == '1') {
                  gender_labels[[sample_id]] <- 1
                } else {
                  gender_labels[[sample_id]] <- 0
                }
              } else {
                gender_labels[[sample_id]] <- as.integer(as.numeric(gender_value))
              }
            }, error = function(e) {
              gender_labels[[sample_id]] <- 0
            })
          }
        }
      } else {
        cat("警告: 未找到健康对照组性别信息列，所有健康对照样本标记为0\n")
        control_ids <- trimws(as.character(na.omit(df_controls$SAMPLE_ID)))
        for (sample_id in control_ids) {
          gender_labels[[sample_id]] <- 0
        }
      }
    }
    
    cat(sprintf("提取了 %d 个样本的性别标签\n", length(gender_labels)))
    
  }, error = function(e) {
    cat(paste("提取性别标签时出错:", e$message, "\n"))
  })
  
  return(gender_labels)
}

# 6. 提取年龄标签
extract_age_labels <- function(df_patients, df_controls) {
  cat("提取年龄标签...\n")
  
  age_labels <- list()
  
  tryCatch({
    # 对于患者样本，从sheet1获取年龄信息
    if (!is.null(df_patients)) {
      # 查找包含'年龄'或'age'的列
      age_cols <- grep('年龄|age', colnames(df_patients), value = TRUE, ignore.case = TRUE)
      
      if (length(age_cols) > 0) {
        # 使用第一个包含'年龄'或'age'的列
        age_col <- age_cols[1]
        cat(sprintf("使用列 '%s' 获取年龄信息\n", age_col))
        
        # 提取患者样本的年龄标签
        for (i in 1:nrow(df_patients)) {
          sample_id <- trimws(as.character(df_patients$SAMPLE_ID[i]))
          age_value <- df_patients[[age_col]][i]
          
          # 处理缺失值
          if (is.na(age_value)) {
            age_labels[[sample_id]] <- 0  # 默认值为0
          } else {
            tryCatch({
              age_labels[[sample_id]] <- as.numeric(age_value)
            }, error = function(e) {
              age_labels[[sample_id]] <- 0
            })
          }
        }
      } else {
        cat("警告: 未找到年龄信息列，所有患者样本标记为0\n")
        patient_ids <- trimws(as.character(na.omit(df_patients$SAMPLE_ID)))
        for (sample_id in patient_ids) {
          age_labels[[sample_id]] <- 0
        }
      }
    }
    
    # 对于健康对照组，从sheet2获取年龄信息
    if (!is.null(df_controls)) {
      # 查找包含'年龄'或'age'的列
      age_cols <- grep('年龄|age', colnames(df_controls), value = TRUE, ignore.case = TRUE)
      
      if (length(age_cols) > 0) {
        # 使用第一个包含'年龄'或'age'的列
        age_col <- age_cols[1]
        cat(sprintf("使用列 '%s' 获取健康对照组年龄信息\n", age_col))
        
        # 提取健康对照样本的年龄标签
        for (i in 1:nrow(df_controls)) {
          sample_id <- trimws(as.character(df_controls$SAMPLE_ID[i]))
          age_value <- df_controls[[age_col]][i]
          
          # 处理缺失值
          if (is.na(age_value)) {
            age_labels[[sample_id]] <- 0  # 默认值为0
          } else {
            tryCatch({
              age_labels[[sample_id]] <- as.numeric(age_value)
            }, error = function(e) {
              age_labels[[sample_id]] <- 0
            })
          }
        }
      } else {
        cat("警告: 未找到健康对照组年龄信息列，所有健康对照样本标记为0\n")
        control_ids <- trimws(as.character(na.omit(df_controls$SAMPLE_ID)))
        for (sample_id in control_ids) {
          age_labels[[sample_id]] <- 0
        }
      }
    }
    
    cat(sprintf("提取了 %d 个样本的年龄标签\n", length(age_labels)))
    
  }, error = function(e) {
    cat(paste("提取年龄标签时出错:", e$message, "\n"))
  })
  
  return(age_labels)
}

# 7. 提取TNM分期标签
extract_tnm_labels <- function(df_patients, df_controls) {
  cat("提取TNM分期标签...\n")
  
  tnm_labels <- list()
  
  tryCatch({
    # 对于患者样本，从sheet1获取TNM分期信息
    if (!is.null(df_patients)) {
      # 查找包含'TNM'或'分期'的列
      tnm_cols <- grep('TNM|分期', colnames(df_patients), value = TRUE)
      
      if (length(tnm_cols) > 0) {
        # 使用第一个包含'TNM'或'分期'的列
        tnm_col <- tnm_cols[1]
        cat(sprintf("使用列 '%s' 获取TNM分期信息\n", tnm_col))
        
        # 提取患者样本的TNM分期标签
        for (i in 1:nrow(df_patients)) {
          sample_id <- trimws(as.character(df_patients$SAMPLE_ID[i]))
          tnm_value <- df_patients[[tnm_col]][i]
          
          # 处理缺失值
          if (is.na(tnm_value)) {
            tnm_labels[[sample_id]] <- 0
          } else {
            tryCatch({
              tnm_labels[[sample_id]] <- as.integer(tnm_value)
            }, error = function(e) {
              tnm_labels[[sample_id]] <- 0
            })
          }
        }
      } else {
        cat("警告: 未找到TNM分期信息列，所有患者样本标记为0\n")
        patient_ids <- trimws(as.character(na.omit(df_patients$SAMPLE_ID)))
        for (sample_id in patient_ids) {
          tnm_labels[[sample_id]] <- 0
        }
      }
    }
    
    # 对于健康对照组，标记为0
    if (!is.null(df_controls)) {
      control_ids <- trimws(as.character(na.omit(df_controls$SAMPLE_ID)))
      for (sample_id in control_ids) {
        tnm_labels[[sample_id]] <- 0
      }
    }
    
    cat(sprintf("提取了 %d 个样本的TNM分期标签\n", length(tnm_labels)))
    
  }, error = function(e) {
    cat(paste("提取TNM分期标签时出错:", e$message, "\n"))
  })
  
  return(tnm_labels)
}

# 8. 提取分化标签
extract_differentiation_labels <- function(df_patients, df_controls) {
  cat("提取分化标签...\n")
  
  differentiation_labels <- list()
  
  tryCatch({
    # 对于患者样本，从sheet1获取分化信息
    if (!is.null(df_patients)) {
      # 查找包含'分化'的列
      diff_cols <- grep('分化|Differentiation', colnames(df_patients), value = TRUE)
      
      if (length(diff_cols) > 0) {
        # 使用第一个包含'分化'的列
        diff_col <- diff_cols[1]
        cat(sprintf("使用列 '%s' 获取分化信息\n", diff_col))
        
        # 提取患者样本的分化标签
        for (i in 1:nrow(df_patients)) {
          sample_id <- trimws(as.character(df_patients$SAMPLE_ID[i]))
          diff_value <- df_patients[[diff_col]][i]
          
          # 处理缺失值
          if (is.na(diff_value)) {
            differentiation_labels[[sample_id]] <- 0  # 默认中-高分化
          } else {
            tryCatch({
              if (is.character(diff_value)) {
                # 转换为整数：0=中-高分化, 1=低分化
                if (grepl('低', diff_value) || 
                    grepl('1', diff_value)) {
                  differentiation_labels[[sample_id]] <- 1
                } else {
                  differentiation_labels[[sample_id]] <- 0
                }
              } else {
                differentiation_labels[[sample_id]] <- as.integer(as.numeric(diff_value))
              }
            }, error = function(e) {
              differentiation_labels[[sample_id]] <- 0
            })
          }
        }
      } else {
        cat("警告: 未找到分化信息列，所有患者样本标记为0\n")
        patient_ids <- trimws(as.character(na.omit(df_patients$SAMPLE_ID)))
        for (sample_id in patient_ids) {
          differentiation_labels[[sample_id]] <- 0
        }
      }
    }
    
    # 对于健康对照组，标记为0（无分化信息）
    if (!is.null(df_controls)) {
      control_ids <- trimws(as.character(na.omit(df_controls$SAMPLE_ID)))
      for (sample_id in control_ids) {
        differentiation_labels[[sample_id]] <- 0
      }
    }
    
    cat(sprintf("提取了 %d 个样本的分化标签\n", length(differentiation_labels)))
    
  }, error = function(e) {
    cat(paste("提取分化标签时出错:", e$message, "\n"))
  })
  
  return(differentiation_labels)
}

# 9. 合并数据为数据框
merge_to_dataframe <- function(taxonomy_data, metabolomics_data, 
                               df_patients, df_controls, sample_ids,
                               smoking_labels, gender_labels, age_labels, tnm_labels,
                               differentiation_labels) {
  cat("合并数据...\n")
  
  # 处理宏基因组数据：Species列为索引，转置为样本为行
  tax_df <- taxonomy_data
  species_col <- colnames(tax_df)[1]
  
  # 设置物种名为行名并转置
  rownames(tax_df) <- tax_df[[species_col]]
  tax_df[[species_col]] <- NULL
  tax_t <- as.data.frame(t(tax_df))
  colnames(tax_t) <- paste0("tax_", colnames(tax_t))
  tax_t$SAMPLE_ID <- rownames(tax_t)
  
  # 处理代谢组数据
  meta <- metabolomics_data
  
  # 检查第一行是否为标签行
  if (!is.na(meta[1, 1]) && grepl("label", tolower(meta[1, 1]))) {
    meta <- meta[-1, ]  # 跳过标签行
  }
  
  # 设置代谢物名为行名并转置
  metab_col <- colnames(meta)[1]
  rownames(meta) <- meta[[metab_col]]
  meta[[metab_col]] <- NULL
  meta_t <- as.data.frame(t(meta))
  colnames(meta_t) <- paste0("met_", colnames(meta_t))
  meta_t$SAMPLE_ID <- rownames(meta_t)
  
  # 按样本ID合并（外连接）
  merged <- full_join(tax_t, meta_t, by = "SAMPLE_ID")
  
  # 如果sample_ids存在，按顺序排列
  if (length(sample_ids) > 0) {
    # 确保样本ID格式一致
    ordered_idx <- sample_ids[sample_ids %in% merged$SAMPLE_ID]
    merged <- merged[match(ordered_idx, merged$SAMPLE_ID), ]
  }
  
  # 添加CRC/CTRL标签
  patient_set <- character()
  control_set <- character()
  
  if (!is.null(df_patients)) {
    patient_set <- trimws(as.character(na.omit(df_patients$SAMPLE_ID)))
  }
  
  if (!is.null(df_controls)) {
    control_set <- trimws(as.character(na.omit(df_controls$SAMPLE_ID)))
  }
  
  merged$crc_label <- sapply(merged$SAMPLE_ID, function(sid) {
    sid <- trimws(as.character(sid))
    if (sid %in% patient_set) return(1)
    if (sid %in% control_set) return(0)
    return(0)
  })
  
  # 添加吸烟标签
  smoking_df <- data.frame(
    SAMPLE_ID = names(smoking_labels),
    smoking_label = unlist(smoking_labels)
  )
  
  # 添加性别标签
  gender_df <- data.frame(
    SAMPLE_ID = names(gender_labels),
    gender_label = unlist(gender_labels)
  )
  
  # 添加年龄标签
  age_df <- data.frame(
    SAMPLE_ID = names(age_labels),
    age = unlist(age_labels)
  )
  
  # 添加TNM分期标签
  tnm_df <- data.frame(
    SAMPLE_ID = names(tnm_labels),
    tnm_stage = unlist(tnm_labels)
  )
  
  # 添加分化标签
  diff_df <- data.frame(
    SAMPLE_ID = names(differentiation_labels),
    differentiation = unlist(differentiation_labels)
  )
  
  # 合并所有标签
  merged <- merged %>%
    left_join(smoking_df, by = "SAMPLE_ID") %>%
    left_join(gender_df, by = "SAMPLE_ID") %>%
    left_join(age_df, by = "SAMPLE_ID") %>%
    left_join(tnm_df, by = "SAMPLE_ID") %>%
    left_join(diff_df, by = "SAMPLE_ID")
  
  # 将缺失值替换为0
  merged[is.na(merged)] <- 0
  
  return(merged)
}

# 9. 预处理宏基因组数据
preprocess_taxonomy_data <- function(df,
                                      remove_duplicates = TRUE,
                                      remove_low_expression = TRUE,
                                      expression_threshold = 0.01,
                                      prevalence_threshold = 0.1,
                                      remove_outliers = FALSE,
                                      outlier_method = "iqr",
                                      outlier_threshold = 3.0,
                                      relative_expression = TRUE,
                                      transform = FALSE,
                                      transform_method = "log") {
  cat("\n预处理宏基因组数据...\n")
  original_shape <- c(nrow(df), ncol(df))
  cat(sprintf("原始数据形状: %d行(物种) x %d列(样本)\n", original_shape[1], original_shape[2]))
  
  id_col <- colnames(df)[1]
  sample_cols <- colnames(df)[-1]
  n_samples <- length(sample_cols)
  
  if (remove_duplicates && nrow(df) > 0) {
    dup_rows <- duplicated(df[[id_col]])
    if (sum(dup_rows) > 0) {
      cat(sprintf("发现重复物种: %d 个\n", sum(dup_rows)))
      df <- df[!dup_rows, ]
    }
  }
  
  if (remove_low_expression && length(sample_cols) > 0) {
    expression_counts <- rowSums(df[, sample_cols, drop = FALSE] > expression_threshold, na.rm = TRUE)
    prevalence <- expression_counts / n_samples
    kept_rows <- prevalence >= prevalence_threshold
    df <- df[kept_rows, , drop = FALSE]

    cat(sprintf("过滤低流行度物种后保留: %d 个物种 (流行度阈值: %.1f%%)\n", nrow(df), prevalence_threshold * 100))
  }
  
  if (remove_outliers && nrow(df) > 0 && length(sample_cols) > 0) {
    outlier_count <- 0
    
    for (i in 1:nrow(df)) {
      row_data <- as.numeric(df[i, sample_cols])
      
      if (outlier_method == "iqr") {
        Q1 <- quantile(row_data, 0.25, na.rm = TRUE)
        Q3 <- quantile(row_data, 0.75, na.rm = TRUE)
        IQR <- Q3 - Q1
        lower_bound <- Q1 - outlier_threshold * IQR
        upper_bound <- Q3 + outlier_threshold * IQR
        
        outliers <- !is.na(row_data) & (row_data < lower_bound | row_data > upper_bound)
        df[i, sample_cols][outliers] <- pmax(pmin(row_data[outliers], upper_bound), lower_bound)
        outlier_count <- outlier_count + sum(outliers, na.rm = TRUE)
        
      } else if (outlier_method == "zscore") {
        z_scores <- abs((row_data - mean(row_data, na.rm = TRUE)) / sd(row_data, na.rm = TRUE))
        outliers <- !is.na(z_scores) & z_scores > outlier_threshold
        median_val <- median(row_data, na.rm = TRUE)
        df[i, sample_cols][outliers] <- median_val
        outlier_count <- outlier_count + sum(outliers, na.rm = TRUE)
      }
    }
    
    cat(sprintf("去除异常值 (%s方法): 处理了 %d 个异常值\n", outlier_method, outlier_count))
  }
  
  if (nrow(df) > 0 && length(sample_cols) > 0) {
    feature_df <- df[, sample_cols, drop = FALSE]
    
    if (relative_expression) {
      col_sums <- colSums(feature_df, na.rm = TRUE)
      col_sums[col_sums == 0] <- 1
      feature_df <- as.data.frame(mapply(function(x, s) x / s * 100, feature_df, col_sums))
      cat("宏基因组数据转换为相对丰度\n")
    }

    if (transform && transform_method == "log") {
      feature_df <- log1p(feature_df)
      cat("宏基因组数据变换 (log)\n")
    }
    
    if (transform && transform_method == "log2") {
      feature_df <- log2(feature_df + 1e-6)
      cat("宏基因组数据变换 (log2)\n")
    }

    if (transform && transform_method == "clr") {
      feature_matrix <- as.matrix(feature_df)
      feature_matrix[feature_matrix <= 0] <- 0.5 * min(feature_matrix[feature_matrix > 0], na.rm = TRUE)
      clr_result <- log(feature_matrix / exp(rowMeans(log(feature_matrix), na.rm = TRUE)))
      feature_df <- as.data.frame(clr_result)
      cat("宏基因组数据变换 (CLR)\n")
    }
    
    df[, sample_cols] <- feature_df
  }
  
  final_shape <- c(nrow(df), ncol(df))
  cat(sprintf("宏基因组数据预处理完成: (%d, %d) -> (%d, %d)\n", 
              original_shape[1], original_shape[2], final_shape[1], final_shape[2]))
  
  return(df)
}

# 10. 预处理代谢组数据
# 输入数据格式: 行为代谢物, 列为样本, 第一列为代谢物名称
preprocess_metabolomics_data <- function(df,
                                         remove_duplicates = TRUE,
                                         remove_low_expression = TRUE,
                                         expression_threshold = 100,
                                         prevalence_threshold = 0.1,
                                         remove_outliers = FALSE,
                                         outlier_method = "iqr",
                                         outlier_threshold = 3.0,
                                         relative_expression = FALSE,
                                         transform = FALSE,
                                         transform_method = "log",
                                         scale = FALSE) {
  cat("\n预处理代谢组数据...\n")
  original_shape <- c(nrow(df), ncol(df))
  cat(sprintf("原始数据形状: %d行(代谢物) x %d列(样本)\n", original_shape[1], original_shape[2]))
  
  id_col <- colnames(df)[1]
  sample_cols <- colnames(df)[-1]
  n_samples <- length(sample_cols)
  
  if (remove_duplicates && nrow(df) > 0) {
    dup_rows <- duplicated(df[[id_col]])
    if (sum(dup_rows) > 0) {
      cat(sprintf("发现重复代谢物: %d 个\n", sum(dup_rows)))
      df <- df[!dup_rows, ]
    }
  }

  if (remove_low_expression && length(sample_cols) > 0) {
    expression_counts <- rowSums(df[, sample_cols, drop = FALSE] > expression_threshold, na.rm = TRUE)
    prevalence <- expression_counts / n_samples
    kept_rows <- prevalence >= prevalence_threshold
    df <- df[kept_rows, , drop = FALSE]

    cat(sprintf("过滤低流行度代谢物后保留: %d 个代谢物 (流行度阈值: %.1f%%)\n", nrow(df), prevalence_threshold * 100))
  }
  
  if (remove_outliers && nrow(df) > 0 && length(sample_cols) > 0) {
    outlier_count <- 0
    
    for (i in 1:nrow(df)) {
      row_data <- as.numeric(df[i, sample_cols])
      
      if (outlier_method == "iqr") {
        Q1 <- quantile(row_data, 0.25, na.rm = TRUE)
        Q3 <- quantile(row_data, 0.75, na.rm = TRUE)
        IQR <- Q3 - Q1
        lower_bound <- Q1 - outlier_threshold * IQR
        upper_bound <- Q3 + outlier_threshold * IQR
        
        outliers <- !is.na(row_data) & (row_data < lower_bound | row_data > upper_bound)
        df[i, sample_cols][outliers] <- pmax(pmin(row_data[outliers], upper_bound), lower_bound)
        outlier_count <- outlier_count + sum(outliers, na.rm = TRUE)
        
      } else if (outlier_method == "zscore") {
        z_scores <- abs((row_data - mean(row_data, na.rm = TRUE)) / sd(row_data, na.rm = TRUE))
        outliers <- !is.na(z_scores) & z_scores > outlier_threshold
        median_val <- median(row_data, na.rm = TRUE)
        df[i, sample_cols][outliers] <- median_val
        outlier_count <- outlier_count + sum(outliers, na.rm = TRUE)
      }
    }
    
    cat(sprintf("去除异常值 (%s方法): 处理了 %d 个异常值\n", outlier_method, outlier_count))
  }
  
  if (nrow(df) > 0 && length(sample_cols) > 0) {
    feature_df <- df[, sample_cols, drop = FALSE]
    
    feature_df <- as.data.frame(lapply(feature_df, function(x) {
      x <- as.numeric(x)
      x[is.na(x) | x == ""] <- median(x, na.rm = TRUE)
      x
    }))
    cat("代谢组数据填充缺失值 (中位数)\n")
    
    if (relative_expression) {
      col_sums <- colSums(feature_df, na.rm = TRUE)
      col_sums[col_sums == 0] <- 1
      feature_df <- as.data.frame(mapply(function(x, s) x / s * 100, feature_df, col_sums))
      cat("代谢组数据转换为相对丰度\n")

    }
    if (transform && transform_method == "log") {
      feature_df <- log1p(feature_df)
      cat("代谢组数据变换 (log)\n")
    }
    
    if (transform && transform_method == "log2") {
      feature_df <- log2(feature_df + 1e-6)
      cat("代谢组数据变换 (log2)\n")
    }
    
    if (scale) {
      feature_df <- scale(feature_df)
      cat("代谢组数据标准化 (Z-score)\n")
    }
    
    df[, sample_cols] <- feature_df
  }
  
  final_shape <- c(nrow(df), ncol(df))
  cat(sprintf("代谢组数据预处理完成: (%d, %d) -> (%d, %d)\n", 
              original_shape[1], original_shape[2], final_shape[1], final_shape[2]))
  
  return(df)
}

# 主程序
cat("开始数据处理流程...\n\n")

# 步骤1：加载样本ID
sample_info <- load_sample_ids()
df_patients <- sample_info$df_patients
df_controls <- sample_info$df_controls
sample_ids <- sample_info$sample_ids

# 步骤2：加载宏基因组数据
taxonomy_data <- load_taxonomy_data()

# 步骤3：加载代谢组数据
metabolomics_data <- load_metabolomics_data()

# 步骤4：提取吸烟标签
smoking_labels <- extract_smoking_labels(df_patients, df_controls)

# 步骤5：提取性别标签
gender_labels <- extract_gender_labels(df_patients, df_controls)

# 步骤6：提取年龄标签
age_labels <- extract_age_labels(df_patients, df_controls)

# 步骤7：提取TNM分期标签
tnm_labels <- extract_tnm_labels(df_patients, df_controls)

# 步骤8：提取分化标签
differentiation_labels <- extract_differentiation_labels(df_patients, df_controls)

# 步骤9：预处理数据
taxonomy_data_relative <- preprocess_taxonomy_data(taxonomy_data)

# 步骤10：预处理代谢组数据
metabolomics_data_relative <- preprocess_metabolomics_data(metabolomics_data)

# 步骤9：预处理数据
taxonomy_data_log <- preprocess_taxonomy_data(taxonomy_data,
                                              remove_duplicates = TRUE,
                                              remove_low_expression = TRUE,
                                              expression_threshold = 0.01,
                                              prevalence_threshold = 0.1,
                                              remove_outliers = FALSE,
                                              outlier_method = "iqr",
                                              outlier_threshold = 3.0,
                                              relative_expression = FALSE,
                                              transform = TRUE,
                                              transform_method = "log")

# 步骤10：预处理代谢组数据
metabolomics_data_log <- preprocess_metabolomics_data(metabolomics_data,
                                                      remove_duplicates = TRUE,
                                                      remove_low_expression = TRUE,
                                                      expression_threshold = 100,
                                                      prevalence_threshold = 0.1,
                                                      remove_outliers = FALSE,
                                                      outlier_method = "iqr",
                                                      outlier_threshold = 3.0,
                                                      relative_expression = FALSE,
                                                      transform = TRUE,
                                                      transform_method = "log",
                                                      scale = FALSE)

# 步骤8：合并数据
merged_data_raw <- merge_to_dataframe(taxonomy_data, metabolomics_data,
                                  df_patients, df_controls, sample_ids,
                                  smoking_labels, gender_labels, age_labels, tnm_labels,
                                  differentiation_labels)

merged_data_relative <- merge_to_dataframe(taxonomy_data_relative, metabolomics_data_relative,
                                  df_patients, df_controls, sample_ids,
                                  smoking_labels, gender_labels, age_labels, tnm_labels,
                                  differentiation_labels)

merged_data_processed <- merge_to_dataframe(taxonomy_data_log, metabolomics_data_log,
                                  df_patients, df_controls, sample_ids,
                                  smoking_labels, gender_labels, age_labels, tnm_labels,
                                  differentiation_labels)

# 显示结果
cat("\n数据处理完成！\n")
cat(sprintf("预处理后数据形状: %d行 x %d列\n", nrow(merged_data_processed), ncol(merged_data_processed)))

# 查看标签分布
cat("\nCRC标签分布:\n")
print(table(merged_data_processed$crc_label))

cat("\n吸烟标签分布:\n")
print(table(merged_data_processed$smoking_label))

cat("\n性别标签分布:\n")
print(table(merged_data_processed$gender_label))

cat("\n年龄分布:\n")
print(summary(merged_data_processed$age))

cat("\nTNM分期标签分布:\n")
print(table(merged_data_processed$tnm_stage))

cat("\n分化标签分布:\n")
print(table(merged_data_processed$differentiation))

# 保存结果到CSV文件（可选）
write.csv(merged_data_raw, "dataset/merged_dataset_raw.csv", row.names = FALSE)
cat("\n数据已保存到 'dataset/merged_dataset_raw.csv'\n")

write.csv(merged_data_relative, "dataset/merged_dataset_relative.csv", row.names = FALSE)
cat("\n数据已保存到 'dataset/merged_dataset_relative.csv'\n")

write.csv(merged_data_processed, "dataset/merged_dataset_processed.csv", row.names = FALSE)
cat("\n数据已保存到 'dataset/merged_dataset_processed.csv'\n")
