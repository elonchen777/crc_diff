suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(purrr)
  library(ggplot2)
})

set.seed(20260325)

input_file <- "dataset/merged_dataset_relative.csv"
output_dir <- "results/R_plots/mediation_species_ko_metabolite_diff_group"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

FIXED_KO_LIST <- c(
  "K02548",
  "K10200",
  "K10201",
  "K02114",
  "K01512",
  "K01647",
  "K02804",
  "K03785"
)

FIXED_SPECIES_LIST <- c(
  "Peptostreptococcus_stomatis",
  "Faecalibacterium_prausnitzii"
)

FIXED_METABOLITES_LIST <- c(
  "SQDG 26:2; SQDG(13:1/13:1)",
  "Chenodeoxycholic acid sulfate",
  "4-Hydroxy-5-(phenyl)-valeric acid-O-sulphate",
  "Cholesterol"
)

USE_CRC_ONLY <- FALSE
N_BOOT <- as.integer(Sys.getenv("N_BOOT", unset = "10"))
MIN_COMPLETE_CASES <- 25

normalize_feature_name <- function(x) {
  x <- trimws(as.character(x))
  x <- tolower(x)
  x <- gsub("^tax_s__", "", x)
  x <- gsub("^kegg_", "", x)
  x <- gsub("^met_", "", x)
  x <- gsub("[^a-z0-9]+", "", x)
  x
}

match_targets_to_cols <- function(target_list, available_cols, remove_prefix = "", strict = TRUE) {
  display_names <- gsub(remove_prefix, "", available_cols)
  norm_display <- normalize_feature_name(display_names)

  matched_cols <- character(0)
  missing_targets <- character(0)

  for (target in target_list) {
    norm_target <- normalize_feature_name(target)

    exact_idx <- which(norm_display == norm_target)
    if (length(exact_idx) > 0) {
      matched_cols <- c(matched_cols, available_cols[exact_idx[1]])
      next
    }

    partial_idx <- which(nzchar(norm_target) & grepl(norm_target, norm_display, fixed = TRUE))
    if (length(partial_idx) == 1) {
      matched_cols <- c(matched_cols, available_cols[partial_idx[1]])
      next
    }

    if (length(partial_idx) > 1) {
      message(sprintf("%s 匹配到多个候选列，跳过模糊匹配", target))
    }
    missing_targets <- c(missing_targets, target)
  }

  if (length(missing_targets) > 0) {
    msg <- paste("Target features not found:", paste(missing_targets, collapse = ", "))
    if (strict) stop(msg)
    warning(msg)
  }

  unique(matched_cols)
}

quote_var <- function(x) paste0("`", x, "`")

coef_by_var <- function(fit, var_name) {
  cc <- coef(fit)
  nm <- gsub("`", "", names(cc), fixed = TRUE)
  idx <- which(nm == var_name)
  if (length(idx) == 0) return(NA_real_)
  as.numeric(cc[idx[1]])
}

safe_scale <- function(x) {
  x <- as.numeric(x)
  sd_x <- suppressWarnings(sd(x, na.rm = TRUE))
  if (!is.finite(sd_x) || sd_x == 0) return(rep(0, length(x)))
  as.numeric(scale(x))
}

make_formula <- function(lhs, rhs_vars) {
  rhs <- paste(vapply(rhs_vars, quote_var, character(1)), collapse = " + ")
  as.formula(paste0(quote_var(lhs), " ~ ", rhs))
}

compute_effects_once <- function(dat, x, m1, m2, y, covars) {
  f_m1 <- make_formula(m1, c(x, covars))
  f_m2 <- make_formula(m2, c(x, m1, covars))
  f_y <- make_formula(y, c(x, m1, m2, covars))

  fit_m1 <- lm(f_m1, data = dat)
  fit_m2 <- lm(f_m2, data = dat)
  fit_y <- lm(f_y, data = dat)

  a1 <- coef_by_var(fit_m1, x)
  d21 <- coef_by_var(fit_m2, m1)
  a2 <- coef_by_var(fit_m2, x)
  b1 <- coef_by_var(fit_y, m1)
  b2 <- coef_by_var(fit_y, m2)
  c_prime <- coef_by_var(fit_y, x)

  ind_x_m1_y <- a1 * b1
  ind_x_m2_y <- a2 * b2
  serial <- a1 * d21 * b2
  total_indirect <- ind_x_m1_y + ind_x_m2_y + serial
  total_effect <- c_prime + total_indirect

  list(
    effects = c(
      ind_x_m1_y = ind_x_m1_y,
      ind_x_m2_y = ind_x_m2_y,
      serial = serial,
      total_indirect = total_indirect,
      direct = c_prime,
      total_effect = total_effect
    ),
    paths = c(a1 = a1, d21 = d21, a2 = a2, b1 = b1, b2 = b2, c_prime = c_prime),
    r2_y = summary(fit_y)$r.squared
  )
}

bootstrap_serial_mediation <- function(dat, x, m1, m2, y, covars, n_boot = 500) {
  used_vars <- unique(c(x, m1, m2, y, covars))
  sub <- dat[, used_vars, drop = FALSE]
  complete_idx <- complete.cases(sub)
  sub <- sub[complete_idx, , drop = FALSE]

  if (nrow(sub) < MIN_COMPLETE_CASES) {
    return(tibble(
      species = x,
      ko = m1,
      metabolite = m2,
      n = nrow(sub),
      status = "insufficient_n"
    ))
  }

  # Scale continuous predictors to improve comparability of path coefficients.
  sub[[x]] <- safe_scale(sub[[x]])
  sub[[m1]] <- safe_scale(sub[[m1]])
  sub[[m2]] <- safe_scale(sub[[m2]])
  for (cv in covars) {
    sub[[cv]] <- safe_scale(sub[[cv]])
  }
  sub[[y]] <- as.numeric(sub[[y]])

  obs <- tryCatch(
    compute_effects_once(sub, x, m1, m2, y, covars),
    error = function(e) NULL
  )
  if (is.null(obs)) {
    return(tibble(
      species = x,
      ko = m1,
      metabolite = m2,
      n = nrow(sub),
      status = "model_failed"
    ))
  }

  boot_mat <- matrix(NA_real_, nrow = n_boot, ncol = 6)
  colnames(boot_mat) <- names(obs$effects)

  for (i in seq_len(n_boot)) {
    idx <- sample.int(nrow(sub), size = nrow(sub), replace = TRUE)
    one <- tryCatch(
      compute_effects_once(sub[idx, , drop = FALSE], x, m1, m2, y, covars)$effects,
      error = function(e) rep(NA_real_, 6)
    )
    boot_mat[i, ] <- one
  }

  ci_fn <- function(v) {
    v <- v[is.finite(v)]
    if (length(v) < 10) return(c(NA_real_, NA_real_))
    as.numeric(quantile(v, probs = c(0.025, 0.975), na.rm = TRUE, type = 6))
  }

  p_fn <- function(v) {
    v <- v[is.finite(v)]
    if (length(v) < 10) return(NA_real_)
    2 * min(mean(v <= 0), mean(v >= 0))
  }

  ci_serial <- ci_fn(boot_mat[, "serial"])
  ci_total_ind <- ci_fn(boot_mat[, "total_indirect"])
  ci_direct <- ci_fn(boot_mat[, "direct"])
  ci_total <- ci_fn(boot_mat[, "total_effect"])

  tibble(
    species = x,
    ko = m1,
    metabolite = m2,
    n = nrow(sub),
    status = "ok",
    a1 = obs$paths[["a1"]],
    d21 = obs$paths[["d21"]],
    a2 = obs$paths[["a2"]],
    b1 = obs$paths[["b1"]],
    b2 = obs$paths[["b2"]],
    c_prime = obs$paths[["c_prime"]],
    ind_x_m1_y = obs$effects[["ind_x_m1_y"]],
    ind_x_m2_y = obs$effects[["ind_x_m2_y"]],
    serial_indirect = obs$effects[["serial"]],
    total_indirect = obs$effects[["total_indirect"]],
    direct_effect = obs$effects[["direct"]],
    total_effect = obs$effects[["total_effect"]],
    serial_ci_low = ci_serial[1],
    serial_ci_high = ci_serial[2],
    total_indirect_ci_low = ci_total_ind[1],
    total_indirect_ci_high = ci_total_ind[2],
    direct_ci_low = ci_direct[1],
    direct_ci_high = ci_direct[2],
    total_effect_ci_low = ci_total[1],
    total_effect_ci_high = ci_total[2],
    p_serial_boot = p_fn(boot_mat[, "serial"]),
    p_total_indirect_boot = p_fn(boot_mat[, "total_indirect"]),
    r2_y = obs$r2_y
  )
}

clean_name <- function(x) {
  x <- gsub("^tax_s__", "", x)
  x <- gsub("^kegg_", "", x)
  x <- gsub("^met_", "", x)
  x
}

make_diff_group <- function(df) {
  crc <- suppressWarnings(as.integer(df$crc_label))
  differentiation <- suppressWarnings(as.integer(df$differentiation))

  out <- rep(NA_integer_, nrow(df))
  out[!is.na(crc) & crc == 0] <- 0
  out[!is.na(crc) & crc == 1 & !is.na(differentiation) & differentiation == 0] <- 1
  out[!is.na(crc) & crc == 1 & !is.na(differentiation) & differentiation == 1] <- 2
  out
}

cat("读取数据...\n")
dat <- fread(input_file, data.table = FALSE)

if (USE_CRC_ONLY) {
  dat <- dat %>% filter(!is.na(crc_label), as.integer(crc_label) == 1)
}

if (!all(c("crc_label", "differentiation") %in% colnames(dat))) {
  stop("数据缺少生成 diff_group 所需列: crc_label 或 differentiation")
}

dat$diff_group <- make_diff_group(dat)
dat <- dat %>% filter(!is.na(diff_group))

species_cols_all <- grep("^tax_s__", colnames(dat), value = TRUE)
ko_cols_all <- grep("^kegg_", colnames(dat), value = TRUE)
met_cols_all <- grep("^met_", colnames(dat), value = TRUE)

species_cols <- match_targets_to_cols(FIXED_SPECIES_LIST, species_cols_all, remove_prefix = "^tax_s__", strict = TRUE)
ko_cols <- match_targets_to_cols(FIXED_KO_LIST, ko_cols_all, remove_prefix = "^kegg_", strict = TRUE)
met_cols <- match_targets_to_cols(FIXED_METABOLITES_LIST, met_cols_all, remove_prefix = "^met_", strict = TRUE)

cat(sprintf("目标数量: species=%d, KO=%d, metabolite=%d\n", length(species_cols), length(ko_cols), length(met_cols)))

covariate_candidates <- c("age", "gender_label")
covariates <- covariate_candidates[covariate_candidates %in% colnames(dat)]
cat(sprintf("纳入协变量: %s\n", ifelse(length(covariates) == 0, "无", paste(covariates, collapse = ", "))))

all_numeric_vars <- unique(c(species_cols, ko_cols, met_cols, covariates, "diff_group"))
for (vv in all_numeric_vars) {
  dat[[vv]] <- suppressWarnings(as.numeric(dat[[vv]]))
}

triplets <- expand.grid(
  species = species_cols,
  ko = ko_cols,
  metabolite = met_cols,
  stringsAsFactors = FALSE
)

cat(sprintf("开始中介分析，共 %d 条链路，bootstrap=%d ...\n", nrow(triplets), N_BOOT))

res <- pmap_dfr(
  triplets,
  function(species, ko, metabolite) {
    bootstrap_serial_mediation(
      dat = dat,
      x = species,
      m1 = ko,
      m2 = metabolite,
      y = "diff_group",
      covars = covariates,
      n_boot = N_BOOT
    )
  }
)

res <- res %>%
  mutate(
    species_name = clean_name(species),
    ko_name = clean_name(ko),
    metabolite_name = clean_name(metabolite),
    path_label = paste0(species_name, " -> ", ko_name, " -> ", metabolite_name, " -> diff_group"),
    serial_sig = !is.na(serial_ci_low) & !is.na(serial_ci_high) & serial_ci_low * serial_ci_high > 0
  ) %>%
  arrange(p_serial_boot, desc(abs(serial_indirect)))

out_all <- file.path(output_dir, "serial_mediation_all_triplets.csv")
out_sig <- file.path(output_dir, "serial_mediation_significant_triplets.csv")
fwrite(res, out_all)
fwrite(res %>% filter(status == "ok", serial_sig), out_sig)

plot_df <- res %>%
  filter(status == "ok", is.finite(serial_indirect)) %>%
  arrange(desc(abs(serial_indirect))) %>%
  slice_head(n = 20)

if (nrow(plot_df) > 0) {
  plot_df$path_label <- factor(plot_df$path_label, levels = rev(plot_df$path_label))
  p1 <- ggplot(plot_df, aes(x = serial_indirect, y = path_label, color = serial_indirect > 0)) +
    geom_vline(xintercept = 0, color = "#7A7A7A", linetype = "dashed", linewidth = 0.4) +
    geom_errorbarh(aes(xmin = serial_ci_low, xmax = serial_ci_high), height = 0.2, linewidth = 0.7) +
    geom_point(size = 2.5) +
    scale_color_manual(values = c("TRUE" = "#B2182B", "FALSE" = "#2166AC"), guide = "none") +
    labs(
      title = "Serial mediation effects (top |indirect|)",
      subtitle = "Path: species -> KO -> metabolite -> diff_group",
      x = "Serial indirect effect (bootstrap 95% CI)",
      y = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold")
    )

  ggsave(file.path(output_dir, "serial_mediation_forest_top20.pdf"), p1, width = 12.5, height = 8.5, device = cairo_pdf)
  ggsave(file.path(output_dir, "serial_mediation_forest_top20.png"), p1, width = 12.5, height = 8.5, dpi = 320)
}

best <- res %>%
  filter(status == "ok", is.finite(serial_indirect)) %>%
  arrange(desc(abs(serial_indirect))) %>%
  slice_head(n = 1)

if (nrow(best) == 1) {
  nd <- tibble(
    node = c("Species", "KO", "Metabolite", "DiffGroup"),
    x = c(1, 2, 3, 4),
    y = c(1, 1, 1, 1),
    label = c(best$species_name, best$ko_name, best$metabolite_name, "diff_group")
  )

  ed <- tibble(
    x = c(1, 2, 3, 1),
    xend = c(2, 3, 4, 4),
    y = c(1, 1, 1, 0.90),
    yend = c(1, 1, 1, 0.90),
    label = c(
      paste0("a1=", round(best$a1, 3)),
      paste0("d21=", round(best$d21, 3)),
      paste0("b2=", round(best$b2, 3)),
      paste0("c' =", round(best$c_prime, 3))
    ),
    edge_type = c("main", "main", "main", "direct")
  )

  p2 <- ggplot() +
    geom_curve(
      data = ed,
      aes(x = x, y = y, xend = xend, yend = yend, linetype = edge_type),
      curvature = 0.06,
      arrow = arrow(length = unit(0.14, "inches")),
      linewidth = 0.7,
      color = "#3A3A3A"
    ) +
    geom_label(
      data = ed,
      aes(x = (x + xend) / 2, y = (y + yend) / 2 + ifelse(edge_type == "direct", -0.08, 0.08), label = label),
      size = 3,
      fill = "white",
      label.size = 0.2
    ) +
    geom_label(
      data = nd,
      aes(x = x, y = y, label = label),
      size = 3.2,
      fill = "#F7F7F7",
      label.size = 0.25
    ) +
    annotate(
      "text",
      x = 2.5,
      y = 1.35,
      label = paste0(
        "Serial indirect = ", round(best$serial_indirect, 4),
        "\n95% CI [", round(best$serial_ci_low, 4), ", ", round(best$serial_ci_high, 4), "]"
      ),
      size = 3.5,
      fontface = "bold"
    ) +
    scale_linetype_manual(values = c(main = "solid", direct = "dashed"), guide = "none") +
    coord_cartesian(xlim = c(0.5, 4.5), ylim = c(0.65, 1.5), clip = "off") +
    labs(
      title = "Best serial mediation chain",
      subtitle = "species -> KO -> metabolite -> diff_group"
    ) +
    theme_void(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5)
    )

  ggsave(file.path(output_dir, "serial_mediation_best_chain.pdf"), p2, width = 10.8, height = 4.6, device = cairo_pdf)
  ggsave(file.path(output_dir, "serial_mediation_best_chain.png"), p2, width = 10.8, height = 4.6, dpi = 320)
}

summary_txt <- file.path(output_dir, "serial_mediation_summary.txt")
ok_n <- sum(res$status == "ok")
sig_n <- sum(res$status == "ok" & res$serial_sig, na.rm = TRUE)

summary_lines <- c(
  "Species -> KO -> Metabolite -> diff_group 串联中介分析完成",
  "diff_group coding: control=0, CRC_Well=1, CRC_Poor=2",
  sprintf("总链路数: %d", nrow(triplets)),
  sprintf("成功拟合: %d", ok_n),
  sprintf("串联中介显著(95%%CI不跨0): %d", sig_n),
  sprintf("结果文件: %s", out_all),
  sprintf("显著结果: %s", out_sig)
)
writeLines(summary_lines, summary_txt)

cat("完成。输出目录: ", output_dir, "\n", sep = "")