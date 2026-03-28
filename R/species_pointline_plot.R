#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
})

GROUP_COLORS <- c(
  CTRL = "#2E86AB",
  CRC_welldiff = "#F18F01",
  CRC_poordiff = "#D7263D"
)
GROUP_ORDER <- c("CTRL", "CRC_welldiff", "CRC_poordiff")
GROUP_TITLE <- c(
  CTRL = "CTRL",
  CRC_welldiff = "CRC-Well",
  CRC_poordiff = "CRC-Poor"
)

normalize_name <- function(x) {
  tolower(gsub("[^a-z0-9]+", "", as.character(x)))
}

canonicalize_label <- function(x) {
  text <- trimws(gsub('"', "", as.character(x)))
  text <- gsub("^(tax_|met_|kegg_)", "", text, ignore.case = TRUE)
  text <- gsub("^[a-z]__", "", text, ignore.case = TRUE)
  normalize_name(text)
}

pretty_species <- function(x) {
  text <- trimws(gsub('"', "", as.character(x)))
  text <- gsub("^(tax_|met_|kegg_)", "", text, ignore.case = TRUE)
  text <- gsub("^[a-z]__", "", text, ignore.case = TRUE)
  gsub("_", " ", text)
}

find_species_column <- function(df, target_name) {
  species_cols <- names(df)[grepl("^tax_", names(df))]
  target_key <- canonicalize_label(target_name)

  exact_matches <- species_cols[vapply(species_cols, function(col) {
    canonicalize_label(col) == target_key
  }, logical(1))]

  if (length(exact_matches) == 1) {
    return(exact_matches[1])
  }
  if (length(exact_matches) > 1) {
    stop(sprintf("Species %s matched multiple columns: %s",
                 target_name, paste(exact_matches, collapse = ", ")))
  }

  fallback_matches <- species_cols[vapply(species_cols, function(col) {
    nzchar(target_key) && grepl(target_key, canonicalize_label(col), fixed = TRUE)
  }, logical(1))]

  if (length(fallback_matches) == 1) {
    return(fallback_matches[1])
  }
  if (length(fallback_matches) > 1) {
    stop(sprintf("Species %s matched multiple columns: %s",
                 target_name, paste(fallback_matches, collapse = ", ")))
  }

  available <- paste(head(species_cols, 20), collapse = ", ")
  stop(sprintf("No column found for species %s. Example available columns: %s",
               target_name, available))
}

pearson_ci <- function(x, y, confidence = 0.95) {
  ct <- suppressWarnings(cor.test(x, y, method = "pearson"))
  r <- unname(ct$estimate)
  p_value <- ct$p.value
  n <- length(x)

  if (n <= 3 || isTRUE(all.equal(abs(r), 1))) {
    return(list(r = r, p = p_value, ci_low = NA_real_, ci_high = NA_real_))
  }

  clipped_r <- min(max(r, -0.999999), 0.999999)
  z <- atanh(clipped_r)
  se <- 1 / sqrt(n - 3)
  z_crit <- qnorm((1 + confidence) / 2)
  ci_low <- tanh(z - z_crit * se)
  ci_high <- tanh(z + z_crit * se)

  list(r = r, p = p_value, ci_low = ci_low, ci_high = ci_high)
}

prepare_crc_diff_groups <- function(df) {
  crc <- suppressWarnings(as.integer(df$crc_label))
  diff <- suppressWarnings(as.integer(df$differentiation))

  grp <- ifelse(is.na(crc), NA_character_, ifelse(
    crc == 1,
    ifelse(diff == 1, "CRC_poordiff", "CRC_welldiff"),
    "CTRL"
  ))

  df$group <- grp
  df <- df[!is.na(df$group), , drop = FALSE]
  df$group <- factor(df$group, levels = GROUP_ORDER, ordered = TRUE)
  df
}

prepare_plot_data <- function(df, species_a, species_b, transform = "log1p") {
  x_col <- find_species_column(df, species_a)
  y_col <- find_species_column(df, species_b)

  plot_df <- df[, c(x_col, y_col, "group"), drop = FALSE]
  names(plot_df) <- c("x", "y", "group")

  plot_df$x <- suppressWarnings(as.numeric(plot_df$x))
  plot_df$y <- suppressWarnings(as.numeric(plot_df$y))
  plot_df$x[is.infinite(plot_df$x)] <- NA_real_
  plot_df$y[is.infinite(plot_df$y)] <- NA_real_
  plot_df <- plot_df[complete.cases(plot_df), , drop = FALSE]

  plot_df <- plot_df[plot_df$group %in% GROUP_ORDER, , drop = FALSE]
  plot_df$group <- factor(as.character(plot_df$group), levels = GROUP_ORDER, ordered = TRUE)

  if (transform == "log1p") {
    plot_df$x <- log1p(pmax(plot_df$x, 0))
    plot_df$y <- log1p(pmax(plot_df$y, 0))
    x_label <- paste0(pretty_species(species_a), " (log1p abundance)")
    y_label <- paste0(pretty_species(species_b), " (log1p abundance)")
  } else if (transform == "none") {
    x_label <- pretty_species(species_a)
    y_label <- pretty_species(species_b)
  } else {
    stop("transform only supports: log1p or none")
  }

  list(plot_df = plot_df, x_col = x_col, y_col = y_col, x_label = x_label, y_label = y_label)
}

build_annotation <- function(df_group) {
  if (nrow(df_group) < 3) {
    return(list(text = "Insufficient data", n = nrow(df_group)))
  }

  stats_out <- pearson_ci(df_group$x, df_group$y, confidence = 0.95)
  ci_text <- if (is.na(stats_out$ci_low) || is.na(stats_out$ci_high)) {
    "95% CI = [NA, NA]"
  } else {
    sprintf("95%% CI = [%.3f, %.3f]", stats_out$ci_low, stats_out$ci_high)
  }

  txt <- paste(
    sprintf("r = %.3f", stats_out$r),
    ci_text,
    sprintf("p = %.2e", stats_out$p),
    sep = "\n"
  )

  list(text = txt, n = nrow(df_group))
}

make_single_plot <- function(group_df, group_name, x_label, y_label, title_text) {
  ann <- build_annotation(group_df)

  p <- ggplot(group_df, aes(x = x, y = y)) +
    geom_point(size = 2.2, alpha = 0.82, color = GROUP_COLORS[[group_name]]) +
    stat_smooth(method = "lm", se = TRUE, level = 0.95,
                color = GROUP_COLORS[[group_name]], linewidth = 0.9) +
    annotate("label", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.05,
             label = ann$text, size = 3.2,
             label.size = 0.3,
             color = "black", fill = "white") +
    labs(
      title = sprintf("%s (n=%d)", GROUP_TITLE[[group_name]], ann$n),
      subtitle = title_text,
      x = x_label,
      y = y_label
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 10),
      axis.title = element_text(size = 11),
      panel.grid = element_blank(),
      panel.border = element_rect(color = "#202020", fill = NA, linewidth = 0.7)
    )

  p
}

build_species_pair_plot <- function(species_a,
                                    species_b,
                                    input_file = "dataset/merged_dataset_processed.csv",
                                    output_dir = "results/species_pointline",
                                    transform = "none") {
  df <- read.csv(input_file, check.names = FALSE, stringsAsFactors = FALSE)
  df <- prepare_crc_diff_groups(df)

  prepared <- prepare_plot_data(df, species_a, species_b, transform = transform)
  plot_df <- prepared$plot_df

  if (nrow(plot_df) < 3) {
    stop("Fewer than 3 points available for plotting.")
  }

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  safe_a <- normalize_name(species_a)
  safe_b <- normalize_name(species_b)

  data_file <- file.path(output_dir, sprintf("%s_vs_%s_pointline_data.csv", safe_a, safe_b))
  combined_file <- file.path(output_dir, sprintf("%s_vs_%s_pointline_grouped.png", safe_a, safe_b))

  write.csv(plot_df, data_file, row.names = FALSE)

  combined_plot <- ggplot(plot_df, aes(x = x, y = y, color = group)) +
    geom_point(size = 2.0, alpha = 0.82) +
    stat_smooth(method = "lm", se = TRUE, level = 0.95, linewidth = 0.85) +
    facet_wrap(~group, nrow = 1, scales = "fixed", labeller = as_labeller(GROUP_TITLE)) +
    scale_color_manual(values = GROUP_COLORS, drop = FALSE) +
    labs(
      title = sprintf("%s vs %s", pretty_species(species_a), pretty_species(species_b)),
      x = prepared$x_label,
      y = prepared$y_label,
      color = "Group"
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
      strip.text = element_text(face = "bold", size = 11),
      panel.grid = element_blank(),
      panel.border = element_rect(color = "#202020", fill = NA, linewidth = 0.7),
      legend.position = "none"
    )

  ann_df <- do.call(rbind, lapply(GROUP_ORDER, function(g) {
    g_df <- plot_df[plot_df$group == g, , drop = FALSE]
    ann <- build_annotation(g_df)
    data.frame(group = g, label = ann$text, stringsAsFactors = FALSE)
  }))
  ann_df$group <- factor(ann_df$group, levels = GROUP_ORDER, ordered = TRUE)

  combined_plot <- combined_plot +
    geom_label(
      data = ann_df,
      aes(x = -Inf, y = Inf, label = label),
      inherit.aes = FALSE,
      hjust = -0.05,
      vjust = 1.05,
      size = 3.1,
      label.size = 0.3,
      fill = "white",
      color = "black"
    )

  ggsave(combined_file, combined_plot, width = 18.5, height = 6.0, dpi = 300)

  single_files <- character(0)
  for (group_name in GROUP_ORDER) {
    group_df <- plot_df[plot_df$group == group_name, , drop = FALSE]
    title_text <- sprintf("%s: %s vs %s",
                          GROUP_TITLE[[group_name]],
                          pretty_species(species_a),
                          pretty_species(species_b))
    p <- make_single_plot(group_df, group_name, prepared$x_label, prepared$y_label, title_text)

    out_file <- file.path(output_dir, sprintf("%s_vs_%s_%s_pointline.png",
                                              safe_a, safe_b, tolower(group_name)))
    ggsave(out_file, p, width = 6.9, height = 6.0, dpi = 300)
    single_files <- c(single_files, out_file)
  }

  message("Saved combined plot: ", combined_file)
  for (f in single_files) {
    message("Saved group plot: ", f)
  }
  message("Saved data: ", data_file)

  invisible(combined_file)
}

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)

  defaults <- list(
    species_a = "Peptostreptococcus_stomatis",
    species_b = "Faecalibacterium_prausnitzii",
    input_file = "dataset/merged_dataset_processed.csv",
    output_dir = "results/R_plots/species_pointline",
    transform = "none"
  )

  if (length(args) == 0) {
    return(defaults)
  }

  i <- 1
  while (i <= length(args)) {
    key <- args[i]
    if (i == length(args)) {
      stop(sprintf("Missing value for argument: %s", key))
    }
    value <- args[i + 1]

    if (key == "--species-a") {
      defaults$species_a <- value
    } else if (key == "--species-b") {
      defaults$species_b <- value
    } else if (key == "--input-file") {
      defaults$input_file <- value
    } else if (key == "--output-dir") {
      defaults$output_dir <- value
    } else if (key == "--transform") {
      defaults$transform <- value
    } else {
      stop(sprintf("Unknown argument: %s", key))
    }
    i <- i + 2
  }

  if (!(defaults$transform %in% c("log1p", "none"))) {
    stop("--transform only supports: log1p or none")
  }

  defaults
}

main <- function() {
  args <- parse_args()
  build_species_pair_plot(
    species_a = args$species_a,
    species_b = args$species_b,
    input_file = args$input_file,
    output_dir = args$output_dir,
    transform = args$transform
  )
}

if (identical(environment(), globalenv())) {
  main()
}
