suppressPackageStartupMessages({
  library(data.table)
  library(KEGGREST)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

# ============================================================
# Pathway lift-up from metabolites, then correlate with KO.
# Step1  metabolite -> KEGG ID mapping
# Step2  KEGG pathway annotation
# Step3  build pathway x sample matrix
# Step4  compute pathway score (z-score)
# Step5  pathway vs Diff correlation
# Step6  pathway vs KO correlation + chessboard visualization
# ============================================================

INPUT_MERGED <- "dataset/merged_dataset_processed.csv"
INPUT_COMBINE <- "dataset/rawdata/combine.intensity.xls"
INPUT_KO_SIG <- "results/maaslin2_ko/significant_results.tsv"
OUT_DIR <- "results/pathway_ko_chessboard"
KEGG_CACHE_DIR <- "dataset/KEGG/cache"

Q_CUTOFF <- 0.01
TOP_N_PATHWAY <- 40
TOP_N_KO <- 80
MIN_MET_PER_PATHWAY <- 2
MIN_KO_PER_PATHWAY <- 3
TOP_N_KO_PATHWAY <- 40

# Fixed pathway level2 list used to control displayed pathways in heatmaps.
# If empty, no level2 filtering is applied.
LEVEL2_FIXED_LIST <- c(
  "Carbohydrate metabolism",
  "Amino acid metabolism",
  "Lipid metabolism",
  "Metabolism of cofactors and vitamins",
  "Nucleotide metabolism",
  "Energy metabolism",
  "Membrane transport",
  "Signal transduction"
)

if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)
if (!dir.exists(KEGG_CACHE_DIR)) dir.create(KEGG_CACHE_DIR, recursive = TRUE)

read_tsv_flexible <- function(path) {
  if (!file.exists(path)) {
    stop("File not found: ", path)
  }
  fread(path, sep = "\t", data.table = FALSE, check.names = FALSE)
}

extract_kegg_ids <- function(x) {
  if (is.na(x) || x == "") return(character(0))
  ids <- regmatches(x, gregexpr("C[0-9]{5}", x, perl = TRUE))[[1]]
  unique(ids)
}

safe_numeric_matrix <- function(df) {
  mat <- as.matrix(df)
  suppressWarnings(mode(mat) <- "numeric")
  mat
}

fetch_compound_pathway_links <- function(compounds, cache_dir = KEGG_CACHE_DIR) {
  link_cache_file <- file.path(cache_dir, "kegg_compound_pathway_links.csv")

  link_df <- data.frame(kegg_id = character(0), pathway_id = character(0), stringsAsFactors = FALSE)
  if (file.exists(link_cache_file)) {
    tmp <- fread(link_cache_file, data.table = FALSE)
    if (all(c("kegg_id", "pathway_id") %in% colnames(tmp))) {
      link_df <- unique(tmp[, c("kegg_id", "pathway_id"), drop = FALSE])
    }
  }

  miss_compounds <- setdiff(compounds, unique(link_df$kegg_id))
  if (length(miss_compounds) > 0) {
    message("Fetching missing KEGG compound-pathway links: ", length(miss_compounds))
    batch_size <- 100
    idx <- seq(1, length(miss_compounds), by = batch_size)

    new_rows <- list()
    for (i in idx) {
      batch <- miss_compounds[i:min(i + batch_size - 1, length(miss_compounds))]
      link_vec <- tryCatch(
        KEGGREST::keggLink("pathway", paste0("cpd:", batch)),
        error = function(e) NULL
      )
      if (is.null(link_vec) || length(link_vec) == 0) next

      new_rows[[length(new_rows) + 1]] <- data.frame(
        kegg_id = sub("^cpd:", "", names(link_vec)),
        pathway_id = sub("^path:", "", unname(link_vec)),
        stringsAsFactors = FALSE
      )
    }

    if (length(new_rows) > 0) {
      link_df <- unique(rbind(link_df, do.call(rbind, new_rows)))
      write.csv(link_df, link_cache_file, row.names = FALSE)
    }
  }

  unique(link_df[link_df$kegg_id %in% compounds, c("kegg_id", "pathway_id"), drop = FALSE])
}

fetch_pathway_metadata <- function(pathway_ids, cache_dir = KEGG_CACHE_DIR) {
  pathway_cache_file <- file.path(cache_dir, "kegg_pathway_metadata.csv")
  meta_cols <- c("pathway_id", "pathway_name", "level1", "level2")

  path_meta <- NULL
  if (file.exists(pathway_cache_file)) {
    tmp <- fread(pathway_cache_file, data.table = FALSE)
    needed <- c("pathway_id", "pathway_name")
    if (all(needed %in% colnames(tmp))) {
      if (!("level1" %in% colnames(tmp))) tmp$level1 <- NA_character_
      if (!("level2" %in% colnames(tmp))) tmp$level2 <- NA_character_
      path_meta <- tmp[, meta_cols, drop = FALSE]
    }
  }

  if (is.null(path_meta)) {
    path_meta <- data.frame(
      pathway_id = character(0),
      pathway_name = character(0),
      level1 = character(0),
      level2 = character(0),
      stringsAsFactors = FALSE
    )
  }

  miss <- setdiff(pathway_ids, path_meta$pathway_id)
  if (length(miss) > 0) {
    name_vec <- tryCatch(KEGGREST::keggList("pathway"), error = function(e) NULL)
    if (!is.null(name_vec) && length(name_vec) > 0) {
      add_df <- data.frame(
        pathway_id = sub("^path:", "", names(name_vec)),
        pathway_name = unname(name_vec),
        level1 = NA_character_,
        level2 = NA_character_,
        stringsAsFactors = FALSE
      )
      add_df <- add_df[, meta_cols, drop = FALSE]
      path_meta <- unique(rbind(path_meta, add_df))
      write.csv(path_meta, pathway_cache_file, row.names = FALSE)
    }
  }

  out <- unique(path_meta[path_meta$pathway_id %in% pathway_ids,
                          c("pathway_id", "pathway_name", "level1", "level2"), drop = FALSE])

  miss_name <- setdiff(pathway_ids, out$pathway_id)
  if (length(miss_name) > 0) {
    out <- rbind(
      out,
      data.frame(
        pathway_id = miss_name,
        pathway_name = miss_name,
        level1 = "Unclassified",
        level2 = "Unclassified",
        stringsAsFactors = FALSE
      )
    )
  }

  out$pathway_name[is.na(out$pathway_name) | out$pathway_name == ""] <- out$pathway_id
  out$level1[is.na(out$level1) | out$level1 == ""] <- "Unclassified"
  out$level2[is.na(out$level2) | out$level2 == ""] <- "Unclassified"
  out
}

corr_with_p <- function(x, y, method = "spearman") {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 6) return(c(rho = NA_real_, p = NA_real_, n = sum(ok)))

  ct <- tryCatch(cor.test(x[ok], y[ok], method = method, exact = FALSE), error = function(e) NULL)
  if (is.null(ct)) return(c(rho = NA_real_, p = NA_real_, n = sum(ok)))
  c(rho = unname(ct$estimate), p = ct$p.value, n = sum(ok))
}

filter_by_level2 <- function(df, level2_col = "level2", fixed_list = LEVEL2_FIXED_LIST) {
  if (nrow(df) == 0) return(df)
  if (!(level2_col %in% colnames(df))) return(df)
  if (length(fixed_list) == 0) return(df)

  out <- df[df[[level2_col]] %in% fixed_list, , drop = FALSE]
  if (nrow(out) == 0) {
    stop("No pathways left after LEVEL2_FIXED_LIST filtering. Please update LEVEL2_FIXED_LIST.")
  }
  out
}

message("[1/6] Load merged data and map metabolite -> KEGG IDs")
merged_data <- fread(INPUT_MERGED, data.table = FALSE, check.names = FALSE)
met_cols <- colnames(merged_data)[grepl("^met_", colnames(merged_data))]
ko_cols_all <- colnames(merged_data)[grepl("^kegg_", colnames(merged_data))]

if (!all(c("crc_label", "differentiation") %in% colnames(merged_data))) {
  stop("merged dataset must include columns: crc_label, differentiation")
}

comb <- read_tsv_flexible(INPUT_COMBINE)
for (nm in c("MS2Metabolite", "MS2kegg")) {
  if (!nm %in% colnames(comb)) comb[[nm]] <- ""
}

ms2_name <- trimws(as.character(comb$MS2Metabolite))
ms2_kegg <- trimws(as.character(comb$MS2kegg))

map_tbl <- data.frame(
  feature = make.names(paste0("met_", ms2_name), unique = FALSE),
  metabolite_name = ms2_name,
  ms2kegg = ms2_kegg,
  stringsAsFactors = FALSE
)
map_tbl <- map_tbl[map_tbl$feature %in% met_cols & map_tbl$metabolite_name != "" & map_tbl$ms2kegg != "", , drop = FALSE]

kegg_list <- lapply(map_tbl$ms2kegg, extract_kegg_ids)
expand_n <- lengths(kegg_list)

map_expand <- map_tbl[rep(seq_len(nrow(map_tbl)), expand_n), c("feature", "metabolite_name"), drop = FALSE]
if (nrow(map_expand) > 0) {
  map_expand$kegg_id <- unlist(kegg_list)
  map_expand <- unique(map_expand)
} else {
  stop("No metabolite KEGG IDs mapped from combine file.")
}

write.csv(map_expand, file.path(OUT_DIR, "step1_metabolite_to_kegg_map.csv"), row.names = FALSE)
message("Mapped metabolites with KEGG IDs: ", nrow(map_expand))

message("[2/6] Annotate KEGG pathways")
cpd2path <- fetch_compound_pathway_links(unique(map_expand$kegg_id), cache_dir = KEGG_CACHE_DIR)
if (nrow(cpd2path) == 0) {
  stop("No compound-pathway mapping available.")
}

path_meta <- fetch_pathway_metadata(unique(cpd2path$pathway_id), cache_dir = KEGG_CACHE_DIR)
cpd2path_anno <- merge(cpd2path, path_meta, by = "pathway_id", all.x = TRUE)
cpd2path_anno <- unique(cpd2path_anno)

write.csv(cpd2path_anno, file.path(OUT_DIR, "step2_compound_to_pathway_annotation.csv"), row.names = FALSE)

message("[3/6] Build pathway x sample matrix")
feat2cpd <- unique(map_expand[, c("feature", "kegg_id"), drop = FALSE])
feat2path <- merge(feat2cpd, cpd2path_anno[, c("kegg_id", "pathway_id", "pathway_name", "level1", "level2"), drop = FALSE],
                   by = "kegg_id", all.x = FALSE)
feat2path <- unique(feat2path)

pathway_feature_list <- split(feat2path$feature, feat2path$pathway_id)
pathway_feature_list <- lapply(pathway_feature_list, unique)
pathway_feature_list <- pathway_feature_list[sapply(pathway_feature_list, length) >= MIN_MET_PER_PATHWAY]

if (length(pathway_feature_list) == 0) {
  stop("No pathways pass minimum metabolite count threshold.")
}

pathway_ids <- names(pathway_feature_list)
pathway_mat <- sapply(pathway_ids, function(pid) {
  cols <- pathway_feature_list[[pid]]
  vals <- merged_data[, cols, drop = FALSE]
  rowMeans(vals, na.rm = TRUE)
})

if (is.vector(pathway_mat)) {
  pathway_mat <- matrix(pathway_mat, ncol = 1)
  colnames(pathway_mat) <- pathway_ids[1]
}

pathway_mat <- t(pathway_mat)
colnames(pathway_mat) <- paste0("S", seq_len(ncol(pathway_mat)))

path_info <- unique(feat2path[, c("pathway_id", "pathway_name", "level1", "level2"), drop = FALSE])
path_info <- path_info[match(pathway_ids, path_info$pathway_id), , drop = FALSE]

rownames(pathway_mat) <- paste0(pathway_ids, " | ", path_info$pathway_name)
write.csv(pathway_mat, file.path(OUT_DIR, "step3_pathway_by_sample_matrix.csv"), row.names = TRUE)
write.csv(path_info, file.path(OUT_DIR, "step3_pathway_metadata.csv"), row.names = FALSE)

message("[4/6] Compute pathway scores")
pathway_score <- t(scale(t(pathway_mat)))
pathway_score[!is.finite(pathway_score)] <- 0
write.csv(pathway_score, file.path(OUT_DIR, "step4_pathway_score_matrix_zscore.csv"), row.names = TRUE)

message("[5/6] Correlate pathway score vs Diff")
crc_idx <- merged_data$crc_label == 1 & is.finite(merged_data$differentiation)
diff_num <- as.numeric(merged_data$differentiation)

path_diff_rows <- lapply(seq_len(nrow(pathway_score)), function(i) {
  cc <- corr_with_p(pathway_score[i, crc_idx], diff_num[crc_idx], method = "spearman")
  data.frame(
    pathway_id = pathway_ids[i],
    pathway_name = path_info$pathway_name[i],
    level1 = path_info$level1[i],
    level2 = path_info$level2[i],
    rho = unname(cc["rho"]),
    pval = unname(cc["p"]),
    n = unname(cc["n"]),
    stringsAsFactors = FALSE
  )
})

path_diff <- do.call(rbind, path_diff_rows)
path_diff$qval <- p.adjust(path_diff$pval, method = "BH")
path_diff <- path_diff[order(path_diff$pval), , drop = FALSE]
write.csv(path_diff, file.path(OUT_DIR, "step5_pathway_vs_diff_spearman.csv"), row.names = FALSE)

message("[6/6] Correlate pathway vs KO and draw chessboard heatmap")
ko_cols <- ko_cols_all
if (file.exists(INPUT_KO_SIG)) {
  ko_sig <- read_tsv_flexible(INPUT_KO_SIG)
  if ("feature" %in% colnames(ko_sig)) {
    ko_sig <- ko_sig[grepl("^kegg_K[0-9]{5}$", ko_sig$feature), , drop = FALSE]
    if (nrow(ko_sig) > 0) {
      if ("qval" %in% colnames(ko_sig)) {
        ko_sig$qval <- suppressWarnings(as.numeric(ko_sig$qval))
        ko_sig <- ko_sig[is.finite(ko_sig$qval) & ko_sig$qval < Q_CUTOFF, , drop = FALSE]
      }
      if ("coef" %in% colnames(ko_sig)) {
        ko_sig$coef <- suppressWarnings(as.numeric(ko_sig$coef))
        ko_sig <- ko_sig[order(-abs(ko_sig$coef)), , drop = FALSE]
      }
      ko_pick <- unique(ko_sig$feature)
      if (length(ko_pick) > 0) {
        ko_cols <- intersect(ko_pick, ko_cols_all)
      }
    }
  }
}
if (length(ko_cols) > TOP_N_KO) ko_cols <- ko_cols[seq_len(TOP_N_KO)]

path_pick <- path_diff
path_pick <- path_pick[is.finite(path_pick$pval), , drop = FALSE]
path_pick <- filter_by_level2(path_pick, level2_col = "level2", fixed_list = LEVEL2_FIXED_LIST)
path_pick <- path_pick[order(path_pick$qval, path_pick$pval), , drop = FALSE]
path_pick <- path_pick[seq_len(min(TOP_N_PATHWAY, nrow(path_pick))), , drop = FALSE]

path_sel <- path_pick$pathway_id
path_idx <- match(path_sel, pathway_ids)
path_score_sel <- pathway_score[path_idx, , drop = FALSE]

ko_mat <- safe_numeric_matrix(merged_data[, ko_cols, drop = FALSE])
if (ncol(ko_mat) == 0 || nrow(path_score_sel) == 0) {
  stop("No KO columns or pathway rows selected for correlation.")
}

rho_mat <- matrix(NA_real_, nrow = nrow(path_score_sel), ncol = ncol(ko_mat))
p_mat <- matrix(NA_real_, nrow = nrow(path_score_sel), ncol = ncol(ko_mat))
rownames(rho_mat) <- rownames(path_score_sel)
colnames(rho_mat) <- gsub("^kegg_", "", colnames(ko_mat))
rownames(p_mat) <- rownames(rho_mat)
colnames(p_mat) <- colnames(rho_mat)

for (i in seq_len(nrow(path_score_sel))) {
  for (j in seq_len(ncol(ko_mat))) {
    cc <- corr_with_p(as.numeric(path_score_sel[i, ]), as.numeric(ko_mat[, j]), method = "spearman")
    rho_mat[i, j] <- unname(cc["rho"])
    p_mat[i, j] <- unname(cc["p"])
  }
}

corr_long <- data.frame(
  pathway = rep(rownames(rho_mat), times = ncol(rho_mat)),
  ko = rep(colnames(rho_mat), each = nrow(rho_mat)),
  rho = as.vector(rho_mat),
  pval = as.vector(p_mat),
  stringsAsFactors = FALSE
)
corr_long$qval <- p.adjust(corr_long$pval, method = "BH")
write.csv(corr_long, file.path(OUT_DIR, "step6_pathway_vs_ko_correlations.csv"), row.names = FALSE)

ko_link_cache <- file.path(KEGG_CACHE_DIR, "kegg_ko_pathway_links.csv")
ko_group <- data.frame(ko = colnames(rho_mat), ko_level2 = "Unclassified", stringsAsFactors = FALSE)
if (file.exists(ko_link_cache)) {
  ko_link <- fread(ko_link_cache, data.table = FALSE)
  if (all(c("ko_id", "pathway_id") %in% colnames(ko_link))) {
    ko_map <- data.frame(ko = gsub("^K", "K", ko_link$ko_id), pathway_id = ko_link$pathway_id, stringsAsFactors = FALSE)
    ko_map <- ko_map[ko_map$ko %in% ko_group$ko, , drop = FALSE]

    if (nrow(ko_map) > 0) {
      ko_meta <- fetch_pathway_metadata(unique(ko_map$pathway_id), cache_dir = KEGG_CACHE_DIR)
      ko_map2 <- merge(ko_map, ko_meta[, c("pathway_id", "level2"), drop = FALSE], by = "pathway_id", all.x = TRUE)
      ko_map2$level2[is.na(ko_map2$level2) | ko_map2$level2 == ""] <- "Unclassified"
      ko_map2 <- ko_map2[order(ko_map2$ko, ko_map2$level2), , drop = FALSE]
      ko_map2 <- ko_map2[!duplicated(ko_map2$ko), c("ko", "level2"), drop = FALSE]
      colnames(ko_map2)[2] <- "ko_level2"
      ko_group <- merge(ko_group, ko_map2, by = "ko", all.x = TRUE, suffixes = c("", ".new"))
      ko_group$ko_level2 <- ifelse(!is.na(ko_group$ko_level2.new), ko_group$ko_level2.new, ko_group$ko_level2)
      ko_group$ko_level2.new <- NULL
    }
  }
}

ko_group$ko_level2[is.na(ko_group$ko_level2) | ko_group$ko_level2 == ""] <- "Unclassified"
ko_group <- ko_group[match(colnames(rho_mat), ko_group$ko), , drop = FALSE]

path_group <- path_info[match(path_sel, path_info$pathway_id), c("pathway_id", "level2"), drop = FALSE]
path_group$level2[is.na(path_group$level2) | path_group$level2 == ""] <- "Unclassified"

sig_mat <- matrix("", nrow = nrow(p_mat), ncol = ncol(p_mat))
sig_mat[p_mat < 0.001] <- "***"
sig_mat[p_mat >= 0.001 & p_mat < 0.01] <- "**"
sig_mat[p_mat >= 0.01 & p_mat < 0.05] <- "*"

ko_levels <- unique(ko_group$ko_level2)
ko_cols_palette <- structure(
  hcl(
    h = seq(15, 375, length.out = length(ko_levels) + 1)[-1],
    c = 85,
    l = 65
  ),
  names = ko_levels
)

path_levels <- unique(path_group$level2)
path_cols_palette <- structure(
  hcl(
    h = seq(30, 390, length.out = length(path_levels) + 1)[-1],
    c = 70,
    l = 55
  ),
  names = path_levels
)

ha_top <- HeatmapAnnotation(
  KO_Group = factor(ko_group$ko_level2, levels = ko_levels),
  col = list(KO_Group = ko_cols_palette),
  annotation_name_gp = gpar(fontsize = 9, fontface = "bold"),
  show_legend = TRUE
)

ha_left <- rowAnnotation(
  Pathway_Group = factor(path_group$level2, levels = path_levels),
  col = list(Pathway_Group = path_cols_palette),
  annotation_name_gp = gpar(fontsize = 9, fontface = "bold"),
  show_legend = TRUE
)

col_fun <- colorRamp2(c(-1, -0.5, 0, 0.5, 1), c("#3B6FB6", "#8FB9E2", "#F7F7F7", "#F2A65A", "#C63D2F"))

ht <- Heatmap(
  rho_mat,
  name = "Spearman r",
  col = col_fun,
  top_annotation = ha_top,
  left_annotation = ha_left,
  column_split = factor(ko_group$ko_level2, levels = ko_levels),
  row_split = factor(path_group$level2, levels = path_levels),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_dend = FALSE,
  show_column_dend = FALSE,
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 7),
  column_names_rot = 45,
  rect_gp = gpar(col = "white", lwd = 0.6),
  border = TRUE,
  cell_fun = function(j, i, x, y, width, height, fill) {
    if (sig_mat[i, j] != "") {
      grid.text(sig_mat[i, j], x, y, gp = gpar(fontsize = 6, col = "black"))
    }
  }
)

png(file.path(OUT_DIR, "pathway_ko_chessboard_heatmap.png"), width = 3600, height = 2200, res = 300)
draw(
  ht,
  column_title = "Pathway vs KO Correlation Chessboard (with KO Group Labels)",
  column_title_gp = gpar(fontsize = 14, fontface = "bold"),
  heatmap_legend_side = "right",
  annotation_legend_side = "right",
  merge_legend = FALSE
)
dev.off()

pdf(file.path(OUT_DIR, "pathway_ko_chessboard_heatmap.pdf"), width = 16, height = 10)
draw(
  ht,
  column_title = "Pathway vs KO Correlation Chessboard (with KO Group Labels)",
  column_title_gp = gpar(fontsize = 14, fontface = "bold"),
  heatmap_legend_side = "right",
  annotation_legend_side = "right",
  merge_legend = FALSE
)
dev.off()

write.csv(ko_group, file.path(OUT_DIR, "ko_group_labels.csv"), row.names = FALSE)
write.csv(path_group, file.path(OUT_DIR, "pathway_group_labels.csv"), row.names = FALSE)

message("[7/7] Compute KO pathway scores and draw KO-pathway vs metabolite-pathway chessboard")
ko_link_cache <- file.path(KEGG_CACHE_DIR, "kegg_ko_pathway_links.csv")
if (!file.exists(ko_link_cache)) {
  stop("KO-pathway cache not found: ", ko_link_cache)
}

ko_link_all <- fread(ko_link_cache, data.table = FALSE)
if (!all(c("ko_id", "pathway_id") %in% colnames(ko_link_all))) {
  stop("KO-pathway cache must include columns: ko_id, pathway_id")
}

ko_link_all$ko_id <- sub("^ko:", "", ko_link_all$ko_id)
ko_link_all$ko_col <- paste0("kegg_", ko_link_all$ko_id)
ko_link_all <- unique(ko_link_all[, c("ko_col", "pathway_id"), drop = FALSE])
ko_link_all <- ko_link_all[ko_link_all$ko_col %in% ko_cols_all, , drop = FALSE]

if (nrow(ko_link_all) == 0) {
  stop("No KO columns in merged data can be mapped to KEGG pathways.")
}

ko_path_meta <- fetch_pathway_metadata(unique(ko_link_all$pathway_id), cache_dir = KEGG_CACHE_DIR)
ko_feat2path <- merge(
  ko_link_all,
  ko_path_meta[, c("pathway_id", "pathway_name", "level1", "level2"), drop = FALSE],
  by = "pathway_id",
  all.x = TRUE
)
ko_feat2path <- unique(ko_feat2path)

ko_pathway_feature_list <- split(ko_feat2path$ko_col, ko_feat2path$pathway_id)
ko_pathway_feature_list <- lapply(ko_pathway_feature_list, unique)
ko_pathway_feature_list <- ko_pathway_feature_list[sapply(ko_pathway_feature_list, length) >= MIN_KO_PER_PATHWAY]

if (length(ko_pathway_feature_list) == 0) {
  stop("No KO pathways pass minimum KO count threshold.")
}

ko_pathway_ids <- names(ko_pathway_feature_list)
ko_pathway_mat <- sapply(ko_pathway_ids, function(pid) {
  cols <- ko_pathway_feature_list[[pid]]
  vals <- merged_data[, cols, drop = FALSE]
  rowMeans(vals, na.rm = TRUE)
})

if (is.vector(ko_pathway_mat)) {
  ko_pathway_mat <- matrix(ko_pathway_mat, ncol = 1)
  colnames(ko_pathway_mat) <- ko_pathway_ids[1]
}

ko_pathway_mat <- t(ko_pathway_mat)
colnames(ko_pathway_mat) <- paste0("S", seq_len(ncol(ko_pathway_mat)))

ko_path_info <- unique(ko_feat2path[, c("pathway_id", "pathway_name", "level1", "level2"), drop = FALSE])
ko_path_info <- ko_path_info[match(ko_pathway_ids, ko_path_info$pathway_id), , drop = FALSE]
rownames(ko_pathway_mat) <- paste0(ko_pathway_ids, " | ", ko_path_info$pathway_name)

write.csv(ko_pathway_mat, file.path(OUT_DIR, "step7_ko_pathway_by_sample_matrix.csv"), row.names = TRUE)
write.csv(ko_path_info, file.path(OUT_DIR, "step7_ko_pathway_metadata.csv"), row.names = FALSE)

ko_pathway_score <- t(scale(t(ko_pathway_mat)))
ko_pathway_score[!is.finite(ko_pathway_score)] <- 0
write.csv(ko_pathway_score, file.path(OUT_DIR, "step7_ko_pathway_score_matrix_zscore.csv"), row.names = TRUE)

ko_path_diff_rows <- lapply(seq_len(nrow(ko_pathway_score)), function(i) {
  cc <- corr_with_p(ko_pathway_score[i, crc_idx], diff_num[crc_idx], method = "spearman")
  data.frame(
    pathway_id = ko_pathway_ids[i],
    pathway_name = ko_path_info$pathway_name[i],
    level1 = ko_path_info$level1[i],
    level2 = ko_path_info$level2[i],
    rho = unname(cc["rho"]),
    pval = unname(cc["p"]),
    n = unname(cc["n"]),
    stringsAsFactors = FALSE
  )
})
ko_path_diff <- do.call(rbind, ko_path_diff_rows)
ko_path_diff$qval <- p.adjust(ko_path_diff$pval, method = "BH")
ko_path_diff <- ko_path_diff[order(ko_path_diff$pval), , drop = FALSE]
write.csv(ko_path_diff, file.path(OUT_DIR, "step7_ko_pathway_vs_diff_spearman.csv"), row.names = FALSE)

ko_path_pick <- ko_path_diff[is.finite(ko_path_diff$pval), , drop = FALSE]
ko_path_pick <- filter_by_level2(ko_path_pick, level2_col = "level2", fixed_list = LEVEL2_FIXED_LIST)
ko_path_pick <- ko_path_pick[order(ko_path_pick$qval, ko_path_pick$pval), , drop = FALSE]
ko_path_pick <- ko_path_pick[seq_len(min(TOP_N_KO_PATHWAY, nrow(ko_path_pick))), , drop = FALSE]

ko_path_sel <- ko_path_pick$pathway_id
ko_path_idx <- match(ko_path_sel, ko_pathway_ids)
ko_path_score_sel <- ko_pathway_score[ko_path_idx, , drop = FALSE]

if (nrow(path_score_sel) == 0 || nrow(ko_path_score_sel) == 0) {
  stop("No pathway rows selected for metabolite-vs-KO pathway chessboard.")
}

rho_path2path <- matrix(NA_real_, nrow = nrow(path_score_sel), ncol = nrow(ko_path_score_sel))
p_path2path <- matrix(NA_real_, nrow = nrow(path_score_sel), ncol = nrow(ko_path_score_sel))
rownames(rho_path2path) <- rownames(path_score_sel)
colnames(rho_path2path) <- rownames(ko_path_score_sel)
rownames(p_path2path) <- rownames(path_score_sel)
colnames(p_path2path) <- rownames(ko_path_score_sel)

for (i in seq_len(nrow(path_score_sel))) {
  for (j in seq_len(nrow(ko_path_score_sel))) {
    cc <- corr_with_p(as.numeric(path_score_sel[i, ]), as.numeric(ko_path_score_sel[j, ]), method = "spearman")
    rho_path2path[i, j] <- unname(cc["rho"])
    p_path2path[i, j] <- unname(cc["p"])
  }
}

cross_path_corr <- data.frame(
  met_pathway = rep(rownames(rho_path2path), times = ncol(rho_path2path)),
  ko_pathway = rep(colnames(rho_path2path), each = nrow(rho_path2path)),
  rho = as.vector(rho_path2path),
  pval = as.vector(p_path2path),
  stringsAsFactors = FALSE
)
cross_path_corr$qval <- p.adjust(cross_path_corr$pval, method = "BH")
write.csv(cross_path_corr, file.path(OUT_DIR, "step7_met_pathway_vs_ko_pathway_correlations.csv"), row.names = FALSE)

met_path_group <- path_info[match(path_sel, path_info$pathway_id), c("pathway_id", "level2"), drop = FALSE]
met_path_group$level2[is.na(met_path_group$level2) | met_path_group$level2 == ""] <- "Unclassified"

ko_path_group <- ko_path_info[match(ko_path_sel, ko_path_info$pathway_id), c("pathway_id", "level2"), drop = FALSE]
ko_path_group$level2[is.na(ko_path_group$level2) | ko_path_group$level2 == ""] <- "Unclassified"

sig_path2path <- matrix("", nrow = nrow(p_path2path), ncol = ncol(p_path2path))
sig_path2path[p_path2path < 0.001] <- "***"
sig_path2path[p_path2path >= 0.001 & p_path2path < 0.01] <- "**"
sig_path2path[p_path2path >= 0.01 & p_path2path < 0.05] <- "*"

met_levels <- unique(met_path_group$level2)
ko_path_levels <- unique(ko_path_group$level2)

met_palette <- structure(
  hcl(h = seq(20, 380, length.out = length(met_levels) + 1)[-1], c = 65, l = 55),
  names = met_levels
)
ko_path_palette <- structure(
  hcl(h = seq(35, 395, length.out = length(ko_path_levels) + 1)[-1], c = 85, l = 68),
  names = ko_path_levels
)

ha_top_path <- HeatmapAnnotation(
  KO_Pathway_Group = factor(ko_path_group$level2, levels = ko_path_levels),
  col = list(KO_Pathway_Group = ko_path_palette),
  annotation_name_gp = gpar(fontsize = 9, fontface = "bold"),
  show_legend = TRUE
)

ha_left_path <- rowAnnotation(
  Met_Pathway_Group = factor(met_path_group$level2, levels = met_levels),
  col = list(Met_Pathway_Group = met_palette),
  annotation_name_gp = gpar(fontsize = 9, fontface = "bold"),
  show_legend = TRUE
)

ht_path2path <- Heatmap(
  rho_path2path,
  name = "Spearman r",
  col = col_fun,
  top_annotation = ha_top_path,
  left_annotation = ha_left_path,
  row_split = factor(met_path_group$level2, levels = met_levels),
  column_split = factor(ko_path_group$level2, levels = ko_path_levels),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_dend = FALSE,
  show_column_dend = FALSE,
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 7),
  column_names_rot = 45,
  rect_gp = gpar(col = "white", lwd = 0.6),
  border = TRUE,
  cell_fun = function(j, i, x, y, width, height, fill) {
    if (sig_path2path[i, j] != "") {
      grid.text(sig_path2path[i, j], x, y, gp = gpar(fontsize = 6, col = "black"))
    }
  }
)

png(file.path(OUT_DIR, "met_pathway_vs_ko_pathway_chessboard_heatmap.png"), width = 3600, height = 2200, res = 300)
draw(
  ht_path2path,
  column_title = "Metabolite Pathway vs KO Pathway Correlation Chessboard",
  column_title_gp = gpar(fontsize = 14, fontface = "bold"),
  heatmap_legend_side = "right",
  annotation_legend_side = "right",
  merge_legend = FALSE
)
dev.off()

pdf(file.path(OUT_DIR, "met_pathway_vs_ko_pathway_chessboard_heatmap.pdf"), width = 16, height = 10)
draw(
  ht_path2path,
  column_title = "Metabolite Pathway vs KO Pathway Correlation Chessboard",
  column_title_gp = gpar(fontsize = 14, fontface = "bold"),
  heatmap_legend_side = "right",
  annotation_legend_side = "right",
  merge_legend = FALSE
)
dev.off()

write.csv(met_path_group, file.path(OUT_DIR, "met_pathway_group_labels.csv"), row.names = FALSE)
write.csv(ko_path_group, file.path(OUT_DIR, "ko_pathway_group_labels.csv"), row.names = FALSE)

summary_txt <- c(
  sprintf("Pathways in matrix: %d", nrow(pathway_mat)),
  sprintf("KOs used for correlation: %d", ncol(ko_mat)),
  sprintf("Pathways used for chessboard: %d", nrow(path_score_sel)),
  sprintf("Significant pathway-vs-Diff (q < %.2f): %d", Q_CUTOFF, sum(path_diff$qval < Q_CUTOFF, na.rm = TRUE)),
  sprintf("Significant pathway-vs-KO cells (p < 0.05): %d", sum(p_mat < 0.05, na.rm = TRUE)),
  sprintf("KO pathways in matrix: %d", nrow(ko_pathway_mat)),
  sprintf("KO pathways used for chessboard: %d", nrow(ko_path_score_sel)),
  sprintf("Significant met-pathway vs KO-pathway cells (p < 0.05): %d", sum(p_path2path < 0.05, na.rm = TRUE))
)
writeLines(summary_txt, con = file.path(OUT_DIR, "run_summary.txt"))

message("Done. Output directory: ", OUT_DIR)
