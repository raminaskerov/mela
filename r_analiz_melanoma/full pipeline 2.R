# =============================================================================
# INTEGRATIVE ANALYSIS: ENCORAFENIB RESISTANCE IN MALIGNANT MELANOMA
# GSE283251 | Colakoglu Bergel et al., Scientific Reports 2025
#
# GOALS:
#   Stage 0  — Data loading & preprocessing
#   Stage 1  — Paper reproduction: NOISeq DEGs + miRNAs
#   Stage 2  — Improved visualizations (volcano, bubble, annotated heatmap)
#   Priority 1 — miRNA–mRNA integration network
#   Priority 2 — Co-expression / WGCNA-style module analysis
#   Priority 3 — TF enrichment analysis
#
# SAMPLES (n=1 per condition, no replicates):
#   A375-RK  → A375_RC  : Resistant control (no drug)
#   A375-R-10 → A375_R10 : Resistant + 10 nM Encorafenib
#   A375-S-K  → A375_SC  : Sensitive control (parental)
#   A375-S-10 → A375_S10 : Sensitive + 10 nM Encorafenib
#
# COMPARISONS (matching paper Fig 4d):
#   Comp1: A375_SC vs A375_RC  ← primary resistance signal
#   Comp2: A375_SC vs A375_S10 ← drug effect in sensitive
#   Comp3: A375_RC vs A375_R10 ← drug effect in resistant
#
# NOTES:
#   - Files are UTF-16LE encoded (GEO submission format)
#   - WGCNA with n=4 is exploratory; interpretation caveats noted inline
#   - multiMiR queries external databases; internet connection required at runtime
# =============================================================================


# =============================================================================
# SECTION 0: SETUP & DATA LOADING
# =============================================================================

# ── 0.1  Install packages (run once) ─────────────────────────────────────────
setwd("C:/bioinformatics")
a <- "na"
print(a)
cran_pkgs <- c(
  "readr",          # file reading with encoding control
  "dplyr",          # data wrangling
  "tidyr",          # reshaping
  "ggplot2",        # plotting
  "ggrepel",        # non-overlapping labels
  "pheatmap",       # heatmaps
  "RColorBrewer",   # color palettes
  "igraph",         # network construction
  "ggraph",         # network visualization
  "tidygraph",      # tidy network manipulation
  "enrichR",        # TF enrichment (ChEA3, ENCODE)
  "WGCNA",          # weighted co-expression
  "flashClust",     # fast hierarchical clustering for #WGCNA aktıv deıl
  "scales",         # axis formatting
  "patchwork",      # combine ggplots
  "viridis",        # color-blind-safe palettes
  "stringr"         # string manipulation
)
for (p in cran_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}

# ── 0.2  Load libraries ───────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(readr);      library(dplyr);     library(tidyr)
  library(ggplot2);    library(ggrepel);   library(pheatmap)
  library(RColorBrewer); library(scales);  library(patchwork)
  library(viridis);    library(stringr)
  library(NOISeq)
  library(igraph);     library(ggraph);    library(tidygraph)
  library(enrichR)
  library(WGCNA)
  library(multiMiR)
  library(org.Hs.eg.db)
})

# WGCNA: allow multi-threading if available
enableWGCNAThreads()

# ── 0.3  Set paths ────────────────────────────────────────────────────────────

DATA_DIR   <- "C:/bioinformatics" # directory containing the .gz files
OUTPUT_DIR <- "output"     # all plots and tables go here
dir.create(OUTPUT_DIR, showWarnings = FALSE)

# =============================================================================
# FIX v2 — readr-based loader (CONFIRMED WORKING after correct DATA_DIR)
#
# This version uses readr::read_tsv() with locale(encoding="UTF-16LE").
# It works once DATA_DIR points directly to the folder containing the .gz files
# (no spaces or special characters in path).
#
# Applied fixes vs original pipeline:
#   1. readr instead of rawToChar() → handles large UTF-16LE gzipped files
#   2. as.data.frame() before NA assignment → tibble doesn't support df[is.na]
#   3. drop=FALSE on [,-1] → guards against single-column edge case
# =============================================================================


# ── 0.4  Helper: read UTF-16LE gzipped count files ───────────────────────────

read_count_file <- function(path) {
  readr::read_tsv(
    path,
    locale         = readr::locale(encoding = "UTF-16LE"),
    show_col_types = FALSE,
    name_repair    = "minimal",
    progress       = FALSE
  )
}


# ── 0.5  Load all mRNA count files ───────────────────────────────────────────


mrna_files <- list(
  A375_RC  = file.path(DATA_DIR, "GSM8658654_A375-RK_mRNA_count.txt.gz"),
  A375_R10 = file.path(DATA_DIR, "GSM8658655_A375-R-10_mRNA_count.txt.gz"),
  A375_SC  = file.path(DATA_DIR, "GSM8658656_A375-S-K_mRNA_count.txt.gz"),
  A375_S10 = file.path(DATA_DIR, "GSM8658657_A375-S-10_mRNA_count.txt.gz")
)

mrna_raw <- list()
for (nm in names(mrna_files)) {
  mrna_raw[[nm]] <- tryCatch(
    read_count_file(mrna_files[[nm]]),
    error = function(e) stop(sprintf("Failed to read %s: %s", nm, e$message))
  )
  message(sprintf("  Loaded %-10s: %d rows x %d cols | first cols: %s",
                  nm,
                  nrow(mrna_raw[[nm]]),
                  ncol(mrna_raw[[nm]]),
                  paste(head(colnames(mrna_raw[[nm]]), 4), collapse = ", ")))
}


# ── extract_counts: Gene_Symbol → count, protein-coding only ─────────────────

extract_counts <- function(df, sample_name) {
  count_col <- grep("Read_Count", colnames(df), value = TRUE)[1]
  
  df |>
    dplyr::select(Gene_Symbol, gene_biotype, all_of(count_col)) |>
    dplyr::rename(count = all_of(count_col)) |>
    dplyr::filter(
      gene_biotype == "protein_coding",
      !is.na(Gene_Symbol),
      Gene_Symbol  != "",
      Gene_Symbol  != "."
    ) |>
    dplyr::group_by(Gene_Symbol) |>
    dplyr::summarise(count = sum(count, na.rm = TRUE), .groups = "drop") |>
    dplyr::rename(!!sample_name := count)
}

mrna_counts_list <- mapply(extract_counts,
                           mrna_raw, names(mrna_raw),
 
                                                     SIMPLIFY = FALSE)
# =  == = ==  =  = = =  == ====== ==  = = == ==   == = = = = = = == = ===    ==   ================================
# protein coding only
# ======= = = = = = = =  = = = = =  = = = = = = = = =  = = = = = = = ==========================================
# Report row counts before joining (useful for diagnosing 0-row issues)
for (nm in names(mrna_counts_list)) {
  message(sprintf("  extract_counts(%s): %d protein-coding genes",
                  nm, nrow(mrna_counts_list[[nm]])))
}

mrna_counts <- Reduce(
  function(x, y) dplyr::full_join(x, y, by = "Gene_Symbol"),
  mrna_counts_list
)

# FIX: tibble → base data.frame before NA assignment
mrna_counts                <- as.data.frame(mrna_counts)
mrna_counts[is.na(mrna_counts)] <- 0L
mrna_mat                   <- as.matrix(mrna_counts[, -1, drop = FALSE])
rownames(mrna_mat)         <- mrna_counts$Gene_Symbol
storage.mode(mrna_mat)     <- "integer"

# ── 0.6  Load all miRNA count files ──────────────────────────────────────────


mirna_files <- list(
  A375_RC  = file.path(DATA_DIR, "GSM8658654_A375-RK_miRNA_count.txt.gz"),
  A375_R10 = file.path(DATA_DIR, "GSM8658655_A375-R-10_miRNA_count.txt.gz"),
  A375_SC  = file.path(DATA_DIR, "GSM8658656_A375-S-K_miRNA_count.txt.gz"),
  A375_S10 = file.path(DATA_DIR, "GSM8658657_A375-S-10_miRNA_count.txt.gz")
)

mirna_raw <- list()
for (nm in names(mirna_files)) {
  mirna_raw[[nm]] <- tryCatch(
    read_count_file(mirna_files[[nm]]),
    error = function(e) stop(sprintf("Failed to read %s: %s", nm, e$message))
  )
  message(sprintf("  Loaded %-10s: %d miRNAs", nm, nrow(mirna_raw[[nm]])))
}

extract_mirna_counts <- function(df, sample_name) {
  count_col <- grep("Read_Count", colnames(df), value = TRUE)[1]
  df |>
    dplyr::select(Mature_ID, all_of(count_col)) |>
    dplyr::rename(!!sample_name := all_of(count_col)) |>
    dplyr::filter(!is.na(Mature_ID), Mature_ID != "")
}

mirna_counts <- Reduce(
  function(x, y) dplyr::full_join(x, y, by = "Mature_ID"),
  mapply(extract_mirna_counts, mirna_raw, names(mirna_raw), SIMPLIFY = FALSE)
)

# FIX: tibble → base data.frame before NA assignment
mirna_counts                <- as.data.frame(mirna_counts)
mirna_counts[is.na(mirna_counts)] <- 0L
mirna_mat                   <- as.matrix(mirna_counts[, -1, drop = FALSE])
rownames(mirna_mat)         <- mirna_counts$Mature_ID
storage.mode(mirna_mat)     <- "integer"


# ── 0.7  Normalization ────────────────────────────────────────────────────────

tmm_normalize <- function(count_mat) {
  lib_sizes <- colSums(count_mat)
  ref_lib   <- exp(mean(log(lib_sizes[lib_sizes > 0])))
  scale_fac <- lib_sizes / ref_lib
  sweep(count_mat, 2, scale_fac, FUN = "/")
}

mrna_norm <- tmm_normalize(mrna_mat)
mirna_cpm <- sweep(mirna_mat, 2, colSums(mirna_mat), FUN = "/") * 1e6

# =============================================================================
# SECTION 1: NOISeq DEG ANALYSIS (PAPER REPRODUCTION)
# =============================================================================

# NOISeq is designed for single-sample (no-replicate) differential expression.
# It uses a noise distribution estimated from within-sample variability.

# ── 1.1  Build NOISeq data object ────────────────────────────────────────────

# Sample metadata
sample_meta <- data.frame(
  condition = c("Resistant", "Resistant", "Sensitive", "Sensitive"),
  treatment = c("Control", "Drug10nM", "Control", "Drug10nM"),
  row.names = c("A375_RC", "A375_R10", "A375_SC", "A375_S10")
)

# Filter low-count genes: keep genes with CPM > 0.5 in at least 1 sample
mrna_cpm_check <- sweep(mrna_mat, 2, colSums(mrna_mat) / 1e6, FUN = "/")
keep_genes <- rowSums(mrna_cpm_check > 0.5) >= 1
mrna_filtered <- mrna_mat[keep_genes, ]
message(sprintf("After CPM filtering: %d genes retained", nrow(mrna_filtered)))
# ── Histone gene filter ───────────────────────────────────────────────────────
# Removes replication-dependent canonical histones (multi-mapping artifacts)
# Run this BEFORE NOISeq, applied to mrna_filtered

histone_pattern <- "^HIST[0-9]|^H[1-4][A-Z]|^HIST|^H2A|^H2B|^H3[^F]|^H4C"

histone_genes <- grep(histone_pattern, rownames(mrna_filtered), value = TRUE)
message(sprintf("Histone genes identified: %d", length(histone_genes)))
print(histone_genes)  # review before removing

mrna_filtered <- mrna_filtered[!rownames(mrna_filtered) %in% histone_genes, ]
message(sprintf("Genes after histone removal: %d", nrow(mrna_filtered)))

# Protect legitimate genes if they appear in histone_genes
protect <- c("H3F3A", "H3F3B", "HMOX1")
histone_genes <- histone_genes[!histone_genes %in% protect]
mrna_filtered <- mrna_filtered[!rownames(mrna_filtered) %in% histone_genes, ]

run_noiseq <- function(count_mat, meta, samples_A, samples_B,
                       label, q = 0.8) {
  sel_cols  <- c(samples_A, samples_B)
  sub_mat   <- count_mat[, sel_cols, drop = FALSE]
  sub_meta  <- data.frame(
    group = c(rep("A", length(samples_A)), rep("B", length(samples_B))),
    row.names = sel_cols
  )
  keep    <- rowSums(sub_mat) > 0
  sub_mat <- sub_mat[keep, ]
  
  noiseq_obj <- NOISeq::readData(data = sub_mat, factors = sub_meta)
  res <- NOISeq::noiseq(
    input      = noiseq_obj,
    factor     = "group",
    conditions = c("A", "B"),
    norm       = "tmm",
    replicates = "no",
    k          = 0.5
  )
  
  full_results                 <- NOISeq::degenes(res, q = 0, M = NULL)
  
  degs <- full_results[full_results$prob >= q, ]
  full_results$M <- full_results$M*-1
  degs$M <- degs$M*-1
  full_results$ranking <- full_results$ranking*-1
  degs$ranking <- degs$ranking*-1
  message(sprintf("  [%s] DEGs (prob>=%.2f): %d",
                  label, q, nrow(degs)))
  
  list(result   = res,
       degs     = degs,
       full     = full_results,
       label    = label,
       samples_A = samples_A,
       samples_B = samples_B)
}

comp1 <- run_noiseq(
  mrna_filtered,
  meta = sample_meta,
  samples_A = "A375_SC",
  samples_B = "A375_RC",
  label = "SC_vs_RC"
)

comp2 <- run_noiseq(
  mrna_filtered,
  meta = sample_meta,
  samples_A = "A375_SC",
  samples_B = "A375_S10",
  label = "SC_vs_S10"
)

comp3 <- run_noiseq(
  mrna_filtered,
  meta = sample_meta,
  samples_A = "A375_RC",
  samples_B = "A375_R10",
  label = "RC_vs_R10"
)

comp4 <- run_noiseq(
  mrna_filtered,
  meta = sample_meta,
  samples_A = "A375_S10",
  samples_B = "A375_R10",
  label = "S10_vs_R10"
)

save_deg_table <- function(comp_obj, filename) {
  # 1. Target the specific dataframe inside the list
  df <- as.data.frame(comp_obj$degs)
  
  # 2. Add rownames as a NEW column (so they aren't lost when row.names=FALSE)
  df$gene <- rownames(df)
  
  # 3. Rename the columns
  df <- df |> dplyr::rename(log2FC = M, meanA = A_mean, meanB = B_mean)
  
  # 4. Add the direction column
  #df$direction <- ifelse(df$log2FC > 0, "Up_in_B", "Down_in_B")
  
  # 5. Save the file
  write.csv(df, file = file.path(OUTPUT_DIR, filename), row.names = FALSE)
  
  message(sprintf("  Saved: %s (%d DEGs)", filename, nrow(df)))
  invisible(df)
}

deg1 <- save_deg_table(comp1, "DEGs_SC_vs_RC_.csv")
deg2 <- save_deg_table(comp2, "DEGs_SC_vs_S10_.csv")
deg3 <- save_deg_table(comp3, "DEGs_RC_vs_R10_.csv")
deg4 <- save_deg_table(comp4, "DEGs_S10_vs_R10_.csv")


# ── 1.3  Run NOISeq for miRNAs ───────────────────────────────────────────────

#––––––––––––––––––––Mirna noiseq––––––––––––––––––––––––#

mirna_cpm <- sweep(mirna_mat, 2, colSums(mirna_mat), FUN = "/") * 1e6

sample_meta <- data.frame(
  condition = c("Resistant", "Resistant", "Sensitive", "Sensitive"),
  treatment = c("Control", "Drug10nM", "Control", "Drug10nM"),
  row.names = c("A375_RC", "A375_R10", "A375_SC", "A375_S10")
)

run_noiseq_comparison <- function(count_mat, meta, samples_A, samples_B,
                                  label, q = 0.8) {
  sel_cols  <- c(samples_A, samples_B)
  sub_mat   <- count_mat[, sel_cols, drop = FALSE]
  sub_meta  <- data.frame(
    group = c(rep("A", length(samples_A)), rep("B", length(samples_B))),
    row.names = sel_cols
  )
  keep    <- rowSums(sub_mat) > 0
  sub_mat <- sub_mat[keep, ]
  
  noiseq_obj <- NOISeq::readData(data = sub_mat, factors = sub_meta)
  res <- NOISeq::noiseq(
    input      = noiseq_obj,
    factor     = "group",
    conditions = c("A", "B"),
    norm       = "tmm",
    replicates = "no",
    k          = 0.5
  )
  
  
  full_results                 <- NOISeq::degenes(res, q = 0, M = NULL)
  full_results$M <- full_results$M*-1
  
  full_results$ranking <- full_results$ranking*-1
  scores <- scale(abs(full_results$ranking)) * full_results$prob^1
  full_results$priority_score <- as.numeric(scores)
  
  threshold     <- quantile(full_results$priority_score, 0.40)
  degs_priority <- full_results[full_results$priority_score >= threshold &
                                  full_results$prob >= q, ]
  
  message(sprintf("  [%s] DEGs (prob>=%.2f): %d",
                  label, q, nrow(degs_priority)))
  
  list(result   = res,
       degs     = degs_priority,
       full     = full_results,
       label    = label,
       samples_A = samples_A,
       samples_B = samples_B)
}



# ── 1.3  Run NOISeq for miRNAs ───────────────────────────────────────────────

# Filter low-count miRNAs
mirna_keep <- rowSums(mirna_mat > 2) >= 1
mirna_filtered <- mirna_mat[mirna_keep, ]
message(sprintf("miRNAs after filtering: %d", nrow(mirna_filtered)))

mirna_comp1 <- run_noiseq_comparison(
  mirna_filtered,
  meta = sample_meta,
  samples_A = "A375_SC",
  samples_B = "A375_RC",
  label = "miRNA_SC_vs_RC"
)

mirna_comp2 <- run_noiseq_comparison(
  mirna_filtered,
  meta = sample_meta,
  samples_A = "A375_SC",
  samples_B = "A375_S10",
  label = "miRNA_SC_vs_S10"
)

mirna_comp3 <- run_noiseq_comparison(
  mirna_filtered,
  meta = sample_meta,
  samples_A = "A375_RC",
  samples_B = "A375_R10",
  label = "miRNA_RC_vs_R10"
)

mirna_comp4 <- run_noiseq_comparison(
  mirna_filtered,
  meta = sample_meta,
  samples_A = "A375_S10",
  samples_B = "A375_R10",
  label = "miRNA_S10_vs_R10"
)
mirna_deg1 <- save_deg_table(mirna_comp1, "DEGs_miRNA_SC_vs_RC_.csv")
mirna_deg2 <- save_deg_table(mirna_comp2, "DEGs_miRNA_SC_vs_S10_.csv")
mirna_deg3 <- save_deg_table(mirna_comp3, "DEGs_miRNA_RC_vs_R10_.csv")
mirna_deg4 <- save_deg_table(mirna_comp4, "DEGs_miRNA_S10_vs_R10_.csv")
# =============================================================================
# SECTION 2: IMPROVED VISUALIZATIONS
# =============================================================================

# ── Gene sets for annotation ──────────────────────────────────────────────────

IRON_GENES <- c("NCOA4", "FTH1", "TFRC", "SLC7A11", "GPX4", "IREB2",
                "FTL", "HAMP", "SLC40A1", "HMOX1", "CYBRD1", "STEAP3")
AUTOPHAGY_GENES <- c("BECN1", "MAP1LC3B", "ATG5", "ATG7", "ATG12",
                     "SQSTM1", "ULK1", "WIPI2", "LAMP2", "RAB7A")
MAPK_AKT_GENES  <- c("BRAF", "KRAS", "NRAS", "MAP2K1", "MAPK1", "MAPK3",
                     "AKT1", "AKT2", "PIK3CA", "PIK3R1", "PTEN",
                     "DUSP6", "DUSP4", "RAF1", "CRAF")
MDR_GENES        <- c("ABCB1", "YBX1", "MDR1")

GENE_CATEGORIES <- c(
  setNames(rep("Iron/Ferroptosis", length(IRON_GENES)), IRON_GENES),
  setNames(rep("Autophagy",        length(AUTOPHAGY_GENES)), AUTOPHAGY_GENES),
  setNames(rep("MAPK/AKT",         length(MAPK_AKT_GENES)), MAPK_AKT_GENES),
  setNames(rep("MDR",              length(MDR_GENES)), MDR_GENES)
)

CATEGORY_COLORS <- c(
  "Iron/Ferroptosis" = "#E64B35",
  "Autophagy"        = "#4DBBD5",
  "MAPK/AKT"         = "#00A087",
  "MDR"              = "#F39B7F",
  "Other DEG"        = "#8491B4",
  "NS"               = "grey75"
)

# ── 2.1  Volcano plots (one per comparison) ───────────────────────────────────

make_volcano <- function(deg_df, title, label_genes = NULL,
                         fc_threshold = 1, prob_threshold = 0.8) {
  
  df <- as.data.frame(deg_df$result@results[[1]])
  df$gene <- rownames(df)
  colnames(df)[colnames(df) == "M"]    <- "log2FC"
  colnames(df)[colnames(df) == "prob"] <- "prob"
  
  # Significance classification
  df$sig <- "NS"
  df$sig[df$prob >= prob_threshold & df$log2FC >  fc_threshold] <- "Up"
  df$sig[df$prob >= prob_threshold & df$log2FC < -fc_threshold] <- "Down"
  
  # Annotate gene categories
  df$category <- dplyr::case_when(
    df$gene %in% IRON_GENES      ~ "Iron/Ferroptosis",
    df$gene %in% AUTOPHAGY_GENES ~ "Autophagy",
    df$gene %in% MAPK_AKT_GENES  ~ "MAPK/AKT",
    df$gene %in% MDR_GENES       ~ "MDR",
    df$sig != "NS"               ~ "Other DEG",
    TRUE                         ~ "NS"
  )
  
  # Genes to label: paper key genes + top DEGs by |FC|
  key_label <- c(IRON_GENES, AUTOPHAGY_GENES, MAPK_AKT_GENES, MDR_GENES)
  if (!is.null(label_genes)) key_label <- c(key_label, label_genes)
  top_degs <- df |>
    dplyr::filter(sig != "NS") |>
    dplyr::arrange(desc(abs(log2FC))) |>
    dplyr::slice_head(n = 20) |>
    dplyr::pull(gene)
  df$label <- ifelse(df$gene %in% c(key_label, top_degs) & df$sig != "NS",
                     df$gene, NA)
  
  n_up   <- sum(df$sig == "Up",   na.rm = TRUE)
  n_down <- sum(df$sig == "Down", na.rm = TRUE)
  
  ggplot(df, aes(x = log2FC, y = prob, color = category, label = label)) +
    geom_point(alpha = 0.5, size = 1.2) +
    geom_point(data = df[!is.na(df$label), ],
               aes(color = category), size = 2.5, alpha = 0.9) +
    geom_hline(yintercept = prob_threshold, linetype = "dashed",
               color = "grey40", linewidth = 0.4) +
    geom_vline(xintercept = c(-fc_threshold, fc_threshold),
               linetype = "dashed", color = "grey40", linewidth = 0.4) +
    ggrepel::geom_text_repel(
      na.rm = TRUE, size = 2.8, max.overlaps = 20,
      segment.size = 0.2, segment.color = "grey50",
      fontface = "italic"
    ) +
    scale_color_manual(values = CATEGORY_COLORS, drop = FALSE) +
    scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
    labs(
      title    = title,
      subtitle = sprintf("↑ %d up  |  ↓ %d down  (prob ≥ %.2f, |log2FC| > %g)",
                         n_up, n_down, prob_threshold, fc_threshold),
      x        = "log₂ Fold Change (B vs A)",
      y        = "NOISeq Probability",
      color    = "Gene Category"
    ) +
    theme_bw(base_size = 11) +
    theme(
      legend.position  = "right",
      plot.title       = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )
}

p_vol1 <- make_volcano(comp1, "Resistant vs Sensitive (Control)\nSC → RC")
p_vol2 <- make_volcano(comp2, "Drug Effect in Sensitive\nSC → S10 (10 nM)")
p_vol3 <- make_volcano(comp3, "Drug Effect in Resistant\nRC → R10 (10 nM)")
p_vol4 <- make_volcano(comp4, "Drug Effect Sensitive vs Resistant\nS10 → R10 (10 nM)")

ggsave(file.path(OUTPUT_DIR, "volcano_SC_vs_RC.pdf"),  p_vol1, width = 9, height = 7)
ggsave(file.path(OUTPUT_DIR, "volcano_SC_vs_S10.pdf"), p_vol2, width = 9, height = 7)
ggsave(file.path(OUTPUT_DIR, "volcano_RC_vs_R10.pdf"), p_vol3, width = 9, height = 7)
ggsave(file.path(OUTPUT_DIR, "volcano_S10_vs_R10.pdf"), p_vol4, width = 9, height = 7)
message("  Volcano plots saved.")

# ── 2.2  Annotated heatmap: key pathway genes ─────────────────────────────────

make_pathway_heatmap <- function(count_mat, gene_sets, title, filename) {
  all_genes <- unlist(gene_sets)
  present   <- all_genes[all_genes %in% rownames(count_mat)]
  
  # log2(CPM+1) normalization for display
  cpm_mat  <- sweep(count_mat[present, ], 2,
                    colSums(count_mat) / 1e6, FUN = "/")
  log_mat  <- log2(cpm_mat + 1)
  
  # Z-score per gene across samples
  z_mat <- t(scale(t(log_mat)))
  z_mat[is.nan(z_mat)] <- 0
  
  # Row annotations: gene category
  row_annot <- data.frame(
    Category = sapply(present, function(g) {
      for (nm in names(gene_sets)) if (g %in% gene_sets[[nm]]) return(nm)
      return("Other")
    }),
    row.names = present
  )
  
  annot_colors <- list(
    Category = c(
      "Iron/Ferroptosis" = "#E64B35",
      "Autophagy"        = "#4DBBD5",
      "MAPK/AKT"         = "#00A087",
      "MDR"              = "#F39B7F"
    )
  )
  
  col_annot <- data.frame(
    Condition  = c("Resistant", "Resistant", "Sensitive", "Sensitive"),
    Treatment  = c("Control",   "Drug",      "Control",   "Drug"),
    row.names  = colnames(count_mat)
  )
  annot_colors$Condition <- c(Resistant = "#CC0000", Sensitive = "#0066CC")
  annot_colors$Treatment <- c(Control = "#999999", Drug = "#FF9900")
  
  pheatmap::pheatmap(
    z_mat,
    color            = colorRampPalette(c("#2166AC","white","#B2182B"))(100),
    annotation_row   = row_annot,
    annotation_col   = col_annot,
    annotation_colors= annot_colors,
    cluster_cols     = FALSE,
    cluster_rows     = TRUE,
    show_rownames    = TRUE,
    fontsize_row     = 8,
    fontsize_col     = 9,
    border_color     = NA,
    main             = title,
    filename         = file.path(OUTPUT_DIR, filename),
    width = 8, height = 10
  )
  message(sprintf("  Heatmap saved: %s (%d genes)", filename, nrow(z_mat)))
}

gene_sets_heatmap <- list(
  "Iron/Ferroptosis" = IRON_GENES,
  "Autophagy"        = AUTOPHAGY_GENES,
  "MAPK/AKT"         = MAPK_AKT_GENES,
  "MDR"              = MDR_GENES
)

make_pathway_heatmap(
  mrna_mat,
  gene_sets_heatmap,
  title    = "Key Pathway Genes Across All Conditions\n(Z-scored log2CPM)",
  filename = "heatmap_key_genes.pdf"
)

# ── 2.3  Iron gene expression barplot (replicates Fig 6a conceptually) ────────

iron_genes_plot <- c("FTH1","NCOA4","TFRC","SLC7A11","GPX4","IREB2",
                     "FTL","HMOX1","STEAP3")
iron_present <- iron_genes_plot[iron_genes_plot %in% rownames(mrna_mat)]

iron_plot_df <- as.data.frame(mrna_mat[iron_present, ]) |>
  tibble::rownames_to_column("gene") |>
  tidyr::pivot_longer(-gene, names_to = "sample", values_to = "count") |>
  dplyr::mutate(
    condition = dplyr::case_when(
      sample == "A375_RC"  ~ "Resistant\nControl",
      sample == "A375_R10" ~ "Resistant\n10 nM",
      sample == "A375_SC"  ~ "Sensitive\nControl",
      sample == "A375_S10" ~ "Sensitive\n10 nM"
    ),
    condition = factor(condition,
                       levels = c("Sensitive\nControl","Sensitive\n10 nM",
                                  "Resistant\nControl","Resistant\n10 nM")),
    log2_cpm  = log2(count / (sum(mrna_mat[, sample]) / 1e6) + 1)
  )

p_iron_bar <- ggplot(iron_plot_df,
                     aes(x = condition, y = log2_cpm, fill = condition)) +
  geom_col(width = 0.7, color = "white", linewidth = 0.3) +
  facet_wrap(~ gene, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = c(
    "Sensitive\nControl"  = "#6BAED6",
    "Sensitive\n10 nM"    = "#2171B5",
    "Resistant\nControl"  = "#FC8D59",
    "Resistant\n10 nM"    = "#D7301F"
  )) +
  labs(
    title = "Iron Metabolism Gene Expression",
    subtitle = "log₂(CPM + 1) across all conditions",
    x = NULL, y = "log₂(CPM + 1)", fill = "Condition"
  ) +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x     = element_text(size = 7),
    strip.background= element_rect(fill = "#f0f0f0"),
    strip.text      = element_text(face = "bold", size = 9),
    legend.position = "none",
    plot.title      = element_text(face = "bold")
  )

ggsave(file.path(OUTPUT_DIR, "iron_gene_barplot.pdf"),
       p_iron_bar, width = 10, height = 8)
message("  Iron gene barplot saved.")

# ── 2.4  miRNA heatmap ────────────────────────────────────────────────────────

# Show all detected miRNAs from the paper + miR-140-3p highlight
paper_mirna_ids <- c(
  "hsa-miR-331-3p", "hsa-let-7c-5p",  "hsa-miR-296-5p", "hsa-let-7b-5p",
  "hsa-miR-31-5p",  "hsa-miR-1260b",  "hsa-let-7a-5p",  "hsa-miR-1229-3p",
  "hsa-miR-6516-3p","hsa-miR-744-5p", "hsa-miR-181a-5p","hsa-miR-34a-5p",
  "hsa-miR-103a-3p","hsa-miR-140-3p"
)
mirna_present <- paper_mirna_ids[paper_mirna_ids %in% rownames(mirna_mat)]

mirna_heatmap_mat <- mirna_cpm[mirna_present, ]
mirna_z <- t(scale(t(log2(mirna_heatmap_mat + 1))))
mirna_z[is.nan(mirna_z)] <- 0

mirna_row_annot <- data.frame(
  Role = c(rep("RC_vs_SC_DEG", length(mirna_present))),
  row.names = mirna_present
)
mirna_row_annot["hsa-miR-140-3p", "Role"] <- "Iron-regulatory"

pheatmap::pheatmap(
  mirna_z,
  color          = colorRampPalette(c("#2166AC","white","#B2182B"))(100),
  annotation_col = data.frame(
    Condition = c("Resistant","Resistant","Sensitive","Sensitive"),
    Treatment = c("Control","Drug","Control","Drug"),
    row.names = colnames(mirna_z)
  ),
  annotation_row = mirna_row_annot,
  cluster_cols   = FALSE,
  cluster_rows   = TRUE,
  fontsize_row   = 9,
  main           = "DEG miRNAs (Paper Fig 5a)\nZ-scored log₂CPM",
  filename       = file.path(OUTPUT_DIR, "heatmap_mirna_paper_DEGs.pdf"),
  width = 7, height = 6
)


#–––––––––––––––––––––––ENRICHMENT–––––––––––––––––––––––––––––––#

gseacomparisons <- list(
  list(obj = comp1, label = "SC_vs_RC"),
  list(obj = comp2, label = "SC_vs_S10"),
  list(obj = comp3, label = "RC_vs_R10"),
  list(obj = comp4, label = "S10_vs_R10")
)

to_entrez <- function(symbols) {
  mapped <- clusterProfiler::bitr(
    symbols,
    fromType = "SYMBOL",
    toType   = "ENTREZID",
    OrgDb    = org.Hs.eg.db
  )
  message(sprintf("    Mapped %d / %d symbols to ENTREZID",
                  nrow(mapped), length(symbols)))
  mapped
}

# ── Helper: build ranked vector for GSEA ─────────────────────────────────────
# Uses full NOISeq results (all genes), ranked by log2FC (column M)

build_ranked_list <- function(full_df, mapped_df) {
  # full_df = comp$full which has rownames as genes, column M = log2FC
  scores <- full_df$M
  names(scores) <- rownames(full_df)
  
  # Convert to ENTREZID
  scores_df <- data.frame(
    SYMBOL   = names(scores),
    log2FC   = scores,
    stringsAsFactors = FALSE
  ) |>
    dplyr::inner_join(mapped_df, by = c("SYMBOL" = "SYMBOL")) |>
    dplyr::arrange(desc(log2FC))
  
  ranked <- scores_df$log2FC
  names(ranked) <- scores_df$ENTREZID
  # Remove duplicates (keep highest |FC| per ENTREZID)
  ranked <- ranked[!duplicated(names(ranked))]
  ranked
}

for (comp in gseacomparisons) {
  
  lbl  <- comp$label
  obj  <- comp$obj
  message(sprintf("\n=== %s ===", lbl))
  
  deg_genes <- rownames(obj$degs)
  message(sprintf("  DEGs: %d", length(deg_genes)))
  
  if (length(deg_genes) < 5) {
    message("  Too few DEGs — skipping.")
    next
  }
  
  deg_mapped <- tryCatch(
    to_entrez(deg_genes),
    error = function(e) { message("  bitr failed: ", e$message); NULL }
  )
  if (is.null(deg_mapped) || nrow(deg_mapped) == 0) next
  
  full_mapped <- tryCatch(
    to_entrez(rownames(obj$full)),
    error = function(e) NULL
  )
  
  if (!is.null(full_mapped)) {
    ranked_vec <- build_ranked_list(obj$full, full_mapped)
    
    
    
    # ── GSEA: KEGG ───────────────────────────────────────────────────────────
    
    message("  Running GSEA KEGG...")
    gsea_kegg <- tryCatch(
      clusterProfiler::gseKEGG(
        geneList      = ranked_vec,
        organism      = "hsa",
        minGSSize     = 10,
        maxGSSize     = 500,
        pAdjustMethod = "BH",
        pvalueCutoff  = 0.05,
        verbose       = FALSE
      ),
      error = function(e) { message("  GSEA KEGG failed: ", e$message); NULL }
    )
    
    if (!is.null(gsea_kegg) && nrow(as.data.frame(gsea_kegg)) > 0) {
      write.csv(as.data.frame(gsea_kegg),
                file.path(OUTPUT_DIR, sprintf("GSEA_KEGG_%s.csv", lbl)),
                row.names = FALSE)
      p_gkegg <- enrichplot::dotplot(gsea_kegg, showCategory = 20,
                                     split = ".sign") +
        facet_grid(. ~ .sign) +
        ggtitle(sprintf("GSEA KEGG | %s", lbl))
      ggsave(file.path(OUTPUT_DIR, sprintf("GSEA_KEGG_%s.pdf", lbl)),
             p_gkegg, width = 12, height = 8)
      assign(paste0("gsea_kegg_", lbl), gsea_kegg)
      message(sprintf("  GSEA KEGG: %d pathways", nrow(as.data.frame(gsea_kegg))))
    } else {
      message("  GSEA KEGG: no significant pathways.")
    }
    
    # ── GSEA: Reactome ───────────────────────────────────────────────────────
    
    message("  Running GSEA Reactome...")
    gsea_reactome <- tryCatch(
      ReactomePA::gsePathway(
        geneList      = ranked_vec,
        organism      = "human",
        minGSSize     = 10,
        maxGSSize     = 500,
        pAdjustMethod = "BH",
        pvalueCutoff  = 0.05,
        verbose       = FALSE
      ),
      error = function(e) { message("  GSEA Reactome failed: ", e$message); NULL }
    )
    
    if (!is.null(gsea_reactome) && nrow(as.data.frame(gsea_reactome)) > 0) {
      write.csv(as.data.frame(gsea_reactome),
                file.path(OUTPUT_DIR, sprintf("GSEA_Reactome_%s.csv", lbl)),
                row.names = FALSE)
      p_greactome <- enrichplot::dotplot(gsea_reactome, showCategory = 20,
                                         split = ".sign") +
        facet_grid(. ~ .sign) +
        ggtitle(sprintf("GSEA Reactome | %s", lbl))
      ggsave(file.path(OUTPUT_DIR, sprintf("GSEA_Reactome_%s.pdf", lbl)),
             p_greactome, width = 12, height = 9)
      assign(paste0("gsea_reactome_", lbl), gsea_reactome)
      message(sprintf("  GSEA Reactome: %d pathways", nrow(as.data.frame(gsea_reactome))))
    } else {
      message("  GSEA Reactome: no significant pathways.")
    }
  }
}

#=========================WEBGESTALT===================================#
comparisons_gsea <- list(
  list(comp=comp1, label = "SC_vs_RC",  A = "A375_SC", B = "A375_RC"),
  list(comp=comp2, label = "SC_vs_S10", A = "A375_SC", B = "A375_S10"),
  list(comp=comp3, label = "RC_vs_R10", A = "A375_RC", B = "A375_R10"),
  list(comp=comp4, label = "S10_vs_R10",A = "A375_S10",B = "A375_R10")
)
rm(comparisons_gsea)
for (c in comparisons_gsea) {
  ranked <- data.frame(
    gene  = rownames(c$comp$full),
    score = c$comp$full$M        # log2FC column from NOISeq full results
  ) |> dplyr::arrange(desc(score))
  
  write.table(ranked,
              file      = file.path(OUTPUT_DIR, 
                                    paste0("webgestalt_GSEA_", c$label, ".rnk")),
              sep       = "\t",
              row.names = FALSE,
              col.names = FALSE,
              quote     = FALSE)
}

# =============================================================================
# PRIORITY 1: miRNA–mRNA INTEGRATION NETWORK
# =============================================================================
# =============================================================================
# METHOD B: PATHWAY-SPECIFIC SUBNETWORKS — REVERSE multiMiR QUERY
# Logic: for each pathway, query multiMiR with pathway genes as TARGETS.
# This returns ALL miRNAs known to regulate those genes, then we filter
# down to your 47 DEG miRNAs (mirna_comp1) + anti-correlation check.
# Advantage over Method A: captures pathway-gene regulations that were
# missed when querying broadly from all DEG miRNAs.
#
# REQUIRED SESSION OBJECTS:
#   pathway_gene_df   — from pathway_subnetworks.R (or reloaded below)
#   mirna_comp1       — NOISeq result, SC vs RC miRNAs (47 DEGs)
#   mirna_cpm         — miRNA CPM matrix
#   mrna_mat          — mRNA count matrix
#   mrna_norm         — TMM-normalized mRNA matrix
#   comp1.gsea / comp1.gsea$full  — full NOISeq mRNA results SC vs RC
#                        (column M = log2FC)
# =============================================================================
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2)
  library(ggraph); library(igraph); library(tidygraph)
  library(ggrepel); library(multiMiR)
  library(org.Hs.eg.db); library(clusterProfiler)
  library(pheatmap); library(RColorBrewer);library(dynamicTreeCut)
})

OUTPUT_DIR <- "output"
dir.create(file.path(OUTPUT_DIR, "method_b"), showWarnings = FALSE)

# =============================================================================
# SECTION 1: PREPARE INPUTS
# =============================================================================

# ── 1.1  Reload pathway_gene_df if not in session ────────────────────────────

if (!exists("pathway_gene_df")) {
  f <- file.path(OUTPUT_DIR, "pathway_gene_assignments.csv")
  if (!file.exists(f)) stop("pathway_gene_assignments.csv not found. Run pathway_subnetworks.R first.")
  pathway_gene_df <- read.csv(f, stringsAsFactors = FALSE)
  message("  Reloaded pathway_gene_df from CSV.")
}

# ── 1.2  Extract 47 DEG miRNAs from mirna_comp1 ──────────────────────────────

#sadece 1. comp assign olub
deg_mirna_df <- as.data.frame(mirna_comp1$degs)
deg_mirna_df$mirna <- rownames(deg_mirna_df)

# Standardize column names (noiseq_adj variant)
if ("M" %in% colnames(deg_mirna_df))
  colnames(deg_mirna_df)[colnames(deg_mirna_df) == "M"] <- "log2FC"
if ("A_mean" %in% colnames(deg_mirna_df))
  colnames(deg_mirna_df)[colnames(deg_mirna_df) == "A_mean"] <- "mean_SC"
if ("B_mean" %in% colnames(deg_mirna_df))
  colnames(deg_mirna_df)[colnames(deg_mirna_df) == "B_mean"] <- "mean_RC"

deg_mirna_ids  <- deg_mirna_df$mirna
message(sprintf("  DEG miRNAs (SC vs RC): %d", length(deg_mirna_ids)))

# ── 1.3  Compute miRNA fold changes (SC → RC) for all DEG miRNAs ─────────────

mirna_fc_b <- data.frame(
  mirna  = rownames(mirna_comp1$full),
  log2FC = mirna_comp1$full$M,   # sign correction applied here
  stringsAsFactors = FALSE
) |>
  dplyr::filter(mirna %in% deg_mirna_ids)

message(sprintf("  miRNA FC computed for: %d miRNAs", nrow(mirna_fc_b)))

# ── 1.4  Compute mRNA fold changes from full NOISeq results ──────────────────
# Using comp1.gsea$full (column M = log2FC, rownames = gene symbols)
# Fallback: compute directly from mrna_norm

if (exists("comp1") && !is.null(comp1$full)) {
  mrna_fc_b <- data.frame(
    gene   = rownames(comp1$full),
    log2FC = comp1$full$M,
    stringsAsFactors = FALSE
  )
  message("  mRNA FC source: comp1$full")
  
} else {
  # Direct computation fallback
  mrna_fc_b <- data.frame(
    gene   = rownames(mrna_norm),
    log2FC = log2((mrna_norm[, "A375_RC"] + 1) /
                    (mrna_norm[, "A375_SC"] + 1)),
    stringsAsFactors = FALSE
  )
  message("  mRNA FC source: computed from mrna_norm (fallback)")
}

# =============================================================================
# SECTION 2: BATCH REVERSE multiMiR QUERY
# =============================================================================
# Query once with ALL unique pathway genes (not per-pathway) to avoid
# redundant API calls. Then split results by pathway afterward.

message("\n--- Section 2: Reverse multiMiR query ---")

all_pathway_genes <- unique(pathway_gene_df$gene)
message(sprintf("  Unique genes across all pathways: %d", length(all_pathway_genes)))
message("  Querying multiMiR (this takes 2-5 minutes)...")

# ── Validated interactions: pathway genes as targets ─────────────────────────

val_result <- tryCatch(
  multiMiR::get_multimir(
    target  = all_pathway_genes,
    table   = "validated",
    summary = FALSE
  ),
  error = function(e) {
    message(sprintf("  Validated query failed: %s", e$message))
    NULL
  }
)

# ── Predicted interactions: pathway genes as targets ─────────────────────────

pred_result <- tryCatch(
  multiMiR::get_multimir(
    target                = all_pathway_genes,
    table                 = "predicted",
    predicted.cutoff      = 30,
    predicted.cutoff.type = "p",
    predicted.site        = "all",
    summary               = FALSE
  ),
  error = function(e) {
    message(sprintf("  Predicted query failed: %s", e$message))
    NULL
  }
)

# ── Combine results ───────────────────────────────────────────────────────────

raw_b <- dplyr::bind_rows(
  if (!is.null(val_result))  val_result@data  |> dplyr::mutate(type = "validated"),
  if (!is.null(pred_result)) pred_result@data |> dplyr::mutate(type = "predicted")
)

message(sprintf("  Raw interactions returned: %d", nrow(raw_b)))

# Standardize to: mirna, gene, type
# multiMiR columns: mature_mirna_id, target_symbol
target_raw <- raw_b |>
  dplyr::select(
    mirna = mature_mirna_id,
    gene  = target_symbol,
    type
  ) |>
  dplyr::filter(
    !is.na(mirna), mirna != "",
    !is.na(gene),  gene  != ""
  ) |>
  # Validation priority: if pair has both, keep as validated
  dplyr::group_by(mirna, gene) |>
  dplyr::summarise(
    type = ifelse("validated" %in% type, "validated", "predicted"),
    .groups = "drop"
  )

message(sprintf("  Unique miRNA-gene pairs after dedup: %d", nrow(target_raw)))

# ── Filter 1: keep only our 47 DEG miRNAs ────────────────────────────────────

target_deg <- target_raw |>
  dplyr::filter(mirna %in% deg_mirna_ids)

message(sprintf("  After DEG miRNA filter: %d pairs", nrow(target_deg)))
message(sprintf("  DEG miRNAs with pathway targets: %d / %d",
                length(unique(target_deg$mirna)), length(deg_mirna_ids)))

# ── Filter 2: anti-correlation ────────────────────────────────────────────────

target_anticor <- target_deg |>
  dplyr::inner_join(mirna_fc_b |> dplyr::select(mirna, log2FC_mirna = log2FC),
                    by = "mirna") |>
  dplyr::inner_join(mrna_fc_b  |> dplyr::select(gene,  log2FC_mrna  = log2FC),
                    by = "gene") |>
  dplyr::filter(sign(log2FC_mirna) != sign(log2FC_mrna)) |>
  dplyr::distinct(mirna, gene, .keep_all = TRUE)

message(sprintf("  After anti-correlation filter: %d pairs", nrow(target_anticor)))
message(sprintf("  Validated pairs: %d", sum(target_anticor$type == "validated")))

# =============================================================================
# SECTION 3: JOIN PATHWAY ANNOTATIONS
# =============================================================================

message("\n--- Section 3: Pathway annotation ---")

# Attach pathway membership (already redundancy-removed from Method A)
pw_lookup <- pathway_gene_df |>
  dplyr::select(gene, pathway, NES, source)

method_b_network <- target_anticor |>
  dplyr::inner_join(pw_lookup, by = "gene") |>
  dplyr::rename(gene_pathway = pathway, pathway_NES = NES)

message(sprintf("  Method B network edges: %d", nrow(method_b_network)))
message(sprintf("  Pathways covered: %d", length(unique(method_b_network$gene_pathway))))
message("  Coverage per pathway:")
print(sort(table(method_b_network$gene_pathway), decreasing = TRUE))



# Save full Method B network
write.csv(method_b_network,
          file.path(OUTPUT_DIR, "method_b", "method_b_full_network.csv"),
          row.names = FALSE)

# =============================================================================
# SECTION 4: PER-PATHWAY SUBNETWORK VISUALIZATION
# =============================================================================

message("\n--- Section 4: Per-pathway plots ---")

pathway_colors <- c(
  "DNA Damage Response"      = "#E64B35",
  "Cell Cycle / G2M"         = "#F39B7F",
  "Homologous Recombination"  = "#D62728",
  "OXPHOS / Complex I"       = "#9467BD",
  "Iron Transport"            = "#8C564B",
  "Interferon Signaling"      = "#4DBBD5",
  "Lysosome / Autophagy"      = "#00A087",
  "Cytokine Signaling"        = "#AEC7E8",
  "Ferroptosis"               = "#C5B0D5",
  "EMT / Mesenchymal"         = "#2CA02C",
  "p53 Signaling"             = "#FF7F0E",
  "PI3K / AKT / mTOR"        = "#1F77B4",
  "Hypoxia / HIF1A"           = "#BCBD22",
  "WNT Signaling"             = "#17BECF",
  "Apoptosis"                 = "#7F7F7F",
  "miRNA"                     = "#7B2D8B"
)

pathways_b <- names(sort(table(method_b_network$gene_pathway), decreasing = TRUE))
pathways_b <- pathways_b[table(method_b_network$gene_pathway)[pathways_b] >= 3]
message(sprintf("  Pathways with ≥3 edges: %d", length(pathways_b)))

for (pw_label in pathways_b) {
  
  pw_edges <- method_b_network |>
    dplyr::filter(gene_pathway == pw_label)
  
  # Node table
  pw_mirna_nodes <- data.frame(
    name      = unique(pw_edges$mirna),
    node_type = "miRNA",
    category  = "miRNA",
    log2FC    = mirna_fc_b$log2FC[match(unique(pw_edges$mirna), mirna_fc_b$mirna)],
    stringsAsFactors = FALSE
  )
  pw_gene_nodes <- pw_edges |>
    dplyr::distinct(gene, log2FC_mrna) |>
    dplyr::rename(name = gene, log2FC = log2FC_mrna) |>
    dplyr::mutate(node_type = "mRNA", category = pw_label)
  
  pw_nodes <- dplyr::bind_rows(pw_mirna_nodes, pw_gene_nodes)
  
  g_pw <- igraph::graph_from_data_frame(
    d        = pw_edges[, c("mirna","gene","type","log2FC_mirna","log2FC_mrna")],
    directed = TRUE,
    vertices = pw_nodes
  )
  
  n_nodes <- igraph::vcount(g_pw)
  n_edges <- igraph::ecount(g_pw)
  
  # Label strategy: all nodes if small, top-degree + all miRNAs if large
  if (n_nodes <= 35) {
    label_set <- pw_nodes$name
  } else {
    top_genes <- names(sort(igraph::degree(g_pw), decreasing = TRUE))[1:20]
    label_set <- unique(c(unique(pw_edges$mirna), top_genes))
  }
  
  tg_pw <- tidygraph::as_tbl_graph(g_pw) |>
    tidygraph::activate(nodes) |>
    dplyr::mutate(
      label  = ifelse(name %in% label_set, name, NA_character_),
      degree = igraph::degree(g_pw)
    )
  
  pw_color <- pathway_colors[pw_label]
  if (is.na(pw_color)) pw_color <- "#8491B4"
  layout_type <- if (n_nodes <= 25) "stress" else "fr"
  
  p_pw <- ggraph::ggraph(tg_pw, layout = layout_type) +
    ggraph::geom_edge_arc(
      aes(color    = log2FC_mirna > 0,
          linetype = type),
      arrow     = arrow(length = unit(1.5, "mm"), type = "closed"),
      end_cap   = ggraph::circle(2.5, "mm"),
      alpha     = 0.55,
      linewidth = 0.45,
      strength  = 0.15
    ) +
    ggraph::geom_node_point(
      aes(size  = degree + 1,
          shape = node_type,
          color = ifelse(node_type == "miRNA", "miRNA", pw_label)),
      alpha = 0.92
    ) +
    ggraph::geom_node_text(
      aes(label = label),
      size         = 2.7,
      repel        = TRUE,
      fontface     = "italic",
      max.overlaps = 30,
      segment.size = 0.2,
      segment.color= "grey50"
    ) +
    scale_color_manual(
      values = setNames(c("#7B2D8B", pw_color), c("miRNA", pw_label)),
      name   = "Node type"
    ) +
    scale_shape_manual(values = c(miRNA = 18, mRNA = 16), guide = "none") +
    scale_size_continuous(range = c(2.5, 9), name = "Degree") +
    ggraph::scale_edge_color_manual(
      values = c(`TRUE`  = "#B2182B", `FALSE` = "#2166AC"),
      labels = c("miRNA ↑ in RC", "miRNA ↓ in RC"),
      name   = "miRNA FC"
    ) +
    ggraph::scale_edge_linetype_manual(
      values = c(validated = "solid", predicted = "dashed"),
      name   = "Evidence"
    ) +
    labs(
      title    = sprintf("[Method B] %s — Reverse Query Subnetwork", pw_label),
      subtitle = sprintf(
        "%d miRNAs → %d genes | %d edges (%d validated, %d predicted) | anti-correlated SC→RC",
        sum(igraph::V(g_pw)$node_type == "miRNA"),
        sum(igraph::V(g_pw)$node_type == "mRNA"),
        n_edges,
        sum(pw_edges$type == "validated"),
        sum(pw_edges$type == "predicted")
      )
    ) +
    ggraph::theme_graph(base_family = "sans") +
    theme(
      plot.title      = element_text(face = "bold", size = 12),
      plot.subtitle   = element_text(size = 8),
      legend.position = "bottom",
      legend.text     = element_text(size = 7)
    )
  
  fname <- gsub("[/ ]", "_", pw_label)
  ggplot2::ggsave(
    file.path(OUTPUT_DIR, "method_b", sprintf("methodB_subnetwork_%s.pdf", fname)),
    p_pw, width = 10, height = 8, device = cairo_pdf
  )
  
  igraph::write_graph(
    g_pw,
    file   = file.path(OUTPUT_DIR, "method_b",
                       sprintf("methodB_subnetwork_%s.graphml", fname)),
    format = "graphml"
  )
  
  message(sprintf("  [%s] %d nodes, %d edges — saved.", pw_label, n_nodes, n_edges))
}

# =============================================================================
# SECTION 5: COMBINED OVERVIEW — METHOD B
# =============================================================================

message("\n--- Section 5: Combined Method B network ---")

all_mirna_nodes_b <- data.frame(
  name      = unique(method_b_network$mirna),
  node_type = "miRNA",
  category  = "miRNA",
  log2FC    = mirna_fc_b$log2FC[match(unique(method_b_network$mirna), mirna_fc_b$mirna)],
  stringsAsFactors = FALSE
)
all_gene_nodes_b <- method_b_network |>
  dplyr::distinct(gene, gene_pathway, log2FC_mrna) |>
  dplyr::rename(name = gene, category = gene_pathway, log2FC = log2FC_mrna) |>
  dplyr::mutate(node_type = "mRNA")

all_nodes_b <- dplyr::bind_rows(all_mirna_nodes_b, all_gene_nodes_b)

g_combined_b <- igraph::graph_from_data_frame(
  d = method_b_network[, c("mirna","gene","type","gene_pathway",
                           "log2FC_mirna","log2FC_mrna")],
  directed = TRUE,
  vertices = all_nodes_b
)

message(sprintf("  Combined Method B graph: %d nodes, %d edges",
                igraph::vcount(g_combined_b), igraph::ecount(g_combined_b)))

igraph::write_graph(
  g_combined_b,
  file   = file.path(OUTPUT_DIR, "method_b", "methodB_combined.graphml"),
  format = "graphml"
)

# ── miRNA × pathway regulatory map ───────────────────────────────────────────

mirna_pw_map_b <- method_b_network |>
  dplyr::group_by(mirna, gene_pathway) |>
  dplyr::summarise(
    n_targets       = n(),
    n_validated     = sum(type == "validated"),
    target_genes    = paste(sort(gene), collapse = "; "),
    mirna_log2FC    = unique(log2FC_mirna),
    mirna_direction = ifelse(unique(log2FC_mirna) > 0, "Up_in_RC", "Down_in_RC"),
    .groups = "drop"
  ) |>
  dplyr::arrange(gene_pathway, desc(n_targets))

write.csv(mirna_pw_map_b,
          file.path(OUTPUT_DIR, "method_b", "methodB_mirna_pathway_map.csv"),
          row.names = FALSE)

# ── Heatmap +──────────────────────────────────────────────────────────────────

map_wide_b <- mirna_pw_map_b |>
  dplyr::select(mirna, gene_pathway, n_targets) |>
  tidyr::pivot_wider(names_from = gene_pathway, values_from = n_targets,
                     values_fill = 0L) |>
  as.data.frame()
rownames(map_wide_b) <- map_wide_b$mirna
map_wide_b$mirna     <- NULL
map_mat_b <- as.matrix(map_wide_b)
map_mat_b <- map_mat_b[rowSums(map_mat_b) > 0, colSums(map_mat_b) > 0, drop = FALSE]

if (nrow(map_mat_b) >= 2 && ncol(map_mat_b) >= 2) {
  
  dir_ann_b <- mirna_pw_map_b |>
    dplyr::distinct(mirna, mirna_direction) |>
    dplyr::filter(mirna %in% rownames(map_mat_b)) |>
    tibble::column_to_rownames("mirna")
  
  pheatmap::pheatmap(
    map_mat_b,
    color             = colorRampPalette(c("white","#FDD0A2","#E64B35"))(50),
    annotation_row    = dir_ann_b,
    annotation_colors = list(
      mirna_direction = c(Up_in_RC = "#B2182B", Down_in_RC = "#2166AC")
    ),
    cluster_rows      = TRUE,
    cluster_cols      = TRUE,
    fontsize_row      = 7,
    fontsize_col      = 8,
    display_numbers   = TRUE,
    number_format     = "%d",
    number_color      = "black",
    main              = "Method B: miRNA → Pathway Regulatory Map\n(n_targets | reverse query | 47 DEG miRNAs | SC vs RC)",
    filename          = file.path(OUTPUT_DIR, "method_b",
                                  "methodB_heatmap_mirna_pathway.pdf"),
    width  = 10,
    height = max(5, nrow(map_mat_b) * 0.28 + 3)
  )
  message("  Method B heatmap saved.")
}

# =============================================================================
# SECTION 6: CONSOLIDATED PRIORITY TABLE
# =============================================================================
# Ranks miRNA-gene pairs by: validated > predicted, then n_pathways regulated,
# then |mirna_log2FC|. This is the shortlist for wet lab follow-up.


priority_table_b <- method_b_network |>
  dplyr::group_by(mirna, gene) |>
  dplyr::summarise(
    type            = unique(type),
    pathways        = paste(sort(unique(gene_pathway)), collapse = "; "),
    n_pathways      = n_distinct(gene_pathway),
    log2FC_mirna    = unique(log2FC_mirna),
    log2FC_mrna     = unique(log2FC_mrna),
    mirna_direction = ifelse(unique(log2FC_mirna) > 0, "Up_in_RC", "Down_in_RC"),
    .groups         = "drop"
  ) |>
  dplyr::mutate(
    evidence_score = ifelse(type == "validated", 2, 1),
    iron_gene      = gene %in% c("FTH1","NCOA4","TFRC","SLC7A11","GPX4",
                                 "IREB2","FTL","HMOX1","STEAP3","CYBRD1"),
    priority_score = evidence_score + n_pathways + (iron_gene * 3)
  ) |>
  dplyr::arrange(desc(priority_score), desc(abs(log2FC_mirna)))

write.csv(priority_table_b,
          file.path(OUTPUT_DIR, "method_b", "methodB_priority_wetlab.csv"),
          row.names = FALSE)

message("\nTop 20 miRNA-gene pairs for wet lab follow-up:")
print(head(priority_table_b[, c("mirna","gene","type","pathways",
                                "log2FC_mirna","log2FC_mrna",
                                "iron_gene","priority_score")], 20))

# =============================================================================
# PRIORITY 2: CO-EXPRESSION / WGCNA MODULE ANALYSIS
# =============================================================================
#
# Design rationale for n = 4:
#   Classical WGCNA requires n ≥ 15 for stable modules. With 4 conditions we
#   adopt a "micro-WGCNA" strategy:
#     (a) Select top 5,000 genes by variance (rank-based, not threshold-based)
#     (b) Signed network; soft power defaults to 6 when undetermined
#     (c) Dynamic tree cut + module merging to reduce fragmentation
#     (d) All outputs labelled EXPLORATORY; interpretation focuses on gene-set
#         overlap (Fisher test) rather than module stability per se.
# =============================================================================
library(dynamicTreeCut)
message("\n=== PRIORITY 2: Co-expression Module Analysis (n=4 Exploratory) ===")

# ── Gene sets ─────────────────────────────────────────────────────────────────
# Extend iron set and use proper HGNC symbols throughout

IRON_EXTENDED <- c(
  IRON_GENES,                                             # from Section 2
  "FTL", "CYBRD1", "TF", "LCN2", "HAMP", "SLC40A1",
  "STEAP3", "HMOX1", "ABCB1", "SLC11A2", "FBXL5",
  "ACO1", "ACO2", "ISCU", "CISD1", "PCBP1"
)

STEM_GENES <- c(
  "POU5F1",  # OCT4 — proper HGNC
  "SOX2", "NANOG", "KLF4", "MYC",
  "CD44",
  "PROM1",   # CD133 — proper HGNC
  "ALDH1A1", "ABCB5", "KDM5B",  # KDM5B = JARID1B
  "ZEB1", "TWIST1", "SNAI1", "SNAI2", "VIM",
  "NOTCH1", "WNT5A", "AXL", "MITF"
)
# Keep "OCT4" and "CD133" as aliases in case they appear in the matrix under
# those names (some annotation pipelines use them):
STEM_GENES <- unique(c(STEM_GENES, "OCT4", "CD133"))

# ── P2.1  Prepare expression matrix ──────────────────────────────────────────

mrna_tmm_log <- log2(mrna_norm + 1)

# Filter: require expression in all 4 samples and non-zero variance
gene_min_expr <- apply(mrna_tmm_log, 1, min)
gene_var_all  <- apply(mrna_tmm_log, 1, var)
keep_base     <- gene_min_expr > 0 & gene_var_all > 0

# Rank-based top-N selection (avoids arbitrary absolute variance thresholds)
# Using 5,000 genes; reduces computational load while retaining strong signal
TOP_NGENES  <- min(5000L, sum(keep_base))
var_rank    <- rank(-gene_var_all[keep_base], ties.method = "first")
wgcna_mat   <- mrna_tmm_log[keep_base, ][var_rank <= TOP_NGENES, ]
message(sprintf("  Genes entering WGCNA: %d (top %d by variance)", nrow(wgcna_mat), TOP_NGENES))

# WGCNA convention: samples in rows
wgcna_t <- t(wgcna_mat)

# ── P2.2  Sample and gene quality check ──────────────────────────────────────

gsg <- WGCNA::goodSamplesGenes(wgcna_t, verbose = 0)
if (!gsg$allOK) {
  n_bad_genes   <- sum(!gsg$goodGenes)
  n_bad_samples <- sum(!gsg$goodSamples)
  message(sprintf("  goodSamplesGenes: removing %d genes and %d samples",
                  n_bad_genes, n_bad_samples))
  wgcna_t   <- wgcna_t[gsg$goodSamples, gsg$goodGenes]
  wgcna_mat <- t(wgcna_t)
}
message(sprintf("  WGCNA matrix after QC: %d genes x %d samples",
                ncol(wgcna_t), nrow(wgcna_t)))

# ── P2.3  Soft-threshold selection ───────────────────────────────────────────

powers  <- c(1:20)
sft_out <- WGCNA::pickSoftThreshold(
  wgcna_t,
  powerVector  = powers,
  networkType  = "signed",
  RsquaredCut  = 0.80,
  verbose      = 0
)

soft_power <- sft_out$powerEstimate
if (is.na(soft_power)) {
  soft_power <- 6L
  message(sprintf(
    "  Soft power undetermined (expected with n=4); using default = %d", soft_power))
} else {
  message(sprintf("  Soft power selected: %d", soft_power))
}
enableWGCNAThreads()
# Scale-free topology diagnostic plot
pdf(file.path(OUTPUT_DIR, "wgcna_soft_threshold.pdf"), width = 9, height = 4)
par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
R2_vec <- -sign(sft_out$fitIndices$slope) * sft_out$fitIndices$SFT.R.sq
plot(powers, R2_vec,
     xlab = "Soft Threshold (power)", ylab = "R² (scale-free topology)",
     main = "Scale-Free Topology Fit", type = "n", ylim = c(0, 1))
text(powers, R2_vec, labels = powers, col = "#E64B35", cex = 0.8)
abline(h = 0.80, col = "#2166AC", lty = 2, lwd = 1.2)
legend("bottomright", legend = "Target R² = 0.80", lty = 2,
       col = "#2166AC", bty = "n", cex = 0.8)
plot(powers, sft_out$fitIndices$mean.k.,
     xlab = "Soft Threshold (power)", ylab = "Mean Connectivity",
     main = "Mean Network Connectivity", type = "n")
text(powers, sft_out$fitIndices$mean.k., labels = powers, col = "#E64B35", cex = 0.8)
dev.off()
message("  Soft threshold diagnostic plot saved.")

# ── P2.4  Adjacency and topological overlap matrix ───────────────────────────

adjacency <- WGCNA::adjacency(wgcna_t, power = soft_power, type = "signed")

TOM     <- WGCNA::TOMsimilarity(adjacency)
dissTOM <- 1 - TOM
rownames(dissTOM) <- colnames(dissTOM) <- colnames(wgcna_t)

# ── P2.5  Hierarchical clustering and dynamic module detection ────────────────

gene_tree <- flashClust::flashClust(as.dist(dissTOM), method = "average")

# deepSplit = 2: moderate sensitivity; minClusterSize = 10 suits n=4
raw_modules   <- cutreeDynamic(
  dendro             = gene_tree,
  distM              = dissTOM,
  deepSplit          = 2,
  pamRespectsDendro  = FALSE,
  minClusterSize     = 10
)
raw_colors <- WGCNA::labels2colors(raw_modules)
message(sprintf("  Raw modules before merging: %d (+grey unassigned)",
                length(unique(raw_colors[raw_colors != "grey"]))))

# ── P2.6  Merge similar modules ───────────────────────────────────────────────
# With n=4, dynamic cut is prone to over-splitting. Merge modules whose
# eigengene Pearson correlation > 0.85 (dissimilarity cutHeight = 0.15).

merge_res     <- WGCNA::mergeCloseModules(
  wgcna_t,
  colors     = raw_colors,
  cutHeight  = 0.15,
  verbose    = 0
)
module_colors <- merge_res$colors
names(module_colors) <- colnames(wgcna_t)  # ADD THIS


MEs           <- WGCNA::orderMEs(merge_res$newMEs)


n_mods <- length(unique(module_colors[module_colors != "grey"]))
message(sprintf("  Modules after merging (cutHeight = 0.15): %d", n_mods))
print(sort(table(module_colors), decreasing = TRUE))

# ── P2.7  Dendrogram + module color bar ──────────────────────────────────────
# This is the canonical WGCNA visualization; absent from the original file.

pdf(file.path(OUTPUT_DIR, "wgcna_dendrogram_modules.pdf"), width = 11, height = 5)
WGCNA::plotDendroAndColors(
  dendro       = gene_tree,
  colors       = cbind(raw_colors, module_colors),
  groupLabels  = c("Before merge", "After merge"),
  dendroLabels = FALSE,
  hang         = 0.03,
  addGuide     = TRUE,
  guideHang    = 0.05,
  marAll       = c(0, 5, 2, 0),
  main         = sprintf(
    "Gene Co-expression Dendrogram | Signed WGCNA | n=4 exploratory\n%d modules | soft power = %d | top %d genes by variance",
    n_mods, soft_power, nrow(wgcna_mat))
)
dev.off()
message("  Dendrogram + module color bar saved.")

# ── P2.8  Module–trait correlation ───────────────────────────────────────────
# Four binary traits capture the experimental design contrasts.

trait_mat <- matrix(
  c(1, 1, 0, 0,   # Resistance:    1 = resistant line
    0, 1, 0, 1,   # DrugTreatment: 1 = drug-treated
    1, 0, 0, 0,   # RC_specific:   1 = resistant control only
    0, 0, 1, 0),  # SC_specific:   1 = sensitive control only
  nrow = 4, ncol = 4,
  dimnames = list(
    c("A375_RC", "A375_R10", "A375_SC", "A375_S10"),
    c("Resistance", "DrugTreatment", "RC_specific", "SC_specific")
  )
)

# Align MEs row order to trait_mat row order (important!)
shared_samples   <- intersect(rownames(MEs), rownames(trait_mat))
MEs_aligned      <- MEs[shared_samples, , drop = FALSE]
trait_aligned    <- trait_mat[shared_samples, , drop = FALSE]

module_trait_cor  <- cor(MEs_aligned, trait_aligned, use = "pairwise.complete.obs")
module_trait_pval <- WGCNA::corPvalueStudent(module_trait_cor, nSamples = nrow(MEs_aligned))

# Build text matrix for the labeled heatmap
text_mat <- matrix(
  paste0(round(module_trait_cor, 2), "\n(", signif(module_trait_pval, 1), ")"),
  nrow = nrow(module_trait_cor),
  ncol = ncol(module_trait_cor),
  dimnames = dimnames(module_trait_cor)
)

pdf(file.path(OUTPUT_DIR, "wgcna_module_trait_correlation.pdf"),
    width = 8,
    height = max(4, nrow(module_trait_cor) * 0.45 + 2))
WGCNA::labeledHeatmap(
  Matrix        = module_trait_cor,
  xLabels       = colnames(trait_aligned),
  yLabels       = rownames(module_trait_cor),
  ySymbols      = rownames(module_trait_cor),
  colorLabels   = FALSE,
  colors        = WGCNA::blueWhiteRed(50),
  textMatrix    = text_mat,
  setStdMargins = FALSE,
  cex.text      = 0.72,
  zlim          = c(-1, 1),
  main          = "Module-Trait Correlation\n(r, p-value from Student t | n=4 exploratory)"
)
dev.off()
message("  Module-trait correlation heatmap saved.")

# ── P2.9  Identify resistance-associated modules ──────────────────────────────

# Primary threshold: |r| > 0.7 with 'Resistance'
resistance_modules <- rownames(module_trait_cor)[
  abs(module_trait_cor[, "Resistance"]) > 0.7
]

if (length(resistance_modules) == 0) {
  # n=4 often yields |r| in the 0.50–0.70 range without reaching 0.7.
  # Fallback: take top 3 modules most correlated with Resistance.
  top3_idx <- order(abs(module_trait_cor[, "Resistance"]), decreasing = TRUE)[
    1:min(3L, nrow(module_trait_cor))
  ]
  resistance_modules <- rownames(module_trait_cor)[top3_idx]
  message(sprintf(
    "  No module passed |r|>0.7; using top %d by |r| with Resistance",
    length(resistance_modules)))
}
message(sprintf("  Resistance-associated modules (%d): %s",
  length(resistance_modules),
  paste(resistance_modules, collapse = ", ")))

# ── P2.10  Module membership (kME) and hub gene identification ─────────────────
# kME = Pearson r between a gene's expression profile and the module eigengene.
# kME > 0.80 is the conventional hub threshold.
# Before kME_all computation, replace:
#kME_all <- WGCNA::signedKME(wgcna_t, MEs_aligned)
# ── P2.10 fix: kME via direct correlation ────────────────────────────────────
# cor(wgcna_t, MEs_aligned) = gene × eigengene Pearson r
# Column names stay as "MEblue", "MEmidnightblue" etc. — no renaming needed

kME_all <- cor(wgcna_t, MEs_aligned, use = "pairwise.complete.obs")

hub_genes  <- list()
hub_all_df <- data.frame(module       = character(),
                         gene         = character(),
                         kME          = numeric(),
                         iron_related = logical(),
                         stem_related = logical(),
                         mapk_related = logical(),
                         stringsAsFactors = FALSE)

for (mod_ME in resistance_modules) {
  
  mod_color <- gsub("^ME", "", mod_ME)
  mod_genes <- names(module_colors)[module_colors == mod_color]
  
  if (length(mod_genes) < 3) {
    message(sprintf("  Module '%s': %d genes — skipping.", mod_color, length(mod_genes)))
    next
  }
  
  # Use mod_ME directly as column name — matches MEs_aligned exactly
  if (!mod_ME %in% colnames(kME_all)) {
    message(sprintf("  Eigengene '%s' not found in kME_all — skipping.", mod_ME))
    next
  }
  
  kME_mod <- sort(kME_all[mod_genes, mod_ME], decreasing = TRUE)
  top_n   <- min(30L, length(kME_mod))
  hub_genes[[mod_color]] <- names(kME_mod)[seq_len(top_n)]
  
  iron_in <- names(kME_mod)[names(kME_mod) %in% IRON_EXTENDED]
  stem_in <- names(kME_mod)[names(kME_mod) %in% STEM_GENES]
  
  message(sprintf("  Module '%s' (%d genes) | top hubs: %s",
                  mod_color, length(mod_genes),
                  paste(names(kME_mod)[seq_len(min(5, length(kME_mod)))], collapse = ", ")))
  if (length(iron_in) > 0)
    message("    Iron genes: ", paste(iron_in[seq_len(min(6,length(iron_in)))], collapse=", "))
  if (length(stem_in) > 0)
    message("    Stem genes: ", paste(stem_in[seq_len(min(6,length(stem_in)))], collapse=", "))
  
  hub_all_df <- dplyr::bind_rows(hub_all_df, data.frame(
    module       = mod_color,
    gene         = names(kME_mod),
    kME          = as.numeric(kME_mod),
    iron_related = names(kME_mod) %in% IRON_EXTENDED,
    stem_related = names(kME_mod) %in% STEM_GENES,
    mapk_related = names(kME_mod) %in% MAPK_AKT_GENES,
    stringsAsFactors = FALSE
  ))
}

message(sprintf("\n  Hub genes populated for %d modules: %s",
                length(hub_genes), paste(names(hub_genes), collapse=", ")))
write.csv(hub_all_df, file.path(OUTPUT_DIR, "wgcna_hub_genes.csv"), row.names=FALSE)
# ── P2.11  Module gene-set overlap (Fisher's exact test) ─────────────────────
# Background: all genes in WGCNA (colnames of wgcna_t)

bg_size <- ncol(wgcna_t)

genesets_for_overlap <- list(
  Iron_Ferroptosis = IRON_EXTENDED,
  Stemness         = STEM_GENES,
  Autophagy        = AUTOPHAGY_GENES,
  MAPK_AKT         = MAPK_AKT_GENES
)

overlap_rows <- list()
for (mod_ME in resistance_modules) {
  mod_color <- gsub("^ME", "", mod_ME)
  mod_genes <- names(module_colors)[module_colors == mod_color]

  for (gs_name in names(genesets_for_overlap)) {
    gs    <- genesets_for_overlap[[gs_name]]
    hits  <- mod_genes[mod_genes %in% gs]
    k     <- length(hits)
    K     <- length(gs)
    N     <- bg_size
    n     <- length(mod_genes)
    pval  <- phyper(k - 1L, K, N - K, n, lower.tail = FALSE)
    OR    <- (k / n) / (K / N)

    overlap_rows[[length(overlap_rows) + 1]] <- data.frame(
      module        = mod_color,
      gene_set      = gs_name,
      module_size   = n,
      geneset_size  = K,
      overlap_n     = k,
      overlap_genes = paste(hits, collapse = "; "),
      p_value       = pval,
      odds_ratio    = OR,
      stringsAsFactors = FALSE
    )
  }
}

if (length(overlap_rows) > 0) {
  overlap_df <- dplyr::bind_rows(overlap_rows) |>
    dplyr::mutate(padj_BH = p.adjust(p_value, method = "BH")) |>
    dplyr::arrange(p_value)

  write.csv(overlap_df,
            file.path(OUTPUT_DIR, "wgcna_module_genesets_overlap.csv"),
            row.names = FALSE)
  message("  Gene-set overlap table saved:")
  print(overlap_df[, c("module", "gene_set", "overlap_n", "p_value", "padj_BH",
                        "odds_ratio")])
} else {
  overlap_df <- data.frame()
}

# ── P2.12  Hub gene expression heatmap per resistance module ─────────────────

for (mod_color in names(hub_genes)) {
  hubs         <- hub_genes[[mod_color]]
  hubs_present <- hubs[hubs %in% rownames(mrna_mat)]
  if (length(hubs_present) < 3) next

  # log2CPM z-scored for display (using mrna_mat, not the WGCNA subset)
  cpm_h <- sweep(mrna_mat[hubs_present, ], 2,
                 colSums(mrna_mat) / 1e6, FUN = "/")
  z_h   <- t(scale(t(log2(cpm_h + 1))))
  z_h[is.nan(z_h)] <- 0

  # Row annotation
  row_ann <- data.frame(
    Category = dplyr::case_when(
      hubs_present %in% IRON_EXTENDED  ~ "Iron/Ferroptosis",
      hubs_present %in% STEM_GENES     ~ "Stemness",
      hubs_present %in% MAPK_AKT_GENES ~ "MAPK/AKT",
      hubs_present %in% AUTOPHAGY_GENES~ "Autophagy",
      TRUE                              ~ "Other"
    ),
    row.names = hubs_present
  )

  col_ann <- data.frame(
    Condition = c("Resistant","Resistant","Sensitive","Sensitive"),
    Treatment = c("Control",  "Drug",     "Control",  "Drug"),
    row.names = colnames(mrna_mat)
  )

  hub_colors <- list(
    Category  = c("Iron/Ferroptosis" = "#E64B35", "Stemness"  = "#7B2D8B",
                  "MAPK/AKT"         = "#00A087", "Autophagy" = "#4DBBD5",
                  "Other"            = "#8491B4"),
    Condition = c(Resistant = "#CC0000", Sensitive = "#0066CC"),
    Treatment = c(Control = "#999999", Drug = "#FF9900")
  )

  pheatmap::pheatmap(
    z_h,
    color             = colorRampPalette(c("#2166AC","white","#B2182B"))(100),
    annotation_row    = row_ann,
    annotation_col    = col_ann,
    annotation_colors = hub_colors,
    cluster_cols      = FALSE,
    cluster_rows      = TRUE,
    show_rownames     = TRUE,
    fontsize_row      = 7,
    main              = sprintf(
      "Module '%s' — Top %d Hub Genes\nZ-scored log2CPM | Resistance-associated | EXPLORATORY",
      mod_color, length(hubs_present)),
    filename = file.path(OUTPUT_DIR,
                         sprintf("wgcna_heatmap_module_%s.pdf", mod_color)),
    width  = 7,
    height = max(5, length(hubs_present) * 0.22 + 2)
  )
  message(sprintf("  Hub gene heatmap saved: module '%s' (%d genes).",
                  mod_color, length(hubs_present)))
}

# ── P2.13  GSEA on kME-ranked module genes (fgsea) ───────────────────────────# 
#This is the step that gives modules biological meaning.


library(fgsea)
library(clusterProfiler)
library(org.Hs.eg.db)

# Build KEGG gene sets once (reused across modules)
kegg_sets <- tryCatch({
  # Download current KEGG pathways for human
  kegg_db <- clusterProfiler::download_KEGG("hsa")
  # kegg_db$KEGGPATHID2EXTID: pathway → entrez
  # Convert to named list of gene symbols
  id2sym <- tryCatch(
    setNames(mapIds(org.Hs.eg.db, keys = unique(unlist(kegg_db$KEGGPATHID2EXTID)),
                    column = "SYMBOL", keytype = "ENTREZID", multiVals = "first"),
             unique(unlist(kegg_db$KEGGPATHID2EXTID))),
    error = function(e) NULL
  )
  if (is.null(id2sym)) stop("symbol mapping failed")
  lapply(split(kegg_db$KEGGPATHID2EXTID$to,
               kegg_db$KEGGPATHID2EXTID$from), function(entrez_ids) {
                 syms <- id2sym[entrez_ids]
                 syms[!is.na(syms)]
               })
}, error = function(e) {
  message("  KEGG download failed: ", conditionMessage(e), " — using MSigDB fallback")
  NULL
})
# Which of the 24 overlap genes appear in TF fgsea leading edges?
target_genes <- c("TFRC","FTH1","TCIRG1","ATP6V1F","NDUFA9","HLA-A","IKBKG","ABL1",
                  "SNAI1","TWIST1","HFE","NDUFA2",
                  "CHEK1","MRE11","EYA1","HIPK2","APAF1",
                  "MDM4","CDC6","LIN52","FAS","ATP6V1C1","ATP5F1E","GSTM4","ETV2")

# From TF_fgsea_significant - find leading edge hits
fgsea_sig <- read.csv("C:\\bioinformatics\\output\\tf_enrich\\TF_fgsea_significant.csv")
lapply(seq_len(nrow(fgsea_sig)), function(i) {
  le <- unlist(strsplit(fgsea_sig$leadingEdge[i], ", "))
  hits <- intersect(le, target_genes)
  if (length(hits) > 0)
    data.frame(TF=fgsea_sig$pathway[i], NES=fgsea_sig$NES[i], 
               direction=fgsea_sig$direction_in_RC[i], hits=paste(hits, collapse=";"))
}) |> dplyr::bind_rows()
# Fallback: use your manually defined gene sets as a sanity check
manual_sets <- list(
  Iron_Ferroptosis = IRON_EXTENDED,
  Autophagy        = AUTOPHAGY_GENES,
  MAPK_AKT         = MAPK_AKT_GENES,
  Stemness         = STEM_GENES
)
# Add this helper function ONCE before the module loop
save_fgsea <- function(result_df, filepath) {
  df <- as.data.frame(result_df)
  # leadingEdge is a list column — collapse to semicolon-separated string
  if ("leadingEdge" %in% colnames(df)) {
    df$leadingEdge <- sapply(df$leadingEdge, paste, collapse = ";")
  }
  write.csv(df, filepath, row.names = FALSE)
}
for (mod_ME in resistance_modules) {
  mod_color <- gsub("^ME", "", mod_ME)
  mod_genes <- names(module_colors)[module_colors == mod_color]
  if (length(mod_genes) < 10) next
  if (!mod_ME %in% colnames(kME_all)) next
  
  # Ranked list: ALL module genes by kME (not just top 30)
  ranked_kme <- sort(kME_all[mod_genes, mod_ME], decreasing = TRUE)
  
  message(sprintf("\n  GSEA for module '%s' (%d genes)...", mod_color, length(ranked_kme)))
  
  # ── Manual gene sets (always works, no internet needed) ──────────────────
  fgsea_manual <- tryCatch(
    fgsea::fgsea(
      pathways   = manual_sets,
      stats      = ranked_kme,
      minSize    = 2,
      maxSize    = 500,
      nPermSimple= 10000
    ),
    error = function(e) { message("  fgsea manual failed: ", e$message); NULL }
  )
  if (!is.null(fgsea_manual) && nrow(fgsea_manual) > 0) {
    fgsea_manual <- fgsea_manual[order(fgsea_manual$pval), ]
    message(sprintf("  Manual gene sets:"))
    print(fgsea_manual[, c("pathway","NES","pval","padj","size")])
    save_fgsea(
      as.data.frame(fgsea_manual),
      file.path(OUTPUT_DIR, sprintf("wgcna_manual_fgsea_%s.csv", mod_color))
    )
    fgsea_kegg <- tryCatch(
      fgsea::fgsea(
        pathways   = kegg_sets,
        stats      = ranked_kme,
        minSize    = 5,
        maxSize    = 500,
        nPermSimple= 10000
      ),
      error = function(e) { message("  fgsea KEGG failed: ", e$message); NULL }
    )
    if (!is.null(fgsea_kegg) && nrow(fgsea_kegg) > 0) {
      fgsea_kegg_sig <- fgsea_kegg[fgsea_kegg$padj < 0.5, ]  # relaxed: n is tiny
      fgsea_kegg_sig <- fgsea_kegg_sig[order(fgsea_kegg_sig$NES, decreasing = TRUE), ]
      if (nrow(fgsea_kegg_sig) > 0) {
        message(sprintf("  KEGG significant (padj<0.25): %d pathways", nrow(fgsea_kegg_sig)))
        print(head(fgsea_kegg_sig[, c("pathway","NES","pval","padj","size")], 15))
      } else {
        message("  KEGG: no pathways at padj<0.25 (expected with small module)")
      }
      save_fgsea(
        as.data.frame(fgsea_kegg),
        file.path(OUTPUT_DIR, sprintf("wgcna_kegg_fgsea_%s.csv", mod_color))
      )
    }
  }
}
# ── P2.14  Export full module gene assignment table ───────────────────────────
# Critical for downstream use: which gene belongs to which module.

module_assign_df <- data.frame(
  gene                  = colnames(wgcna_t),
  module                = module_colors,
  in_resistance_module  = module_colors %in% gsub("^ME", "", resistance_modules),
  iron_related          = colnames(wgcna_t) %in% IRON_EXTENDED,
  stem_related          = colnames(wgcna_t) %in% STEM_GENES,
  mapk_related          = colnames(wgcna_t) %in% MAPK_AKT_GENES,
  stringsAsFactors = FALSE
) |>
  dplyr::arrange(module, gene)

write.csv(module_assign_df,
          file.path(OUTPUT_DIR, "wgcna_gene_module_assignments.csv"),
          row.names = FALSE)
message(sprintf("  Module assignments saved: %d genes across %d modules (+grey).",
                nrow(module_assign_df), n_mods))

# =============================================================================
# PRIORITY 3 (REWRITTEN): TRANSCRIPTION FACTOR ENRICHMENT ANALYSIS
# DESIGN RATIONALE:
#   Two conceptually distinct tools, each used for what it is suited for:
#
#   Layer 1 — TF activity:  decoupleR ULM on CollecTRI signed regulons
#             ULM fits expression ~ TF_activity × mor, so activation and
#             repression targets are handled algebraically.  This is the
#             correct method for signed networks.
#
#   Layer 2 — Pathway enrichment:  fgsea on Hallmark (MSigDB) gene sets.
#             Hallmark gene sets are curated, unsigned, and conceptually
#             homogeneous — exactly the use case fgsea was designed for.
#             Fallback: KEGG gene sets via msigdbr or clusterProfiler.
#
#   Layer 3 — Convergence:  TFs that are (a) differentially active by ULM
#             AND (b) known to regulate genes that drive the top enriched
#             pathways (regulon overlap test).  These are the strongest
#             mechanistic hypotheses for wet-lab follow-up.
#
# REQUIRED SESSION OBJECTS (from Sections 0-1 + Priority 2):
#   mrna_norm, mrna_mat, comp1, deg1
#   IRON_EXTENDED, STEM_GENES, AUTOPHAGY_GENES, MAPK_AKT_GENES
#   hub_genes (from WGCNA; optional — gracefully handled if absent)
#   OUTPUT_DIR
# =============================================================================

message("\n=== PRIORITY 3 (REWRITTEN): TF Enrichment Analysis ===")

# ── Package checks ─────────────────────────────────────────────────────────

needed_pkgs <- c("decoupleR", "fgsea", "msigdbr", "dplyr", "tidyr",
                 "ggplot2", "ggrepel", "pheatmap")
for (p in needed_pkgs) {
  if (!requireNamespace(p, quietly = TRUE))
    BiocManager::install(p, ask = FALSE)
}

suppressPackageStartupMessages({
  library(decoupleR)
  library(fgsea)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
})


# =============================================================================
# P3.0  Build ranked gene list (SC vs RC, full expressed genes)
# =============================================================================
# M column from NOISeq = log2(A/B) = log2(SC/RC)
# We negate so that positive values = higher in RC (resistance signal).
# This ranked list is used only for fgsea pathway enrichment (Layer 2).

full_sc_rc <- as.data.frame(comp1$result@results[[1]])
full_sc_rc$gene <- rownames(full_sc_rc)

ranked_genes <- full_sc_rc |>
  dplyr::filter(!is.na(M)) |>
  dplyr::arrange(dplyr::desc(M)) |>
  dplyr::distinct(gene, .keep_all = TRUE)

# Named vector: positive = higher in RC (resistance direction)
ranked_vec <- setNames(-ranked_genes$M, ranked_genes$gene)
message(sprintf("  Ranked gene list: %d genes | range: [%.2f, %.2f]",
                length(ranked_vec), min(ranked_vec), max(ranked_vec)))


# =============================================================================
# P3.1  LAYER 1: TF activity via decoupleR ULM (signed regulons)
# =============================================================================
# ULM regresses each gene's expression against the TF's signed regulon
# (mor = +1 activation, -1 repression).  The t-statistic of the TF coefficient
# is the activity score.  Activation and repression targets contribute in
# opposite directions automatically — no manual splitting needed.

message("\n--- Layer 1: TF activity (decoupleR ULM) ---")

# Expression matrix: log2(TMM + 1), genes × samples
tf_expr_mat <- log2(mrna_norm + 1)

# Load regulon: CollecTRI preferred, dorothea fallback
regulon_net  <- NULL
regulon_src  <- "none"

# Attempt 1: CollecTRI
regulon_net <- tryCatch({
  net <- decoupleR::get_collectri(organism = "human", split_complexes = FALSE)
  regulon_src <- "CollecTRI"
  message(sprintf("  Regulon: CollecTRI — %d TF-target pairs", nrow(net)))
  net
}, error = function(e) {
  message("  CollecTRI unavailable: ", conditionMessage(e))
  NULL
})

# Attempt 2: dorothea A+B+C
if (is.null(regulon_net)) {
  regulon_net <- tryCatch({
    net <- decoupleR::get_dorothea(organism = "human", levels = c("A","B","C"))
    regulon_src <- "DoRothEA_ABC"
    message(sprintf("  Regulon: DoRothEA A+B+C — %d TF-target pairs", nrow(net)))
    net
  }, error = function(e) {
    message("  DoRothEA also unavailable: ", conditionMessage(e))
    NULL
  })
}

if (is.null(regulon_net))
  stop("  No regulon database available. Check internet connection.")

# Run ULM
ulm_res <- tryCatch(
  decoupleR::run_ulm(
    mat     = as.matrix(tf_expr_mat),
    net     = regulon_net,
    .source = "source",
    .target = "target",
    .mor    = "mor",
    minsize = 5L
  ),
  error = function(e) {
    stop("decoupleR::run_ulm failed: ", conditionMessage(e))
  }
)

# Pivot to TF × sample matrix
ulm_wide <- ulm_res |>
  dplyr::filter(statistic == "ulm") |>
  dplyr::select(source, condition, score) |>
  tidyr::pivot_wider(names_from = "condition", values_from = "score") |>
  as.data.frame()

rownames(ulm_wide) <- ulm_wide$source
message(sprintf("  ULM: %d TF activity profiles computed.", nrow(ulm_wide)))

# Compute resistance contrast: delta = activity(RC) - activity(SC)
# Positive delta = TF more active in resistant cells
if (all(c("A375_RC", "A375_SC") %in% colnames(ulm_wide))) {
  ulm_wide$delta_RC_vs_SC <- ulm_wide$A375_RC - ulm_wide$A375_SC
} else {
  # Fallback if sample names differ
  sample_cols <- intersect(c("A375_RC","A375_R10","A375_SC","A375_S10"),
                           colnames(ulm_wide))
  message("  Available sample columns: ", paste(sample_cols, collapse = ", "))
  ulm_wide$delta_RC_vs_SC <- ulm_wide[[sample_cols[grep("RC|R_K|RK",
                                                         sample_cols)[1]]]] -
    ulm_wide[[sample_cols[grep("SC|S_K|SK",  sample_cols)[1]]]]
}

# Pathway/category flags
hub_genes_vec <- if (exists("hub_genes") && length(hub_genes) > 0)
  unlist(hub_genes, use.names = FALSE) else character(0)

ulm_wide <- ulm_wide |>
  dplyr::rename(TF = source) |>
  dplyr::mutate(
    is_DEG       = TF %in% deg1$gene,
    is_hub       = TF %in% hub_genes_vec,
    iron_related = TF %in% c(IRON_EXTENDED,
                             "NFE2L2","NRF2","BACH1","HIF1A","MTF1",
                             "TFEB", "STAT3", "SP1", "HMOX1"),
    mapk_akt     = TF %in% MAPK_AKT_GENES,
    stem_related = TF %in% STEM_GENES,
    # Activity direction label
    activity_direction = ifelse(delta_RC_vs_SC > 0,
                                "More active in RC", "Less active in RC")
  ) |>
  dplyr::arrange(dplyr::desc(abs(delta_RC_vs_SC)))

# Save full activity table
write.csv(ulm_wide,
          file.path(OUTPUT_DIR, "TF_activity_scores.csv"),
          row.names = FALSE)
message(sprintf("  TF activity scores saved (%d TFs).", nrow(ulm_wide)))


# =============================================================================
# P3.2  Prioritize TFs from ULM results
# =============================================================================
# Priority score: |delta| (effect size) + pathway relevance weights.
# Unlike the old approach, there is no p-value from ULM to use here —
# the t-score magnitude is the evidence.  We rank by |delta| and annotate.

# Define top-activity threshold: |delta| in top 20% (quantile 0.80)
delta_q80 <- quantile(abs(ulm_wide$delta_RC_vs_SC), 0.80, na.rm = TRUE) #islenmedi
tf_priority <- ulm_wide |> #ful idi deyisdim
  dplyr::mutate(
    priority_score = abs(delta_RC_vs_SC) +
      (is_DEG    * 1.5) +
      (is_hub    * 1.5) +
      (iron_related * 2.5) +
      (mapk_akt  * 1.5) +
      (stem_related * 1.0)
  ) |>
  dplyr::arrange(dplyr::desc(priority_score))

write.csv(tf_priority,
          file.path(OUTPUT_DIR, "TF_ULM_priority.csv"),
          row.names = FALSE)

message("\n  Top 15 TFs by priority score:")
print(head(tf_priority[, c("TF","delta_RC_vs_SC","activity_direction",
                            "priority_score","iron_related","is_DEG")], 15))


# =============================================================================
# P3.3  TF activity heatmap (top 40 by |delta|)
# =============================================================================

sample_cols_present <- intersect(
  c("A375_RC","A375_R10","A375_SC","A375_S10"),
  colnames(ulm_wide)
)

top40_tfs <- head(
  ulm_wide[order(abs(ulm_wide$delta_RC_vs_SC), decreasing = TRUE), ],
  40
)

if (nrow(top40_tfs) >= 5 && length(sample_cols_present) >= 2) {

  act_mat <- as.matrix(top40_tfs[, sample_cols_present])
  rownames(act_mat) <- top40_tfs$TF

  act_row_ann <- data.frame(
    Pathway   = dplyr::case_when(
      top40_tfs$iron_related ~ "Iron/Ferroptosis",
      top40_tfs$mapk_akt     ~ "MAPK/AKT",
      top40_tfs$stem_related ~ "Stemness",
      TRUE                   ~ "Other"
    ),
    Direction = top40_tfs$activity_direction,
    row.names = top40_tfs$TF
  )

  cond_map <- c(A375_RC="Resistant", A375_R10="Resistant",
                A375_SC="Sensitive", A375_S10="Sensitive")
  trt_map  <- c(A375_RC="Control",   A375_R10="Drug",
                A375_SC="Control",   A375_S10="Drug")
  act_col_ann <- data.frame(
    Condition = cond_map[sample_cols_present],
    Treatment = trt_map[sample_cols_present],
    row.names = sample_cols_present
  )

  ann_colors <- list(
    Pathway   = c("Iron/Ferroptosis" = "#E64B35", "MAPK/AKT" = "#00A087",
                  "Stemness"         = "#7B2D8B", "Other"    = "#8491B4"),
    Direction = c("More active in RC" = "#B2182B", "Less active in RC" = "#2166AC"),
    Condition = c(Resistant = "#CC0000", Sensitive = "#0066CC"),
    Treatment = c(Control = "#999999",   Drug      = "#FF9900")
  )

  pheatmap::pheatmap(
    act_mat,
    color             = colorRampPalette(c("#2166AC","white","#B2182B"))(100),
    annotation_row    = act_row_ann,
    annotation_col    = act_col_ann,
    annotation_colors = ann_colors,
    cluster_cols      = FALSE,
    cluster_rows      = TRUE,
    scale             = "row",
    fontsize_row      = 7,
    main              = sprintf(
      "TF Activity Scores (ULM | %s | z-scaled per TF)\nTop 40 by |RC - SC|",
      regulon_src),
    filename = file.path(OUTPUT_DIR, "TF_activity_heatmap.pdf"),
    width = 7, height = 11
  )
  message("  TF activity heatmap saved.")
}

# Bubble plot: top 30 by priority
top30 <- head(tf_priority, 15)
top30$pathway_cat <- dplyr::case_when(
  top30$iron_related ~ "Iron/Ferroptosis",
  top30$mapk_akt     ~ "MAPK/AKT",
  top30$stem_related ~ "Stemness",
  TRUE               ~ "Other"
)
top30$abs_delta <- abs(top30$delta_RC_vs_SC)

p_ulm <- ggplot2::ggplot(top30,
  ggplot2::aes(
    x     = delta_RC_vs_SC,
    y     = reorder(TF, delta_RC_vs_SC),
    size  = abs_delta,
    color = pathway_cat,
    shape = activity_direction
  )) +
  ggplot2::geom_point(alpha = 0.85) +
  ggplot2::scale_color_manual(
    values = c("Iron/Ferroptosis" = "#E64B35", "MAPK/AKT" = "#00A087",
               "Stemness" = "#7B2D8B", "Other" = "#8491B4"),
    name = "Pathway"
  ) +
  ggplot2::scale_shape_manual(
    values = c("More active in RC" = 17L, "Less active in RC" = 16L),
    name = "Activity in RC"
  ) +
  ggplot2::scale_size_continuous(range = c(3, 10), name = "|delta activity|") +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
                      color = "grey40", linewidth = 0.4) +
  ggplot2::labs(
    title    = "TF Activity: Resistant vs Sensitive Cells",
    subtitle = sprintf(
      "decoupleR ULM | %s | delta = activity(RC) - activity(SC)\n▲ = more active in RC (resistance driver candidate)",
      regulon_src),
    x = "Activity delta (RC - SC)",
    y = NULL
  ) +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"),
                 legend.position = "right")

ggplot2::ggsave(file.path(OUTPUT_DIR, "TF_ULM_bubble.pdf"),
                p_ulm, width = 11, height = 9, device = cairo_pdf)
message("  TF ULM bubble plot saved.")


# =============================================================================
# P3.4  LAYER 2: Pathway enrichment via fgsea on Hallmark gene sets
# =============================================================================
# Hallmark gene sets are unsigned and internally coherent — the correct use
# case for fgsea.  We use msigdbr to fetch them without downloading files.
# Fallback: KEGG gene sets from msigdbr (C2:KEGG subcollection).

message("\n--- Layer 2: Pathway enrichment (fgsea on Hallmark) ---")

pathway_sets <- NULL
pathway_src  <- "none"

if (requireNamespace("msigdbr", quietly = TRUE)) {

  # Attempt Hallmark
  pathway_sets <- tryCatch({
    h <- msigdbr::msigdbr(species = "Homo sapiens", category = "H")
    sets <- split(h$gene_symbol, h$gs_name)
    pathway_src <- "MSigDB_Hallmark"
    message(sprintf("  Hallmark gene sets: %d pathways", length(sets)))
    sets
  }, error = function(e) {
    message("  Hallmark unavailable: ", conditionMessage(e))
    NULL
  })

  # Fallback: KEGG (C2 subcollection)
  if (is.null(pathway_sets)) {
    pathway_sets <- tryCatch({
      k <- msigdbr::msigdbr(species = "Homo sapiens",
                            category = "C2", subcategory = "CP:KEGG")
      sets <- split(k$gene_symbol, k$gs_name)
      pathway_src <- "MSigDB_KEGG"
      message(sprintf("  KEGG gene sets: %d pathways", length(sets)))
      sets
    }, error = function(e) NULL)
  }

} else {
  message("  msigdbr not available — install with install.packages('msigdbr')")
  message("  Skipping pathway enrichment layer.")
}

fgsea_pathway <- NULL

if (!is.null(pathway_sets)) {

  fgsea_pathway <- tryCatch(
    fgsea::fgsea(
      pathways    = pathway_sets,
      stats       = ranked_vec,
      minSize     = 10,
      maxSize     = 500,
      nPermSimple = 10000,
      eps         = 0
    ),
    error = function(e) {
      message("  fgsea pathway enrichment failed: ", conditionMessage(e))
      NULL
    }
  )

  if (!is.null(fgsea_pathway) && nrow(fgsea_pathway) > 0) {

    fgsea_pathway <- as.data.frame(fgsea_pathway)
    fgsea_pathway$direction_in_RC <- ifelse(
      fgsea_pathway$NES > 0,
      "Enriched in RC", "Enriched in SC"
    )

    fgsea_sig <- fgsea_pathway[fgsea_pathway$padj < 0.05, ]
    message(sprintf("  Significant pathways (padj < 0.05): %d", nrow(fgsea_sig)))

    # Save: collapse leadingEdge list column
    save_df <- fgsea_pathway
    save_df$leadingEdge <- sapply(save_df$leadingEdge, paste, collapse = ";")
    write.csv(save_df,
              file.path(OUTPUT_DIR, "pathway_fgsea_all.csv"),
              row.names = FALSE)

    sig_df <- fgsea_sig
    sig_df$leadingEdge <- sapply(sig_df$leadingEdge, paste, collapse = ";")
    write.csv(sig_df,
              file.path(OUTPUT_DIR, "pathway_fgsea_significant.csv"),
              row.names = FALSE)
    message("  Pathway fgsea results saved.")

    # Pathway bubble plot
    if (nrow(fgsea_sig) >= 3) {
      top_paths <- fgsea_sig[order(fgsea_sig$NES), ]
      n_show    <- min(30, nrow(top_paths))
      top_paths <- rbind(head(top_paths, ceiling(n_show / 2)),
                         tail(top_paths, floor(n_show / 2)))
      top_paths$leadingEdge_str <- sapply(
        top_paths$leadingEdge,
        function(x) if (is.list(x)) paste(x[[1]], collapse = ";") else x
      )

      p_path <- ggplot2::ggplot(top_paths,
        ggplot2::aes(
          x     = NES,
          y     = reorder(pathway, NES),
          size  = size,
          color = direction_in_RC,
          fill  = direction_in_RC
        )) +
        ggplot2::geom_point(shape = 21, alpha = 0.85) +
        ggplot2::scale_color_manual(
          values = c("Enriched in RC" = "#B2182B",
                     "Enriched in SC" = "#2166AC"),
          name = ""
        ) +
        ggplot2::scale_fill_manual(
          values = c("Enriched in RC" = "#B2182B",
                     "Enriched in SC" = "#2166AC"),
          name = ""
        ) +
        ggplot2::scale_size_continuous(range = c(3, 9), name = "Gene set size") +
        ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
                            color = "grey40", linewidth = 0.4) +
        ggplot2::scale_x_continuous(
          labels = function(x) sprintf("%.1f", x)
        ) +
        ggplot2::labs(
          title    = sprintf("Pathway Enrichment (%s)", pathway_src),
          subtitle = "fgsea | ranked by log2FC(RC/SC) | positive NES = enriched in RC",
          x = "Normalized Enrichment Score (NES)",
          y = NULL
        ) +
        ggplot2::theme_bw(base_size = 10) +
        ggplot2::theme(
          plot.title  = ggplot2::element_text(face = "bold"),
          axis.text.y = ggplot2::element_text(size = 8)
        )

      ggplot2::ggsave(file.path(OUTPUT_DIR, "pathway_fgsea_bubble.pdf"),
                      p_path, width = 12, height = 9, device = cairo_pdf)
      message("  Pathway enrichment bubble plot saved.")
    }
  }
}


# =============================================================================
# P3.5  LAYER 3: Convergence — active TFs × enriched pathways
# =============================================================================
# For each significantly enriched pathway, ask: which of our top-active TFs
# has a regulon that overlaps the pathway's leading-edge genes?
# A TF that (a) is differentially active by ULM AND (b) has many targets among
# the genes driving pathway enrichment is the strongest mechanistic candidate.


conv_df <- NULL

if (!is.null(fgsea_pathway) && nrow(fgsea_sig) > 0) {

  # Top active TFs: |delta| in top 20%
  top_active_tfs <- tf_priority |> # deyisdim
    dplyr::filter(abs(delta_RC_vs_SC) >= delta_q80) |>
    dplyr::pull(TF)
  message(sprintf("  Top-active TFs (|delta| >= %.2f): %d",
                  delta_q80, length(top_active_tfs)))

  # Build regulon as a named list: TF → target genes
  tf_regulon_list <- split(regulon_net$target, regulon_net$source) ##duzelmelidi

  # For each significant pathway, find overlapping TFs
  conv_rows <- list()
  for (i in seq_len(nrow(fgsea_sig))) {
    pw_name    <- fgsea_sig$pathway[i]
    pw_NES     <- fgsea_sig$NES[i]
    pw_dir     <- fgsea_sig$direction_in_RC[i]

    pw_leading <- fgsea_sig$leadingEdge[[i]]
    if (is.list(pw_leading)) {
      pw_leading <- unlist(pw_leading, use.names = FALSE)
    } else if (is.character(pw_leading) && length(pw_leading) == 1L) {
      pw_leading <- strsplit(pw_leading, ";\\s*")[[1]]
    }
    pw_leading <- as.character(pw_leading)
    pw_leading <- pw_leading[!is.na(pw_leading) & pw_leading != ""]

    if (length(pw_leading) == 0L) next

    for (tf in top_active_tfs) {
      if (!tf %in% names(tf_regulon_list)) next
      tf_targets <- tf_regulon_list[[tf]]
      overlap    <- intersect(tf_targets, pw_leading)
      if (length(overlap) < 3) next ## duzelmelidi

      tf_info <- tf_priority[tf_priority$TF == tf, ]
      if (nrow(tf_info) == 0) next

      conv_rows[[length(conv_rows) + 1]] <- data.frame(
        TF                   = tf,
        pathway              = pw_name,
        pathway_NES          = pw_NES,
        pathway_direction    = pw_dir,
        pathway_padj         = fgsea_sig$padj[i],
        TF_delta_RC_vs_SC    = tf_info$delta_RC_vs_SC[1],
        TF_activity_direction= tf_info$activity_direction[1],
        TF_priority_score    = tf_info$priority_score[1],
        overlap_n            = length(overlap), ##duzelmelidi
        overlap_genes        = paste(overlap, collapse = ";"),
        iron_related         = tf_info$iron_related[1],
        mapk_akt             = tf_info$mapk_akt[1],
        stem_related         = tf_info$stem_related[1],
        is_DEG               = tf_info$is_DEG[1],
        is_hub               = tf_info$is_hub[1],
        stringsAsFactors     = FALSE
      )
    }
  }

  if (length(conv_rows) > 0) {
    conv_df <- dplyr::bind_rows(conv_rows) |>
      # Flag concordant: TF active in RC AND pathway enriched in RC (or both in SC)
      dplyr::mutate(
        concordant = (TF_delta_RC_vs_SC > 0 & pathway_NES > 0) |
                     (TF_delta_RC_vs_SC < 0 & pathway_NES < 0),
        # Combined score for ranking
        combined_score = abs(TF_delta_RC_vs_SC) * abs(pathway_NES) *
                         log1p(overlap_n) * (1 + concordant)
      ) |>
      dplyr::arrange(dplyr::desc(combined_score))

    write.csv(conv_df,
              file.path(OUTPUT_DIR, "TF_pathway_convergence.csv"),
              row.names = FALSE)
    message(sprintf("  Convergence table: %d TF-pathway pairs (%d concordant).",
                    nrow(conv_df), sum(conv_df$concordant)))

    # Also save a summary: unique TFs that appear in convergence, ranked
    tf_conv_summary <- conv_df |>
      dplyr::group_by(TF) |>
      dplyr::summarise(
        n_pathways       = dplyr::n(),
        n_concordant     = sum(concordant),
        top_pathway      = pathway[which.max(combined_score)],
        max_combined     = max(combined_score),
        delta_RC_vs_SC   = TF_delta_RC_vs_SC[1],
        activity_dir     = TF_activity_direction[1],
        priority_score   = TF_priority_score[1],
        iron_related     = iron_related[1],
        mapk_akt         = mapk_akt[1],
        stem_related     = stem_related[1],
        is_DEG           = is_DEG[1],
        is_hub           = is_hub[1],
        .groups          = "drop"
      ) |>
      dplyr::arrange(dplyr::desc(max_combined))

    write.csv(tf_conv_summary,
              file.path(OUTPUT_DIR, "TF_converged_final_candidates.csv"),
              row.names = FALSE)
    message(sprintf("  Final candidates: %d converged TFs saved.",
                    nrow(tf_conv_summary)))

    message("\n  Top 15 converged TF candidates:")
    print(head(tf_conv_summary[, c("TF","n_pathways","n_concordant",
                                   "delta_RC_vs_SC","activity_dir",
                                   "top_pathway","iron_related")], 15))

    # ── Convergence scatter plot ──────────────────────────────────────────────
    # x = TF activity delta, y = number of enriched pathways regulated,
    # size = max overlap genes, color = concordance rate

    conv_scatter <- tf_conv_summary |>
      dplyr::mutate(
        concordance_rate = n_concordant / n_pathways,
        label_tf         = ifelse(
          rank(-max_combined) <= 20 | iron_related | mapk_akt,
          TF, NA_character_)
      )

    p_conv <- ggplot2::ggplot(conv_scatter,
      ggplot2::aes(
        x     = delta_RC_vs_SC,
        y     = n_pathways,
        size  = max_combined,
        color = concordance_rate,
        label = label_tf
      )) +
      ggplot2::geom_point(alpha = 0.75) +
      ggrepel::geom_text_repel(na.rm = TRUE, size = 2.8, fontface = "italic",
                               max.overlaps = 25,
                               segment.size = 0.2, segment.color = "grey50") +
      ggplot2::scale_color_gradient2(
        low      = "#2166AC",
        mid      = "#F7F7F7",
        high     = "#B2182B",
        midpoint = 0.5,
        name     = "Concordance\nrate"
      ) +
      ggplot2::scale_size_continuous(range = c(2, 9), name = "Combined\nscore") +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
                          color = "grey40", linewidth = 0.4) +
      ggplot2::labs(
        title    = "Convergence: TF Activity × Pathway Enrichment",
        subtitle = paste0(
          "x-axis: TF activity delta (RC - SC, ULM)\n",
          "y-axis: number of enriched pathways with ≥3 regulon targets\n",
          "color: fraction of pathways where TF activity and pathway direction agree"
        ),
        x = "TF activity delta (RC - SC)",
        y = "Number of enriched pathways with regulon overlap"
      ) +
      ggplot2::theme_bw(base_size = 11) +
      ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"),
                     plot.subtitle = ggplot2::element_text(size = 8))

    ggplot2::ggsave(file.path(OUTPUT_DIR, "TF_convergence_scatter.pdf"),
                    p_conv, width = 10, height = 8, device = cairo_pdf)
    message("  Convergence scatter plot saved.")

  } else {
    message("  No TF-pathway pairs found with overlap >= 3 genes.")
    message("  Consider reducing the overlap threshold or the |delta| cutoff.")
  }

} else {
  message("  Pathway fgsea results not available — skipping convergence layer.")
  message("  Saving ULM-only candidates as final output.")

  # If no pathway layer, use ULM-only top TFs as the final candidate table
  final_candidates <- tf_priority |>
    dplyr::filter(abs(delta_RC_vs_SC) >= delta_q80) |>
    dplyr::select(TF, delta_RC_vs_SC, activity_direction, priority_score,
                  iron_related, mapk_akt, stem_related, is_DEG, is_hub) |>
    dplyr::arrange(dplyr::desc(priority_score))

  write.csv(final_candidates,
            file.path(OUTPUT_DIR, "TF_converged_final_candidates.csv"),
            row.names = FALSE)
  message(sprintf("  %d top-active TFs saved as final candidates.",
                  nrow(final_candidates)))
}


# =============================================================================
# SUMMARY
# =============================================================================

message("\n", paste(rep("=", 60), collapse = ""))
message("PRIORITY 3 COMPLETE")
message(paste(rep("=", 60), collapse = ""))
message(sprintf("  Regulon source:       %s", regulon_src))
message(sprintf("  TFs scored by ULM:   %d", nrow(ulm_wide)))
message(sprintf("  Pathway gene sets:   %s", pathway_src))
if (!is.null(fgsea_pathway))
  message(sprintf("  Significant pathways: %d (padj < 0.05)",
                  sum(fgsea_pathway$padj < 0.05, na.rm = TRUE)))
if (!is.null(conv_df))
  message(sprintf("  TF-pathway pairs:     %d (%d concordant)",
                  nrow(conv_df), sum(conv_df$concordant)))
message("\n  Output files:")
message("    TF_activity_scores.csv          — full ULM scores for all TFs")
message("    TF_ULM_priority.csv             — ULM results with priority scoring")
message("    TF_activity_heatmap.pdf         — top 40 TFs by |delta|")
message("    TF_ULM_bubble.pdf               — bubble plot top 30 TFs")
message("    pathway_fgsea_all.csv           — all pathway fgsea results")
message("    pathway_fgsea_significant.csv   — significant pathways (padj<0.05)")
message("    pathway_fgsea_bubble.pdf        — pathway enrichment plot")
message("    TF_pathway_convergence.csv      — all TF-pathway pairs")
message("    TF_converged_final_candidates.csv — final ranked TF candidates\n")
# =============================================================================
# FINAL SUMMARY
# =============================================================================

message("\n", paste(rep("=", 60), collapse = ""))
message("ANALYSIS COMPLETE")
message(paste(rep("=", 60), collapse = ""))
message(sprintf("\nOutput files in: %s/", OUTPUT_DIR))
message("\n[Stage 0-1] Data & DEGs:")
message(sprintf("  DEGs SC→RC:   %d genes", nrow(deg1)))
message(sprintf("  DEGs SC→S10:  %d genes", nrow(deg2)))
message(sprintf("  DEGs RC→R10:  %d genes", nrow(deg3)))
message(sprintf("  miRNA DEGs:   %d (SC→RC)", nrow(mirna_deg1)))
message("\n[Stage 2] Visualizations:")
message("  - Volcano plots (3 comparisons)")
message("  - Annotated pathway heatmap")
message("  - Iron gene expression barplot")
message("  - miRNA DEG heatmap")
# Example for GSE45558
gse <- getGEO("GSE45558")[[1]]
metadata <- pData(gse)
braf_only <- metadata[grep("Vemurafenib", metadata$title, ignore.case = TRUE), ]
library("GEOquery")
# Function to get BRAFi-only samples from GSE45558
get_gse45558_samples <- function() {
  gse <- getGEO("GSE45558", getGPL = FALSE)[[1]]
  meta <- pData(gse)
  
  # Filter for Parental (non-resistant) lines treated with Vemurafenib
  # Note: Resistance studies often include 'R' (Resistant) and 'P' (Parental)
  braf_samples <- meta[grepl("vemurafenib", meta$title, ignore.case = TRUE) & 
                        !grepl("-R", meta$title), ]
  return(list(expr = exprs(gse)[, rownames(braf_samples)], meta = braf_samples))
}
# Function to get BRAFi-only samples from GSE186108
get_gse186108_samples <- function() {
  gse <- getGEO("GSE186108", getGPL = FALSE)[[1]]
  meta <- pData(gse)
  
  # Filtering logic: Include Encorafenib, Exclude Binimetinib (Combo therapy)
  # In GSE186108, treatment is often in 'characteristics_ch1'
  braf_only18 <- meta[grepl("encorafenib", meta$characteristics_ch1, ignore.case = TRUE) & 
                    !grepl("binimetinib", meta$characteristics_ch1, ignore.case = TRUE), ]
  
  return(list(expr = exprs(gse)[, rownames(braf_only18)], meta = braf_only18))
}