# =============================================================================
# 04_visualizations.R
# INTEGRATIVE ANALYSIS: ENCORAFENIB RESISTANCE IN MALIGNANT MELANOMA
# GSE283251 | Colakoglu Bergel et al., Scientific Reports 2025
#
# Input:   checkpoint_03.RData
# Output:  figures in newoutput/
#
# This script keeps the original stage-2 visualization block with only the
# minimum fixes needed for a clean standalone run.
# =============================================================================


# ── Configuration ─────────────────────────────────────────────────────────────

DATA_DIR   <- "/home/ramin/mela/data/ilk"  # <-- UPDATE THIS PATH
OUTPUT_DIR <- file.path(DATA_DIR, "newoutput")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)


# ── Packages ──────────────────────────────────────────────────────────────────

for (pkg in c("dplyr", "tidyr", "ggplot2", "ggrepel", "pheatmap",
              "RColorBrewer", "scales", "patchwork")) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(RColorBrewer)
  library(scales)
  library(patchwork)
})


# ── Load checkpoint ───────────────────────────────────────────────────────────

load(file.path(DATA_DIR, "checkpoint_02.RData"))


# ── Gene categories ──────────────────────────────────────────────────────────

IRON_GENES <- c("NCOA4", "FTH1", "TFRC", "SLC7A11", "GPX4", "IREB2",
                "FTL", "HAMP", "SLC40A1", "HMOX1", "CYBRD1", "STEAP3")
AUTOPHAGY_GENES <- c("BECN1", "MAP1LC3B", "ATG5", "ATG7", "ATG12",
                     "SQSTM1", "ULK1", "WIPI2", "LAMP2", "RAB7A")
MAPK_AKT_GENES  <- c("BRAF", "KRAS", "NRAS", "MAP2K1", "MAPK1", "MAPK3",
                     "AKT1", "AKT2", "PIK3CA", "PIK3R1", "PTEN",
                     "DUSP6", "DUSP4", "RAF1", "CRAF")
MDR_GENES       <- c("ABCB1", "YBX1", "MDR1")

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


# =============================================================================
# 2.1 Volcano plots
# =============================================================================

make_volcano <- function(deg_df, title, label_genes = NULL,
                         fc_threshold = 1, prob_threshold = 0.8) {
  df <- as.data.frame(deg_df$full)
  df$gene <- rownames(df)
  df$log2FC <- df$M

  df$sig <- "NS"
  df$sig[df$prob >= prob_threshold & df$log2FC >  fc_threshold] <- "Up"
  df$sig[df$prob >= prob_threshold & df$log2FC < -fc_threshold] <- "Down"

  df$category <- dplyr::case_when(
    df$gene %in% IRON_GENES      ~ "Iron/Ferroptosis",
    df$gene %in% AUTOPHAGY_GENES ~ "Autophagy",
    df$gene %in% MAPK_AKT_GENES  ~ "MAPK/AKT",
    df$gene %in% MDR_GENES       ~ "MDR",
    df$sig != "NS"               ~ "Other DEG",
    TRUE                         ~ "NS"
  )

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
      x        = "log2 Fold Change (B vs A)",
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

p_vol1 <- make_volcano(comp1, "Resistant vs Sensitive (Control)\nSC -> RC")
p_vol2 <- make_volcano(comp2, "Drug Effect in Sensitive\nSC -> S10 (10 nM)")
p_vol3 <- make_volcano(comp3, "Drug Effect in Resistant\nRC -> R10 (10 nM)")
p_vol4 <- make_volcano(comp4, "Drug Effect Sensitive vs Resistant\nS10 -> R10 (10 nM)")

ggsave(file.path(OUTPUT_DIR, "volcano_SC_vs_RC.pdf"),  p_vol1, width = 9, height = 7)
ggsave(file.path(OUTPUT_DIR, "volcano_SC_vs_S10.pdf"), p_vol2, width = 9, height = 7)
ggsave(file.path(OUTPUT_DIR, "volcano_RC_vs_R10.pdf"), p_vol3, width = 9, height = 7)
ggsave(file.path(OUTPUT_DIR, "volcano_S10_vs_R10.pdf"), p_vol4, width = 9, height = 7)
message("  Volcano plots saved.")


# =============================================================================
# 2.2 Annotated heatmap: key pathway genes
# =============================================================================

make_pathway_heatmap <- function(count_mat, gene_sets, title, filename) {
  all_genes <- unlist(gene_sets)
  present   <- all_genes[all_genes %in% rownames(count_mat)]

  if (length(present) == 0) {
    message(sprintf("  No genes found for heatmap: %s", filename))
    return(invisible(NULL))
  }

  cpm_mat <- sweep(count_mat[present, , drop = FALSE], 2,
                   colSums(count_mat) / 1e6, FUN = "/")
  log_mat <- log2(cpm_mat + 1)
  z_mat   <- t(scale(t(log_mat)))
  z_mat[is.nan(z_mat)] <- 0

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
    Condition = c("Resistant", "Resistant", "Sensitive", "Sensitive"),
    Treatment = c("Control", "Drug", "Control", "Drug"),
    row.names = colnames(count_mat)
  )
  annot_colors$Condition <- c(Resistant = "#CC0000", Sensitive = "#0066CC")
  annot_colors$Treatment <- c(Control = "#999999", Drug = "#FF9900")

  pheatmap::pheatmap(
    z_mat,
    color             = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
    annotation_row    = row_annot,
    annotation_col    = col_annot,
    annotation_colors = annot_colors,
    cluster_cols      = FALSE,
    cluster_rows      = TRUE,
    show_rownames     = TRUE,
    fontsize_row      = 8,
    fontsize_col      = 9,
    border_color      = NA,
    main              = title,
    filename          = file.path(OUTPUT_DIR, filename),
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


# =============================================================================
# 2.3 Iron gene expression barplot
# =============================================================================

iron_genes_plot <- c("FTH1", "NCOA4", "TFRC", "SLC7A11", "GPX4", "IREB2",
                     "FTL", "HMOX1", "STEAP3")
iron_present <- iron_genes_plot[iron_genes_plot %in% rownames(mrna_mat)]

iron_plot_df <- as.data.frame(mrna_mat[iron_present, , drop = FALSE]) |>
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
                       levels = c("Sensitive\nControl", "Sensitive\n10 nM",
                                  "Resistant\nControl", "Resistant\n10 nM")),
    log2_cpm = log2(count / (sum(mrna_mat[, sample]) / 1e6) + 1)
  )

p_iron_bar <- ggplot(iron_plot_df, aes(x = condition, y = log2_cpm, fill = condition)) +
  geom_col(width = 0.7, color = "white", linewidth = 0.3) +
  facet_wrap(~ gene, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = c(
    "Sensitive\nControl" = "#6BAED6",
    "Sensitive\n10 nM"   = "#2171B5",
    "Resistant\nControl"  = "#FC8D59",
    "Resistant\n10 nM"    = "#D7301F"
  )) +
  labs(
    title = "Iron Metabolism Gene Expression",
    subtitle = "log2(CPM + 1) across all conditions",
    x = NULL, y = "log2(CPM + 1)", fill = "Condition"
  ) +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x      = element_text(size = 7),
    strip.background = element_rect(fill = "#f0f0f0"),
    strip.text       = element_text(face = "bold", size = 9),
    legend.position  = "none",
    plot.title       = element_text(face = "bold")
  )

ggsave(file.path(OUTPUT_DIR, "iron_gene_barplot.pdf"),
       p_iron_bar, width = 10, height = 8)
message("  Iron gene barplot saved.")


# =============================================================================
# 2.4 miRNA heatmap
# =============================================================================

paper_mirna_ids <- c(
  "hsa-miR-331-3p", "hsa-let-7c-5p",  "hsa-miR-296-5p", "hsa-let-7b-5p",
  "hsa-miR-31-5p",  "hsa-miR-1260b",  "hsa-let-7a-5p",  "hsa-miR-1229-3p",
  "hsa-miR-6516-3p", "hsa-miR-744-5p", "hsa-miR-181a-5p", "hsa-miR-34a-5p",
  "hsa-miR-103a-3p", "hsa-miR-140-3p"
)
mirna_present <- paper_mirna_ids[paper_mirna_ids %in% rownames(mirna_mat)]

if (length(mirna_present) > 0) {
  mirna_heatmap_mat <- mirna_cpm[mirna_present, , drop = FALSE]
  mirna_z <- t(scale(t(log2(mirna_heatmap_mat + 1))))
  mirna_z[is.nan(mirna_z)] <- 0

  mirna_row_annot <- data.frame(
    Role = rep("RC_vs_SC_DEG", length(mirna_present)),
    row.names = mirna_present
  )
  if ("hsa-miR-140-3p" %in% rownames(mirna_row_annot)) {
    mirna_row_annot["hsa-miR-140-3p", "Role"] <- "Iron-regulatory"
  }

  pheatmap::pheatmap(
    mirna_z,
    color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
    annotation_col = data.frame(
      Condition = c("Resistant", "Resistant", "Sensitive", "Sensitive"),
      Treatment = c("Control", "Drug", "Control", "Drug"),
      row.names = colnames(mirna_z)
    ),
    annotation_row = mirna_row_annot,
    cluster_cols   = FALSE,
    cluster_rows   = TRUE,
    fontsize_row   = 9,
    main           = "DEG miRNAs (Paper Fig 5a)\nZ-scored log2CPM",
    filename       = file.path(OUTPUT_DIR, "heatmap_mirna_paper_DEGs.pdf"),
    width = 7, height = 6
  )
  message("  miRNA heatmap saved.")
}

