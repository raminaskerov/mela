# =============================================================================
# 06_TF_enrichment.R
# INTEGRATIVE ANALYSIS: ENCORAFENIB RESISTANCE IN MALIGNANT MELANOMA
#
# Input:   checkpoint_03.RData
# Output:  TF fgsea tables, TF activity tables, checkpoint_06.RData
#
# This stage is independent of WGCNA outputs.
# =============================================================================


# ── Configuration ─────────────────────────────────────────────────────────────

DATA_DIR   <- "/home/ramin/mela/data/ilk"  # <-- UPDATE THIS PATH
OUTPUT_DIR <- file.path(DATA_DIR, "newoutput")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)


# ── Packages ──────────────────────────────────────────────────────────────────

cran_pkgs <- c("dplyr", "ggplot2", "ggrepel", "pheatmap", "tidyr")
bioc_pkgs  <- c("fgsea", "decoupleR", "org.Hs.eg.db")

for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
for (pkg in bioc_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
}

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(tidyr)
  library(fgsea)
  library(decoupleR)
  library(org.Hs.eg.db)
})


# ── Load checkpoint ───────────────────────────────────────────────────────────

load(file.path(DATA_DIR, "checkpoint_03.RData"))


# ── Helpers ───────────────────────────────────────────────────────────────────

save_fgsea_table <- function(result_df, filepath) {
  df <- as.data.frame(result_df)
  if ("leadingEdge" %in% colnames(df)) {
    df$leadingEdge <- vapply(
      df$leadingEdge,
      function(x) paste(x, collapse = ";"),
      character(1)
    )
  }
  write.csv(df, filepath, row.names = FALSE)
  invisible(df)
}

build_ranked_vec <- function(full_df) {
  ranked_df <- as.data.frame(full_df)
  ranked_df$gene <- rownames(ranked_df)
  ranked_df <- ranked_df |>
    dplyr::filter(!is.na(M)) |>
    dplyr::arrange(dplyr::desc(M)) |>
    dplyr::distinct(gene, .keep_all = TRUE)
  ranked_vec <- ranked_df$M
  names(ranked_vec) <- ranked_df$gene
  ranked_vec
}


# =============================================================================
# TF REGULON GSEA
# =============================================================================

message("\n=== PRIORITY 3: TF enrichment ===")

ranked_vec <- build_ranked_vec(comp1$full)
message(sprintf("  Ranked gene list: %d genes | range [%.2f, %.2f]",
                length(ranked_vec), min(ranked_vec), max(ranked_vec)))

tf_regulons <- NULL
regulon_source <- "none"

tf_regulons <- tryCatch({
  net <- decoupleR::get_collectri(organism = "human", split_complexes = FALSE)
  regulon_source <- "CollecTRI"
  message(sprintf("  CollecTRI loaded: %d TF-target pairs", nrow(net)))
  net
}, error = function(e) {
  message("  CollecTRI failed: ", conditionMessage(e))
  NULL
})

if (is.null(tf_regulons)) {
  tf_regulons <- tryCatch({
    net <- decoupleR::get_dorothea(organism = "human", levels = c("A", "B", "C"))
    regulon_source <- "DoRothEA_ABC"
    message(sprintf("  DoRothEA loaded: %d TF-target pairs", nrow(net)))
    net
  }, error = function(e) {
    message("  DoRothEA failed: ", conditionMessage(e))
    NULL
  })
}

if (is.null(tf_regulons)) {
  message("  No regulon database available; falling back to a manual TF set.")
  regulon_source <- "manual"
  tf_regulons <- data.frame(
    source = rep(c("NFE2L2", "STAT3", "HIF1A", "SP1", "MYC", "TP53",
                   "MITF", "ZEB1", "BACH1", "TFEB"), each = 6),
    target = c(
      "HMOX1", "NQO1", "GCLC", "GCLM", "SLC7A11", "FTH1",
      "BCL2", "MCL1", "CCND1", "MYC", "VEGFA", "TWIST1",
      "LDHA", "VEGFA", "SLC2A1", "HMOX1", "BNIP3", "SLC7A11",
      "TFRC", "FTH1", "FTL", "HMOX1", "SLC11A2", "NCOA4",
      "CDK4", "CCND2", "E2F1", "LDHA", "PKM", "MCM2",
      "CDKN1A", "BAX", "MDM2", "GADD45A", "TIGAR", "SCO2",
      "DCT", "TYRP1", "TYR", "CDK2", "BCL2", "MLANA",
      "VIM", "FN1", "CDH2", "SNAI1", "SNAI2", "TWIST1",
      "HMOX1", "FTH1", "FTL", "BLVRB", "SLC7A11", "SQSTM1",
      "LAMP1", "LAMP2", "CTSD", "SQSTM1", "BECN1", "ATG5"
    ),
    mor = 1,
    stringsAsFactors = FALSE
  )
}

tf_sets_all <- split(tf_regulons$target, tf_regulons$source)
tf_sets_all <- tf_sets_all[lengths(tf_sets_all) >= 5]

message(sprintf("  TF gene sets: %d TFs with >= 5 targets", length(tf_sets_all)))

fgsea_tf <- tryCatch(
  fgsea::fgsea(
    pathways    = tf_sets_all,
    stats       = ranked_vec,
    minSize     = 5,
    maxSize     = 500,
    nPermSimple = 10000,
    eps         = 0
  ),
  error = function(e) {
    message("  fgsea failed: ", conditionMessage(e))
    NULL
  }
)

fgsea_tf_sig <- NULL
tf_summary <- NULL
conv_df <- NULL

if (!is.null(fgsea_tf) && nrow(fgsea_tf) > 0) {
  fgsea_tf <- as.data.frame(fgsea_tf)
  fgsea_tf$direction_in_RC <- ifelse(
    fgsea_tf$NES > 0, "More active in RC", "Less active in RC"
  )

  fgsea_tf_sig <- fgsea_tf[fgsea_tf$padj < 0.05, ]
  if (nrow(fgsea_tf_sig) == 0) {
    fgsea_tf_sig <- fgsea_tf[fgsea_tf$pval < 0.05, ]
    message("  No TFs at padj<0.05; using nominal p<0.05")
  }
  fgsea_tf_sig <- fgsea_tf_sig[order(fgsea_tf_sig$NES, decreasing = TRUE), ]

  save_fgsea_table(fgsea_tf, file.path(OUTPUT_DIR, "TF_fgsea_all.csv"))
  save_fgsea_table(fgsea_tf_sig, file.path(OUTPUT_DIR, "TF_fgsea_significant.csv"))
  message(sprintf("  Significant TFs: %d", nrow(fgsea_tf_sig)))

  hub_genes_vec <- if (exists("hub_genes") && length(hub_genes) > 0) {
    unlist(hub_genes, use.names = FALSE)
  } else {
    character(0)
  }

  tf_summary <- fgsea_tf_sig |>
    dplyr::rename(TF_name = pathway) |>
    dplyr::mutate(
      is_DEG       = TF_name %in% deg1$gene,
      is_hub       = TF_name %in% hub_genes_vec,
      iron_related = TF_name %in% c("NFE2L2", "NRF2", "BACH1", "HIF1A", "MTF1",
                                    "TFEB", "STAT3", "SP1", "HMOX1", "TF"),
      mapk_akt     = TF_name %in% c("BRAF", "KRAS", "NRAS", "MAP2K1", "MAPK1",
                                    "MAPK3", "AKT1", "AKT2", "PIK3CA", "PIK3R1",
                                    "PTEN", "DUSP6", "DUSP4", "RAF1", "CRAF"),
      stem_related = TF_name %in% c("POU5F1", "SOX2", "NANOG", "KLF4", "MYC",
                                    "CD44", "PROM1", "ALDH1A1", "ABCB5", "KDM5B",
                                    "ZEB1", "TWIST1", "SNAI1", "SNAI2", "VIM",
                                    "NOTCH1", "WNT5A", "AXL", "MITF"),
      priority_score = abs(NES) +
        (direction_in_RC == "More active in RC") * 1.0 +
        (is_DEG * 1.5) +
        (is_hub * 1.5) +
        (iron_related * 2.5) +
        (mapk_akt * 1.5) +
        (stem_related * 1.0)
    ) |>
    dplyr::arrange(dplyr::desc(priority_score))

  save_fgsea_table(tf_summary, file.path(OUTPUT_DIR, "TF_ULM_priority.csv"))

  top25 <- head(tf_summary, 25)
  top25$pathway_cat <- dplyr::case_when(
    top25$iron_related ~ "Iron/Ferroptosis",
    top25$mapk_akt     ~ "MAPK/AKT",
    top25$stem_related ~ "Stemness",
    TRUE               ~ "Other"
  )

  p_tf <- ggplot(top25, aes(
    x = NES,
    y = reorder(TF_name, NES),
    size = size,
    color = pathway_cat,
    shape = direction_in_RC
  )) +
    geom_point(alpha = 0.85) +
    scale_color_manual(
      values = c("Iron/Ferroptosis" = "#E64B35", "MAPK/AKT" = "#00A087",
                 "Stemness" = "#7B2D8B", "Other" = "#8491B4"),
      name = "Pathway"
    ) +
    scale_shape_manual(
      values = c("More active in RC" = 17L, "Less active in RC" = 16L),
      name = "Activity in RC"
    ) +
    scale_size_continuous(range = c(3, 10), name = "Regulon size") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.4) +
    labs(
      title = "TF Activity in Resistant vs Sensitive Cells",
      subtitle = sprintf("fgsea on %s regulons | ranked by log2FC(RC/SC)", regulon_source),
      x = "Normalized Enrichment Score (NES)",
      y = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(face = "bold"), legend.position = "right")

  ggsave(file.path(OUTPUT_DIR, "TF_fgsea_bubble.pdf"), p_tf, width = 11, height = 8)
} else {
  message("  fgsea returned no TF results.")
}


# =============================================================================
# TF ACTIVITY SCORING
# =============================================================================

tf_expr_mat <- log2(mrna_norm + 1)
tf_act_df <- NULL

regulon_net <- tryCatch(
  decoupleR::get_collectri(organism = "human", split_complexes = FALSE),
  error = function(e) {
    message("  CollecTRI activity regulon unavailable: ", conditionMessage(e))
    tryCatch(
      decoupleR::get_dorothea(organism = "human", levels = c("A", "B", "C")),
      error = function(e2) {
        message("  DoRothEA activity regulon unavailable: ", conditionMessage(e2))
        NULL
      }
    )
  }
)

if (!is.null(regulon_net)) {
  message(sprintf("  Activity regulon loaded: %d TF-target interactions", nrow(regulon_net)))

  tf_acts_long <- tryCatch(
    decoupleR::run_ulm(
      mat = as.matrix(tf_expr_mat),
      net = regulon_net,
      .source = "source",
      .target = "target",
      .mor = "mor",
      minsize = 5L
    ),
    error = function(e) {
      message("  decoupleR::run_ulm failed: ", conditionMessage(e))
      NULL
    }
  )

  if (!is.null(tf_acts_long)) {
    tf_act_wide <- tf_acts_long |>
      dplyr::filter(statistic == "ulm") |>
      dplyr::select(source, condition, score) |>
      tidyr::pivot_wider(names_from = condition, values_from = score)

    if (is.data.frame(tf_act_wide) && nrow(tf_act_wide) > 0) {
      tf_act_df <- as.data.frame(tf_act_wide)
      rownames(tf_act_df) <- tf_act_df$source
      tf_act_df$TF <- tf_act_df$source
      tf_act_df$source <- NULL
      message(sprintf("  TF activity profiles computed: %d TFs", nrow(tf_act_df)))
    }
  }
}

if (is.null(tf_act_df) || nrow(tf_act_df) == 0) {
  message("  TF activity scoring unavailable.")
} else {
  sample_cols_present <- intersect(
    c("A375_RC", "A375_R10", "A375_SC", "A375_S10"),
    colnames(tf_act_df)
  )

  if (all(c("A375_RC", "A375_SC") %in% sample_cols_present)) {
    tf_act_df$delta_RC_vs_SC <- tf_act_df$A375_RC - tf_act_df$A375_SC
  } else if (length(sample_cols_present) >= 2) {
    tf_act_df$delta_RC_vs_SC <- tf_act_df[[sample_cols_present[1]]] - tf_act_df[[sample_cols_present[2]]]
  }

  tf_act_df <- tf_act_df |>
    dplyr::mutate(
      iron_related = TF %in% c("NFE2L2", "NRF2", "BACH1", "HIF1A", "MTF1",
                               "TFEB", "STAT3", "SP1", "HMOX1"),
      mapk_akt     = TF %in% c("BRAF", "KRAS", "NRAS", "MAP2K1", "MAPK1",
                               "MAPK3", "AKT1", "AKT2", "PIK3CA", "PIK3R1",
                               "PTEN", "DUSP6", "DUSP4", "RAF1", "CRAF"),
      stem_related = TF %in% c("POU5F1", "SOX2", "NANOG", "KLF4", "MYC",
                               "CD44", "PROM1", "ALDH1A1", "ABCB5", "KDM5B",
                               "ZEB1", "TWIST1", "SNAI1", "SNAI2", "VIM",
                               "NOTCH1", "WNT5A", "AXL", "MITF")
    )

  if ("delta_RC_vs_SC" %in% colnames(tf_act_df)) {
    tf_act_df <- dplyr::arrange(tf_act_df, dplyr::desc(abs(delta_RC_vs_SC)))
  }

  write.csv(tf_act_df, file.path(OUTPUT_DIR, "TF_activity_scores.csv"), row.names = FALSE)
  message(sprintf("  TF activity scores saved: %d TFs", nrow(tf_act_df)))

  if ("delta_RC_vs_SC" %in% colnames(tf_act_df) && length(sample_cols_present) >= 2) {
    q90 <- quantile(abs(tf_act_df$delta_RC_vs_SC), 0.90, na.rm = TRUE)
    top_act <- tf_act_df[abs(tf_act_df$delta_RC_vs_SC) >= q90, ]
    top_act <- head(top_act[order(abs(top_act$delta_RC_vs_SC), decreasing = TRUE), ], 40)

    if (nrow(top_act) >= 5) {
      act_mat <- as.matrix(top_act[, sample_cols_present, drop = FALSE])
      rownames(act_mat) <- top_act$TF

      act_row_ann <- data.frame(
        Pathway = dplyr::case_when(
          top_act$iron_related ~ "Iron/Ferroptosis",
          top_act$mapk_akt     ~ "MAPK/AKT",
          top_act$stem_related ~ "Stemness",
          TRUE                 ~ "Other"
        ),
        Direction = ifelse(top_act$delta_RC_vs_SC > 0, "More active RC", "Less active RC"),
        row.names = top_act$TF
      )

      cond_lookup <- c(A375_RC = "Resistant", A375_R10 = "Resistant",
                       A375_SC = "Sensitive", A375_S10 = "Sensitive")
      trt_lookup <- c(A375_RC = "Control", A375_R10 = "Drug",
                      A375_SC = "Control", A375_S10 = "Drug")
      act_col_ann <- data.frame(
        Condition = cond_lookup[sample_cols_present],
        Treatment = trt_lookup[sample_cols_present],
        row.names = sample_cols_present
      )

      act_annot_col <- list(
        Pathway = c("Iron/Ferroptosis" = "#E64B35", "MAPK/AKT" = "#00A087",
                    "Stemness" = "#7B2D8B", "Other" = "#8491B4"),
        Direction = c("More active RC" = "#B2182B", "Less active RC" = "#2166AC"),
        Condition = c(Resistant = "#CC0000", Sensitive = "#0066CC"),
        Treatment = c(Control = "#999999", Drug = "#FF9900")
      )

      pheatmap::pheatmap(
        act_mat,
        color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
        annotation_row = act_row_ann,
        annotation_col = act_col_ann,
        annotation_colors = act_annot_col,
        cluster_cols = FALSE,
        cluster_rows = TRUE,
        scale = "row",
        fontsize_row = 7,
        main = "TF Activity Scores (ULM | z-scaled per TF)\nTop differentially active TFs: RC vs SC",
        filename = file.path(OUTPUT_DIR, "TF_activity_heatmap.pdf"),
        width = 7,
        height = 10
      )
    }
  }

  if (!is.null(tf_summary) && nrow(tf_summary) > 0 && "delta_RC_vs_SC" %in% colnames(tf_act_df)) {
    q90_act <- quantile(abs(tf_act_df$delta_RC_vs_SC), 0.90, na.rm = TRUE)
    top_act_tfs <- tf_act_df$TF[abs(tf_act_df$delta_RC_vs_SC) >= q90_act]
    converged <- intersect(tf_summary$TF_name, top_act_tfs)

    if (length(converged) > 0) {
      conv_df <- dplyr::left_join(
        tf_summary[tf_summary$TF_name %in% converged, ],
        tf_act_df[tf_act_df$TF %in% converged, c("TF", "delta_RC_vs_SC")],
        by = c("TF_name" = "TF")
      ) |>
        dplyr::mutate(
          activity_direction = dplyr::case_when(
            !is.na(delta_RC_vs_SC) & delta_RC_vs_SC > 0 ~ "More active in RC",
            !is.na(delta_RC_vs_SC) & delta_RC_vs_SC < 0 ~ "Less active in RC",
            TRUE                                         ~ "Not computed"
          )
        ) |>
        dplyr::arrange(dplyr::desc(priority_score))

      write.csv(conv_df, file.path(OUTPUT_DIR, "TF_converged_final_candidates.csv"), row.names = FALSE)

      if (length(converged) >= 3) {
        plot_conv <- dplyr::inner_join(
          tf_summary[, c("TF_name", "priority_score", "iron_related", "stem_related")],
          tf_act_df[, c("TF", "delta_RC_vs_SC")],
          by = c("TF_name" = "TF")
        ) |>
          dplyr::filter(!is.na(delta_RC_vs_SC)) |>
          dplyr::mutate(
            is_converged = TF_name %in% converged,
            label_gene = ifelse(is_converged, TF_name, NA_character_)
          )

        p_conv <- ggplot(plot_conv, aes(
          x = delta_RC_vs_SC,
          y = priority_score,
          color = is_converged,
          label = label_gene
        )) +
          geom_point(alpha = 0.6, size = 2) +
          ggrepel::geom_text_repel(na.rm = TRUE, size = 2.8, fontface = "italic",
                                   max.overlaps = 20) +
          scale_color_manual(
            values = c(`TRUE` = "#E64B35", `FALSE` = "grey60"),
            labels = c("Other", "Converged"),
            name = ""
          ) +
          geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.4) +
          labs(
            title = "TF Activity vs fgsea Priority Score",
            subtitle = "Red = TFs significant in both layers",
            x = "TF activity delta (RC - SC | NES)",
            y = "fgsea priority score"
          ) +
          theme_bw(base_size = 11) +
          theme(plot.title = element_text(face = "bold"))

        ggsave(file.path(OUTPUT_DIR, "TF_convergence_scatter.pdf"), p_conv, width = 8, height = 6)
      }
    }
  }
}


# =============================================================================
# CHECKPOINT SAVE
# =============================================================================

save(
  mrna_mat, mrna_filtered, mrna_norm,
  mirna_mat, mirna_filtered, mirna_cpm,
  sample_meta,
  comp1, comp2, comp3, comp4,
  mirna_comp1, mirna_comp2, mirna_comp3, mirna_comp4,
  deg1, deg2, deg3, deg4,
  mirna_deg1, mirna_deg2, mirna_deg3, mirna_deg4,
  ranked_vec, tf_regulons, regulon_source,
  fgsea_tf, fgsea_tf_sig, tf_summary,
  tf_act_df, conv_df,
  file = file.path(DATA_DIR, "checkpoint_06.RData")
)

message("\ncheckpoint_06.RData saved.")
