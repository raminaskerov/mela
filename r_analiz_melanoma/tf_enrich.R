# =============================================================================
# PRIORITY 3: TRANSCRIPTION FACTOR ENRICHMENT ANALYSIS
# =============================================================================
#
# =============================================================================
# PRIORITY 3 (revised): TF enrichment via GSEA on full ranked gene list
# Replaces the enrichR ORA approach with fgsea against TF regulon databases
# =============================================================================


if (!requireNamespace("fgsea", quietly = TRUE)) BiocManager::install("fgsea")
library(fgsea)
library(org.Hs.eg.db)

# ── Helper: save fgsea result (handles list column) ──────────────────────────

save_fgsea <- function(result_df, filepath) {
  df <- as.data.frame(result_df)
  if ("leadingEdge" %in% colnames(df))
    df$leadingEdge <- sapply(df$leadingEdge, paste, collapse = ";")
  write.csv(df, filepath, row.names = FALSE)
}

# ── P3.0  Build ranked gene list from SC vs RC full NOISeq results ────────────
# Uses ALL expressed genes (not just DEGs) — required for valid GSEA
# M column = log2FC(A/B) = log2(SC/RC); positive = higher in SC, negative = higher in RC

full_results_comp1 <- as.data.frame(comp1$result@results[[1]])
full_results_comp1$gene <- rownames(full_results_comp1)

# Score: M (log2FC) — keeps direction and magnitude
# Remove genes with NA scores
ranked_genes <- full_results_comp1 |>
  dplyr::filter(!is.na(M)) |>
  dplyr::arrange(dplyr::desc(M)) |>
  dplyr::distinct(gene, .keep_all = TRUE)

# Named numeric vector: value = log2FC(SC/RC), name = gene symbol
# Negative values = higher in RC (resistant)
ranked_vec <- setNames(ranked_genes$M, ranked_genes$gene)
message(sprintf("  Ranked gene list: %d genes", length(ranked_vec)))
message(sprintf("  Range: %.2f to %.2f", min(ranked_vec), max(ranked_vec)))

# ── P3.1  Build TF regulon gene sets ─────────────────────────────────────────
# Three sources, each tried in order:
#   (A) CollecTRI via decoupleR   — broadest, most current
#   (B) DoRothEA via decoupleR    — confidence A+B+C
#   (C) Manual curated set        — offline fallback

tf_regulons <- NULL
regulon_source <- "none"

# Path A: CollecTRI
if (requireNamespace("decoupleR", quietly = TRUE)) {
  tf_regulons <- tryCatch({
    net <- decoupleR::get_collectri(organism = "human", split_complexes = FALSE)
    message(sprintf("  CollecTRI loaded: %d TF-target pairs", nrow(net)))
    regulon_source <- "CollecTRI"
    net
  }, error = function(e) {
    message("  CollecTRI failed: ", conditionMessage(e))
    NULL
  })
}

# Path B: DoRothEA via decoupleR
if (is.null(tf_regulons) && requireNamespace("decoupleR", quietly = TRUE)) {
  tf_regulons <- tryCatch({
    net <- decoupleR::get_dorothea(organism = "human", levels = c("A","B","C"))
    message(sprintf("  DoRothEA loaded: %d TF-target pairs", nrow(net)))
    regulon_source <- "DoRothEA"
    net
  }, error = function(e) {
    message("  DoRothEA failed: ", conditionMessage(e))
    NULL
  })
}

# Convert to named list for fgsea: each TF → character vector of target genes
# For CollecTRI/DoRothEA: source = TF, target = gene, mor = mode of regulation (+1/-1)
# We split into activation and repression sets for directional interpretation

if (!is.null(tf_regulons)) {
  
  # Full regulon (all targets regardless of direction)
  tf_sets_all <- split(tf_regulons$target, tf_regulons$source)
  # Keep only TFs with enough targets (fgsea minSize handles this but pre-filter saves time)
  tf_sets_all <- tf_sets_all[lengths(tf_sets_all) >= 5]
  message(sprintf("  TF gene sets: %d TFs with >= 5 targets", length(tf_sets_all)))
  
  # Activation targets only (mor > 0) — for interpreting NES direction
  tf_sets_act <- split(
    tf_regulons$target[tf_regulons$mor > 0],
    tf_regulons$source[tf_regulons$mor > 0]
  )
  tf_sets_act <- tf_sets_act[lengths(tf_sets_act) >= 5]
  
  # Repression targets only (mor < 0)
  tf_sets_rep <- split(
    tf_regulons$target[tf_regulons$mor < 0],
    tf_regulons$source[tf_regulons$mor < 0]
  )
  tf_sets_rep <- tf_sets_rep[lengths(tf_sets_rep) >= 5]
  
} else {
  # Path C: manual fallback using known iron/resistance TFs
  message("  No regulon database available. Using manual curated TF sets.")
  regulon_source <- "manual"
  tf_sets_all <- list(
    NFE2L2  = c("HMOX1","NQO1","GCLC","GCLM","SLC7A11","FTH1","FTL","TXNRD1"),
    STAT3   = c("BCL2","MCL1","CCND1","MYC","VEGFA","FGF2","TWIST1","VIM"),
    HIF1A   = c("LDHA","VEGFA","SLC2A1","HMOX1","BNIP3","BNIP3L","SLC7A11"),
    SP1     = c("TFRC","FTH1","FTL","HMOX1","SLC11A2","NCOA4"),
    MYC     = c("CDK4","CCND2","E2F1","LDHA","PKM","MCM2","MCM4","MCM5"),
    TP53    = c("CDKN1A","BAX","PUMA","MDM2","GADD45A","TIGAR","SCO2"),
    MITF    = c("DCT","TYRP1","TYR","CDK2","BCL2","MLANA","GPR143"),
    ZEB1    = c("VIM","FN1","CDH2","SNAI1","SNAI2","TWIST1","MMP2","MMP9"),
    BACH1   = c("HMOX1","FTH1","FTL","BLVRB","SLC7A11","SQSTM1"),
    TFEB    = c("LAMP1","LAMP2","CTSD","SQSTM1","BECN1","ATG5","MAP1LC3B","NCOA4")
  )
  tf_sets_act <- tf_sets_all
  tf_sets_rep <- list()
}

# ── P3.2  Run fgsea: full TF regulon ─────────────────────────────────────────

message(sprintf("\n  Running fgsea against %s regulons...", regulon_source))

fgsea_tf <- tryCatch(
  fgsea::fgsea(
    pathways    = tf_sets_all,
    stats       = ranked_vec,
    minSize     = 5,
    maxSize     = 500,
    nPermSimple = 10000,
    eps         = 0       # use exact p-values (slower but better for reporting)
  ),
  error = function(e) {
    message("  fgsea failed: ", conditionMessage(e))
    NULL
  }
)

if (!is.null(fgsea_tf) && nrow(fgsea_tf) > 0) {
  
  fgsea_tf <- as.data.frame(fgsea_tf)
  
  # NES interpretation with M = log2(SC/RC):
  #   Positive NES → TF targets enriched among genes HIGH in SC (low in RC)
  #                → TF is LESS active in RC
  #   Negative NES → TF targets enriched among genes LOW in SC (high in RC)
  #                → TF is MORE active in RC (resistance driver candidate)
  fgsea_tf$direction_in_RC <- ifelse(
    fgsea_tf$NES < 0, "More active in RC", "Less active in RC"
  )
  
  # Significance filter
  fgsea_tf_sig <- fgsea_tf[fgsea_tf$padj < 0.05, ]
  if (nrow(fgsea_tf_sig) == 0) {
    fgsea_tf_sig <- fgsea_tf[fgsea_tf$pval < 0.05, ]
    message("  No TFs at padj<0.05; using nominal p<0.05")
  }
  fgsea_tf_sig <- fgsea_tf_sig[order(fgsea_tf_sig$NES), ]  # most negative first
  
  message(sprintf("  Significant TFs: %d", nrow(fgsea_tf_sig)))
  message("  Top RC-activated TFs (NES most negative):")
  print(head(fgsea_tf_sig[, c("pathway","NES","pval","padj","size",
                              "direction_in_RC")], 15))
  
  save_fgsea(fgsea_tf, file.path(OUTPUT_DIR, "TF_fgsea_all.csv"))
  save_fgsea(fgsea_tf_sig, file.path(OUTPUT_DIR, "TF_fgsea_significant.csv"))
  message("  Full and significant TF GSEA results saved.")
  
  # ── P3.3  Priority scoring ──────────────────────────────────────────────────
  # Same logic as before but now NES replaces enrichR adjusted p-value
  
  tf_summary <- fgsea_tf_sig |>
    dplyr::rename(TF_name = pathway) |>
    dplyr::mutate(
      is_DEG       = TF_name %in% deg1$gene,
      is_hub       = TF_name %in% unlist(hub_genes, use.names = FALSE),
      iron_related = TF_name %in% c(IRON_EXTENDED, "NFE2L2","NRF2","BACH1",
                                    "HIF1A","MTF1","TFEB","STAT3","SP1","HMOX1"),
      mapk_akt     = TF_name %in% MAPK_AKT_GENES,
      stem_related = TF_name %in% STEM_GENES,
      # Score: |NES| (activity strength) + pathway weights
      # Negative NES = more active in RC = resistance driver, so weight accordingly
      priority_score = abs(NES) +
        (direction_in_RC == "More active in RC") * 1.0 +
        (is_DEG    * 1.5) +
        (is_hub    * 1.5) +
        (iron_related * 2.5) +
        (mapk_akt  * 1.5) +
        (stem_related * 1.0)
    ) |>
    dplyr::arrange(dplyr::desc(priority_score))
  
  save_fgsea(tf_summary,
            OUTPUT_DIR
            #row.names = FALSE
  )
  
  message("\n  Top 10 prioritized TFs:")
  print(head(tf_summary[, c("TF_name","NES","padj","priority_score",
                            "iron_related","direction_in_RC")], 10))
  
  # ── P3.4  TF GSEA bubble plot ───────────────────────────────────────────────
  
  top25 <- head(tf_summary, 25)
  top25$pathway_cat <- dplyr::case_when(
    top25$iron_related ~ "Iron/Ferroptosis",
    top25$mapk_akt     ~ "MAPK/AKT",
    top25$stem_related ~ "Stemness",
    TRUE               ~ "Other"
  )
  
  p_tf <- ggplot2::ggplot(top25,
                          ggplot2::aes(
                            x     = NES,
                            y     = reorder(TF_name, NES),
                            size  = size,
                            color = pathway_cat,
                            shape = direction_in_RC
                          )) +
    ggplot2::geom_point(alpha = 0.85) +
    ggplot2::scale_color_manual(
      values = c("Iron/Ferroptosis" = "#E64B35", "MAPK/AKT" = "#00A087",
                 "Stemness" = "#7B2D8B", "Other" = "#8491B4"),
      name = "Pathway"
    ) +
    ggplot2::scale_shape_manual(
      values = c("More active in RC" = 17L, "Less active in RC" = 16L),
      name   = "Activity in RC"
    ) +
    ggplot2::scale_size_continuous(range = c(3, 10), name = "Regulon size") +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
                        color = "grey40", linewidth = 0.4) +
    ggplot2::labs(
      title    = "TF Activity in Resistant vs Sensitive Cells",
      subtitle = sprintf("fgsea on %s regulons | ranked by log2FC(SC/RC)\nNegative NES = more active in RC",
                         regulon_source),
      x        = "Normalized Enrichment Score (NES)",
      y        = NULL
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      plot.title  = ggplot2::element_text(face = "bold"),
      legend.position = "right"
    )
  
  ggplot2::ggsave(file.path(OUTPUT_DIR, "TF_fgsea_bubble.pdf"),
                  p_tf, width = 11, height = 8, device = cairo_pdf)
  message("  TF GSEA bubble plot saved.")
  
} else {
  message("  fgsea returned no results.")
}


# ── P3.5  TF activity estimation ─────────────────────────────────────────────
#
# Preferred path:  decoupleR  (Bioc 3.20 standard; no VIPER dep needed)
# Fallback path:   viper + dorothea
#
# Both estimate per-sample TF activity scores (NES) from expression data.
# decoupleR ULM method: fits a univariate linear model of expression ~ regulon.
# We then compute delta = activity(RC) - activity(SC) as the resistance signal.

# Expression input for TF activity: log2(TMM+1), genes × samples
tf_expr_mat <- log2(mrna_norm + 1)

tf_act_df <- NULL   # will hold the result; NULL if neither method succeeds

# ─── Path A: decoupleR ───────────────────────────────────────────────────────

if (requireNamespace("decoupleR", quietly = TRUE)) {
  library(decoupleR)
  message("  TF activity via decoupleR (ULM)...")
  
  # Prefer CollecTRI (broader coverage) over dorothea via decoupleR
  regulon_net <- tryCatch(
    decoupleR::get_collectri(organism = "human", split_complexes = FALSE),
    error = function(e) {
      message(sprintf("  CollecTRI unavailable (%s); trying dorothea...",
                      conditionMessage(e)))
      tryCatch(
        decoupleR::get_dorothea(organism = "human", levels = c("A", "B", "C")),
        error = function(e2) {
          message(sprintf("  dorothea via decoupleR also failed: %s",
                          conditionMessage(e2)))
          NULL
        }
      )
    }
  )
  
  if (!is.null(regulon_net)) {
    message(sprintf("  Regulon loaded: %d TF-target interactions", nrow(regulon_net)))
    
    tf_acts_long <- tryCatch(
      decoupleR::run_ulm(
        mat     = as.matrix(tf_expr_mat),
        net     = regulon_net,
        .source = "source",
        .target = "target",
        .mor    = "mor",
        minsize = 5L
      ),
      error = function(e) {
        message(sprintf("  decoupleR::run_ulm failed: %s", conditionMessage(e)))
        NULL
      }
    )
    
    if (!is.null(tf_acts_long)) {
      # Pivot to wide format: TF (rows) × sample (columns)
      tf_act_wide <- tf_acts_long |>
        dplyr::filter(statistic == "ulm") |>
        dplyr::select(source, condition, score) |>
        tidyr::pivot_wider(names_from = "condition", values_from = "score")
      
      if (!is.data.frame(tf_act_wide) || nrow(tf_act_wide) == 0) {
        message("  decoupleR returned empty activity table.")
      } else {
        tf_act_df <- as.data.frame(tf_act_wide)
        rownames(tf_act_df) <- tf_act_df$source
        tf_act_df$TF <- tf_act_df$source
        tf_act_df$source <- NULL
        message(sprintf("  decoupleR: %d TF activity profiles computed.", nrow(tf_act_df)))
      }
    }
  }
}

# ─── Path B: viper + dorothea fallback ───────────────────────────────────────

if (is.null(tf_act_df) &&
    requireNamespace("viper",    quietly = TRUE) &&
    requireNamespace("dorothea", quietly = TRUE)) {
  
  message("  TF activity via viper + dorothea (fallback)...")
  library(viper); library(dorothea)
  
  # dorothea API varies by version; try multiple load approaches
  regulon_df <- tryCatch({
    # Modern dorothea: direct data access
    dorothea::dorothea_hs |>
      dplyr::filter(confidence %in% c("A", "B"))
  }, error = function(e) {
    tryCatch({
      # Older dorothea: data() loading
      tmp_env <- new.env()
      data("dorothea_hs", package = "dorothea", envir = tmp_env)
      get("dorothea_hs", envir = tmp_env) |>
        dplyr::filter(confidence %in% c("A", "B"))
    }, error = function(e2) {
      message(sprintf("  Could not load dorothea regulon: %s", conditionMessage(e2)))
      NULL
    })
  })
  
  if (!is.null(regulon_df) && nrow(regulon_df) > 0) {
    # Restrict expression matrix to regulated targets
    viper_mat <- tf_expr_mat[rownames(tf_expr_mat) %in% unique(regulon_df$target), ]
    
    viper_res <- tryCatch(
      viper::viper(
        eset        = viper_mat,
        regulon     = viper::df2regulon(regulon_df),
        minsize     = 4L,
        eset.filter = FALSE,
        verbose     = FALSE
      ),
      error = function(e) {
        message(sprintf("  viper::viper failed: %s", conditionMessage(e)))
        NULL
      }
    )
    
    if (!is.null(viper_res)) {
      sample_cols <- intersect(
        colnames(viper_res),
        c("A375_RC","A375_R10","A375_SC","A375_S10")
      )
      tf_act_df <- as.data.frame(viper_res[, sample_cols, drop = FALSE]) |>
        tibble::rownames_to_column("TF")
      message(sprintf("  VIPER: %d TF activity profiles computed.", nrow(tf_act_df)))
    }
  }
}

# ─── Post-process activity results ────────────────────────────────────────────

if (!is.null(tf_act_df) && nrow(tf_act_df) > 0) {
  
  # Ensure correct column names for the four samples
  sample_cols_present <- intersect(
    c("A375_RC","A375_R10","A375_SC","A375_S10"),
    colnames(tf_act_df)
  )
  message(sprintf("  Activity samples available: %s",
                  paste(sample_cols_present, collapse = ", ")))
  
  # Resistance signal: delta = activity(RC) - activity(SC)
  if (all(c("A375_RC","A375_SC") %in% sample_cols_present)) {
    tf_act_df$delta_RC_vs_SC <- tf_act_df$A375_RC - tf_act_df$A375_SC
  } else if (length(sample_cols_present) >= 2) {
    # If sample names don't match exactly, use first two columns
    tf_act_df$delta_RC_vs_SC <- tf_act_df[[sample_cols_present[1]]] -
      tf_act_df[[sample_cols_present[2]]]
    message("  Note: delta_RC_vs_SC computed from first two available samples.")
  }
  
  tf_act_df <- tf_act_df |>
    dplyr::mutate(
      iron_related = TF %in% c(IRON_EXTENDED, "NFE2L2","NRF2","BACH1",
                               "HIF1A","MTF1","TFEB","STAT3","SP1"),
      mapk_akt     = TF %in% MAPK_AKT_GENES,
      stem_related = TF %in% STEM_GENES
    )
  
  if ("delta_RC_vs_SC" %in% colnames(tf_act_df)) {
    tf_act_df <- dplyr::arrange(tf_act_df, desc(abs(delta_RC_vs_SC)))
  }
  
  write.csv(tf_act_df,
            file.path(OUTPUT_DIR, "TF_activity_scores.csv"),
            row.names = FALSE)
  message(sprintf("  TF activity scores saved: %d TFs.", nrow(tf_act_df)))
  
  # ── P3.6  TF activity heatmap ──────────────────────────────────────────────
  
  if ("delta_RC_vs_SC" %in% colnames(tf_act_df) &&
      length(sample_cols_present) >= 2) {
    
    # Top 40 TFs by |delta| across resistance contrast
    q90 <- quantile(abs(tf_act_df$delta_RC_vs_SC), 0.90, na.rm = TRUE)
    top_act <- tf_act_df[abs(tf_act_df$delta_RC_vs_SC) >= q90, ]
    top_act <- head(top_act[order(abs(top_act$delta_RC_vs_SC), decreasing = TRUE), ], 40)
    
    if (nrow(top_act) >= 5) {
      act_mat <- as.matrix(top_act[, sample_cols_present])
      rownames(act_mat) <- top_act$TF
      
      # Row annotations
      act_row_ann <- data.frame(
        Pathway   = dplyr::case_when(
          top_act$iron_related ~ "Iron/Ferroptosis",
          top_act$mapk_akt     ~ "MAPK/AKT",
          top_act$stem_related ~ "Stemness",
          TRUE                 ~ "Other"
        ),
        Direction = ifelse(top_act$delta_RC_vs_SC > 0,
                           "More active RC", "Less active RC"),
        row.names = top_act$TF
      )
      
      # Column annotations — map sample names to condition/treatment
      cond_lookup <- c(A375_RC="Resistant",A375_R10="Resistant",
                       A375_SC="Sensitive",A375_S10="Sensitive")
      trt_lookup  <- c(A375_RC="Control",  A375_R10="Drug",
                       A375_SC="Control",  A375_S10="Drug")
      act_col_ann <- data.frame(
        Condition = cond_lookup[sample_cols_present],
        Treatment = trt_lookup[sample_cols_present],
        row.names = sample_cols_present
      )
      
      act_annot_col <- list(
        Pathway   = c("Iron/Ferroptosis"="#E64B35","MAPK/AKT"="#00A087",
                      "Stemness"="#7B2D8B","Other"="#8491B4"),
        Direction = c("More active RC"="#B2182B","Less active RC"="#2166AC"),
        Condition = c(Resistant="#CC0000",Sensitive="#0066CC"),
        Treatment = c(Control="#999999",Drug="#FF9900")
      )
      
      pheatmap::pheatmap(
        act_mat,
        color             = colorRampPalette(c("#2166AC","white","#B2182B"))(100),
        annotation_row    = act_row_ann,
        annotation_col    = act_col_ann,
        annotation_colors = act_annot_col,
        cluster_cols      = FALSE,
        cluster_rows      = TRUE,
        scale             = "row",
        fontsize_row      = 7,
        main              = "TF Activity Scores (NES | z-scaled per TF)\nTop differentially active TFs: RC vs SC",
        filename          = file.path(OUTPUT_DIR, "TF_activity_heatmap.pdf"),
        width = 7, height = 10
      )
      message("  TF activity heatmap saved.")
    }
  }
  
  # ── P3.7  enrichR + activity convergence analysis ─────────────────────────
  # TFs appearing in BOTH enrichR (overrepresentation) AND activity scoring
  # (differential NES) are the strongest mechanistic candidates.
  
  if (!is.null(tf_summary) && nrow(tf_summary) > 0 &&
      "delta_RC_vs_SC" %in% colnames(tf_act_df)) {
    
    # Top-10% most differentially active TFs
    q90_act <- quantile(abs(tf_act_df$delta_RC_vs_SC), 0.90, na.rm = TRUE)
    top_act_tfs <- tf_act_df$TF[abs(tf_act_df$delta_RC_vs_SC) >= q90_act]
    
    enrichr_tfs  <- tf_summary$TF_name
    converged    <- intersect(enrichr_tfs, top_act_tfs)
    
    message(sprintf("\n  === TF CONVERGENCE SUMMARY ==="))
    message(sprintf("  enrichR significant TFs:           %d", length(enrichr_tfs)))
    message(sprintf("  Top-10%% activity TFs (|delta NES|):%d", length(top_act_tfs)))
    message(sprintf("  CONVERGED (both layers):           %d", length(converged)))
    
    if (length(converged) > 0) {
      message("  Converged TFs: ", paste(converged, collapse = ", "))
      
      # Build convergence table: merge enrichR priority + activity delta
      conv_df <- dplyr::left_join(
        tf_summary[tf_summary$TF_name %in% converged, ],
        tf_act_df[tf_act_df$TF %in% converged,
                  c("TF","delta_RC_vs_SC")],
        by = c("TF_name" = "TF")
      ) |>
        dplyr::mutate(
          activity_direction = dplyr::case_when(
            !is.na(delta_RC_vs_SC) & delta_RC_vs_SC > 0 ~ "More active in RC",
            !is.na(delta_RC_vs_SC) & delta_RC_vs_SC < 0 ~ "Less active in RC",
            TRUE                                         ~ "Not computed"
          )
        ) |>
        dplyr::arrange(desc(priority_score))
      
      write.csv(conv_df,
                file.path(OUTPUT_DIR, "TF_converged_final_candidates.csv"),
                row.names = FALSE)
      message("  Converged TF table saved.")
      
      # Quick visual: convergence scatter
      if (length(converged) >= 3 && !is.null(tf_summary)) {
        plot_conv <- dplyr::inner_join(
          tf_summary[, c("TF_name","priority_score","iron_related","stem_related")],
          tf_act_df[, c("TF","delta_RC_vs_SC")],
          by = c("TF_name" = "TF")
        ) |>
          dplyr::filter(!is.na(delta_RC_vs_SC)) |>
          dplyr::mutate(
            is_converged = TF_name %in% converged,
            label_gene   = ifelse(is_converged, TF_name, NA_character_)
          )
        
        p_conv <- ggplot2::ggplot(plot_conv,
                                  ggplot2::aes(
                                    x     = delta_RC_vs_SC,
                                    y     = priority_score,
                                    color = is_converged,
                                    label = label_gene
                                  )) +
          ggplot2::geom_point(alpha = 0.6, size = 2) +
          ggrepel::geom_text_repel(na.rm = TRUE, size = 2.8, fontface = "italic",
                                   max.overlaps = 20) +
          ggplot2::scale_color_manual(
            values = c(`TRUE` = "#E64B35", `FALSE` = "grey60"),
            labels = c("Other", "Converged"),
            name   = ""
          ) +
          ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
                              color = "grey40", linewidth = 0.4) +
          ggplot2::labs(
            title    = "TF Activity vs enrichR Priority Score",
            subtitle = "Red = TFs significant in BOTH layers (top candidates)",
            x        = "TF activity delta (RC - SC | NES)",
            y        = "enrichR priority score"
          ) +
          ggplot2::theme_bw(base_size = 11) +
          ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))
        
        ggplot2::ggsave(
          file.path(OUTPUT_DIR, "TF_convergence_scatter.pdf"),
          p_conv, width = 8, height = 6)
        message("  Convergence scatter plot saved.")
      }
    } else {
      message("  No TFs converged across both layers (small n is expected cause).")
    }
  }
  
} else {
  message("  TF activity scoring unavailable.")
  message("  Install decoupleR (primary): BiocManager::install('decoupleR')")
  message("  Or install viper + dorothea (fallback): BiocManager::install(c('viper','dorothea'))")
}
