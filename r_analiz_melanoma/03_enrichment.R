# =============================================================================
# 03_enrichment.R
# INTEGRATIVE ANALYSIS: ENCORAFENIB RESISTANCE IN MALIGNANT MELANOMA
# GSE283251 | Colakoglu Bergel et al., Scientific Reports 2025
#
# Input:   checkpoint_02.RData
# Output:  pathway enrichment tables, .rnk files, checkpoint_03.RData
#
# This stage is pathway-only. It does not depend on TF or WGCNA outputs.
# =============================================================================


# ── Configuration ─────────────────────────────────────────────────────────────

DATA_DIR   <- "/home/ramin/mela/data/ilk"  # <-- UPDATE THIS PATH
OUTPUT_DIR <- file.path(DATA_DIR, "newoutput")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)


# ── Packages ──────────────────────────────────────────────────────────────────

cran_pkgs <- c("dplyr", "tidyr", "ggplot2", "ggrepel", "pheatmap")
bioc_pkgs  <- c("clusterProfiler", "ReactomePA", "msigdbr", "enrichplot",
                "org.Hs.eg.db")

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
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(clusterProfiler)
  library(ReactomePA)
  library(msigdbr)
  library(enrichplot)
  library(org.Hs.eg.db)
})


# ── Load checkpoint ───────────────────────────────────────────────────────────

load(file.path(DATA_DIR, "checkpoint_02.RData"))


# ── Helpers ───────────────────────────────────────────────────────────────────

to_entrez <- function(symbols) {
  symbols <- unique(symbols)
  symbols <- symbols[!is.na(symbols) & nzchar(symbols)]
  if (length(symbols) == 0) {
    return(data.frame(SYMBOL = character(), ENTREZID = character()))
  }

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

build_ranked_list <- function(full_df, mapped_df) {
  ranked_df <- as.data.frame(full_df)
  ranked_df$gene <- rownames(ranked_df)

  ranked_df <- ranked_df |>
    dplyr::filter(!is.na(M)) |>
    dplyr::select(gene, M) |>
    dplyr::inner_join(mapped_df, by = c("gene" = "SYMBOL")) |>
    dplyr::arrange(dplyr::desc(M)) |>
    dplyr::distinct(ENTREZID, .keep_all = TRUE)

  ranked_vec <- ranked_df$M
  names(ranked_vec) <- ranked_df$ENTREZID
  ranked_vec
}

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

run_pathway_gsea <- function(full_df, label) {
  symbols <- rownames(full_df)
  mapped  <- to_entrez(symbols)
  if (nrow(mapped) == 0) {
    message("  No ENTREZ mappings available; skipping ", label)
    return(list(
      label = label,
      ranked_vec = numeric(),
      kegg = NULL,
      reactome = NULL,
      kegg_df = NULL,
      reactome_df = NULL
    ))
  }

  ranked_vec <- build_ranked_list(full_df, mapped)
  message(sprintf("  %s ranked list: %d genes | range [%.2f, %.2f]",
                  label, length(ranked_vec), min(ranked_vec), max(ranked_vec)))

  kegg_res <- tryCatch(
    clusterProfiler::gseKEGG(
      geneList      = ranked_vec,
      organism      = "hsa",
      minGSSize     = 10,
      maxGSSize     = 500,
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      verbose       = FALSE
    ),
    error = function(e) {
      message("  GSEA KEGG failed for ", label, ": ", conditionMessage(e))
      NULL
    }
  )

  reactome_res <- tryCatch(
    ReactomePA::gsePathway(
      geneList      = ranked_vec,
      organism      = "human",
      minGSSize     = 10,
      maxGSSize     = 500,
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      verbose       = FALSE
    ),
    error = function(e) {
      message("  GSEA Reactome failed for ", label, ": ", conditionMessage(e))
      NULL
    }
  )

  kegg_df <- if (!is.null(kegg_res) && nrow(as.data.frame(kegg_res)) > 0) {
    as.data.frame(kegg_res)
  } else {
    NULL
  }
  reactome_df <- if (!is.null(reactome_res) && nrow(as.data.frame(reactome_res)) > 0) {
    as.data.frame(reactome_res)
  } else {
    NULL
  }

  list(
    label = label,
    ranked_vec = ranked_vec,
    kegg = kegg_res,
    reactome = reactome_res,
    kegg_df = kegg_df,
    reactome_df = reactome_df
  )
}


# =============================================================================
# PATHWAY ENRICHMENT
# =============================================================================

message("\n=== PATHWAY ENRICHMENT ===")

comparison_list <- list(
  list(obj = comp1, label = "SC_vs_RC"),
  list(obj = comp2, label = "SC_vs_S10"),
  list(obj = comp3, label = "RC_vs_R10"),
  list(obj = comp4, label = "S10_vs_R10")
)

pathway_results <- list()

for (comp in comparison_list) {
  lbl <- comp$label
  message(sprintf("\n--- %s ---", lbl))

  full_df <- as.data.frame(comp$obj$full)
  full_df$gene <- rownames(full_df)

  pathway_results[[lbl]] <- run_pathway_gsea(full_df, lbl)

  if (length(pathway_results[[lbl]]$ranked_vec) == 0) {
    next
  }

  write.table(
    data.frame(
      gene = names(pathway_results[[lbl]]$ranked_vec),
      score = unname(pathway_results[[lbl]]$ranked_vec)
    ) |>
      dplyr::arrange(dplyr::desc(score)),
    file      = file.path(OUTPUT_DIR, paste0("webgestalt_GSEA_", lbl, ".rnk")),
    sep       = "\t",
    row.names = FALSE,
    col.names = FALSE,
    quote     = FALSE
  )

  if (!is.null(pathway_results[[lbl]]$kegg_df) &&
      nrow(pathway_results[[lbl]]$kegg_df) > 0) {
    save_fgsea_table(
      pathway_results[[lbl]]$kegg_df,
      file.path(OUTPUT_DIR, paste0("GSEA_KEGG_", lbl, ".csv"))
    )

    p_kegg <- tryCatch(
      enrichplot::dotplot(pathway_results[[lbl]]$kegg, showCategory = 20, split = ".sign") +
        facet_grid(. ~ .sign) +
        ggtitle(sprintf("GSEA KEGG | %s", lbl)),
      error = function(e) NULL
    )
    if (!is.null(p_kegg)) {
      ggsave(
        file.path(OUTPUT_DIR, paste0("GSEA_KEGG_", lbl, ".pdf")),
        p_kegg,
        width = 12,
        height = 8
      )
    }
    message(sprintf("  KEGG pathways: %d", nrow(pathway_results[[lbl]]$kegg_df)))
  } else {
    message("  KEGG: no significant pathways.")
  }

  if (!is.null(pathway_results[[lbl]]$reactome_df) &&
      nrow(pathway_results[[lbl]]$reactome_df) > 0) {
    save_fgsea_table(
      pathway_results[[lbl]]$reactome_df,
      file.path(OUTPUT_DIR, paste0("GSEA_Reactome_", lbl, ".csv"))
    )

    p_reactome <- tryCatch(
      enrichplot::dotplot(pathway_results[[lbl]]$reactome, showCategory = 20, split = ".sign") +
        facet_grid(. ~ .sign) +
        ggtitle(sprintf("GSEA Reactome | %s", lbl)),
      error = function(e) NULL
    )
    if (!is.null(p_reactome)) {
      ggsave(
        file.path(OUTPUT_DIR, paste0("GSEA_Reactome_", lbl, ".pdf")),
        p_reactome,
        width = 12,
        height = 9
      )
    }
    message(sprintf("  Reactome pathways: %d", nrow(pathway_results[[lbl]]$reactome_df)))
  } else {
    message("  Reactome: no significant pathways.")
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
  pathway_results,
  file = file.path(DATA_DIR, "checkpoint_03.RData")
)

message("\ncheckpoint_03.RData saved.")
