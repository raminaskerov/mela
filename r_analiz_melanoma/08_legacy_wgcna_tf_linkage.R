# =============================================================================
# 08_legacy_wgcna_tf_linkage.R
# INTEGRATIVE ANALYSIS: ENCORAFENIB RESISTANCE IN MALIGNANT MELANOMA
#
# Input:   checkpoint_06.RData + checkpoint_07.RData
# Output:  legacy TF/WGCNA linkage tables, checkpoint_08.RData
#
# This is an optional compatibility stage for the old 2.13 logic.
# It is intentionally separate from the core WGCNA stage.
# =============================================================================


# ── Configuration ─────────────────────────────────────────────────────────────

DATA_DIR   <- "/home/ramin/mela/data/ilk"  # <-- UPDATE THIS PATH
OUTPUT_DIR <- file.path(DATA_DIR, "newoutput")
LEGACY_DIR <- file.path(OUTPUT_DIR, "legacy_linkage")
dir.create(LEGACY_DIR, showWarnings = FALSE, recursive = TRUE)


# ── Packages ──────────────────────────────────────────────────────────────────

for (pkg in c("dplyr", "ggplot2", "ggrepel")) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
})


# ── Load checkpoints ──────────────────────────────────────────────────────────

load(file.path(DATA_DIR, "checkpoint_06.RData"))
load(file.path(DATA_DIR, "checkpoint_07.RData"))

if (is.null(fgsea_tf_sig) || nrow(fgsea_tf_sig) == 0) {
  stop("checkpoint_06.RData does not contain TF fgsea significant results.")
}
if (is.null(module_assign_df) || nrow(module_assign_df) == 0) {
  stop("checkpoint_07.RData does not contain WGCNA module assignments.")
}


# ── Legacy target gene list ───────────────────────────────────────────────────

target_genes <- c(
  "TFRC", "FTH1", "TCIRG1", "ATP6V1F", "NDUFA9", "HLA-A", "IKBKG", "ABL1",
  "SNAI1", "TWIST1", "HFE", "NDUFA2",
  "CHEK1", "MRE11", "EYA1", "HIPK2", "APAF1",
  "MDM4", "CDC6", "LIN52", "FAS", "ATP6V1C1", "ATP5F1E", "GSTM4", "ETV2"
)


# =============================================================================
# Legacy TF leading-edge overlap
# =============================================================================

message("\n=== LEGACY TF/WGCNA LINKAGE ===")

tf_leading_rows <- list()

for (i in seq_len(nrow(fgsea_tf_sig))) {
  tf_name <- fgsea_tf_sig$pathway[i]
  tf_nes  <- fgsea_tf_sig$NES[i]
  tf_dir  <- if ("direction_in_RC" %in% colnames(fgsea_tf_sig)) {
    fgsea_tf_sig$direction_in_RC[i]
  } else {
    ifelse(tf_nes > 0, "More active in RC", "Less active in RC")
  }

  leading <- fgsea_tf_sig$leadingEdge[[i]]
  if (is.character(leading) && length(leading) == 1L) {
    leading <- unlist(strsplit(leading, ";", fixed = TRUE), use.names = FALSE)
  }
  leading <- unique(trimws(as.character(leading)))
  leading <- leading[nzchar(leading)]
  if (length(leading) == 0) {
    next
  }

  hits <- intersect(leading, target_genes)
  if (length(hits) == 0) {
    next
  }

  for (gene in hits) {
    tf_leading_rows[[length(tf_leading_rows) + 1]] <- data.frame(
      TF = tf_name,
      TF_NES = tf_nes,
      TF_direction = tf_dir,
      target_gene = gene,
      target_in_legacy_list = TRUE,
      leading_edge_contains_target = TRUE,
      stringsAsFactors = FALSE
    )
  }
}

if (length(tf_leading_rows) == 0) {
  stop("No TF leading-edge hits matched the legacy target_genes list.")
}

tf_target_hits <- dplyr::bind_rows(tf_leading_rows)


# =============================================================================
# WGCNA annotation for the legacy target list
# =============================================================================

target_context <- data.frame(
  gene = target_genes,
  stringsAsFactors = FALSE
) |>
  dplyr::left_join(
    module_assign_df |>
      dplyr::select(gene, module, in_resistance_module, iron_related, stem_related, mapk_related),
    by = "gene"
  )

if (exists("hub_all_df") && nrow(hub_all_df) > 0) {
  hub_lookup <- hub_all_df |>
    dplyr::group_by(gene) |>
    dplyr::summarise(
      max_kME = max(kME, na.rm = TRUE),
      best_module = module[which.max(kME)],
      .groups = "drop"
    )
  target_context <- target_context |>
    dplyr::left_join(hub_lookup, by = "gene")
}


# =============================================================================
# Combined legacy table
# =============================================================================

legacy_linkage_df <- tf_target_hits |>
  dplyr::left_join(
    target_context,
    by = c("target_gene" = "gene")
  ) |>
  dplyr::mutate(
    is_wgcna_hub = !is.na(max_kME),
    legacy_priority_score = abs(TF_NES) +
      (target_in_legacy_list * 1.0) +
      (is_wgcna_hub * 1.5) +
      (dplyr::coalesce(in_resistance_module, FALSE) * 1.0) +
      (dplyr::coalesce(iron_related, FALSE) * 1.5)
  ) |>
  dplyr::arrange(dplyr::desc(legacy_priority_score), dplyr::desc(abs(TF_NES)))

write.csv(
  legacy_linkage_df,
  file.path(LEGACY_DIR, "TF_target_legacy_linkage.csv"),
  row.names = FALSE
)

write.csv(
  target_context,
  file.path(LEGACY_DIR, "legacy_target_gene_context.csv"),
  row.names = FALSE
)

tf_summary <- legacy_linkage_df |>
  dplyr::group_by(TF) |>
  dplyr::summarise(
    n_hits = n(),
    n_wgcna_hub_hits = sum(is_wgcna_hub, na.rm = TRUE),
    n_resistance_module_hits = sum(in_resistance_module, na.rm = TRUE),
    top_hit_gene = target_gene[which.max(legacy_priority_score)],
    max_legacy_score = max(legacy_priority_score, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::arrange(dplyr::desc(max_legacy_score))

write.csv(
  tf_summary,
  file.path(LEGACY_DIR, "TF_target_legacy_summary.csv"),
  row.names = FALSE
)

if (nrow(tf_summary) >= 3) {
  p_legacy <- ggplot(tf_summary, aes(
    x = reorder(TF, max_legacy_score),
    y = max_legacy_score,
    size = n_hits,
    color = n_wgcna_hub_hits
  )) +
    geom_point(alpha = 0.85) +
    coord_flip() +
    scale_color_gradient(low = "#2166AC", high = "#B2182B", name = "Hub hits") +
    scale_size_continuous(range = c(2.5, 8), name = "Target hits") +
    labs(
      title = "Legacy TF/WGCNA linkage",
      subtitle = "TF leading-edge hits against the old target gene list",
      x = NULL,
      y = "Legacy priority score"
    ) +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(face = "bold"))

  ggsave(file.path(LEGACY_DIR, "TF_target_legacy_linkage.pdf"), p_legacy,
         width = 10, height = max(6, nrow(tf_summary) * 0.25 + 2))
}


# =============================================================================
# CHECKPOINT SAVE
# =============================================================================

save(
  target_genes,
  tf_target_hits, target_context, legacy_linkage_df, tf_summary,
  file = file.path(DATA_DIR, "checkpoint_08.RData")
)

message("\ncheckpoint_08.RData saved.")
