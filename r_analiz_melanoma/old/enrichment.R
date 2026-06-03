# =============================================================================
# ENRICHMENT ANALYSIS: ORA + GSEA + STRINGdb PPI
# Requires: comp1.1–comp4.1 in session (from noiseq_adj.R)
# =============================================================================

# ── Packages ──────────────────────────────────────────────────────────────────

bioc_needed <- c("clusterProfiler", "org.Hs.eg.db", "enrichplot", "STRINGdb")
for (p in bioc_needed) {
  if (!requireNamespace(p, quietly = TRUE)) BiocManager::install(p)
}

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(STRINGdb)
  library(ggplot2)
  library(dplyr)
})

# ── Comparison registry ───────────────────────────────────────────────────────

comparisons <- list(
  list(obj = comp1.2, label = "SC_vs_RC"),
  list(obj = comp2.2, label = "SC_vs_S10"),
  list(obj = comp3.2, label = "RC_vs_R10"),
  list(obj = comp4.2, label = "S10_vs_R10")
)

# ── Helper: symbol → ENTREZID ────────────────────────────────────────────────

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

# =============================================================================
# MAIN LOOP: ORA + GSEA per comparison
# =============================================================================

for (comp in comparisons) {
  
  lbl  <- comp$label
  obj  <- comp$obj
  message(sprintf("\n=== %s ===", lbl))
  
  # ── Gene lists ─────────────────────────────────────────────────────────────
  
  deg_genes <- rownames(obj$degs)         # filtered DEGs (priority score)
  message(sprintf("  DEGs: %d", length(deg_genes)))
  
  if (length(deg_genes) < 5) {
    message("  Too few DEGs — skipping ORA/STRINGdb for this comparison.")
    next
  }
  
  # Symbol → ENTREZID for DEGs
  deg_mapped <- tryCatch(
    to_entrez(deg_genes),
    error = function(e) { message("  bitr failed: ", e$message); NULL }
  )
  if (is.null(deg_mapped) || nrow(deg_mapped) == 0) next
  
  # Full ranked list for GSEA (all expressed genes)
  full_mapped <- tryCatch(
    to_entrez(rownames(obj$full)),
    error = function(e) NULL
  )
  
  # ── ORA: GO Biological Process ─────────────────────────────────────────────
  
  message("  Running ORA GO-BP...")
  ora_go <- tryCatch(
    clusterProfiler::enrichGO(
      gene          = deg_mapped$ENTREZID,
      OrgDb         = org.Hs.eg.db,
      ont           = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      qvalueCutoff  = 0.20,
      readable      = TRUE
    ),
    error = function(e) { message("  ORA GO failed: ", e$message); NULL }
  )
  
  if (!is.null(ora_go) && nrow(as.data.frame(ora_go)) > 0) {
    write.csv(as.data.frame(ora_go),
              file.path(OUTPUT_DIR, sprintf("ORA_GO_BP_%s.csv", lbl)),
              row.names = FALSE)
    p_ora <- enrichplot::dotplot(ora_go, showCategory = 20) +
      ggtitle(sprintf("ORA GO-BP | %s", lbl))
    ggsave(file.path(OUTPUT_DIR, sprintf("ORA_GO_BP_%s.pdf", lbl)),
           p_ora, width = 10, height = 9)
    message(sprintf("  ORA GO-BP: %d terms", nrow(as.data.frame(ora_go))))
  } else {
    message("  ORA GO-BP: no significant terms.")
  }
  
  # ── ORA: KEGG ──────────────────────────────────────────────────────────────
  
  message("  Running ORA KEGG...")
  ora_kegg <- tryCatch(
    clusterProfiler::enrichKEGG(
      gene          = deg_mapped$ENTREZID,
      organism      = "hsa",
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      qvalueCutoff  = 0.20
    ),
    error = function(e) { message("  ORA KEGG failed: ", e$message); NULL }
  )
  
  if (!is.null(ora_kegg) && nrow(as.data.frame(ora_kegg)) > 0) {
    write.csv(as.data.frame(ora_kegg),
              file.path(OUTPUT_DIR, sprintf("ORA_KEGG_%s.csv", lbl)),
              row.names = FALSE)
    p_kegg <- enrichplot::dotplot(ora_kegg, showCategory = 20) +
      ggtitle(sprintf("ORA KEGG | %s", lbl))
    ggsave(file.path(OUTPUT_DIR, sprintf("ORA_KEGG_%s.pdf", lbl)),
           p_kegg, width = 10, height = 8)
    message(sprintf("  ORA KEGG: %d pathways", nrow(as.data.frame(ora_kegg))))
  } else {
    message("  ORA KEGG: no significant pathways.")
  }
  
 
}
# ORA input — filtered DEG gene lists
write.table(rownames(comp1.2$degs), 
            file.path(OUTPUT_DIR, "webgestalt_ORA_SC_vs_RC_filtered.txt"))
write_lines(rownames(comp2.2$degs),
            file.path(OUTPUT_DIR, "webgestalt_ORA_SC_vs_S10_filtered.txt"))
write_lines(rownames(comp3.2$degs),
            file.path(OUTPUT_DIR, "webgestalt_ORA_RC_vs_R10_filtered.txt"))
write_lines(rownames(comp4.2$degs),
            file.path(OUTPUT_DIR, "webgestalt_ORA_S10_vs_R10_filtered.txt"))

comparisons_gsea <- list(
  list(comp=comp1, label = "SC_vs_RC",  A = "A375_SC", B = "A375_RC"),
  list(comp=comp2, label = "SC_vs_S10", A = "A375_SC", B = "A375_S10"),
  list(comp=comp3, label = "RC_vs_R10", A = "A375_RC", B = "A375_R10"),
  list(comp=comp4, label = "S10_vs_R10",A = "A375_S10",B = "A375_R10")
)

write.table((meqale$"Gene"),
            file      = file.path(OUTPUT_DIR, "meqale.txt"),
            row.names = FALSE,
            col.names = FALSE,
            quote     = FALSE)
write.table(rownames(comp2.2$degs),
            file      = file.path(OUTPUT_DIR, "webg_ora2.txt"),
            row.names = FALSE,
            col.names = FALSE,
            quote     = FALSE)

for (c in comparisons_gsea) {
  ranked <- data.frame(
    gene  = rownames(c$comp$degs),
    score = c$comp$degs$M        # log2FC column from NOISeq full results
  ) |> dplyr::arrange(desc(score))
  
  write.table(ranked,
              file      = file.path(OUTPUT_DIR, 
                                    paste0("webgestalt_GSEA_*", c$label, ".rnk")),
              sep       = "\t",
              row.names = FALSE,
              col.names = FALSE,
              quote     = FALSE)
}
ranked <- data.frame(
  gene  = rownames(comp1ş2$degs),
  score = comp4$degs$M        # log2FC column from NOISeq full results
) |> dplyr::arrange(desc(score))

write.table(ranked,
            file      = file.path(OUTPUT_DIR, 
                                  paste0("webgestalt_GSEA.2_","S10_vs_R10" , ".rnk")),
            sep       = "\t",
            row.names = FALSE,
            col.names = FALSE,
            quote     = FALSE)
message("\n=== Enrichment analysis complete ===")
message(sprintf("All outputs in: %s/", OUTPUT_DIR))


for (comp in comparisons){
  lbl  <- comp$label
  obj  <- comp$obj
  message(sprintf("\n=== %s ===", lbl))
}
# ── Gene lists ─────────────────────────────────────────────────────────────

  deg_genes <- rownames(obj$degs)         # filtered DEGs (priority score)
  message(sprintf("  DEGs: %d", length(deg_genes)))

if (length(deg_genes) < 5) {
  message("  Too few DEGs — skipping ORA/STRINGdb for this comparison.")
  next
}

# Symbol → ENTREZID for DEGs
deg_mapped <- tryCatch(
  to_entrez(deg_genes),
  error = function(e) { message("  bitr failed: ", e$message); NULL }
)
if (is.null(deg_mapped) || nrow(deg_mapped) == 0) next

# Full ranked list for GSEA (all expressed genes)
full_mapped <- tryCatch(
  to_entrez(rownames(obj$full)),
  error = function(e) NULL
)

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
  message(sprintf("  GSEA KEGG: %d pathways", nrow(as.data.frame(gsea_kegg))))
} else {
  message("  GSEA KEGG: no significant pathways.")
}
# ── GSEA: GO Biological Process ────────────────────────────────────────────

if (!is.null(full_mapped)) {
  message("  Running GSEA GO-BP...")
  ranked_vec <- build_ranked_list(obj$full, full_mapped)
  
  gsea_go <- tryCatch(
    clusterProfiler::gseGO(
      geneList      = ranked_vec,
      OrgDb         = org.Hs.eg.db,
      ont           = "BP",
      minGSSize     = 10,
      maxGSSize     = 500,
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      verbose       = FALSE
    ),
    error = function(e) { message("  GSEA GO failed: ", e$message); NULL }
  )
  
  if (!is.null(gsea_go) && nrow(as.data.frame(gsea_go)) > 0) {
    write.csv(as.data.frame(gsea_go),
              file.path(OUTPUT_DIR, sprintf("GSEA_GO_BP_%s.csv", lbl)),
              row.names = FALSE)
    p_gsea <- enrichplot::dotplot(gsea_go, showCategory = 20,
                                  split = ".sign") +
      facet_grid(. ~ .sign) +
      ggtitle(sprintf("GSEA GO-BP | %s", lbl))
    ggsave(file.path(OUTPUT_DIR, sprintf("GSEA_GO_BP_%s.pdf", lbl)),
           p_gsea, width = 12, height = 9)
    message(sprintf("  GSEA GO-BP: %d terms", nrow(as.data.frame(gsea_go))))
  } else {
    message("  GSEA GO-BP: no significant terms.")
  }
  
  # ── GSEA: KEGG ─────────────────────────────────────────────────────────
  
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
    message(sprintf("  GSEA KEGG: %d pathways", nrow(as.data.frame(gsea_kegg))))
  } else {
    message("  GSEA KEGG: no significant pathways.")
  }
}

# ── STRINGdb PPI network ────────────────────────────────────────────────────

message("  Building STRINGdb PPI network...")
string_db <- STRINGdb$new(
  version      = "11.5",
  species      = 9606,          # human
  score_threshold = 400,        # medium confidence
  network_type = "functional",
  input_directory = OUTPUT_DIR
)

# Map DEG gene symbols
deg_df <- data.frame(gene = deg_genes, stringsAsFactors = FALSE)
deg_mapped_string <- tryCatch(
  string_db$map(deg_df, "gene", removeUnmappedRows = TRUE),
  error = function(e) { message("  STRINGdb mapping failed: ", e$message); NULL }
)

if (!is.null(deg_mapped_string) && nrow(deg_mapped_string) > 2) {
  # Save network image
  pdf(file.path(OUTPUT_DIR, sprintf("STRINGdb_PPI_%s.pdf", lbl)),
      width = 10, height = 10)
  string_db$plot_network(deg_mapped_string$STRING_id)
  dev.off()
  
  # Save interaction table
  interactions <- tryCatch(
    string_db$get_interactions(deg_mapped_string$STRING_id),
    error = function(e) NULL
  )
  if (!is.null(interactions)) {
    write.csv(interactions,
              file.path(OUTPUT_DIR, sprintf("STRINGdb_interactions_%s.csv", lbl)),
              row.names = FALSE)
    message(sprintf("  STRINGdb: %d nodes, %d interactions",
                    nrow(deg_mapped_string), nrow(interactions)))
  }
} else {
  message("  STRINGdb: too few mapped genes for network.")
}

