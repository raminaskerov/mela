tmm_normalize <- function(count_mat) {
  lib_sizes <- colSums(count_mat)
  ref_lib   <- exp(mean(log(lib_sizes[lib_sizes > 0])))
  scale_fac <- lib_sizes / ref_lib
  sweep(count_mat, 2, scale_fac, FUN = "/")
}

mrna_norm <- tmm_normalize(mrna_mat)
mirna_cpm <- sweep(mirna_mat, 2, colSums(mirna_mat), FUN = "/") * 1e6


sample_meta <- data.frame(
  condition = c("Resistant", "Resistant", "Sensitive", "Sensitive"),
  treatment = c("Control", "Drug10nM", "Control", "Drug10nM"),
  row.names = c("A375_RC", "A375_R10", "A375_SC", "A375_S10")
)

# Filter low-count genes: keep genes with CPM > 0.5 in at least 1 sample
mrna_cpm_check <- sweep(mrna_mat, 2, colSums(mrna_mat) / 1e6, FUN = "/")
keep_genes <- rowSums(mrna_cpm_check > 0.5) >= 1
mrna_filtered <- mrna_mat[keep_genes, ]




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
  scores <- scale(abs(full_results$ranking)) * full_results$prob^3
  full_results$priority_score <- as.numeric(scores)
  
  threshold     <- quantile(full_results$priority_score, 0.90)
  degs_priority <- full_results[full_results$priority_score >= threshold &
                                  full_results$prob >= q, ]
  
  message(sprintf("  [%s] DEGs (prob>=%.2f + top 10%% priority): %d",
                  label, q, nrow(degs_priority)))
  
  list(result   = res,
       degs     = degs_priority,
       full     = full_results,
       label    = label,
       samples_A = samples_A,
       samples_B = samples_B)
}

mean_counts <- rowMeans(mrna_mat)
mrna_filtered2 <- mrna_mat[mean_counts >= 10 & keep_genes, ]


# Comp1: SC vs RC (primary resistance comparison)

comp1.ora <- run_noiseq_comparison(
  mrna_filtered2,
  meta = sample_meta,
  samples_A = "A375_SC",
  samples_B = "A375_RC",
  label = "SC_vs_RC"
)

# Comp2: SC vs S10 (drug effect in sensitive)
comp2.ora <- run_noiseq_comparison(
  mrna_filtered2,
  meta = sample_meta,
  samples_A = "A375_SC",
  samples_B = "A375_S10",
  label = "SC_vs_S10"
)

# Comp3: RC vs R10 (drug effect in resistant)
comp3.ora <- run_noiseq_comparison(
  mrna_filtered2,
  meta = sample_meta,
  samples_A = "A375_RC",
  samples_B = "A375_R10",
  label = "RC_vs_R10"
)

comp4.ora <- run_noiseq_comparison(
  mrna_filtered2,
  meta = sample_meta,
  samples_A = "A375_S10",
  samples_B = "A375_R10",
  label = "S10_vs_R10"
)

write.table(rownames(comp1.ora$degs),
            file      = file.path(OUTPUT_DIR, "webg_ora1.txt"),
            row.names = FALSE,
            col.names = FALSE,
            quote     = FALSE)

#----------------------------ORA----KEGG------GOBP-----------------
bioc_needed <- c("clusterProfiler", "org.Hs.eg.db", "enrichplot", "STRINGdb", "ReactomePA")
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
  library(ReactomePA)
})

# ── Comparison registry ───────────────────────────────────────────────────────

oracomparisons <- list(
  list(obj = comp1.ora, label = "SC_vs_RC"),
  list(obj = comp2.ora, label = "SC_vs_S10"),
  list(obj = comp3.ora, label = "RC_vs_R10"),
  list(obj = comp4.ora, label = "S10_vs_R10")
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

for (comp in oracomparisons) {
  
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
  
  # ── ORA: Reactome ──────────────────────────────────────────────────────────
  
  message("  Running ORA Reactome...")
  ora_reactome <- tryCatch(
    ReactomePA::enrichPathway(
      gene          = deg_mapped$ENTREZID,
      organism      = "human",
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      qvalueCutoff  = 0.20,
      readable      = TRUE
    ),
    error = function(e) { message("  ORA Reactome failed: ", e$message); NULL }
  )
  
  if (!is.null(ora_reactome) && nrow(as.data.frame(ora_reactome)) > 0) {
    write.csv(as.data.frame(ora_reactome),
              file.path(OUTPUT_DIR, sprintf("ORA_Reactome_%s.csv", lbl)),
              row.names = FALSE)
    p_reactome <- enrichplot::dotplot(ora_reactome, showCategory = 20) +
      ggtitle(sprintf("ORA Reactome | %s", lbl))
    ggsave(file.path(OUTPUT_DIR, sprintf("ORA_Reactome_%s.pdf", lbl)),
           p_reactome, width = 10, height = 9)
    #assign(paste0("ora_reactome_", lbl), ora_reactome)
    message(sprintf("  ORA Reactome: %d pathways", nrow(as.data.frame(ora_reactome))))
  } else {
    message("  ORA Reactome: no significant pathways.")
  }
  
}