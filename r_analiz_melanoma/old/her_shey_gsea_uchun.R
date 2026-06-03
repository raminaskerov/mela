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
  
  message(sprintf("  [%s] DEGs (prob>=%.2f): %d",
                  label, q, nrow(degs)))
  
  list(result   = res,
       degs     = degs,
       full     = full_results,
       label    = label,
       samples_A = samples_A,
       samples_B = samples_B)
}


comp1.gsea <- run_noiseq(
  mrna_filtered,
  meta = sample_meta,
  samples_A = "A375_SC",
  samples_B = "A375_RC",
  label = "SC_vs_RC"
)

# Comp2: SC vs S10 (drug effect in sensitive)
comp2.gsea <- run_noiseq(
  mrna_filtered,
  meta = sample_meta,
  samples_A = "A375_SC",
  samples_B = "A375_S10",
  label = "SC_vs_S10"
)

# Comp3: RC vs R10 (drug effect in resistant)
comp3.gsea <- run_noiseq(
  mrna_filtered,
  meta = sample_meta,
  samples_A = "A375_RC",
  samples_B = "A375_R10",
  label = "RC_vs_R10"
)

comp4.gsea <- run_noiseq(
  mrna_filtered,
  meta = sample_meta,
  samples_A = "A375_S10",
  samples_B = "A375_R10",
  label = "S10_vs_R10"
)



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

gseacomparisons <- list(
  list(obj = comp1.gsea, label = "SC_vs_RC"),
  list(obj = comp2.gsea, label = "SC_vs_S10"),
  list(obj = comp3.gsea, label = "RC_vs_R10"),
  list(obj = comp4.gsea, label = "S10_vs_R10")
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
comparisons_gsea <- list(
  list(comp=comp1.gsea, label = "SC_vs_RC",  A = "A375_SC", B = "A375_RC"),
  list(comp=comp2.gsea, label = "SC_vs_S10", A = "A375_SC", B = "A375_S10"),
  list(comp=comp3.gsea, label = "RC_vs_R10", A = "A375_RC", B = "A375_R10"),
  list(comp=comp4.gsea, label = "S10_vs_R10",A = "A375_S10",B = "A375_R10")
)

edger_webgsea <- list(
  list(comp=sim_comp4, label = "SC_vs_RC",  A = "A375_SC", B = "A375_RC")
  #list(comp=edger_comp2, label = "SC_vs_S10", A = "A375_SC", B = "A375_S10"),
  #list(comp=edger_comp3, label = "RC_vs_R10", A = "A375_RC", B = "A375_R10"),
  #list(comp=edger_comp4, label = "S10_vs_R10",A = "A375_S10",B = "A375_R10")
)

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
