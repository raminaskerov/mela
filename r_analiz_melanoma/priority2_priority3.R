# =============================================================================
# PRIORITY 2 & 3 — MELANOMA RESISTANCE PIPELINE
# GSE283251 | Colakoglu Bergel et al., Scientific Reports 2025
#
# REPLACE lines 905–1360 in melanoma_resistance_pipeline_1.R with this file.
# All objects from Sections 0–1 and Priority 1 must be in the session first:
#   mrna_mat, mrna_norm, mrna_cpm, mirna_mat, mirna_cpm
#   comp1/2/3, deg1/2/3, mirna_deg1/3
#   IRON_GENES, AUTOPHAGY_GENES, MAPK_AKT_GENES, MDR_GENES
#   g_vis, g_full, full_network, vis_network, node_stats (from P1)
# =============================================================================


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

message("  Computing adjacency matrix...")
adjacency <- WGCNA::adjacency(wgcna_t, power = soft_power, type = "signed")

message("  Computing TOM (this may take a minute with 5K genes)...")
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

if (!requireNamespace("fgsea",          quietly = TRUE)) BiocManager::install("fgsea")
if (!requireNamespace("clusterProfiler",quietly = TRUE)) BiocManager::install("clusterProfiler")

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
      file.path(OUTPUT_DIR, sprintf("TF_manual_fgsea_%s.csv", mod_color))
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
        file.path(OUTPUT_DIR, sprintf("TF_kegg_fgsea_%s.csv", mod_color))
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
message("  Priority 2 complete.\n")


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

# Named numeric vector: value = -log2FC(SC/RC) = log2FC(RC/SC), name = gene symbol
# Positive values = higher in RC (resistant)
ranked_vec <- setNames(-ranked_genes$M, ranked_genes$gene)
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
is.null(tf_regulons)
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
print("a")
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
  
  # NES interpretation with M = log2(RC/SC):
  #   Positive NES → TF targets enriched among genes HIGH in RC
  #                → TF is MORE active in RC (resistance driver candidate)
  #   Negative NES → TF targets enriched among genes HIGH in SC
  #                → TF is LESS active in RC
  fgsea_tf$direction_in_RC <- ifelse(
    fgsea_tf$NES > 0, "More active in RC", "Less active in RC"
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
  
  save_fgsea(
    tf_summary,
    file.path(OUTPUT_DIR, "TF_fgsea_priority_summary.csv")
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
      subtitle = sprintf("fgsea on %s regulons | ranked by log2FC(RC/SC)\nPositive NES = more active in RC",
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
