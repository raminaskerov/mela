# =============================================================================
# 07_WGCNA.R
# INTEGRATIVE ANALYSIS: ENCORAFENIB RESISTANCE IN MALIGNANT MELANOMA
#
# Input:   checkpoint_03.RData
# Output:  WGCNA module tables/plots, checkpoint_07.RData
#
# This stage is intentionally independent of TF enrichment.
# The legacy TF linkage step has been moved to 08.
# =============================================================================


# ── Configuration ─────────────────────────────────────────────────────────────

DATA_DIR   <- "/home/ramin/mela/data/ilk"  # <-- UPDATE THIS PATH
OUTPUT_DIR <- file.path(DATA_DIR, "newoutput")
WGCNA_DIR  <- file.path(OUTPUT_DIR, "wgcna")
dir.create(WGCNA_DIR, showWarnings = FALSE, recursive = TRUE)


# ── Packages ──────────────────────────────────────────────────────────────────

cran_pkgs <- c("dplyr", "ggplot2", "pheatmap", "flashClust", "stringr")
bioc_pkgs  <- c("WGCNA")

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
  library(pheatmap)
  library(WGCNA)
  library(flashClust)
  library(stringr)
})


# ── Load checkpoint ───────────────────────────────────────────────────────────

load(file.path(DATA_DIR, "checkpoint_03.RData"))


# ── Gene sets ─────────────────────────────────────────────────────────────────

IRON_GENES <- c("NCOA4", "FTH1", "TFRC", "SLC7A11", "GPX4", "IREB2",
                "FTL", "HAMP", "SLC40A1", "HMOX1", "CYBRD1", "STEAP3")
AUTOPHAGY_GENES <- c("BECN1", "MAP1LC3B", "ATG5", "ATG7", "ATG12",
                     "SQSTM1", "ULK1", "WIPI2", "LAMP2", "RAB7A")
MAPK_AKT_GENES <- c("BRAF", "KRAS", "NRAS", "MAP2K1", "MAPK1", "MAPK3",
                    "AKT1", "AKT2", "PIK3CA", "PIK3R1", "PTEN",
                    "DUSP6", "DUSP4", "RAF1", "CRAF")
STEM_GENES <- unique(c(
  "POU5F1", "SOX2", "NANOG", "KLF4", "MYC", "CD44",
  "PROM1", "ALDH1A1", "ABCB5", "KDM5B",
  "ZEB1", "TWIST1", "SNAI1", "SNAI2", "VIM",
  "NOTCH1", "WNT5A", "AXL", "MITF", "OCT4", "CD133"
))


# =============================================================================
# PRIORITY 2: CO-EXPRESSION / WGCNA MODULE ANALYSIS
# =============================================================================

message("\n=== PRIORITY 2: Co-expression / WGCNA ===")

mrna_tmm_log <- log2(mrna_norm + 1)

gene_min_expr <- apply(mrna_tmm_log, 1, min)
gene_var_all  <- apply(mrna_tmm_log, 1, var)
keep_base     <- gene_min_expr > 0 & gene_var_all > 0

top_n_genes <- min(5000L, sum(keep_base))
var_rank    <- rank(-gene_var_all[keep_base], ties.method = "first")
wgcna_mat   <- mrna_tmm_log[keep_base, ][var_rank <= top_n_genes, , drop = FALSE]
message(sprintf("  Genes entering WGCNA: %d (top %d by variance)",
                nrow(wgcna_mat), top_n_genes))

wgcna_t <- t(wgcna_mat)

gsg <- WGCNA::goodSamplesGenes(wgcna_t, verbose = 0)
if (!gsg$allOK) {
  wgcna_t <- wgcna_t[gsg$goodSamples, gsg$goodGenes, drop = FALSE]
  wgcna_mat <- t(wgcna_t)
  message(sprintf("  WGCNA QC removed %d genes and %d samples",
                  sum(!gsg$goodGenes), sum(!gsg$goodSamples)))
}
message(sprintf("  WGCNA matrix after QC: %d genes x %d samples",
                ncol(wgcna_t), nrow(wgcna_t)))

powers  <- 1:20
sft_out <- WGCNA::pickSoftThreshold(
  wgcna_t,
  powerVector = powers,
  networkType = "signed",
  RsquaredCut = 0.80,
  verbose = 0
)

soft_power <- sft_out$powerEstimate
if (is.na(soft_power)) {
  soft_power <- 6L
  message("  Soft power undetermined; using default = 6")
} else {
  message(sprintf("  Soft power selected: %d", soft_power))
}

if ("enableWGCNAThreads" %in% getNamespaceExports("WGCNA")) {
  WGCNA::enableWGCNAThreads()
}

pdf(file.path(WGCNA_DIR, "wgcna_soft_threshold.pdf"), width = 9, height = 4)
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

adjacency <- WGCNA::adjacency(wgcna_t, power = soft_power, type = "signed")
TOM     <- WGCNA::TOMsimilarity(adjacency)
dissTOM <- 1 - TOM
rownames(dissTOM) <- colnames(dissTOM) <- colnames(wgcna_t)

gene_tree <- flashClust::flashClust(as.dist(dissTOM), method = "average")

raw_modules <- cutreeDynamic(
  dendro = gene_tree,
  distM = dissTOM,
  deepSplit = 2,
  pamRespectsDendro = FALSE,
  minClusterSize = 10
)
raw_colors <- WGCNA::labels2colors(raw_modules)
message(sprintf("  Raw modules before merging: %d (+grey)",
                length(unique(raw_colors[raw_colors != "grey"]))))

merge_res <- WGCNA::mergeCloseModules(
  wgcna_t,
  colors = raw_colors,
  cutHeight = 0.15,
  verbose = 0
)
module_colors <- merge_res$colors
names(module_colors) <- colnames(wgcna_t)
MEs <- WGCNA::orderMEs(merge_res$newMEs)

n_mods <- length(unique(module_colors[module_colors != "grey"]))
message(sprintf("  Modules after merging: %d", n_mods))

pdf(file.path(WGCNA_DIR, "wgcna_dendrogram_modules.pdf"), width = 11, height = 5)
WGCNA::plotDendroAndColors(
  dendro = gene_tree,
  colors = cbind(raw_colors, module_colors),
  groupLabels = c("Before merge", "After merge"),
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05,
  marAll = c(0, 5, 2, 0),
  main = sprintf(
    "Gene Co-expression Dendrogram | Signed WGCNA | n=4 exploratory\n%d modules | soft power = %d | top %d genes by variance",
    n_mods, soft_power, nrow(wgcna_mat)
  )
)
dev.off()

trait_mat <- matrix(
  c(1, 1, 0, 0,
    0, 1, 0, 1,
    1, 0, 0, 0,
    0, 0, 1, 0),
  nrow = 4, ncol = 4,
  dimnames = list(
    c("A375_RC", "A375_R10", "A375_SC", "A375_S10"),
    c("Resistance", "DrugTreatment", "RC_specific", "SC_specific")
  )
)

shared_samples <- intersect(rownames(MEs), rownames(trait_mat))
MEs_aligned   <- MEs[shared_samples, , drop = FALSE]
trait_aligned <- trait_mat[shared_samples, , drop = FALSE]

module_trait_cor  <- cor(MEs_aligned, trait_aligned, use = "pairwise.complete.obs")
module_trait_pval <- WGCNA::corPvalueStudent(module_trait_cor, nSamples = nrow(MEs_aligned))

text_mat <- matrix(
  paste0(round(module_trait_cor, 2), "\n(", signif(module_trait_pval, 1), ")"),
  nrow = nrow(module_trait_cor),
  ncol = ncol(module_trait_cor),
  dimnames = dimnames(module_trait_cor)
)

pdf(file.path(WGCNA_DIR, "wgcna_module_trait_correlation.pdf"),
    width = 8,
    height = max(4, nrow(module_trait_cor) * 0.45 + 2))
WGCNA::labeledHeatmap(
  Matrix = module_trait_cor,
  xLabels = colnames(trait_aligned),
  yLabels = rownames(module_trait_cor),
  ySymbols = rownames(module_trait_cor),
  colorLabels = FALSE,
  colors = WGCNA::blueWhiteRed(50),
  textMatrix = text_mat,
  setStdMargins = FALSE,
  cex.text = 0.72,
  zlim = c(-1, 1),
  main = "Module-Trait Correlation\n(r, p-value from Student t | n=4 exploratory)"
)
dev.off()

resistance_modules <- rownames(module_trait_cor)[abs(module_trait_cor[, "Resistance"]) > 0.7]
if (length(resistance_modules) == 0) {
  top3_idx <- order(abs(module_trait_cor[, "Resistance"]), decreasing = TRUE)[1:min(3L, nrow(module_trait_cor))]
  resistance_modules <- rownames(module_trait_cor)[top3_idx]
  message(sprintf("  No module passed |r| > 0.7; using top %d by |r| with Resistance",
                  length(resistance_modules)))
}
message(sprintf("  Resistance-associated modules: %s",
                paste(resistance_modules, collapse = ", ")))

kME_all <- cor(wgcna_t, MEs_aligned, use = "pairwise.complete.obs")

hub_genes  <- list()
hub_all_df <- data.frame(
  module = character(),
  gene = character(),
  kME = numeric(),
  iron_related = logical(),
  stem_related = logical(),
  mapk_related = logical(),
  stringsAsFactors = FALSE
)

for (mod_ME in resistance_modules) {
  mod_color <- gsub("^ME", "", mod_ME)
  mod_genes <- names(module_colors)[module_colors == mod_color]
  if (length(mod_genes) < 3 || !mod_ME %in% colnames(kME_all)) {
    next
  }

  kME_mod <- sort(kME_all[mod_genes, mod_ME], decreasing = TRUE)
  top_n <- min(30L, length(kME_mod))
  hub_genes[[mod_color]] <- names(kME_mod)[seq_len(top_n)]

  hub_all_df <- dplyr::bind_rows(hub_all_df, data.frame(
    module = mod_color,
    gene = names(kME_mod),
    kME = as.numeric(kME_mod),
    iron_related = names(kME_mod) %in% IRON_GENES,
    stem_related = names(kME_mod) %in% STEM_GENES,
    mapk_related = names(kME_mod) %in% MAPK_AKT_GENES,
    stringsAsFactors = FALSE
  ))
}

write.csv(hub_all_df, file.path(WGCNA_DIR, "wgcna_hub_genes.csv"), row.names = FALSE)
message(sprintf("  Hub genes populated for %d modules", length(hub_genes)))

bg_size <- ncol(wgcna_t)
genesets_for_overlap <- list(
  Iron_Ferroptosis = IRON_GENES,
  Stemness = STEM_GENES,
  Autophagy = AUTOPHAGY_GENES,
  MAPK_AKT = MAPK_AKT_GENES
)

overlap_rows <- list()
for (mod_ME in resistance_modules) {
  mod_color <- gsub("^ME", "", mod_ME)
  mod_genes <- names(module_colors)[module_colors == mod_color]

  for (gs_name in names(genesets_for_overlap)) {
    gs <- genesets_for_overlap[[gs_name]]
    hits <- mod_genes[mod_genes %in% gs]
    k <- length(hits)
    K <- length(gs)
    N <- bg_size
    n <- length(mod_genes)
    pval <- phyper(k - 1L, K, N - K, n, lower.tail = FALSE)
    OR <- (k / n) / (K / N)

    overlap_rows[[length(overlap_rows) + 1]] <- data.frame(
      module = mod_color,
      gene_set = gs_name,
      module_size = n,
      geneset_size = K,
      overlap_n = k,
      overlap_genes = paste(hits, collapse = "; "),
      p_value = pval,
      odds_ratio = OR,
      stringsAsFactors = FALSE
    )
  }
}

if (length(overlap_rows) > 0) {
  overlap_df <- dplyr::bind_rows(overlap_rows) |>
    dplyr::mutate(padj_BH = p.adjust(p_value, method = "BH")) |>
    dplyr::arrange(p_value)
  write.csv(overlap_df, file.path(WGCNA_DIR, "wgcna_module_genesets_overlap.csv"), row.names = FALSE)
} else {
  overlap_df <- data.frame()
}

for (mod_color in names(hub_genes)) {
  hubs <- hub_genes[[mod_color]]
  hubs_present <- hubs[hubs %in% rownames(mrna_mat)]
  if (length(hubs_present) < 3) next

  cpm_h <- sweep(mrna_mat[hubs_present, , drop = FALSE], 2, colSums(mrna_mat) / 1e6, FUN = "/")
  z_h   <- t(scale(t(log2(cpm_h + 1))))
  z_h[is.nan(z_h)] <- 0

  row_ann <- data.frame(
    Category = dplyr::case_when(
      hubs_present %in% IRON_GENES      ~ "Iron/Ferroptosis",
      hubs_present %in% STEM_GENES      ~ "Stemness",
      hubs_present %in% MAPK_AKT_GENES  ~ "MAPK/AKT",
      hubs_present %in% AUTOPHAGY_GENES ~ "Autophagy",
      TRUE                              ~ "Other"
    ),
    row.names = hubs_present
  )

  col_ann <- data.frame(
    Condition = c("Resistant", "Resistant", "Sensitive", "Sensitive"),
    Treatment = c("Control", "Drug", "Control", "Drug"),
    row.names = colnames(mrna_mat)
  )

  hub_colors <- list(
    Category = c("Iron/Ferroptosis" = "#E64B35", "Stemness" = "#7B2D8B",
                 "MAPK/AKT" = "#00A087", "Autophagy" = "#4DBBD5",
                 "Other" = "#8491B4"),
    Condition = c(Resistant = "#CC0000", Sensitive = "#0066CC"),
    Treatment = c(Control = "#999999", Drug = "#FF9900")
  )

  pheatmap::pheatmap(
    z_h,
    color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
    annotation_row = row_ann,
    annotation_col = col_ann,
    annotation_colors = hub_colors,
    cluster_cols = FALSE,
    cluster_rows = TRUE,
    show_rownames = TRUE,
    fontsize_row = 7,
    main = sprintf("Module '%s' - Top Hub Genes\nZ-scored log2CPM | Resistance-associated", mod_color),
    filename = file.path(WGCNA_DIR, sprintf("wgcna_heatmap_module_%s.pdf", mod_color)),
    width = 7,
    height = max(5, length(hubs_present) * 0.22 + 2)
  )
}

module_assign_df <- data.frame(
  gene = colnames(wgcna_t),
  module = module_colors,
  in_resistance_module = module_colors %in% gsub("^ME", "", resistance_modules),
  iron_related = colnames(wgcna_t) %in% IRON_GENES,
  stem_related = colnames(wgcna_t) %in% STEM_GENES,
  mapk_related = colnames(wgcna_t) %in% MAPK_AKT_GENES,
  stringsAsFactors = FALSE
) |>
  dplyr::arrange(module, gene)

write.csv(module_assign_df, file.path(WGCNA_DIR, "wgcna_gene_module_assignments.csv"), row.names = FALSE)


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
  wgcna_mat, wgcna_t, soft_power, module_colors, MEs,
  resistance_modules, hub_genes, hub_all_df, overlap_df, module_assign_df,
  file = file.path(DATA_DIR, "checkpoint_07.RData")
)

message("\ncheckpoint_07.RData saved.")
