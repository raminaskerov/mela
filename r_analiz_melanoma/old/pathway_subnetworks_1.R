# =============================================================================
# PATHWAY-GUIDED SUBNETWORK RECONSTRUCTION
# Subset approach: filter existing full_network by GSEA leading edge genes
#
# REQUIRED SESSION OBJECTS:
#   full_network  — all anti-correlated miRNA→mRNA pairs (from Priority 1)
#   mirna_fc      — miRNA fold changes (SC→RC)
#
# REQUIRED FILES (in OUTPUT_DIR):
#   GSEA_KEGG_SC_vs_RC.csv
#   GSEA_Reactome_SC_vs_RC.csv
#   GSEA_KEGG_S10_vs_R10.csv
#   GSEA_Reactome_S10_vs_R10.csv
# =============================================================================
#Ilk koddan qalan
get_deg_df <- function(comp_obj, q = 0.8) {
  df <- as.data.frame(comp_obj$result@results[[1]])
  df$mirna <- rownames(df)
  colnames(df)[colnames(df) == "M"]    <- "log2FC"
  colnames(df)[colnames(df) == "prob"] <- "prob"
  df |> dplyr::filter(prob >= q) |>
    dplyr::mutate(direction = ifelse(log2FC > 0, "Up_in_B", "Down_in_B"))
}
# Get DEG fold changes for the primary resistance comparison (SC vs RC)
deg_mrna_SC_vs_RC <- get_deg_df(comp1)   # mRNA DEGs
mrna_fc <- deg_mrna_SC_vs_RC |>
  dplyr::select(gene = mirna, log2FC_mrna = log2FC, prob_mrna = prob)
all_deg_mirnas <- unique(c(mirna_deg1$mirna,
                           mirna_deg4$mirna))

# Get miRNA fold changes — use all queried miRNAs, not just NOISeq DEGs,
# because miR-140-3p borderline cases matter here
mirna_fc <- data.frame(
  mirna = rownames(mirna_mat),
  log2FC_mirna = log2((mirna_cpm[, "A375_RC"] + 1) /
                        (mirna_cpm[, "A375_SC"] + 1)),
  stringsAsFactors = FALSE
) |>
  dplyr::filter(mirna %in% all_deg_mirnas)  # keep only our DEG miRNAs

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2)
  library(ggraph); library(igraph); library(tidygraph)
  library(ggrepel); library(patchwork)
  library(org.Hs.eg.db); library(clusterProfiler)
  library(RColorBrewer); library(viridis)
})

OUTPUT_DIR <- "output"
dir.create(file.path(OUTPUT_DIR, "subnetworks"), showWarnings = FALSE)

# =============================================================================
# SECTION 1: LOAD & PARSE GSEA RESULTS
# =============================================================================

message("Loading GSEA results...")

gsea_files <- list(
  KEGG_SC_vs_RC      = file.path(OUTPUT_DIR, "GSEA_KEGG_SC_vs_RC.csv"),
  Reactome_SC_vs_RC  = file.path(OUTPUT_DIR, "GSEA_Reactome_SC_vs_RC.csv"),
  KEGG_S10_vs_R10    = file.path(OUTPUT_DIR, "GSEA_KEGG_S10_vs_R10.csv"),
  Reactome_S10_vs_R10= file.path(OUTPUT_DIR, "GSEA_Reactome_S10_vs_R10.csv")
)

gsea_all <- lapply(names(gsea_files), function(nm) {
  f <- gsea_files[[nm]]
  if (!file.exists(f)) {
    message(sprintf("  Missing: %s", f)); return(NULL)
  }
  df <- read.csv(f, stringsAsFactors = FALSE)
  df$source <- nm
  parts <- strsplit(nm, "_", fixed = TRUE)[[1]]
  df$db         <- parts[1]                          # KEGG / Reactome
  df$comparison <- paste(parts[-1], collapse = "_")  # SC_vs_RC etc.
  df
})
gsea_all <- do.call(rbind, Filter(Negate(is.null), gsea_all))
message(sprintf("  Total GSEA terms loaded: %d", nrow(gsea_all)))

# =============================================================================
# SECTION 2: PATHWAY SELECTION
# =============================================================================
# Hand-curated list covering:
#   (a) Top GSEA hits from SC_vs_RC
#   (b) Key cancer pathways absent from GSEA but biologically essential
# Each entry: label, ID pattern to match, source comparison, forced_include flag

pathway_catalog <- data.frame(
  label = c(
    # ── From your GSEA SC_vs_RC ──────────────────────────────────────────────
    "DNA Damage Response",
    "Cell Cycle / G2M",
    "Homologous Recombination",
    "OXPHOS / Complex I",
    "Iron Transport",
    "Interferon Signaling",
    "Lysosome / Autophagy",
    "Cytokine Signaling",
    # ── From S10 vs R10 ──────────────────────────────────────────────────────
    "Ferroptosis",
    # ── Canonical cancer pathways (add regardless of GSEA) ──────────────────
    "EMT / Mesenchymal",
    "p53 Signaling",
    "PI3K / AKT / mTOR",
    "Hypoxia / HIF1A",
    "WNT Signaling",
    "Apoptosis"
  ),
  # Search terms matched against Description column (case-insensitive)
  search_term = c(
    "DNA.*break|DNA.*repair|double.*strand|checkpoint",
    "cell cycle|G2.*M|mitotic|G1.*S",
    "homolog.*recomb|HDR|HRR",
    "oxidative phosphoryl|complex I|OXPHOS|electron transport",
    "iron|transferrin|ferritin|ferrop",
    "interferon|IFN|innate immune",
    "lysosom|autopha|vacuol",
    "cytokine|interleukin|IL-",
    "ferroptosis",
    "epithelial.*mesench|EMT|mesench.*transit",
    "p53|TP53|apoptosis.*p53",
    "PI3K|AKT|mTOR|PTEN",
    "hypoxia|HIF|oxygen",
    "Wnt|WNT|beta.catenin",
    "apoptosis|caspase"
  ),
  # Which comparison to pull leading edges from
  comparison = c(
    "SC_vs_RC", "SC_vs_RC", "SC_vs_RC", "SC_vs_RC",
    "SC_vs_RC", "SC_vs_RC", "SC_vs_RC", "SC_vs_RC",
    "S10_vs_R10",
    rep("SC_vs_RC", 6)   # canonical pathways: pull from SC_vs_RC if present
  ),
  # If TRUE: include even if not found in GSEA (use curated gene set instead)
  force_include = c(
    FALSE, FALSE, FALSE, FALSE,
    FALSE, FALSE, FALSE, FALSE,
    FALSE,
    TRUE, TRUE, TRUE, TRUE, TRUE, TRUE
  ),
  stringsAsFactors = FALSE
)

message(sprintf("  Pathway catalog: %d entries", nrow(pathway_catalog)))

# =============================================================================
# SECTION 3: ENTREZ → SYMBOL CONVERSION
# =============================================================================

# Build a global ENTREZ → SYMBOL lookup from all GSEA core_enrichment fields
all_entrez <- unique(unlist(
  strsplit(gsea_all$core_enrichment, "/", fixed = TRUE)
))
all_entrez <- all_entrez[all_entrez != "" & !is.na(all_entrez)]
message(sprintf("  Converting %d unique ENTREZ IDs to symbols...", length(all_entrez)))

entrez_map <- suppressMessages(
  clusterProfiler::bitr(
    all_entrez,
    fromType = "ENTREZID",
    toType   = "SYMBOL",
    OrgDb    = org.Hs.eg.db
  )
)
# entrez_map: ENTREZID, SYMBOL
message(sprintf("  Mapped: %d / %d IDs", nrow(entrez_map), length(all_entrez)))

# Helper: convert slash-separated ENTREZ string → vector of gene symbols
entrez_to_symbols <- function(entrez_str) {
  ids <- strsplit(entrez_str, "/", fixed = TRUE)[[1]]
  syms <- entrez_map$SYMBOL[entrez_map$ENTREZID %in% ids]
  unique(syms[!is.na(syms)])
}

# =============================================================================
# SECTION 4: CURATED GENE SETS FOR FORCED-INCLUDE PATHWAYS
# =============================================================================
# These are used when a pathway has no GSEA hit (force_include = TRUE).
# Gene symbols only — restricted to what's in your expression data.

canonical_genesets <- list(
  "EMT / Mesenchymal" = c(
    "ZEB1","ZEB2","TWIST1","TWIST2","SNAI1","SNAI2",
    "VIM","CDH1","CDH2","FN1","MMP2","MMP9",
    "TGFB1","TGFB2","TGFBR1","TGFBR2",
    "AXL","EGFR","NGFR","SOX10","MITF"
  ),
  "p53 Signaling" = c(
    "TP53","CDKN1A","MDM2","MDM4","BBC3","PUMA",
    "BAX","GADD45A","GADD45B","PERP","FAS",
    "CASP3","CASP7","CASP8","CASP9","CYCS",
    "RB1","E2F1","ATM","CHEK1","CHEK2"
  ),
  "PI3K / AKT / mTOR" = c(
    "PIK3CA","PIK3CB","PIK3CD","PIK3R1","PIK3R2",
    "AKT1","AKT2","AKT3","PTEN",
    "MTOR","RPTOR","RICTOR","RPS6KB1","RPS6KB2",
    "EIF4EBP1","TSC1","TSC2","DEPTOR",
    "PDK1","GSK3B","FOXO1","FOXO3"
  ),
  "Hypoxia / HIF1A" = c(
    "HIF1A","HIF1B","EPAS1","ARNT",
    "VEGFA","VEGFC","SLC2A1","LDHA","PKM",
    "PFKL","PFKP","ENO1","ENO2",
    "ADM","BNIP3","BNIP3L","CA9",
    "HMOX1","SLC7A11","NDRG1","P4HA1"
  ),
  "WNT Signaling" = c(
    "CTNNB1","APC","AXIN1","AXIN2",
    "GSK3B","CK1A","DVL1","DVL2","DVL3",
    "FZD1","FZD2","FZD3","FZD7",
    "TCF7","TCF7L2","LEF1","MYC",
    "CCND1","CD44","SURVIVIN","BIRC5"
  ),
  "Apoptosis" = c(
    "BCL2","BCL2L1","MCL1","BCL2L11",
    "BAX","BAK1","BAD","BID","BIM",
    "CASP3","CASP7","CASP8","CASP9","CASP6",
    "CYCS","APAF1","DIABLO","XIAP",
    "PARP1","FAS","FASLG","TRAIL","TNFRSF10A"
  )
)

# =============================================================================
# SECTION 5: BUILD PATHWAY GENE TABLE WITH REDUNDANCY REMOVAL
# =============================================================================
# For each pathway: find best GSEA hit (highest |NES|) + parse leading edges.
# Redundancy removal: if a gene appears in multiple pathways, assign it to the
# pathway with the highest |NES|. Ties broken by comparison priority (SC_vs_RC > S10_vs_R10).

message("Building pathway gene table with redundancy removal...")

pathway_gene_table <- list()

#pseudo full_network
#full_network <- target_df |>
  #dplyr::rename(mirna = mature_mirna_id, gene = target_symbol) |>
  #dplyr::inner_join(mirna_fc, by = "mirna") |>
  #dplyr::inner_join(mrna_fc,  by = "gene") |>
  #dplyr::filter(sign(log2FC_mirna) != sign(log2FC_mrna)) |>  # anti-correlation
  #dplyr::distinct(mirna, gene, .keep_all = TRUE)

for (i in seq_len(nrow(pathway_catalog))) {

  pw    <- pathway_catalog[i, ]
  label <- pw$label

  # ── Find GSEA hits matching this pathway ──────────────────────────────────
  hits <- gsea_all |>
    dplyr::filter(
      comparison == pw$comparison,
      grepl(pw$search_term, Description, ignore.case = TRUE, perl = TRUE)
    ) |>
    dplyr::arrange(desc(abs(NES)))

  if (nrow(hits) > 0) {
    # Best hit: highest |NES|
    best     <- hits[1, ]
    genes    <- entrez_to_symbols(best$core_enrichment)
    nes      <- best$NES
    padj     <- best$p.adjust
    term_id  <- best$ID
    term_desc<- best$Description
    source   <- "GSEA_leading_edge"
    message(sprintf("  [%s] GSEA hit: '%s' | NES=%.2f | %d genes",
                    label, term_desc, nes, length(genes)))
  } else if (pw$force_include && label %in% names(canonical_genesets)) {
    # Fallback: use curated gene set
    genes    <- canonical_genesets[[label]]
    nes      <- 0          # unknown from GSEA
    padj     <- NA_real_
    term_id  <- "curated"
    term_desc<- label
    source   <- "curated_geneset"
    message(sprintf("  [%s] No GSEA hit; using curated gene set (%d genes)",
                    label, length(genes)))
  } else {
    message(sprintf("  [%s] No GSEA hit and no curated set — skipping.", label))
    next
  }

  # Restrict to genes present in full_network targets
  #genes_in_network <- genes[genes %in% unique(full_network$gene)]
  #message(sprintf("    -> %d / %d genes present in full_network",
  #                length(genes_in_network), length(genes)))

  pathway_gene_table[[label]] <- data.frame(
    pathway      = label,
    gene         = genes,
    #in_network   = genes %in% unique(full_network$gene),
    NES          = nes,
    padj         = padj,
    term_id      = term_id,
    term_desc    = term_desc,
    source       = source,
    comparison   = pw$comparison,
    stringsAsFactors = FALSE
  )
}

pathway_gene_df <- dplyr::bind_rows(pathway_gene_table)

# ── Redundancy removal ────────────────────────────────────────────────────────
# For genes appearing in multiple pathways: keep assignment with highest |NES|.
# Curated sets (NES=0) lose to any GSEA hit.

pathway_gene_df <- pathway_gene_df |>
  dplyr::group_by(gene) |>
  dplyr::arrange(desc(abs(NES))) |>        # highest |NES| first
  dplyr::slice(1) |>                       # keep only one pathway per gene
  dplyr::ungroup()

message(sprintf("\nAfter redundancy removal: %d unique genes across %d pathways",
                nrow(pathway_gene_df),
                length(unique(pathway_gene_df$pathway))))
#print(table(pathway_gene_df$pathway[pathway_gene_df$in_network == TRUE]))

write.csv(pathway_gene_df,
          file.path(OUTPUT_DIR, "pathway_gene_assignments2.csv"),
          row.names = FALSE)

# =============================================================================
# SECTION 6: BUILD PATHWAY-ANNOTATED SUBNETWORKS
# =============================================================================

message("\nBuilding pathway-annotated subnetworks...")

# Genes assigned to a pathway and present in full_network
pathway_genes_network <- pathway_gene_df |>
  dplyr::filter(in_network) |>
  dplyr::select(gene, pathway, NES, source)

# Join pathway annotation onto full_network
pathway_network <- full_network |>
  dplyr::inner_join(
    pathway_genes_network |> dplyr::rename(gene_pathway = pathway,
                                            pathway_NES  = NES,
                                            gene_source  = source),
    by = "gene"
  )

message(sprintf("  Pathway-annotated edges: %d (from %d total)",
                nrow(pathway_network), nrow(full_network)))
message(sprintf("  miRNAs with pathway targets: %d",
                length(unique(pathway_network$mirna))))
message(sprintf("  Pathway coverage:"))
print(sort(table(pathway_network$gene_pathway), decreasing = TRUE))

# ── Coverage report: what fraction of each pathway is covered ────────────────
coverage_df <- pathway_gene_df |>
  dplyr::filter(in_network) |>
  dplyr::group_by(pathway) |>
  dplyr::summarise(
    n_genes_in_network = n(),
    .groups = "drop"
  ) |>
  dplyr::left_join(
    pathway_gene_df |>
      dplyr::group_by(pathway) |>
      dplyr::summarise(n_genes_total = n(), .groups = "drop"),
    by = "pathway"
  ) |>
  dplyr::mutate(
    coverage_pct = round(100 * n_genes_in_network / n_genes_total, 1)
  ) |>
  dplyr::left_join(
    pathway_gene_df |>
      dplyr::distinct(pathway, NES, padj, term_desc, source) |>
      dplyr::group_by(pathway) |>
      dplyr::slice(1),
    by = "pathway"
  ) |>
  dplyr::arrange(desc(n_genes_in_network))

write.csv(coverage_df,
          file.path(OUTPUT_DIR, "subnetworks", "pathway_network_coverage.csv"),
          row.names = FALSE)
message("\nPathway network coverage:")
print(coverage_df[, c("pathway","n_genes_total","n_genes_in_network","coverage_pct","NES")])

# =============================================================================
# SECTION 7: VISUALIZATION — COMBINED OVERVIEW PLOT
# =============================================================================

message("\nBuilding combined overview visualization...")

# ── Color palette ─────────────────────────────────────────────────────────────
pathway_colors <- c(
  "DNA Damage Response"   = "#E64B35",
  "Cell Cycle / G2M"      = "#F39B7F",
  "Homologous Recombination"= "#D62728",
  "OXPHOS / Complex I"    = "#9467BD",
  "Iron Transport"        = "#8C564B",
  "Interferon Signaling"  = "#4DBBD5",
  "Lysosome / Autophagy"  = "#00A087",
  "Cytokine Signaling"    = "#AEC7E8",
  "Ferroptosis"           = "#C5B0D5",
  "EMT / Mesenchymal"     = "#2CA02C",
  "p53 Signaling"         = "#FF7F0E",
  "PI3K / AKT / mTOR"    = "#1F77B4",
  "Hypoxia / HIF1A"       = "#BCBD22",
  "WNT Signaling"         = "#17BECF",
  "Apoptosis"             = "#7F7F7F",
  "miRNA"                 = "#7B2D8B"
)

# ── Build igraph ──────────────────────────────────────────────────────────────

# Node table — miRNAs + pathway-annotated genes
mirna_nodes <- data.frame(
  name      = unique(pathway_network$mirna),
  node_type = "miRNA",
  category  = "miRNA",
  log2FC    = mirna_fc$log2FC_mirna[match(unique(pathway_network$mirna), mirna_fc$mirna)],
  stringsAsFactors = FALSE
)

gene_nodes <- pathway_network |>
  dplyr::distinct(gene, gene_pathway, log2FC_mrna) |>
  dplyr::rename(name = gene, category = gene_pathway, log2FC = log2FC_mrna) |>
  dplyr::mutate(node_type = "mRNA")

all_nodes <- dplyr::bind_rows(mirna_nodes, gene_nodes)

g_pathway <- igraph::graph_from_data_frame(
  d = pathway_network[, c("mirna","gene","type","gene_pathway",
                           "log2FC_mirna","log2FC_mrna")],
  directed = TRUE,
  vertices = all_nodes
)

message(sprintf("  Combined pathway graph: %d nodes, %d edges",
                igraph::vcount(g_pathway), igraph::ecount(g_pathway)))

# ── Export combined graphml ───────────────────────────────────────────────────
igraph::write_graph(
  g_pathway,
  file   = file.path(OUTPUT_DIR, "subnetworks", "network_pathway_combined.graphml"),
  format = "graphml"
)
message("  Combined graphml exported.")

# ── ggraph plot ───────────────────────────────────────────────────────────────
# Label: all miRNAs + top-degree gene per pathway

gene_degree <- igraph::degree(g_pathway)
top_per_pathway <- pathway_network |>
  dplyr::mutate(degree = gene_degree[gene]) |>
  dplyr::group_by(gene_pathway) |>
  dplyr::slice_max(degree, n = 2, with_ties = FALSE) |>
  dplyr::pull(gene)

label_nodes_combined <- unique(c(
  unique(pathway_network$mirna),
  top_per_pathway,
  "FTH1", "NCOA4", "TFRC", "IREB2", "SLC7A11", "GPX4"  # always label iron genes
))

tg_combined <- tidygraph::as_tbl_graph(g_pathway) |>
  tidygraph::activate(nodes) |>
  dplyr::mutate(
    label  = ifelse(name %in% label_nodes_combined, name, NA_character_),
    degree = igraph::degree(g_pathway)
  )

p_combined <- ggraph::ggraph(tg_combined, layout = "stress") +
  ggraph::geom_edge_arc(
    aes(color = log2FC_mirna > 0),
    arrow       = arrow(length = unit(1.5, "mm"), type = "closed"),
    end_cap     = ggraph::circle(2.5, "mm"),
    alpha       = 0.35,
    linewidth   = 0.35,
    strength    = 0.12
  ) +
  ggraph::geom_node_point(
    aes(size  = ifelse(node_type == "miRNA", abs(log2FC) + 2, abs(log2FC) + 0.5),
        color = category,
        shape = node_type),
    alpha = 0.9
  ) +
  ggraph::geom_node_text(
    aes(label = label),
    size         = 2.5,
    repel        = TRUE,
    fontface     = "italic",
    max.overlaps = 30,
    segment.size = 0.2,
    segment.color= "grey50"
  ) +
  scale_color_manual(values = pathway_colors, name = "Pathway") +
  scale_shape_manual(values = c(miRNA = 18, mRNA = 16), name = "Node type") +
  scale_size_continuous(range = c(1.5, 7), name = "|log₂FC|") +
  ggraph::scale_edge_color_manual(
    values = c(`TRUE`  = "#B2182B", `FALSE` = "#2166AC"),
    labels = c("miRNA ↑ in RC", "miRNA ↓ in RC"),
    name   = "miRNA direction"
  ) +
  labs(
    title    = "Pathway-Guided miRNA–mRNA Network",
    subtitle = sprintf("%d miRNAs | %d pathway-annotated target genes | %d pathways",
                       sum(igraph::V(g_pathway)$node_type == "miRNA"),
                       sum(igraph::V(g_pathway)$node_type == "mRNA"),
                       length(unique(pathway_network$gene_pathway))),
    caption  = "Gene colors = pathway assignment (redundancy-removed by max|NES|) | Comparison: SC vs RC"
  ) +
  ggraph::theme_graph(base_family = "sans") +
  theme(
    plot.title    = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size  = 9),
    plot.caption  = element_text(size  = 7, color = "grey50"),
    legend.position = "right",
    legend.text   = element_text(size = 7),
    legend.title  = element_text(size = 8, face = "bold")
  )

ggplot2::ggsave(
  file.path(OUTPUT_DIR, "subnetworks", "network_pathway_combined.pdf"),
  p_combined, width = 16, height = 13,
  device = cairo_pdf
)
message("  Combined pathway network plot saved.")

# =============================================================================
# SECTION 8: PER-PATHWAY SUBNETWORK PLOTS
# =============================================================================
# Only plot pathways with ≥ 3 genes in the network (otherwise uninformative).

message("\nBuilding per-pathway subnetwork plots...")

pathways_to_plot <- coverage_df |>
  dplyr::filter(n_genes_in_network >= 3) |>
  dplyr::pull(pathway)

message(sprintf("  Pathways with ≥3 genes in network: %d", length(pathways_to_plot)))

for (pw_label in pathways_to_plot) {

  pw_edges <- pathway_network |>
    dplyr::filter(gene_pathway == pw_label)

  if (nrow(pw_edges) == 0) next

  # Node table for this subnetwork
  pw_mirnas <- data.frame(
    name      = unique(pw_edges$mirna),
    node_type = "miRNA",
    category  = "miRNA",
    log2FC    = mirna_fc$log2FC_mirna[match(unique(pw_edges$mirna), mirna_fc$mirna)],
    stringsAsFactors = FALSE
  )
  pw_genes <- pw_edges |>
    dplyr::distinct(gene, log2FC_mrna) |>
    dplyr::rename(name = gene, log2FC = log2FC_mrna) |>
    dplyr::mutate(node_type = "mRNA", category = pw_label)

  pw_nodes <- dplyr::bind_rows(pw_mirnas, pw_genes)

  g_pw <- igraph::graph_from_data_frame(
    d        = pw_edges[, c("mirna","gene","type","log2FC_mirna","log2FC_mrna")],
    directed = TRUE,
    vertices = pw_nodes
  )

  # Label all nodes in small subnetworks, top nodes in large ones
  n_nodes <- igraph::vcount(g_pw)
  label_all <- if (n_nodes <= 30) pw_nodes$name else {
    c(unique(pw_edges$mirna),
      names(sort(igraph::degree(g_pw), decreasing = TRUE))[1:min(15, n_nodes)])
  }

  tg_pw <- tidygraph::as_tbl_graph(g_pw) |>
    tidygraph::activate(nodes) |>
    dplyr::mutate(
      label  = ifelse(name %in% label_all, name, NA_character_),
      degree = igraph::degree(g_pw)
    )

  # Choose color for this pathway's gene nodes
  pw_color <- pathway_colors[pw_label]
  if (is.na(pw_color)) pw_color <- "#8491B4"

  p_pw <- ggraph::ggraph(tg_pw, layout = ifelse(n_nodes <= 20, "stress", "fr")) +
    ggraph::geom_edge_arc(
      aes(color = log2FC_mirna > 0,
          linetype = type),
      arrow     = arrow(length = unit(1.5, "mm"), type = "closed"),
      end_cap   = ggraph::circle(2.5, "mm"),
      alpha     = 0.55,
      linewidth = 0.45,
      strength  = 0.15
    ) +
    ggraph::geom_node_point(
      aes(size  = degree + 1,
          shape = node_type,
          color = ifelse(node_type == "miRNA", "miRNA", pw_label)),
      alpha = 0.92
    ) +
    ggraph::geom_node_text(
      aes(label = label),
      size         = 2.7,
      repel        = TRUE,
      fontface     = "italic",
      max.overlaps = 25,
      segment.size = 0.2
    ) +
    scale_color_manual(
      values = setNames(
        c("#7B2D8B", pw_color),
        c("miRNA",   pw_label)
      ),
      name = "Node type"
    ) +
    scale_shape_manual(values = c(miRNA = 18, mRNA = 16), guide = "none") +
    scale_size_continuous(range = c(2.5, 9), name = "Degree") +
    ggraph::scale_edge_color_manual(
      values = c(`TRUE`  = "#B2182B", `FALSE` = "#2166AC"),
      labels = c("miRNA ↑ in RC", "miRNA ↓ in RC"),
      name   = "miRNA FC"
    ) +
    ggraph::scale_edge_linetype_manual(
      values = c(validated = "solid", predicted = "dashed"),
      name   = "Evidence"
    ) +
    labs(
      title    = sprintf("%s — miRNA Regulatory Subnetwork", pw_label),
      subtitle = sprintf("%d miRNAs → %d genes | %d edges (%d validated, %d predicted)",
                         sum(igraph::V(g_pw)$node_type == "miRNA"),
                         sum(igraph::V(g_pw)$node_type == "mRNA"),
                         nrow(pw_edges),
                         sum(pw_edges$type == "validated"),
                         sum(pw_edges$type == "predicted")),
      caption  = sprintf("Leading edge from GSEA: SC vs RC | Comparison: %s",
                         unique(pw_edges$gene_pathway)[1])
    ) +
    ggraph::theme_graph(base_family = "sans") +
    theme(
      plot.title    = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 8),
      plot.caption  = element_text(size = 7, color = "grey50"),
      legend.position = "bottom",
      legend.text   = element_text(size = 7)
    )

  fname <- gsub("[/ ]", "_", pw_label)
  ggplot2::ggsave(
    file.path(OUTPUT_DIR, "subnetworks", sprintf("subnetwork_%s.pdf", fname)),
    p_pw,
    width  = 10,
    height = 8,
    device = cairo_pdf
  )

  # Export per-pathway graphml
  igraph::write_graph(
    g_pw,
    file   = file.path(OUTPUT_DIR, "subnetworks",
                       sprintf("subnetwork_%s.graphml", fname)),
    format = "graphml"
  )

  message(sprintf("  [%s] %d nodes, %d edges — saved.", pw_label, n_nodes, nrow(pw_edges)))
}

# =============================================================================
# SECTION 9: SUMMARY TABLE — miRNA → PATHWAY REGULATORY MAP
# =============================================================================
# One row per miRNA×pathway pair: how many genes does each miRNA regulate
# in each pathway? This is the key table for biological interpretation.

message("\nBuilding miRNA × pathway regulatory map...")

mirna_pathway_map <- pathway_network |>
  dplyr::group_by(mirna, gene_pathway) |>
  dplyr::summarise(
    n_targets         = n(),
    n_validated       = sum(type == "validated"),
    target_genes      = paste(sort(gene), collapse = "; "),
    mirna_log2FC      = unique(log2FC_mirna),
    mirna_direction   = ifelse(unique(log2FC_mirna) > 0, "Up_in_RC", "Down_in_RC"),
    .groups           = "drop"
  ) |>
  dplyr::arrange(gene_pathway, desc(n_targets))

write.csv(mirna_pathway_map,
          file.path(OUTPUT_DIR, "subnetworks", "mirna_pathway_regulatory_map.csv"),
          row.names = FALSE)

message("\n  miRNA × Pathway regulatory map (top entries):")
print(head(mirna_pathway_map[, c("mirna","gene_pathway","n_targets",
                                  "n_validated","mirna_direction")], 20))

# ── Heatmap: miRNAs (rows) × pathways (cols), value = n_targets ──────────────

map_wide <- mirna_pathway_map |>
  dplyr::select(mirna, gene_pathway, n_targets) |>
  tidyr::pivot_wider(names_from = gene_pathway, values_from = n_targets,
                     values_fill = 0L) |>
  as.data.frame()
rownames(map_wide) <- map_wide$mirna
map_wide$mirna     <- NULL

# Only keep pathways with ≥1 target for ≥1 miRNA
map_mat <- as.matrix(map_wide)
keep_cols <- colSums(map_mat) > 0
keep_rows <- rowSums(map_mat) > 0
map_mat   <- map_mat[keep_rows, keep_cols, drop = FALSE]

if (nrow(map_mat) >= 2 && ncol(map_mat) >= 2) {
  # Row annotation: miRNA direction
  mirna_dir_df <- mirna_pathway_map |>
    dplyr::distinct(mirna, mirna_direction) |>
    dplyr::filter(mirna %in% rownames(map_mat)) |>
    tibble::column_to_rownames("mirna")

  pheatmap::pheatmap(
    map_mat,
    color             = colorRampPalette(c("white","#FDD0A2","#E64B35"))(50),
    annotation_row    = mirna_dir_df,
    annotation_colors = list(
      mirna_direction = c(Up_in_RC = "#B2182B", Down_in_RC = "#2166AC")
    ),
    cluster_rows      = TRUE,
    cluster_cols      = TRUE,
    fontsize_row      = 7,
    fontsize_col      = 8,
    display_numbers   = TRUE,
    number_format     = "%d",
    number_color      = "black",
    main              = "miRNA → Pathway Regulatory Map\n(n_targets per pathway)",
    filename          = file.path(OUTPUT_DIR, "subnetworks",
                                   "heatmap_mirna_pathway_map.pdf"),
    width = 10, height = max(5, nrow(map_mat) * 0.3 + 3)
  )
  message("  miRNA × pathway heatmap saved.")
}

# =============================================================================
# FINAL SUMMARY
# =============================================================================

message("\n", paste(rep("=", 60), collapse = ""))
message("PATHWAY SUBNETWORK ANALYSIS COMPLETE")
message(paste(rep("=", 60), collapse = ""))
message(sprintf("  Pathways analyzed:        %d",
                length(unique(pathway_network$gene_pathway))))
message(sprintf("  Total pathway edges:       %d", nrow(pathway_network)))
message(sprintf("  miRNAs with pw targets:    %d",
                length(unique(pathway_network$mirna))))
message(sprintf("  Pathway-annotated genes:   %d",
                length(unique(pathway_network$gene))))
message(sprintf("\nOutputs in: %s/subnetworks/", OUTPUT_DIR))
message("  network_pathway_combined.pdf/.graphml   — all pathways overlay")
message("  subnetwork_[pathway].pdf/.graphml        — per-pathway plots")
message("  mirna_pathway_regulatory_map.csv         — miRNA × pathway table")
message("  heatmap_mirna_pathway_map.pdf            — regulatory map heatmap")
message("  pathway_gene_assignments.csv             — gene→pathway (redundancy-removed)")
message("  pathway_network_coverage.csv             — GSEA coverage per pathway")
