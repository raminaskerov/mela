

# ── P1.2  Retrieve targets via multiMiR ──────────────────────────────────────

message("Querying multiMiR for targets...")

# Check your multiMiR version first — API differs between versions
packageVersion("multiMiR")

# Get validated interactions only first (highest confidence)
mirna_targets_validated <- multiMiR::get_multimir(
  mirna   = all_deg_mirnas,
  table   = "validated",    # miRTarBase: experimentally confirmed
  summary = FALSE
)

# Get predicted interactions separately (TargetScan + miRanda + other)
mirna_targets_predicted <- multiMiR::get_multimir(
  mirna                 = all_deg_mirnas,
  table                 = "predicted",
  predicted.cutoff      = 30,       # top 30 percentile by score
  predicted.cutoff.type = "p",      # "p" = percentile, "n" = top N per miRNA
  predicted.site        = "all",    # seed + non-seed sites
  summary               = FALSE
)

# Access the data slot directly (works in all multiMiR versions)
val_df  <- mirna_targets_validated@data
pred_df <- mirna_targets_predicted@data

message(sprintf("  Validated interactions:  %d", nrow(val_df)))
message(sprintf("  Predicted interactions:  %d", nrow(pred_df)))

# Combine and standardize column names
# multiMiR returns: mature_mirna_id, target_symbol, target_entrez,
#                   target_ensembl, experiment, support_type, pubmed_id
target_df <- dplyr::bind_rows(
  val_df  |> dplyr::mutate(type = "validated"),
  pred_df |> dplyr::mutate(type = "predicted")
) |>
  dplyr::select(
    mature_mirna_id, 
    target_symbol,
    type
  ) |>
  dplyr::filter(
    !is.na(target_symbol),
    target_symbol != ""
  ) |>
  dplyr::distinct(mature_mirna_id, target_symbol, .keep_all = TRUE)

message(sprintf("  Total unique pairs after dedup: %d", nrow(target_df)))
message(sprintf("  Targets per miRNA (median): %.0f",
                median(table(target_df$mature_mirna_id))))


# ── Add this AFTER target_df is built, BEFORE P1.3 ──────────────────────────

# Pre-filter 1: keep only targets actually expressed in your data
expressed_genes <- rownames(mrna_mat)[rowSums(mrna_mat) > 0]
target_df <- target_df |>
  dplyr::filter(target_symbol %in% expressed_genes)
message(sprintf("  After expression filter: %d pairs", nrow(target_df)))

# Pre-filter 2: validated interactions get priority —
# if a pair has both validated + predicted, keep as validated
target_df <- target_df |>
  dplyr::group_by(mature_mirna_id, target_symbol) |>
  dplyr::summarise(
    type = ifelse("validated" %in% type, "validated", "predicted"),
    .groups = "drop"
  )
message(sprintf("  After dedup with validation priority: %d pairs", nrow(target_df)))

# ── P1.3  Anti-correlation filter ─────────────────────────────────────────────

message("Applying anti-correlation filter...")

# Get DEG fold changes for the primary resistance comparison (SC vs RC)
deg_mrna_SC_vs_RC <- get_deg_df(comp1)   # mRNA DEGs
mrna_fc <- deg_mrna_SC_vs_RC |>
  dplyr::select(gene = mirna, log2FC_mrna = log2FC, prob_mrna = prob)

# Get miRNA fold changes — use all queried miRNAs, not just NOISeq DEGs,
# because miR-140-3p borderline cases matter here
mirna_fc <- data.frame(
  mirna = rownames(mirna_mat),
  log2FC_mirna = log2((mirna_cpm[, "A375_RC"] + 1) /
                        (mirna_cpm[, "A375_SC"] + 1)),
  stringsAsFactors = FALSE
) |>
  dplyr::filter(mirna %in% all_deg_mirnas)  # keep only our DEG miRNAs

# Join everything and apply anti-correlation
full_network <- target_df |>
  dplyr::rename(mirna = mature_mirna_id, gene = target_symbol) |>
  dplyr::inner_join(mirna_fc, by = "mirna") |>
  dplyr::inner_join(mrna_fc,  by = "gene") |>
  dplyr::filter(sign(log2FC_mirna) != sign(log2FC_mrna)) |>  # anti-correlation
  dplyr::distinct(mirna, gene, .keep_all = TRUE)

message(sprintf("  Full anti-correlated network: %d pairs", nrow(full_network)))
message(sprintf("  Validated pairs:              %d", sum(full_network$type == "validated")))
message(sprintf("  miRNAs with targets:          %d", length(unique(full_network$mirna))))
message(sprintf("  Unique target genes:          %d", length(unique(full_network$gene))))


# ── P1.4  Build two-tier network ──────────────────────────────────────────────

# Define key gene sets for biological focus

# ── Define gene sets (run before ALL_KEY_GENES) ───────────────────────────────

IRON_GENES <- c("NCOA4", "FTH1", "TFRC", "SLC7A11", "GPX4", "IREB2",
                "FTL", "HAMP", "SLC40A1", "HMOX1", "CYBRD1", "STEAP3")

AUTOPHAGY_GENES <- c("BECN1", "MAP1LC3B", "ATG5", "ATG7", "ATG12",
                     "SQSTM1", "ULK1", "WIPI2", "LAMP2", "RAB7A")

MAPK_AKT_GENES <- c("BRAF", "KRAS", "NRAS", "MAP2K1", "MAPK1", "MAPK3",
                    "AKT1", "AKT2", "PIK3CA", "PIK3R1", "PTEN",
                    "DUSP6", "DUSP4", "RAF1")

MDR_GENES <- c("ABCB1", "YBX1")

# Melanoma cancer stem cell markers (not mesenchymal stem cell)
STEM_GENES <- c("SOX2", "POU5F1", "NANOG", "KLF4", "MYC",
                "CD44", "PROM1", "ALDH1A1", "ABCB5", "JARID1B",
                "ZEB1", "TWIST1", "SNAI1", "SNAI2", "VIM",
                "NGFR", "EGFR", "AXL", "NES", "SOX10")

# ── get_deg_df helper (needed before deg_mrna_SC_vs_RC) ───────────────────────

get_deg_df <- function(comp_obj, q = 0.8) {
  df <- as.data.frame(comp_obj$degs)
  df$mirna <- rownames(df)                          # column named "mirna" for both
  colnames(df)[colnames(df) == "M"]      <- "log2FC"
  colnames(df)[colnames(df) == "A_mean"] <- "mean_A"
  colnames(df)[colnames(df) == "B_mean"] <- "mean_B"
  df |> dplyr::filter(prob >= q) |>
    dplyr::mutate(direction = ifelse(log2FC > 0, "Up_in_B", "Down_in_B"))
}

# ── Primary resistance DEGs (SC vs RC) ────────────────────────────────────────

deg_mrna_SC_vs_RC <- get_deg_df(comp1)
message(sprintf("DEGs SC vs RC (q>=0.8): %d", nrow(deg_mrna_SC_vs_RC)))

# ── Now ALL_KEY_GENES is safe to run ─────────────────────────────────────────

ALL_KEY_GENES <- unique(c(
  IRON_GENES, AUTOPHAGY_GENES, MAPK_AKT_GENES, MDR_GENES, STEM_GENES,
  deg_mrna_SC_vs_RC$mirna[deg_mrna_SC_vs_RC$prob >= 0.9]
))

message(sprintf("ALL_KEY_GENES: %d genes", length(ALL_KEY_GENES)))
ALL_KEY_GENES <- unique(c(
  IRON_GENES, AUTOPHAGY_GENES, MAPK_AKT_GENES, MDR_GENES, STEM_GENES,
  # Add any additional genes that appeared significant in your NOISeq results
  deg_mrna_SC_vs_RC$gene[deg_mrna_SC_vs_RC$prob >= 0.9]  # highest-confidence DEGs
))

# ── TIER 1: Visualization network ────────────────────────────────────────────
# Rule: validated interactions regardless of target
#     + predicted interactions only if target is in key gene sets
#     This is what gets plotted in R

vis_network <- full_network |>
  dplyr::filter(
    type == "validated" |
      gene %in% ALL_KEY_GENES
  ) |>
  dplyr::mutate(
    gene_category = dplyr::case_when(
      gene %in% IRON_GENES      ~ "Iron/Ferroptosis",
      gene %in% AUTOPHAGY_GENES ~ "Autophagy",
      gene %in% MAPK_AKT_GENES  ~ "MAPK/AKT",
      gene %in% MDR_GENES       ~ "MDR",
      gene %in% STEM_GENES      ~ "Stemness",
      TRUE                      ~ "Other DEG"
    )
  )

message(sprintf("\n  Visualization network: %d pairs", nrow(vis_network)))
message(sprintf("  Breakdown by category:"))
print(table(vis_network$gene_category))

# ── TIER 2: Full network → Cytoscape only ────────────────────────────────────
# Keep everything from full_network; don't try to plot this in R

message(sprintf("\n  Full network (Cytoscape): %d pairs", nrow(full_network)))

# Annotate full network gene categories too
full_network <- full_network |>
  dplyr::mutate(
    gene_category = dplyr::case_when(
      gene %in% IRON_GENES      ~ "Iron/Ferroptosis",
      gene %in% AUTOPHAGY_GENES ~ "Autophagy",
      gene %in% MAPK_AKT_GENES  ~ "MAPK/AKT",
      gene %in% MDR_GENES       ~ "MDR",
      gene %in% STEM_GENES      ~ "Stemness",
      TRUE                      ~ "Other DEG"
    )
  )


# ── P1.5  Build igraph objects ────────────────────────────────────────────────

build_graph <- function(network_df, mirna_fc_df, mrna_fc_df) {
  # Node table
  mirna_nodes <- data.frame(
    name      = unique(network_df$mirna),
    node_type = "miRNA",
    log2FC    = mirna_fc_df$log2FC_mirna[
      match(unique(network_df$mirna), mirna_fc_df$mirna)],
    category  = "miRNA",
    stringsAsFactors = FALSE
  )
  gene_nodes <- data.frame(
    name      = unique(network_df$gene),
    node_type = "mRNA",
    log2FC    = mrna_fc_df$log2FC_mrna[
      match(unique(network_df$gene), mrna_fc_df$gene)],
    category  = network_df$gene_category[
      match(unique(network_df$gene), network_df$gene)],
    stringsAsFactors = FALSE
  )
  all_nodes <- dplyr::bind_rows(mirna_nodes, gene_nodes)
  
  igraph::graph_from_data_frame(
    d        = network_df[, c("mirna","gene","type",
                              "log2FC_mirna","log2FC_mrna","gene_category")],
    directed = TRUE,
    vertices = all_nodes
  )
}

g_vis  <- build_graph(vis_network,  mirna_fc, mrna_fc)
g_full <- build_graph(full_network, mirna_fc, mrna_fc)

message(sprintf("\n  Visualization graph: %d nodes, %d edges",
                igraph::vcount(g_vis), igraph::ecount(g_vis)))
message(sprintf("  Full graph:          %d nodes, %d edges",
                igraph::vcount(g_full), igraph::ecount(g_full)))


# ── P1.6  Network visualization (vis_network only) ────────────────────────────
all_de
NODE_COLORS <- c(
  "miRNA"            = "#7B2D8B",
  "Iron/Ferroptosis" = "#E64B35",
  "Autophagy"        = "#4DBBD5",
  "MAPK/AKT"         = "#00A087",
  "MDR"              = "#F39B7F",
  "Stemness"         = "#3C5488",
  "Other DEG"        = "#8491B4"
)

# Label: all miRNAs + iron/autophagy genes + top-degree mRNA nodes
top_degree_genes <- names(sort(igraph::degree(g_vis), decreasing = TRUE))[1:15]
label_nodes <- unique(c(
  unique(vis_network$mirna),
  IRON_GENES, AUTOPHAGY_GENES,
  top_degree_genes
))

tg_vis <- tidygraph::as_tbl_graph(g_vis) |>
  tidygraph::activate(nodes) |>
  dplyr::mutate(
    label     = ifelse(name %in% label_nodes, name, NA_character_),
    degree    = igraph::degree(g_vis)
  )

p_network <- ggraph::ggraph(tg_vis, layout = "stress") +
  ggraph::geom_edge_arc(
    aes(color = log2FC_mirna > 0),
    arrow       = arrow(length = unit(1.5, "mm"), type = "closed"),
    end_cap     = circle(2.5, "mm"),
    alpha       = 0.4,
    linewidth   = 0.4,
    strength    = 0.1
  ) +
  ggraph::geom_node_point(
    aes(size  = abs(log2FC),
        color = category,
        shape = node_type),
    alpha = 0.9
  ) +
  ggraph::geom_node_text(
    aes(label = label),
    size      = 2.6,
    repel     = TRUE,
    fontface  = "italic",
    max.overlaps = 25
  ) +
  scale_color_manual(values = NODE_COLORS, name = "Category") +
  scale_shape_manual(
    values = c(miRNA = 18, mRNA = 16),
    name   = "Node type"
  ) +
  scale_size_continuous(range = c(2, 8), name = "|log2FC|") +
  ggraph::scale_edge_color_manual(
    values = c(`TRUE`  = "#B2182B",
               `FALSE` = "#2166AC"),
    labels = c("miRNA up in RC", "miRNA down in RC"),
    name   = "miRNA direction"
  ) +
  labs(
    title    = "miRNA-mRNA Integration Network (Visualization tier)",
    subtitle = sprintf("%d miRNAs -> %d target genes | validated + key pathway predicted",
                       sum(igraph::V(g_vis)$node_type == "miRNA"),
                       sum(igraph::V(g_vis)$node_type == "mRNA")),
    caption  = "Full network exported to Cytoscape (.graphml)"
  ) +
  theme_graph(base_family = "sans") +
  theme(
    plot.title    = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 9),
    legend.position = "right"
  )

ggsave(file.path(OUTPUT_DIR, "network_mirna_mrna.pdf"),
       p_network, width = 14, height = 11, device = cairo_pdf)
message("  Network plot saved.")


# ── P1.7  Export both graphs for Cytoscape ────────────────────────────────────

igraph::write_graph(g_vis,
                    file   = file.path(OUTPUT_DIR, "network_vis_tier.graphml"),
                    format = "graphml")

igraph::write_graph(g_full,
                    file   = file.path(OUTPUT_DIR, "network_full_tier.graphml"),
                    format = "graphml")

message("  Both networks exported as .graphml for Cytoscape.")


# ── P1.8  Node statistics (run on full network for complete picture) ──────────

node_stats <- data.frame(
  gene        = igraph::V(g_full)$name,
  node_type   = igraph::V(g_full)$node_type,
  category    = igraph::V(g_full)$category,
  log2FC      = igraph::V(g_full)$log2FC,
  degree      = igraph::degree(g_full),
  betweenness = igraph::betweenness(g_full, normalized = TRUE)
) |>
  dplyr::arrange(desc(betweenness))

write.csv(node_stats,
          file.path(OUTPUT_DIR, "network_node_statistics.csv"),
          row.names = FALSE)

# Report hub nodes — high betweenness means they bridge miRNA and mRNA layers
top_hubs <- node_stats |>
  dplyr::filter(node_type == "mRNA") |>
  dplyr::slice_head(n = 10)

message("\n  Top hub target genes (by betweenness in full network):")
print(top_hubs[, c("gene","category","log2FC","degree","betweenness")])

BiocManager::install("RCy3")