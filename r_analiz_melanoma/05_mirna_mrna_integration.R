# =============================================================================
# 05_mirna_mrna_integration.R
# INTEGRATIVE ANALYSIS: ENCORAFENIB RESISTANCE IN MALIGNANT MELANOMA
#
# Input:   checkpoint_03.RData + pathway_gene_df artifact
# Output:  method_b network tables/plots, checkpoint_05.RData
#
# This stage is independent of TF and WGCNA outputs.
# =============================================================================


# ── Configuration ─────────────────────────────────────────────────────────────

DATA_DIR   <- "/home/ramin/mela/data/ilk"  # <-- UPDATE THIS PATH
OUTPUT_DIR <- file.path(DATA_DIR, "newoutput")
METHOD_B_DIR <- file.path(OUTPUT_DIR, "method_b")
dir.create(METHOD_B_DIR, showWarnings = FALSE, recursive = TRUE)


# ── Packages ──────────────────────────────────────────────────────────────────

cran_pkgs <- c(
  "dplyr", "tidyr", "ggplot2", "ggraph", "igraph", "tidygraph",
  "ggrepel", "pheatmap", "RColorBrewer", "multiMiR", "tibble"
)
for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggraph)
  library(igraph)
  library(tidygraph)
  library(ggrepel)
  library(pheatmap)
  library(RColorBrewer)
  library(multiMiR)
  library(tibble)
})


# ── Load checkpoint ───────────────────────────────────────────────────────────

load(file.path(DATA_DIR, "checkpoint_03.RData"))


# ── Load pathway_gene_df ─────────────────────────────────────────────────────

resolve_pathway_gene_file <- function() {
  candidates <- c(
    file.path(DATA_DIR, "newoutput", "pathway_gene_assignments2.csv"),
    file.path(DATA_DIR, "newoutput", "pathway_gene_assignments.csv2"),
    file.path(DATA_DIR, "output", "pathway_gene_assignments2.csv"),
    file.path(DATA_DIR, "output", "pathway_gene_assignments.csv2"),
    file.path(OUTPUT_DIR, "pathway_gene_assignments2.csv"),
    file.path(OUTPUT_DIR, "pathway_gene_assignments.csv2")
  )

  existing <- candidates[file.exists(candidates)]
  if (length(existing) > 0) {
    return(existing[1])
  }
  NULL
}

if (!exists("pathway_gene_df")) {
  pathway_file <- resolve_pathway_gene_file()
  if (is.null(pathway_file)) {
    stop("pathway_gene_df not found. Provide the object in-session or place a pathway assignment CSV in newoutput/ or output/.")
  }
  pathway_gene_df <- read.csv(pathway_file, stringsAsFactors = FALSE, check.names = FALSE)
  message("  Reloaded pathway_gene_df from: ", pathway_file)
}

if ("gene_symbol" %in% colnames(pathway_gene_df) && !"gene" %in% colnames(pathway_gene_df)) {
  pathway_gene_df <- pathway_gene_df |> dplyr::rename(gene = gene_symbol)
}

required_pathway_cols <- c("gene", "pathway")
missing_cols <- setdiff(required_pathway_cols, colnames(pathway_gene_df))
if (length(missing_cols) > 0) {
  stop("pathway_gene_df is missing required columns: ", paste(missing_cols, collapse = ", "))
}


# ── Helpers ───────────────────────────────────────────────────────────────────

save_graph_if_possible <- function(graph_obj, path) {
  tryCatch(
    igraph::write_graph(graph_obj, path, format = "graphml"),
    error = function(e) message("  Skipped graph export (", basename(path), "): ", conditionMessage(e))
  )
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

load_or_compute_mrna_fc <- function() {
  if (exists("comp1") && !is.null(comp1$full)) {
    df <- data.frame(
      gene = rownames(comp1$full),
      log2FC = comp1$full$M,
      stringsAsFactors = FALSE
    )
    message("  mRNA FC source: comp1$full")
    return(df)
  }

  df <- data.frame(
    gene = rownames(mrna_norm),
    log2FC = log2((mrna_norm[, "A375_RC"] + 1) / (mrna_norm[, "A375_SC"] + 1)),
    stringsAsFactors = FALSE
  )
  message("  mRNA FC source: computed from mrna_norm")
  df
}


# =============================================================================
# METHOD B: reverse multiMiR query using pathway genes as targets
# =============================================================================

message("\n=== PRIORITY 1: miRNA-mRNA integration network ===")

deg_mirna_df <- as.data.frame(mirna_comp1$degs)
deg_mirna_df$mirna <- rownames(deg_mirna_df)
if ("M" %in% colnames(deg_mirna_df)) {
  deg_mirna_df <- deg_mirna_df |> dplyr::rename(log2FC = M)
}
if ("A_mean" %in% colnames(deg_mirna_df)) {
  deg_mirna_df <- deg_mirna_df |> dplyr::rename(mean_SC = A_mean)
}
if ("B_mean" %in% colnames(deg_mirna_df)) {
  deg_mirna_df <- deg_mirna_df |> dplyr::rename(mean_RC = B_mean)
}
deg_mirna_ids <- unique(deg_mirna_df$mirna)
message(sprintf("  DEG miRNAs (SC vs RC): %d", length(deg_mirna_ids)))

mirna_fc_b <- data.frame(
  mirna = rownames(mirna_comp1$full),
  log2FC = mirna_comp1$full$M,
  stringsAsFactors = FALSE
) |>
  dplyr::filter(mirna %in% deg_mirna_ids)

mrna_fc_b <- load_or_compute_mrna_fc()

all_pathway_genes <- unique(pathway_gene_df$gene)
message(sprintf("  Unique pathway genes: %d", length(all_pathway_genes)))

val_result <- tryCatch(
  multiMiR::get_multimir(
    target = all_pathway_genes,
    table = "validated",
    summary = FALSE
  ),
  error = function(e) {
    message("  Validated query failed: ", conditionMessage(e))
    NULL
  }
)

pred_result <- tryCatch(
  multiMiR::get_multimir(
    target = all_pathway_genes,
    table = "predicted",
    predicted.cutoff = 30,
    predicted.cutoff.type = "p",
    predicted.site = "all",
    summary = FALSE
  ),
  error = function(e) {
    message("  Predicted query failed: ", conditionMessage(e))
    NULL
  }
)

raw_b <- dplyr::bind_rows(
  if (!is.null(val_result))  val_result@data  |> dplyr::mutate(type = "validated"),
  if (!is.null(pred_result)) pred_result@data |> dplyr::mutate(type = "predicted")
)
message(sprintf("  Raw multiMiR interactions: %d", nrow(raw_b)))

target_raw <- raw_b |>
  dplyr::select(
    mirna = mature_mirna_id,
    gene = target_symbol,
    type
  ) |>
  dplyr::filter(
    !is.na(mirna), nzchar(mirna),
    !is.na(gene), nzchar(gene)
  ) |>
  dplyr::group_by(mirna, gene) |>
  dplyr::summarise(
    type = ifelse("validated" %in% type, "validated", "predicted"),
    .groups = "drop"
  )

target_deg <- target_raw |>
  dplyr::filter(mirna %in% deg_mirna_ids)

target_anticor <- target_deg |>
  dplyr::inner_join(mirna_fc_b |> dplyr::select(mirna, log2FC_mirna = log2FC), by = "mirna") |>
  dplyr::inner_join(mrna_fc_b  |> dplyr::select(gene, log2FC_mrna = log2FC), by = "gene") |>
  dplyr::filter(sign(log2FC_mirna) != sign(log2FC_mrna)) |>
  dplyr::distinct(mirna, gene, .keep_all = TRUE)

message(sprintf("  Anti-correlated miRNA-mRNA pairs: %d", nrow(target_anticor)))

pw_lookup <- pathway_gene_df |>
  dplyr::select(gene, pathway, dplyr::any_of(c("NES", "source")))

method_b_network <- target_anticor |>
  dplyr::inner_join(pw_lookup, by = "gene") |>
  dplyr::rename(gene_pathway = pathway)

if (nrow(method_b_network) == 0) {
  stop("No pathway-annotated miRNA-mRNA pairs were found.")
}

if (!"pathway_NES" %in% colnames(method_b_network)) {
  if ("NES" %in% colnames(method_b_network)) {
    method_b_network$pathway_NES <- method_b_network$NES
  } else {
    method_b_network$pathway_NES <- NA_real_
  }
}

write.csv(
  method_b_network,
  file.path(METHOD_B_DIR, "method_b_full_network.csv"),
  row.names = FALSE
)


# =============================================================================
# Per-pathway subnetworks
# =============================================================================

pathway_colors <- c(
  "DNA Damage Response"     = "#E64B35",
  "Cell Cycle / G2M"        = "#F39B7F",
  "Homologous Recombination" = "#D62728",
  "OXPHOS / Complex I"      = "#9467BD",
  "Iron Transport"         = "#8C564B",
  "Interferon Signaling"   = "#4DBBD5",
  "Lysosome / Autophagy"    = "#00A087",
  "Cytokine Signaling"     = "#AEC7E8",
  "Ferroptosis"            = "#C5B0D5",
  "EMT / Mesenchymal"      = "#2CA02C",
  "p53 Signaling"          = "#FF7F0E",
  "PI3K / AKT / mTOR"      = "#1F77B4",
  "Hypoxia / HIF1A"        = "#BCBD22",
  "WNT Signaling"          = "#17BECF",
  "Apoptosis"              = "#7F7F7F",
  "miRNA"                  = "#7B2D8B"
)

paths_b <- names(sort(table(method_b_network$gene_pathway), decreasing = TRUE))
paths_b <- paths_b[table(method_b_network$gene_pathway)[paths_b] >= 3]
message(sprintf("  Pathways with >=3 edges: %d", length(paths_b)))

for (pw_label in paths_b) {
  pw_edges <- method_b_network |> dplyr::filter(gene_pathway == pw_label)

  pw_mirna_nodes <- data.frame(
    name = unique(pw_edges$mirna),
    node_type = "miRNA",
    category = "miRNA",
    log2FC = mirna_fc_b$log2FC[match(unique(pw_edges$mirna), mirna_fc_b$mirna)],
    stringsAsFactors = FALSE
  )

  pw_gene_nodes <- pw_edges |>
    dplyr::distinct(gene, log2FC_mrna) |>
    dplyr::rename(name = gene, log2FC = log2FC_mrna) |>
    dplyr::mutate(node_type = "mRNA", category = pw_label)

  pw_nodes <- dplyr::bind_rows(pw_mirna_nodes, pw_gene_nodes)
  g_pw <- igraph::graph_from_data_frame(
    d = pw_edges[, c("mirna", "gene", "type", "log2FC_mirna", "log2FC_mrna")],
    directed = TRUE,
    vertices = pw_nodes
  )

  n_nodes <- igraph::vcount(g_pw)
  label_set <- if (n_nodes <= 35) {
    pw_nodes$name
  } else {
    top_genes <- names(sort(igraph::degree(g_pw), decreasing = TRUE))[1:min(20, n_nodes)]
    unique(c(unique(pw_edges$mirna), top_genes))
  }

  tg_pw <- tidygraph::as_tbl_graph(g_pw) |>
    tidygraph::activate(nodes) |>
    dplyr::mutate(
      label = ifelse(name %in% label_set, name, NA_character_),
      degree = igraph::degree(g_pw)
    )

  pw_color <- pathway_colors[[pw_label]]
  if (is.na(pw_color) || is.null(pw_color)) {
    pw_color <- "#8491B4"
  }

  p_pw <- ggraph::ggraph(tg_pw, layout = if (n_nodes <= 25) "stress" else "fr") +
    ggraph::geom_edge_arc(
      aes(color = log2FC_mirna > 0, linetype = type),
      arrow = arrow(length = grid::unit(1.5, "mm"), type = "closed"),
      end_cap = ggraph::circle(2.5, "mm"),
      alpha = 0.55,
      linewidth = 0.45,
      strength = 0.15
    ) +
    ggraph::geom_node_point(
      aes(size = degree + 1, shape = node_type, color = ifelse(node_type == "miRNA", "miRNA", pw_label)),
      alpha = 0.92
    ) +
    ggraph::geom_node_text(
      aes(label = label),
      size = 2.7,
      repel = TRUE,
      fontface = "italic",
      max.overlaps = 30,
      segment.size = 0.2,
      segment.color = "grey50"
    ) +
    scale_color_manual(
      values = setNames(c("#7B2D8B", pw_color), c("miRNA", pw_label)),
      name = "Node type"
    ) +
    scale_shape_manual(values = c(miRNA = 18, mRNA = 16), guide = "none") +
    scale_size_continuous(range = c(2.5, 9), name = "Degree") +
    ggraph::scale_edge_color_manual(
      values = c(`TRUE` = "#B2182B", `FALSE` = "#2166AC"),
      labels = c("miRNA up in RC", "miRNA down in RC"),
      name = "miRNA FC"
    ) +
    ggraph::scale_edge_linetype_manual(
      values = c(validated = "solid", predicted = "dashed"),
      name = "Evidence"
    ) +
    labs(
      title = sprintf("[Method B] %s - Reverse Query Subnetwork", pw_label),
      subtitle = sprintf(
        "%d miRNAs -> %d genes | %d edges (%d validated, %d predicted)",
        sum(igraph::V(g_pw)$node_type == "miRNA"),
        sum(igraph::V(g_pw)$node_type == "mRNA"),
        igraph::ecount(g_pw),
        sum(pw_edges$type == "validated"),
        sum(pw_edges$type == "predicted")
      )
    ) +
    ggraph::theme_graph(base_family = "sans") +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 8),
      legend.position = "bottom",
      legend.text = element_text(size = 7)
    )

  fname <- gsub("[/ ]", "_", pw_label)
  ggsave(
    file.path(METHOD_B_DIR, sprintf("methodB_subnetwork_%s.pdf", fname)),
    p_pw, width = 10, height = 8
  )
  save_graph_if_possible(
    g_pw,
    file.path(METHOD_B_DIR, sprintf("methodB_subnetwork_%s.graphml", fname))
  )
}


# =============================================================================
# Combined outputs
# =============================================================================

all_mirna_nodes_b <- data.frame(
  name = unique(method_b_network$mirna),
  node_type = "miRNA",
  category = "miRNA",
  log2FC = mirna_fc_b$log2FC[match(unique(method_b_network$mirna), mirna_fc_b$mirna)],
  stringsAsFactors = FALSE
)

all_gene_nodes_b <- method_b_network |>
  dplyr::distinct(gene, gene_pathway, log2FC_mrna) |>
  dplyr::rename(name = gene, category = gene_pathway, log2FC = log2FC_mrna) |>
  dplyr::mutate(node_type = "mRNA")

all_nodes_b <- dplyr::bind_rows(all_mirna_nodes_b, all_gene_nodes_b)
g_combined_b <- igraph::graph_from_data_frame(
  d = method_b_network[, c("mirna", "gene", "type", "gene_pathway", "log2FC_mirna", "log2FC_mrna")],
  directed = TRUE,
  vertices = all_nodes_b
)
save_graph_if_possible(g_combined_b, file.path(METHOD_B_DIR, "methodB_combined.graphml"))

mirna_pw_map_b <- method_b_network |>
  dplyr::group_by(mirna, gene_pathway) |>
  dplyr::summarise(
    n_targets = n(),
    n_validated = sum(type == "validated"),
    target_genes = paste(sort(gene), collapse = "; "),
    mirna_log2FC = unique(log2FC_mirna),
    mirna_direction = ifelse(unique(log2FC_mirna) > 0, "Up_in_RC", "Down_in_RC"),
    .groups = "drop"
  ) |>
  dplyr::arrange(gene_pathway, desc(n_targets))

write.csv(
  mirna_pw_map_b,
  file.path(METHOD_B_DIR, "methodB_mirna_pathway_map.csv"),
  row.names = FALSE
)

map_wide_b <- mirna_pw_map_b |>
  dplyr::select(mirna, gene_pathway, n_targets) |>
  tidyr::pivot_wider(names_from = gene_pathway, values_from = n_targets, values_fill = 0L) |>
  as.data.frame()
rownames(map_wide_b) <- map_wide_b$mirna
map_wide_b$mirna <- NULL
map_mat_b <- as.matrix(map_wide_b)
map_mat_b <- map_mat_b[rowSums(map_mat_b) > 0, colSums(map_mat_b) > 0, drop = FALSE]

if (nrow(map_mat_b) >= 2 && ncol(map_mat_b) >= 2) {
  dir_ann_b <- mirna_pw_map_b |>
    dplyr::distinct(mirna, mirna_direction) |>
    dplyr::filter(mirna %in% rownames(map_mat_b)) |>
    tibble::column_to_rownames("mirna")

  pheatmap::pheatmap(
    map_mat_b,
    color = colorRampPalette(c("white", "#FDD0A2", "#E64B35"))(50),
    annotation_row = dir_ann_b,
    annotation_colors = list(
      mirna_direction = c(Up_in_RC = "#B2182B", Down_in_RC = "#2166AC")
    ),
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    fontsize_row = 7,
    fontsize_col = 8,
    display_numbers = TRUE,
    number_format = "%d",
    number_color = "black",
    main = "Method B: miRNA -> Pathway Regulatory Map",
    filename = file.path(METHOD_B_DIR, "methodB_heatmap_mirna_pathway.pdf"),
    width = 10,
    height = max(5, nrow(map_mat_b) * 0.28 + 3)
  )
}

priority_table_b <- method_b_network |>
  dplyr::group_by(mirna, gene) |>
  dplyr::summarise(
    type = unique(type),
    pathways = paste(sort(unique(gene_pathway)), collapse = "; "),
    n_pathways = n_distinct(gene_pathway),
    log2FC_mirna = unique(log2FC_mirna),
    log2FC_mrna = unique(log2FC_mrna),
    mirna_direction = ifelse(unique(log2FC_mirna) > 0, "Up_in_RC", "Down_in_RC"),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    evidence_score = ifelse(type == "validated", 2, 1),
    iron_gene = gene %in% c("FTH1", "NCOA4", "TFRC", "SLC7A11", "GPX4",
                            "IREB2", "FTL", "HMOX1", "STEAP3", "CYBRD1"),
    priority_score = evidence_score + n_pathways + (iron_gene * 3)
  ) |>
  dplyr::arrange(dplyr::desc(priority_score), dplyr::desc(abs(log2FC_mirna)))

write.csv(
  priority_table_b,
  file.path(METHOD_B_DIR, "methodB_priority_wetlab.csv"),
  row.names = FALSE
)


# =============================================================================
# CHECKPOINT SAVE
# =============================================================================

save(
  pathway_gene_df,
  deg_mirna_df, deg_mirna_ids,
  mirna_fc_b, mrna_fc_b,
  method_b_network, mirna_pw_map_b, priority_table_b,
  file = file.path(DATA_DIR, "checkpoint_05.RData")
)

message("\ncheckpoint_05.RData saved.")
