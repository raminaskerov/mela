#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(igraph)
})

ROOT_DIR <- "C:/bioinformatics"
OUTPUT_DIR <- file.path(ROOT_DIR, "output")
CYTO_DIR <- file.path(OUTPUT_DIR, "tf_cytoscape")

tf_file <- file.path(OUTPUT_DIR, "TF_converged_final_candidates.csv")
pathway_file <- file.path(OUTPUT_DIR, "TF_pathway_convergence.csv")
deg_file <- file.path(OUTPUT_DIR, "deg", "DEGs_SC_vs_RC_.csv")
wgcna_hub_file <- file.path(OUTPUT_DIR, "wgcna", "wgcna_hub_genes.csv")

edge_file <- file.path(CYTO_DIR, "TF_mRNA_edges.csv")
graph_file <- file.path(CYTO_DIR, "TF_mRNA_network.graphml")

stop_if_missing <- function(path) {
  if (!file.exists(path)) {
    stop("Missing required input file: ", path)
  }
}

clean_text <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  trimws(x)
}

as_flag <- function(x) {
  x <- toupper(clean_text(x))
  x %in% c("TRUE", "T", "1", "YES", "Y")
}

split_genes <- function(x) {
  genes <- unlist(strsplit(clean_text(x), ";", fixed = TRUE), use.names = FALSE)
  genes <- trimws(genes)
  unique(genes[!is.na(genes) & nzchar(genes)])
}

require_columns <- function(data, columns, label) {
  missing <- setdiff(columns, colnames(data))
  if (length(missing) > 0) {
    stop(label, " is missing required columns: ", paste(missing, collapse = ", "))
  }
}

collapse_unique <- function(x) {
  x <- unique(clean_text(x))
  x <- x[nzchar(x)]
  if (length(x) == 0) "" else paste(x, collapse = ";")
}

derive_target_log2fc <- function(nes) {
  nes <- as.numeric(nes)
  signs <- unique(sign(nes[!is.na(nes) & nes != 0]))
  if (length(signs) == 1) {
    return(signs)
  }
  0
}

for (path in c(tf_file, pathway_file)) {
  stop_if_missing(path)
}
dir.create(CYTO_DIR, recursive = TRUE, showWarnings = FALSE)

tf_candidates <- read.csv(tf_file, stringsAsFactors = FALSE, check.names = FALSE)
tf_pathways <- read.csv(pathway_file, stringsAsFactors = FALSE, check.names = FALSE)
deg_df <- NULL
if (file.exists(deg_file)) {
  deg_df <- read.csv(deg_file, stringsAsFactors = FALSE, check.names = FALSE)
}

require_columns(
  tf_candidates,
  c(
    "TF", "delta_RC_vs_SC", "activity_dir", "priority_score", "iron_related",
    "mapk_akt", "stem_related", "is_DEG", "is_hub"
  ),
  "TF candidate table"
)

require_columns(
  tf_pathways,
  c(
    "TF", "pathway", "pathway_NES", "pathway_padj", "overlap_genes",
    "combined_score", "concordant"
  ),
  "TF pathway convergence table"
)

if (!is.null(deg_df)) {
  require_columns(
    deg_df,
    c("gene", "log2FC"),
    "DEG table"
  )
}

tf_candidates$TF <- clean_text(tf_candidates$TF)
candidate_tfs <- unique(tf_candidates$TF[nzchar(tf_candidates$TF)])

tf_pathways$TF <- clean_text(tf_pathways$TF)
tf_pathways <- tf_pathways[
  tf_pathways$TF %in% candidate_tfs & as_flag(tf_pathways$concordant),
  ,
  drop = FALSE
]

edge_rows <- vector("list", nrow(tf_pathways))
for (i in seq_len(nrow(tf_pathways))) {
  genes <- split_genes(tf_pathways$overlap_genes[i])
  if (length(genes) == 0) {
    next
  }

  edge_rows[[i]] <- data.frame(
    from = tf_pathways$TF[i],
    to = genes,
    TF = tf_pathways$TF[i],
    target_gene = genes,
    pathway = tf_pathways$pathway[i],
    pathway_NES = as.numeric(tf_pathways$pathway_NES[i]),
    pathway_padj = as.numeric(tf_pathways$pathway_padj[i]),
    combined_score = as.numeric(tf_pathways$combined_score[i]),
    concordant = TRUE,
    stringsAsFactors = FALSE
  )
}

edges <- do.call(rbind, edge_rows[!vapply(edge_rows, is.null, logical(1))])
if (is.null(edges) || nrow(edges) == 0) {
  stop("No concordant TF-target edges were found after filtering to candidate TFs.")
}

graph_edges <- data.frame(
  from = paste0("TF:", edges$TF),
  to = paste0("mRNA:", edges$target_gene),
  TF = edges$TF,
  target_gene = edges$target_gene,
  pathway = edges$pathway,
  pathway_NES = edges$pathway_NES,
  pathway_padj = edges$pathway_padj,
  combined_score = edges$combined_score,
  concordant = edges$concordant,
  stringsAsFactors = FALSE
)

tf_lookup <- tf_candidates[match(candidate_tfs, tf_candidates$TF), , drop = FALSE]
tf_pathway_membership <- tapply(tf_pathways$pathway, tf_pathways$TF, collapse_unique)
target_pathway_membership <- tapply(edges$pathway, edges$target_gene, collapse_unique)

edges$TF_pathway_membership <- unname(tf_pathway_membership[edges$TF])
edges$target_pathway_membership <- unname(target_pathway_membership[edges$target_gene])
edges$TF_pathway_membership[is.na(edges$TF_pathway_membership)] <- ""
edges$target_pathway_membership[is.na(edges$target_pathway_membership)] <- ""

target_log2fc_lookup <- numeric()
if (!is.null(deg_df)) {
  deg_df$gene <- clean_text(deg_df$gene)
  deg_df <- deg_df[nzchar(deg_df$gene), , drop = FALSE]
  deg_df <- deg_df[!duplicated(deg_df$gene), , drop = FALSE]
  target_log2fc_lookup <- setNames(as.numeric(deg_df$log2FC), deg_df$gene)
}

hub_module <- setNames(rep("", length(candidate_tfs)), candidate_tfs)
if (file.exists(wgcna_hub_file)) {
  wgcna_hubs <- read.csv(wgcna_hub_file, stringsAsFactors = FALSE, check.names = FALSE)
  if (all(c("gene", "module") %in% colnames(wgcna_hubs))) {
    module_lookup <- tapply(wgcna_hubs$module, wgcna_hubs$gene, collapse_unique)
    matched_modules <- unname(module_lookup[candidate_tfs])
    matched_modules[is.na(matched_modules)] <- ""
    hub_module <- setNames(matched_modules, candidate_tfs)
  }
}

tf_nodes <- data.frame(
  name = paste0("TF:", candidate_tfs),
  label = candidate_tfs,
  gene_symbol = candidate_tfs,
  node_type = "TF",
  pathway_membership = unname(tf_pathway_membership[candidate_tfs]),
  delta_RC_vs_SC = as.numeric(tf_lookup$delta_RC_vs_SC),
  activity_dir = tf_lookup$activity_dir,
  priority_score = as.numeric(tf_lookup$priority_score),
  iron_related = as_flag(tf_lookup$iron_related),
  mapk_akt = as_flag(tf_lookup$mapk_akt),
  stem_related = as_flag(tf_lookup$stem_related),
  is_DEG = as_flag(tf_lookup$is_DEG),
  is_hub = as_flag(tf_lookup$is_hub),
  wgcna_hub_module = unname(hub_module[candidate_tfs]),
  log2FC = NA_real_,
  stringsAsFactors = FALSE
)
tf_nodes$pathway_membership[is.na(tf_nodes$pathway_membership)] <- ""

target_genes <- sort(unique(edges$to))
target_nodes <- data.frame(
  name = paste0("mRNA:", target_genes),
  label = target_genes,
  gene_symbol = target_genes,
  node_type = "mRNA",
  pathway_membership = unname(target_pathway_membership[target_genes]),
  delta_RC_vs_SC = NA_real_,
  activity_dir = "",
  priority_score = NA_real_,
  iron_related = FALSE,
  mapk_akt = FALSE,
  stem_related = FALSE,
  is_DEG = target_genes %in% names(target_log2fc_lookup),
  is_hub = FALSE,
  wgcna_hub_module = "",
  log2FC = as.numeric(target_log2fc_lookup[target_genes]),
  stringsAsFactors = FALSE
)
target_nodes$pathway_membership[is.na(target_nodes$pathway_membership)] <- ""

if (all(is.na(target_nodes$log2FC))) {
  target_log2fc_fallback <- tapply(edges$pathway_NES, edges$to, derive_target_log2fc)
  target_nodes$log2FC <- as.numeric(target_log2fc_fallback[target_genes])
  target_nodes$is_DEG <- FALSE
}

nodes <- rbind(tf_nodes, target_nodes)
nodes <- nodes[!duplicated(nodes$name), , drop = FALSE]

g <- graph_from_data_frame(graph_edges, directed = TRUE, vertices = nodes)

write.csv(edges, edge_file, row.names = FALSE)
write_graph(g, graph_file, format = "graphml")

message("Wrote: ", edge_file)
message("Wrote: ", graph_file)
message("Nodes: ", vcount(g), " (TF: ", sum(V(g)$node_type == "TF"),
        ", mRNA: ", sum(V(g)$node_type == "mRNA"), ")")
message("Edges: ", ecount(g))
