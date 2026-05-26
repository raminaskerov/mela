#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(STRINGdb)
  library(igraph)
})

ROOT_DIR <- "C:/bioinformatics"
PPI_OUTPUT_DIR <- file.path(ROOT_DIR, "output", "ppi")
CACHE_DIR <- file.path(ROOT_DIR, "cash", "stringdb_cache")
dir.create(PPI_OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(CACHE_DIR, recursive = TRUE, showWarnings = FALSE)

session_candidates <- c(
  file.path(ROOT_DIR, "melanoma_session.RData"),
  file.path(ROOT_DIR, "r_analiz_melanoma", "melanoma_session.RData"),
  file.path(ROOT_DIR, "r_analiz_melanoma", "finished(supheli).RData"),
  file.path(ROOT_DIR, "r_analiz_melanoma", "fresh.RData"),
  file.path(ROOT_DIR, "r_analiz_melanoma", "wgcna.RData")
)

session_file <- session_candidates[file.exists(session_candidates)][1]
if (is.na(session_file)) {
  stop("No melanoma session RData file found. Expected melanoma_session.RData or a compatible saved session.")
}

loaded_objects <- load(session_file)
if (!exists("comp1")) {
  stop("Loaded session does not contain comp1: ", session_file)
}
if (is.null(comp1$degs) || is.null(comp1$full) || !"M" %in% colnames(as.data.frame(comp1$full))) {
  stop("comp1 must contain degs and full$M.")
}

deg_genes <- unique(rownames(as.data.frame(comp1$degs)))
deg_genes <- deg_genes[!is.na(deg_genes) & nzchar(deg_genes)]
if (length(deg_genes) < 2) {
  stop("Fewer than two DEG gene symbols found in comp1$degs.")
}

full_df <- as.data.frame(comp1$full)
log2fc_lookup <- setNames(1 * as.numeric(full_df$M), rownames(full_df))

string_db <- STRINGdb$new(
  version = "11.5",
  species = 9606,
  score_threshold = 400,
  input_directory = CACHE_DIR
)

deg_df <- data.frame(gene = deg_genes, stringsAsFactors = FALSE)
mapped <- string_db$map(deg_df, "gene", removeUnmappedRows = TRUE)
mapped <- mapped[!duplicated(mapped$STRING_id), c("gene", "STRING_id")]
if (nrow(mapped) < 2) {
  stop("Fewer than two DEG genes mapped to STRING.")
}

interactions <- string_db$get_interactions(mapped$STRING_id)
if (is.null(interactions) || nrow(interactions) == 0) {
  stop("STRING returned no interactions for mapped SC vs RC DEGs.")
}

id_to_gene <- setNames(mapped$gene, mapped$STRING_id)
edge_df <- data.frame(
  from = id_to_gene[interactions$from],
  to = id_to_gene[interactions$to],
  combined_score = interactions$combined_score,
  stringsAsFactors = FALSE
)
edge_df <- edge_df[!is.na(edge_df$from) & !is.na(edge_df$to) & edge_df$from != edge_df$to, ]
edge_df <- unique(edge_df)
if (nrow(edge_df) == 0) {
  stop("No intra-DEG STRING interactions remained after mapping STRING IDs back to gene symbols.")
}

network_genes <- sort(unique(c(edge_df$from, edge_df$to)))
g <- graph_from_data_frame(edge_df, directed = FALSE, vertices = data.frame(name = network_genes))
g <- simplify(g, remove.multiple = TRUE, remove.loops = TRUE, edge.attr.comb = "max")

V(g)$gene_symbol <- V(g)$name
V(g)$log2FC <- unname(log2fc_lookup[V(g)$name])
V(g)$degree_centrality <- degree(g, mode = "all", normalized = FALSE)

degree_df <- data.frame(
  gene = V(g)$gene_symbol,
  ppi_degree = V(g)$degree_centrality,
  log2FC = V(g)$log2FC,
  stringsAsFactors = FALSE
)
degree_df <- degree_df[order(-degree_df$ppi_degree, degree_df$gene), ]
hub_n <- max(1L, ceiling(0.10 * nrow(degree_df)))
hub_df <- degree_df[seq_len(hub_n), ]

get_gene_column <- function(df) {
  lower_names <- tolower(colnames(df))
  gene_col <- which(lower_names %in% c("gene", "symbol", "genes", "tf"))[1]
  if (is.na(gene_col)) {
    return(character())
  }
  unique(as.character(df[[gene_col]]))
}

wgcna_hubs <- character()
if (exists("wgcna_hub_genes")) {
  obj <- get("wgcna_hub_genes")
  if (!is.null(rownames(obj))) {
    wgcna_hubs <- rownames(obj)
  }
  if (is.data.frame(obj)) {
    wgcna_hubs <- unique(c(wgcna_hubs, get_gene_column(obj)))
  }
}
wgcna_csv <- file.path(ROOT_DIR, "output", "wgcna", "wgcna_hub_genes.csv")
if (length(wgcna_hubs) == 0 && file.exists(wgcna_csv)) {
  wgcna_hubs <- get_gene_column(read.csv(wgcna_csv, stringsAsFactors = FALSE, check.names = FALSE))
}
wgcna_hubs <- unique(wgcna_hubs[!is.na(wgcna_hubs) & nzchar(wgcna_hubs)])

tf_candidates <- character()
tf_csv_candidates <- c(
  file.path(ROOT_DIR, "output", "TF_converged_final_candidates.csv"),
  file.path(ROOT_DIR, "output", "tf_enrich", "TF_converged_final_candidates.csv")
)
tf_csv <- tf_csv_candidates[file.exists(tf_csv_candidates)][1]
if (!is.na(tf_csv)) {
  tf_candidates <- get_gene_column(read.csv(tf_csv, stringsAsFactors = FALSE, check.names = FALSE))
}
tf_candidates <- unique(tf_candidates[!is.na(tf_candidates) & nzchar(tf_candidates)])

iron_genes <- c("NCOA4", "FTH1", "TFRC", "SLC7A11", "GPX4", "IREB2", "FTL", "HMOX1")
hub_summary <- data.frame(
  gene = hub_df$gene,
  ppi_degree = hub_df$ppi_degree,
  in_wgcna_hub = hub_df$gene %in% wgcna_hubs,
  in_tf_candidates = hub_df$gene %in% tf_candidates,
  log2FC = hub_df$log2FC,
  iron_related = hub_df$gene %in% iron_genes,
  stringsAsFactors = FALSE
)

write_graph(g, file.path(PPI_OUTPUT_DIR, "ppi_SC_vs_RC.graphml"), format = "graphml")
write.csv(hub_summary, file.path(PPI_OUTPUT_DIR, "ppi_hub_convergence.csv"), row.names = FALSE)

stats_lines <- c(
  "PPI network: SC vs RC DEGs",
  paste0("session_file: ", "finished(supheli).RData"),
  paste0("requested_session_found: ", file.exists(file.path(ROOT_DIR, "melanoma_session.RData")) ||
           file.exists(file.path(ROOT_DIR, "r_analiz_melanoma", "melanoma_session.RData"))),
  paste0("stringdb_version: 11.5"),
  paste0("species: 9606"),
  paste0("score_threshold: 400"),
  paste0("deg_genes_input: ", length(deg_genes)),
  paste0("mapped_deg_genes: ", nrow(mapped)),
  paste0("total_nodes: ", vcount(g)),
  paste0("total_edges: ", ecount(g)),
  paste0("hub_gene_count: ", nrow(hub_summary)),
  paste0("hubs_in_wgcna_hub: ", sum(hub_summary$in_wgcna_hub)),
  paste0("hubs_in_tf_candidates: ", sum(hub_summary$in_tf_candidates)),
  paste0("hubs_in_both_wgcna_and_tf: ", sum(hub_summary$in_wgcna_hub & hub_summary$in_tf_candidates)),
  paste0("hubs_iron_related: ", sum(hub_summary$iron_related))
)
writeLines(stats_lines, file.path(PPI_OUTPUT_DIR, "ppi_network_stats.txt"))

message("Wrote: ", file.path(PPI_OUTPUT_DIR, "ppi_SC_vs_RC.graphml"))
message("Wrote: ", file.path(PPI_OUTPUT_DIR, "ppi_hub_convergence.csv"))
message("Wrote: ", file.path(PPI_OUTPUT_DIR, "ppi_network_stats.txt"))

gse45558 <- getGEO("GSE45558", getGPL = FALSE)[[1]]
meta45558 <- pData(gse45558)

sort(table(df$gene_biotype), decreasing = TRUE)