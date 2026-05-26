#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(igraph)
})

ROOT_DIR <- "C:/bioinformatics"
PPI_DIR <- file.path(ROOT_DIR, "output", "ppi")

hub_file <- file.path(PPI_DIR, "ppi_hub_convergence.csv")
graph_file <- file.path(PPI_DIR, "ppi_SC_vs_RC.graphml")
reactome_file <- file.path(ROOT_DIR, "output", "enrichment", "GSEA_Reactome_SC_vs_RC.csv")
hallmark_file <- file.path(ROOT_DIR, "output", "pathway_fgsea_significant.csv")
wgcna_file <- file.path(ROOT_DIR, "output", "wgcna", "wgcna_hub_genes.csv")

stop_if_missing <- function(path) {
  if (!file.exists(path)) {
    stop("Missing required input file: ", path)
  }
}

for (path in c(hub_file, graph_file, reactome_file, hallmark_file, wgcna_file)) {
  stop_if_missing(path)
}

clean_text <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  trimws(x)
}

split_genes <- function(x, sep_pattern) {
  x <- clean_text(x)
  if (!nzchar(x)) {
    return(character())
  }
  genes <- unlist(strsplit(x, sep_pattern, perl = TRUE), use.names = FALSE)
  genes <- trimws(genes)
  genes[!is.na(genes) & nzchar(genes)]
}

collapse_unique <- function(x) {
  x <- unique(clean_text(x))
  x <- x[nzchar(x)]
  if (length(x) == 0) "" else paste(x, collapse = ";")
}

nes_range <- function(x) {
  x <- as.numeric(x)
  x <- x[!is.na(x)]
  if (length(x) == 0) {
    return("")
  }
  if (length(x) == 1 || min(x) == max(x)) {
    return(sprintf("%.3f", x[which.max(abs(x))[1]]))
  }
  sprintf("%.3f..%.3f", min(x), max(x))
}

direction_from_nes <- function(nes) {
  ifelse(as.numeric(nes) >= 0, "Enriched_in_RC", "Enriched_in_SC")
}

direction_from_log2fc <- function(log2fc) {
  ifelse(as.numeric(log2fc) >= 0, "Enriched_in_RC", "Enriched_in_SC")
}

make_long_reactome <- function(path) {
  reactome <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  needed <- c("Description", "NES", "core_enrichment")
  if (!all(needed %in% colnames(reactome))) {
    stop("Reactome file must contain columns: ", paste(needed, collapse = ", "))
  }

  rows <- vector("list", nrow(reactome))
  for (i in seq_len(nrow(reactome))) {
    entrez <- split_genes(reactome$core_enrichment[i], "/")
    if (length(entrez) == 0) {
      next
    }
    symbols <- suppressMessages(
      AnnotationDbi::mapIds(
        org.Hs.eg.db,
        keys = entrez,
        keytype = "ENTREZID",
        column = "SYMBOL",
        multiVals = "first"
      )
    )
    symbols <- unique(unname(symbols))
    symbols <- symbols[!is.na(symbols) & nzchar(symbols)]
    if (length(symbols) == 0) {
      next
    }
    rows[[i]] <- data.frame(
      gene = symbols,
      source = "Reactome",
      pathway = reactome$Description[i],
      NES = as.numeric(reactome$NES[i]),
      direction = direction_from_nes(reactome$NES[i]),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows[!vapply(rows, is.null, logical(1))])
}

make_long_hallmark <- function(path) {
  hallmark <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  needed <- c("pathway", "NES", "leadingEdge")
  if (!all(needed %in% colnames(hallmark))) {
    stop("Hallmark file must contain columns: ", paste(needed, collapse = ", "))
  }

  rows <- vector("list", nrow(hallmark))
  for (i in seq_len(nrow(hallmark))) {
    genes <- split_genes(hallmark$leadingEdge[i], ";")
    if (length(genes) == 0) {
      next
    }
    rows[[i]] <- data.frame(
      gene = unique(genes),
      source = "Hallmark",
      pathway = hallmark$pathway[i],
      NES = as.numeric(hallmark$NES[i]),
      direction = direction_from_nes(hallmark$NES[i]),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows[!vapply(rows, is.null, logical(1))])
}

primary_for_gene <- function(gene, memberships) {
  hit <- memberships[memberships$gene == gene, , drop = FALSE]
  if (nrow(hit) == 0) {
    return(c(primary_pathway = "", primary_direction = ""))
  }
  hit <- hit[order(-abs(hit$NES), hit$source, hit$pathway), , drop = FALSE]
  c(primary_pathway = hit$pathway[1], primary_direction = hit$direction[1])
}

annotate_genes <- function(genes, reactome_long, hallmark_long, memberships) {
  rows <- lapply(genes, function(gene) {
    r_hit <- reactome_long[reactome_long$gene == gene, , drop = FALSE]
    h_hit <- hallmark_long[hallmark_long$gene == gene, , drop = FALSE]
    primary <- primary_for_gene(gene, memberships)

    data.frame(
      gene = gene,
      reactome_pathways = collapse_unique(r_hit$pathway),
      reactome_NES_range = nes_range(r_hit$NES),
      hallmark_pathways = collapse_unique(h_hit$pathway),
      hallmark_NES_range = nes_range(h_hit$NES),
      primary_pathway = primary[["primary_pathway"]],
      primary_direction = primary[["primary_direction"]],
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

hubs <- read.csv(hub_file, stringsAsFactors = FALSE, check.names = FALSE)
reactome_long <- make_long_reactome(reactome_file)
hallmark_long <- make_long_hallmark(hallmark_file)
memberships <- rbind(reactome_long, hallmark_long)

wgcna <- read.csv(wgcna_file, stringsAsFactors = FALSE, check.names = FALSE)
module_lookup <- tapply(wgcna$module, wgcna$gene, function(x) paste(unique(x), collapse = ";"))
hubs$wgcna_hub_module <- unname(module_lookup[hubs$gene])
hubs$wgcna_hub_module[is.na(hubs$wgcna_hub_module)] <- ""

annotated <- merge(hubs, annotate_genes(hubs$gene, reactome_long, hallmark_long, memberships),
                   by = "gene", sort = FALSE)
annotated <- annotated[match(hubs$gene, annotated$gene), ]
hub_annotated_file <- file.path(PPI_DIR, "ppi_hub_annotated.csv")
hub_annotated_written <- tryCatch({
  write.csv(annotated, hub_annotated_file, row.names = FALSE)
  TRUE
}, error = function(e) {
  message("Skipped locked file: ", hub_annotated_file, " (", e$message, ")")
  FALSE
})

g <- read_graph(graph_file, format = "graphml")
node_genes <- if ("gene_symbol" %in% vertex_attr_names(g)) V(g)$gene_symbol else V(g)$name
full_ann <- annotate_genes(node_genes, reactome_long, hallmark_long, memberships)

hub_lookup <- hubs[match(node_genes, hubs$gene), , drop = FALSE]
V(g)$ppi_degree <- ifelse(is.na(hub_lookup$ppi_degree), V(g)$degree_centrality, hub_lookup$ppi_degree)
V(g)$is_ppi_hub <- node_genes %in% hubs$gene
V(g)$in_wgcna_hub <- ifelse(is.na(hub_lookup$in_wgcna_hub), FALSE, hub_lookup$in_wgcna_hub)
V(g)$in_tf_candidates <- ifelse(is.na(hub_lookup$in_tf_candidates), FALSE, hub_lookup$in_tf_candidates)
V(g)$iron_related <- ifelse(is.na(hub_lookup$iron_related), FALSE, hub_lookup$iron_related)
V(g)$wgcna_hub_module <- unname(module_lookup[node_genes])
V(g)$wgcna_hub_module[is.na(V(g)$wgcna_hub_module)] <- ""
V(g)$reactome_pathways <- full_ann$reactome_pathways
V(g)$reactome_NES_range <- full_ann$reactome_NES_range
V(g)$hallmark_pathways <- full_ann$hallmark_pathways
V(g)$hallmark_NES_range <- full_ann$hallmark_NES_range
V(g)$primary_pathway <- full_ann$primary_pathway
V(g)$primary_direction <- full_ann$primary_direction

write_graph(g, file.path(PPI_DIR, "ppi_SC_vs_RC_annotated.graphml"), format = "graphml")

if (hub_annotated_written) {
  message("Wrote: ", hub_annotated_file)
}
message("Wrote: ", file.path(PPI_DIR, "ppi_SC_vs_RC_annotated.graphml"))
