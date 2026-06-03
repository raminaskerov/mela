#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(pathview)
})

ROOT_DIR <- "C:/bioinformatics"
SESSION_FILE <- file.path(ROOT_DIR, "r_analiz_melanoma", "finished(supheli).RData")
PATHVIEW_DIR <- file.path(ROOT_DIR, "output", "pathview")

dir.create(PATHVIEW_DIR, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(SESSION_FILE)) {
  stop("Missing required session file: ", SESSION_FILE)
}

load(SESSION_FILE)

if (!exists("comp1")) {
  stop("Loaded session does not contain comp1: ", SESSION_FILE)
}

full_df <- as.data.frame(comp1$full)
if (!"M" %in% colnames(full_df)) {
  stop("comp1$full must contain an M column with SC vs RC log2FC values.")
}

gene_fc <- as.numeric(full_df$M)
names(gene_fc) <- rownames(full_df)
gene_fc <- gene_fc[!is.na(gene_fc) & !is.na(names(gene_fc)) & nzchar(names(gene_fc))]

mapped <- suppressMessages(
  clusterProfiler::bitr(
    names(gene_fc),
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  )
)

fc_df <- data.frame(
  SYMBOL = names(gene_fc),
  log2FC = unname(gene_fc),
  stringsAsFactors = FALSE
)
fc_df <- merge(mapped, fc_df, by = "SYMBOL", all.x = FALSE, all.y = FALSE)
fc_df <- fc_df[order(fc_df$SYMBOL, -abs(fc_df$log2FC), fc_df$ENTREZID), ]
fc_df <- fc_df[!duplicated(fc_df$SYMBOL), ]
fc_df <- fc_df[order(fc_df$ENTREZID, -abs(fc_df$log2FC), fc_df$SYMBOL), ]
fc_df <- fc_df[!duplicated(fc_df$ENTREZID), ]

gene_fc_entrez <- fc_df$log2FC
names(gene_fc_entrez) <- fc_df$ENTREZID

pathway_ids <- c(
  "hsa04978",
  "hsa00190",
  "hsa04110",
  "hsa04210",
  "hsa04668",
  "hsa04630"
)

message("Mapped ", length(gene_fc_entrez), " unique Entrez IDs from ", length(gene_fc), " symbols.")
message("Writing Pathview output to: ", PATHVIEW_DIR)

old_files <- list.files(getwd(), pattern = "\\.(xml|png)$", full.names = TRUE)

for (pathway_id in pathway_ids) {
  message("Running Pathview for ", pathway_id, "...")
  tryCatch(
    {
      pathview::pathview(
        gene.data = gene_fc_entrez,
        pathway.id = pathway_id,
        species = "hsa",
        kegg.dir = PATHVIEW_DIR,
        out.suffix = "SC_vs_RC",
        low = "blue",
        mid = "white",
        high = "red",
        limit = list(gene = 3, cpd = 1)
      )
      message("Finished: ", pathway_id)
    },
    error = function(e) {
      message("Skipping ", pathway_id, " after Pathview failure: ", e$message)
    }
  )
}

new_pngs <- setdiff(list.files(getwd(), pattern = "\\.png$", full.names = TRUE), old_files)
if (length(new_pngs) > 0) {
  file.copy(new_pngs, file.path(PATHVIEW_DIR, basename(new_pngs)), overwrite = TRUE)
  unlink(new_pngs)
}

new_xmls <- setdiff(list.files(getwd(), pattern = "\\.xml$", full.names = TRUE), old_files)
if (length(new_xmls) > 0) {
  unlink(new_xmls)
}

message("PNG files in output directory:")
print(list.files(PATHVIEW_DIR, pattern = "\\.png$", full.names = TRUE))
