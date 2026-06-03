# =============================================================================
# 01_preprocessing.R
# INTEGRATIVE ANALYSIS: ENCORAFENIB RESISTANCE IN MALIGNANT MELANOMA
# GSE283251 | Colakoglu Bergel et al., Scientific Reports 2025
#
# Inputs:  8 GEO count files (mRNA + miRNA, 4 samples each)
# Outputs: checkpoint_01.RData
#          Objects: mrna_mat, mrna_filtered, mrna_norm,
#                   mirna_mat, mirna_filtered, mirna_cpm
# =============================================================================


# ── Configuration ─────────────────────────────────────────────────────────────

DATA_DIR   <- "/home/ramin/mela/data/ilk"  # <-- UPDATE THIS PATH
OUTPUT_DIR <- file.path(DATA_DIR, "newoutput")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)


# ── Packages ──────────────────────────────────────────────────────────────────

for (pkg in c("readr", "dplyr")) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})


# =============================================================================
# SECTION 1: LOAD COUNT FILES
# =============================================================================

# GEO files are UTF-16LE encoded gzip TSVs

read_count_file <- function(path) {
  readr::read_tsv(
    path,
    locale         = readr::locale(encoding = "UTF-16LE"),
    show_col_types = FALSE,
    name_repair    = "minimal",
    progress       = FALSE
  )
}


# ── mRNA ──────────────────────────────────────────────────────────────────────

mrna_files <- list(
  A375_RC  = file.path(DATA_DIR, "GSM8658654_A375-RK_mRNA_count.txt.gz"),
  A375_R10 = file.path(DATA_DIR, "GSM8658655_A375-R-10_mRNA_count.txt.gz"),
  A375_SC  = file.path(DATA_DIR, "GSM8658656_A375-S-K_mRNA_count.txt.gz"),
  A375_S10 = file.path(DATA_DIR, "GSM8658657_A375-S-10_mRNA_count.txt.gz")
)

mrna_raw <- list()
for (nm in names(mrna_files)) {
  mrna_raw[[nm]] <- tryCatch(
    read_count_file(mrna_files[[nm]]),
    error = function(e) stop(sprintf("Failed to read %s: %s", nm, e$message))
  )
  message(sprintf("  Loaded %-10s: %d rows", nm, nrow(mrna_raw[[nm]])))
}

# Keep protein-coding genes only; sum counts for duplicate gene symbols
extract_mrna_counts <- function(df, sample_name) {
  count_col <- grep("Read_Count", colnames(df), value = TRUE)[1]
  df |>
    dplyr::select(Gene_Symbol, gene_biotype, all_of(count_col)) |>
    dplyr::rename(count = all_of(count_col)) |>
    dplyr::filter(
      gene_biotype == "protein_coding",
      !is.na(Gene_Symbol),
      Gene_Symbol  != "",
      Gene_Symbol  != "."
    ) |>
    dplyr::group_by(Gene_Symbol) |>
    dplyr::summarise(count = sum(count, na.rm = TRUE), .groups = "drop") |>
    dplyr::rename(!!sample_name := count)
}

mrna_counts_list <- mapply(extract_mrna_counts, mrna_raw, names(mrna_raw),
                           SIMPLIFY = FALSE)

for (nm in names(mrna_counts_list)) {
  message(sprintf("  %s: %d protein-coding genes", nm, nrow(mrna_counts_list[[nm]])))
}

mrna_counts <- Reduce(
  function(x, y) dplyr::full_join(x, y, by = "Gene_Symbol"),
  mrna_counts_list
)

mrna_counts                   <- as.data.frame(mrna_counts)
mrna_counts[is.na(mrna_counts)] <- 0L
mrna_mat                      <- as.matrix(mrna_counts[, -1, drop = FALSE])
rownames(mrna_mat)             <- mrna_counts$Gene_Symbol
storage.mode(mrna_mat)         <- "integer"


# ── miRNA ─────────────────────────────────────────────────────────────────────

mirna_files <- list(
  A375_RC  = file.path(DATA_DIR, "GSM8658654_A375-RK_miRNA_count.txt.gz"),
  A375_R10 = file.path(DATA_DIR, "GSM8658655_A375-R-10_miRNA_count.txt.gz"),
  A375_SC  = file.path(DATA_DIR, "GSM8658656_A375-S-K_miRNA_count.txt.gz"),
  A375_S10 = file.path(DATA_DIR, "GSM8658657_A375-S-10_miRNA_count.txt.gz")
)

mirna_raw <- list()
for (nm in names(mirna_files)) {
  mirna_raw[[nm]] <- tryCatch(
    read_count_file(mirna_files[[nm]]),
    error = function(e) stop(sprintf("Failed to read %s: %s", nm, e$message))
  )
  message(sprintf("  Loaded %-10s: %d rows", nm, nrow(mirna_raw[[nm]])))
}

extract_mirna_counts <- function(df, sample_name) {
  count_col <- grep("Read_Count", colnames(df), value = TRUE)[1]
  df |>
    dplyr::select(Mature_ID, all_of(count_col)) |>
    dplyr::rename(!!sample_name := all_of(count_col)) |>
    dplyr::filter(!is.na(Mature_ID), Mature_ID != "")
}

mirna_counts <- Reduce(
  function(x, y) dplyr::full_join(x, y, by = "Mature_ID"),
  mapply(extract_mirna_counts, mirna_raw, names(mirna_raw), SIMPLIFY = FALSE)
)

mirna_counts                    <- as.data.frame(mirna_counts)
mirna_counts[is.na(mirna_counts)] <- 0L
mirna_mat                       <- as.matrix(mirna_counts[, -1, drop = FALSE])
rownames(mirna_mat)              <- mirna_counts$Mature_ID
storage.mode(mirna_mat)          <- "integer"


# =============================================================================
# SECTION 2: FILTERING
# =============================================================================

# ── mRNA: keep genes with CPM > 0.5 in at least one sample ───────────────────

mrna_cpm_check <- sweep(mrna_mat, 2, colSums(mrna_mat) / 1e6, FUN = "/")
keep_genes     <- rowSums(mrna_cpm_check > 0.5) >= 1
mrna_filtered  <- mrna_mat[keep_genes, ]
message(sprintf("mRNA CPM filter: %d / %d genes retained",
                sum(keep_genes), nrow(mrna_mat)))

# ── Histone filter ────────────────────────────────────────────────────────────
# Replication-dependent canonical histones are removed before enrichment
# analysis to prevent multi-mapping artifacts from inflating pathway terms.
# Genes in protect_genes are retained even if they match the pattern.

histone_pattern <- "^HIST[0-9]|^H[1-4][A-Z]|^HIST|^H2A|^H2B|^H3[^F]|^H4C"
protect_genes   <- c("H3F3A", "H3F3B", "HMOX1")

histone_candidates <- grep(histone_pattern, rownames(mrna_filtered), value = TRUE)
histone_remove     <- setdiff(histone_candidates, protect_genes)
mrna_filtered      <- mrna_filtered[!rownames(mrna_filtered) %in% histone_remove, ]

message(sprintf("Histone filter: %d genes removed", length(histone_remove)))
message(sprintf("mRNA after histone filter: %d genes", nrow(mrna_filtered)))

# ── miRNA: keep features with count > 2 in at least one sample ───────────────

mirna_keep     <- rowSums(mirna_mat > 2) >= 1
mirna_filtered <- mirna_mat[mirna_keep, ]
message(sprintf("miRNA filter: %d / %d retained",
                sum(mirna_keep), nrow(mirna_mat)))


# =============================================================================
# SECTION 3: NORMALIZATION
# =============================================================================

# Geometric mean library-size scaling for visualization and WGCNA input.
# NOISeq applies its own internal TMM normalization; these objects are
# not passed to NOISeq.

libsize_normalize <- function(count_mat) {
  lib_sizes <- colSums(count_mat)
  ref_lib   <- exp(mean(log(lib_sizes[lib_sizes > 0])))
  sweep(count_mat, 2, lib_sizes / ref_lib, FUN = "/")
}

mrna_norm <- libsize_normalize(mrna_mat)
mirna_cpm <- sweep(mirna_mat, 2, colSums(mirna_mat), FUN = "/") * 1e6


# =============================================================================
# CHECKPOINT SAVE
# =============================================================================

save(mrna_mat, mrna_filtered, mrna_norm,
     mirna_mat, mirna_filtered, mirna_cpm,
     file = file.path(DATA_DIR, "checkpoint_01.RData"))

message("\ncheckpoint_01.RData saved.")
message(sprintf("  mrna_filtered:  %d genes  x %d samples",
                nrow(mrna_filtered), ncol(mrna_filtered)))
message(sprintf("  mirna_filtered: %d miRNAs x %d samples",
                nrow(mirna_filtered), ncol(mirna_filtered)))
