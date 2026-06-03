# =============================================================================
# FIX v2 — readr-based loader (CONFIRMED WORKING after correct DATA_DIR)
#
# This version uses readr::read_tsv() with locale(encoding="UTF-16LE").
# It works once DATA_DIR points directly to the folder containing the .gz files
# (no spaces or special characters in path).
#
# Applied fixes vs original pipeline:
#   1. readr instead of rawToChar() → handles large UTF-16LE gzipped files
#   2. as.data.frame() before NA assignment → tibble doesn't support df[is.na]
#   3. drop=FALSE on [,-1] → guards against single-column edge case
# =============================================================================


# ── 0.4  Helper: read UTF-16LE gzipped count files ───────────────────────────

read_count_file <- function(path) {
  readr::read_tsv(
    path,
    locale         = readr::locale(encoding = "UTF-16LE"),
    show_col_types = FALSE,
    name_repair    = "minimal",
    progress       = FALSE
  )
}


# ── 0.5  Load all mRNA count files ───────────────────────────────────────────

message("Loading mRNA count files...")

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
  message(sprintf("  Loaded %-10s: %d rows x %d cols | first cols: %s",
                  nm,
                  nrow(mrna_raw[[nm]]),
                  ncol(mrna_raw[[nm]]),
                  paste(head(colnames(mrna_raw[[nm]]), 4), collapse = ", ")))
}


# ── extract_counts: Gene_Symbol → count, protein-coding only ─────────────────

extract_counts <- function(df, sample_name) {
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

mrna_counts_list <- mapply(extract_counts,
                           mrna_raw, names(mrna_raw),
                           SIMPLIFY = FALSE)

# Report row counts before joining (useful for diagnosing 0-row issues)
for (nm in names(mrna_counts_list)) {
  message(sprintf("  extract_counts(%s): %d protein-coding genes",
                  nm, nrow(mrna_counts_list[[nm]])))
}

mrna_counts <- Reduce(
  function(x, y) dplyr::full_join(x, y, by = "Gene_Symbol"),
  mrna_counts_list
)

# FIX: tibble → base data.frame before NA assignment
mrna_counts                <- as.data.frame(mrna_counts)
mrna_counts[is.na(mrna_counts)] <- 0L
mrna_mat                   <- as.matrix(mrna_counts[, -1, drop = FALSE])
rownames(mrna_mat)         <- mrna_counts$Gene_Symbol
storage.mode(mrna_mat)     <- "integer"

message(sprintf("mRNA matrix: %d genes x %d samples", nrow(mrna_mat), ncol(mrna_mat)))

# Validation against Python audit values
cat("\n  Validation (expected from Python audit):\n")
cat("  Gene    | RC   | R10  | SC   | S10\n")
cat("  --------|------|------|------|------\n")
for (g in c("FTH1","NCOA4","DUSP6","PIK3R1")) {
  if (g %in% rownames(mrna_mat)) {
    cat(sprintf("  %-7s | %4d | %4d | %4d | %4d\n", g,
                mrna_mat[g, "A375_RC"],  mrna_mat[g, "A375_R10"],
                mrna_mat[g, "A375_SC"],  mrna_mat[g, "A375_S10"]))
  } else {
    cat(sprintf("  %-7s | MISSING\n", g))
  }
}


# Expected:
#   FTH1    |  627 | 1361 |  162 |  263
#   NCOA4   |   33 |   28 |   31 |   19
#   DUSP6   |  127 |   24 |  117 |   12
#   PIK3R1  |   22 |  108 |   25 |   67


# ── 0.6  Load all miRNA count files ──────────────────────────────────────────

message("\nLoading miRNA count files...")

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
  message(sprintf("  Loaded %-10s: %d miRNAs", nm, nrow(mirna_raw[[nm]])))
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

# FIX: tibble → base data.frame before NA assignment
mirna_counts                <- as.data.frame(mirna_counts)
mirna_counts[is.na(mirna_counts)] <- 0L
mirna_mat                   <- as.matrix(mirna_counts[, -1, drop = FALSE])
rownames(mirna_mat)         <- mirna_counts$Mature_ID
storage.mode(mirna_mat)     <- "integer"

message(sprintf("miRNA matrix: %d miRNAs x %d samples", nrow(mirna_mat), ncol(mirna_mat)))

# Validate hsa-miR-140-3p (expected: RC=58, R10=703, SC=22, S10=10)
if ("hsa-miR-140-3p" %in% rownames(mirna_mat)) {
  cat(sprintf("  miR-140-3p: RC=%d  R10=%d  SC=%d  S10=%d\n",
              mirna_mat["hsa-miR-140-3p", "A375_RC"],
              mirna_mat["hsa-miR-140-3p", "A375_R10"],
              mirna_mat["hsa-miR-140-3p", "A375_SC"],
              mirna_mat["hsa-miR-140-3p", "A375_S10"]))
  # Expected: RC=58  R10=703  SC=22  S10=10
}


# ── 0.7  Normalization ────────────────────────────────────────────────────────

tmm_normalize <- function(count_mat) {
  lib_sizes <- colSums(count_mat)
  ref_lib   <- exp(mean(log(lib_sizes[lib_sizes > 0])))
  scale_fac <- lib_sizes / ref_lib
  sweep(count_mat, 2, scale_fac, FUN = "/")
}

mrna_norm <- tmm_normalize(mrna_mat)
mirna_cpm <- sweep(mirna_mat, 2, colSums(mirna_mat), FUN = "/") * 1e6

message("\nData loading complete. Proceeding to Section 1 (NOISeq)...")

# =============================================================================
# CONTINUE: paste Section 1 onward from melanoma_resistance_pipeline.R
# (everything from "SECTION 1: NOISeq DEG ANALYSIS" is unchanged)
# =============================================================================