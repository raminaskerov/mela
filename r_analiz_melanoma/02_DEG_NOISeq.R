# =============================================================================
# 02_DEG_NOISeq.R
# INTEGRATIVE ANALYSIS: ENCORAFENIB RESISTANCE IN MALIGNANT MELANOMA
# GSE283251 | Colakoglu Bergel et al., Scientific Reports 2025
#
# Inputs:  checkpoint_01.RData
# Outputs: 8 DEG tables (4 mRNA + 4 miRNA comparisons)
#          checkpoint_02.RData
#
# NOISeq notes:
#   - norm = "tmm" is applied internally; input must be raw integer counts
#   - NOISeq M = log2(A/B); sign is inverted so positive log2FC = up in B
#   - Threshold: prob >= 0.8 (BH correction is not applicable to NOISeq prob)
#   - miRNA: additional priority score filter applied (see run_noiseq comments)
# =============================================================================


# ── Configuration ─────────────────────────────────────────────────────────────

DATA_DIR   <- "/home/ramin/mela/data/ilk"  # <-- UPDATE THIS PATH
OUTPUT_DIR <- file.path(DATA_DIR, "newoutput")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)


# ── Packages ──────────────────────────────────────────────────────────────────

for (pkg in c("dplyr")) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}
if (!requireNamespace("NOISeq", quietly = TRUE)) BiocManager::install("NOISeq")

suppressPackageStartupMessages({
  library(NOISeq)
  library(dplyr)
})


# ── Load checkpoint ───────────────────────────────────────────────────────────

load(file.path(DATA_DIR, "checkpoint_01.RData"))
# Loaded: mrna_mat, mrna_filtered, mrna_norm,
#         mirna_mat, mirna_filtered, mirna_cpm


# ── Sample metadata ───────────────────────────────────────────────────────────

sample_meta <- data.frame(
  condition = c("Resistant", "Resistant", "Sensitive", "Sensitive"),
  treatment = c("Control", "Drug10nM", "Control", "Drug10nM"),
  row.names = c("A375_RC", "A375_R10", "A375_SC", "A375_S10")
)


# =============================================================================
# NOISeq WRAPPER
# =============================================================================

# Runs a single NOISeq comparison between two sample groups.
#
# Sign correction: NOISeq reports M = log2(A/B). We invert the sign so that
# positive log2FC = higher in B (the second/test group).
#
# priority_filter: when TRUE, computes a priority score
#   (scaled |ranking| * prob) and applies a 40th-percentile cutoff in addition
#   to prob >= q. Used for miRNA to reduce false positives from low-count
#   features that can reach high fold-change with minimal expression.

run_noiseq <- function(count_mat, samples_A, samples_B, label,
                       q = 0.8, priority_filter = FALSE) {

  sel_cols <- c(samples_A, samples_B)
  sub_mat  <- count_mat[, sel_cols, drop = FALSE]
  sub_meta <- data.frame(
    group     = c(rep("A", length(samples_A)), rep("B", length(samples_B))),
    row.names = sel_cols
  )
  sub_mat <- sub_mat[rowSums(sub_mat) > 0, ]

  noiseq_obj <- NOISeq::readData(data = sub_mat, factors = sub_meta)
  res <- NOISeq::noiseq(
    input      = noiseq_obj,
    factor     = "group",
    conditions = c("A", "B"),
    norm       = "tmm",
    replicates = "no",
    k          = 0.5
  )

  full_results         <- NOISeq::degenes(res, q = 0, M = NULL)
  full_results$M       <- full_results$M * -1
  full_results$ranking <- full_results$ranking * -1

  if (priority_filter) {
    scores <- as.numeric(scale(abs(full_results$ranking)) * full_results$prob)
    full_results$priority_score <- scores
    threshold <- quantile(scores, 0.40)
    degs <- full_results[full_results$prob >= q &
                           full_results$priority_score >= threshold, ]
  } else {
    degs <- full_results[full_results$prob >= q, ]
  }

  message(sprintf("  [%s] DEGs (prob >= %.2f): %d", label, q, nrow(degs)))

  list(result    = res,
       full      = full_results,
       degs      = degs,
       label     = label,
       samples_A = samples_A,
       samples_B = samples_B)
}


# ── Save DEG table ────────────────────────────────────────────────────────────

save_deg_table <- function(comp_obj, filename) {
  df          <- as.data.frame(comp_obj$degs)
  df$gene     <- rownames(df)
  df          <- df |> dplyr::rename(log2FC = M, meanA = A_mean, meanB = B_mean)
  write.csv(df, file = file.path(OUTPUT_DIR, filename), row.names = FALSE)
  message(sprintf("  Saved: %s (%d DEGs)", filename, nrow(df)))
  invisible(df)
}


# =============================================================================
# mRNA COMPARISONS
# =============================================================================

message("\n--- mRNA NOISeq ---")

comp1 <- run_noiseq(mrna_filtered, "A375_SC",  "A375_RC",  "SC_vs_RC")
comp2 <- run_noiseq(mrna_filtered, "A375_SC",  "A375_S10", "SC_vs_S10")
comp3 <- run_noiseq(mrna_filtered, "A375_RC",  "A375_R10", "RC_vs_R10")
comp4 <- run_noiseq(mrna_filtered, "A375_S10", "A375_R10", "S10_vs_R10")

deg1 <- save_deg_table(comp1, "DEGs_SC_vs_RC.csv")
deg2 <- save_deg_table(comp2, "DEGs_SC_vs_S10.csv")
deg3 <- save_deg_table(comp3, "DEGs_RC_vs_R10.csv")
deg4 <- save_deg_table(comp4, "DEGs_S10_vs_R10.csv")


# =============================================================================
# miRNA COMPARISONS
# =============================================================================

message("\n--- miRNA NOISeq ---")

mirna_comp1 <- run_noiseq(mirna_filtered, "A375_SC",  "A375_RC",
                          "miRNA_SC_vs_RC",   priority_filter = TRUE)
mirna_comp2 <- run_noiseq(mirna_filtered, "A375_SC",  "A375_S10",
                          "miRNA_SC_vs_S10",  priority_filter = TRUE)
mirna_comp3 <- run_noiseq(mirna_filtered, "A375_RC",  "A375_R10",
                          "miRNA_RC_vs_R10",  priority_filter = TRUE)
mirna_comp4 <- run_noiseq(mirna_filtered, "A375_S10", "A375_R10",
                          "miRNA_S10_vs_R10", priority_filter = TRUE)

mirna_deg1 <- save_deg_table(mirna_comp1, "DEGs_miRNA_SC_vs_RC.csv")
mirna_deg2 <- save_deg_table(mirna_comp2, "DEGs_miRNA_SC_vs_S10.csv")
mirna_deg3 <- save_deg_table(mirna_comp3, "DEGs_miRNA_RC_vs_R10.csv")
mirna_deg4 <- save_deg_table(mirna_comp4, "DEGs_miRNA_S10_vs_R10.csv")


# =============================================================================
# CHECKPOINT SAVE
# =============================================================================

save(mrna_mat, mrna_filtered, mrna_norm,
     mirna_mat, mirna_filtered, mirna_cpm,
     sample_meta,
     comp1, comp2, comp3, comp4,
     mirna_comp1, mirna_comp2, mirna_comp3, mirna_comp4,
     deg1, deg2, deg3, deg4,
     mirna_deg1, mirna_deg2, mirna_deg3, mirna_deg4,
     file = file.path(DATA_DIR, "checkpoint_02.RData"))

message("\ncheckpoint_02.RData saved.")
message(sprintf("  mRNA  DEGs  SC vs RC:   %d", nrow(deg1)))
message(sprintf("  mRNA  DEGs  SC vs S10:  %d", nrow(deg2)))
message(sprintf("  mRNA  DEGs  RC vs R10:  %d", nrow(deg3)))
message(sprintf("  mRNA  DEGs  S10 vs R10: %d", nrow(deg4)))
message(sprintf("  miRNA DEGs  SC vs RC:   %d", nrow(mirna_deg1)))
message(sprintf("  miRNA DEGs  SC vs S10:  %d", nrow(mirna_deg2)))
message(sprintf("  miRNA DEGs  RC vs R10:  %d", nrow(mirna_deg3)))
message(sprintf("  miRNA DEGs  S10 vs R10: %d", nrow(mirna_deg4)))
