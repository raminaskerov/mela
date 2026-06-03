tmm_normalize <- function(count_mat) {
  lib_sizes <- colSums(count_mat)
  ref_lib   <- exp(mean(log(lib_sizes[lib_sizes > 0])))
  scale_fac <- lib_sizes / ref_lib
  sweep(count_mat, 2, scale_fac, FUN = "/")
}

mirna_cpm <- sweep(mirna_mat, 2, colSums(mirna_mat), FUN = "/") * 1e6

message("\nData loading complete. Proceeding to Section 1 (NOISeq)...")


# =============================================================================
# SECTION 1: NOISeq DEG ANALYSIS (PAPER REPRODUCTION)
# =============================================================================

# NOISeq is designed for single-sample (no-replicate) differential expression.
# It uses a noise distribution estimated from within-sample variability.

message("\n=== SECTION 1: NOISeq DEG Analysis ===")

# ── 1.1  Build NOISeq data object ────────────────────────────────────────────

# Sample metadata
sample_meta <- data.frame(
  condition = c("Resistant", "Resistant", "Sensitive", "Sensitive"),
  treatment = c("Control", "Drug10nM", "Control", "Drug10nM"),
  row.names = c("A375_RC", "A375_R10", "A375_SC", "A375_S10")
)



# ── 1.2  Run NOISeq comparisons ───────────────────────────────────────────────


# The three comparisons from the paper
# Note: NOISeq factor levels must match your sample_meta 'condition' column
# We'll restructure for each comparison by subsetting

run_noiseq_comparison <- function(count_mat, meta, samples_A, samples_B,
                                  label, q = 0.8) {
  sel_cols  <- c(samples_A, samples_B)
  sub_mat   <- count_mat[, sel_cols, drop = FALSE]
  sub_meta  <- data.frame(
    group = c(rep("A", length(samples_A)), rep("B", length(samples_B))),
    row.names = sel_cols
  )
  keep    <- rowSums(sub_mat) > 0
  sub_mat <- sub_mat[keep, ]
  
  noiseq_obj <- NOISeq::readData(data = sub_mat, factors = sub_meta)
  res <- NOISeq::noiseq(
    input      = noiseq_obj,
    factor     = "group",
    conditions = c("A", "B"),
    norm       = "tmm",
    replicates = "no",
    k          = 0.5
  )
  
  full_results                 <- NOISeq::degenes(res, q = 0, M = NULL)
  scores <- scale(abs(full_results$ranking)) * full_results$prob^1
  full_results$priority_score <- as.numeric(scores)
  
  threshold     <- quantile(full_results$priority_score, 0.40)
  degs_priority <- full_results[full_results$priority_score >= threshold &
                                  full_results$prob >= q, ]
  
  message(sprintf("  [%s] DEGs (prob>=%.2f): %d",
                  label, q, nrow(degs_priority)))
  
  list(result   = res,
       degs     = degs_priority,
       full     = full_results,
       label    = label,
       samples_A = samples_A,
       samples_B = samples_B)
}



# ── 1.3  Run NOISeq for miRNAs ───────────────────────────────────────────────

# Filter low-count miRNAs
mirna_keep <- rowSums(mirna_mat > 2) >= 1
mirna_filtered <- mirna_mat[mirna_keep, ]
message(sprintf("miRNAs after filtering: %d", nrow(mirna_filtered)))

mirna_comp1 <- run_noiseq_comparison(
  mirna_filtered,
  meta = sample_meta,
  samples_A = "A375_SC",
  samples_B = "A375_RC",
  label = "miRNA_SC_vs_RC"
)

mirna_comp2 <- run_noiseq_comparison(
  mirna_filtered,
  meta = sample_meta,
  samples_A = "A375_SC",
  samples_B = "A375_S10",
  label = "miRNA_SC_vs_S10"
)

mirna_comp3 <- run_noiseq_comparison(
  mirna_filtered,
  meta = sample_meta,
  samples_A = "A375_RC",
  samples_B = "A375_R10",
  label = "miRNA_RC_vs_R10"
)

mirna_comp4 <- run_noiseq_comparison(
  mirna_filtered,
  meta = sample_meta,
  samples_A = "A375_S10",
  samples_B = "A375_R10",
  label = "miRNA_S10_vs_R10"
)
