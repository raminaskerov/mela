# =============================================================================
# main pipeline 2.R
# INTEGRATIVE ANALYSIS: ENCORAFENIB RESISTANCE IN MALIGNANT MELANOMA
#
# This file is now a thin runner only.
#
# Stage graph:
#   01 -> 02 -> 03 -> 04
#   03 -> 05
#   03 -> 06
#   03 -> 07
#   06 + 07 -> 08
#
# Each stage is a standalone script that reads the latest checkpoint and writes
# its own output checkpoint or analysis tables.
# =============================================================================


get_script_dir <- function() {
  cmd <- commandArgs(trailingOnly = FALSE)
  file_arg <- sub("^--file=", "", cmd[grep("^--file=", cmd)])
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(file_arg[1])))
  }
  getwd()
}

SCRIPT_DIR <- get_script_dir()

stage_files <- c(
  "01_preprocessing.R",
  "02_DEG_NOISeq.R",
  "03_enrichment.R",
  "04_visualizations.R",
  "05_mirna_mrna_integration.R",
  "06_TF_enrichment.R",
  "07_WGCNA.R",
  "08_legacy_wgcna_tf_linkage.R"
)

message("Running melanoma resistance pipeline from: ", SCRIPT_DIR)

for (stage in stage_files) {
  stage_path <- file.path(SCRIPT_DIR, stage)
  if (!file.exists(stage_path)) {
    stop("Missing stage script: ", stage_path)
  }
  message("\n=== ", stage, " ===")
  source(stage_path, local = new.env(parent = globalenv()))
}

message("\nPipeline finished.")
