# =============================================================================
# [NEW] run_all.R
# Purpose : Run the whole pipeline (scripts 01-05) in order in one R session.
# Inputs  : DBH CSV (default data/raw/dbh_sample.csv, override with DBH_CSV)
# Outputs : files in outputs/tables/ and outputs/figures/
# Author  : Oyediran Abiodun
# Date    : 2026-10-04
# Usage   : run from the repository root:  Rscript run_all.R
# =============================================================================

pipeline <- c("R/01_data_cleaning.R",
              "R/02_diameter_classes.R",
              "R/03_distribution_fitting.R",
              "R/04_goodness_of_fit.R",
              "R/05_visualization.R")

missing_files <- pipeline[!file.exists(pipeline)]
if(length(missing_files)) {
  stop("Run this script from the repository root. Missing: ", paste(missing_files, collapse = ", "))
}

for(script in pipeline) {
  cat("\n>>>>>>>>>> Running", script, "\n")
  source(script, echo = FALSE)
}
