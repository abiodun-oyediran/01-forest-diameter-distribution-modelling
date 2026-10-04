# =============================================================================
# [ADAPTED] R/01_data_cleaning.R
# Purpose : Load packages, read the DBH CSV, apply minimal cleaning and compute
#           the summary values used by the later scripts.
# Inputs  : DBH CSV with a column "Dbh(cm)" (read by read.csv as "Dbh.cm.").
#           Default: data/raw/dbh_sample.csv. Override with env var DBH_CSV.
# Outputs : objects in the R session: data, n, xmin, xmax, mean_data, sd_data,
#           out_fig_dir, out_tab_dir
# Author  : Oyediran Abiodun
# Date    : 2026-10-04
# Notes   : [adapted from my script: package loading, input data, summary values]
#           Input path and output folders are now configurable (previously
#           file.choose()). The cleaning block is [NEW].
# =============================================================================

# ---------- REQUIRED PACKAGES ----------
required <- c("fitdistrplus", "goftest")
to_install <- required[!(required %in% installed.packages()[,1])]
if(length(to_install)) install.packages(to_install, repos = "https://cloud.r-project.org")

library(fitdistrplus)
library(goftest)

# ---------- PATHS ----------
dbh_csv     <- Sys.getenv("DBH_CSV", unset = file.path("data", "raw", "dbh_sample.csv"))
out_fig_dir <- Sys.getenv("OUT_FIG_DIR", unset = file.path("outputs", "figures"))
out_tab_dir <- Sys.getenv("OUT_TAB_DIR", unset = file.path("outputs", "tables"))
for(dir_i in c(out_fig_dir, out_tab_dir)) dir.create(dir_i, recursive = TRUE, showWarnings = FALSE)

# ---------- INPUT DATA ----------
if(!file.exists(dbh_csv)) {
  stop("Input file not found: ", dbh_csv,
       "\nRun 'Rscript scripts/make_sample_data.R' or set the DBH_CSV environment variable.")
}
HD<- read.csv(dbh_csv)   # [adapted from my script: file.choose() replaced by dbh_csv]
names(HD)
if(!"Dbh.cm." %in% names(HD)) {
  stop("Column 'Dbh.cm.' (CSV header 'Dbh(cm)') not found in ", dbh_csv)
}
data<-HD$Dbh.cm.

# ---------- MINIMAL CLEANING ----------
# Keep finite, positive numeric DBH values only; report what was removed.
n_raw <- length(data)
data <- suppressWarnings(as.numeric(data))
data <- data[is.finite(data) & data > 0]
cat("DBH records read:", n_raw, "| removed (NA, non-numeric or <= 0):",
    n_raw - length(data), "| kept:", length(data), "\n")
if(length(data) < 5) stop("Fewer than 5 valid DBH values remain; cannot fit distributions.")

n <- length(data)
xmin <- min(data)
xmax <- max(data)
mean_data <- mean(data)
sd_data <- sd(data)
