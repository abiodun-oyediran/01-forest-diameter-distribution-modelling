# =============================================================================
# [NEW] scripts/make_sample_data.R
# Purpose : Generate a reproducible SYNTHETIC tropical-forest DBH sample.
#           The data are simulated; they are not field measurements.
# Inputs  : none (optional environment variables: DBH_N, DBH_SEED, DBH_CSV)
# Outputs : data/raw/dbh_sample.csv  (one column, header "Dbh(cm)")
# Author  : Oyediran Abiodun
# Date    : 2026-10-04
# Usage   : run from the repository root:  Rscript scripts/make_sample_data.R
# =============================================================================

n        <- as.integer(Sys.getenv("DBH_N", unset = "500"))
seed     <- as.integer(Sys.getenv("DBH_SEED", unset = "2026"))
out_path <- Sys.getenv("DBH_CSV", unset = file.path("data", "raw", "dbh_sample.csv"))

if (is.na(n) || n < 30) stop("DBH_N must be an integer >= 30.")
if (is.na(seed)) stop("DBH_SEED must be an integer.")

set.seed(seed)

# Right-skewed (reverse-J) size structure above a 10 cm inventory threshold,
# plus a small share of large canopy / emergent trees.
n_large <- round(0.06 * n)
small   <- 10 + rweibull(n - n_large, shape = 1.2, scale = 24)
large   <- rlnorm(n_large, meanlog = log(85), sdlog = 0.3)

dbh <- round(c(small, large), 1)
dbh <- pmin(pmax(dbh, 10), 180)
dbh <- sample(dbh)

# Header "Dbh(cm)" is read by read.csv() as the column name "Dbh.cm."
out <- data.frame(`Dbh(cm)` = dbh, check.names = FALSE)

dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
res <- tryCatch(
  write.csv(out, out_path, row.names = FALSE),
  error = function(e) stop("Could not write ", out_path, ": ", conditionMessage(e))
)

cat("Wrote", nrow(out), "synthetic DBH records to", out_path, "\n")
cat("Summary (cm): min =", min(dbh), "| mean =", round(mean(dbh), 1),
    "| median =", median(dbh), "| max =", max(dbh), "\n")
