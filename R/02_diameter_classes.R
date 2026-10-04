# =============================================================================
# [NEW] R/02_diameter_classes.R
# Purpose : Assign trees to diameter classes and tabulate class frequencies.
#           Minimal working code written because the source script had no
#           diameter-class step.
# Inputs  : objects from R/01_data_cleaning.R (data, xmin, xmax, out_tab_dir)
#           optional env var DBH_CLASS_WIDTH (cm, default 10)
# Outputs : outputs/tables/diameter_class_table.csv ; object dbh_class_table
# Author  : Oyediran Abiodun
# Date    : 2026-10-04
# =============================================================================

if(!exists("data") || !exists("xmin") || !exists("out_tab_dir")) {
  stop("Run R/01_data_cleaning.R first (objects 'data', 'xmin', 'out_tab_dir' not found).")
}

class_width <- as.numeric(Sys.getenv("DBH_CLASS_WIDTH", unset = "10"))
if(is.na(class_width) || class_width <= 0) stop("DBH_CLASS_WIDTH must be a positive number.")

# Classes are [lower, upper); the last break lies above the largest tree
class_breaks <- seq(floor(xmin / class_width) * class_width,
                    (floor(xmax / class_width) + 1) * class_width,
                    by = class_width)

dbh_class <- cut(data, breaks = class_breaks, right = FALSE)
class_freq <- as.vector(table(dbh_class))

dbh_class_table <- data.frame(
  Class       = levels(dbh_class),
  Lower_cm    = head(class_breaks, -1),
  Upper_cm    = class_breaks[-1],
  Midpoint_cm = head(class_breaks, -1) + class_width / 2,
  Trees       = class_freq,
  Rel_Freq    = class_freq / length(data),
  stringsAsFactors = FALSE
)

cat("\nDiameter class table (class width =", class_width, "cm):\n")
print(dbh_class_table)

write.csv(dbh_class_table, file.path(out_tab_dir, "diameter_class_table.csv"), row.names = FALSE)
cat("\nDiameter class table saved to", file.path(out_tab_dir, "diameter_class_table.csv"), "\n")
