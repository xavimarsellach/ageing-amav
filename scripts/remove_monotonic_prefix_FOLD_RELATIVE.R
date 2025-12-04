#| label: remove_monotonic_prefix_FOLD_RELATIVE.R
#| echo: true
#| eval: true
#| warning: false
#| message: true

# remove_monotonic_prefix_FOLD_RELATIVE.R
# Goal:
#   Read AMAV_DATA.xlsx, detect for each disease the initial stretch
#   where Fold_increase = RelYear + 1, replace those values by NA,
#   and write a new file:
#       AMAV_DATA-no-monotonic_AUTO.xlsx
#
# Requirements:
#   - Packages: readxl, writexl
#   - Input file: "AMAV_DATA.xlsx"
#   - Expected sheets: "AMAV", "FOLD_YEARLY", "FOLD_RELATIVE", "LOG_used_column"

library(readxl)
library(writexl)

# --------- 1. Basic configuration ---------------------------------

# Ensure output directory exists
if (!dir.exists("output")) {
  dir.create("output")
}

input_file  <- "output/AMAV_DATA.xlsx"
output_file <- "output/AMAV_DATA-no-monotonic_AUTO.xlsx"

# Read all sheets we want to preserve
amav_p   <- read_excel(input_file, sheet = "AMAV")
amav_fy  <- read_excel(input_file, sheet = "FOLD_YEARLY")
amav_fr  <- read_excel(input_file, sheet = "FOLD_RELATIVE")
log_used <- read_excel(input_file, sheet = "LOG_used_column")

# Check that the first column is RelYear
if (names(amav_fr)[1] != "RelYear") {
  stop("Expected first column of 'FOLD_RELATIVE' to be 'RelYear'. Please check the file.")
}

# --------- 2. Function to remove monotonic prefix -----------------

# Rule:
#   For each disease column, starting from the smallest RelYear,
#   while value == RelYear + 1 (within a tiny tolerance),
#   set that value to NA. Stop as soon as the condition fails.

remove_monotonic_prefix <- function(df) {
  out <- df
  rel <- out$RelYear
  disease_cols <- names(out)[-1]  # all columns except RelYear
  
  tol <- 1e-9  # numerical tolerance for equality
  
  removed_counts <- numeric(length(disease_cols))
  names(removed_counts) <- disease_cols
  
  for (col in disease_cols) {
    vals <- out[[col]]
    removed <- 0L
    
    for (i in seq_along(vals)) {
      y <- vals[i]
      t <- rel[i]
      
      # Skip if NA
      if (is.na(y)) next
      
      # Check if this value follows the artefactual rule y = RelYear + 1
      if (abs(y - (t + 1)) < tol) {
        vals[i] <- NA_real_
        removed <- removed + 1L
      } else {
        # First value that does not follow y = RelYear + 1:
        # stop removing for this disease
        break
      }
    }
    
    out[[col]] <- vals
    removed_counts[col] <- removed
  }
  
  message("Monotonic prefix removal summary (number of points removed per disease):")
  print(removed_counts)
  
  return(out)
}

# --------- 3. Apply function to FOLD_RELATIVE --------------------------

amav_fr_no_mono <- remove_monotonic_prefix(amav_fr)

# --------- 4. Write new Excel file --------------------------------

write_xlsx(
  list(
    "AMAV"          = amav_p,
    "FOLD_YEARLY"        = amav_fy,
    "FOLD_RELATIVE"        = amav_fr_no_mono,
    "LOG_used_column" = log_used
  ),
  path = output_file
)

message("Done. Written file: ", output_file)