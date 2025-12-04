# Figure 2 and Figure 3 from AMAV_DATA.xlsx (sheet "AMAV")
# Physical vs mental disease prevalence (AMAV)
# ------------------------------------------------
# This script reads the master file `AMAV_DATA.xlsx` and generates:
#   - Figure 2: physical diseases
#   - Figure 3: mental/behavioural traits
#
# Adjust `file_path` below if your AMAV_DATA.xlsx lives in another folder.

library(readxl)
library(ggplot2)
library(reshape2)
library(scales)

# ---------- Config ----------

file_path <- "data/AMAV_DATA.xlsx"

# Physical diseases for Figure 2
physical_diseases <- c(
  "Allergic Rhinitis",
  "Alopecia Areata",
  "Anaemia",
  "Asthma",
  "Atopic Dermatitis",
  "Cancer",
  "Coeliac Disease",
  "Food Allergies",
  "Hypertension",
  "Lupus",
  "Myopia",
  "Obesity",
  "Osteoarthritis",
  "Rheumatoid Arthritis",
  "Type 1 Diabetes",
  "Type 2 Diabetes"
)

# Mental / behavioural traits for Figure 3
mental_diseases <- c(
  "ADHD",
  "Alzheimer",
  "Anxiety",
  "Autism",
  "Bipolar Disorder",
  "Dementia",
  "Depression",
  "Eating Disorders",
  "Handedness",
  "Homosexuality - Men",
  "Homosexuality - Women",
  "Personality Disorders",
  "Schizophrenia",
  "Suicide",
  "Transgender & Gender Dysphoria"
)

# Common aesthetics (same as old Figure_2.R / Figure_3.R)
manual_colours <- c(
  "#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e",
  "#e6ab02", "#a6761d", "#666666", "#8c564b", "#d62728",
  "#9467bd", "#2ca02c", "#1f77b4", "#ff7f0e", "#bcbd22",
  "#17becf"
)

manual_shapes <- c(
  16, 17, 18, 19, 15, 0, 1, 2,
  3, 4, 5, 6, 7, 8, 9, 10
)

# ---------- Helper to build the plot ----------

make_amav_plot <- function(df_wide, plot_title) {
  # df_wide must have a "Year" column and one column per disease/trait
  
  # Drop rows where *all* disease columns are NA (keep rows with at least one value)
  if (!"Year" %in% names(df_wide)) {
    stop("Input data frame must contain a 'Year' column.")
  }
  
  if (ncol(df_wide) <= 1) {
    stop("Input data frame must contain at least one disease column besides 'Year'.")
  }
  
  disease_cols <- df_wide[, -1, drop = FALSE]
  keep_rows <- rowSums(is.na(disease_cols)) < ncol(disease_cols)
  data_filtered <- df_wide[keep_rows, ]
  
  # Reshape to long format
  data_melted <- melt(
    data_filtered,
    id.vars = "Year",
    variable.name = "Disease",
    value.name = "value"
  )
  
  # Build the plot (same structure as original scripts)
  ggplot(
    data_melted,
    aes(x = Year, y = value, colour = Disease, shape = Disease)
  ) +
    geom_smooth(method = "loess", se = FALSE, linewidth = 1, na.rm = TRUE) +
    geom_point(size = 3, na.rm = TRUE) +
    scale_y_log10(labels = scales::percent_format(accuracy = 0.01)) +
    labs(
      title = plot_title,
      x = "Year",
      y = "Accumulated Mean Annual Variation (AMAV)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
      axis.title.x = element_text(size = 16),
      axis.title.y = element_text(size = 16),
      axis.text.x  = element_text(size = 14),
      axis.text.y  = element_text(size = 14)
    ) +
    scale_colour_manual(values = manual_colours) +
    scale_shape_manual(values = manual_shapes)
}

# ---------- Read master file and subset ----------

amav_p <- read_excel(file_path, sheet = "AMAV")

# Subset for Figure 2 (physical diseases)
fig2_data <- amav_p[, c("Year", physical_diseases)]

# Subset for Figure 3 (mental / behavioural traits)
fig3_data <- amav_p[, c("Year", mental_diseases)]

# ---------- Build plots ----------

figure_2 <- make_amav_plot(
  df_wide    = fig2_data,
  plot_title = "Epigenetic Degeneration Examples. Set Number One"
)

figure_3 <- make_amav_plot(
  df_wide    = fig3_data,
  plot_title = "Epigenetic Degeneration Examples. Set Number Two. Mental Health Issues"
)

# Print static versions (for scripts / small screens)
print(figure_2)
print(figure_3)

# Save to files in output/ (like the GitHub pipeline)
dir.create("output", showWarnings = FALSE)

ggsave(file.path("output", "Figure_2_AMAV_physical.png"),
       figure_2, width = 8, height = 6, dpi = 300)

ggsave(file.path("output", "Figure_3_AMAV_mental.png"),
       figure_3, width = 8, height = 6, dpi = 300)