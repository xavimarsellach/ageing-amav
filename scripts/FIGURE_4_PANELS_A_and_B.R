#| label: FIGURE_4_PANELS_A_and_B.R
#| echo: true
#| eval: true
#| warning: false
#| message: false

# FIGURE_4_PANELS_A_and_B.R
# Panel Figure 4A–B: per-disease trajectories + global geometric means

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(patchwork)

# --------- 1. CONFIG ---------------------------------------------------

file_path <- "output/AMAV_DATA.xlsx"   # adjust path if needed

physical_diseases <- c(
  "Allergic Rhinitis","Alopecia Aereata","Anaemia","Asthma","Atopic Dermatitis",
  "Cancer","Coeliac Disease","Food Allergies","Hypertension","Lupus",
  "Myopia","Obesity","Osteoarthritis","Rheumatoid Arthritis",
  "Type 1 Diabetes","Type 2 Diabetes"
)

mental_diseases <- c(
  "ADHD","Alzheimer","Anxiety","Autism","Bipolar Disorder","Dementia",
  "Depression","Eating Disorders","Handedness","Homosexuality (Men)",
  "Homosexuality (Women)","Personality Disorders","Schizophrenia",
  "Suicide","Transgender & Gender Dysphoria"
)

# Common colours for all panels
group_cols <- c(Mental = "red", Physical = "blue")

# Ensure output directory exists
if (!dir.exists("output")) {
  dir.create("output")
}

# --------- 2. READ + LONG FORMAT (COMMON OBJECT) ----------------------

amav_fr <- read_excel(file_path, sheet = "FOLD_RELATIVE")

rel_col <- "RelYear"
disease_cols <- setdiff(names(amav_fr), rel_col)

long_all <- amav_fr %>%
  pivot_longer(
    cols      = all_of(disease_cols),
    names_to  = "Disease",
    values_to = "Fold_increase"
  ) %>%
  mutate(
    Group = case_when(
      Disease %in% physical_diseases ~ "Physical",
      Disease %in% mental_diseases   ~ "Mental",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Group), is.finite(Fold_increase)) %>%
  mutate(
    Group = factor(Group, levels = c("Mental", "Physical"))
  )

# --------- 3. FIGURE 4A: PER-DISEASE TRAJECTORIES --------------------

p4A <- ggplot(
  long_all,
  aes(x = .data[[rel_col]],
      y = Fold_increase,
      group = Disease,
      colour = Group)
) +
  geom_line(size = 0.9, alpha = 0.85) +
  geom_point(size = 1.6, alpha = 0.85) +
  scale_y_log10() +
  scale_colour_manual(values = group_cols) +
  labs(
    title  = "Per-disease trajectories",
    x      = "Years since first datapoint",
    y      = "Fold increase",
    colour = "Group"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title   = element_text(hjust = 0, face = "bold"),
    legend.position = "none"
  )

# --------- 4. FIGURE 4B: GLOBAL GEOMETRIC MEANS (ALL DATA) -----------

geom_means_all <- long_all %>%
  group_by(Group, RelYear = .data[[rel_col]]) %>%
  summarise(
    GeomMean   = exp(mean(log(Fold_increase), na.rm = TRUE)),
    n_diseases = dplyr::n(),
    .groups    = "drop"
  )

p4B <- ggplot(
  geom_means_all,
  aes(x = RelYear, y = GeomMean, colour = Group)
) +
  geom_line(size = 1.4) +
  geom_point(size = 2) +
  scale_y_log10() +
  scale_colour_manual(values = group_cols) +
  labs(
    title  = "Global geometric means",
    x      = "Years since first datapoint",
    y      = "Fold increase (Geometric mean)",
    colour = "Group"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0, face = "bold"),
    legend.position = "right"
  )

# --------- 5. COMBINED PANEL & SAVE ----------------------------------

fig4_panel <- p4A + p4B + plot_layout(widths = c(1.1, 1))

print(fig4_panel)

# Save inside output/ directory
ggsave(
  filename = "output/FIGURE_4A-B.png",
  plot     = fig4_panel,
  width    = 14, height = 6, dpi = 300
)

ggsave("output/FIGURE_4A.png",
       p4A, width = 7, height = 6, dpi = 300)

ggsave("output/FIGURE_4B.png",
       p4B, width = 7, height = 6, dpi = 300)