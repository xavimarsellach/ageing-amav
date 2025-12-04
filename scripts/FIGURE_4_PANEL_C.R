#| label: FIGURE_4_PANEL_C.R
#| echo: true
#| eval: true
#| warning: false
#| message: false

# FIGURE_4_PANEL_C.R
# Geometric means of fold increases (Mental vs Physical)
# + progressive thresholds k = 0..10
# + ANCOVA per threshold
# Uses AMAV_DATA.xlsx, sheet "FOLD_RELATIVE"

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(writexl)

# ------------------ 0. Ensure output/ exists ----------------------------

if (!dir.exists("output")) {
  dir.create("output")
}

# ------------------ 1. Read FOLD_RELATIVE -----------------------------------

file_path <- "data/AMAV_DATA.xlsx"  # same as for panels A–B
df <- read_excel(file_path, sheet = "FOLD_RELATIVE")

rel_col <- "RelYear"
disease_cols <- setdiff(names(df), rel_col)

# ------------------ 2. Mental vs Physical sets --------------------------

physical <- c(
  "Allergic Rhinitis", "Alopecia Aereata", "Anaemia", "Asthma", "Atopic Dermatitis",
  "Cancer", "Coeliac Disease", "Food Allergies", "Hypertension", "Lupus",
  "Myopia", "Obesity", "Osteoarthritis", "Rheumatoid Arthritis",
  "Type 1 Diabetes", "Type 2 Diabetes"
)

mental <- c(
  "ADHD", "Alzheimer", "Anxiety", "Autism", "Bipolar Disorder",
  "Dementia", "Depression", "Eating Disorders", "Handedness",
  "Homosexuality (Men)", "Homosexuality (Women)",
  "Personality Disorders", "Schizophrenia",
  "Suicide", "Transgender & Gender Dysphoria"
)

group_cols <- c(Mental = "red", Physical = "blue")

# ------------------ 3. Long format --------------------------------------

long_df <- df %>%
  pivot_longer(
    cols = all_of(disease_cols),
    names_to  = "Disease",
    values_to = "Fold_increase"
  ) %>%
  mutate(
    GROUP = case_when(
      Disease %in% physical ~ "Physical",
      Disease %in% mental   ~ "Mental",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(GROUP)) %>%
  mutate(GROUP = factor(GROUP, levels = c("Mental", "Physical")))

# ------------------ 4. Geometric means (all data) ------------------------

geom_means_all <- long_df %>%
  group_by(GROUP, !!sym(rel_col)) %>%
  summarise(
    GeomMean    = exp(mean(log(Fold_increase), na.rm = TRUE)),
    n_diseases  = sum(is.finite(Fold_increase)),
    .groups     = "drop"
  ) %>%
  rename(RelYear = !!sym(rel_col))

# ------------------ 5. Threshold loop k = 0..10 --------------------------

thresholds <- 0:10
all_geom_list  <- list()
ancova_results <- list()

for (k in thresholds) {
  
  gm_tmp <- geom_means_all
  
  if (k > 0) {
    gm_tmp <- gm_tmp %>%
      group_by(GROUP) %>%
      filter(n_diseases >= k) %>%
      ungroup()
  }
  
  years_keep <- gm_tmp %>%
    count(RelYear, GROUP) %>%
    count(RelYear, name = "n_groups") %>%
    filter(n_groups == 2) %>%
    pull(RelYear)
  
  gm_k <- gm_tmp %>%
    filter(RelYear %in% years_keep) %>%
    arrange(GROUP, RelYear) %>%
    mutate(k_threshold = k)
  
  all_geom_list[[as.character(k)]] <- gm_k
  
  n_years <- length(unique(gm_k$RelYear))
  
  if (n_years < 3) {
    
    ancova_results[[as.character(k)]] <- tibble(
      k_threshold    = k,
      n_years_shared = n_years,
      n_points       = nrow(gm_k),
      p_interaction  = NA_real_,
      p_group        = NA_real_,
      slope_common   = NA_real_,
      intercept      = NA_real_,
      R2_adj         = NA_real_
    )
    
  } else {
    
    data_k <- gm_k %>% mutate(logGeom = log10(GeomMean))
    
    model_int    <- lm(logGeom ~ RelYear * GROUP, data = data_k)
    model_no_int <- lm(logGeom ~ RelYear + GROUP, data = data_k)
    
    aov_int <- anova(model_no_int, model_int)
    p_int   <- aov_int$`Pr(>F)`[2]
    
    sm <- summary(model_no_int)
    coef_tab <- sm$coefficients
    
    ancova_results[[as.character(k)]] <- tibble(
      k_threshold    = k,
      n_years_shared = n_years,
      n_points       = nrow(gm_k),
      p_interaction  = p_int,
      p_group        = coef_tab["GROUPPhysical", "Pr(>|t|)"],
      slope_common   = coef_tab["RelYear", "Estimate"],
      intercept      = coef_tab["(Intercept)", "Estimate"],
      R2_adj         = sm$adj.r.squared
    )
  }
}

geom_all_thresholds <- bind_rows(all_geom_list)
ancova_table        <- bind_rows(ancova_results)

# ------------------ 6. Export to Excel (in output/) ---------------------

write_xlsx(
  list(
    GeomMeans_by_threshold = geom_all_thresholds,
    ANCOVA_results         = ancova_table
  ),
  "output/GEOMETRIC_MEANS_MENTAL_vs_PHYSICAL_ANCOVA_thresholds_k0_to_k10.xlsx"
)

# ------------------ 7. Facet labels with p-values -----------------------

facet_labels_df <- ancova_table %>%
  mutate(
    p_lab = case_when(
      is.na(p_interaction) ~ paste0("k = ", k_threshold, " (no ANCOVA)"),
      TRUE ~ paste0(
        "k = ", k_threshold,
        " (p = ", sprintf("%.6f", p_interaction), ")"
      )
    )
  )

facet_label_vec <- setNames(facet_labels_df$p_lab,
                            facet_labels_df$k_threshold)

# ------------------ 8. Plot facets (Mental vs Physical) -----------------

p_facets <- ggplot(
  geom_all_thresholds,
  aes(x = RelYear, y = GeomMean, colour = GROUP)
) +
  geom_line(size = 0.9) +
  geom_point(size = 1.3) +
  scale_y_log10() +
  scale_colour_manual(values = group_cols) +
  facet_wrap(
    ~ k_threshold,
    ncol = 3,
    labeller = labeller(k_threshold = function(z) facet_label_vec[as.character(z)])
  ) +
  labs(
    title = "Global Fold Increase: Mental vs Physical\nGeometric Means by Threshold (common years only)",
    x     = "Years since first datapoint",
    y     = "Fold Increase (Geometric Mean)",
    colour = "Group"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title      = element_text(hjust = 0.5, face = "bold"),
    legend.position = "bottom"
  )

# Mostrar la figura al Quarto
print(p_facets)

# Guardar figura a output/
ggsave(
  filename = "output/FIGURE_4C.png",
  plot     = p_facets,
  width    = 12,
  height   = 10,
  dpi      = 300
)