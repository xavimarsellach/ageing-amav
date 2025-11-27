#| label: FIGURE_S1B.R
#| echo: true
#| eval: true
#| warning: false
#| message: false

# FIGURE_S1B.R
# Control analysis: remove early monotonic stretches
# Geometric means (Mental vs Physical) + ANCOVA for k = 0..10

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(writexl)
library(forcats)

# ---------------- CONFIG ----------------

if (!dir.exists("output")) {
  dir.create("output")
}

file_path    <- "output/AMAV_DATA-no-monotonic_AUTO.xlsx"  # adjust path if needed
sheet_name   <- "FOLD_RELATIVE"
min_rel_year <- 0                  # set to 10 for RelYear >= 10

thresholds   <- 0:10

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

# ---------------- READ & LONG FORMAT ----------------

df <- read_excel(file_path, sheet = sheet_name)

rel_col      <- "RelYear"
disease_cols <- setdiff(names(df), rel_col)

long_df <- df %>%
  pivot_longer(
    cols      = all_of(disease_cols),
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
  filter(.data[[rel_col]] >= min_rel_year)

# ---------------- GEOMETRIC MEANS (BASE) ----------------

geom_means_all <- long_df %>%
  group_by(GROUP, !!sym(rel_col)) %>%
  summarise(
    GeomMean   = exp(mean(log(Fold_increase), na.rm = TRUE)),
    n_diseases = sum(is.finite(Fold_increase)),
    .groups    = "drop"
  ) %>%
  rename(RelYear = !!sym(rel_col))

# ---------------- THRESHOLD LOOP ----------------

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
      p_group        = NA_real_
    )
  } else {
    data_k <- gm_k %>% mutate(logGeom = log10(GeomMean))
    
    model_int    <- lm(logGeom ~ RelYear * GROUP, data = data_k)
    model_no_int <- lm(logGeom ~ RelYear + GROUP, data = data_k)
    
    aov_int <- anova(model_no_int, model_int)
    p_int   <- aov_int$`Pr(>F)`[2]
    
    sm       <- summary(model_no_int)
    coef_tab <- sm$coefficients
    p_group  <- coef_tab["GROUPPhysical", "Pr(>|t|)"]
    
    ancova_results[[as.character(k)]] <- tibble(
      k_threshold    = k,
      n_years_shared = n_years,
      n_points       = nrow(gm_k),
      p_interaction  = p_int,
      p_group        = p_group
    )
  }
}

geom_all_thresholds <- bind_rows(all_geom_list)
ancova_table        <- bind_rows(ancova_results)

# ---------------- ORDER FACETS (0..10) ----------------

geom_all_thresholds <- geom_all_thresholds %>%
  mutate(
    k_factor = factor(k_threshold, levels = thresholds)
  )

ancova_table <- ancova_table %>%
  mutate(
    k_factor = factor(k_threshold, levels = thresholds)
  )

labels_df <- ancova_table %>%
  group_by(k_factor) %>%
  summarise(
    p_int = unique(p_interaction),
    label = ifelse(
      is.na(p_int),
      paste0("k = ", as.character(k_factor), " (p = NA)"),
      paste0(
        "k = ", as.character(k_factor),
        " (p = ", formatC(p_int, format = "f", digits = 4), ")"
      )
    ),
    .groups = "drop"
  )

# ---------------- PLOT ----------------

title_prefix <- if (min_rel_year == 0) {
  "Global Fold Increase: FULL DATA (no monotonic early stretches)"
} else {
  paste0("Global Fold Increase (no monotonic early stretches): RelYear \u2265 ", min_rel_year)
}

p_facets <- ggplot(
  geom_all_thresholds,
  aes(x = RelYear, y = GeomMean, colour = GROUP)
) +
  geom_line(size = 0.9) +
  geom_point(size = 1.3) +
  scale_y_log10() +
  facet_wrap(
    ~ k_factor,
    ncol = 3,
    labeller = labeller(
      k_factor = function(z) {
        labels_df$label[match(z, labels_df$k_factor)]
      }
    )
  ) +
  labs(
    title  = paste0(title_prefix, "\nGeometric Means by Threshold"),
    x      = "Years since first datapoint (RelYear)",
    y      = "Fold Increase (Geometric Mean)",
    colour = "GROUP"
  ) +
  scale_colour_manual(values = c(Mental = "red", Physical = "blue")) +
  theme_bw(base_size = 11) +
  theme(
    plot.title      = element_text(hjust = 0.5, face = "bold"),
    legend.position = "bottom"
  )

print(p_facets)

ggsave(
  filename = if (min_rel_year == 0)
    "output/FIGURE_S1B.png"
  else
    paste0("output/FIGURE_S1B_RelYear", min_rel_year, "_k0_to_k10_ordered.png"),
  plot  = p_facets,
  width = 12,
  height = 10,
  dpi   = 300
)

write_xlsx(
  list(
    GeomMeans_by_threshold = geom_all_thresholds,
    ANCOVA_results         = ancova_table
  ),
  if (min_rel_year == 0)
    "output/GEOMETRIC_MEANS_thresholds_FULL_no_monotonic_ordered.xlsx"
  else
    paste0("output/GEOMETRIC_MEANS_thresholds_RelYear", min_rel_year, "_no_monotonic_ordered.xlsx")
)