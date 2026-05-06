# ============================================================================
# Figure 2F: ISR Activating Drugs — with individual replicate data points
# Journal requirement: show individual data points behind mean ± SEM
# Generates Factor1 + GDF15 plots (not EIF2AK3) for each Parental line
# ============================================================================

library(here)

# Source the original script to get all data objects
source(here("Main Figure Scripts", "Helper Scripts", "Figure_2F_Base.R"), local = TRUE)

# Re-define output directory AFTER sourcing (source script clears environment)
out_dir <- here("Results", "Figures", "Figure_2F")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# After sourcing we have:
#   merged_df  — raw replicates with cols: Timepoint, fa_vs_EIF2AK3, value, Parental, Condition, Hours, Rep
#   summ_df    — summary (mean/sem) with baseline-extended 0h for non-UT conditions
#   exprs_ISR  — full expression matrix (Timepoint x genes), includes GDF15

# ============================================================================
# Build GDF15 versions of merged_df and summ_df
# ============================================================================

# 1. Get Factor1 scores per timepoint (one row per replicate)
factor1_raw <- merged_df %>%
  filter(fa_vs_EIF2AK3 == "Factor1") %>%
  select(Timepoint, Parental, Condition, Hours, Rep, Factor1_value = value)

# 2. Get GDF15 expression per timepoint from exprs_ISR
gdf15_expr <- exprs_ISR %>%
  select(Timepoint, GDF15)

# 3. Merge and pivot into long format matching merged_df structure
merged_gdf15 <- factor1_raw %>%
  left_join(gdf15_expr, by = "Timepoint") %>%
  pivot_longer(
    cols = c("Factor1_value", "GDF15"),
    names_to = "fa_vs_GDF15",
    values_to = "value"
  ) %>%
  mutate(fa_vs_GDF15 = ifelse(fa_vs_GDF15 == "Factor1_value", "Factor1", "GDF15"))

# 4. Normalize conditions to match original
normalize_condition <- function(x) {
  case_when(
    x %in% c("UT") ~ "UT",
    x %in% c("Dim") ~ "Dim",
    x == "Ars" ~ "Ars",
    x == "Tg"  ~ "Tg",
    TRUE ~ x
  )
}
merged_gdf15 <- merged_gdf15 %>%
  mutate(Condition = normalize_condition(Condition))

conditions_wanted <- c("UT", "Dim", "Ars", "Tg")
merged_gdf15 <- merged_gdf15 %>%
  filter(Condition %in% conditions_wanted) %>%
  mutate(Condition = factor(Condition, levels = conditions_wanted))

# 5. Create summary stats for GDF15 data
summ_gdf15 <- merged_gdf15 %>%
  group_by(Parental, Condition, Hours, fa_vs_GDF15) %>%
  summarise(
    n    = dplyr::n(),
    mean = mean(value, na.rm = TRUE),
    sd   = sd(value, na.rm = TRUE),
    sem  = sd / sqrt(n),
    .groups = "drop"
  )

# 6. Extend 0h baseline: use UT 0h for all non-UT conditions
zero_ut_gdf <- summ_gdf15 %>%
  filter(Condition == "UT", Hours == 0) %>%
  select(Parental, fa_vs_GDF15, n, mean, sd, sem)

zero_rows_gdf <- summ_gdf15 %>%
  filter(Condition != "UT") %>%
  distinct(Parental, Condition, fa_vs_GDF15) %>%
  left_join(zero_ut_gdf, by = c("Parental", "fa_vs_GDF15")) %>%
  mutate(Hours = 0L)

summ_gdf15 <- bind_rows(summ_gdf15, zero_rows_gdf) %>%
  arrange(Parental, Condition, fa_vs_GDF15, Hours) %>%
  distinct(Parental, Condition, Hours, fa_vs_GDF15, .keep_all = TRUE)

# Also need baseline raw replicates for non-UT conditions at hour 0
raw_baselines_gdf <- purrr::map_dfr(
  setdiff(conditions_wanted, "UT"),
  function(cond) {
    merged_gdf15 %>%
      filter(Condition == "UT", Hours == 0) %>%
      mutate(Condition = cond)
  }
)
merged_gdf15 <- bind_rows(merged_gdf15, raw_baselines_gdf)

# ============================================================================
# Color palettes for conditions
# ============================================================================
# Factor1 plots: blue palette (matching manuscript Figure 2F)
condition_colors_factor1 <- c(
  "UT"  = "#DAEBFA",  # palest blue  (untreated)
  "Ars" = "#AFD2F3",  # light blue   (arsenite)
  "Tg"  = "#74A8DF",  # medium blue  (thapsigargin)
  "Dim" = "#476CA9"   # dark blue    (dimerizable PERK)
)

# GDF15 plots: pink palette (matching 2E GDF15 style)
condition_colors_gdf15 <- c(
  "UT"  = "#FAD0E8",  # light pink
  "Dim" = "#F48FB1",  # medium pink
  "Ars" = "#C2185B",  # dark pink
  "Tg"  = "#AD1457"   # very dark pink
)

# ============================================================================
# Plotting function
# ============================================================================
make_enhanced_plot <- function(summ_data, raw_data, parental_val, analyte_val,
                               analyte_col, title_suffix, color_palette) {
  # Filter summary data
  summ_sub <- summ_data %>%
    filter(Parental == parental_val, !!sym(analyte_col) == analyte_val) %>%
    mutate(ymin = mean - sem, ymax = mean + sem)

  # Filter raw data
  raw_sub <- raw_data %>%
    filter(Parental == parental_val, !!sym(analyte_col) == analyte_val)

  p <- ggplot() +
    # Error bars
    geom_errorbar(
      data = summ_sub,
      aes(x = Hours, ymin = ymin, ymax = ymax, color = Condition),
      width = 0.6, linewidth = 0.5
    ) +
    # Mean lines
    geom_line(
      data = summ_sub,
      aes(x = Hours, y = mean, group = Condition, color = Condition),
      linewidth = 1, alpha = 0.9
    ) +
    # Mean points — translucent fill
    geom_point(
      data = summ_sub,
      aes(x = Hours, y = mean, color = Condition),
      size = 3, alpha = 0.3
    ) +
    # Mean points — solid outline
    geom_point(
      data = summ_sub,
      aes(x = Hours, y = mean, color = Condition),
      shape = 1, size = 3, stroke = 1
    ) +
    # Individual replicates ON TOP so they are always visible
    geom_point(
      data = raw_sub,
      aes(x = Hours, y = value, color = Condition),
      size = 2.8, alpha = 0.4,
      position = position_identity()
    ) +
    scale_color_manual(values = color_palette) +
    scale_x_continuous(name = "Time (hours)", breaks = seq(0, 24, 4)) +
    ylab(analyte_val) +
    labs(title = paste0(parental_val, " — ", analyte_val, " ", title_suffix)) +
    theme_minimal(base_size = 17) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.text  = element_text(colour = "black"),
      axis.title = element_text(colour = "black"),
      plot.margin = margin(t = 5.5, r = 15, b = 5.5, l = 15)
    )

  return(p)
}

# ============================================================================
# Generate Factor1 plots (from original merged_df/summ_df)
# ============================================================================

# Also need baseline raw replicates for non-UT conditions at hour 0 in original data
raw_baselines_orig <- purrr::map_dfr(
  setdiff(conditions_wanted, "UT"),
  function(cond) {
    merged_df %>%
      filter(Condition == "UT", Hours == 0) %>%
      mutate(Condition = cond)
  }
)
merged_df_extended <- bind_rows(merged_df, raw_baselines_orig)

parentals <- unique(summ_df$Parental)

for (par in parentals) {
  # Factor1 plot
  p <- make_enhanced_plot(summ_df, merged_df_extended, par, "Factor1",
                          "fa_vs_EIF2AK3", "— Individual Replicates Shown",
                          condition_colors_factor1)
  fname <- paste0("Figure_2F_", gsub(" ", "_", par), "_Factor1_with_individual_points.png")
  ggsave(file.path(out_dir, fname), plot = p, width = 8, height = 5.5, dpi = 300)
  message("Saved: ", fname)
}

# ============================================================================
# Generate GDF15 plots (from new merged_gdf15/summ_gdf15)
# ============================================================================

for (par in parentals) {
  # GDF15 plot
  p <- make_enhanced_plot(summ_gdf15, merged_gdf15, par, "GDF15",
                          "fa_vs_GDF15", "— Individual Replicates Shown",
                          condition_colors_gdf15)
  fname <- paste0("Figure_2F_", gsub(" ", "_", par), "_GDF15_with_individual_points.png")
  ggsave(file.path(out_dir, fname), plot = p, width = 8, height = 5.5, dpi = 300)
  message("Saved: ", fname)
}

# ============================================================================
# EXPORT SOURCE DATA CSV — individual datapoints
# ============================================================================
# Factor1 replicates (from merged_df_extended which includes baseline copies)
factor1_individual <- merged_df_extended %>%
  filter(fa_vs_EIF2AK3 == "Factor1",
         Condition %in% conditions_wanted) %>%
  select(Parental, Condition, Hours, Rep, value) %>%
  mutate(Measure = "Factor1")

# GDF15 replicates (from merged_gdf15 which includes baseline copies)
gdf15_individual <- merged_gdf15 %>%
  filter(fa_vs_GDF15 == "GDF15") %>%
  select(Parental, Condition, Hours, Rep, value) %>%
  mutate(Measure = "GDF15")

# Join summary stats (mean, sem) onto individual datapoints
summ_f1_join <- summ_df %>%
  filter(fa_vs_EIF2AK3 == "Factor1") %>%
  select(Parental, Condition, Hours, mean, sem) %>%
  mutate(Measure = "Factor1")

summ_gdf_join <- summ_gdf15 %>%
  filter(fa_vs_GDF15 == "GDF15") %>%
  select(Parental, Condition, Hours, mean, sem) %>%
  mutate(Measure = "GDF15")

summ_join <- bind_rows(summ_f1_join, summ_gdf_join)

source_data_individual <- bind_rows(factor1_individual, gdf15_individual) %>%
  left_join(summ_join, by = c("Parental", "Condition", "Hours", "Measure")) %>%
  arrange(Measure, Parental, Condition, Hours, Rep)

write.csv(source_data_individual,
          file.path(out_dir, "Figure_2F_source_data_individual.csv"),
          row.names = FALSE)
message("  Saved: Figure_2F_source_data_individual.csv (",
        nrow(source_data_individual), " rows)")

message("\n=== Figure 2F journal-compliant plots saved ===")
message("  Factor1 and GDF15 plots for each Parental line in: ", out_dir)
