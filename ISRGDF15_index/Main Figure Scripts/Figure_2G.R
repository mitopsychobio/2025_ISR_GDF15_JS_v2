# ============================================================================
# Figure 2G: DmrPERK vs ATF4KO — with individual replicate data points
# Journal requirement: show individual data points behind mean ± SEM
# ============================================================================

library(here)

# Source the original script to get all data objects
source(here("Main Figure Scripts", "Helper Scripts", "Figure_2G_Base.R"), local = TRUE)

# Re-define output directory AFTER sourcing (source script clears environment)
out_dir <- here("Results", "Figures", "Figure_2G")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# merged_df has raw replicates: Hours, value, Condition, Kinase, fa_vs_GDF15, Rep
# summ_df has summarized: Hours, mean, sem, Condition, Kinase, fa_vs_GDF15

# Colors per analyte (matching original)
analyte_colors <- c("Factor1" = "#0b3d91", "GDF15" = "#cc6699")

# Function to create enhanced plot for a given Kinase/Analyte combination
make_enhanced_plot <- function(kinase_val, analyte_val) {
  # Filter summary data
  summ_sub <- summ_df %>%
    filter(Kinase == kinase_val, fa_vs_GDF15 == analyte_val) %>%
    mutate(ymin = mean - sem, ymax = mean + sem)

  # Filter raw data
  raw_sub <- merged_df %>%
    filter(Kinase == kinase_val, fa_vs_GDF15 == analyte_val)

  base_color <- analyte_colors[[analyte_val]]

  # Split summary and raw data by condition for layered styling
  summ_wt  <- summ_sub %>% filter(Condition == "WT")
  summ_ko  <- summ_sub %>% filter(Condition == "ATF4KO")
  raw_wt   <- raw_sub %>% filter(Condition == "WT")
  raw_ko   <- raw_sub %>% filter(Condition == "ATF4KO")

  p <- ggplot() +
    # Error bars
    geom_errorbar(
      data = summ_sub,
      aes(x = Hours, ymin = ymin, ymax = ymax),
      color = base_color,
      width = 0.6, linewidth = 0.5
    ) +
    # Mean lines (solid for WT, dashed for ATF4KO)
    geom_line(
      data = summ_sub,
      aes(x = Hours, y = mean, group = Condition, linetype = Condition),
      color = base_color,
      linewidth = 1, alpha = 0.9
    ) +
    # WT mean points — filled circle with translucent fill + solid outline
    geom_point(
      data = summ_wt,
      aes(x = Hours, y = mean),
      shape = 21,
      fill = scales::alpha(base_color, 0.3),
      color = base_color,
      size = 3, stroke = 1
    ) +
    # ATF4KO mean points — open circle (ring) to match dashed line
    geom_point(
      data = summ_ko,
      aes(x = Hours, y = mean),
      shape = 1,
      color = base_color,
      size = 3, stroke = 1
    ) +
    # WT individual replicates — solid dots ON TOP
    geom_point(
      data = raw_wt,
      aes(x = Hours, y = value),
      color = base_color, shape = 16,
      size = 2.8, alpha = 0.4,
      position = position_identity()
    ) +
    # ATF4KO individual replicates — open circles ON TOP
    geom_point(
      data = raw_ko,
      aes(x = Hours, y = value),
      color = base_color, shape = 1,
      size = 2.8, alpha = 0.4,
      position = position_identity()
    ) +
    scale_linetype_manual(values = c("WT" = "solid", "ATF4KO" = "dashed")) +
    scale_x_continuous(name = "Time (hours)", breaks = seq(0, 24, 4)) +
    ylab(analyte_val) +
    labs(title = paste0(kinase_val, " — ", analyte_val, " — Individual Replicates Shown")) +
    theme_minimal(base_size = 17) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.text  = element_text(colour = "black"),
      axis.title = element_text(colour = "black"),
      plot.margin = margin(t = 5.5, r = 15, b = 5.5, l = 15)
    )

  return(p)
}

# Generate plots for all Kinase/Analyte combinations
kinases  <- unique(summ_df$Kinase)
analytes <- unique(summ_df$fa_vs_GDF15)

for (kin in kinases) {
  for (ana in analytes) {
    p <- make_enhanced_plot(kin, ana)
    fname <- paste0("Figure_2G_", gsub(" ", "_", kin), "_", ana, "_with_individual_points.png")
    ggsave(file.path(out_dir, fname), plot = p, width = 8, height = 5.5, dpi = 300)
    message("Saved: ", fname)
  }
}

# ============================================================================
# EXPORT SOURCE DATA CSV — individual datapoints
# ============================================================================
# Join summary stats (mean, sem) onto individual datapoints
summ_join <- summ_df %>%
  select(Condition, Kinase, Hours, fa_vs_GDF15, mean, sem) %>%
  rename(Measure = fa_vs_GDF15)

source_data_individual <- merged_df %>%
  rename(Measure = fa_vs_GDF15) %>%
  select(Condition, Kinase, Hours, Rep, Measure, value) %>%
  left_join(summ_join, by = c("Condition", "Kinase", "Hours", "Measure")) %>%
  arrange(Measure, Kinase, Condition, Hours, Rep)

write.csv(source_data_individual,
          file.path(out_dir, "Figure_2G_source_data_individual.csv"),
          row.names = FALSE)
message("  Saved: Figure_2G_source_data_individual.csv (",
        nrow(source_data_individual), " rows)")

message("\n=== Figure 2G journal-compliant plots saved ===")
