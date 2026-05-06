# ============================================================================
# Figure 4B: Factor1/GDF15 vs Place of Death - DATA PROCESSING AND FIGURE
# ============================================================================
# This script creates plots showing:
#   - Factor1 (ISR) score vs place of death for brain frontal cortex and heart left ventricle
#   - GDF15 expression vs place of death for the same tissues
#
# Output:
#   - 4 plots (2 tissues x 2 measures)
#   - CSV files for each plot for Prism import
# ============================================================================

library(here)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(effsize)   # cohen.d() with hedges.correction

# ============================================================================
# CONFIGURATION
# ============================================================================
# Tissues to include in Figure 4B (matching manuscript)
target_tissues <- c("brain_frontal_cortex", "heart_left_ventricle")

# Place of death categories
control <- "Hospital inpatient"
condition <- "Emergency room"

# Output directories
output_dir <- here("Results", "Figures", "Figure_4B")
csv_output_dir <- file.path(output_dir, "CSV_for_Prism")

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
if (!dir.exists(csv_output_dir)) dir.create(csv_output_dir, recursive = TRUE)

# ============================================================================
# LOAD DATA
# ============================================================================
# Build the per-donor / per-tissue long-format dataframe by reading the
# per-tissue Plot_Data_<tissue>.csv files written by
# Figure_3B_4A_Tissue_Scatter.R. The helper returns columns
#   SAMPID, Tissue, DTHPLCE, fa_vs_gdf15, value, group, ...
# with NO DTHPLCE filtering. The Hospital-inpatient / Emergency-room
# restriction for Figure 4B is applied below via the case_when on `group`.
# ============================================================================
message("Loading GTEx place of death data...")

source(here("Main Figure Scripts", "Helper Scripts",
            "Load_All_Tissue_DTHPLCE.R"))
data <- load_all_tissue_dthplce()

# Restrict to the two place-of-death categories Figure 4B compares.
# (Helper passes through every DTHPLCE category; small ones like "At MVA
# scene" have N < 20 and would cause the all(N >= 20) check below to drop
# every tissue if left in.)
data <- data %>%
  filter(DTHPLCE %in% c("Hospital inpatient", "Emergency room"))

message("  Data loaded: ", nrow(data), " rows")

# ============================================================================
# FILTER AND PROCESS DATA
# ============================================================================
message("\nProcessing data for target tissues...")

# Check what tissues are available
available_tissues <- unique(data$Tissue)
message("  Available tissues: ", length(available_tissues))

# Find matching tissue names (handle case differences)
tissue_matches <- sapply(target_tissues, function(t) {
  matches <- grep(gsub("_", ".", t), available_tissues, ignore.case = TRUE, value = TRUE)
  if (length(matches) == 0) {
    matches <- grep(gsub("_", " ", t), available_tissues, ignore.case = TRUE, value = TRUE)
  }
  if (length(matches) == 0) matches <- NA
  return(matches[1])
})

message("  Tissue mapping:")
for (i in seq_along(target_tissues)) {
  message("    ", target_tissues[i], " -> ", tissue_matches[i])
}

# Filter for target tissues
data_filtered <- data %>%
  filter(Tissue %in% na.omit(tissue_matches))

# Add sample counts per group
data_filtered <- data_filtered %>%
  group_by(Tissue, group) %>%
  mutate(N = n()) %>%
  ungroup()

# Filter for groups with N >= 20
data_filtered <- data_filtered %>%
  group_by(Tissue) %>%
  filter(all(N >= 20)) %>%
  ungroup()

message("  Samples after filtering: ", nrow(data_filtered))

# ============================================================================
# SEPARATE FACTOR1 AND GDF15 DATA
# ============================================================================
message("\nSeparating Factor1 and GDF15 data...")

# Parse the group column to extract measure (Factor1/GDF15) and death place
data_filtered <- data_filtered %>%
  mutate(
    death_place = case_when(
      grepl("Hospital inpatient", group) ~ "Hospital inpatient",
      grepl("Emergency room", group) ~ "Emergency room",
      TRUE ~ NA_character_
    ),
    measure = case_when(
      grepl("Factor1", group) ~ "Factor1",
      grepl("GDF15", group) ~ "GDF15",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(death_place) & !is.na(measure))

# Split into Factor1 and GDF15 datasets
factor1_data <- data_filtered %>% filter(measure == "Factor1")
gdf15_data <- data_filtered %>% filter(measure == "GDF15")

message("  Factor1 samples: ", nrow(factor1_data))
message("  GDF15 samples: ", nrow(gdf15_data))

# ============================================================================
# CREATE PLOTS AND SAVE CSVs
# ============================================================================
message("\nGenerating plots and CSVs...")

# Color palette
death_colors <- c("Hospital inpatient" = "#E69F00", "Emergency room" = "#56B4E9")

# Function to create plot and CSV for one tissue/measure combination
create_plot_and_csv <- function(plot_data, tissue_name, measure_name) {
  tissue_clean <- gsub("_", " ", tissue_name)
  tissue_clean <- tools::toTitleCase(tissue_clean)

  # Calculate statistics
  stats <- plot_data %>%
    group_by(death_place) %>%
    summarise(
      n = n(),
      mean = mean(value, na.rm = TRUE),
      median = median(value, na.rm = TRUE),
      sd = sd(value, na.rm = TRUE),
      .groups = "drop"
    )

  # Wilcoxon test
  wilcox_result <- wilcox.test(value ~ death_place, data = plot_data)
  p_value <- wilcox_result$p.value
  p_label <- ifelse(p_value < 0.001, "***",
                    ifelse(p_value < 0.01, "**",
                           ifelse(p_value < 0.05, "*", "ns")))

  # Hedges' g for the H vs E comparison.
  # Factor levels are set explicitly so the sign convention matches Figure 4C:
  # Hospital inpatient first, Emergency room second  =>  g > 0 means HI > ER.
  death_factor <- factor(plot_data$death_place,
                         levels = c("Hospital inpatient", "Emergency room"))
  hedges <- tryCatch(
    effsize::cohen.d(plot_data$value, death_factor,
                     hedges.correction = TRUE),
    error = function(e) NULL
  )
  hedges_g <- if (!is.null(hedges)) as.numeric(hedges$estimate) else NA_real_
  comparison_str <- paste0("Hospital inpatient_", measure_name,
                           " vs Emergency room_", measure_name)
  n_hi <- sum(plot_data$death_place == "Hospital inpatient", na.rm = TRUE)
  n_er <- sum(plot_data$death_place == "Emergency room",     na.rm = TRUE)

  hedges_row <- data.frame(
    Tissue               = tissue_name,
    Measure              = measure_name,
    Comparison           = comparison_str,
    Hedges_G             = hedges_g,
    Wilcoxon_P           = p_value,
    N_Hospital_inpatient = n_hi,
    N_Emergency_room     = n_er,
    stringsAsFactors     = FALSE
  )

  # Y-axis label
  y_label <- ifelse(measure_name == "Factor1", "ISR Score (Factor1)", "GDF15 Expression")

  # Create plot
  p <- ggplot(plot_data, aes(x = death_place, y = value, fill = death_place)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.5, size = 1.5) +
    scale_fill_manual(values = death_colors) +
    stat_compare_means(
      method = "wilcox.test",
      label = "p.signif",
      comparisons = list(c("Hospital inpatient", "Emergency room")),
      label.y = max(plot_data$value, na.rm = TRUE) * 1.1
    ) +
    labs(
      title = paste0(measure_name, " vs Place of Death - ", tissue_clean),
      x = "Place of Death",
      y = y_label,
      caption = paste0("Wilcoxon p = ", format(p_value, digits = 3))
    ) +
    theme_minimal() +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 12, face = "bold"),
      axis.text.x = element_text(size = 10),
      axis.title = element_text(size = 11)
    )

  # Save CSV for Prism (wide format)
  csv_wide <- plot_data %>%
    select(death_place, value) %>%
    group_by(death_place) %>%
    mutate(row_id = row_number()) %>%
    pivot_wider(names_from = death_place, values_from = value) %>%
    select(-row_id)

  csv_filename <- paste0("Figure_4B_", measure_name, "_", gsub(" ", "_", tissue_name), ".csv")
  write.csv(csv_wide, file.path(csv_output_dir, csv_filename), row.names = FALSE, na = "")

  # Also save stats
  stats_filename <- paste0("Figure_4B_", measure_name, "_", gsub(" ", "_", tissue_name), "_stats.csv")
  write.csv(stats, file.path(csv_output_dir, stats_filename), row.names = FALSE)

  return(list(plot = p, csv_file = csv_filename, stats = stats,
              hedges = hedges_row))
}

# Generate all 4 plots
plots <- list()
hedges_results <- list()
for (tissue in unique(data_filtered$Tissue)) {
  for (measure in c("Factor1", "GDF15")) {
    plot_data <- data_filtered %>%
      filter(Tissue == tissue, measure == !!measure)

    if (nrow(plot_data) > 0) {
      result <- create_plot_and_csv(plot_data, tissue, measure)
      plot_key <- paste0(tissue, "_", measure)
      plots[[plot_key]]          <- result$plot
      hedges_results[[plot_key]] <- result$hedges
      message("  Created: ", result$csv_file)
    }
  }
}

# ============================================================================
# SAVE HEDGES G RESULTS
# ============================================================================
# One row per (tissue, measure) — Hedges g for the Hospital inpatient vs
# Emergency room comparison. Sign matches Figure 4C convention:
# positive g  =>  Hospital inpatient > Emergency room.
hedges_df <- bind_rows(hedges_results)
hedges_csv_file <- file.path(output_dir, "Figure_4B_Hedges_G_results.csv")
write.csv(hedges_df, hedges_csv_file, row.names = FALSE)
message("  Saved: Figure_4B_Hedges_G_results.csv (",
        nrow(hedges_df), " comparisons)")

# ============================================================================
# SAVE COMBINED FIGURE
# ============================================================================
message("\nSaving combined figure...")

if (length(plots) >= 4) {
  combined_plot <- ggarrange(
    plotlist = plots,
    ncol = 2, nrow = 2,
    labels = c("A", "B", "C", "D")
  )

  ggsave(
    file.path(output_dir, "Figure_4B_Place_of_Death_Combined.png"),
    combined_plot,
    width = 12, height = 10, dpi = 300
  )
  message("  Saved: Figure_4B_Place_of_Death_Combined.png")

  ggsave(
    file.path(output_dir, "Figure_4B_Place_of_Death_Combined.pdf"),
    combined_plot,
    width = 12, height = 10
  )
  message("  Saved: Figure_4B_Place_of_Death_Combined.pdf")
}

# Save individual plots
for (name in names(plots)) {
  ggsave(
    file.path(output_dir, paste0("Figure_4B_", name, ".png")),
    plots[[name]],
    width = 6, height = 5, dpi = 300
  )
}
message("  Saved individual plots")

# ============================================================================
# COPY TO SHARED OUTPUT FOLDER
# ============================================================================
shared_out_dir <- here("Results", "Figures")
if (!dir.exists(shared_out_dir)) dir.create(shared_out_dir, recursive = TRUE)
if (file.exists(file.path(output_dir, "Figure_4B_Place_of_Death_Combined.png"))) {
  file.copy(
    file.path(output_dir, "Figure_4B_Place_of_Death_Combined.png"),
    file.path(shared_out_dir, "Figure_4B.png"),
    overwrite = TRUE
  )
  message("\nCopied to shared Output folder: ", shared_out_dir)
}

message("\n============ COMPLETE ============")
message("CSVs for Prism saved to: ", csv_output_dir)

# ============================================================================
# SAVE SOURCE DATA FOR FIGURE 4B
# ============================================================================
message("\nExporting source data for Figure 4B...")

# Prepare source data containing only the plotted data (filtered by target tissues)
source_data <- data_filtered %>%
  select(Tissue, death_place, value, measure)

# Save to CSV
csv_source_file <- file.path(output_dir, "Figure_4B_source_data.csv")
write.csv(source_data, csv_source_file, row.names = FALSE)

message("  Saved: Figure_4B_source_data.csv")
message("  Location: ", output_dir)
