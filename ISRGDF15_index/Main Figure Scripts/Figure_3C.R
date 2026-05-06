################################################################################
# Figure 3C with Individual Points
# Creates multi-panel grid of all tissue scatter plots (AGE vs Factor1)
# with individual sample points colored by DTHPLCE and copies individual
# tissue plots to the journal guidelines folder
################################################################################

library(here)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(tidyr)

# Source the original Figure 3C script to get all data objects
source(here("Main Figure Scripts", "Helper Scripts", "Figure_3C_Base.R"), local = TRUE)

# After sourcing, we have:
# - all_tissue_scores: per-sample data with columns SAMPID, Factor1, AGE, SEX, GDF15, Tissue
# - results_df: per-tissue summary (when using cached data)
# We need to add DTHPLCE from the attributes

# Set up output directories
journal_out_dir <- here("Results", "Figures", "Figure_3C")
if (!dir.exists(journal_out_dir)) dir.create(journal_out_dir, recursive = TRUE)

individual_tissues_dir <- file.path(journal_out_dir, "Individual_Tissues")
if (!dir.exists(individual_tissues_dir)) dir.create(individual_tissues_dir, recursive = TRUE)

# Define color mapping for DTHPLCE
color_mapping <- c(
  "Hospital inpatient" = "orange",
  "Emergency room" = "maroon",
  "Other" = "gray"
)

################################################################################
# Merge all_tissue_scores with DTHPLCE from attributes
################################################################################

# Get unique SAMPID-DTHPLCE mapping from Atts_merged (created by sourcing the Atts_COD script)
# Atts_merged should be available from sourcing the 01_Generate_Figure.R script
# which sources Attibutes_Phenos_Merged_plus_COD.R

if (!exists("Atts_merged")) {
  # If Atts_merged isn't available, source the attributes script directly
  Atts_COD_readin <- here("Main Figure Scripts", "Helper Scripts", "gtex", "Attibutes_Phenos_Merged_plus_COD.R")
  source(Atts_COD_readin, local = TRUE)
}

# Extract SAMPID and DTHPLCE from Atts_merged
dthplce_mapping <- Atts_merged %>%
  select(SAMPID, DTHPLCE) %>%
  distinct() %>%
  mutate(DTHPLCE = replace_na(DTHPLCE, "Other"))

# Merge with all_tissue_scores
all_tissue_scores <- all_tissue_scores %>%
  left_join(dthplce_mapping, by = "SAMPID") %>%
  mutate(
    DTHPLCE = ifelse(is.na(DTHPLCE), "Other", DTHPLCE),
    # Ensure DTHPLCE levels are consistent
    DTHPLCE = case_when(
      DTHPLCE %in% c("Hospital inpatient", "Emergency room") ~ DTHPLCE,
      TRUE ~ "Other"
    )
  )

# Rename Factor1 to chosen_Factor for consistency
boxplot_data <- all_tissue_scores %>%
  rename(chosen_Factor = Factor1)

# Ensure DTHPLCE is a factor with correct levels
boxplot_data$DTHPLCE <- factor(boxplot_data$DTHPLCE,
                                levels = c("Hospital inpatient", "Emergency room", "Other"),
                                exclude = NULL)

# Convert columns to numeric
boxplot_data <- boxplot_data %>%
  mutate(
    AGE = as.numeric(AGE),
    chosen_Factor = as.numeric(as.character(chosen_Factor)),
    GDF15 = as.numeric(as.character(GDF15))
  )

# Remove any rows with NA in key columns
boxplot_data_clean <- boxplot_data %>%
  filter(!is.na(AGE), !is.na(chosen_Factor), !is.na(DTHPLCE))

################################################################################
# PART 1: Create multi-panel faceted grid of all tissues
################################################################################

# Get unique tissues and ensure they're ordered
tissues <- unique(boxplot_data_clean$Tissue)
n_tissues <- length(tissues)

cat("Creating multi-panel figure for", n_tissues, "tissues\n")

# Create the faceted plot
faceted_plot <- ggplot(boxplot_data_clean, aes(x = AGE, y = chosen_Factor, color = DTHPLCE)) +
  geom_point(alpha = 0.3, size = 1.5) +  # Individual sample points
  geom_smooth(method = "lm", se = TRUE, color = "red", fill = "red", alpha = 0.1) +
  facet_wrap(~ Tissue, scales = "free_y", ncol = 6) +
  scale_color_manual(
    name = "Place of Death",
    values = color_mapping,
    na.value = "gray"
  ) +
  labs(
    title = "Figure 3C: AGE vs Factor1 — Individual Samples per Tissue",
    x = "Age (years)",
    y = "Index Score (Factor1)",
    caption = "Individual sample points colored by place of death. Red line = linear regression fit with 95% CI."
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, hjust = 0.5, margin = margin(b = 10)),
    plot.caption = element_text(size = 10, hjust = 0, vjust = 1),
    strip.text = element_text(size = 9, face = "bold"),
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 8),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    legend.position = "bottom",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  )

# Save the faceted plot as a large PNG
ggsave(
  filename = "Figure_3C_All_Tissues_Faceted.png",
  plot = faceted_plot,
  path = journal_out_dir,
  width = 24,
  height = 30,
  dpi = 300,
  units = "in"
)

cat("Saved faceted plot: Figure_3C_All_Tissues_Faceted.png\n")

################################################################################
# PART 2: Save individual tissue scatter plots to subfolder
################################################################################

cat("Saving individual tissue plots to subfolder...\n")

# Additional output folder with larger fonts
individual_spearman_dir <- file.path(journal_out_dir, "Individual Tissue Spearman")
if (!dir.exists(individual_spearman_dir)) dir.create(individual_spearman_dir, recursive = TRUE)

# Create a tissue-level summary for annotations (if available)
tissue_summary <- NULL
if (exists("results_df")) {
  # Filter for Factor1 comparisons only
  tissue_summary <- results_df %>%
    filter(Comparison == "chosen_Factor") %>%
    select(Tissue, Spearman_Rho, Asterisks)
}

for (tissue in tissues) {
  # Filter data for this tissue
  tissue_data <- boxplot_data_clean %>%
    filter(Tissue == tissue)

  # Get correlation stats for this tissue (if available)
  rho_val <- NA
  asterisks <- ""

  if (!is.null(tissue_summary)) {
    tissue_stats <- tissue_summary %>%
      filter(Tissue == tissue)

    if (nrow(tissue_stats) > 0) {
      rho_val <- tissue_stats$Spearman_Rho[1]
      asterisks <- tissue_stats$Asterisks[1]
    }
  }

  # Title is just the tissue name (no rho/p); stats are overlaid on the plot
  title_text <- tissue

  # Build annotation label for Spearman stats
  spearman_label <- if (!is.na(rho_val)) {
    # Compute p-value for this tissue
    spearman_test <- cor.test(tissue_data$AGE, tissue_data$chosen_Factor, method = "spearman", exact = FALSE)
    p_val <- spearman_test$p.value
    paste0("Spearman rho: ", round(rho_val, 4), "\n",
           "p-value: ", formatC(p_val, format = "e", digits = 2))
  } else {
    NULL
  }

  tissue_filename <- gsub(" ", "_", tissue)

  # --- Original size version (saved to Individual_Tissues) ---
  individual_plot <- ggplot(tissue_data, aes(x = AGE, y = chosen_Factor)) +
    geom_point(alpha = 0.3, size = 1.5, color = "gray50") +
    geom_smooth(method = "lm", se = TRUE, color = "red", fill = "red", alpha = 0.1) +
    labs(
      title = title_text,
      x = "Age (years)",
      y = "Index Score (Factor1)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 12, hjust = 0.5),
      axis.title = element_text(size = 11),
      axis.text = element_text(size = 10)
    )

  # Add Spearman annotation to the plot (upper-left)
  if (!is.null(spearman_label)) {
    individual_plot <- individual_plot +
      annotate("text", x = -Inf, y = Inf,
               label = spearman_label,
               hjust = -0.05, vjust = 1.2, size = 3.5,
               color = "black", fontface = "italic")
  }

  ggsave(
    filename = paste0("Figure_3C_", tissue_filename, ".png"),
    plot = individual_plot,
    path = individual_tissues_dir,
    width = 8,
    height = 6,
    dpi = 300,
    units = "in"
  )

  # --- Large font version (3x font sizes, gray points, saved to Individual Tissue Spearman) ---
  individual_plot_lg <- ggplot(tissue_data, aes(x = AGE, y = chosen_Factor)) +
    geom_point(alpha = 0.3, size = 3, color = "gray50") +
    geom_smooth(method = "lm", se = TRUE, color = "red", fill = "red", alpha = 0.1) +
    labs(
      title = title_text,
      x = "Age (years)",
      y = "Index Score (Factor1)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 36, hjust = 0.5),
      axis.title = element_text(size = 33),
      axis.text = element_text(size = 30)
    )

  # Add Spearman annotation (larger font for the large version)
  if (!is.null(spearman_label)) {
    individual_plot_lg <- individual_plot_lg +
      annotate("text", x = -Inf, y = Inf,
               label = spearman_label,
               hjust = -0.05, vjust = 1.2, size = 9,
               color = "black", fontface = "italic")
  }

  ggsave(
    filename = paste0("Figure_3C_", tissue_filename, ".png"),
    plot = individual_plot_lg,
    path = individual_spearman_dir,
    width = 12,
    height = 9,
    dpi = 300,
    units = "in"
  )
}

cat("Saved", length(tissues), "individual tissue plots to:", individual_tissues_dir, "\n")
cat("Saved", length(tissues), "large-font tissue plots to:", individual_spearman_dir, "\n")

################################################################################
# PART 3: Summary information
################################################################################

cat("\n========== Figure 3C Summary ==========\n")
cat("Total tissues analyzed:", n_tissues, "\n")
cat("Total individual samples:", nrow(boxplot_data_clean), "\n")

if (exists("results_df")) {
  cat("\nPer-tissue summary statistics (Factor1 vs AGE):\n")
  results_chosen_factor <- results_df %>%
    filter(Comparison == "chosen_Factor")
  print(results_chosen_factor)
}

cat("\n========== Output Files ==========\n")
cat("Multi-panel faceted figure: ", file.path(journal_out_dir, "Figure_3C_All_Tissues_Faceted.png"), "\n")
cat("Individual tissue plots folder: ", individual_tissues_dir, "\n")

################################################################################
