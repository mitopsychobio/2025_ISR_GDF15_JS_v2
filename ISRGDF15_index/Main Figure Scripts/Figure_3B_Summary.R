################################################################################
# Figure 3B with Individual Points
# Creates multi-panel grid of all tissue scatter plots (Factor1 vs GDF15)
# with individual sample points colored by DTHPLCE and copies individual
# tissue plots to the journal guidelines folder
################################################################################

library(here)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(tidyr)

# Source the original Figure 3B script to get all data objects
source(here("Main Figure Scripts", "Helper Scripts", "Figure_3B_4A_Base.R"), local = TRUE)

# After sourcing, we have:
# - boxplot_data: per-tissue data with columns SAMPID, Factor1, GDF15, DTHPLCE, Tissue, rho, p_value, p_adjust, asterisks
# - results_df: per-tissue summary with Tissue, rho, p_adjust, asterisks

# Set up output directories
journal_out_dir <- here("Results", "Figures", "Figure_3B_Summary")
if (!dir.exists(journal_out_dir)) dir.create(journal_out_dir, recursive = TRUE)

individual_tissues_dir <- file.path(journal_out_dir, "Individual_Tissues")
if (!dir.exists(individual_tissues_dir)) dir.create(individual_tissues_dir, recursive = TRUE)

# Define color mapping for DTHPLCE
color_mapping <- c(
  "Hospital inpatient" = "orange",
  "Emergency room" = "maroon",
  "Other" = "gray"
)

# Ensure DTHPLCE is a factor with correct levels
boxplot_data$DTHPLCE <- factor(boxplot_data$DTHPLCE,
                                levels = c("Hospital inpatient", "Emergency room", "Other"),
                                exclude = NULL)

# Replace NA values in DTHPLCE with "Other"
boxplot_data <- boxplot_data %>%
  mutate(DTHPLCE = replace_na(DTHPLCE, "Other"))

# Convert Factor1 and GDF15 to numeric
boxplot_data <- boxplot_data %>%
  mutate(
    Factor1 = as.numeric(as.character(Factor1)),
    GDF15 = as.numeric(as.character(GDF15))
  )

# Remove any rows with NA in Factor1 or GDF15
boxplot_data_clean <- boxplot_data %>%
  filter(!is.na(Factor1), !is.na(GDF15))

################################################################################
# Create multi-panel faceted grid of all tissues
################################################################################

# Get unique tissues and ensure they're ordered
tissues <- unique(boxplot_data_clean$Tissue)
n_tissues <- length(tissues)

cat("Creating multi-panel figure for", n_tissues, "tissues\n")

# Create the faceted plot
faceted_plot <- ggplot(boxplot_data_clean, aes(x = Factor1, y = GDF15, color = DTHPLCE)) +
  geom_point(alpha = 0.3, size = 1.5) +  # Individual sample points
  geom_smooth(method = "lm", se = TRUE, color = "red", fill = "red", alpha = 0.1) +
  facet_wrap(~ Tissue, scales = "free", ncol = 6) +
  scale_color_manual(
    name = "Place of Death",
    values = color_mapping,
    na.value = "gray"
  ) +
  labs(
    title = "Figure 3B: Factor1 vs GDF15 — Individual Samples per Tissue",
    x = "Index Score (Factor1)",
    y = "GDF15 (scaled)",
    caption = "Individual sample points colored by place of death. Red line = linear regression fit with 95% CI."
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 10)),
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
  filename = "Figure_3B_All_Tissues_Faceted.png",
  plot = faceted_plot,
  path = journal_out_dir,
  width = 24,
  height = 30,
  dpi = 300,
  units = "in"
)

cat("Saved faceted plot: Figure_3B_All_Tissues_Faceted.png\n")

################################################################################
# PART 2: Save individual tissue scatter plots to subfolder
################################################################################

cat("Saving individual tissue plots to subfolder...\n")

for (tissue in tissues) {
  # Filter data for this tissue
  tissue_data <- boxplot_data_clean %>%
    filter(Tissue == tissue)

  # Get correlation stats for this tissue
  tissue_stats <- boxplot_data_clean %>%
    filter(Tissue == tissue) %>%
    distinct(Tissue, rho, p_adjust, asterisks) %>%
    slice(1)

  rho_val <- tissue_stats$rho
  asterisks <- tissue_stats$asterisks

  # Create individual tissue plot
  individual_plot <- ggplot(tissue_data, aes(x = Factor1, y = GDF15, color = DTHPLCE)) +
    geom_point(alpha = 0.3, size = 1.5) +
    geom_smooth(method = "lm", se = TRUE, color = "red", fill = "red", alpha = 0.1) +
    scale_color_manual(
      name = "Place of Death",
      values = color_mapping,
      na.value = "gray"
    ) +
    labs(
      title = paste0(tissue, " (ρ = ", round(rho_val, 3), " ", asterisks, ")"),
      x = "Index Score (Factor1)",
      y = "GDF15 (scaled)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 11),
      axis.text = element_text(size = 10),
      legend.position = "bottom",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9)
    )

  # Save individual tissue plot
  tissue_filename <- gsub(" ", "_", tissue)
  ggsave(
    filename = paste0("Figure_3B_", tissue_filename, ".png"),
    plot = individual_plot,
    path = individual_tissues_dir,
    width = 8,
    height = 6,
    dpi = 300,
    units = "in"
  )
}

cat("Saved", length(tissues), "individual tissue plots to:", individual_tissues_dir, "\n")

################################################################################
# PART 3: Summary information
################################################################################

cat("\n========== Figure 3B Summary ==========\n")
cat("Total tissues analyzed:", n_tissues, "\n")
cat("Total individual samples:", nrow(boxplot_data_clean), "\n")
cat("\nPer-tissue summary statistics:\n")
print(results_df)

cat("\n========== Output Files ==========\n")
cat("Multi-panel faceted figure: ", file.path(journal_out_dir, "Figure_3B_All_Tissues_Faceted.png"), "\n")
cat("Individual tissue plots folder: ", individual_tissues_dir, "\n")

################################################################################
