# script to get the average of the spearman rhos of each Factor# vs AGE. Then to compare against GDF15 vs Age, CHOP, ATF4, and others.

#This code currently gets the log2 TMM expression of GDF15, then centers and scales it, and then compared the pc2 score (which also gets centerd and scaled) for each tissue

# Uncommented out many of the "print" plots

rm(list = ls())

# BiocManager::install("edgeR")

library(tidyverse)
library(edgeR)
library(corrr)
library(dplyr)
library(tibble)
library(ggpubr)
library(readr)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(beepr)
library(ggplot2)
library(ggpubr)
library(FSA)      # For Dunn's test
library(dplyr)
library(here)

# Output directories
local_out_dir <- here("Results", "Figures", "Figure_3D")
shared_out_dir <- here("Results", "Figures")
if (!dir.exists(local_out_dir)) dir.create(local_out_dir, recursive = TRUE)
if (!dir.exists(shared_out_dir)) dir.create(shared_out_dir, recursive = TRUE)
folder_path <- local_out_dir  # Keep folder_path for compatibility with rest of script

############################################################################################

script_path <- here("Main Figure Scripts", "Helper Scripts", "gtex", "Attibutes_Phenos_Merged_plus_COD.R")
source(script_path, local = TRUE)


# Set the path to the main folder
CSV_folder_path <- here("Data", "gtex", "Age_comparisons", "Comparing_AGE_Other_Factors_BH")

# List all subdirectories (12 folders)
subfolders <- list.dirs(CSV_folder_path, full.names = TRUE, recursive = FALSE)

# Initialize an empty list to store the dataframes
data_list <- list()

# Loop through each subfolder to read the specific CSV file
for (subfolder in subfolders) {
  # Find the CSV file in the current subfolder
  csv_file <- list.files(subfolder, pattern = "\\.csv$", full.names = TRUE)

  # Check if there's exactly one CSV file, then read it
  if (length(csv_file) == 1) {
    # Extract the file name without the extension to use as the dataframe name
    file_name <- tools::file_path_sans_ext(basename(csv_file))

    # Read the CSV file and store it in the list with the file name as key
    data_list[[file_name]] <- read.csv(csv_file)
  } else {
    message(paste("Skipping folder", subfolder, "due to multiple or no CSV files."))
  }
}

# Initialize an empty dataframe for the result, starting with the "Tissue" column
All_Spearmans <- data.frame(Tissue = unique(data_list[[1]]$Tissue))

# Loop through each dataframe in data_list
for (name in names(data_list)) {
  # Get the current dataframe
  df <- data_list[[name]]

  # Find the column with "Spearman_Rho_chosen" in its name
  p_value_col <- grep("Spearman_Rho_chosen", colnames(df), value = TRUE)

  # Check if we found the exact column and "Tissue" column exists
  if (length(p_value_col) == 1 && "Tissue" %in% colnames(df)) {
    # Select the relevant columns (Tissue and Spearman_Rho_chosen column)
    df_selected <- df[, c("Tissue", p_value_col)]

    # Rename the Spearman_Rho_chosen column to the CSV file name
    colnames(df_selected)[2] <- name

    # Merge this dataframe with All_Spearmans by "Tissue"
    All_Spearmans <- merge(All_Spearmans, df_selected, by = "Tissue", all = TRUE)
  } else {
    message(paste("Skipping file", name, "due to missing or multiple 'Spearman_Rho_chosen' columns."))
  }
}


# Set the path to the main folder
ATF4_csv <- here("Data", "gtex", "Age_comparisons", "Comparing_AGE_Other_Factors_BH_ATF4","2024_10_30_All_Factors_AGE_ATF4", "Factor1_vs_AGE_SpearmanRhos.csv")
ATF4_df <- read.csv(ATF4_csv)
# Merge the Spearman_Rho_ATF4 column from ATF4_df into All_Spearmans by Tissue
All_Spearmans <- merge(All_Spearmans, ATF4_df[, c("Tissue", "Spearman_Rho_ATF4")], by = "Tissue", all.x = TRUE)



# Set the path to the main folder
DDIT3_csv <- here("Data", "gtex", "Age_comparisons", "Comparing_AGE_Other_Factors_BH_DDIT3", "2024_10_30_All_Factors_AGE_DDIT3", "Factor1_vs_AGE_SpearmanRhos.csv")
DDIT3_df <- read.csv(DDIT3_csv)
# Merge the Spearman_Rho_DDIT3 column from DDIT3_df into All_Spearmans by Tissue
All_Spearmans <- merge(All_Spearmans, DDIT3_df[, c("Tissue", "Spearman_Rho_DDIT3")], by = "Tissue", all.x = TRUE)


# Set the path to the main folder
ATF5_csv <-  here("Data", "gtex", "Age_comparisons", "Comparing_AGE_Other_Factors_BH_ATF5", "2024_10_30_All_Factors_AGE_ATF5", "Factor1_vs_AGE_SpearmanRhos.csv")
ATF5_df <- read.csv(ATF5_csv)
# Merge the Spearman_Rho_ATF5 column from ATF5_df into All_Spearmans by Tissue
All_Spearmans <- merge(All_Spearmans, ATF5_df[, c("Tissue", "Spearman_Rho_ATF5")], by = "Tissue", all.x = TRUE)

GDF15_csv <-  here("Data", "gtex", "Age_comparisons", "Comparing_AGE_Other_Factors_BH", "2024_10_28_All_Factors_AGE_Factor1", "Factor1_vs_AGE_SpearmanRhos.csv")
GDF15_df <- read.csv(GDF15_csv)
All_Spearmans <- merge(All_Spearmans, GDF15_df[, c("Tissue", "Spearman_Rho_GDF15")], by = "Tissue", all.x = TRUE)

Proliferation_csv <-  here("Data", "gtex", "Age_comparisons", "2024_12_02_All_Factors_AGE_AvgOf3ProliferationGenes", "Factor1_vs_AGE_SpearmanRhos.csv")
Proliferation_df <- read.csv(Proliferation_csv)
Proliferation_df <- Proliferation_df %>%
  rename(Spearman_Rho_Proliferation = Spearman_Rho)
All_Spearmans <- merge(All_Spearmans, Proliferation_df[, c("Tissue", "Spearman_Rho_Proliferation")], by = "Tissue", all.x = TRUE)

Senescence_csv <-  here("Data", "gtex", "Age_comparisons", "2024_12_02_All_Factors_AGE_AvgOf3SenescentGenes", "Factor1_vs_AGE_SpearmanRhos.csv")
Senescence_df <- read.csv(Senescence_csv)
Senescence_df <- Senescence_df %>%
  rename(Spearman_Rho_Senescence = Spearman_Rho)
All_Spearmans <- merge(All_Spearmans, Senescence_df[, c("Tissue", "Spearman_Rho_Senescence")], by = "Tissue", all.x = TRUE)


# # Calculate mean and SEM for each specified column, assuming All_Spearmans is your dataframe
# plot_columns <- c('Spearman_Rho_GDF15', 'Factor1_vs_AGE_SpearmanRhos', 'Factor2_vs_AGE_SpearmanRhos',
#                   'Factor3_vs_AGE_SpearmanRhos', 'Factor4_vs_AGE_SpearmanRhos', 'Factor5_vs_AGE_SpearmanRhos',
#                   'Factor6_vs_AGE_SpearmanRhos', 'Factor7_vs_AGE_SpearmanRhos', 'Factor8_vs_AGE_SpearmanRhos',
#                   'Factor9_vs_AGE_SpearmanRhos', 'Factor10_vs_AGE_SpearmanRhos', 'Factor11_vs_AGE_SpearmanRhos',
#                   'Factor12_vs_AGE_SpearmanRhos', 'Spearman_Rho_ATF4', 'Spearman_Rho_DDIT3', 'Spearman_Rho_ATF5')

# Calculate mean and SEM for each specified column, assuming All_Spearmans is your dataframe
plot_columns <- c('Spearman_Rho_GDF15', 'Spearman_Rho_Senescence', 'Spearman_Rho_Proliferation',
                  'Spearman_Rho_ATF4', 'Spearman_Rho_ATF5', 'Spearman_Rho_DDIT3', 'Factor1_vs_AGE_SpearmanRhos',
                  'Factor12_vs_AGE_SpearmanRhos', 'Factor6_vs_AGE_SpearmanRhos', 'Factor2_vs_AGE_SpearmanRhos',
                  'Factor11_vs_AGE_SpearmanRhos', 'Factor4_vs_AGE_SpearmanRhos', 'Factor9_vs_AGE_SpearmanRhos',
                  'Factor8_vs_AGE_SpearmanRhos', 'Factor3_vs_AGE_SpearmanRhos', 'Factor10_vs_AGE_SpearmanRhos',
                  'Factor7_vs_AGE_SpearmanRhos', 'Factor5_vs_AGE_SpearmanRhos'
                  )

# Compute means and SEM for each column
means <- sapply(All_Spearmans[plot_columns], mean, na.rm = TRUE)
sems <- sapply(All_Spearmans[plot_columns], function(x) sd(x, na.rm = TRUE) / sqrt(length(na.omit(x))))

# Convert to a dataframe for plotting
plot_data <- data.frame(
  Factor = plot_columns,
  Mean = means,
  SEM = sems
)




# Define custom colors for each specific column
color_mapping <- c(
  "Spearman_Rho_GDF15" = "pink",
  "Factor1_vs_AGE_SpearmanRhos" = "red",
  "Factor2_vs_AGE_SpearmanRhos" = "#808080",
  "Factor3_vs_AGE_SpearmanRhos" = "#808080",
  "Factor4_vs_AGE_SpearmanRhos" = "#CCCCCC",
  "Factor5_vs_AGE_SpearmanRhos" = "#CCCCCC",
  "Factor6_vs_AGE_SpearmanRhos" = "#CCCCCC",
  "Factor7_vs_AGE_SpearmanRhos" = "#808080",
  "Factor8_vs_AGE_SpearmanRhos" = "#CCCCCC",
  "Factor9_vs_AGE_SpearmanRhos" = "#808080",
  "Factor10_vs_AGE_SpearmanRhos" = "#CCCCCC",
  "Factor11_vs_AGE_SpearmanRhos" = "#808080",
  "Factor12_vs_AGE_SpearmanRhos" = "#808080",
  "Spearman_Rho_ATF4" = "orange",  # tomato red
  "Spearman_Rho_DDIT3" = "lightblue",  # firebrick red
  "Spearman_Rho_ATF5" = "brown",    # dark red
  "Spearman_Rho_Senescence" = "blue",
  "Spearman_Rho_Proliferation" = "purple"
)

# plot_data$Factor <- factor(plot_data$Factor, levels = plot_columns)

# Ensure the levels of Factor match the keys in color_mapping
plot_data$Factor <- factor(plot_data$Factor, levels = names(color_mapping))


plot_columns_csv <- c('Tissue', 'Spearman_Rho_GDF15', 'Spearman_Rho_Senescence', 'Spearman_Rho_Proliferation',
                  'Spearman_Rho_ATF4', 'Spearman_Rho_ATF5', 'Spearman_Rho_DDIT3', 'Factor1_vs_AGE_SpearmanRhos',
                  'Factor12_vs_AGE_SpearmanRhos', 'Factor6_vs_AGE_SpearmanRhos', 'Factor2_vs_AGE_SpearmanRhos',
                 'Factor4_vs_AGE_SpearmanRhos', 'Factor9_vs_AGE_SpearmanRhos',
                  'Factor8_vs_AGE_SpearmanRhos', 'Factor3_vs_AGE_SpearmanRhos', 'Factor10_vs_AGE_SpearmanRhos',
                  'Factor7_vs_AGE_SpearmanRhos', 'Factor5_vs_AGE_SpearmanRhos',  'Factor11_vs_AGE_SpearmanRhos'
)



# Ensure Factor column has the desired order
# Reverse the order
plot_columns_csv_reversed <- rev(plot_columns_csv)

# Ensure Factor column has the reversed order
plot_data$Factor <- factor(plot_data$Factor, levels = plot_columns_csv_reversed)


# NOTE: All GTEx plotting has been moved to Figure_3D.R.
# This script is now strictly a data-prep step that builds All_Spearmans /
# plot_data and writes the source-data CSVs that Figure_3D.R consumes.

# ===================== Save Figure 3D GTEx source data =====================
write.csv(
  plot_data,
  file = file.path(folder_path, "Figure_3D_GTEx_source_data.csv"),
  row.names = FALSE
)
message("Saved Figure 3D GTEx source data to: ", file.path(folder_path, "Figure_3D_GTEx_source_data.csv"))

# ===================== Save per-tissue Spearman rho CSVs for each metric =====================
tissue_csv_dir <- file.path(folder_path, "Per_Tissue_Spearman_Rhos")
if (!dir.exists(tissue_csv_dir)) dir.create(tissue_csv_dir, recursive = TRUE)

for (col_name in plot_columns) {
  if (col_name %in% colnames(All_Spearmans)) {
    tissue_data <- All_Spearmans[, c("Tissue", col_name), drop = FALSE]
    tissue_data <- tissue_data[order(tissue_data[[col_name]], decreasing = TRUE), ]
    csv_name <- paste0(col_name, "_per_tissue.csv")
    write.csv(tissue_data, file = file.path(tissue_csv_dir, csv_name), row.names = FALSE)
  }
}
message("Saved per-tissue Spearman rho CSVs to: ", tissue_csv_dir)

# ============================================================================
# FIGURE 3D - FIBROBLAST VERSION
# ============================================================================
message("\n========== Generating Figure 3D Fibroblast version ==========")

# Read fibroblast Spearman Rho summary with confidence intervals
fibroblast_csv <- here("Results", "Fibroblast_lifespan", "Factor_Analysis",
                       "Generate_All_Figures", "Spearman_Rhos_Age_vs_Genes",
                       "Spearman_Rho_Summary_With_CI.csv")

if (!file.exists(fibroblast_csv)) {
  message("WARNING: Fibroblast Spearman Rho data not found at: ", fibroblast_csv)
  message("Skipping Fibroblast version of Figure 3D.")
} else {
  fb_data <- read.csv(fibroblast_csv, stringsAsFactors = FALSE)

  # Use original Column names as labels (GDF15, sen_mean, Factor1, etc.)
  fb_plot_data <- fb_data %>%
    select(Column, Spearman_Rho, CI_Lower, CI_Upper, P_Value, Adjusted_P_Value, Significant)

  # Order to match GTEx plot (top to bottom: GDF15, sen_mean, prolif_mean, ATF4, ATF5, DDIT3, Factor1, Factor12, ...)
  fb_order <- c('GDF15', 'sen_mean', 'prolif_mean',
                'ATF4', 'ATF5', 'DDIT3',
                'Factor1', 'Factor12', 'Factor6', 'Factor2',
                'Factor11', 'Factor4', 'Factor9', 'Factor8',
                'Factor3', 'Factor10', 'Factor7', 'Factor5')

  fb_plot_data <- fb_plot_data %>%
    filter(Column %in% fb_order)

  fb_plot_data$Column <- factor(fb_plot_data$Column, levels = rev(fb_order))

  # Color mapping using original names
  fb_color_mapping <- c(
    "GDF15" = "pink",
    "sen_mean" = "blue",
    "prolif_mean" = "purple",
    "ATF4" = "orange",
    "ATF5" = "brown",
    "DDIT3" = "lightblue",
    "Factor1" = "red",
    "Factor2" = "#808080",
    "Factor3" = "#808080",
    "Factor4" = "#CCCCCC",
    "Factor5" = "#CCCCCC",
    "Factor6" = "#CCCCCC",
    "Factor7" = "#808080",
    "Factor8" = "#CCCCCC",
    "Factor9" = "#808080",
    "Factor10" = "#CCCCCC",
    "Factor11" = "#808080",
    "Factor12" = "#808080"
  )

  # NOTE: Fibroblast plotting has been moved to Figure_3D.R.
  # This block now only writes the source-data CSV that Figure_3D.R consumes.

  # Save Fibroblast source data CSV (rename Column back to Factor for consistency)
  fb_export <- fb_plot_data %>% rename(Factor = Column)
  write.csv(
    fb_export,
    file = file.path(folder_path, "Figure_3D_Fibroblast_source_data.csv"),
    row.names = FALSE
  )

  message("Saved Figure 3D Fibroblast source data to: ", file.path(folder_path, "Figure_3D_Fibroblast_source_data.csv"))
}
