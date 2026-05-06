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
local_out_dir <- here("Results", "Figures", "Figure_3C")
shared_out_dir <- here("Results", "Figures")
if (!dir.exists(local_out_dir)) dir.create(local_out_dir, recursive = TRUE)
if (!dir.exists(shared_out_dir)) dir.create(shared_out_dir, recursive = TRUE)
folder_path <- local_out_dir  # Keep folder_path for compatibility

Plot_Save = "ON"
chosen_Factor <- "Factor1"
#------------------------------------------------------
Folder_name <- paste0("2024_10_28_All_Factors_AGE_", chosen_Factor)
#------------------------------------------------------


folder_path_1stfolder <- folder_path

# Check if the folder exists
if (!dir.exists(folder_path_1stfolder)) {
  # If it doesn't exist, create the folder
  dir.create(folder_path_1stfolder, recursive = TRUE)
  cat("Folder created at:", folder_path_1stfolder, "\n")
} else {
  cat("Folder already exists at:", folder_path_1stfolder, "\n")
}

folder_path <- folder_path_1stfolder


# Read necessary input files
Atts_COD_readin <- here("Main Figure Scripts", "Helper Scripts", "gtex", "Attibutes_Phenos_Merged_plus_COD.R")
source(Atts_COD_readin, local = TRUE)

getwd()

# correlations_readin <- file.path(paste0(gtex_input, "/2024_09_23_correlations_for_scaling.csv"))
# correlations <- read.csv(correlations_readin)


# Path to the folder containing individual tissue folders
#########################################################################################################
individual_tissues_folder <- here("Data", "gtex", "All_Tissues_Indiv_Folders")
# individual_tissues_folder <- here("Data", "gtex", "Single_Tissue", "brain") #using brain only for testing
#########################################################################################################

# Debug: print the path to check if it is correct
print(individual_tissues_folder)

# ============================================================================
# CACHING: Check if pre-computed Factor1 scores exist
# ============================================================================
cached_scores_file <- file.path(folder_path, "Tissue_Factor1scores_per_SAMPID.csv")
USE_CACHED_DATA <- file.exists(cached_scores_file)

if (USE_CACHED_DATA) {
  message("\n========== USING CACHED DATA ==========")
  message("Loading pre-computed Factor1 scores from: ", cached_scores_file)

  # Load cached data
  all_tissue_scores <- read.csv(cached_scores_file)
  message("Loaded ", nrow(all_tissue_scores), " samples from ",
          length(unique(all_tissue_scores$Tissue)), " tissues")

  # Create results_df with both raw and adjusted p-values
  results_df <- data.frame(Tissue = character(),
                           Comparison = character(),
                           Spearman_Rho = numeric(),
                           P_Value_Raw = numeric(),
                           P_Value_Adjusted = numeric(),
                           Asterisks = character(),
                           N_obs = integer(),
                           CI_lower = numeric(),
                           CI_upper = numeric(),
                           stringsAsFactors = FALSE)

  # Compute correlations from cached data
  tissues <- unique(all_tissue_scores$Tissue)

  for (tissue in tissues) {
    tissue_data <- all_tissue_scores %>% filter(Tissue == tissue)

    # Spearman correlation for Factor1 vs AGE
    rho_chosen_Factor <- cor.test(tissue_data$AGE, tissue_data$Factor1, method = "spearman")

    # Spearman correlation for GDF15 vs AGE
    rho_gdf15 <- cor.test(tissue_data$AGE, tissue_data$GDF15, method = "spearman")

    # Compute n_obs and 95% CIs using Fisher z-transformation
    n_obs_tissue <- nrow(tissue_data)

    # For chosen_Factor
    z_cf <- atanh(rho_chosen_Factor$estimate)
    se_z_cf <- 1 / sqrt(max(n_obs_tissue - 3, 1))
    ci_lower_cf <- tanh(z_cf - 1.96 * se_z_cf)
    ci_upper_cf <- tanh(z_cf + 1.96 * se_z_cf)

    # For GDF15
    z_gdf15 <- atanh(rho_gdf15$estimate)
    se_z_gdf15 <- 1 / sqrt(max(n_obs_tissue - 3, 1))
    ci_lower_gdf15 <- tanh(z_gdf15 - 1.96 * se_z_gdf15)
    ci_upper_gdf15 <- tanh(z_gdf15 + 1.96 * se_z_gdf15)

    # Add to results (p-value adjustment will be done after all tissues)
    results_df <- rbind(results_df, data.frame(
      Tissue = tissue,
      Comparison = "chosen_Factor",
      Spearman_Rho = rho_chosen_Factor$estimate,
      P_Value_Raw = rho_chosen_Factor$p.value,
      P_Value_Adjusted = NA,  # Will be filled in after loop
      Asterisks = "",
      N_obs = n_obs_tissue,
      CI_lower = ci_lower_cf,
      CI_upper = ci_upper_cf
    ))

    results_df <- rbind(results_df, data.frame(
      Tissue = tissue,
      Comparison = "GDF15",
      Spearman_Rho = rho_gdf15$estimate,
      P_Value_Raw = rho_gdf15$p.value,
      P_Value_Adjusted = NA,
      Asterisks = "",
      N_obs = n_obs_tissue,
      CI_lower = ci_lower_gdf15,
      CI_upper = ci_upper_gdf15
    ))
  }

  # Apply Bonferroni correction across all tests
  results_df$P_Value_Adjusted <- p.adjust(results_df$P_Value_Raw, method = "bonferroni")

  # Add asterisks based on adjusted p-values
  results_df$Asterisks <- ifelse(results_df$P_Value_Adjusted < 0.0001, "****",
                                  ifelse(results_df$P_Value_Adjusted < 0.001, "***",
                                         ifelse(results_df$P_Value_Adjusted < 0.01, "**",
                                                ifelse(results_df$P_Value_Adjusted < 0.05, "*", "ns"))))

  message("Computed correlations for ", length(tissues), " tissues")
  message("========================================\n")

} else {
  message("\n========== PROCESSING ALL TISSUES ==========")
  message("No cached data found. Processing all tissue folders...")
  message("This may take a while. Results will be cached for future runs.\n")

# Create empty dataframe to store results (with both raw and adjusted p-values)
results_df <- data.frame(Tissue = character(),
                         Comparison = character(),
                         Spearman_Rho = numeric(),
                         P_Value_Raw = numeric(),
                         P_Value_Adjusted = numeric(),
                         Asterisks = character(),
                         N_obs = integer(),
                         CI_lower = numeric(),
                         CI_upper = numeric(),
                         stringsAsFactors = FALSE)


# Initialize an empty list to store results from each tissue
all_results <- list()

all_dunn_res <- data.frame()
# Initialize lists to collect data needed for plotting
boxplot_data_list <- list()
cohens_d_results_list <- list()

# Initialize list to collect Factor1 scores for all samples
all_tissue_scores_list <- list()


# Get a list of all subdirectories within the main folder
subfolders <- list.dirs(individual_tissues_folder, recursive = FALSE, full.names = TRUE)

# Loop through each subfolder
for (folder in subfolders) {
  # Print the folder name
  print(paste("Processing folder:", basename(folder)))

  last_component <- basename(folder)

  # Set the path to the current folder
  current_folder <- folder

  correlations_readin <- here("Data", "Fibroblast_lifespan", "Processed", "Fibroblast_FactorAnalysis_correlations_for_scaling.csv")
  correlations <- read.csv(correlations_readin)



  # Print the last component
  print(last_component)


  # new_folder_path <- paste0(folder_path, "/3_Outputs/gtex/", Folder_name, "/", last_component)
  #
  # # Check if the folder exists
  # if (!dir.exists(new_folder_path)) {
  #   # If it doesn't exist, create the folder
  #   dir.create(new_folder_path, recursive = TRUE)
  #   cat("Folder created at:", new_folder_path, "\n")
  # } else {
  #   cat("Folder already exists at:", new_folder_path, "\n")
  # }
  #

  # Get a list of all .gz files in the tissue folder
  gz_files <- list.files(folder, pattern = "\\.gz$", full.names = TRUE) # CHANGE BACK TO THIS AFTER DONE WITH BRAIN
  # gz_files <- list.files(individual_tissues_folder, pattern = "\\.gz$", full.names = TRUE)

  # Initialize dataframes for each tissue
  df_tissue_list <- list()
  gene_symbols <- NULL

  # Loop through each .gz file and read data
  for (gz_file in gz_files) {
    file_name <- basename(gz_file)
    tissue <- sub(".*_v8_(.*?)\\.gct\\.gz", "\\1", file_name)

    # Read data from the current .gz file
    raw_gene_counts <- read.delim(gz_file, skip = 2)

    # Extract gene symbols if not already done
    if (is.null(gene_symbols)) {
      gene_symbols <- raw_gene_counts[, 2:3]
    }

    # Extract sample IDs and gene counts
    gene_counts <- raw_gene_counts %>%
      select(-id, -Description) %>%
      as.data.frame()

    # Store the gene counts in the dataframe list with tissue name
    df_tissue_list[[tissue]] <- gene_counts
  }

  # Combine dataframes if there are multiple files
  combined_gene_counts <- gene_symbols
  for (tissue in names(df_tissue_list)) {
    combined_gene_counts <- combined_gene_counts %>%
      left_join(df_tissue_list[[tissue]], by = "Name")
  }

  # Remove the "Description" column if it exists
  if ("Description" %in% colnames(combined_gene_counts)) {
    combined_gene_counts <- combined_gene_counts %>%
      select(-Description)
  }

  combined_gene_counts <- combined_gene_counts %>%
    column_to_rownames("Name")

  annotations_sample_attributes_filtered <- annotations_merged %>%
    filter(SAMPID %in% colnames(combined_gene_counts))

  # TMM normalization
  dge <- DGEList(combined_gene_counts[, annotations_sample_attributes_filtered$SAMPID], group = factor(annotations_sample_attributes_filtered$COD))
  keep <- filterByExpr(dge)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge, method = "TMM")
  exprs <- cpm(dge, log = TRUE) # changed from FALST to TRUE 9/27/2024

  gene_symbols <- gene_symbols %>%
    mutate(Description = ifelse(Description == "KIAA0141", "DELE1", Description),  # changing so that DELE1 is the gene symbol
           Description = ifelse(Description == "WARS", "WARS1", Description),
           Description = ifelse(Description == "NARS", "NARS1", Description))


  # Join with gene symbols
  exprs <- exprs %>%
    as.data.frame() %>%
    rownames_to_column("Name") %>%
    full_join(gene_symbols, by = "Name") %>%
    na.omit() %>%
    select(-Name) %>%
    rename("Gene" = Description)

  # Filter genes of the ISR only into our exprs
  exprs_ISR <- exprs %>%
    filter(Gene %in% Total_ISR_List$Gene) %>%
    column_to_rownames("Gene") %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("SAMPID") %>%
    unique()


  #------------------ finding missing genes ----------------------

  # Step 1: Extract gene names
  exprs_genes <- colnames(exprs_ISR)[-1]
  list_genes <- Total_ISR_List$Gene

  # Step 2: Standardize gene names
  exprs_genes_clean <- trimws(toupper(exprs_genes))
  list_genes_clean <- trimws(toupper(list_genes))

  # Remove duplicates if necessary
  exprs_genes_clean <- unique(exprs_genes_clean)
  list_genes_clean <- unique(list_genes_clean)

  # Step 3: Find missing genes
  genes_in_exprs_not_in_list <- setdiff(exprs_genes_clean, list_genes_clean)
  genes_in_list_not_in_exprs <- setdiff(list_genes_clean, exprs_genes_clean)

  # Step 4: Print results
  if (length(genes_in_exprs_not_in_list) > 0) {
    cat("Genes in exprs_ISR but not in Total_ISR_List:\n")
    print(genes_in_exprs_not_in_list)
  } else {
    cat("All genes in exprs_ISR are present in Total_ISR_List.\n")
  }

  if (length(genes_in_list_not_in_exprs) > 0) {
    cat("Genes in Total_ISR_List but not in exprs_ISR:\n")
    print(genes_in_list_not_in_exprs)
  } else {
    cat("All genes in Total_ISR_List are present in exprs_ISR.\n")
  }

  # Check if there are missing genes
  if (length(genes_in_list_not_in_exprs) > 0 || length(genes_in_exprs_not_in_list) > 0 ) {
    cat(paste("Missing genes detected in folder:", last_component, "- Skipping this folder.\n"))
    # Skip to the next iteration of the loop
    next
  }

  #-------------------------------------------------------------------------------------------------------

  # CHANGING THIS mean centering and scaling needs to be performed using fibroblast data
  ###################################################################################################
  pathtoMEANofX <- here("Data", "Fibroblast_lifespan", "Processed", "Fibroblast_FactorAnalysis_meanx_for_scaling.csv")
  mean_x_fbdata <- read.csv(pathtoMEANofX)

  pathtoSDofX <- here("Data", "Fibroblast_lifespan", "Processed", "Fibroblast_FactorAnalysis_sdx_for_scaling.csv")
  sd_x_fbdata <- read.csv(pathtoSDofX)
  ###################################################################################################

  # Rename columns using dplyr
  sd_x_fbdata <- sd_x_fbdata %>%
    rename(
      Gene = X,
      sd = x
    )

  # Rename columns using dplyr
  mean_x_fbdata <- mean_x_fbdata %>%
    rename(
      Gene = X,
      mean = x
    )

  sd_x_fbdata <- sd_x_fbdata %>%
    left_join(mean_x_fbdata %>% select(Gene, mean), by = "Gene")


  # Step 1: Extract gene names from exprs_ISR (excluding SampID)
  gene_columns <- setdiff(names(exprs_ISR), "SAMPID")

  # Step 2: Ensure that gene names match between exprs_ISR and sd_x_fbdata
  common_genes <- intersect(gene_columns, sd_x_fbdata$Gene)
  # view(common_genes)


  #Finding any missing genes
  # Step 1: Extract gene names
  exprs_genes <- setdiff(names(exprs_ISR), "SAMPID")
  sd_genes <- sd_x_fbdata$Gene

  # Step 2: Standardize gene names
  exprs_genes_clean <- trimws(toupper(exprs_genes))
  sd_genes_clean <- trimws(toupper(sd_genes))

  # Step 3: Find missing genes
  genes_in_exprs_not_in_sd <- setdiff(exprs_genes_clean, sd_genes_clean)
  genes_in_sd_not_in_exprs <- setdiff(sd_genes_clean, exprs_genes_clean)

  # Step 4: Print missing genes
  if (length(genes_in_exprs_not_in_sd) > 0) {
    cat("Genes in exprs_ISR but not in sd_x_fbdata:\n")
    print(genes_in_exprs_not_in_sd)
  } else {
    cat("All genes in exprs_ISR are present in sd_x_fbdata.\n")
  }

  if (length(genes_in_sd_not_in_exprs) > 0) {
    cat("Genes in sd_x_fbdata but not in exprs_ISR:\n")
    print(genes_in_sd_not_in_exprs)
  } else {
    cat("All genes in sd_x_fbdata are present in exprs_ISR.\n")
  }




  #Scaling the data
  # Step 3: Subset exprs_ISR to include only the common genes
  exprs_ISR_subset <- exprs_ISR[, c("SAMPID", common_genes)]

  # Step 4: Create named vectors for mean and sd
  mean_values <- sd_x_fbdata$mean
  names(mean_values) <- sd_x_fbdata$Gene

  sd_values <- sd_x_fbdata$sd
  names(sd_values) <- sd_x_fbdata$Gene

  # Subset mean and sd to include only common genes
  mean_values <- mean_values[common_genes]
  sd_values <- sd_values[common_genes]

  # Check for zero standard deviations
  if (any(sd_values == 0)) {
    stop("Standard deviation is zero for some genes. Cannot perform scaling.")
  }

  # Step 5: Manually scale the data
  exprs_ISR_scaled <- exprs_ISR_subset

  for (gene in common_genes) {
    exprs_ISR_scaled[[gene]] <- (exprs_ISR_subset[[gene]] - mean_values[gene]) / sd_values[gene]
  }


  exprs_gdf15 <- exprs_ISR_scaled %>%
    select(SAMPID, GDF15)

  # # Join with subject phenotypes
  exprs_gdf15 <- exprs_gdf15 %>%
    full_join(annotations_sample_attributes_filtered, by = "SAMPID") %>%
    select(GDF15, SAMPID, AGE, SEX) %>%
    unique() %>%
    na.omit()
  #
  # # Join with subject phenotypes
  exprs_gdf15_tissue <- exprs_gdf15 %>%
    full_join(annotations_sample_attributes_filtered, by = "SAMPID") %>%
    select(GDF15, SMTSD, SAMPID) %>%
    unique() %>%
    na.omit()

  ########################################################################################

  scaled_exprs_ISR <- exprs_ISR_scaled %>%
    full_join(annotations_sample_attributes_filtered, by = "SAMPID") %>%
    unique()
  #
  # scaled_exprs_ISR <- scaled_exprs_ISR %>% # not sure how this is different
  #   select(SMTSD, SUBJID, everything()) %>%
  #   unique()

  tissue_results <- data.frame(PC = character(), Spearman_Rho = numeric(), P_Value = numeric(), regression_coefficient = numeric(), reg_p_value = numeric(), N = numeric(), stringsAsFactors = FALSE)


  loadings_readin <- here("Data", "Fibroblast_lifespan", "Processed", "Fibroblast_FactorAnalysis_FacLoads_12_factors.csv") # produced by Figure_1E_1F.R
  loadings <- read.csv(loadings_readin)
  loadings <- loadings %>%
    select(X, chosen_Factor)


  selected_loadings <- loadings %>%
    rename(Gene = X, Loading = chosen_Factor) %>%
    as.data.frame()

  common_genes <- intersect(selected_loadings$Gene, colnames(scaled_exprs_ISR))
  selected_loadings <- selected_loadings %>%
    filter(Gene %in% common_genes)

  # Step 1: Prepare the correlations matrix
  rownames(loadings) <- loadings$X
  loadings$X <- NULL
  loadings_matrix <- as.matrix(loadings)

  scaled_data <- scaled_exprs_ISR %>%
    select(SAMPID, all_of(common_genes))

  rownames(scaled_data) <- scaled_data$SAMPID
  scaled_data$SAMPID <- NULL
  scaled_data_matrix <- as.matrix(scaled_data)


  # t(loadings) %*% solve(correlations)
  # Step 1: Prepare the correlations matrix
  rownames(correlations) <- correlations$X
  correlations$X <- NULL
  correlations_matrix <- as.matrix(correlations)

  W <- t(loadings_matrix) %*% solve(correlations_matrix)
  Wtransposed <- t(W)
  final <- scaled_data_matrix %*% Wtransposed

  final_df <- as.data.frame(final)

  final_df <- final_df %>%
    tibble::rownames_to_column(var = "SAMPID")

  final_data <- final_df %>%
    full_join(annotations_sample_attributes_filtered, by = "SAMPID") %>%
    unique()

  ########################################################################################
  ########################################################################################


  boxplot_data <- final_data %>%
    select(SAMPID, chosen_Factor, AGE, SEX, DTHPLCE)

  colnames(boxplot_data)[2] <- "chosen_Factor"

  exprs_GDF15 <- exprs_ISR_scaled %>%
    select(SAMPID, GDF15)

  boxplot_data <- boxplot_data %>%
    left_join(exprs_GDF15, by = "SAMPID")

  # boxplot_data <- boxplot_data %>%
  #   filter(DTHPLCE %in% c("Hospital inpatient", "Emergency room"))


  # SCATTERPLOTS

  # Scatterplot for AGE vs chosen_Factor
  rho_chosen_Factor <- cor.test(boxplot_data$AGE, boxplot_data$chosen_Factor, method = "spearman")
  adjusted_p_chosen_Factor <- p.adjust(rho_chosen_Factor$p.value, method = "BH")

  asterisks_chosen_Factor <- ifelse(adjusted_p_chosen_Factor < 0.0001, "****",
                                    ifelse(adjusted_p_chosen_Factor < 0.001, "***",
                                           ifelse(adjusted_p_chosen_Factor < 0.01, "**",
                                                  ifelse(adjusted_p_chosen_Factor < 0.05, "*", "ns"))))

  # Scatterplot for AGE vs GDF15
  rho_gdf15 <- cor.test(boxplot_data$AGE, boxplot_data$GDF15, method = "spearman")
  adjusted_p_gdf15 <- p.adjust(rho_gdf15$p.value, method = "BH")

  asterisks_gdf15 <- ifelse(adjusted_p_gdf15 < 0.0001, "****",
                            ifelse(adjusted_p_gdf15 < 0.001, "***",
                                   ifelse(adjusted_p_gdf15 < 0.01, "**",
                                          ifelse(adjusted_p_gdf15 < 0.05, "*", "ns"))))

  # Compute n_obs and 95% CIs using Fisher z-transformation
  n_obs_tissue <- nrow(boxplot_data)

  # For chosen_Factor
  z_cf <- atanh(rho_chosen_Factor$estimate)
  se_z_cf <- 1 / sqrt(max(n_obs_tissue - 3, 1))
  ci_lower_cf <- tanh(z_cf - 1.96 * se_z_cf)
  ci_upper_cf <- tanh(z_cf + 1.96 * se_z_cf)

  # For GDF15
  z_gdf15 <- atanh(rho_gdf15$estimate)
  se_z_gdf15 <- 1 / sqrt(max(n_obs_tissue - 3, 1))
  ci_lower_gdf15 <- tanh(z_gdf15 - 1.96 * se_z_gdf15)
  ci_upper_gdf15 <- tanh(z_gdf15 + 1.96 * se_z_gdf15)

  # Add results to results dataframe (storing RAW p-values - adjustment done after loop)
  results_df <- rbind(results_df, data.frame(Tissue = tissue,
                                             Comparison = "chosen_Factor",
                                             Spearman_Rho = rho_chosen_Factor$estimate,
                                             P_Value_Raw = rho_chosen_Factor$p.value,
                                             P_Value_Adjusted = NA,
                                             Asterisks = asterisks_chosen_Factor,
                                             N_obs = n_obs_tissue,
                                             CI_lower = ci_lower_cf,
                                             CI_upper = ci_upper_cf))

  results_df <- rbind(results_df, data.frame(Tissue = tissue,
                                             Comparison = "GDF15",
                                             Spearman_Rho = rho_gdf15$estimate,
                                             P_Value_Raw = rho_gdf15$p.value,
                                             P_Value_Adjusted = NA,
                                             Asterisks = asterisks_gdf15,
                                             N_obs = n_obs_tissue,
                                             CI_lower = ci_lower_gdf15,
                                             CI_upper = ci_upper_gdf15))

  # Collect Factor1 scores for caching
  tissue_scores <- boxplot_data %>%
    select(SAMPID, chosen_Factor, AGE, SEX, GDF15) %>%
    mutate(Tissue = tissue) %>%
    rename(Factor1 = chosen_Factor)

  all_tissue_scores_list[[tissue]] <- tissue_scores

  # Generate scatterplots
  p1 <- ggplot(boxplot_data, aes(x = AGE, y = chosen_Factor)) +
    geom_point(color = "gray") +
    geom_smooth(method = "lm", color = "red", se = TRUE) +
    annotate("text", x = Inf, y = Inf, label = paste0("rho: ", round(rho_chosen_Factor$estimate, 3),
                                                      "\nP: ", asterisks_chosen_Factor),
             hjust = 1.1, vjust = 2, size = 5) +
    ggtitle(paste("AGE vs ", chosen_Factor, "in", tissue)) +
    theme_minimal()

  p2 <- ggplot(boxplot_data, aes(x = AGE, y = GDF15)) +
    geom_point(color = "violet") +
    geom_smooth(method = "lm", color = "red", se = TRUE) +
    annotate("text", x = Inf, y = Inf, label = paste0("rho: ", round(rho_gdf15$estimate, 3),
                                                      "\nP: ", asterisks_gdf15),
             hjust = 1.1, vjust = 2, size = 5) +
    ggtitle(paste("AGE vs GDF15 in", tissue)) +
    theme_minimal()

  print(p1)
  print(p2)


  ggsave(filename = paste0(folder_path, "/", last_component, "_AGEvs", chosen_Factor, "_Scatterplot.png"), plot = p1, width = 10, height = 8, dpi = 300)
  #ggsave(filename = paste0(folder_path, "/", last_component, "_AGEvsGDF15_Scatterplot.png"), plot = p2, width = 10, height = 8, dpi = 300)



  folder_path_colorcodeDTHPLCE <- paste0(folder_path, "/DTHPLCE_colorcode")


  # Check if the folder exists
  if (!dir.exists(folder_path_colorcodeDTHPLCE)) {
    # If it doesn't exist, create the folder
    dir.create(folder_path_colorcodeDTHPLCE, recursive = TRUE)
    cat("Folder created at:", folder_path_colorcodeDTHPLCE, "\n")
  } else {
    cat("Folder already exists at:", folder_path_colorcodeDTHPLCE, "\n")
  }




  # Scatterplot with color coding, Spearman rho, and p-value asterisks
  plot_with_spearman_and_dthplce <- function(data, x_var, y_var, rho_test, title) {
    ggplot(data, aes_string(x = x_var, y = y_var)) +
      geom_point(aes(color = DTHPLCE), alpha = 0.6, size = 3) +  # Adjust transparency with alpha
      scale_color_manual(values = c("Hospital inpatient" = "orange",
                                    "Emergency room" = "maroon",
                                    "Other" = "gray")) +
      geom_smooth(method = "lm", color = "red", se = TRUE) +  # Add red linear regression line with SEM
      annotate("text", x = Inf, y = Inf, label = paste0("rho: ", round(rho_test$estimate, 3),
                                                        "\nP: ", rho_test$asterisks),
               hjust = 1.1, vjust = 2, size = 5) +  # Spearman rho and p-value asterisks
      ggtitle(title) +
      theme_minimal()
  }

  # Recode DTHPLCE column to categorize other values as "Other"
  boxplot_data <- boxplot_data %>%
    mutate(DTHPLCE = ifelse(DTHPLCE %in% c("Hospital inpatient", "Emergency room"), DTHPLCE, "Other"))

  # Spearman tests for chosen_Factor and GDF15
  rho_chosen_Factor <- cor.test(boxplot_data$AGE, boxplot_data$chosen_Factor, method = "spearman")
  adjusted_p_chosen_Factor <- p.adjust(rho_chosen_Factor$p.value, method = "BH")
  rho_chosen_Factor$asterisks <- ifelse(adjusted_p_chosen_Factor < 0.0001, "****",
                                        ifelse(adjusted_p_chosen_Factor < 0.001, "***",
                                               ifelse(adjusted_p_chosen_Factor < 0.01, "**",
                                                      ifelse(adjusted_p_chosen_Factor < 0.05, "*", "ns"))))

  rho_gdf15 <- cor.test(boxplot_data$AGE, boxplot_data$GDF15, method = "spearman")
  adjusted_p_gdf15 <- p.adjust(rho_gdf15$p.value, method = "BH")
  rho_gdf15$asterisks <- ifelse(adjusted_p_gdf15 < 0.0001, "****",
                                ifelse(adjusted_p_gdf15 < 0.001, "***",
                                       ifelse(adjusted_p_gdf15 < 0.01, "**",
                                              ifelse(adjusted_p_gdf15 < 0.05, "*", "ns"))))

  # Generate scatterplots
  p1_dthplce <- plot_with_spearman_and_dthplce(
    boxplot_data,
    "AGE",
    "chosen_Factor",  # Use the column name as a string
    rho_chosen_Factor,
    paste0(last_component, "AGE vs ", chosen_Factor, " (color by DTHPLCE)")
  )

  p2_dthplce <- plot_with_spearman_and_dthplce(
    boxplot_data,
    "AGE",
    "GDF15",
    rho_gdf15,
    paste(last_component, "AGE vs GDF15 (color by DTHPLCE)")
  )

  # print(p1_dthplce)
  # print(p2_dthplce)



  ggsave(filename = paste0(folder_path_colorcodeDTHPLCE, "/", last_component, "_AGEvs", chosen_Factor, "_Scatterplot_DTHPLCE.png"), plot = p1_dthplce, width = 10, height = 8, dpi = 300)
  #ggsave(filename = paste0(folder_path_colorcodeDTHPLCE, "/", last_component, "_AGEvsGDF15_Scatterplot_DTHPLCE.png"), plot = p2_dthplce, width = 10, height = 8, dpi = 300)


}

# ============================================================================
# SAVE CACHED DATA: Factor1 scores per SAMPID
# ============================================================================
message("\n========== SAVING CACHED DATA ==========")

# Combine all tissue scores into one dataframe
all_tissue_scores <- bind_rows(all_tissue_scores_list)

# Save to CSV for future runs
write.csv(all_tissue_scores, cached_scores_file, row.names = FALSE)
message("Saved Factor1 scores to: ", cached_scores_file)
message("Total samples: ", nrow(all_tissue_scores))
message("Total tissues: ", length(unique(all_tissue_scores$Tissue)))

# ============================================================================
# APPLY BONFERRONI CORRECTION ACROSS ALL TESTS
# ============================================================================
message("\n========== APPLYING BONFERRONI CORRECTION ==========")

# Apply Bonferroni correction across all tests
results_df$P_Value_Adjusted <- p.adjust(results_df$P_Value_Raw, method = "bonferroni")

# Update asterisks based on adjusted p-values
results_df$Asterisks <- ifelse(results_df$P_Value_Adjusted < 0.0001, "****",
                                ifelse(results_df$P_Value_Adjusted < 0.001, "***",
                                       ifelse(results_df$P_Value_Adjusted < 0.01, "**",
                                              ifelse(results_df$P_Value_Adjusted < 0.05, "*", "ns"))))

message("Applied Bonferroni correction to ", nrow(results_df), " tests")
message("========================================\n")

} # End of if (!USE_CACHED_DATA)

# Final plot for Spearman rho
# One plot for chosen_Factor and one for GDF15
results_chosen_Factor <- results_df[results_df$Comparison == "chosen_Factor", ]
results_gdf15 <- results_df[results_df$Comparison == "GDF15", ]

p_final_chosen_Factor <- ggplot(results_chosen_Factor, aes(x = Spearman_Rho, y = reorder(Tissue, Spearman_Rho))) +
  geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper), height = 0.3, color = "gray50", linewidth = 0.5) +
  geom_point() +
  ggtitle("Spearman Rho for AGE vs", chosen_Factor) +
  theme_minimal()

p_final_gdf15 <- ggplot(results_gdf15, aes(x = Spearman_Rho, y = reorder(Tissue, Spearman_Rho))) +
  geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper), height = 0.3, color = "gray50", linewidth = 0.5) +
  geom_point() +
  ggtitle("Spearman Rho for AGE vs GDF15") +
  theme_minimal()

# print(p_final_chosen_Factor)
# print(p_final_gdf15)

ggsave(filename = paste0(folder_path, "/SpearmanR_AGEvs", chosen_Factor, ".png"), plot = p_final_chosen_Factor, width = 10, height = 8, dpi = 300)
ggsave(filename = paste0(folder_path, "/SpearmanR_AGEvsGDF15.png"), plot = p_final_gdf15, width = 10, height = 8, dpi = 300)



# Plot 1: Significant Spearman Rhos in red
plot_significant_rhos <- function(data, comparison) {
  ggplot(data, aes(x = Spearman_Rho, y = reorder(Tissue, Spearman_Rho), color = P_Value_Adjusted < 0.05)) +
    geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper), height = 0.3, color = "gray50", linewidth = 0.5) +
    geom_point(size = 3) +
    scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
    scale_x_continuous(limits = c(-0.6, 0.6)) +
    ggtitle(paste("Spearman Rho for AGE vs", comparison, "- Significant in Red (Bonferroni)")) +
    theme_minimal()
}

# chosen_Factor
p_significant_chosen_Factor <- plot_significant_rhos(results_chosen_Factor, chosen_Factor)
# GDF15
p_significant_gdf15 <- plot_significant_rhos(results_gdf15, "GDF15")

# print(p_significant_chosen_Factor)
# print(p_significant_gdf15)

# Plot 2: Asterisks next to datapoints for significant rhos
plot_with_asterisks <- function(data, comparison) {
  ggplot(data, aes(x = Spearman_Rho, y = reorder(Tissue, Spearman_Rho))) +
    geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper), height = 0.3, color = "gray50", linewidth = 0.5) +
    geom_point(size = 3) +
    geom_text(aes(label = Asterisks), vjust = -1, size = 5) +
    ggtitle(paste("Spearman Rho for AGE vs", comparison, "- Asterisks")) +
    theme_minimal()
}

# chosen_Factor with Asterisks
p_asterisks_chosen_Factor <- plot_with_asterisks(results_chosen_Factor, chosen_Factor)
# GDF15 with Asterisks
p_asterisks_gdf15 <- plot_with_asterisks(results_gdf15, "GDF15")

# print(p_asterisks_chosen_Factor)
# print(p_asterisks_gdf15)

# Plot 3: Combined Spearman Rho plot for chosen_Factor and GDF15
combined_data <- merge(results_chosen_Factor, results_gdf15, by = "Tissue", suffixes = c(chosen_Factor, "_GDF15"))

colnames(combined_data)[3] <- "Spearman_Rho_chosen_Factor"

p_combined <- ggplot(combined_data) +
  geom_errorbarh(aes(xmin = CI_lower_GDF15, xmax = CI_upper_GDF15, y = reorder(Tissue, Spearman_Rho_GDF15)), height = 0.2, color = "pink", alpha = 0.5, linewidth = 0.5) +
  geom_errorbarh(aes(xmin = CI_lowerFactor1, xmax = CI_upperFactor1, y = reorder(Tissue, Spearman_Rho_GDF15)), height = 0.2, color = "gray", alpha = 0.5, linewidth = 0.5) +
  geom_point(aes(x = Spearman_Rho_GDF15, y = reorder(Tissue, Spearman_Rho_GDF15)), color = "pink", size = 3) +
  geom_point(aes(x = Spearman_Rho_chosen_Factor, y = reorder(Tissue, Spearman_Rho_GDF15)), color = "gray", size = 3) +
  ggtitle(paste("Combined Spearman Rho for AGE vs GDF15 (pink) and", chosen_Factor, "(gray)")) +
  theme_minimal()

print(p_combined)

# Plot 4: Combined Spearman Rho plot for chosen_Factor and GDF15
combined_data <- merge(results_chosen_Factor, results_gdf15, by = "Tissue", suffixes = c(chosen_Factor, "_GDF15"))

colnames(combined_data)[3] <- "Spearman_Rho_chosen_Factor"

p_combined2 <- ggplot(combined_data) +
  geom_errorbarh(aes(xmin = CI_lower_GDF15, xmax = CI_upper_GDF15, y = reorder(Tissue, Spearman_Rho_chosen_Factor)), height = 0.2, color = "pink", alpha = 0.5, linewidth = 0.5) +
  geom_errorbarh(aes(xmin = CI_lowerFactor1, xmax = CI_upperFactor1, y = reorder(Tissue, Spearman_Rho_chosen_Factor)), height = 0.2, color = "gray", alpha = 0.5, linewidth = 0.5) +
  geom_point(aes(x = Spearman_Rho_GDF15, y = reorder(Tissue, Spearman_Rho_chosen_Factor)), color = "pink", size = 3) +
  geom_point(aes(x = Spearman_Rho_chosen_Factor, y = reorder(Tissue, Spearman_Rho_chosen_Factor)), color = "gray", size = 3) +
  ggtitle(paste("Combined Spearman Rho for AGE vs GDF15 (pink) and", chosen_Factor, "(gray)")) +
  theme_minimal()

print(p_combined2)


ggsave(filename = paste0(folder_path, "/SpearmanR_p_significant_", chosen_Factor, ".png"), plot = p_significant_chosen_Factor, width = 10, height = 8, dpi = 300)
# ggsave(filename = paste0(folder_path, "/SpearmanR_p_significant_gdf15.png"), plot = p_significant_gdf15, width = 10, height = 8, dpi = 300)
# ggsave(filename = paste0(folder_path, "/SpearmanR_p_asterisks_", chosen_Factor, ".png"), plot = p_asterisks_chosen_Factor, width = 10, height = 8, dpi = 300)
# ggsave(filename = paste0(folder_path, "/SpearmanR_p_asterisks_gdf15.png"), plot = p_asterisks_gdf15, width = 10, height = 8, dpi = 300)
# ggsave(filename = paste0(folder_path, "/SpearmanR_p_combined.png"), plot = p_combined, width = 10, height = 8, dpi = 300)
# ggsave(filename = paste0(folder_path, "/SpearmanR_p_combined_factorOrder.png"), plot = p_combined2, width = 10, height = 8, dpi = 300)

# Save main figure to shared Output directory as Figure_3C.png
# Using the significant Spearman Rho plot (red = significant after Bonferroni)
ggsave(filename = paste0(shared_out_dir, "/Figure_3C.png"), plot = p_significant_chosen_Factor, width = 10, height = 8, dpi = 300)

# ===================== Save Figure 3C source data (Factor1 only) =====================
write.csv(
  results_chosen_Factor,
  file = file.path(folder_path, "Figure_3C_source_data.csv"),
  row.names = FALSE
)
message("Saved Figure 3C source data to: ", file.path(folder_path, "Figure_3C_source_data.csv"))

# ===================== Extended Data Figure 14: GDF15 vs AGE (separate CSV) =====================
write.csv(
  results_gdf15,
  file = file.path(folder_path, "Extended_Data_Figure_14.csv"),
  row.names = FALSE
)
message("Saved Extended Data Figure 14 (GDF15 vs AGE) to: ", file.path(folder_path, "Extended_Data_Figure_14.csv"))
