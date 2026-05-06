
#This code currently gets the log2 TMM expression of GDF15, then centers and scales it, and then compared the pc2 score (which also gets centerd and scaled) for each tissue

rm(list = ls())
#
# BiocManager::install("ComplexHeatmap")
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



Plot_Save = "ON"
#------------------------------------------------------
#------------------------------------------------------ NAMING FOLDER AFTER THE NAME OF THE SCRIPT

# Output directories
local_out_dir <- here("Results", "Figures", "Figure_3B_4A")
shared_out_dir <- here("Results", "Figures")
if (!dir.exists(local_out_dir)) dir.create(local_out_dir, recursive = TRUE)
if (!dir.exists(shared_out_dir)) dir.create(shared_out_dir, recursive = TRUE)
folder_path <- local_out_dir  # Keep folder_path for compatibility with rest of script

# ===================== Output folder for plot CSVs =====================
plot_csv_dir <- file.path(folder_path, "Plot_Data_CSVs")
if (!dir.exists(plot_csv_dir)) dir.create(plot_csv_dir, recursive = TRUE)

# ===================== Output folder for Understanding Sample Sizes =====================
sample_size_dir <- file.path(folder_path, "Understanding_Sample_Sizes")
if (!dir.exists(sample_size_dir)) dir.create(sample_size_dir, recursive = TRUE)

# Store sample sizes across tissues (filled inside loop)
all_tissue_sample_sizes <- data.frame(
  Tissue = character(),
  N_total_SAMPID = integer(),
  N_used_complete = integer(),
  stringsAsFactors = FALSE
)

# ===================== Tracking sample counts at each filtering stage =====================
# Stage 1: No filtering (raw from .gz files)
counts_no_filtering <- data.frame(
  Tissue = character(),
  N_Samples_Raw = integer(),
  stringsAsFactors = FALSE
)

# Stage 2: After RIN >= 6 filtering (but before N >= 20)
counts_after_RIN <- data.frame(
  Tissue = character(),
  N_Samples_After_RIN = integer(),
  stringsAsFactors = FALSE
)

# Stage 3: After both RIN >= 6 AND N >= 20 filtering
counts_final <- data.frame(
  Tissue = character(),
  N_Samples_Final = integer(),
  stringsAsFactors = FALSE
)



# Read necessary input files

script_path <- here("Main Figure Scripts", "Helper Scripts", "gtex", "Attibutes_Phenos_Merged_plus_COD.R")
source(script_path, local = TRUE)



# Path to the folder containing individual tissue folders
#########################################################################################################
individual_tissues_folder <- here("Data", "gtex", "All_Tissues_Indiv_Folders")
#########################################################################################################

# Debug: print the path to check if it is correct
print(individual_tissues_folder)

correlation_results <- data.frame(
  Tissue = character(),
  rho = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

# Initialize an empty dataframe to store the results
results_df <- data.frame(
  Tissue = character(),
  rho = numeric(),
  p_adjust = numeric(),
  asterisks = character(),
  stringsAsFactors = FALSE
)



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



# new_folder_path <- here("Results", "gtex", Folder_name, last_component)

# # Check if the folder exists
# if (!dir.exists(new_folder_path)) {
#   # If it doesn't exist, create the folder
#   dir.create(new_folder_path, recursive = TRUE)
#   cat("Folder created at:", new_folder_path, "\n")
# } else {
#   cat("Folder already exists at:", new_folder_path, "\n")
# }


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

  # ===================== STAGE 1: Count raw samples (NO filtering) =====================
  # This counts all samples in the .gz file before any filtering
  N_raw <- ncol(combined_gene_counts)  # All sample columns
  counts_no_filtering <- rbind(
    counts_no_filtering,
    data.frame(Tissue = last_component, N_Samples_Raw = N_raw, stringsAsFactors = FALSE)
  )

  annotations_sample_attributes_filtered <- annotations_merged %>%
    filter(SAMPID %in% colnames(combined_gene_counts))

  # ===================== STAGE 2: Count after RIN >= 6 filtering =====================
  # annotations_merged already has RIN >= 6 filter applied (from Attibutes_Phenos_Merged_plus_COD.R)
  N_after_RIN <- nrow(annotations_sample_attributes_filtered)
  counts_after_RIN <- rbind(
    counts_after_RIN,
    data.frame(Tissue = last_component, N_Samples_After_RIN = N_after_RIN, stringsAsFactors = FALSE)
  )

  #### TMM normalization
  ##########       Option 1

  dge <- DGEList(combined_gene_counts[, annotations_sample_attributes_filtered$SAMPID], group = factor(annotations_sample_attributes_filtered$COD))
  keep <- filterByExpr(dge)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge, method = "TMM")
  exprs <- cpm(dge, log = TRUE)


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


  loadings_readin <- here("Data", "Fibroblast_lifespan", "Processed", "Fibroblast_FactorAnalysis_FacLoads_12_factors.csv") # produced by Figure_1E_1F.R; includes gdf15
  loadings <- read.csv(loadings_readin)
  loadings <- loadings %>%
    select(X, Factor1)


    selected_loadings <- loadings %>%
      rename(Gene = X, Loading = Factor1) %>%
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



     boxplot_data <- final_data %>%
       select(SAMPID, Factor1, AGE, SEX, DTHPLCE)

     exprs_GDF15 <- exprs_ISR_scaled %>%
       select(SAMPID, GDF15)

     boxplot_data <- boxplot_data %>%
       left_join(exprs_GDF15, by = "SAMPID")

     ## 1. Get the tissue name
     Tissue <- last_component  # or use 'last_component' if that's the correct variable

     boxplot_data$Tissue <- Tissue

     # boxplot_data <- boxplot_data %>%
     #   filter(DTHPLCE %in% c("Hospital inpatient", "Emergency room"))

     # Reshape the dataframe
     boxplot_data <- pivot_longer(
       boxplot_data,
       cols = c(Factor1, GDF15),
       names_to = "fa_vs_gdf15",
       values_to = "value"
     )

     # Reshape data from long to wide format
     boxplot_data_wide <- boxplot_data %>%
       pivot_wider(
         names_from = fa_vs_gdf15,  # Columns to create from unique values in fa_vs_gdf15
         values_from = value         # Values to fill in the new columns
       )

     boxplot_data <- boxplot_data_wide

     # Convert Factor1 and GDF15 to numeric
     boxplot_data <- boxplot_data %>%
       mutate(
         Factor1 = as.numeric(as.character(Factor1)),
         GDF15 = as.numeric(as.character(GDF15))
       )

     # ===================== Ensure one row per SAMPID =====================
     boxplot_data <- boxplot_data %>% distinct(SAMPID, .keep_all = TRUE)

     # ===================== STAGE 3: Count final samples used (after all processing) =====================
     # This is what actually gets used in the analysis
     N_final <- nrow(boxplot_data)
     counts_final <- rbind(
       counts_final,
       data.frame(Tissue = Tissue, N_Samples_Final = N_final, stringsAsFactors = FALSE)
     )

     # ===================== Sample size counts (NOT by deathplace) =====================
     N_total_SAMPID <- boxplot_data %>% distinct(SAMPID) %>% nrow()

     N_used_complete <- boxplot_data %>%
       filter(!is.na(Factor1), !is.na(GDF15)) %>%
       distinct(SAMPID) %>%
       nrow()

     # store counts for this tissue
     all_tissue_sample_sizes <- rbind(
       all_tissue_sample_sizes,
       data.frame(
         Tissue = Tissue,
         N_total_SAMPID = N_total_SAMPID,
         N_used_complete = N_used_complete,
         stringsAsFactors = FALSE
       )
     )



     # boxplot_data$value <- as.numeric(boxplot_data$value)
     # class(boxplot_data$value)





  ##### where I pasted

     #
     # # Reshape data from long to wide format
     # boxplot_data_wide <- boxplot_data %>%
     #   pivot_wider(
     #     names_from = fa_vs_gdf15,  # Columns to create from unique values in fa_vs_gdf15
     #     values_from = value         # Values to fill in the new columns
     #   )
     #
     # boxplot_data <- boxplot_data_wide

     # Calculate Spearman correlations per Tissue
     correlation_results <- boxplot_data %>%
       group_by(Tissue) %>%
       summarise(
         rho = cor(Factor1, GDF15, method = "spearman", use = "complete.obs"),
         p_value = cor.test(Factor1, GDF15, method = "spearman", exact = FALSE)$p.value,
         n_obs = sum(!is.na(Factor1) & !is.na(GDF15))
       ) %>%
       ungroup() %>%
       mutate(
         z = atanh(rho),
         se_z = 1 / sqrt(pmax(n_obs - 3, 1)),
         CI_lower = tanh(z - 1.96 * se_z),
         CI_upper = tanh(z + 1.96 * se_z)
       ) %>%
       select(-z, -se_z)

     # Adjust p-values for multiple comparisons (Bonferroni correction with 44 tissues)
     correlation_results <- correlation_results %>%
       mutate(
         p_adjust = p_value * 44,          # Bonferroni correction
         p_adjust = ifelse(p_adjust > 1, 1, p_adjust),  # Cap at 1
         asterisks = case_when(
           p_adjust < 0.0001 ~ "****",
           p_adjust < 0.001  ~ "***",
           p_adjust < 0.01   ~ "**",
           p_adjust < 0.05   ~ "*",
           TRUE              ~ ""
         )
       )

     # Merge correlation results back to boxplot_data
     boxplot_data <- boxplot_data %>%
       left_join(correlation_results, by = "Tissue")




     # ===================== Save the exact data that goes into the plot =====================
     plot_data_out <- boxplot_data %>%
       mutate(
         Tissue = Tissue,
         N_total_SAMPID = N_total_SAMPID,
         N_used_complete = N_used_complete
       ) %>%
       select(SAMPID, Tissue, DTHPLCE, Factor1, GDF15,
              rho, p_value, p_adjust, asterisks,
              N_total_SAMPID, N_used_complete) %>%
       arrange(DTHPLCE)

     write.csv(
       plot_data_out,
       file = file.path(plot_csv_dir, paste0("Plot_Data_", gsub(" ", "_", Tissue), ".csv")),
       row.names = FALSE
     )





     # Define color mapping for DTHPLCE
     color_mapping <- c(
       "Hospital inpatient" = "orange",
       "Emergency room" = "maroon",
       "Other" = "gray"   # Assign specific color if needed
     )

     # Ensure 'DTHPLCE' is a factor with the correct levels
     boxplot_data$DTHPLCE <- factor(boxplot_data$DTHPLCE, levels = names(color_mapping))

     # Replace NA values in DTHPLCE with "Other"
     boxplot_data <- boxplot_data %>%
       mutate(
         DTHPLCE = replace_na(DTHPLCE, "Other")
       )
     sum(is.na(boxplot_data$DTHPLCE))  # Should be 0 or handled in color mapping

     sum(is.na(boxplot_data$Factor1))  # Count of NA in Factor1
     sum(is.na(boxplot_data$GDF15))    # Count of NA in GDF15
     boxplot_data$Factor1 <- as.numeric(as.character(boxplot_data$Factor1))
     boxplot_data$GDF15 <- as.numeric(as.character(boxplot_data$GDF15))


       # Retrieve correlation results
       rho_value <- unique(boxplot_data$rho)
       asterisks <- unique(boxplot_data$asterisks)

       # Handle potential multiple entries
       rho_value <- ifelse(length(rho_value) > 1, rho_value[1], rho_value)
       asterisks <- ifelse(length(asterisks) > 1, asterisks[1], asterisks)

       # Base plot (no annotation)
       p_base <- ggplot(boxplot_data, aes(x = Factor1, y = GDF15, color = DTHPLCE)) +
         geom_point(size = 3, alpha = 0.7) +  # Scatter points
         geom_smooth(method = "lm", se = TRUE, color = "red") +
         ggtitle(Tissue) +
         theme_minimal() +
         theme(
           plot.title = element_text(size = 35, face = "bold", hjust = 0.5),
           axis.title = element_text(size = 30),
           axis.text  = element_text(size = 28),
           legend.position = "none"
         ) +
         labs(
           x = "Index Score",
           y = "GDF15"
         ) +
         # Apply custom color mapping
         scale_color_manual(values = color_mapping)

       # Version WITH Spearman rho and p-value (top right)
       p <- p_base +
         annotate("text",
                  x = Inf, y = Inf,
                  label = paste0(
                    "\u03c1 = ", round(unique(boxplot_data$rho)[1], 3),
                    "\np = ",    signif(unique(boxplot_data$p_value)[1], 3)
                  ),
                  hjust = 1.1, vjust = 1.5,
                  size = 10, color = "black")

       # Display the annotated plot
       print(p)

       save_path <- paste0(folder_path, "/")
       # Create the directory if it doesn't exist
       if (!dir.exists(save_path)) {
         dir.create(save_path, recursive = TRUE)
         cat("Directory created at:", save_path, "\n")
       } else {
         cat("Directory already exists at:", save_path, "\n")
       }

       # Subfolder for plots without annotation
       save_path_no_rho <- file.path(save_path, "No_Rho_Label")
       if (!dir.exists(save_path_no_rho)) dir.create(save_path_no_rho, recursive = TRUE)

       # Save annotated version (with rho)
       ggsave(
         filename = paste0("Spearman_Correlation_", gsub(" ", "_", Tissue), ".png"),
         plot = p,
         path = save_path,
         width = 8,
         height = 6,
         dpi = 300
       )

       # Save version WITHOUT rho annotation
       ggsave(
         filename = paste0("Spearman_Correlation_", gsub(" ", "_", Tissue), ".png"),
         plot = p_base,
         path = save_path_no_rho,
         width = 8,
         height = 6,
         dpi = 300
       )



       # Extract rho and p_adjust from boxplot_data
       rho_value <- unique(boxplot_data$rho)
       p_adjust_value <- unique(boxplot_data$p_adjust)

       # Assign asterisks based on adjusted p-value
       asterisks <- unique(boxplot_data$asterisks)

       # Create a temporary dataframe for the current Tissue
       temp_df <- data.frame(
         Tissue = Tissue,
         rho = rho_value,
         p_adjust = p_adjust_value,
         asterisks = asterisks,
         n_obs = unique(boxplot_data$n_obs)[1],
         CI_lower = unique(boxplot_data$CI_lower)[1],
         CI_upper = unique(boxplot_data$CI_upper)[1],
         stringsAsFactors = FALSE
       )

       # Append the temporary dataframe to the results dataframe
       results_df <- rbind(results_df, temp_df)



       cat("Plot saved for tissue:", Tissue, "\n")
     }




# Optional: Order tissues by rho for better visualization
results_df$Tissue <- factor(results_df$Tissue, levels = results_df$Tissue[order(results_df$rho)])

# Create the correlation plot
correlation_plot <- ggplot(results_df, aes(x = rho, y = Tissue)) +
  geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper), height = 0.3, color = "gray50", linewidth = 0.5) +
  geom_point(color = "black",
             size = 4,
             alpha = 0.7) +  # Semi-transparent points
  # geom_vline(xintercept = 0, linetype = "dashed", color = "black") +  # Reference line at rho = 0
  # geom_text(aes(label = asterisks),
  #           hjust = -0.3,  # Adjust horizontal position
  #           vjust = 0.5,   # Adjust vertical position
  #           size = 6,
  #           color = "black") +  # Add asterisks next to significant points
  expand_limits(x = 0) +  # Ensure x-axis includes
  labs(
    title = "Spearman Correlation between Factor1 and GDF15 by Tissue",
    x = "Spearman Rho",
    y = "Tissue"
    # color = "Tissue"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.position = "none"  # Hide legend if colors are self-explanatory
  )

# Display the plot
print(correlation_plot)

# Define the directory where the summary plot will be saved
save_path2 <- paste0(folder_path, "/")

# Save the summary plot
ggsave("Spearman_Correlation_Summary_x0.png", plot = correlation_plot, path = save_path2, width = 10, height = 8, dpi = 300)

# Also save to shared output directory
ggsave("Figure_3B.png", plot = correlation_plot, path = shared_out_dir, width = 10, height = 8, dpi = 300)


# ===================== Save all-tissue sample sizes =====================
all_tissue_sample_sizes <- all_tissue_sample_sizes %>%
  distinct(Tissue, .keep_all = TRUE) %>%
  arrange(Tissue)

write.csv(
  all_tissue_sample_sizes,
  file = file.path(plot_csv_dir, "All_Tissues_Sample_sizes.csv"),
  row.names = FALSE
)

# ===================== Save Understanding Sample Sizes CSVs =====================

# Stage 1: No filtering
write.csv(
  counts_no_filtering,
  file = file.path(sample_size_dir, "No_Filtering_Sample_Sizes_Per_Tissue.csv"),
  row.names = FALSE
)
message("Saved: No_Filtering_Sample_Sizes_Per_Tissue.csv")

# Stage 2: After RIN >= 6 filtering (but before N >= 20)
write.csv(
  counts_after_RIN,
  file = file.path(sample_size_dir, "Filtered_RIN_but_not_samp_size_20.csv"),
  row.names = FALSE
)
message("Saved: Filtered_RIN_but_not_samp_size_20.csv")

# Stage 3: Final samples used in analysis
write.csv(
  counts_final,
  file = file.path(sample_size_dir, "Final_Samples_After_All_Processing.csv"),
  row.names = FALSE
)
message("Saved: Final_Samples_After_All_Processing.csv")

# ===================== Create a combined summary table =====================
sample_size_summary <- counts_no_filtering %>%
  full_join(counts_after_RIN, by = "Tissue") %>%
  full_join(counts_final, by = "Tissue") %>%
  mutate(
    Lost_to_RIN_Filter = N_Samples_Raw - N_Samples_After_RIN,
    Lost_to_Processing = N_Samples_After_RIN - N_Samples_Final,
    Pct_Retained = round(100 * N_Samples_Final / N_Samples_Raw, 1)
  ) %>%
  arrange(Tissue)

write.csv(
  sample_size_summary,
  file = file.path(sample_size_dir, "Sample_Size_Summary_All_Stages.csv"),
  row.names = FALSE
)
message("Saved: Sample_Size_Summary_All_Stages.csv")

# Spearman_Rho_All_Tissues.csv removed (duplicate of Figure_3B_source_data.csv below)

# Print summary for kidney_cortex specifically
kidney_info <- sample_size_summary %>% filter(Tissue == "kidney_cortex")
if (nrow(kidney_info) > 0) {
  message("\n=== Kidney Cortex Sample Counts ===")
  message("Raw (no filtering): ", kidney_info$N_Samples_Raw)
  message("After RIN >= 6: ", kidney_info$N_Samples_After_RIN)
  message("Final (used in analysis): ", kidney_info$N_Samples_Final)
  message("Lost to RIN filter: ", kidney_info$Lost_to_RIN_Filter)
  message("Lost to other processing: ", kidney_info$Lost_to_Processing)
}

# ===================== Save Figure 3B source data =====================
write.csv(
  results_df,
  file = file.path(folder_path, "Figure_3B_source_data.csv"),
  row.names = FALSE
)
message("Saved Figure 3B source data to: ", file.path(folder_path, "Figure_3B_source_data.csv"))

# ============================================================================
# EXPORT FIGURE 4A SOURCE DATA: 4 tissues with DTHPLCE (incl. "Other")
# ============================================================================
message("\n=== Exporting Figure 4A source data (4 tissues with DTHPLCE) ===")

fig4a_tissues <- c("skin_lower_leg", "liver", "heart_left_ventricle", "brain_frontal_cortex")
fig4a_dthplce_keep <- c("Hospital inpatient", "Emergency room", "Other",
                         "Hospital Inpatient", "Emergency Room")  # handle case variants

plot_csv_dir <- file.path(folder_path, "Plot_Data_CSVs")

for (tissue_name in fig4a_tissues) {
  tissue_csv <- file.path(plot_csv_dir, paste0("Plot_Data_", tissue_name, ".csv"))
  if (file.exists(tissue_csv)) {
    tissue_df <- read.csv(tissue_csv, check.names = FALSE)
    # Filter to Hospital Inpatient, Emergency Room, and Other
    tissue_df_filtered <- tissue_df[tissue_df$DTHPLCE %in% fig4a_dthplce_keep, ]
    # Drop SAMPID for privacy (keep Tissue, DTHPLCE, Factor1, GDF15, rho, etc.)
    cols_to_drop <- intersect(c("SAMPID", "SUBJID"), colnames(tissue_df_filtered))
    if (length(cols_to_drop) > 0) {
      tissue_df_filtered <- tissue_df_filtered[, !(colnames(tissue_df_filtered) %in% cols_to_drop), drop = FALSE]
    }
    out_file <- file.path(folder_path, paste0("Figure_4A_", tissue_name, "_source_data.csv"))
    write.csv(tissue_df_filtered, out_file, row.names = FALSE)
    message("  Saved: Figure_4A_", tissue_name, "_source_data.csv (", nrow(tissue_df_filtered), " rows)")
  } else {
    message("  WARNING: Plot_Data_", tissue_name, ".csv not found in Plot_Data_CSVs/")
  }
}
