

#This code currently gets the log2 TMM expression of GDF15, then centers and scales it, and then compared the pc2 score (which also gets centerd and scaled) for each tissue

# Control F and then replace "Factor1" with Factor_
rm(list = ls())

# BiocManager::install("edgeR")
# install.packages("readxl")

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
library(readxl)
library(here)

# Output directories
local_out_dir <- here("Results", "Figures", "Figure_2C")
shared_out_dir <- here("Results", "Figures")
if (!dir.exists(local_out_dir)) dir.create(local_out_dir, recursive = TRUE)
if (!dir.exists(shared_out_dir)) dir.create(shared_out_dir, recursive = TRUE)
folder_path <- local_out_dir  # Keep folder_path for compatibility with rest of script

Opto_Readin <- here("Data", "Opto_PKR", "Taivan_Normalized_counts_text.xlsx")
Opto_data <- read_excel(Opto_Readin)

correlations_readin <- here("Data", "Fibroblast_lifespan", "Processed", "Fibroblast_FactorAnalysis_correlations_for_scaling.csv")
correlations <- read.csv(correlations_readin)

ISR_read_in <- here("Data", "Opto_PKR", "Total_ISR_Gene_List_plus_gdf15.csv")
Total_ISR_List <- read.csv(ISR_read_in) %>%
  select(-X)



### Converting FPKM to TPM
colnames(Opto_data)
# Rename gene name column for clarity
colnames(Opto_data)[1] <- "X"


# Convert all columns *except* the first to numeric
Opto_data[, -1] <- lapply(Opto_data[, -1], as.numeric)


# Function to convert FPKM to TPM for a single sample (column)

gene_names <- Opto_data$X # 1. Save gene names
fpkm_only <- Opto_data[, -1] # 2. Convert all numeric columns to TPM (skip column 1)
fpkm_only <- as.data.frame(lapply(fpkm_only, as.numeric)) # Ensure all numeric values

# Define function
fpkm_to_tpm <- function(fpkm_vector) {
  tpm <- (fpkm_vector / sum(fpkm_vector, na.rm = TRUE)) * 1e6
  return(tpm)
}
tpm_only <- apply(fpkm_only, 2, fpkm_to_tpm) # Apply TPM conversion to each column
tpm_only <- as.data.frame(tpm_only)# Convert back to data frame
tpm_log <- log2(tpm_only + 1) # Log2-transform TPM
Opto_data_TPMlog <- cbind(X = gene_names, tpm_log) # 3. Add gene names back as the first column
Opto_data <- Opto_data_TPMlog


# Get the list of genes from Opto_data
opto_genes <- Opto_data$X

# Get the list of genes from Total_ISR_List
isr_genes <- Total_ISR_List$Gene

# Check which genes are missing in Opto_data
missing_genes <- isr_genes[!isr_genes %in% opto_genes]

# Print the missing genes
if (length(missing_genes) > 0) {
  cat("Missing genes:\n")
  print(missing_genes)
} else {
  cat("All genes from Total_ISR_List are present in Opto_data.\n")
}

Opto_data <- Opto_data %>%
  mutate(X = case_when(
    X == "KIAA0141" ~ "DELE1",
    X == "WARS"     ~ "WARS1",
    X == "NARS"     ~ "NARS1",
    X == "ZAK"      ~ "MAP3K20",
    X == "ERO1L"    ~ "ERO1A",
    X == "GCN1L1"    ~ "GCN1",
    TRUE ~ X
  ))

# Extract the gene column from Opto_data
opto_genes <- Opto_data$X

# Find the genes still missing after renaming
missing_genes <- Total_ISR_List$Gene[!Total_ISR_List$Gene %in% opto_genes]

# Create a dataframe of 0s for missing genes
missing_df <- data.frame(
  X = missing_genes,
  matrix(0, nrow = length(missing_genes), ncol = ncol(Opto_data) - 1)
)

# Set column names to match Opto_data
colnames(missing_df) <- colnames(Opto_data)

# Combine the original and new data
Opto_data_updated <- bind_rows(Opto_data, missing_df)

Opto_data_filtered <- Opto_data_updated %>%
  filter(X %in% Total_ISR_List$Gene)

Opto_data_filtered <- Opto_data_filtered %>%
  rename(
    Gene = X)

# Make sure the gene column is character (not factor)
Opto_data_filtered$Gene <- as.character(Opto_data_filtered$Gene)

# Set row names to gene names and remove the Gene column
Opto_matrix <- Opto_data_filtered %>%
  column_to_rownames(var = "Gene")

# Transpose the dataframe
Opto_transposed <- as.data.frame(t(Opto_matrix))

# Optional: move rownames (timepoints) into a column
Opto_transposed <- Opto_transposed %>%
  rownames_to_column(var = "Timepoint")

exprs_ISR <- Opto_transposed
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


# Step 1: Extract gene names from exprs_ISR (excluding Timepoint)
gene_columns <- setdiff(names(exprs_ISR), "Timepoint")

# Step 2: Ensure that gene names match between exprs_ISR and sd_x_fbdata
common_genes <- intersect(gene_columns, sd_x_fbdata$Gene)
# view(common_genes)


#Finding any missing genes
# Step 1: Extract gene names
exprs_genes <- setdiff(names(exprs_ISR), "Timepoint")
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
exprs_ISR_subset <- exprs_ISR[, c("Timepoint", common_genes)]

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


# Step 5: Manually scale the data # # # HERE IS THE DATA!!!# # #

exprs_ISR_scaled <- exprs_ISR_subset

for (gene in common_genes) {
  exprs_ISR_scaled[[gene]] <- (exprs_ISR_subset[[gene]] - mean_values[gene]) / sd_values[gene]
}


exprs_gdf15 <- exprs_ISR %>%
  select(Timepoint, GDF15)

# exprs_GADD34 <- exprs_ISR_scaled %>%
#   select(Timepoint, PPP1R15A)

########################################################################################

scaled_exprs_ISR <- exprs_ISR_scaled

loadings_readin <- here("Data", "Fibroblast_lifespan", "Processed", "Fibroblast_FactorAnalysis_FacLoads_12_factors.csv") # produced by Figure_1E_1F.R; includes gdf15
loadings <- read.csv(loadings_readin)

loadings <- loadings %>%
  select(X, Factor1)


# Get all factor columns from the loadings dataframe
factor_columns <- grep("^Factor\\d+$", colnames(loadings), value = TRUE)

# Prepare correlation matrix once
rownames(correlations) <- correlations$X
correlations$X <- NULL
correlations_matrix <- as.matrix(correlations)

# Set up scaled expression matrix
common_genes <- intersect(loadings$X, colnames(scaled_exprs_ISR))
scaled_data <- scaled_exprs_ISR %>% select(Timepoint, all_of(common_genes))
rownames(scaled_data) <- scaled_data$Timepoint
scaled_data$Timepoint <- NULL
scaled_data_matrix <- as.matrix(scaled_data)

###### ONLY FACTOR1 RIGHT NOw!!!!!!!!!

# loadings <- loadings %>%
#   select(X, Factor1)

# Loop through each factor
for (factor_name in factor_columns) {

  message("Processing ", factor_name, "...")

  # Prepare selected loadings
  selected_loadings <- loadings %>%
    select(Gene = X, Loading = all_of(factor_name)) %>%
    filter(Gene %in% common_genes)

  # Create loading matrix
  loadings_matrix <- as.matrix(selected_loadings$Loading)
  rownames(loadings_matrix) <- selected_loadings$Gene

  # Compute factor scores
  W <- t(loadings_matrix) %*% solve(correlations_matrix[selected_loadings$Gene, selected_loadings$Gene])
  Wtransposed <- t(W)
  final <- scaled_data_matrix[, selected_loadings$Gene] %*% Wtransposed

  # Final dataframe
  final_df <- data.frame(Timepoint = rownames(scaled_data_matrix), Score = final[,1])
  colnames(final_df)[2] <- factor_name

  # Merge with GDF15
  merged_df <- final_df %>%
    left_join(exprs_gdf15 %>% select(Timepoint, GDF15), by = "Timepoint") %>%
    pivot_longer(cols = c(all_of(factor_name), "GDF15"), names_to = "fa_vs_gdf15", values_to = "value") %>%
    mutate(
      type = ifelse(fa_vs_gdf15 == factor_name, "fa_Scores", "GDF15_scaled"),
      Time_numeric = as.numeric(gsub("X(\\d+)_.*", "\\1", Timepoint))
    )
}


library(ggplot2)
library(dplyr)

# Step 1: Filter Factor1 and GDF15 rows
df_factor <- merged_df %>% filter(fa_vs_gdf15 == "Factor1")
df_gdf15 <- merged_df %>% filter(fa_vs_gdf15 == "GDF15")

# Step 2: Convert time from minutes to hours
df_factor <- df_factor %>% mutate(Time_hours = Time_numeric / 60)
df_gdf15 <- df_gdf15 %>% mutate(Time_hours = Time_numeric / 60)

# Step 3: Calculate scale factor for GDF15 to match Factor1 visually
scale_factor <- max(df_gdf15$value, na.rm = TRUE) / max(df_factor$value, na.rm = TRUE)


# Step 4: Create the dual-axis plot with points and lines
combined_plot <- ggplot() +
  geom_line(data = df_factor, aes(x = Time_hours, y = value, color = "Factor1"), size = 1.2) +
  geom_point(data = df_factor, aes(x = Time_hours, y = value, color = "Factor1"), size = 4) +
  geom_line(data = df_gdf15, aes(x = Time_hours, y = value / scale_factor, color = "GDF15"), linewidth = 1.2) +
  geom_point(data = df_gdf15, aes(x = Time_hours, y = value / scale_factor, color = "GDF15"), size = 4) +
  scale_x_continuous(
    name = "Time (hours)",
    limits = c(0, 12),
    breaks = seq(0, 12, by = 2)
  ) +
  scale_y_continuous(
    name = "Factor1 value",
    sec.axis = sec_axis(~ . * scale_factor, name = "GDF15 value")
  ) +
  scale_color_manual(
    name = "Legend",
    values = c("Factor1" = "black", "GDF15" = "pink")
  ) +
  labs(title = "Factor1 and GDF15 Over Time") +
  theme_minimal() +
  theme(
    axis.title.y.left = element_text(color = "black"),
    axis.text.y.left = element_text(color = "black"),
    axis.title.y.right = element_text(color = "pink"),
    axis.text.y.right = element_text(color = "pink"),
    legend.position = "top"
  )
#
# print(combined_plot)
#
#
# ggsave(file.path(folder_path, paste0("plot_combined_Factor1_GDF15.png")),
#        plot = combined_plot, width = 7, height = 5, dpi = 300)




library(ggplot2)
library(dplyr)

# Ranges
FA_min  <- -1.5
FA_max  <- 1.5
GDF_min <- 8
GDF_max <- 11

# Compute slope (a) and intercept (b) for GDF15 vs Factor1
a <- (GDF_max - GDF_min) / (FA_max - FA_min)  # should be 1.6
b <- GDF_min - a * FA_min                    # should be 8400

combined_plot <- ggplot() +
  # Plot Factor1 in its actual values
  geom_line(
    data = df_factor,
    aes(x = Time_hours, y = value, color = "Factor1"),
    size = 1.2
  ) +
  geom_point(
    data = df_factor,
    aes(x = Time_hours, y = value, color = "Factor1"),
    size = 4
  ) +
  # Plot GDF15 in Factor1-scale => (GDF15 - b) / a
  geom_line(
    data = df_gdf15,
    aes(x = Time_hours, y = (value - b) / a, color = "GDF15"),
    linewidth = 1.2
  ) +
  geom_point(
    data = df_gdf15,
    aes(x = Time_hours, y = (value - b) / a, color = "GDF15"),
    size = 4
  ) +
  # X-axis from 0 to 12 hours
  scale_x_continuous(
    name = "Time (hours)",
    limits = c(0, 12),
    breaks = seq(0, 12, by = 2)
  ) +
  # Y-axis for Factor1 from -4000 to 6000
  # Right axis transforms Factor1-range to GDF15-range
  scale_y_continuous(
    name = "Factor1 value",
    limits = c(FA_min, FA_max),
    # You can adjust the left-axis breaks to suit
    breaks = seq(FA_min, FA_max, by = 0.5),
    sec.axis = sec_axis(
      trans = ~ . * a + b,
      name = "GDF15 value",
      # Similarly, pick breaks for the right axis if you'd like
      breaks = seq(GDF_min, GDF_max, by = 0.5)
    )
  ) +
  scale_color_manual(
    name = "Legend",
    values = c("Factor1" = "black", "GDF15" = "pink")
  ) +
  labs(title = "ISRgdf15 kinetics with PKR activation") +
  theme_minimal() +
  theme(
    axis.title.y.left   = element_text(color = "black"),
    axis.text.y.left    = element_text(color = "black"),
    axis.title.y.right  = element_text(color = "pink"),
    axis.text.y.right   = element_text(color = "pink"),
    legend.position     = "top"
  )

print(combined_plot)



ggsave(file.path(folder_path, paste0("Figure_2C_ISRgdf15_kinetics_with_PKR_activation.png")),
       plot = combined_plot, width = 7, height = 5, dpi = 300)

# Save to shared output directory
ggsave(file.path(shared_out_dir, paste0("Figure_2C.png")),
       plot = combined_plot, width = 7, height = 5, dpi = 300)




###### Checking other genes ######

chosen_gene <- "EIF2AK2" # PKR
# chosen_gene <- "GDF15"

Gene_check <- Opto_data

Gene_check <- Gene_check %>%
dplyr::filter(X == chosen_gene)

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)

Gene_check_num <- Gene_check %>%
  rename_with(~ as.character(readr::parse_number(.x)), .cols = -X)

plot_df <- Gene_check_num %>%
  pivot_longer(cols = -X, names_to = "Time_min", values_to = "value") %>%
  mutate(Time_min = as.numeric(Time_min)) %>%
  arrange(Time_min)

plot_df <- plot_df %>%
  mutate(Time_hours = Time_min/60)

p_individual <- ggplot(plot_df, aes(Time_hours, value)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 4) +
  scale_x_continuous("Time (hours)", breaks = pretty(plot_df$Time_hours)) +
  labs(title = unique(Gene_check$X), y = "Gene Expression") +
  theme_minimal(base_size = 18)

print(p_individual)


ggsave(file.path(folder_path, paste0(chosen_gene, "_score.png")),
       plot = p_individual, width = 7, height = 5, dpi = 300)

# ================================================================================================
# Export source data for Figure 2C
# ================================================================================================

# Combine both data frames for the plotted data (Factor1 and GDF15 kinetics)
source_data_2C <- bind_rows(
  df_factor %>% select(Time_hours, value, fa_vs_gdf15),
  df_gdf15 %>% select(Time_hours, value, fa_vs_gdf15)
)

# Save to CSV
write.csv(source_data_2C,
          file = file.path(folder_path, "Figure_2C_source_data.csv"),
          row.names = FALSE)

message("Figure 2C source data saved to: ", file.path(folder_path, "Figure_2C_source_data.csv"))
