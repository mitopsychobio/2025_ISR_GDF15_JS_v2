
# Control F and then replace "Factor1" with Factor_
#Control F and replace the chosen gene name but make sure that the gene list csv stays with EIF2AK3 or it wont run
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

chosen_gene_to_compare <- "EIF2AK3"

# Output directories
local_out_dir <- here("Results", "Figures", "Figure_2F")
shared_out_dir <- here("Results", "Figures")
if (!dir.exists(local_out_dir)) dir.create(local_out_dir, recursive = TRUE)
if (!dir.exists(shared_out_dir)) dir.create(shared_out_dir, recursive = TRUE)
folder_path <- local_out_dir  # Keep folder_path for compatibility with rest of script

# Readin <- here("Data", "NatComms_K_Labbe_et_al", "Taivan_Normalized_counts_text.xlsx")
# TPM_data <- read_excel(Readin)
tpm_path <- here("Data", "NatComms_K_Labbe_et_al", "GSE273600_Control_Conditions_TPMs.txt.gz")
TPM_data <- readr::read_tsv(
  file = tpm_path,
  guess_max = 200000,          # increase if needed
  show_col_types = FALSE,
  progress = TRUE
)

correlations_readin <- here("Data", "Fibroblast_lifespan", "Processed", "Fibroblast_FactorAnalysis_correlations_for_scaling.csv")
correlations <- read.csv(correlations_readin)

ISR_read_in <- here("Data", "NatComms_K_Labbe_et_al", "Total_ISR_Gene_List_plus_GDF15.csv")
Total_ISR_List <- read.csv(ISR_read_in) %>%
  select(-X)

# THIS FILE ALREADY CONTAINS TPM #

# ### Converting FPKM to TPM
# colnames(TPM_data)
# # Rename gene name column for clarity
# colnames(TPM_data)[1] <- "X"
#
#
# # Convert all columns *except* the first to numeric
# TPM_data[, -1] <- lapply(TPM_data[, -1], as.numeric)
#
#
# # Function to convert FPKM to TPM for a single sample (column)
#
# gene_names <- TPM_data$X # 1. Save gene names
# fpkm_only <- TPM_data[, -1] # 2. Convert all numeric columns to TPM (skip column 1)
# fpkm_only <- as.data.frame(lapply(fpkm_only, as.numeric)) # Ensure all numeric values
#
# # Define function
# fpkm_to_tpm <- function(fpkm_vector) {
#   tpm <- (fpkm_vector / sum(fpkm_vector, na.rm = TRUE)) * 1e6
#   return(tpm)
# }
# tpm_only <- apply(fpkm_only, 2, fpkm_to_tpm) # Apply TPM conversion to each column
# tpm_only <- as.data.frame(tpm_only)# Convert back to data frame
# tpm_log <- log2(tpm_only + 1) # Log2-transform TPM
# TPM_data_TPMlog <- cbind(X = gene_names, tpm_log) # 3. Add gene names back as the first column
# TPM_data <- TPM_data_TPMlog

TPM_log2 <- TPM_data %>%
  mutate(
    across(
      -c(ID, gene_name),          # all columns except the first two
      ~ log2(as.numeric(.x) + 1)       # +1 to handle zeros safely
    )
  )

TPM_data <-TPM_log2

TPM_data <- TPM_data %>%
  select(-ID)

TPM_data <- TPM_data %>% rename(X = gene_name)
# Get the list of genes from TPM_data
all_genes <- TPM_data$X

# Get the list of genes from Total_ISR_List
isr_genes <- Total_ISR_List$Gene

# Check which genes are missing in TPM_data
missing_genes <- isr_genes[!isr_genes %in% all_genes]

# Print the missing genes
if (length(missing_genes) > 0) {
  cat("Missing genes:\n")
  print(missing_genes)
} else {
  cat("All genes from Total_ISR_List are present in TPM_data.\n")
}

TPM_data <- TPM_data %>%
  mutate(X = case_when(
    X == "KIAA0141" ~ "DELE1",
    X == "WARS"     ~ "WARS1",
    X == "NARS"     ~ "NARS1",
    X == "ZAK"      ~ "MAP3K20",
    X == "ERO1L"    ~ "ERO1A",
    X == "GCN1L1"    ~ "GCN1",
    X == "MLTK"      ~ "MAP3K20",
    X == "MLT"      ~ "MAP3K20",
    X == "MLK"      ~ "MAP3K20",
    TRUE ~ X
  ))

# Extract the gene column from TPM_data
all_genes <- TPM_data$X

# Find the genes still missing after renaming
missing_genes <- Total_ISR_List$Gene[!Total_ISR_List$Gene %in% all_genes]

# Create a dataframe of 0s for missing genes
missing_df <- data.frame(
  X = missing_genes,
  matrix(0, nrow = length(missing_genes), ncol = ncol(TPM_data) -1)
)

# Set column names to match TPM_data
colnames(missing_df) <- colnames(TPM_data)

# Combine the original and new data
TPM_data_updated <- bind_rows(TPM_data, missing_df)

TPM_data_filtered <- TPM_data_updated %>%
  filter(X %in% Total_ISR_List$Gene)

TPM_data_filtered <- TPM_data_filtered %>%
  rename(
    Gene = X)

# Make sure the gene column is character (not factor)
TPM_data_filtered$Gene <- as.character(TPM_data_filtered$Gene)

# Set row names to gene names and remove the Gene column
Opto_matrix <- TPM_data_filtered %>%
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



exprs_EIF2AK3 <- exprs_ISR %>%
  select(Timepoint, EIF2AK3)

# exprs_GADD34 <- exprs_ISR_scaled %>%
#   select(Timepoint, EIF2AK3)

########################################################################################

scaled_exprs_ISR <- exprs_ISR_scaled

loadings_readin <- here("Data", "Fibroblast_lifespan", "Processed", "Fibroblast_FactorAnalysis_FacLoads_12_factors.csv") # produced by Figure_1E_1F.R; includes EIF2AK3
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

  # Merge with EIF2AK3
  merged_df <- final_df %>%
    left_join(exprs_EIF2AK3 %>% select(Timepoint, EIF2AK3), by = "Timepoint") %>%
    pivot_longer(cols = c(all_of(factor_name), "EIF2AK3"), names_to = "fa_vs_EIF2AK3", values_to = "value") %>%
    mutate(
      type = ifelse(fa_vs_EIF2AK3 == factor_name, "fa_Scores", "EIF2AK3_scaled"),
      Time_numeric = as.numeric(gsub("X(\\d+)_.*", "\\1", Timepoint))
    )
}

# =============================== Setup =======================================
suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
})

# ====================== Parse columns from Timepoint ==========================
merged_df <- merged_df %>%
  tidyr::extract(
    col   = Timepoint,
    into  = c("Parental", "Condition", "Hours", "Rep"),
    regex = "^([^_]+)_([^_]+)_(\\d+)h_(rep\\d+)$",
    remove = FALSE
  ) %>%
  mutate(Hours = as.numeric(Hours))

# ====================== Summary stats (mean/sem) ==============================
summ_df <- merged_df %>%
  group_by(Parental, Condition, Hours, fa_vs_EIF2AK3) %>%
  summarise(
    n   = dplyr::n(),
    mean = mean(value, na.rm = TRUE),
    sd   = sd(value,   na.rm = TRUE),
    sem  = sd / sqrt(n),
    .groups = "drop"
  )

# ====================== Harmonize condition names =============================
normalize_condition <- function(x) {
  case_when(
    x %in% c("UT") ~ "UT",
    x %in% c("Dim") ~ "Dim",
    x == "Ars" ~ "Ars",
    x == "Tg"  ~ "Tg",
    TRUE ~ x
  )
}

summ_df <- summ_df %>%
  mutate(Condition = normalize_condition(Condition))

# Keep only the 4 conditions of interest (if present)
conditions_wanted <- c("UT","Dim","Ars","Tg")
summ_df <- summ_df %>%
  filter(Condition %in% conditions_wanted) %>%
  mutate(Condition = factor(Condition, levels = conditions_wanted[conditions_wanted %in% Condition]))



# ── Add 0h for all Conditions based on UT 0h within each Parental x analyte ──
# Get the UT 0h baselines per Parental x analyte
zero_ut <- summ_df %>%
  dplyr::filter(Condition == "UT", Hours == 0) %>%
  dplyr::select(Parental, fa_vs_EIF2AK3, n, mean, sd, sem)

# For every non-UT Condition, create a 0h row using the matching UT 0h baseline
zero_rows_all_cond <- summ_df %>%
  dplyr::filter(Condition != "UT") %>%                     # only non-UT conditions
  dplyr::distinct(Parental, Condition, fa_vs_EIF2AK3) %>%    # one row per group
  dplyr::left_join(zero_ut, by = c("Parental", "fa_vs_EIF2AK3")) %>%
  dplyr::mutate(Hours = 0L)

# Bind back and make sure we don't duplicate any existing rows
summ_df <- dplyr::bind_rows(summ_df, zero_rows_all_cond) %>%
  dplyr::arrange(Parental, Condition, fa_vs_EIF2AK3, Hours) %>%
  dplyr::distinct(Parental, Condition, Hours, fa_vs_EIF2AK3, .keep_all = TRUE)




# ====================== Global y-limits (shared per analyte) ==================
ylims_factor1 <- summ_df %>%
  filter(fa_vs_EIF2AK3 == "Factor1") %>%
  summarise(lo = min(mean, na.rm = TRUE), hi = max(mean, na.rm = TRUE)) %>%
  as.list()

ylims_EIF2AK3 <- summ_df %>%
  filter(fa_vs_EIF2AK3 == "EIF2AK3") %>%
  summarise(lo = 9, hi = 11) %>%
  as.list()

# Add a little padding to the limits
pad_range <- function(lo, hi, mult = 0.05) {
  rng <- hi - lo
  c(lo - mult * rng, hi + mult * rng)
}
yl_factor1 <- pad_range(ylims_factor1$lo, ylims_factor1$hi)
yl_EIF2AK3   <- pad_range(ylims_EIF2AK3$lo,   ylims_EIF2AK3$hi)


# ====================== Plotting helper ======================================
plot_parental_analyte <- function(parental_level, analyte, ylims) {
  df <- summ_df %>%
    dplyr::filter(Parental == parental_level, fa_vs_EIF2AK3 == analyte) %>%
    dplyr::arrange(Condition, Hours)

  # Line type by Parental
  line_type <- if (parental_level == "Parental") "dashed" else "solid"

  # Color by analyte — shades of blue matching manuscript Figure 2F
  cond_cols <- c(
    "UT"  = "#DAEBFA",  # palest blue  (untreated)
    "Ars" = "#AFD2F3",  # light blue   (arsenite)
    "Tg"  = "#74A8DF",  # medium blue  (thapsigargin)
    "Dim" = "#476CA9"   # dark blue    (dimerizable PERK)
  )

  ggplot(df, aes(x = Hours, y = mean, group = Condition, color = Condition)) +
    geom_errorbar(
      aes(ymin = mean - sem, ymax = mean + sem),
      width     = 1,
      linewidth = 0.6,
      alpha     = 0.7
    ) +
    geom_line(
      linewidth = 1.2,
      linetype  = line_type,
      alpha     = 0.7
    ) +
    geom_point(size = 3.5,
               alpha     = 0.7) +
    scale_color_manual(name = "Condition", values = cond_cols) +
    # scale_x_continuous(
    #   name         = "Time (hours)",
    #   breaks       = seq(0, 24, 4),
    #   minor_breaks = NULL,
    #   expand = expansion(add = c(0.5, 0.5))
    # ) +
    scale_x_continuous(
      name   = "Time (hours)",
      breaks = seq(0, 24, 4),
      expand = c(0, 0)
    ) +
    coord_cartesian(xlim = c(0, 24), clip = "off") +
    labs(title = paste0(parental_level, " – ", analyte)) +
    theme_minimal(base_size = 19) +
    theme(
      axis.text   = element_text(colour = "black"),
      axis.title  = element_text(colour = "black"),
      plot.margin = margin(t = 5.5, r = 25, b = 5.5, l = 25)
    )
}



# ====================== Build the four plots =================================
p_parental_factor1 <- plot_parental_analyte("Parental", "Factor1", yl_factor1)
p_dmrperk_factor1  <- plot_parental_analyte("DmrPERK",  "Factor1", yl_factor1)
p_parental_EIF2AK3   <- plot_parental_analyte("Parental", "EIF2AK3",   yl_EIF2AK3)
p_dmrperk_EIF2AK3    <- plot_parental_analyte("DmrPERK",  "EIF2AK3",   yl_EIF2AK3)

# Print to the device (you'll see 4 plots)
print(p_parental_factor1)
print(p_dmrperk_factor1)
print(p_parental_EIF2AK3)
print(p_dmrperk_EIF2AK3)

# ====================== Optional: save all four ===============================
out_dir <- file.path(local_out_dir, "Parental_vs_DmrPERK_byAnalyte")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

base_name <- if (exists("chosen_gene_to_compare")) chosen_gene_to_compare else "EIF2AK3_vs_Factor1"

ggsave(file.path(out_dir, paste0(base_name, "_Parental_Factor1.png")), p_parental_factor1, width = 7, height = 5, dpi = 300)
ggsave(file.path(out_dir, paste0(base_name, "_DmrPERK_Factor1.png")),  p_dmrperk_factor1,  width = 7, height = 5, dpi = 300)
ggsave(file.path(out_dir, paste0(base_name, "_Parental_EIF2AK3.png")),   p_parental_EIF2AK3,   width = 7, height = 5, dpi = 300)
ggsave(file.path(out_dir, paste0(base_name, "_DmrPERK_EIF2AK3.png")),    p_dmrperk_EIF2AK3,    width = 7, height = 5, dpi = 300)

# Save to shared output directory
ggsave(file.path(shared_out_dir, paste0("Figure_2F.png")), p_dmrperk_factor1, width = 7, height = 5, dpi = 300)

# ================================================================================================
# Export source data for Figure 2F
# ================================================================================================

# Use the summary data that was plotted (summ_df contains all four plots data)
# Include the columns relevant to the plots
source_data_2F <- summ_df %>%
  filter(fa_vs_EIF2AK3 != "EIF2AK3", Parental == "DmrPERK") %>%
  select(Hours, mean, sem, Condition, Parental, fa_vs_EIF2AK3) %>%
  rename(Measure = fa_vs_EIF2AK3)

# Save to CSV
write.csv(source_data_2F,
          file = file.path(folder_path, "Figure_2F_source_data.csv"),
          row.names = FALSE)

message("Figure 2F source data saved to: ", file.path(folder_path, "Figure_2F_source_data.csv"))
