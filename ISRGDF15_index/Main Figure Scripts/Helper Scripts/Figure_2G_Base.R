

# Control F and then replace "Factor1" with Factor_
#Control F and replace the chosen gene name but make sure that the gene list csv stays with GDF15 or it wont run
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

chosen_gene_to_compare <- "GDF15"

# Output directories
local_out_dir <- here("Results", "Figures", "Figure_2G")
shared_out_dir <- here("Results", "Figures")
if (!dir.exists(local_out_dir)) dir.create(local_out_dir, recursive = TRUE)
if (!dir.exists(shared_out_dir)) dir.create(shared_out_dir, recursive = TRUE)
folder_path <- local_out_dir  # Keep folder_path for compatibility with rest of script

# Readin <- here("Data", "NatComms_K_Labbe_et_al", "Taivan_Normalized_counts_text.xlsx")
# TPM_data <- read_excel(Readin)
tpm_path <- here("Data", "NatComms_K_Labbe_et_al", "GSE273599_WT_ATF4KO_TPMs.txt.gz")
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


exprs_GDF15 <- exprs_ISR %>%
  select(Timepoint, GDF15)

# exprs_GADD34 <- exprs_ISR_scaled %>%
#   select(Timepoint, GDF15)

########################################################################################

scaled_exprs_ISR <- exprs_ISR_scaled

loadings_readin <- here("Data", "Fibroblast_lifespan", "Processed", "Fibroblast_FactorAnalysis_FacLoads_12_factors.csv") # produced by Figure_1E_1F.R; includes GDF15
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
    left_join(exprs_GDF15 %>% select(Timepoint, GDF15), by = "Timepoint") %>%
    pivot_longer(cols = c(all_of(factor_name), "GDF15"), names_to = "fa_vs_GDF15", values_to = "value") %>%
    mutate(
      type = ifelse(fa_vs_GDF15 == factor_name, "fa_Scores", "GDF15_scaled"),
      Time_numeric = as.numeric(gsub("X(\\d+)_.*", "\\1", Timepoint))
    )
}


merged_df <- merged_df %>%
  tidyr::extract(
    Timepoint,
    into  = c("Condition", "Kinase", "Agent", "Hours", "Rep"),
    regex = "^([^_]+)_([^_]+)_([^_]+)_(\\d+)h_(rep\\d+)$",  # ← fixed
    remove = FALSE
  ) %>%
  mutate(Hours = as.numeric(Hours)) %>%
  select(-any_of("Time_numeric"))   # drop if present


# ── 2) Summary stats per Kinase x Agent x Hours x fa_vs_GDF15 ─────────────────
# sem = sd/sqrt(n)
# cv% = (sd/mean)*100; protect against mean==0
summ_df <- merged_df %>%
  group_by(Condition, Kinase, Hours, fa_vs_GDF15) %>%
  summarise(
    n     = dplyr::n(),
    mean  = mean(value, na.rm = TRUE),
    sd    = sd(value,   na.rm = TRUE),
    sem   = sd / sqrt(n),
    cv    = dplyr::if_else(is.finite(mean) & mean != 0, (sd / abs(mean)) * 100, NA_real_),
    .groups = "drop"
  )



# ------------------------------------------------------------------------------------------------
# Plots: Factor1 and GDF15
# ------------------------------------------------------------------------------------------------

# ---- Aesthetics (keep yours) ----
line_size   <- 1.2
point_size  <- 4
base_size   <- 12
error_bar_width <- 1   # your example uses 1; set to 2 if you really want wider caps

plot_agent_analyte_split <- function(df_agent,
                                     analyte = c("Factor1", "GDF15"),
                                     agent_label = unique(df_agent$Kinase),
                                     ylims = NULL) {
  analyte <- match.arg(analyte)

  df <- df_agent %>%
    dplyr::filter(fa_vs_GDF15 == analyte) %>%
    dplyr::arrange(Condition, Hours)

  # Color by analyte (Factor1 black, GDF15 pink)
  color_val <- if (analyte == "Factor1") "#0b3d91" else "pink"

  # Linetype map: WT solid, everything else dashed (generalized)
  cond_levels <- sort(unique(df$Condition))
  lt_values   <- setNames(rep("dashed", length(cond_levels)), cond_levels)
  if ("WT" %in% names(lt_values)) lt_values["WT"] <- "solid"

  # Shape map (auto, but sensible defaults if present)
  shape_pool <- c(16, 17, 15, 3, 7, 8, 0, 1, 2, 5, 6)
  shape_vals <- setNames(shape_pool[seq_along(cond_levels)], cond_levels)
  if ("WT" %in% names(shape_vals)) shape_vals["WT"] <- 16
  if ("ATF4KO" %in% names(shape_vals)) shape_vals["ATF4KO"] <- 16

  p <- ggplot(df, aes(x = Hours, y = mean, group = Condition)) +
    # SEM error bars (match your example)
    geom_errorbar(
      aes(ymin = mean - sem, ymax = mean + sem),
      width     = error_bar_width,
      linewidth = 0.6,
      color     = color_val,
      alpha     = 0.7
    ) +
    # Lines (slightly translucent)
    geom_line(
      aes(linetype = Condition),
      linewidth = line_size,
      color     = color_val,
      alpha     = 0.7
    ) +
    # Points (solid)
    geom_point(
      aes(shape = Condition),
      size  = point_size,
      color = color_val
    ) +
    scale_shape_manual(
      name   = "Condition",
      values = shape_vals
    ) +
    scale_linetype_manual(
      name   = "Condition",
      values = lt_values
    ) +
    scale_x_continuous(
      name   = "Time (hours)",
      breaks = seq(0, 24, 4),
      expand = c(0, 0)
    ) +
    coord_cartesian(xlim = c(0, 24), clip = "off") +
    theme_minimal(base_size = 19) +
    labs(title = paste0("Timecourse: ", agent_label, " – ", analyte)) +
    theme(
      axis.text   = element_text(colour = "black"),
      axis.title  = element_text(colour = "black"),
      plot.margin = margin(t = 5.5, r = 25, b = 5.5, l = 25)
    )

  p
}

# ---- Build all Kinase plots: TWO per Kinase (Factor1 + GDF15) ----
agent_list  <- summ_df %>% dplyr::group_split(Kinase)
agent_names <- summ_df %>% dplyr::distinct(Kinase) %>% dplyr::pull(Kinase)

plots_split_tbl <- purrr::map2_dfr(agent_list, agent_names, ~{
  tibble::tibble(
    Kinase  = .y,
    analyte = c("Factor1", "GDF15"),
    plot    = list(
      plot_agent_analyte_split(.x, analyte = "Factor1", agent_label = .y),
      plot_agent_analyte_split(.x, analyte = "GDF15",   agent_label = .y)
    )
  )
})

# Print all
invisible(purrr::walk(plots_split_tbl$plot, print))

# ---- Optional: save ----
out_dir <- file.path(local_out_dir, "Agent_Plots_split")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

purrr::pwalk(
  list(plots_split_tbl$plot, plots_split_tbl$Kinase, plots_split_tbl$analyte),
  ~ ggsave(
    filename = file.path(out_dir, paste0(chosen_gene_to_compare, "_Timecourse_", ..2, "_", ..3, ".png")),
    plot     = ..1,
    width    = 7,
    height   = 5,
    dpi      = 300
  )
)

# Save main figure to shared output directory
# Get the first Kinase's Factor1 plot as the main figure
first_kinase <- plots_split_tbl$Kinase[1]
main_plot <- plots_split_tbl %>%
  filter(Kinase == first_kinase, analyte == "Factor1") %>%
  pull(plot) %>%
  .[[1]]

ggsave(
  file.path(shared_out_dir, paste0("Figure_2G.png")),
  plot   = main_plot,
  width  = 7,
  height = 5,
  dpi    = 300
)

# ================================================================================================
# Export source data for Figure 2G
# ================================================================================================

# Use the summary data that was plotted (summ_df contains all plot data)
# Include the columns relevant to the plots
source_data_2G <- summ_df %>%
  select(Hours, mean, sem, Condition, Kinase, fa_vs_GDF15) %>%
  rename(Measure = fa_vs_GDF15)

# Main figure CSV: Factor1 only (GDF15 is in extended data)
source_data_2G_main <- source_data_2G %>%
  filter(Measure != "GDF15")

write.csv(source_data_2G_main,
          file = file.path(folder_path, "Figure_2G_source_data.csv"),
          row.names = FALSE)
message("Figure 2G source data saved to: ", file.path(folder_path, "Figure_2G_source_data.csv"))

# Extended Data Figure 11C: GDF15 only (separate CSV, not included in master Excel)
source_data_ext_11C <- source_data_2G %>%
  filter(Measure == "GDF15")

write.csv(source_data_ext_11C,
          file = file.path(folder_path, "Extended_Data_Figure_11C.csv"),
          row.names = FALSE)
message("Extended Data Figure 11C source data saved to: ", file.path(folder_path, "Extended_Data_Figure_11C.csv"))
