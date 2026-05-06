# legend.position = "top",          # match dual-axis THIS IS WHAT WAS CHANGED IN THIS CODE

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
local_out_dir <- here("Results", "Figures", "Figure_2E")
shared_out_dir <- here("Results", "Figures")
if (!dir.exists(local_out_dir)) dir.create(local_out_dir, recursive = TRUE)
if (!dir.exists(shared_out_dir)) dir.create(shared_out_dir, recursive = TRUE)
folder_path <- local_out_dir  # Keep folder_path for compatibility with rest of script

# Readin <- here("Data", "NatComms_K_Labbe_et_al", "Taivan_Normalized_counts_text.xlsx")
# TPM_data <- read_excel(Readin)
tpm_path <- here("Data", "NatComms_K_Labbe_et_al", "GSE273601_DmrPERK_Dose_Response_TPMs.txt.gz")
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
      -c(gene_id, gene_name),          # all columns except the first two
      ~ log2(as.numeric(.x) + 1)       # +1 to handle zeros safely
    )
  )

TPM_data <-TPM_log2

TPM_data <- TPM_data %>%
  select(-gene_id)

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
    col = Timepoint,
    into = c("Kinase", "Agent", "Hours", "Rep"),
    regex = "^([^_]+)_([^_]+)_(\\d+)h_(rep\\d+)$",
    remove = FALSE
  ) %>%
  mutate(Hours = as.numeric(Hours))%>%
  select(-Time_numeric)


# ── 2) Summary stats per Kinase x Agent x Hours x fa_vs_GDF15 ─────────────────
# sem = sd/sqrt(n)
# cv% = (sd/mean)*100; protect against mean==0
summ_df <- merged_df %>%
  group_by(Kinase, Agent, Hours, fa_vs_GDF15) %>%
  summarise(
    n     = dplyr::n(),
    mean  = mean(value, na.rm = TRUE),
    sd    = sd(value,   na.rm = TRUE),
    sem   = sd / sqrt(n),
    cv    = dplyr::if_else(is.finite(mean) & mean != 0, (sd / abs(mean)) * 100, NA_real_),
    .groups = "drop"
  )

# ── 1) Helper: dual-axis plot for a single Agent ──────────────────────────────
library(dplyr)
library(ggplot2)
library(purrr)
library(tidyr)

# If needed, rebuild summary (same as before)
if (!exists("summ_df")) {
  summ_df <- merged_df %>%
    tidyr::extract(Timepoint, c("Kinase","Agent","Hours","Rep"),
                   "^([^_]+)_([^_]+)_(\\d+)h_(rep\\d+)$", remove = FALSE) %>%
    mutate(Hours = as.numeric(Hours)) %>%
    group_by(Kinase, Agent, Hours, fa_vs_GDF15) %>%
    summarise(n = dplyr::n(),
              mean = mean(value, na.rm = TRUE),
              sd   = sd(value,   na.rm = TRUE),
              sem  = sd/sqrt(n),
              .groups = "drop")
}

# ---- Aesthetics to mirror your example ----
col_factor1      <- "black"
col_GDF15        <- "pink"
line_size        <- 1.2
point_size       <- 4
base_size        <- 12
error_bar_width  <- 1   # or 2 if you prefer the wide bars you used before

plot_dual_axis_agent_match <- function(df_agent, agent_label = unique(df_agent$Agent)) {
  df_fa  <- df_agent %>% dplyr::filter(fa_vs_GDF15 == "Factor1")
  df_gdf <- df_agent %>% dplyr::filter(fa_vs_GDF15 == "GDF15")

  # Map GDF15 (right) ↔ Factor1 (left)
  FA_min  <- min(df_fa$mean,  na.rm = TRUE)
  FA_max  <- max(df_fa$mean,  na.rm = TRUE)
  GDF_min <- min(df_gdf$mean, na.rm = TRUE)
  GDF_max <- max(df_gdf$mean, na.rm = TRUE)

  # Protect against divide-by-zero (all means identical)
  if (FA_max == FA_min) {
    FA_max <- FA_min + 1e-6
  }
  if (GDF_max == GDF_min) {
    GDF_max <- GDF_min + 1e-6
  }

  a <- (GDF_max - GDF_min) / (FA_max - FA_min)
  b <- GDF_min - a * FA_min

  # Add SEM-based error bar columns
  df_fa <- df_fa %>%
    mutate(
      ymin = mean - sem,
      ymax = mean + sem
    )

  df_gdf <- df_gdf %>%
    mutate(
      mean_L = (mean - b) / a,
      ymin_L = (mean - sem - b) / a,
      ymax_L = (mean + sem - b) / a
    )

  ggplot() +
    # ---- Factor1 (left axis) ----
  geom_errorbar(
    data = df_fa,
    aes(x = Hours, ymin = ymin, ymax = ymax, color = "Factor1"),
    width = error_bar_width
  ) +
    geom_line(
      data = df_fa,
      aes(x = Hours, y = mean, color = "Factor1"),
      linewidth = line_size
    ) +
    geom_point(
      data = df_fa,
      aes(x = Hours, y = mean, color = "Factor1"),
      size = point_size
    ) +

    # ---- GDF15 mapped to left axis (with SEM error bars) ----
  geom_errorbar(
    data = df_gdf,
    aes(x = Hours, ymin = ymin_L, ymax = ymax_L, color = "GDF15"),
    width = error_bar_width
  ) +
    geom_line(
      data = df_gdf,
      aes(x = Hours, y = mean_L, color = "GDF15"),
      linewidth = line_size
    ) +
    geom_point(
      data = df_gdf,
      aes(x = Hours, y = mean_L, color = "GDF15"),
      size = point_size
    ) +

    # Axes & scales
    scale_x_continuous(
      name = "Time (hours)",
      breaks = sort(unique(df_agent$Hours))
    ) +
    scale_y_continuous(
      name = "Factor1 value",
      sec.axis = sec_axis(~ . * a + b, name = "GDF15 expression (log2)")
    ) +
    scale_color_manual(
      name = "Legend",
      values = c("Factor1" = col_factor1, "GDF15" = col_GDF15)
    ) +

    labs(title = paste0("Timecourse: ", agent_label)) +
    theme_minimal(base_size = base_size) +
    theme(
      plot.title         = element_text(face = "bold"),
      legend.position    = "top",
      axis.title.y.left  = element_text(color = col_factor1),
      axis.text.y.left   = element_text(color = col_factor1),
      axis.title.y.right = element_text(color = col_GDF15),
      axis.text.y.right  = element_text(color = col_GDF15)
    )
}

# Map: dose agent -> which "control" agent provides the 0h baseline
baseline_map <- c(
  "Dim01"  = "UT",
  "Dim195" = "UT",
  "Dim3"   = "UT2"
)

library(dplyr)
library(purrr)

make_df_agent_with_baseline <- function(agent_name, summ_df, baseline_map) {
  df_agent <- summ_df %>% filter(Agent == agent_name)

  base_agent <- baseline_map[agent_name]

  if (!is.na(base_agent)) {
    df_base <- summ_df %>%
      filter(Agent == base_agent, Hours == 0)

    if (nrow(df_base) > 0) {
      # Relabel baseline Agent so it plots as part of the same series
      df_base <- df_base %>% mutate(Agent = agent_name)

      df_agent <- bind_rows(df_base, df_agent) %>%
        arrange(fa_vs_GDF15, Hours)
    }
  }

  df_agent
}


agents_to_plot <- c("Dim01", "Dim195", "Dim3")

plots_dual_matched <- lapply(agents_to_plot, function(agent) {
  df_agent2 <- make_df_agent_with_baseline(agent, summ_df, baseline_map)
  plot_dual_axis_agent_match(df_agent2, agent_label = agent)
})

names(plots_dual_matched) <- agents_to_plot

invisible(lapply(plots_dual_matched, print))

out_dir <- file.path(local_out_dir, "Agent_Plots_dual_axis_matched")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# pwalk(
#   list(plots_dual_matched, agents_to_plot),
#   ~ ggsave(
#     file.path(out_dir, paste0(chosen_gene_to_compare, "_Timecourse_", ..2, "_dual_axis.png")),
#     plot = ..1, width = 7, height = 5, dpi = 300
#   )
# )




# Color maps for Agents
agent_cols_Factor1 <- c(
  "Dim01"  = "grey80",  # light gray
  "Dim195" = "grey50",  # medium gray
  "Dim3"   = "black"    # black
)

agent_cols_GDF15 <- c(
  "Dim01"  = "#FAD0E8", # light pink
  "Dim195" = "#F48FB1", # medium pink
  "Dim3"   = "#AD1457"  # dark pink
)


df_all_agents <- summ_df


# Reuse the same axis breaks as the dual-axis plots
x_breaks <- sort(unique(df_all_agents$Hours))




make_df_all_agents_with_baselines <- function(summ_df, agents_to_plot, baseline_map) {

  # Keep only the agents you want in the 3-curve plots
  df_main <- summ_df %>%
    dplyr::filter(Agent %in% agents_to_plot) %>%
    dplyr::filter(Hours != 0)  # drop each Dim's own 0h so we replace it with mapped baseline

  df_baselines <- purrr::map_dfr(agents_to_plot, function(agent_name) {

    base_agent <- unname(baseline_map[agent_name])
    if (is.na(base_agent) || !nzchar(base_agent)) return(NULL)

    summ_df %>%
      dplyr::filter(Agent == base_agent, Hours == 0) %>%
      dplyr::mutate(
        Agent = agent_name
        # optional: baseline_source = base_agent
      )
  })

  dplyr::bind_rows(df_baselines, df_main) %>%
    dplyr::arrange(Agent, fa_vs_GDF15, Hours)
}


agents_to_plot <- c("Dim01", "Dim195", "Dim3")

df_all_agents <- make_df_all_agents_with_baselines(
  summ_df       = summ_df,
  agents_to_plot = agents_to_plot,
  baseline_map  = baseline_map
)




# Factor1-only plot (3 curves) — shades of gray
df_fa_all <- df_all_agents %>%
  dplyr::filter(fa_vs_GDF15 == "Factor1") %>%
  dplyr::mutate(
    ymin = mean - sem,
    ymax = mean + sem
  )

p_Factor1_three <- ggplot(
  df_fa_all,
  aes(x = Hours,
      y = mean,
      group = Agent,
      color = Agent)
) +
  # solid error bars, same color as the line
  geom_errorbar(
    aes(ymin = ymin, ymax = ymax),
    width = error_bar_width
  ) +
  geom_line(
    linewidth = line_size,
    alpha     = 0.7    # slight translucency
  ) +
  geom_point(
    size  = point_size,
    alpha = 0.7
  ) +
  scale_color_manual(
    name   = "Agent",
    values = agent_cols_Factor1
  ) +
  scale_x_continuous(
    name   = "Time (hours)",
    breaks = seq(0, 24, 4),
    expand = c(0, 0)
  ) +
  ylab("Factor1 value") +
  labs(
    title = paste0("Timecourse: Factor1 (", chosen_gene_to_compare, ")")
  )  +
  coord_cartesian(xlim = c(0, 24), clip = "off") +
  theme_minimal(base_size = 19) +
  theme(
    plot.title      = element_text(face = "bold"),
    axis.text       = element_text(colour = "black"),
    axis.title      = element_text(colour = "black"),
    # legend.position = "top",          # match dual-axis
    # panel.grid.minor = element_blank(),
    plot.margin = margin(t = 5.5, r = 25, b = 5.5, l = 25)
  )




print(p_Factor1_three)

ggsave(
  file.path(out_dir, paste0(chosen_gene_to_compare,
                            "_Timecourse_Factor1_three_agents_graph_consitstent.png")),
  plot   = p_Factor1_three,
  width  = 7,
  height = 5,
  dpi    = 300
)

# Save to shared output directory
ggsave(
  file.path(shared_out_dir, paste0("Figure_2E.png")),
  plot   = p_Factor1_three,
  width  = 7,
  height = 5,
  dpi    = 300
)


# GDF15-only plot (3 curves) — shades of pink
df_gdf_all <- df_all_agents %>%
  dplyr::filter(fa_vs_GDF15 == "GDF15") %>%
  dplyr::mutate(
    ymin = mean - sem,
    ymax = mean + sem
  )

p_GDF15_three <- ggplot(
  df_gdf_all,
  aes(x = Hours,
      y = mean,
      group = Agent,
      color = Agent)
) +
  geom_errorbar(
    aes(ymin = ymin, ymax = ymax),
    width = error_bar_width
  ) +
  geom_line(
    linewidth = line_size,
    alpha     = 0.7
  ) +
  geom_point(
    size  = point_size,
    alpha = 0.7
  ) +
  scale_color_manual(
    name   = "Agent",
    values = agent_cols_GDF15
  ) +
  scale_x_continuous(
    name   = "Time (hours)",
    breaks = seq(0, 24, 4),
    expand = c(0, 0)
  ) +
  ylab("GDF15 expression (log2)") +
  labs(
    title = paste0("Timecourse: GDF15 (", chosen_gene_to_compare, ")")
  ) +
  coord_cartesian(xlim = c(0, 24), clip = "off") +
  theme_minimal(base_size = 19) +
  theme(
    plot.title      = element_text(face = "bold"),
    axis.text       = element_text(colour = "black"),
    axis.title      = element_text(colour = "black"),
    # legend.position = "top",
    # panel.grid.minor = element_blank(),
    plot.margin = margin(t = 5.5, r = 25, b = 5.5, l = 25)
  )


print(p_GDF15_three)

ggsave(
  file.path(out_dir, paste0(chosen_gene_to_compare,
                            "_Timecourse_GDF15_three_agents_graph_consitstent.png")),
  plot   = p_GDF15_three,
  width  = 7,
  height = 5,
  dpi    = 300
)

# ================================================================================================
# Export source data for Figure 2E
# ================================================================================================

# Use the summary data that was plotted (df_all_agents contains all data)
# Include only the columns relevant to the plots
source_data_2E <- df_all_agents %>%
  filter(fa_vs_GDF15 != "GDF15") %>%
  select(Hours, mean, sem, Agent, fa_vs_GDF15) %>%
  rename(Measure = fa_vs_GDF15)

# Save to CSV
write.csv(source_data_2E,
          file = file.path(local_out_dir, "Figure_2E_source_data.csv"),
          row.names = FALSE)

message("Figure 2E source data saved to: ", file.path(local_out_dir, "Figure_2E_source_data.csv"))
