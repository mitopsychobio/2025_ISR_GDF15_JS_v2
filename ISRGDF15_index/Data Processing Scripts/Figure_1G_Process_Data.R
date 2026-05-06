# ============================================================================
# Figure 1G: ISR Score Boxplots - DATA PROCESSING SCRIPT
# ============================================================================
# This script prepares data for boxplots of ISR (Factor 1) scores by condition.
# Generates CSV output for verification in Prism.
#
# Input: Pre-processed factor scores from Figure 1D OR reads data directly
# Output:
#   - Figure_1G_data.RData (for figure generation)
#   - CSV files for Prism verification
# ============================================================================

# Clear any conflicting variables from previous scripts
rm(list = intersect(ls(), c("data_sub", "factor_scores", "manifest", "boxplot_data")))

library(here)
library(tidyverse)
library(dunn.test)
library(psych)  # for factor analysis

# ============================================================================
# CONFIGURATION
# ============================================================================
output_dir <- here("Data", "Fibroblast_lifespan", "Processed")

# ============================================================================
# ALWAYS READ FRESH DATA (pre-processed files have mismatched IDs)
# ============================================================================
message("Reading data fresh...")

data_loaded <- FALSE

# ============================================================================
# READ DATA USING MEMORY-EFFICIENT APPROACH (same as Intro1 script)
# ============================================================================
if (!data_loaded) {
  message("  Reading raw data directly (memory-efficient approach)...")

  # Read manifest
  manifest <- read.csv(here("Data", "Fibroblast_lifespan", "Lifespan_Study_selected_data.csv"))

  # Clean manifest - same as Intro1
  manifest <- manifest %>%
    filter(!is.na(RNAseq_sampleID)) %>%
    mutate(
      SampleID = paste0("Sample_", stringr::str_trim(as.character(RNAseq_sampleID)))
    ) %>%
    rename(DaysGrown = Days_grown_Udays) %>%
    rename(Cell_Line = Cell_line_inhouse)

  # Create Group column
  manifest$Group <- ifelse(manifest$Cell_Line %in% c("hFB12", "hFB13", "hFB14", "hFB11"), "Control",
                    ifelse(manifest$Cell_Line %in% c("hFB6", "hFB7", "hFB8"), "SURF1", "Ctrl_tech_rep"))

  # Create Experiment column
  manifest <- manifest %>%
    mutate(Experiment = case_when(
      grepl("Mutation_Control", Cell_line_group) ~ "No_Tx",
      grepl("Normal_Control_ox21", Cell_line_group) ~ "No_Tx",
      grepl("Mutation_DEX", Cell_line_group) ~ "DEX",
      grepl("Normal_DEX", Cell_line_group) ~ "DEX",
      grepl("Normal_Oligomycin_", Cell_line_group) ~ "Oligo",
      grepl("Normal_Oligomycin+", Cell_line_group) ~ "DEX_Oligo",
      grepl("Normal_mitoNUITs_", Cell_line_group) ~ "mitoNUITs",
      grepl("Normal_mitoNUITs+", Cell_line_group) ~ "DEX_mitoNUITs",
      grepl("Galactose", Cell_line_group) ~ "Galactose",
      grepl("2DG", Cell_line_group) ~ "2DG",
      grepl("betahydroxybutyrate", Cell_line_group) ~ "betahydroxybutyrate",
      grepl("Control_ox3", Cell_line_group) ~ "ox3",
      grepl("Inhibition_ox21", Cell_line_group) ~ "Contact_Inhibition",
      grepl("Inhibition_ox3", Cell_line_group) ~ "Contact_Inhibition_ox3",
      TRUE ~ NA_character_
    ))

  # Read expression data (wide format - genes as columns)
  message("  Reading expression data...")
  RNAseq_exp_log2 <- read.csv(here("Data", "Fibroblast_lifespan", "GSE179848_processed_cell_lifespan_RNAseq_data.csv"))
  colnames(RNAseq_exp_log2)[1] <- "gene"

  # Transpose: samples as rows, genes as columns (memory efficient)
  message("  Transposing expression matrix...")
  gene_names <- RNAseq_exp_log2$gene
  sample_names <- colnames(RNAseq_exp_log2)[-1]

  exprs_matrix <- t(RNAseq_exp_log2[, -1])
  colnames(exprs_matrix) <- gene_names
  exprs_df <- as.data.frame(exprs_matrix)
  exprs_df$SampleID <- sample_names

  # Merge with manifest
  message("  Merging with manifest...")
  data_sub <- merge(exprs_df, manifest, by = "SampleID", all.x = TRUE)

  # Create interaction column
  data_sub <- data_sub %>%
    mutate(intx = paste(Group, Experiment, sep = "_"))

  # Read ISR gene list
  ISR <- read.csv(here("Data", "Fibroblast_lifespan", "ISR_Gene_Lists_updated.csv"))
  long_ISR <- ISR %>%
    pivot_longer(cols = everything(), names_to = "Source", values_to = "Gene") %>%
    distinct()
  long_ISR$Gene <- toupper(long_ISR$Gene)
  genes_of_interest <- unique(long_ISR$Gene[long_ISR$Gene != ""])
  genes_of_interest <- c(genes_of_interest, "GDF15")

  # Filter to available genes
  valid_genes <- genes_of_interest[genes_of_interest %in% colnames(data_sub)]
  message("  Valid ISR genes found: ", length(valid_genes))

  # Get gene data for factor analysis
  gene_data <- data_sub[, valid_genes, drop = FALSE]

  # Convert to numeric (in case they're character after transpose)
  gene_data <- as.data.frame(lapply(gene_data, as.numeric))

  complete_rows <- complete.cases(gene_data)
  gene_data <- gene_data[complete_rows, ]
  data_sub <- data_sub[complete_rows, ]

  message("  Samples with complete data: ", nrow(data_sub))

  # Run factor analysis
  message("  Running factor analysis...")
  data_scaled <- scale(gene_data)
  Nfacs <- min(12, ncol(data_scaled) - 1)
  fit <- factanal(data_scaled, Nfacs, rotation = "varimax", scores = "regression")
  factor_scores <- as.data.frame(fit$scores)
  rownames(factor_scores) <- data_sub$SampleID

  # Clean up large objects
  rm(exprs_matrix, exprs_df, gene_data, data_scaled)
  gc()

  data_loaded <- TRUE
}

if (!data_loaded) {
  stop("ERROR: Could not load data. Check that data files exist.")
}

# ============================================================================
# PROCESS DATA FOR BOXPLOTS
# ============================================================================
message("\nProcessing data for boxplots...")

# Debug: check what columns exist
message("  Columns in data_sub: ", paste(head(names(data_sub), 10), collapse = ", "), "...")

# Check if intx column exists, create if not
if (!"intx" %in% names(data_sub)) {
  message("  Creating intx column...")
  data_sub <- data_sub %>%
    mutate(intx = paste(Group, Experiment, sep = "_"))
}

boxplot_data <- data.frame(

  SampleID = data_sub$SampleID,
  intx = data_sub$intx,
  Group = data_sub$Group,
  Experiment = data_sub$Experiment
)

message("  Initial boxplot_data rows: ", nrow(boxplot_data))
message("  Unique intx values: ", length(unique(boxplot_data$intx)))

# Match factor scores
matched_idx <- match(boxplot_data$SampleID, rownames(factor_scores))
message("  Matched samples: ", sum(!is.na(matched_idx)), " of ", length(matched_idx))

boxplot_data$ISR_Score <- factor_scores$Factor1[matched_idx]

# Remove NA and filter out conditions with "NA" in name
boxplot_data <- boxplot_data %>%
  filter(!is.na(ISR_Score)) %>%
  filter(!grepl("_NA$", intx)) %>%
  filter(!grepl("^NA_", intx))

message("  After filtering - Samples: ", nrow(boxplot_data))
message("  After filtering - Conditions: ", length(unique(boxplot_data$intx)))

if (nrow(boxplot_data) == 0) {
  stop("ERROR: No data remaining after filtering. Check SampleID matching.")
}

if (length(unique(boxplot_data$intx)) < 2) {
  message("  WARNING: Only one condition found. Listing unique intx values:")
  print(unique(boxplot_data$intx))
  stop("ERROR: Need at least 2 groups for statistical tests.")
}

# ============================================================================
# CALCULATE CONDITION STATISTICS (sorted by ascending median)
# ============================================================================
message("\nCalculating condition statistics...")

condition_stats <- boxplot_data %>%
  group_by(intx) %>%
  summarise(
    n = n(),
    median = median(ISR_Score),
    mean = mean(ISR_Score),
    sd = sd(ISR_Score),
    min = min(ISR_Score),
    max = max(ISR_Score),
    Q1 = quantile(ISR_Score, 0.25),
    Q3 = quantile(ISR_Score, 0.75),
    .groups = "drop"
  ) %>%
  arrange(median)

message("  Conditions by median ISR score (ascending):")
for (i in 1:nrow(condition_stats)) {
  message("    ", i, ". ", condition_stats$intx[i], ": median = ", round(condition_stats$median[i], 3))
}

boxplot_data$intx <- factor(boxplot_data$intx, levels = condition_stats$intx)

# ============================================================================
# STATISTICAL TESTS
# ============================================================================
message("\nRunning Kruskal-Wallis test...")

kw_test <- kruskal.test(ISR_Score ~ intx, data = boxplot_data)
message("  Kruskal-Wallis chi-squared: ", round(kw_test$statistic, 2))
message("  P-value: ", format(kw_test$p.value, scientific = TRUE, digits = 3))

message("\nRunning Dunn's post-hoc test...")
dunn_result <- dunn.test(boxplot_data$ISR_Score, boxplot_data$intx, method = "bonferroni", kw = FALSE)

pairwise_comparisons <- data.frame(
  comparison = dunn_result$comparisons,
  Z = dunn_result$Z,
  P.adjusted = dunn_result$P.adjusted
)

significant_pairs <- pairwise_comparisons %>% filter(P.adjusted < 0.05)
message("  Significant pairwise comparisons: ", nrow(significant_pairs))

# # ============================================================================
# # SAVE CSV FILES FOR PRISM
# # ============================================================================
# message("\nSaving CSV files for Prism...")
# 
# # Wide format
# prism_wide <- boxplot_data %>%
#   select(intx, ISR_Score) %>%
#   group_by(intx) %>%
#   mutate(row_id = row_number()) %>%
#   pivot_wider(names_from = intx, values_from = ISR_Score) %>%
#   select(-row_id)
# prism_wide <- prism_wide[, as.character(condition_stats$intx)]
# 
# write.csv(prism_wide, file.path(csv_output_dir, "Figure_1G_ISR_Scores_by_Condition_WIDE.csv"), row.names = FALSE, na = "")
# message("  Saved: Figure_1G_ISR_Scores_by_Condition_WIDE.csv")
# 
# # Long format
# prism_long <- boxplot_data %>%
#   select(SampleID, intx, Group, Experiment, ISR_Score) %>%
#   arrange(factor(intx, levels = condition_stats$intx), ISR_Score)
# 
# write.csv(prism_long, file.path(csv_output_dir, "Figure_1G_ISR_Scores_by_Condition_LONG.csv"), row.names = FALSE)
# message("  Saved: Figure_1G_ISR_Scores_by_Condition_LONG.csv")
# 
# # Statistics
# write.csv(condition_stats, file.path(csv_output_dir, "Figure_1G_Condition_Statistics.csv"), row.names = FALSE)
# message("  Saved: Figure_1G_Condition_Statistics.csv")
# 
# write.csv(pairwise_comparisons, file.path(csv_output_dir, "Figure_1G_Pairwise_Comparisons.csv"), row.names = FALSE)
# message("  Saved: Figure_1G_Pairwise_Comparisons.csv")

# ============================================================================
# SAVE PROCESSED DATA FOR FIGURE GENERATION
# ============================================================================
message("\nSaving processed data...")

figure_1g_data <- list(
  boxplot_data = boxplot_data,
  condition_stats = condition_stats,
  kw_test = kw_test,
  pairwise_comparisons = pairwise_comparisons,
  significant_pairs = significant_pairs
)

save(figure_1g_data, file = file.path(output_dir, "Figure_1G_data.RData"))
message("  Saved: Figure_1G_data.RData")

message("\n============ DATA PROCESSING COMPLETE ============")
