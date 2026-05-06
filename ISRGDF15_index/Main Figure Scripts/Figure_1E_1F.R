# ============================================================================
# Figure 1E & 1F: Factor Analysis Plots + Source Data CSVs
# ============================================================================
# Standalone script extracted from Figure_1E_1F_Factor_Analysis.Rmd
#
# Figure 1E: (a) Heatmap of factor loadings (genes × factors)
#            (b) Factor 1 loadings bar chart (sorted)
# Figure 1F: ISR Factor1 Score vs scaled GDF15, colored by condition
#
# This script sources the same helper scripts as the Rmd, runs the factor
# analysis, generates all plots, exports source-data CSVs, and copies
# key PNGs to the shared Results/Figures folder.
# ============================================================================

rm(list = ls())

library(here)
library(tidyverse)
library(ggplot2)
library(psych)
library(pheatmap)
library(RColorBrewer)
library(plotly)
library(tools)

Plot_Save <- "ON"

# ============================================================================
# 1) SOURCE HELPER SCRIPTS  (loads manifest, exprs, data_sub, gene lists)
# ============================================================================
message("\n=== Sourcing helper scripts ===")

source(here("Main Figure Scripts", "Helper Scripts", "Fibroblast_lifespan",
            "Intro1_FB_Lifespan_exprs_manifest_ISR_list.R"))
source(here("Main Figure Scripts", "Helper Scripts", "Fibroblast_lifespan",
            "GO_vs_AnyGenes_vs_JacksonLabs_ListsofGenes.R"))

# ============================================================================
# 2) OUTPUT DIRECTORIES
# ============================================================================
fig_1E_dir <- here("Results", "Figures", "Figure_1E")
fig_1F_dir <- here("Results", "Figures", "Figure_1F")
shared_out_dir <- here("Results", "Figures")
if (!dir.exists(fig_1E_dir)) dir.create(fig_1E_dir, recursive = TRUE)
if (!dir.exists(fig_1F_dir)) dir.create(fig_1F_dir, recursive = TRUE)
if (!dir.exists(shared_out_dir)) dir.create(shared_out_dir, recursive = TRUE)

# Processed data folder (for saving intermediate CSVs)
processed_path <- here("Data", "Fibroblast_lifespan", "Processed")
if (!dir.exists(processed_path)) dir.create(processed_path, recursive = TRUE)

# ============================================================================
# 3) PREPARE DATASETS
# ============================================================================
message("\n=== Preparing datasets ===")

genes_of_interest <- genes_of_interest_plus_gdf15

data_sub <- data_sub %>%
  mutate(intx = paste(Group, Experiment, sep = "_"))

datasets <- process_datasets(data_sub, genes_of_interest)
WT_SURF1_all_Txs <- datasets$WT_SURF1_all_Txs

message("  Samples: ", nrow(WT_SURF1_all_Txs))

# ============================================================================
# 4) FACTOR ANALYSIS (scaled data, varimax rotation)
# ============================================================================
message("\n=== Running Factor Analysis ===")

valid_genes <- genes_of_interest[genes_of_interest %in% colnames(WT_SURF1_all_Txs)]
subset_df   <- WT_SURF1_all_Txs[, valid_genes]
data_scaled <- scale(subset_df)

# Determine number of factors via parallel analysis
message("  Running parallel analysis...")
get_NFacs <- fa.parallel(data_scaled, fa = "fa", plot = FALSE)
Nfacs     <- get_NFacs$nfact
message("  Number of factors: ", Nfacs)

# Run factor analysis
fit <- factanal(data_scaled, Nfacs, rotation = "varimax", scores = "regression")

num_loadings   <- nrow(fit$loadings)
FactorLoadings <- round(fit$loadings[1:num_loadings, ], 3)

message("  Factor loadings: ", nrow(FactorLoadings), " genes × ", ncol(FactorLoadings), " factors")

# ============================================================================
# 4b) BOOTSTRAP CSVs — write the scaling parameters and factor outputs that
#     downstream scripts (Figures 2C/2E/2F/2G/3B/3C/4A) need.
# ============================================================================
# These five CSVs were historically produced by the legacy
# Figure_1E_1F_Factor_Analysis.Rmd under a date-prefixed name
# (e.g. 20250526_Spearman_rho_Gene_vs_DaysGrown_*).
#
# We now write them here under the stable prefix "Fibroblast_FactorAnalysis_"
# so that running the master script from an empty Processed/ folder bootstraps
# every input the downstream figure scripts need — no manual seeding required.
#
# Contents:
#   meanx_for_scaling.csv         — column means of the fibroblast gene matrix
#   sdx_for_scaling.csv           — column SDs  of the fibroblast gene matrix
#   correlations_for_scaling.csv  — gene-gene correlation matrix from factanal
#   FacLoads_12_factors.csv       — factor loadings (genes × factors)
#   factor_scores.csv             — per-sample factor scores
# ============================================================================
message("\n=== Writing bootstrap CSVs (scaling params + factor outputs) ===")

fa_prefix <- "Fibroblast_FactorAnalysis"

mean_x <- attr(data_scaled, "scaled:center")
sd_x   <- attr(data_scaled, "scaled:scale")

write.csv(mean_x,
          file = file.path(processed_path, paste0(fa_prefix, "_meanx_for_scaling.csv")),
          row.names = TRUE)
message("  Saved: ", paste0(fa_prefix, "_meanx_for_scaling.csv"))

write.csv(sd_x,
          file = file.path(processed_path, paste0(fa_prefix, "_sdx_for_scaling.csv")),
          row.names = TRUE)
message("  Saved: ", paste0(fa_prefix, "_sdx_for_scaling.csv"))

write.csv(fit$correlation,
          file = file.path(processed_path, paste0(fa_prefix, "_correlations_for_scaling.csv")),
          row.names = TRUE)
message("  Saved: ", paste0(fa_prefix, "_correlations_for_scaling.csv"))

write.csv(FactorLoadings,
          file = file.path(processed_path, paste0(fa_prefix, "_FacLoads_12_factors.csv")),
          row.names = TRUE)
message("  Saved: ", paste0(fa_prefix, "_FacLoads_12_factors.csv"))

write.csv(fit$scores,
          file = file.path(processed_path, paste0(fa_prefix, "_factor_scores.csv")),
          row.names = TRUE)
message("  Saved: ", paste0(fa_prefix, "_factor_scores.csv"))

# ============================================================================
# 5) FIGURE 1E — HEATMAP
# ============================================================================
message("\n=== Figure 1E: Heatmap ===")

# IMPORTANT: namespace-qualify the call as `pheatmap::pheatmap` so we always
# get the real pheatmap, even if a helper script (or knit-from-Rmd session)
# loaded ComplexHeatmap. ComplexHeatmap exports its own pheatmap() that masks
# the original and translates to ComplexHeatmap::Heatmap(), which silently
# drops `filename`/`width`/`height` and renders with different defaults.
heatmap_file <- file.path(shared_out_dir, "Figure_1E_Heatmap.png")
pheatmap::pheatmap(FactorLoadings,
         cluster_rows  = TRUE,
         cluster_cols  = FALSE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         filename = heatmap_file,
         width = 10, height = 15)
message("  Saved: ", heatmap_file)

# ============================================================================
# 6) FIGURE 1E — FACTOR 1 LOADINGS BAR CHART
# ============================================================================
message("\n=== Figure 1E: Factor 1 Loadings Bar Chart ===")

pheatmap_default_colors <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)

factor1_col    <- FactorLoadings[, "Factor1", drop = FALSE]
factor1_bar_df <- data.frame(
  variable = rownames(FactorLoadings),
  loading  = factor1_col[, 1]
)
factor1_bar_df <- factor1_bar_df[order(factor1_bar_df$loading), ]
factor1_bar_df$variable <- factor(factor1_bar_df$variable, levels = factor1_bar_df$variable)

p_factor1_bar <- ggplot(factor1_bar_df, aes(x = variable, y = loading, fill = loading)) +
  geom_col() +
  coord_flip() +
  labs(title = "Factor 1 Loadings",
       x = "Variable",
       y = "Loading") +
  scale_fill_gradientn(colors = pheatmap_default_colors) +
  theme_minimal()

# print(p_factor1_bar)

bar_file <- file.path(fig_1E_dir, "Figure_1E_Factor1_Loadings.png")
ggsave(bar_file, plot = p_factor1_bar, width = 14, height = 14, dpi = 300)
message("  Saved: ", bar_file)

# # Also copy to shared output
file.copy(bar_file, file.path(shared_out_dir, "Figure_1E_Factor1_Loadings.png"), overwrite = TRUE)

# ---- All factor bar charts (Factor 1 through Nfacs) ----
for (i in 1:Nfacs) {
  factor_column_name <- paste0("Factor", i)
  chosen_columns <- FactorLoadings[, factor_column_name, drop = FALSE]
  ordered_df <- data.frame(
    variable = rownames(FactorLoadings),
    loading  = chosen_columns[, 1]
  )
  ordered_df <- ordered_df[order(ordered_df$loading), ]
  ordered_df$variable <- factor(ordered_df$variable, levels = ordered_df$variable)

  p <- ggplot(ordered_df, aes(x = variable, y = loading, fill = loading)) +
    geom_col() +
    coord_flip() +
    labs(title = paste("Factor", i, "Loadings"), x = "Variable", y = "Loading") +
    scale_fill_gradientn(colors = pheatmap_default_colors) +
    theme_minimal()

  ggsave(file.path(fig_1E_dir, paste0("Factor_", i, ".png")),
         plot = p, width = 14, height = 14, dpi = 300)
}
message("  Saved bar charts for all ", Nfacs, " factors")

# ============================================================================
# 7) FIGURE 1F — FACTOR1 vs GDF15 SCATTER (colored by condition)
# ============================================================================
message("\n=== Figure 1F: Factor1 vs GDF15 ===")

# Re-run factanal to get scores aligned with WT_SURF1_all_Txs
fa_result <- factanal(data_scaled, Nfacs, rotation = "varimax", scores = "regression")

# Build the Figure 1F dataset
WT_SURF1_all_Txs_GDF15 <- data_sub %>%
  select(SampleID, Experiment, Group, DaysGrown, intx, any_of("GDF15")) %>%
  filter(Group %in% c("Control", "SURF1"))

WT_SURF1_all_Txs_GDF15$GDF15_scaled <- scale(WT_SURF1_all_Txs_GDF15$GDF15)[, 1]
WT_SURF1_all_Txs_GDF15$fa_Scores    <- fa_result$scores[, "Factor1"]

# Spearman correlation
spearman_test <- cor.test(WT_SURF1_all_Txs_GDF15$GDF15_scaled,
                          WT_SURF1_all_Txs_GDF15$fa_Scores,
                          method = "spearman")
spearman_rho <- spearman_test$estimate
spearman_p   <- spearman_test$p.value

if (spearman_p < 0.0001) {
  p_value_text <- "p < 0.0001"
} else {
  p_value_text <- paste("p-value:", signif(spearman_p, 3))
}

message("  Spearman rho = ", round(spearman_rho, 3), ", ", p_value_text)

color_palette <- c(
  "Control_No_Tx"                  = "gray",
  "Control_DEX"                    = "lightcoral",
  "SURF1_No_Tx"                    = "purple",
  "SURF1_DEX"                      = "maroon",
  "Control_Oligo"                  = "violet",
  "Control_mitoNUITs"              = "orange",
  "Control_DEX_mitoNUITs"          = "lightpink",
  "Control_DEX_Oligo"              = "blue",
  "Control_ox3"                    = "lightgreen",
  "Control_Contact_Inhibition"     = "green",
  "Control_Contact_Inhibition_ox3" = "darkgreen",
  "Control_Galactose"              = "lightblue",
  "Control_betahydroxybutyrate"    = "skyblue",
  "Control_2DG"                    = "yellow"
)

p_1F <- ggplot(WT_SURF1_all_Txs_GDF15,
               aes(x = GDF15_scaled, y = fa_Scores, color = intx)) +
  geom_point(alpha = 0.5, size = 5) +
  geom_smooth(method = "lm", color = "red") +
  scale_color_manual(values = color_palette) +
  labs(x = "Standardized GDF15 Expression",
       y = "Factor1 Scores",
       title = paste("GDF15 and Factor1 Score",
                     " Spearman's rho:", round(spearman_rho, 3),
                     p_value_text)) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 14),
    axis.text  = element_text(size = 14),
    plot.title = element_text(size = 16, face = "bold")
  )


fig1F_file <- file.path(fig_1F_dir, "Figure_1F.png")
ggsave(fig1F_file, plot = p_1F, width = 14, height = 10, dpi = 300)
message("  Saved: ", fig1F_file)

# Copy to shared output
file.copy(fig1F_file, file.path(shared_out_dir, "Figure_1F.png"), overwrite = TRUE)

# ============================================================================
# 8) EXPORT SOURCE DATA CSVs
# ============================================================================
message("\n=== Exporting source data CSVs ===")

# Figure 1E Tab 1: Full heatmap data
heatmap_export <- as.data.frame(FactorLoadings)
heatmap_export$Gene <- rownames(heatmap_export)
heatmap_export <- heatmap_export[, c("Gene", setdiff(colnames(heatmap_export), "Gene"))]
write.csv(heatmap_export,
          file.path(fig_1E_dir, "Figure_1E_heatmap_source_data.csv"),
          row.names = FALSE)
message("  Saved: Figure_1E_heatmap_source_data.csv")

# Figure 1E Tab 2: Factor 1 bar chart data (sorted)
factor1_export <- data.frame(
  Gene            = rownames(factor1_col),
  Factor1_Loading = factor1_col[, 1]
)
factor1_export <- factor1_export[order(factor1_export$Factor1_Loading), ]
write.csv(factor1_export,
          file.path(fig_1E_dir, "Figure_1E_Factor1_bar_source_data.csv"),
          row.names = FALSE)
message("  Saved: Figure_1E_Factor1_bar_source_data.csv")

# Figure 1F: ISR Factor1 Score vs GDF15
fig1F_data <- WT_SURF1_all_Txs %>%
  mutate(ISR_Factor1_Score = fa_result$scores[, "Factor1"],
         GDF15_scaled = scale(GDF15)[, 1]) %>%
  select(SampleID, intx, ISR_Factor1_Score, GDF15_scaled)
write.csv(fig1F_data,
          file.path(fig_1F_dir, "Figure_1F_source_data.csv"),
          row.names = FALSE)
message("  Saved: Figure_1F_source_data.csv")

# ============================================================================
# 9) PER-FACTOR vs DAYSGROWN — write CSVs needed by Figure 3D
# ============================================================================
# Ported from Figure_1E_1F_Factor_Analysis.Rmd (Correlation FA score vs AGE
# chunk, ~lines 1006-1170). For each factor, compute the per-sample factor
# score, attach DaysGrown / Group / intx metadata, build a scatter plot of
# factor score vs DaysGrown (colored by intx) with Spearman rho in the title,
# and write a per-sample CSV that Figure 3D loads to compute per-condition
# Spearman rhos.
#
# Output directory MUST match fb_data_dir in Figure_3D.R:
#   Results/Fibroblast_lifespan/Factor_Analysis/Generate_All_Figures/
#       Spearman_Rhos_Age_vs_Genes/
# ============================================================================
message("\n=== Writing per-Factor DaysGrown CSVs (for Figure 3D) ===")

age_csv_dir <- here("Results", "Fibroblast_lifespan", "Factor_Analysis",
                    "Generate_All_Figures", "Spearman_Rhos_Age_vs_Genes")
age_plot_dir <- file.path(age_csv_dir, "FA_All_Factors_vs_AgeDaysGrown")
if (!dir.exists(age_csv_dir))  dir.create(age_csv_dir,  recursive = TRUE)
if (!dir.exists(age_plot_dir)) dir.create(age_plot_dir, recursive = TRUE)

# Same color palette the Rmd uses for these scatter plots
age_color_palette <- c(
  "Control_No_Tx"                  = "gray",
  "Control_DEX"                    = "lightcoral",
  "SURF1_No_Tx"                    = "purple",
  "SURF1_DEX"                      = "maroon",
  "Control_Oligo"                  = "violet",
  "Control_mitoNUITs"              = "orange",
  "Control_DEX_mitoNUITs"          = "lightpink",
  "Control_DEX_Oligo"              = "blue",
  "Control_ox3"                    = "lightgreen",
  "Control_Contact_Inhibition"     = "green",
  "Control_Contact_Inhibition_ox3" = "darkgreen",
  "Control_Galactose"              = "lightblue",
  "Control_betahydroxybutyrate"    = "skyblue",
  "Control_2DG"                    = "yellow"
)

# fa_result$scores rows are aligned with WT_SURF1_all_Txs rows (same order
# from factor analysis). Loop over every Factor* column.
factor_cols <- grep("^Factor", colnames(fa_result$scores), value = TRUE)
message("  Writing CSVs for ", length(factor_cols), " factors")

for (chosen_Factor in factor_cols) {
  AgeCorr <- WT_SURF1_all_Txs
  AgeCorr$fa_Scores <- fa_result$scores[, chosen_Factor]

  # Per-sample CSV in the exact shape Figure_3D.R expects:
  #   columns: SampleID, Group, DaysGrown, intx, fa_Scores
  FA_Index_Scores <- AgeCorr %>%
    select(SampleID, Group, DaysGrown, intx, fa_Scores)

  csv_path <- file.path(age_csv_dir,
                        paste0(chosen_Factor,
                               "_ISR_Scores_per_Sample_Age_DaysGrown.csv"))
  write.csv(FA_Index_Scores, csv_path)

  # Spearman rho across all samples (for the plot title only)
  AgeCorr_clean <- AgeCorr %>% filter(!is.na(DaysGrown), !is.na(fa_Scores))
  if (nrow(AgeCorr_clean) >= 3) {
    spearman_test <- suppressWarnings(
      cor.test(AgeCorr_clean$DaysGrown, AgeCorr_clean$fa_Scores,
               method = "spearman")
    )
    spearman_rho <- spearman_test$estimate
    spearman_p   <- spearman_test$p.value
    p_value_text <- if (!is.na(spearman_p) && spearman_p < 0.0001) {
      "p < 0.0001"
    } else {
      paste("p-value:", signif(spearman_p, 3))
    }

    p_age <- ggplot(AgeCorr_clean,
                    aes(x = DaysGrown, y = fa_Scores, color = intx)) +
      geom_point(alpha = 0.5, size = 5) +
      geom_smooth(method = "lm", color = "red") +
      scale_color_manual(values = age_color_palette) +
      labs(x = "DaysGrown",
           y = paste(chosen_Factor, " Scores"),
           title = paste("Age and", chosen_Factor, " Score",
                         " Spearman's rho:", round(spearman_rho, 3),
                         p_value_text)) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(size = 24),
        axis.text.y = element_text(size = 24),
        axis.title  = element_text(size = 14),
        plot.title  = element_text(size = 16, face = "bold")
      )

    plot_path <- file.path(age_plot_dir,
                           paste0(chosen_Factor, "scoresvsAge_biggerdots.png"))
    ggsave(plot_path, plot = p_age, width = 14, height = 10, dpi = 300)
  }
}
message("  Saved per-Factor DaysGrown CSVs to: ", age_csv_dir)

# ============================================================================
# 10) SPEARMAN RHO SUMMARY (POOLED, w/ BOOTSTRAP 95% CI) — for Figure 3D
# ============================================================================
# Ported from Figure_1E_1F_Factor_Analysis.Rmd (Chosen Genes vs age chunk,
# ~lines 1273-1412). For each factor and each gene metric, compute the
# Spearman rho vs DaysGrown across all samples (Control + SURF1, all txs),
# plus a 1000-rep percentile bootstrap 95% CI. Figure_3D.R reads this CSV
# to draw the black mean point + horizontal error bar for each row.
#
# Output file MUST match Figure_3D.R (line 465):
#   <Spearman_Rhos_Age_vs_Genes>/Spearman_Rho_Summary_With_CI.csv
# ============================================================================
message("\n=== Computing Spearman summary with bootstrap CI ===")

if (!requireNamespace("boot", quietly = TRUE)) install.packages("boot")
library(boot)

# ---- Build the same `data` frame the Rmd builds ----
other_genes <- c("GDF15", "ATF4", "ATF5", "DDIT3",
                 "CDKN1A", "CDKN2A", "CCND2",
                 "MKI67",  "RRM2",   "TOP2A")

# Re-run process_datasets() with the gene panel (NOT the factor analysis genes)
gene_datasets <- process_datasets(data_sub, other_genes)
gene_data     <- gene_datasets$WT_SURF1_all_Txs

# Scale the gene columns (z-score) — matches the Rmd
present_genes <- intersect(other_genes, colnames(gene_data))
gene_data[, present_genes] <- scale(gene_data[, present_genes])

# Per-row sen_mean (CDKN1A, CDKN2A, CCND2) and prolif_mean (MKI67, RRM2, TOP2A)
sen_genes    <- intersect(c("CDKN1A", "CDKN2A", "CCND2"), present_genes)
prolif_genes <- intersect(c("MKI67",  "RRM2",   "TOP2A"), present_genes)
gene_data <- gene_data %>%
  rowwise() %>%
  mutate(
    sen_mean    = mean(c_across(all_of(sen_genes)),    na.rm = TRUE),
    prolif_mean = mean(c_across(all_of(prolif_genes)), na.rm = TRUE)
  ) %>%
  ungroup()

# Per-sample factor scores (Factor1..FactorNfacs aligned with WT_SURF1_all_Txs)
FA_Scores <- data.frame(
  SampleID = WT_SURF1_all_Txs$SampleID,
  fa_result$scores,
  check.names = FALSE
)

# Merge gene metrics with factor scores by SampleID
summary_data <- merge(gene_data, FA_Scores, by = "SampleID")

# ---- Loop over each column and compute Spearman rho + bootstrap CI ----
for_spearman_order <- c("GDF15", "sen_mean", "prolif_mean",
                        "ATF4", "ATF5", "DDIT3",
                        "Factor1", "Factor12", "Factor6",  "Factor2",
                        "Factor4", "Factor9",  "Factor8",
                        "Factor3", "Factor10", "Factor7",  "Factor5",
                        "Factor11")

# Restrict to columns actually present in summary_data (in case Nfacs < 12)
for_spearman_order <- intersect(for_spearman_order, colnames(summary_data))
n_comparisons <- length(for_spearman_order)
message("  Computing rho + 1000-rep bootstrap CI for ", n_comparisons, " columns")

spearman_boot_fn <- function(data, indices, col) {
  d <- data[indices, ]
  suppressWarnings(cor(d[[col]], d$DaysGrown, method = "spearman"))
}

spearman_results <- data.frame(
  Column           = character(),
  Spearman_Rho     = numeric(),
  P_Value          = numeric(),
  Adjusted_P_Value = numeric(),
  CI_Lower         = numeric(),
  CI_Upper         = numeric(),
  Significant      = character(),
  stringsAsFactors = FALSE
)

# Set seed for reproducible bootstrap CIs
set.seed(42)

for (col in for_spearman_order) {
  test <- suppressWarnings(
    cor.test(summary_data[[col]], summary_data$DaysGrown, method = "spearman")
  )
  rho         <- unname(test$estimate)
  p_value     <- test$p.value
  adj_p_value <- p_value * n_comparisons  # Bonferroni
  significant <- ifelse(!is.na(adj_p_value) & adj_p_value < 0.05, "Yes", "No")

  boot_result <- boot(summary_data, statistic = spearman_boot_fn,
                      R = 1000, col = col)
  ci          <- boot.ci(boot_result, type = "perc")$percent[4:5]

  spearman_results <- rbind(
    spearman_results,
    data.frame(
      Column           = col,
      Spearman_Rho     = rho,
      P_Value          = p_value,
      Adjusted_P_Value = adj_p_value,
      CI_Lower         = ci[1],
      CI_Upper         = ci[2],
      Significant      = significant,
      stringsAsFactors = FALSE
    )
  )
}

summary_csv <- file.path(age_csv_dir, "Spearman_Rho_Summary_With_CI.csv")
write.csv(spearman_results, file = summary_csv, row.names = FALSE)
message("  Saved: ", summary_csv)

# ============================================================================
# DONE
# ============================================================================
message("\n", strrep("=", 60))
message("FIGURES 1E AND 1F COMPLETE")
message(strrep("=", 60))
message("\nFigure 1E outputs:")
message("  ", heatmap_file)
# message("  ", bar_file)
message("  ", file.path(fig_1E_dir, "Figure_1E_heatmap_source_data.csv"))
message("  ", file.path(fig_1E_dir, "Figure_1E_Factor1_bar_source_data.csv"))
message("\nFigure 1F outputs:")
message("  ", fig1F_file)
message("  ", file.path(fig_1F_dir, "Figure_1F_source_data.csv"))
message("\nPer-Factor DaysGrown CSVs (consumed by Figure 3D):")
message("  ", age_csv_dir)
message("\nSpearman summary CSV (consumed by Figure 3D):")
message("  ", summary_csv)
