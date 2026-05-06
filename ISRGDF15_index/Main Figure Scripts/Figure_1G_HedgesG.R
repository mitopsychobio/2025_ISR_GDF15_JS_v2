# ============================================================================
# Figure 1G: Hedges' g Effect Sizes — Each Condition vs Control_No_Tx
# ============================================================================
# Computes Hedges' g (bias-corrected Cohen's d) for each experimental
# condition compared to the Control_No_Tx reference group.
# Also runs Wilcoxon rank-sum tests with BH correction for significance.
#
# PREREQUISITE: Run Figure_1G_Process_Data.R first (creates .RData file)
# ============================================================================

library(here)
library(dplyr)
library(effsize)   # for cohen.d() with hedges.correction

# ============================================================================
# LOAD DATA
# ============================================================================
data_file <- here("Data", "Fibroblast_lifespan", "Processed", "Figure_1G_data.RData")

if (!file.exists(data_file)) {
  stop(
    "\nProcessed data not found: ", data_file,
    "\nRun Data Processing Scripts/Figure_1G_Process_Data.R first.\n"
  )
}

load(data_file)
boxplot_data <- figure_1g_data$boxplot_data
message("Loaded ", nrow(boxplot_data), " samples across ",
        length(unique(boxplot_data$intx)), " conditions")

# ============================================================================
# CONFIGURATION
# ============================================================================
reference_group <- "Control_No_Tx"

# Output directory
output_dir <- here("Results", "Figures", "Figure_1G")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ============================================================================
# COMPUTE HEDGES' g AND WILCOXON P-VALUES
# ============================================================================
ref_scores <- boxplot_data %>%
  filter(intx == reference_group) %>%
  pull(ISR_Score)

conditions <- setdiff(unique(boxplot_data$intx), reference_group)

results <- lapply(conditions, function(cond) {
  cond_scores <- boxplot_data %>%
    filter(intx == cond) %>%
    pull(ISR_Score)

  # Hedges' g (bias-corrected effect size)
  hg <- cohen.d(cond_scores, ref_scores, hedges.correction = TRUE)

  # Wilcoxon rank-sum test (non-parametric, matches existing Kruskal approach)
  wt <- wilcox.test(cond_scores, ref_scores, exact = FALSE)

  data.frame(
    Condition        = cond,
    Reference        = reference_group,
    N_Condition      = length(cond_scores),
    N_Reference      = length(ref_scores),
    Mean_Condition    = mean(cond_scores, na.rm = TRUE),
    Mean_Reference    = mean(ref_scores, na.rm = TRUE),
    Hedges_g         = hg$estimate,
    CI_Lower         = hg$conf.int[1],
    CI_Upper         = hg$conf.int[2],
    Magnitude        = as.character(hg$magnitude),
    Wilcoxon_W       = wt$statistic,
    P_Value_Raw      = wt$p.value,
    stringsAsFactors = FALSE
  )
})

results_df <- bind_rows(results)

# BH-adjusted p-values
results_df$P_Value_Bonferroni <- p.adjust(results_df$P_Value_Raw, method = "bonferroni")

# Significance stars
results_df <- results_df %>%
  mutate(Significance = case_when(
    P_Value_Bonferroni < 0.001  ~ "***",
    P_Value_Bonferroni < 0.01   ~ "**",
    P_Value_Bonferroni < 0.05   ~ "*",
    TRUE                 ~ "ns"
  ))

# Sort by Hedges' g
results_df <- results_df %>% arrange(Hedges_g)

# ============================================================================
# PRINT SUMMARY
# ============================================================================
message("\n", strrep("=", 70))
message("HEDGES' g: Each Condition vs ", reference_group)
message(strrep("=", 70))
for (i in seq_len(nrow(results_df))) {
  r <- results_df[i, ]
  message(sprintf("  %-35s  g = %6.3f  [%6.3f, %6.3f]  %s  (p_adj = %.4g)",
                  r$Condition, r$Hedges_g, r$CI_Lower, r$CI_Upper,
                  r$Significance, r$P_Value_Bonferroni))
}
message(strrep("=", 70))

# ============================================================================
# SAVE CSV
# ============================================================================
csv_path <- file.path(output_dir, "Figure_1G_HedgesG_vs_Control.csv")
write.csv(results_df, csv_path, row.names = FALSE)
message("\nSaved: ", csv_path)
