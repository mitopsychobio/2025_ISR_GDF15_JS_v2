# ============================================================================
# Figure 1G: ISR Score Boxplots - FIGURE GENERATION SCRIPT
# ============================================================================
# This script generates boxplots of ISR (Factor 1) scores for each condition,
# sorted by ascending median score, with Kruskal-Wallis test and Dunn's
# post-hoc comparisons.
#
# PREREQUISITE: Run 01_Process_Data.R first
# ============================================================================

# Load required packages
library(here)
library(ggplot2)
library(ggpubr)
library(dplyr)

# ============================================================================
# CONFIGURATION
# ============================================================================
# Input
data_file <- here("Data", "Fibroblast_lifespan", "Processed", "Figure_1G_data.RData")

# Output
output_dir <- here("Results", "Figures", "Figure_1G")
output_file <- file.path(output_dir, "Figure_1G_ISR_Boxplots.png")

# ============================================================================
# CHECK PREREQUISITE
# ============================================================================
if (!file.exists(data_file)) {
  stop(
    "\n",
    "========== ERROR: PROCESSED DATA NOT FOUND ==========\n",
    "File missing: ", data_file, "\n",
    "\n",
    "You must run 01_Process_Data.R first.\n",
    "========================================================\n"
  )
}

# ============================================================================
# LOAD DATA
# ============================================================================
message("Loading processed data...")
load(data_file)

boxplot_data <- figure_1g_data$boxplot_data
condition_stats <- figure_1g_data$condition_stats
kw_test <- figure_1g_data$kw_test
pairwise_comparisons <- figure_1g_data$pairwise_comparisons

message("  Samples: ", nrow(boxplot_data))
message("  Conditions: ", length(levels(boxplot_data$intx)))

# ============================================================================
# CONFIGURE COLORS
# ============================================================================
condition_colors <- c(
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

# Derive slightly darker shades for box outlines, median lines, and whiskers
darken_hex <- function(color, factor = 0.65) {
  v <- col2rgb(color) / 255
  rgb(v[1, ] * factor, v[2, ] * factor, v[3, ] * factor)
}
condition_colors_dark <- setNames(
  sapply(condition_colors, darken_hex, factor = 0.65),
  names(condition_colors)
)

# ============================================================================
# CREATE BOXPLOT
# ============================================================================
message("\nGenerating boxplot...")

# Format condition labels for x-axis
boxplot_data <- boxplot_data %>%
  mutate(intx_label = gsub("Control_", "Control + ", intx)) %>%
  mutate(intx_label = gsub("SURF1_", "SURF1 + ", intx_label)) %>%
  mutate(intx_label = gsub("_", " ", intx_label))

# Reorder labels
label_order <- condition_stats$intx
label_order_formatted <- gsub("Control_", "Control + ", label_order)
label_order_formatted <- gsub("SURF1_", "SURF1 + ", label_order_formatted)
label_order_formatted <- gsub("_", " ", label_order_formatted)

boxplot_data$intx_label <- factor(boxplot_data$intx_label, levels = label_order_formatted)

# Create base boxplot
p <- ggplot(boxplot_data, aes(x = intx_label, y = ISR_Score, fill = intx, color = intx)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.4, size = 1.5) +
  scale_fill_manual(values = condition_colors, guide = "none") +
  scale_color_manual(values = condition_colors_dark, guide = "none") +
  coord_cartesian(ylim = c(-2, 6)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  labs(
    title = "ISR Score by Condition",
    subtitle = paste0("Kruskal-Wallis p ", if(kw_test$p.value < 0.0001) "< 0.0001" else paste("=", format(kw_test$p.value, digits = 3))),
    x = "",
    y = "ISR Score"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 10, face = "italic"),
    axis.title.y = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    axis.text.y = element_text(size = 10),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )

# Snapshot before significance brackets are added (used for no-stars variant)
p_no_stars <- p

# ============================================================================
# Add significance brackets — Kruskal-Wallis + Dunn's test (all 91 pairs,
# Bonferroni corrected). Only DISPLAY brackets for significant comparisons
# vs Control_No_Tx.
# ============================================================================
library(ggsignif)

# Helper: convert Dunn's Bonferroni-adjusted p-value to asterisks
p_to_stars <- function(p) {
  dplyr::case_when(
    p < 0.0001 ~ "****",
    p < 0.001  ~ "***",
    p < 0.01   ~ "**",
    p < 0.05   ~ "*",
    TRUE       ~ "ns"
  )
}

# Format group names to match x-axis labels
format_label <- function(x) {
  x <- gsub("Control_", "Control + ", x)
  x <- gsub("SURF1_", "SURF1 + ", x)
  x <- gsub("_", " ", x)
  return(x)
}

# --- From the full Dunn's test (91 pairs, Bonferroni), extract only
#     Control_No_Tx comparisons that are significant (P.adjusted < 0.05) ---
control_ref <- "Control_No_Tx"

control_comparisons <- pairwise_comparisons %>%
  mutate(
    group1 = trimws(sub(" - .*", "", comparison)),
    group2 = trimws(sub(".* - ", "", comparison))
  ) %>%
  filter(group1 == control_ref | group2 == control_ref) %>%
  filter(P.adjusted < 0.05) %>%
  mutate(stars = p_to_stars(P.adjusted))

message("  Dunn's test (Bonferroni, all 91 pairs): ",
        nrow(control_comparisons), " significant vs ", control_ref)

if (nrow(control_comparisons) > 0) {
  # Get x-axis positions from factor levels
  condition_levels <- levels(boxplot_data$intx_label)

  control_comparisons <- control_comparisons %>%
    mutate(
      label1 = format_label(group1),
      label2 = format_label(group2),
      x1 = match(label1, condition_levels),
      x2 = match(label2, condition_levels),
      span = abs(x2 - x1)
    ) %>%
    # Order: shortest span (bottom) to longest span (top)
    arrange(span)

  # Assign y-positions: shortest bracket at the bottom, longest at the top
  y_max <- max(boxplot_data$ISR_Score, na.rm = TRUE)
  y_start <- y_max + 0.5
  y_step  <- y_max * 0.14

  control_comparisons <- control_comparisons %>%
    mutate(y_pos = y_start + (row_number() - 1) * y_step)

  # Print comparison details for verification
  message("\n  Bracket details (Dunn's test, Bonferroni across all 91 pairs):")
  for (i in seq_len(nrow(control_comparisons))) {
    message("    ", control_comparisons$comparison[i],
            "  Z=", round(control_comparisons$Z[i], 3),
            "  P.adj=", format(control_comparisons$P.adjusted[i], digits = 3),
            "  ", control_comparisons$stars[i])
  }

  # Add each bracket individually (more reliable than single data frame)
  for (i in seq_len(nrow(control_comparisons))) {
    bracket_df <- data.frame(
      xmin        = control_comparisons$x1[i],
      xmax        = control_comparisons$x2[i],
      y_position  = control_comparisons$y_pos[i],
      annotations = control_comparisons$stars[i],
      stringsAsFactors = FALSE
    )
    p <- p + geom_signif(
      data         = bracket_df,
      aes(xmin = xmin, xmax = xmax, annotations = annotations, y_position = y_position),
      tip_length   = 0.01,
      textsize     = 4.5,
      manual       = TRUE,
      inherit.aes  = FALSE
    )
  }

  message("\n  Total brackets added to plot: ", nrow(control_comparisons))
}

print(p)


# ============================================================================
# SAVE NO-STARS VARIANT
# ============================================================================
output_file_no_stars     <- file.path(output_dir, "FB_ISR_Scores_No_Stars.png")

ggsave(output_file_no_stars,     plot = p_no_stars, width = 14, height = 8, dpi = 300)
message("Saved: ", output_file_no_stars)


# ============================================================================
# EXPORT SOURCE DATA
# ============================================================================
message("\nExporting source data for Figure 1G...")

# Export wide-format CSV for Prism: intx as column headers, ISR_Score as rows
boxplot_export <- boxplot_data %>%
  select(intx, ISR_Score) %>%
  group_by(intx) %>%
  mutate(row_id = row_number()) %>%
  ungroup() %>%
  pivot_wider(names_from = intx, values_from = ISR_Score) %>%
  select(-row_id)

# Reorder columns to match the boxplot order (ascending median)
boxplot_export <- boxplot_export[, as.character(condition_stats$intx)]

csv_file <- file.path(output_dir, "FB_Factor1_for_Prism.csv")
write.csv(boxplot_export, file = csv_file, row.names = FALSE, na = "")
message("Saved source data: ", csv_file)

csv_file2 <- file.path(output_dir, "significant_comparisons.csv")
write.csv(control_comparisons, file = csv_file2, row.names = FALSE, na = "")
message("Saved source data: ", csv_file2)
