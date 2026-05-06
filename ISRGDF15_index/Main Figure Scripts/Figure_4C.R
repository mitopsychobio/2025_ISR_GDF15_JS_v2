#
# Control F and then replace "Factor1" with Factor_
rm(list = ls())

# BiocManager::install("edgeR")
?effsize::cohen.d
library(tidyverse)
library(edgeR)
library(corrr)
library(tibble)
library(ggpubr)
library(readr)
library(circlize)
library(beepr)
library(ggplot2)
library(dplyr)
library(effsize)
library(here)


condition_column_names <- c("DTHPLCE")
control <- "Hospital inpatient"
condition <- "Emergency room"
condition_for_titles <- "deathplace"
method_multiple_comparisons <- "BH"
date_of_csv_folder <- "2024_12_13"
CohenFalse_HedgesTrue <- "TRUE"


#------------------------------------------------------
#------------------------------------------------------ NAMING FOLDER AFTER THE NAME OF THE SCRIPT

# Output directories
local_out_dir <- here("Results", "Figures", "Figure_4C")
shared_out_dir <- here("Results", "Figures")
if (!dir.exists(local_out_dir)) dir.create(local_out_dir, recursive = TRUE)
if (!dir.exists(shared_out_dir)) dir.create(shared_out_dir, recursive = TRUE)
folder_path <- local_out_dir  # Keep folder_path for compatibility with rest of script

#------------------------------------------------------

Plot_Save = "ON"
#------------------------------------------------------


# Build the per-donor / per-tissue long-format dataframe by reading the
# per-tissue Plot_Data_<tissue>.csv files written by
# Figure_3B_4A_Tissue_Scatter.R. The helper returns columns
#   SAMPID, Tissue, DTHPLCE, fa_vs_gdf15, value, group, ...
# with NO DTHPLCE filtering. The Hospital-inpatient / Emergency-room
# restriction for Figure 4C is applied below via the `comparisons` list.
source(here("Main Figure Scripts", "Helper Scripts",
            "Load_All_Tissue_DTHPLCE.R"))
data <- load_all_tissue_dthplce()

# Restrict to the two place-of-death categories Figure 4C compares.
# (Helper passes through every DTHPLCE category; small ones like "At MVA
# scene" have N < 20 and would cause the all(N >= 20) check below to drop
# every tissue if left in, leaving final_p_values empty and crashing the
# downstream `mutate(Group = paste(Tissue, Comparison, ...))`.)
data <- data %>%
  filter(DTHPLCE %in% c("Hospital inpatient", "Emergency room"))

#----------------------------------------------------------------------------------------------------------

# Assuming your dataframe is named 'data'
data <- data %>%
  group_by(Tissue, group) %>%
  mutate(N = n())

data <- data %>%
  group_by(Tissue) %>%
  filter(all(N >= 20))


# Initialize result storage
results_list <- list()
plots_list <- list()
cohen_d_results <- data.frame()


# Define the correct group names for comparisons
comparisons <- list(
  c("Hospital inpatient_Factor1", "Emergency room_Factor1"),
  c("Hospital inpatient_GDF15", "Emergency room_GDF15")
)


# Iterate through each tissue
for (tissue in unique(data$Tissue)) {
  tissue_data <- data %>% filter(Tissue == tissue)

  # Initialize a vector to store unadjusted p-values
  tissue_p_values <- c()
  comparison_names <- c()

  for (comp in comparisons) {
    group1 <- comp[1]
    group2 <- comp[2]

    # Subset data for the comparison
    subset_data <- tissue_data %>% filter(group %in% c(group1, group2))

    if (nrow(subset_data) > 0) {
      # Set group levels explicitly to ensure order
      subset_data$group <- factor(subset_data$group, levels = c(group1, group2))

      # Perform Wilcoxon test
      wilcox_test <- wilcox.test(value ~ group, data = subset_data, exact = FALSE)

      # Store the unadjusted p-value and corresponding comparison
      tissue_p_values <- c(tissue_p_values, wilcox_test$p.value)
      comparison_names <- c(comparison_names, paste(group1, "vs", group2))



    } else {
      print(paste("No data for comparison:", group1, "vs", group2, "in tissue", tissue))
    }
  }

  # Apply p-value adjustment for multiple comparisons within this tissue
  if (length(tissue_p_values) > 0) {
    adjusted_p_values <- p.adjust(tissue_p_values, method = method_multiple_comparisons)

    # Store results for each comparison with adjusted p-values
    for (i in seq_along(comparison_names)) {
      comp <- strsplit(comparison_names[i], " vs ")[[1]]
      group1 <- comp[1]
      group2 <- comp[2]

      subset_data <- tissue_data %>% filter(group %in% c(group1, group2))

      if (nrow(subset_data) > 0) {
        # Set group levels explicitly to ensure order
        subset_data$group <- factor(subset_data$group, levels = c(group1, group2))

        # Calculate Cohen's D
        cohen_d <- cohen.d(subset_data$value, subset_data$group, hedges.correction = CohenFalse_HedgesTrue)

        # Save results
        p_values <- data.frame(
          Tissue = tissue,
          Comparison = comparison_names[i],
          Unadjusted_P = tissue_p_values[i],
          Adjusted_P = adjusted_p_values[i]
        )
        results_list[[paste(tissue, comparison_names[i], sep = "_")]] <- p_values

        cohen_d_results <- rbind(cohen_d_results, data.frame(
          Tissue = tissue,
          Comparison = comparison_names[i],
          Cohens_D = cohen_d$estimate,
          Adjusted_P = adjusted_p_values[i]
        ))

        # Create plot
        plot <- ggplot(subset_data, aes(x = group, y = value)) +
          geom_boxplot() +
          geom_jitter(width = 0.2, alpha = 0.5) +
          ggtitle(paste("Tissue:", tissue, "-", group1, "vs", group2)) +
          annotate("text", x = 1.5, y = max(subset_data$value, na.rm = TRUE),
                   label = paste("Adj P:", signif(adjusted_p_values[i], 4), "\nMethod: ", method_multiple_comparisons),
                   hjust = 0.5) +
          theme_minimal()

        plots_list[[paste(tissue, group1, group2, sep = "_")]] <- plot
      }
    }
  }
}


# Combine all p-values and adjusted p-values into a single dataframe
final_p_values <- bind_rows(results_list)

if (length(final_p_values) > 0) {
  final_p_values$Adjusted_P <- p.adjust(final_p_values$Unadjusted_P, method = method_multiple_comparisons)
}


cohen_d_results <- cohen_d_results %>%
  mutate(Group = paste(Tissue, Comparison, sep = "_"))

cohen_d_results <- cohen_d_results %>%
  select(-Adjusted_P, -Comparison)

final_p_values <- final_p_values %>%
  mutate(Group = paste(Tissue, Comparison, sep = "_"))

final_p_values <- final_p_values %>%
  select(-Comparison, -Tissue)

cohen_d_results <- cohen_d_results %>%
  left_join(final_p_values %>% select(Group, Adjusted_P), by = "Group")

cohen_d_results <- cohen_d_results %>%
  mutate(Comparison = str_extract(Group, "Hospital.*"))


# Apply p-value adjustment for multiple comparisons within this tissue
if (length(tissue_p_values) > 0) {
  adjusted_p_values <- p.adjust(tissue_p_values, method = method_multiple_comparisons)

  # Store results for each comparison with adjusted p-values
  for (i in seq_along(comparison_names)) {
    comp <- strsplit(comparison_names[i], " vs ")[[1]]
    group1 <- comp[1]
    group2 <- comp[2]

    subset_data <- tissue_data %>% filter(group %in% c(group1, group2))

    if (nrow(subset_data) > 0) {

     # Create plot
      plot <- ggplot(subset_data, aes(x = group, y = value)) +
        geom_boxplot() +
        geom_jitter(width = 0.2, alpha = 0.5) +
        ggtitle(paste("Tissue:", tissue, "-", group1, "vs", group2)) +
        annotate("text", x = 1.5, y = max(subset_data$value, na.rm = TRUE),
                 label = paste("Adj P:", signif(adjusted_p_values[i], 4), "\nMethod: ", method_multiple_comparisons),
                 hjust = 0.5) +
        theme_minimal()

      plots_list[[paste(tissue, group1, group2, sep = "_")]] <- plot
    }
  }
}




# Save each plot in the specified folder
for (name in names(plots_list)) {
  ggsave(
    filename = file.path(folder_path, paste0(name, ".png")),
    plot = plots_list[[name]],
    width = 6,
    height = 4
  )
}

# Create Cohen's D plots with significance coloring
cohen_d_results <- cohen_d_results %>%
  mutate(Significant = ifelse(Adjusted_P < 0.05, "Significant", "Not Significant"))

# Plot for Factor1 comparisons
cohen_d_factor1 <- cohen_d_results %>%
  filter(grepl("Factor1", Group)) %>%
  arrange(Cohens_D)

plot_factor1 <- ggplot(cohen_d_factor1, aes(x = Cohens_D, y = reorder(Tissue, Cohens_D), color = Significant)) +
  geom_point(size = 3) +
  scale_color_manual(values = c("Significant" = "red", "Not Significant" = "black")) +
  ggtitle("Hedge's G for Factor1 Comparisons") +
  xlab("Hedge's G") + ylab("Tissue") +
  theme_minimal()

# Plot for GDF15 comparisons
cohen_d_gdf15 <- cohen_d_results %>%
  filter(grepl("GDF15", Group)) %>%
  arrange(Cohens_D)

plot_gdf15 <- ggplot(cohen_d_gdf15, aes(x = Cohens_D, y = reorder(Tissue, Cohens_D), color = Significant)) +
  geom_point(size = 3) +
  scale_color_manual(values = c("Significant" = "red", "Not Significant" = "black")) +
  ggtitle("Hedge's G for GDF15 Comparisons") +
  xlab("Hedge's G") + ylab("Tissue") +
  theme_minimal()

# Display Cohen's D plots
print(plot_factor1)
print(plot_gdf15)

# Save Cohen's D plots
ggsave(filename = file.path(paste0(folder_path, "/1_Hedges_G_Factor1_", method_multiple_comparisons, ".png")), plot = plot_factor1, width = 6, height = 8)
ggsave(filename = file.path(paste0(folder_path, "/1_Hedges_G_GDF15_", method_multiple_comparisons, ".png")), plot = plot_gdf15, width = 6, height = 8)


# Display or save final results
# p_values_results removed (redundant with Hedges_G_results / source_data)
write.csv(cohen_d_results, file = file.path(folder_path, "1_Hedges_G_results.csv"), row.names = FALSE)

# Combine Factor1 and GDF15 data for combined plot
combined_data <- cohen_d_results %>%
  mutate(
    Type = ifelse(grepl("Factor1", Group), "Factor1", "GDF15"),
    Shape = ifelse(Type == "Factor1", "Circle", "Square"),
    Color = ifelse(Significant == "Significant", "Red", "Gray")
  ) %>%
  arrange(Cohens_D) %>%
  mutate(Tissue = factor(Tissue, levels = unique(cohen_d_factor1$Tissue)))

# ---- Diagnostic: tell us what's missing per shape so we can spot why
# squares were silently dropped before. ----
message("\n=== combined_data diagnostics ===")
message("Total rows: ", nrow(combined_data))
for (sh in c("Circle", "Square")) {
  sub <- combined_data %>% filter(Shape == sh)
  message("  Shape '", sh, "': ", nrow(sub), " rows | ",
          "NA Cohens_D: ", sum(is.na(sub$Cohens_D)), " | ",
          "NA Adjusted_P: ", sum(is.na(sub$Adjusted_P)), " | ",
          "NA Significant: ", sum(is.na(sub$Significant)), " | ",
          "NA Color: ", sum(is.na(sub$Color)), " | ",
          "NA Tissue (after factor): ", sum(is.na(sub$Tissue)))
  bad <- sub %>% filter(is.na(Cohens_D) | is.na(Color) | is.na(Tissue))
  if (nrow(bad) > 0) {
    message("    First few problematic rows:")
    print(utils::head(as.data.frame(bad), 5))
  }
}

# Defensive filter: drop rows whose plotted aesthetics are NA so geom_point
# doesn't silently swallow the entire layer.
n_before <- nrow(combined_data)
combined_data <- combined_data %>%
  filter(!is.na(Cohens_D), !is.na(Tissue), !is.na(Color))
n_after <- nrow(combined_data)
if (n_before != n_after) {
  message("Dropped ", n_before - n_after,
          " row(s) with NA in Cohens_D / Tissue / Color before plotting.")
}

# Create combined plot
# Outline rule: Factor1 circles get a black outline; GDF15 squares get a
# transparent outline (so the square fill is visible without a contrasting
# rim). The outline colour is encoded as a literal column on the data and
# pulled in with scale_color_identity() — that avoids the ggplot2 quirk
# where a discrete colour scale containing NA silently drops every point
# with NA outline.
combined_data <- combined_data %>%
  mutate(OutlineColor = ifelse(Shape == "Circle", "black", "transparent"))

combined_plot <- ggplot(combined_data, aes(x = Cohens_D, y = Tissue)) +
  geom_point(
    aes(shape = Shape, fill = Color, color = OutlineColor),
    size = 6, alpha = 0.6, stroke = 1
  ) +
  scale_shape_manual(values = c("Circle" = 21, "Square" = 22)) +
  scale_fill_manual(values  = c("Red" = "red", "Gray" = "gray")) +
  scale_color_identity() +   # use the literal colours in OutlineColor
  geom_vline(xintercept = 0, linetype = "dotted", color = "black", size = 0.5) +
  ggtitle("Hedge's G for Factor1 and GDF15 Comparisons") +
  xlab("Hedge's G") +
  ylab("Tissue") +
  theme_minimal() +
  theme(
    axis.line = element_line(size = 0.25),  # Thinner axis lines
    panel.grid.major = element_line(size = 0.25),  # Thinner major grid lines
    panel.grid.minor = element_line(size = 0.1)  # Thinner minor grid lines
  )

# Display the plot
print(combined_plot)

# Save the combined plot
ggsave(filename = file.path(paste0(folder_path, "/1_Hedges_G_Combined_", method_multiple_comparisons, ".png")), plot = combined_plot, width = 14, height = 10)

# Also save to shared output directory
ggsave("Figure_4C.png", plot = combined_plot, path = shared_out_dir, width = 14, height = 10)


# CHI SQUARE
# Step 1: Create the contingency table
contingency_table <- cohen_d_results %>%
  group_by(Comparison, Significant) %>%
  summarise(Count = n(), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = Significant, values_from = Count, values_fill = 0)

# Rename columns for clarity
contingency_table <- contingency_table %>%
  rename(
    Not_Significant = `Not Significant`,
    Significant = Significant
  )

# Step 2: Perform the chi-square test
chisq_result <- chisq.test(contingency_table[, c("Not_Significant", "Significant")])

# Print the results
print("Contingency Table:")
print(contingency_table)

print("Chi-Square Test Result:")
print(chisq_result)


# NOTE: An exploratory AGE/SEX block (Averages by tissue and per-DTHPLCE
# sex_counts) lived here previously. It only printed to the console — no
# saved figure or CSV depended on it. It has been removed because the
# per-tissue Plot_Data CSVs that now feed this script don't carry AGE / SEX
# columns. If those summaries are ever needed again, join AGE/SEX from the
# GTEx phenotype annotations file before calling load_all_tissue_dthplce().

# ============================================================================
# SAVE SOURCE DATA FOR FIGURE 4C
# ============================================================================
message("\nExporting source data for Figure 4C...")

# Source data already saved as 1_Hedges_G_results.csv above
message("  Figure 4C source data: 1_Hedges_G_results.csv")
message("  Location: ", folder_path)
