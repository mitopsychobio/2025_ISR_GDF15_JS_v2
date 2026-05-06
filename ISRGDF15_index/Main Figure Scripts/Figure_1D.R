# Figure 1D: PCA Biplot of Fibroblast Lifespan Data
# Shows 339 samples colored by treatment condition

library(here)
library(ggplot2)
library(dplyr)
library(tidyr)

# Source the data processing script to get the datasets
# NOTE: Must source BEFORE defining out_dir because Intro1 cleans up the environment
source(here("Main Figure Scripts", "Helper Scripts", "Fibroblast_lifespan", "Intro1_FB_Lifespan_exprs_manifest_ISR_list.R"), local = TRUE)

# Load gene list
source(here("Main Figure Scripts", "Helper Scripts", "Fibroblast_lifespan", "GO_vs_AnyGenes_vs_JacksonLabs_ListsofGenes.R"), local = TRUE)

# Create output directories (defined AFTER source() calls due to Intro1 cleanup)
local_out_dir <- here("Results", "Figures", "Figure_1D")
shared_out_dir <- here("Results", "Figures")
if (!dir.exists(local_out_dir)) dir.create(local_out_dir, recursive = TRUE)
if (!dir.exists(shared_out_dir)) dir.create(shared_out_dir, recursive = TRUE)

# Process datasets
genes_of_interest <- genes_of_interest_plus_gdf15

data_sub <- data_sub %>%
  mutate(intx = paste(Group, Experiment, sep = "_"))

datasets <- process_datasets(data_sub, genes_of_interest)
WT_SURF1_all_Txs <- datasets$WT_SURF1_all_Txs

# Prepare data for PCA
valid_genes <- genes_of_interest[genes_of_interest %in% colnames(WT_SURF1_all_Txs)]
gene_data <- WT_SURF1_all_Txs[, valid_genes]

# Scale data
data_scaled <- scale(gene_data)

# Perform PCA
pca_result <- prcomp(data_scaled, center = FALSE, scale. = FALSE)  # already scaled

# Get PC scores
pc_scores <- as.data.frame(pca_result$x[, 1:2])
pc_scores$intx <- WT_SURF1_all_Txs$intx

# Get variance explained
var_explained <- summary(pca_result)$importance[2, 1:2] * 100

# Color palette - EXACT SAME as Figure 1F in 20251126_Spearman_rho_Gene_vs_DaysGrown.Rmd
color_palette <- c(
  "Control_No_Tx" = "gray",
  "Control_DEX" = "lightcoral",
  "SURF1_No_Tx" = "purple",
  "SURF1_DEX" = "maroon",
  "Control_Oligo" = "violet",
  "Control_mitoNUITs" = "orange",
  "Control_DEX_mitoNUITs" = "lightpink",
  "Control_DEX_Oligo" = "blue",
  "Control_ox3" = "lightgreen",
  "Control_Contact_Inhibition" = "green",
  "Control_Contact_Inhibition_ox3" = "darkgreen",
  "Control_Galactose" = "lightblue",
  "Control_betahydroxybutyrate" = "skyblue",
  "Control_2DG" = "yellow"
)

# Create simple PCA scatter plot (no arrows, no GDF15 label)
p <- ggplot(pc_scores, aes(x = PC1, y = PC2, color = intx)) +
  geom_point(alpha = 0.7, size = 4) +
  scale_color_manual(values = color_palette, name = "Intx") +
  labs(
    title = paste0("PCA Biplot: All data n= ", nrow(pc_scores),
                   " , gene input ", length(valid_genes), " ISR genes - WT_SURF1_all_Txs"),
    x = paste0("PC1 (", round(var_explained[1], 2), "%)"),
    y = paste0("PC2 (", round(var_explained[2], 2), "%)")
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 8),
    plot.title = element_text(size = 10)
  ) +
  coord_fixed()

# Save to local folder
ggsave(
  filename = file.path(local_out_dir, "Figure_1D_PCA_Biplot.png"),
  plot = p,
  width = 14,
  height = 10,
  dpi = 300
)

# Save to shared Output folder
ggsave(
  filename = file.path(shared_out_dir, "Figure_1D.png"),
  plot = p,
  width = 14,
  height = 10,
  dpi = 300
)

# ============================================================================
# EXPORT SOURCE DATA CSV
# ============================================================================
# Save PCA scores (PC1, PC2, treatment condition) and variance explained
source_data_1D <- pc_scores %>%
  select(PC1, PC2, intx)

write.csv(source_data_1D,
          file.path(local_out_dir, "Figure_1D_source_data.csv"),
          row.names = FALSE)
message("  Saved: Figure_1D_source_data.csv (", nrow(source_data_1D), " rows)")

# Also save variance explained as a separate small CSV
var_df <- data.frame(
  Component = c("PC1", "PC2"),
  Variance_Explained_Pct = round(var_explained, 2)
)
write.csv(var_df,
          file.path(local_out_dir, "Figure_1D_variance_explained.csv"),
          row.names = FALSE)
message("  Saved: Figure_1D_variance_explained.csv")

message("Figure 1D saved to: ", local_out_dir)
message("Figure 1D also saved to: ", shared_out_dir)
