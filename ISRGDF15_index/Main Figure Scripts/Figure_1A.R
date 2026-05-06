# ============================================================================
# Figure 1A: Venn Diagram - FIGURE GENERATION SCRIPT
# ============================================================================
# This script generates the Venn diagram showing overlap between three
# ISR gene list sources: GO Consortium, AnyGenes, and JacksonLabs
#
# PREREQUISITE: Run 01_Process_Data.R first to generate processed gene lists
# ============================================================================

# Load required packages
library(here)
library(VennDiagram)
library(grid)
library(futile.logger)  # Suppress VennDiagram logging

# Suppress VennDiagram log messages
flog.threshold(ERROR)

# ============================================================================
# CONFIGURATION
# ============================================================================
# Input: processed gene lists
data_file <- here("Data", "ISR_Gene_Lists", "processed_gene_lists.RData")

# Output
output_dir <- here("Results", "Figures", "Figure_1A")
output_file <- file.path(output_dir, "Figure_1A_Venn_Diagram.png")

# ============================================================================
# CHECK PREREQUISITE - FAIL LOUD IF DATA NOT PROCESSED
# ============================================================================
if (!file.exists(data_file)) {
  stop(
    "\n",
    "========== ERROR: PROCESSED DATA NOT FOUND ==========\n",
    "File missing: ", data_file, "\n",
    "\n",
    "You must run 01_Process_Data.R first to generate the processed gene lists.\n",
    "========================================================\n"
  )
}

# ============================================================================
# LOAD PROCESSED DATA
# ============================================================================
message("Loading processed gene lists...")
load(data_file)

# Extract gene lists
genes_GO <- gene_lists$GO_consortium
genes_AG <- gene_lists$AnyGenes
genes_JL <- gene_lists$JacksonLabs

message("Gene counts:")
message("  GO Consortium: ", length(genes_GO))
message("  AnyGenes: ", length(genes_AG))
message("  JacksonLabs: ", length(genes_JL))

# ============================================================================
# CALCULATE VENN REGIONS
# ============================================================================
# Calculate overlaps for verification
only_GO <- setdiff(setdiff(genes_GO, genes_AG), genes_JL)
only_AG <- setdiff(setdiff(genes_AG, genes_GO), genes_JL)
only_JL <- setdiff(setdiff(genes_JL, genes_GO), genes_AG)
GO_AG <- setdiff(intersect(genes_GO, genes_AG), genes_JL)
GO_JL <- setdiff(intersect(genes_GO, genes_JL), genes_AG)
AG_JL <- setdiff(intersect(genes_AG, genes_JL), genes_GO)
all_three <- intersect(intersect(genes_GO, genes_AG), genes_JL)

total_genes <- length(unique(c(genes_GO, genes_AG, genes_JL)))

message("\nVenn diagram regions:")
message("  Only GO: ", length(only_GO))
message("  Only AnyGenes: ", length(only_AG))
message("  Only JacksonLabs: ", length(only_JL))
message("  GO & AnyGenes only: ", length(GO_AG))
message("  GO & JacksonLabs only: ", length(GO_JL))
message("  AnyGenes & JacksonLabs only: ", length(AG_JL))
message("  All three: ", length(all_three))
message("  Total unique: ", total_genes)

# ============================================================================
# GENERATE VENN DIAGRAM
# ============================================================================
message("\nGenerating Venn diagram...")

# Colors matching the reference figure (translucent red, yellow, blue)
# Based on Figure 1A in figures_smaller.pdf
venn_colors <- c("#E74C3C", "#F1C40F", "#3498DB")  # red, yellow, blue

# Create PNG
png(output_file, width = 8, height = 8, units = "in", res = 300)

venn_plot <- draw.triple.venn(
  area1 = length(genes_GO),
  area2 = length(genes_AG),
  area3 = length(genes_JL),
  n12 = length(intersect(genes_GO, genes_AG)),
  n23 = length(intersect(genes_AG, genes_JL)),
  n13 = length(intersect(genes_GO, genes_JL)),
  n123 = length(all_three),
  category = c("GO Consortium", "AnyGenes", "JacksonLabs"),
  fill = venn_colors,
  alpha = 0.5,
  cex = 2,
  cat.cex = 1.5,
  cat.fontface = "bold",
  cat.dist = c(0.05, 0.05, 0.05),
  cat.pos = c(-30, 30, 180),
  margin = 0.1,
  scaled = FALSE,
  euler.d = FALSE,
  print.mode = c("raw", "percent"),
  sigdigs = 3
)

# Add title
grid.text("Consensus ISR gene lists",
          x = 0.5, y = 0.95,
          gp = gpar(fontsize = 16, fontface = "bold"))

dev.off()
message("Saved: ", output_file)


# ============================================================================
# EXPORT SOURCE DATA
# ============================================================================
message("\nExporting source data for Figure 1A...")

# Create a membership matrix for the Venn diagram
all_genes <- unique(c(genes_GO, genes_AG, genes_JL))
venn_matrix <- data.frame(
  Gene = all_genes,
  GO_Consortium = all_genes %in% genes_GO,
  AnyGenes = all_genes %in% genes_AG,
  JacksonLabs = all_genes %in% genes_JL
)

# Sort by gene name
venn_matrix <- venn_matrix[order(venn_matrix$Gene), ]

# Save CSV
csv_file <- file.path(output_dir, "Figure_1A_source_data.csv")
write.csv(venn_matrix, file = csv_file, row.names = FALSE)
message("Saved source data: ", csv_file)
