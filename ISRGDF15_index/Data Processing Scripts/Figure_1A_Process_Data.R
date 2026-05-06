# ============================================================================
# Figure 1A: Venn Diagram - DATA PROCESSING SCRIPT
# ============================================================================
# This script processes the ISR gene lists from three sources:
# - GO Consortium
# - AnyGenes
# - JacksonLabs (mouse genes - converted to human)
#
# Output: Processed gene lists saved to Figures/Data/ISR_Gene_Lists/
# ============================================================================

# Load required packages
library(here)
library(tidyverse)

# ============================================================================
# CONFIGURATION
# ============================================================================
# Input file
input_file <- here("Data", "ISR_Gene_Lists", "ISR_Gene_Lists_updated.csv")

# Output directory
output_dir <- here("Data", "ISR_Gene_Lists")

# ============================================================================
# MOUSE TO HUMAN GENE NAME MAPPING
# ============================================================================
# Some mouse genes have different naming conventions than human genes
# This mapping handles known differences
mouse_to_human_mapping <- c(
  "IRGM" = "IRGM",      # igtp in mouse -> IRGM in human (already handled by uppercase)
  "NARS1" = "NARS1",    # Nars -> NARS1
  "WARS1" = "WARS1"     # Wars -> WARS1
)

convert_mouse_to_human <- function(gene_name) {
  # First convert to uppercase (standard human gene format)
  human_gene <- toupper(gene_name)

 # Check if there's a specific mapping needed
  if (human_gene %in% names(mouse_to_human_mapping)) {
    return(mouse_to_human_mapping[human_gene])
  }

  return(human_gene)
}

# ============================================================================
# DATA PROCESSING
# ============================================================================
message("Reading ISR gene lists from: ", input_file)

# Check input file exists - FAIL LOUD if not
if (!file.exists(input_file)) {
  stop("INPUT FILE NOT FOUND: ", input_file, "\nCannot proceed without gene list data.")
}

ISR <- read.csv(input_file)

message("Raw data dimensions: ", nrow(ISR), " rows x ", ncol(ISR), " columns")
message("Columns: ", paste(colnames(ISR), collapse = ", "))

# Convert to long format and clean
long_ISR <- ISR %>%
  pivot_longer(cols = everything(), names_to = "Source", values_to = "Gene") %>%
  filter(!is.na(Gene) & Gene != "") %>%
  distinct()

# Convert all genes to human format (uppercase)
# JacksonLabs genes are mouse format - need conversion
long_ISR <- long_ISR %>%
  mutate(Gene_Human = sapply(Gene, convert_mouse_to_human))

message("After conversion to human gene names: ", nrow(long_ISR), " gene entries")

# Create wide format with presence indicators
wide_ISR <- long_ISR %>%
  mutate(Presence = "Yes") %>%
  pivot_wider(
    id_cols = Gene_Human,
    names_from = Source,
    values_from = Presence,
    values_fill = list(Presence = "No")
  ) %>%
  rename(Gene = Gene_Human)

message("Unique genes after deduplication: ", nrow(wide_ISR))

# Extract gene lists for each source
genes_GO <- wide_ISR %>% filter(GO_consortium == "Yes") %>% pull(Gene)
genes_AnyGenes <- wide_ISR %>% filter(AnyGenes == "Yes") %>% pull(Gene)
genes_JacksonLabs <- wide_ISR %>% filter(JacksonLabs == "Yes") %>% pull(Gene)

message("\nGene counts by source:")
message("  GO Consortium: ", length(genes_GO))
message("  AnyGenes: ", length(genes_AnyGenes))
message("  JacksonLabs: ", length(genes_JacksonLabs))

# Create total gene list (union of all three)
total_genes <- unique(wide_ISR$Gene)
message("  Total unique genes: ", length(total_genes))

# Add GDF15 if not already present
if (!"GDF15" %in% total_genes) {
  total_genes <- c(total_genes, "GDF15")
  message("  Added GDF15 to total list")
}

# ============================================================================
# SAVE PROCESSED DATA
# ============================================================================
# Save wide format with all information
write.csv(wide_ISR, file.path(output_dir, "ISR_Gene_Lists_Processed.csv"), row.names = FALSE)
message("\nSaved: ISR_Gene_Lists_Processed.csv")

# Save individual gene lists as simple text files (one gene per line)
writeLines(genes_GO, file.path(output_dir, "genes_GO_consortium.txt"))
writeLines(genes_AnyGenes, file.path(output_dir, "genes_AnyGenes.txt"))
writeLines(genes_JacksonLabs, file.path(output_dir, "genes_JacksonLabs_human.txt"))
writeLines(total_genes, file.path(output_dir, "genes_Total_plus_GDF15.txt"))

message("Saved individual gene list files")

# Save as RData for faster loading in figure script
gene_lists <- list(
  GO_consortium = genes_GO,
  AnyGenes = genes_AnyGenes,
  JacksonLabs = genes_JacksonLabs,
  Total = total_genes,
  wide_data = wide_ISR
)
save(gene_lists, file = file.path(output_dir, "processed_gene_lists.RData"))
message("Saved: processed_gene_lists.RData")

message("\n============ DATA PROCESSING COMPLETE ============")
message("Run 02_Generate_Figure.R to create the Venn diagram")
