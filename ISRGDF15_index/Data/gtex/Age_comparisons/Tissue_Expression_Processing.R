
# Tissue_Expression_Processing.R
# ============================================================================
# Purpose: For each GTEx tissue, perform TMM normalization once and save
#          raw log2 TMM expression for genes needed by any downstream comparison.
#
# Inputs:
#   - GTEx .gz expression files in Data/gtex/All_Tissues_Indiv_Folders
#   - Sample annotations and Total_ISR_List from helper script
#
# Outputs:
#   - Per-tissue CSV files in Data/gtex/Age_comparisons/tissue_expression_data/
#   - Columns: SAMPID, AGE, SEX, DTHPLCE, SMTSD, plus all genes needed
# ============================================================================

library(here)
library(tidyverse)
library(edgeR)

# Source the helper to get annotations_merged and Total_ISR_List
source(here("Main Figure Scripts", "Helper Scripts", "gtex", "Attibutes_Phenos_Merged_plus_COD.R"))

# Define the union of genes any comparison will need:
# - Total_ISR_List: used for all Factor1-12 projections (includes ATF4, ATF5, DDIT3, GDF15)
# - proliferation_genes: TOP2A, RRM2, MKI67
# - senescence_genes: CDKN1A, CDKN2A, CCND2
genes_needed <- unique(c(
  Total_ISR_List$Gene,
  c("TOP2A", "RRM2", "MKI67"),         # proliferation
  c("CDKN1A", "CDKN2A", "CCND2")        # senescence
))

# Output directory
out_dir <- here("Data", "gtex", "Age_comparisons", "tissue_expression_data")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Loop over tissues
individual_tissues_folder <- here("Data", "gtex", "All_Tissues_Indiv_Folders")
subfolders <- list.dirs(individual_tissues_folder, recursive = FALSE, full.names = TRUE)

for (folder in subfolders) {
  tissue_name <- basename(folder)
  message("Processing: ", tissue_name)

  # Get a list of all .gz files in the tissue folder
  gz_files <- list.files(folder, pattern = "\\.gz$", full.names = TRUE)

  # Initialize dataframes for each tissue
  df_tissue_list <- list()
  gene_symbols <- NULL

  # Loop through each .gz file and read data
  for (gz_file in gz_files) {
    file_name <- basename(gz_file)
    tissue <- sub(".*_v8_(.*?)\\.gct\\.gz", "\\1", file_name)

    # Read data from the current .gz file
    raw_gene_counts <- read.delim(gz_file, skip = 2)

    # Extract gene symbols if not already done
    if (is.null(gene_symbols)) {
      gene_symbols <- raw_gene_counts[, 2:3]
    }

    # Extract sample IDs and gene counts
    gene_counts <- raw_gene_counts %>%
      select(-id, -Description) %>%
      as.data.frame()

    # Store the gene counts in the dataframe list with tissue name
    df_tissue_list[[tissue]] <- gene_counts
  }

  # Combine dataframes if there are multiple files
  combined_gene_counts <- gene_symbols
  for (tissue in names(df_tissue_list)) {
    combined_gene_counts <- combined_gene_counts %>%
      left_join(df_tissue_list[[tissue]], by = "Name")
  }

  # Remove the "Description" column if it exists
  if ("Description" %in% colnames(combined_gene_counts)) {
    combined_gene_counts <- combined_gene_counts %>%
      select(-Description)
  }

  combined_gene_counts <- combined_gene_counts %>%
    column_to_rownames("Name")

  # Filter samples to those in annotations_merged
  annotations_sample_attributes_filtered <- annotations_merged %>%
    filter(SAMPID %in% colnames(combined_gene_counts))

  # TMM normalization
  dge <- DGEList(combined_gene_counts[, annotations_sample_attributes_filtered$SAMPID],
                 group = factor(annotations_sample_attributes_filtered$COD))
  keep <- filterByExpr(dge)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge, method = "TMM")
  exprs <- cpm(dge, log = TRUE)

  # Apply gene symbol mappings
  gene_symbols <- gene_symbols %>%
    mutate(Description = ifelse(Description == "KIAA0141", "DELE1", Description),
           Description = ifelse(Description == "WARS", "WARS1", Description),
           Description = ifelse(Description == "NARS", "NARS1", Description))

  # Join with gene symbols
  exprs <- exprs %>%
    as.data.frame() %>%
    rownames_to_column("Name") %>%
    full_join(gene_symbols, by = "Name") %>%
    na.omit() %>%
    select(-Name) %>%
    rename("Gene" = Description)

  # Filter to genes_needed (only those present in this tissue)
  # Note: NOT restricting to Total_ISR_List, but to genes_needed union
  exprs_filtered <- exprs %>%
    filter(Gene %in% genes_needed) %>%
    column_to_rownames("Gene") %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("SAMPID") %>%
    unique()

  # Join with metadata (AGE, SEX, DTHPLCE, SMTSD)
  # IMPORTANT: use inner_join + filter(AGE not NA) instead of full_join + na.omit().
  # The original single-gene scripts used full_join + na.omit, which was fine when
  # the per-tissue dataframe had ~1 gene column. Here we have ~100 gene columns,
  # so na.omit() would drop any sample missing ANY single ISR/proliferation/senescence
  # gene — which would discard most samples. Script 2 handles per-comparison NA filtering
  # correctly, so here we only require SAMPID match + non-NA AGE (needed for Spearman).
  per_tissue_df <- exprs_filtered %>%
    inner_join(
      annotations_sample_attributes_filtered %>%
        select(SAMPID, AGE, SEX, DTHPLCE, SMTSD),
      by = "SAMPID"
    ) %>%
    filter(!is.na(AGE))

  # Reorder columns: metadata first, then genes (for human readability)
  metadata_cols <- intersect(c("SAMPID", "AGE", "SEX", "DTHPLCE", "SMTSD"), colnames(per_tissue_df))
  gene_cols <- setdiff(colnames(per_tissue_df), metadata_cols)
  per_tissue_df <- per_tissue_df[, c(metadata_cols, gene_cols)]

  # Write CSV
  out_path <- file.path(out_dir, paste0(tissue_name, ".csv"))
  write.csv(per_tissue_df, out_path, row.names = FALSE)

  message("  Wrote: ", out_path, " (", nrow(per_tissue_df), " samples)")
}

message("Done. Per-tissue CSVs written to: ", out_dir)
