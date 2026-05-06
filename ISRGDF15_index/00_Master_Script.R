# ============================================================================
# MASTER SCRIPT: Generate All Main Figures for ISRGDF15 Paper
# ============================================================================
#
# This script orchestrates the full pipeline to reproduce all main figures:
#   Step 1: Process raw data (Data Processing Scripts)
#   Step 2: Generate figures (Main Figure Scripts)
#   Step 3: Compile all source data CSVs into a single master Excel file
#
# HOW TO RUN:
#   1. Open ISRGDF15_Paper_Figures.Rproj in RStudio (this anchors all here() paths)
#   2. Run this script (source it or run line by line)
#
# OUTPUT:
#   - Figures saved to:       Results/Figures/<Figure_Name>/
#   - Individual source CSVs: saved alongside figures (Results/Figures/<Figure_Name>/)
#   - Master Excel:           Results/CSVs/All_Figure_Source_Data.xlsx
#                             (one tab per figure with all numerical source data)
#
# DATA REQUIREMENTS:
#   All required data files must be present in the Data/ folder.
#   See README.md for the complete list of required data files.
#
# RUNTIME NOTE:
#   - Figures 3B/4A and 3C (GTEx tissue scatter plots) are slow — they iterate
#     over 49 tissues and compute Spearman correlations.
#   - Figure 3D also requires loading GTEx TPM data (large file).
#   - Total expected runtime: 0.5–4 hours depending on system.
# ============================================================================

library(here)

# ============================================================================
# CONFIGURATION
# ============================================================================
# Set to TRUE to run the full pipeline (including slow GTEx steps)
# Set to FALSE to skip data processing (if processed data already exists)
RUN_DATA_PROCESSING <- TRUE

# Set to TRUE to run all figure scripts; FALSE to run individual ones below
RUN_ALL_FIGURES <- TRUE

# ============================================================================
# HELPER FUNCTION: Create output directories
# ============================================================================
create_output_dirs <- function() {
  dirs <- c(
    # Figure output folders (figures + per-figure source data CSVs saved here)
    here("Results", "Figures", "Figure_1A"),
    here("Results", "Figures", "Figure_1D"),
    here("Results", "Figures", "Figure_1E"),
    here("Results", "Figures", "Figure_1F"),
    here("Results", "Figures", "Figure_1G"),
    here("Results", "Figures", "Figure_2C"),
    here("Results", "Figures", "Figure_2E"),
    here("Results", "Figures", "Figure_2F"),
    here("Results", "Figures", "Figure_2G"),
    here("Results", "Figures", "Figure_3B_4A"),
    here("Results", "Figures", "Figure_3B_Summary"),
    here("Results", "Figures", "Figure_3C"),
    here("Results", "Figures", "Figure_3D"),
    here("Results", "Figures", "Figure_4B"),
    here("Results", "Figures", "Figure_4C"),
    here("Results", "Figures", "Figure_4D_4E"),
    # Data processing output folders
    here("Data", "ISR_Gene_Lists"),
    here("Data", "Fibroblast_lifespan", "Processed")
  )
  for (d in dirs) {
    if (!dir.exists(d)) dir.create(d, recursive = TRUE)
  }
  message("Output directories created.")
}

create_output_dirs()

# ============================================================================
# STEP 1: DATA PROCESSING
# ============================================================================
if (RUN_DATA_PROCESSING) {
  message("\n", strrep("=", 70))
  message("STEP 1: DATA PROCESSING")
  message(strrep("=", 70))

  # Process ISR gene lists (required for Figure 1A)
  message("\n[1/4] Processing ISR gene lists (Figure 1A)...")
  local(source(here("Data Processing Scripts", "Figure_1A_Process_Data.R"), local = TRUE))

  # Process fibroblast ISR boxplot data (required for Figure 1G)
  message("\n[2/4] Processing fibroblast ISR boxplot data (Figure 1G)...")
  local(source(here("Data Processing Scripts", "Figure_1G_Process_Data.R"), local = TRUE))

  # Process GTEx tissue proliferation index (required for Figures 4D, 4E)
  # NOTE: This step reads the ~2.5 GB GTEx TPM file and requires significant memory.
  # It is skipped if the cached output files already exist. To force regeneration,
  # delete the files in Data/gtex/Processed/ and run Tissue_proliferation_index.R
  # in a fresh R session with sufficient memory (R_MAX_VSIZE=32Gb in .Renviron).
  prolif_index_file <- here("Data", "gtex", "Processed", "proliferation_index_Jack_Devine_namesChanged.csv")
  prolif_donor_file <- here("Data", "gtex", "Processed", "proliferation_per_donor.csv")
  if (file.exists(prolif_index_file) && file.exists(prolif_donor_file)) {
    message("\n[3/4] Tissue proliferation index: cached files found, skipping.")
    message("       To regenerate, delete files in Data/gtex/Processed/ and run")
    message("       Tissue_proliferation_index.R in a fresh R session.")
  } else {
    message("\n[3/4] Processing tissue proliferation index (Figures 4D, 4E)...")
    message("       WARNING: This reads the large GTEx TPM file and needs ~20+ GB RAM.")
    local(source(here("Data Processing Scripts", "Tissue_proliferation_index.R"), local = TRUE))
  }

  # Process per-tissue GTEx expression for Age comparisons (required for
  # Age_vs_Chose_Comparison.R, which writes the CSVs that Figure_3D reads).
  # NOTE: This step TMM-normalizes every GTEx tissue and writes one CSV per
  # tissue to Data/gtex/Age_comparisons/tissue_expression_data/. It is the
  # slowest data-processing step (~30+ min). It is skipped if that folder
  # already exists and contains CSVs. To force regeneration, delete the
  # tissue_expression_data/ folder and rerun.
  tissue_expr_dir <- here("Data", "gtex", "Age_comparisons", "tissue_expression_data")
  tissue_expr_existing <- if (dir.exists(tissue_expr_dir)) {
    list.files(tissue_expr_dir, pattern = "\\.csv$", full.names = FALSE)
  } else {
    character(0)
  }
  if (length(tissue_expr_existing) > 0) {
    message("\n[4/4] Per-tissue GTEx expression: cached files found (",
            length(tissue_expr_existing), " tissues), skipping.")
    message("       To regenerate, delete Data/gtex/Age_comparisons/tissue_expression_data/ and rerun.")
  } else {
    message("\n[4/4] Processing per-tissue GTEx expression for Age comparisons...")
    message("       This is the slow step (TMM-normalizes every GTEx tissue, ~30+ min).")
    local(source(here("Data", "gtex", "Age_comparisons", "Tissue_Expression_Processing.R"), local = TRUE))
  }

  message("\nData processing complete.")
} else {
  message("\nSkipping data processing (RUN_DATA_PROCESSING = FALSE)")
  message("Assuming processed data already exists in Data/Fibroblast_lifespan/Processed/")
}

# ============================================================================
# STEP 2: GENERATE FIGURES
# ============================================================================
if (RUN_ALL_FIGURES) {
  message("\n", strrep("=", 70))
  message("STEP 2: GENERATING FIGURES")
  message(strrep("=", 70))

  # --- FIGURE 1A: ISR Gene Venn Diagram ---
  message("\n[Figure 1A] ISR Gene Venn Diagram...")
  local(source(here("Main Figure Scripts", "Figure_1A.R"), local = TRUE))

  # --- FIGURE 1D: PCA Biplot ---
  message("\n[Figure 1D] PCA Biplot (requires Fibroblast_lifespan data)...")
  local(source(here("Main Figure Scripts", "Figure_1D.R"), local = TRUE))

  # --- FIGURES 1E/1F: Factor Analysis Plots + Source Data CSVs ---
  message("\n[Figures 1E/1F] Factor Analysis (heatmap, bar chart, GDF15 scatter)...")
  local(source(here("Main Figure Scripts", "Figure_1E_1F.R"), local = TRUE))

  # --- AGE COMPARISONS: per-tissue Spearman vs AGE for Factor1-12 + 4 single
  # genes + 2 gene-average panels. Writes 17 CSVs into
  # Data/gtex/Age_comparisons/<various subfolders>/Factor*_vs_AGE_SpearmanRhos.csv
  # which Figure_3D_Base.R (sourced by Figure_3D.R) consumes downstream.
  #
  # Depends on:
  #   - Step 1 [4/4]: per-tissue CSVs in tissue_expression_data/
  #   - Figure_1E_1F.R above: Fibroblast_FactorAnalysis_*.csv scaling/loadings
  # Hence this slot — after Figure 1E/1F, before Figure 3D.
  message("\n[Age Comparisons] Computing GTEx Spearman vs AGE per tissue...")
  local(source(here("Data", "gtex", "Age_comparisons", "Age_vs_Chose_Comparison.R"), local = TRUE))

  # --- FIGURE 1G: ISR Boxplots ---
  message("\n[Figure 1G] ISR Boxplots...")
  local(source(here("Main Figure Scripts", "Figure_1G.R"), local = TRUE))

  # --- FIGURE 2C: Opto-PKR Kinetics ---
  message("\n[Figure 2C] Opto-PKR Kinetics...")
  local(source(here("Main Figure Scripts", "Figure_2C.R"), local = TRUE))

  # --- FIGURE 2E: DmrPERK Dose Response ---
  message("\n[Figure 2E] DmrPERK Dose Response (journal-compliant)...")
  local(source(here("Main Figure Scripts", "Figure_2E.R"), local = TRUE))

  # --- FIGURE 2F: ISR-Activating Drugs ---
  message("\n[Figure 2F] ISR-Activating Drugs (journal-compliant)...")
  local(source(here("Main Figure Scripts", "Figure_2F.R"), local = TRUE))

  # --- FIGURE 2G: DmrPERK vs ATF4KO ---
  message("\n[Figure 2G] DmrPERK vs ATF4KO (journal-compliant)...")
  local(source(here("Main Figure Scripts", "Figure_2G.R"), local = TRUE))

  # --- FIGURE 3B/4A: GTEx Tissue Scatter Plots ---
  message("\n[Figures 3B/4A] GTEx Tissue Scatter Plots (SLOW - ~1 hour)...")
  message("  This generates individual Spearman correlation scatter plots for all tissues.")
  local(source(here("Main Figure Scripts", "Figure_3B_4A_Tissue_Scatter.R"), local = TRUE))

  # --- FIGURE 3B SUMMARY: GTEx Summary Plot ---
  message("\n[Figure 3B Summary] GTEx Summary Overview Plot (journal-compliant)...")
  local(source(here("Main Figure Scripts", "Figure_3B_Summary.R"), local = TRUE))

  # --- FIGURE 3C: GTEx Age vs ISR Scatter Plots ---
  message("\n[Figure 3C] GTEx Age vs ISR Scatter Plots (journal-compliant)...")
  local(source(here("Main Figure Scripts", "Figure_3C.R"), local = TRUE))

  # --- FIGURE 3D: GTEx vs Fibroblast Comparison ---
  message("\n[Figure 3D] GTEx vs Fibroblast Comparison (journal-compliant, loads large TPM file)...")
  local(source(here("Main Figure Scripts", "Figure_3D.R"), local = TRUE))

  # --- FIGURE 4B: Place of Death ---
  message("\n[Figure 4B] Place of Death...")
  local(source(here("Main Figure Scripts", "Figure_4B.R"), local = TRUE))

  # --- FIGURE 4C: Hedges' G Effect Sizes ---
  message("\n[Figure 4C] Hedges' G Effect Sizes...")
  local(source(here("Main Figure Scripts", "Figure_4C.R"), local = TRUE))

  # --- FIGURES 4D/4E: Proliferation vs ISR ---
  message("\n[Figures 4D/4E] Proliferation vs ISR (journal-compliant)...")
  local(source(here("Main Figure Scripts", "Figure_4D_4E.R"), local = TRUE))

  message("\n", strrep("=", 70))
  message("ALL FIGURES GENERATED SUCCESSFULLY")
  message("Results saved to: ", here("Results"))
  message(strrep("=", 70))
}



message("\n", strrep("=", 70))
message("PIPELINE COMPLETE")
message("  Figures:       ", here("Results", "Figures"))
message("  Master Excel:  ", here("Results", "CSVs", "All_Figure_Source_Data.xlsx"))
message(strrep("=", 70))
