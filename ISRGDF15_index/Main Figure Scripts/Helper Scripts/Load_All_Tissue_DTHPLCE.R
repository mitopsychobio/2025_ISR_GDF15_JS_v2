# ============================================================================
# Load_All_Tissue_DTHPLCE.R
# ============================================================================
# Shared loader for Figures 4B, 4C, and 4D/4E.
#
# Replaces the legacy single-file All_Tissue_data_DTHPLCE.csv (which no
# script in this repo actually generates). Instead, this helper assembles
# the same per-donor / per-tissue / per-metric long-format dataframe by
# reading the per-tissue Plot_Data_<tissue>.csv files that
# Figure_3B_4A_Tissue_Scatter.R writes to:
#   Results/Figures/Figure_3B_4A/Plot_Data_CSVs/
#
# Each per-tissue CSV has columns:
#   SAMPID, Tissue, DTHPLCE, Factor1, GDF15,
#   rho, p_value, p_adjust, asterisks, N_total_SAMPID, N_used_complete
#
# This helper pivots Factor1 and GDF15 from wide to long, producing rows
# with a `fa_vs_gdf15` column ("Factor1" or "GDF15") and a `value` column.
# It also adds a `group` column = paste(DTHPLCE, fa_vs_gdf15, sep = "_"),
# which is the form Figure_4C expects (e.g. "Hospital inpatient_Factor1").
#
# USAGE:
#   source(here("Main Figure Scripts", "Helper Scripts",
#               "Load_All_Tissue_DTHPLCE.R"))
#   data <- load_all_tissue_dthplce()
#
# RETURNS a data.frame with columns:
#   SAMPID, Tissue, DTHPLCE, fa_vs_gdf15, value, group,
#   rho, p_value, p_adjust, asterisks, N_total_SAMPID, N_used_complete
#
# NO DTHPLCE FILTERING is applied here. Rows include every DTHPLCE category
# present in the per-tissue files ("Hospital inpatient", "Emergency room",
# "Dead on arrival at hospital", "At MVA scene", NAs, etc.). Each consuming
# figure script applies its own filter to the H/E pair as needed.
# ============================================================================

library(here)
library(dplyr)
library(tidyr)

load_all_tissue_dthplce <- function(plot_data_dir = NULL) {

  if (is.null(plot_data_dir)) {
    plot_data_dir <- here("Results", "Figures", "Figure_3B_4A",
                          "Plot_Data_CSVs")
  }

  if (!dir.exists(plot_data_dir)) {
    stop(
      "\n========== ERROR: PER-TISSUE DATA NOT FOUND ==========\n",
      "Could not find Plot_Data_CSVs directory at:\n",
      "  ", plot_data_dir, "\n",
      "These per-tissue CSVs are produced by\n",
      "Figure_3B_4A_Tissue_Scatter.R. Run that script first.\n",
      "======================================================\n"
    )
  }

  csv_files <- list.files(plot_data_dir,
                          pattern    = "^Plot_Data_.*\\.csv$",
                          full.names = TRUE)

  # Defensive: drop any non-tissue files that happen to match the pattern
  csv_files <- csv_files[!grepl("All_Tissues_Sample_sizes",
                                basename(csv_files))]

  if (length(csv_files) == 0) {
    stop("No Plot_Data_*.csv files found in: ", plot_data_dir)
  }

  message("Loading ", length(csv_files),
          " per-tissue CSVs from: ", plot_data_dir)

  per_tissue <- lapply(csv_files, function(f) {
    df <- read.csv(f, stringsAsFactors = FALSE)
    if (!all(c("Factor1", "GDF15") %in% colnames(df))) {
      message("  WARNING: skipping ", basename(f),
              " (missing Factor1 or GDF15 column)")
      return(NULL)
    }
    df %>%
      tidyr::pivot_longer(
        cols      = c("Factor1", "GDF15"),
        names_to  = "fa_vs_gdf15",
        values_to = "value"
      ) %>%
      dplyr::mutate(
        group = paste(DTHPLCE, fa_vs_gdf15, sep = "_")
      )
  })

  combined <- dplyr::bind_rows(per_tissue)

  message("  Combined rows: ", nrow(combined),
          " (", dplyr::n_distinct(combined$Tissue), " tissues, ",
          dplyr::n_distinct(combined$SAMPID), " unique samples)")

  return(combined)
}
