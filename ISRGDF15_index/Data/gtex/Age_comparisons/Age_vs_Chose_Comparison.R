
# Age_vs_Chose_Comparison.R
# ============================================================================
# Purpose: Read per-tissue CSV files, compute all 17 Spearman comparisons
#          (12 Factors + 4 single genes + 2 gene averages) vs AGE, write CSVs.
#
# Inputs:
#   - Per-tissue CSVs from tissue_expression_data/
#   - Fibroblast scaling parameters (mean, sd, loadings, correlations)
#
# Outputs:
#   - 12 Factor outputs in Comparing_AGE_Other_Factors_BH/.../Factor<N>_vs_AGE_SpearmanRhos.csv
#     Format: Tissue, ComparisonFactor<N>, Spearman_Rho_chosen_Factor, P_ValueFactor<N>,
#             AsterisksFactor<N>, Comparison_GDF15, Spearman_Rho_GDF15, P_Value_GDF15,
#             Asterisks_GDF15
#   - 4 single gene outputs (ATF4, ATF5, DDIT3, GDF15) in respective Comparing_AGE_Other_Factors_BH_<gene>
#     Format: Tissue, ComparisonFactor1, Spearman_Rho_chosen_Factor, P_ValueFactor1,
#             AsterisksFactor1, Comparison_<gene>, Spearman_Rho_<gene>, P_Value_<gene>,
#             Asterisks_<gene>
#   - 2 gene average outputs (proliferation, senescence) in 2024_12_02 folders
#     Format: Tissue, ComparisonFactor1, Spearman_Rho_chosen_Factor, P_ValueFactor1,
#             AsterisksFactor1, Comparison, Spearman_Rho, P_Value, Asterisks
#     (note: avg side uses unsuffixed Spearman_Rho/P_Value/Asterisks because
#      Figure_3D_Base.R reads "Spearman_Rho" directly from these files)
#
# Column naming contract matches what Figure_3D_Base.R reads:
#   - Spearman_Rho_chosen_Factor (Factor* files, single gene files, prolif/senes files)
#   - Spearman_Rho_<gene> (single gene files; also GDF15 in BH/Factor* files)
#   - Spearman_Rho (prolif/senes files)
# ============================================================================

library(here)
library(tidyverse)

source(here("Main Figure Scripts", "Helper Scripts", "gtex", "Attibutes_Phenos_Merged_plus_COD.R"))
# Loads Total_ISR_List

# ============================================================================
# Read fibroblast scaling parameters
# ============================================================================

processed <- here("Data", "Fibroblast_lifespan", "Processed")

mean_x <- read.csv(file.path(processed, "Fibroblast_FactorAnalysis_meanx_for_scaling.csv")) %>%
  rename(Gene = X, mean = x)

sd_x <- read.csv(file.path(processed, "Fibroblast_FactorAnalysis_sdx_for_scaling.csv")) %>%
  rename(Gene = X, sd = x)

correlations <- read.csv(file.path(processed, "Fibroblast_FactorAnalysis_correlations_for_scaling.csv"))

loadings_raw <- read.csv(file.path(processed, "Fibroblast_FactorAnalysis_FacLoads_12_factors.csv"))

# Prepare matrices for factor projection
rownames(loadings_raw) <- loadings_raw$X
loadings_matrix <- as.matrix(loadings_raw[, paste0("Factor", 1:12)])

rownames(correlations) <- correlations$X
correlations_matrix <- as.matrix(correlations[, -1])  # drop the X column

# ============================================================================
# Gene panels — proliferation and senescence are NOT in the fibroblast ISR
# factor analysis gene list, so they don't have fibroblast-derived scaling
# parameters. We z-score them within-tissue instead (matches what
# Figure_1E_1F.R does for these same genes via base R scale()).
# ============================================================================

prolif_genes <- c("TOP2A", "RRM2", "MKI67")
senes_genes  <- c("CDKN1A", "CDKN2A", "CCND2")
within_tissue_z_genes <- c(prolif_genes, senes_genes)

# ============================================================================
# Helper: asterisks from p-value
# ============================================================================

compute_asterisks <- function(p) {
  ifelse(is.na(p), "ns",
    ifelse(p < 0.0001, "****",
      ifelse(p < 0.001, "***",
        ifelse(p < 0.01, "**",
          ifelse(p < 0.05, "*", "ns")))))
}

# ============================================================================
# Initialize storage. Each entry is a long-format dataframe with generic
# columns (Tissue, Comparison, Spearman_Rho, P_Value, Asterisks) so we can
# reuse a single merge helper at write time.
# ============================================================================

empty_long_df <- function() {
  data.frame(Tissue = character(),
             Comparison = character(),
             Spearman_Rho = numeric(),
             P_Value = numeric(),
             Asterisks = character(),
             stringsAsFactors = FALSE)
}

factor_results <- setNames(
  lapply(1:12, function(k) empty_long_df()),
  paste0("Factor", 1:12)
)

gene_results <- list(
  ATF4  = empty_long_df(),
  ATF5  = empty_long_df(),
  DDIT3 = empty_long_df(),
  GDF15 = empty_long_df()
)

prolif_results <- empty_long_df()
senes_results  <- empty_long_df()

# ============================================================================
# Per-tissue main loop
# ============================================================================

tissue_files <- list.files(here("Data", "gtex", "Age_comparisons", "tissue_expression_data"),
                           pattern = "\\.csv$", full.names = TRUE)

if (length(tissue_files) == 0) {
  stop("No per-tissue CSVs found in Data/gtex/Age_comparisons/tissue_expression_data/. ",
       "Run Tissue_Expression_Processing.R first.")
}

for (tf in tissue_files) {
  tissue <- tools::file_path_sans_ext(basename(tf))
  message("Processing tissue: ", tissue)

  per_tissue <- read.csv(tf, check.names = FALSE)

  metadata_cols <- c("SAMPID", "AGE", "SEX", "DTHPLCE", "SMTSD")
  gene_cols <- setdiff(colnames(per_tissue), metadata_cols)

  # ----------------------------------------------------------------------
  # MISSING-GENE FILTER (mirrors original BH/ATF4/ATF5/DDIT3 scripts).
  # If any ISR gene from Total_ISR_List$Gene is absent from this tissue's
  # expression matrix (typically because edgeR::filterByExpr in Script 1
  # dropped it for low expression), skip this tissue entirely. This is what
  # reduced the original GTEx tissue count from 51 -> 44 (cell lines and
  # very-small-N tissues like cervix/bladder/fallopian_tube/kidney_medulla
  # tend to lose ISR genes via filterByExpr and get filtered out here).
  #
  # The skip applies to ALL comparisons for this tissue (Factor1-12, single
  # genes, prolif/senesc), matching the original per-script `next` behavior.
  # ----------------------------------------------------------------------
  isr_genes_clean   <- unique(trimws(toupper(Total_ISR_List$Gene)))
  exprs_genes_clean <- unique(trimws(toupper(gene_cols)))
  missing_isr_in_exprs <- setdiff(isr_genes_clean, exprs_genes_clean)
  extras_in_exprs      <- setdiff(exprs_genes_clean, isr_genes_clean)
  # `extras_in_exprs` will normally include the prolif/senesc genes (TOP2A,
  # RRM2, MKI67, CDKN1A, CDKN2A, CCND2) that we deliberately added in Script
  # 1 — those are EXPECTED extras, not a reason to skip. So we only check
  # the missing-from-exprs side, which is what indicates filterByExpr dropped
  # an ISR gene.
  if (length(missing_isr_in_exprs) > 0) {
    message("  Skipping ", tissue, ": ", length(missing_isr_in_exprs),
            " ISR gene(s) missing from expression matrix (",
            paste(head(missing_isr_in_exprs, 5), collapse = ", "),
            if (length(missing_isr_in_exprs) > 5)
              paste0(", + ", length(missing_isr_in_exprs) - 5, " more")
            else "",
            ")")
    next
  }

  # ----------------------------------------------------------------------
  # Apply fibroblast scaling per ISR gene (genes that are in mean_x/sd_x)
  # ----------------------------------------------------------------------
  per_tissue_scaled <- per_tissue
  for (g in gene_cols) {
    mu <- mean_x$mean[mean_x$Gene == g]
    sg <- sd_x$sd[sd_x$Gene == g]
    if (length(mu) == 1 && length(sg) == 1 && sg > 0) {
      per_tissue_scaled[[g]] <- (per_tissue[[g]] - mu) / sg
    } else {
      # Gene not in fibroblast scaling parameters — set to NA for now.
      # The within-tissue z-score block below will overwrite this for
      # the proliferation/senescence panel genes.
      per_tissue_scaled[[g]] <- NA_real_
    }
  }

  # ----------------------------------------------------------------------
  # Within-tissue z-score for proliferation + senescence genes
  # (these are NOT in the fibroblast ISR panel so they have no fibroblast
  # scaling parameters). Mirrors Figure_1E_1F.R's `scale()` of these genes
  # within the dataframe being analyzed.
  # ----------------------------------------------------------------------
  for (g in within_tissue_z_genes) {
    if (g %in% gene_cols) {
      raw_vals <- per_tissue[[g]]
      v <- var(raw_vals, na.rm = TRUE)
      if (!is.na(v) && v > 0) {
        per_tissue_scaled[[g]] <- as.numeric(scale(raw_vals))
      }
    }
  }

  # ----------------------------------------------------------------------
  # Factor projection (only ISR genes participate; prolif/senes are not in
  # loadings_matrix rownames so they're naturally excluded by intersect)
  # ----------------------------------------------------------------------
  common_genes <- intersect(rownames(loadings_matrix), gene_cols)
  common_genes <- intersect(common_genes, rownames(correlations_matrix))

  # Drop genes whose scaled column is fully NA
  common_genes <- common_genes[
    sapply(common_genes, function(g) !all(is.na(per_tissue_scaled[[g]])))
  ]

  factor_scores <- NULL
  if (length(common_genes) >= 2) {
    L_sub <- loadings_matrix[common_genes, , drop = FALSE]
    C_sub <- correlations_matrix[common_genes, common_genes, drop = FALSE]
    W  <- t(L_sub) %*% solve(C_sub)   # 12 x g
    Wt <- t(W)                        # g x 12
    X_sub <- as.matrix(per_tissue_scaled[, common_genes, drop = FALSE])
    factor_scores <- X_sub %*% Wt     # n x 12
    colnames(factor_scores) <- paste0("Factor", 1:12)
  } else {
    message("  Tissue ", tissue, ": only ", length(common_genes),
            " ISR genes available — skipping factor comparisons.")
  }

  # Helper: append a (Spearman vs AGE) result row to a long-format df
  append_spearman <- function(target_df, comparison_label, age_vec, score_vec) {
    df <- data.frame(AGE = age_vec, score = score_vec)
    df <- df[complete.cases(df), ]
    if (nrow(df) < 3) return(target_df)
    rho <- suppressWarnings(cor.test(df$AGE, df$score, method = "spearman"))
    adj_p <- p.adjust(rho$p.value, method = "BH")  # no-op on a single p; matches original script behavior
    rbind(target_df,
          data.frame(Tissue = tissue,
                     Comparison = comparison_label,
                     Spearman_Rho = unname(rho$estimate),
                     P_Value = adj_p,
                     Asterisks = compute_asterisks(adj_p),
                     stringsAsFactors = FALSE))
  }

  # ----- Factor1..Factor12 vs AGE -----
  if (!is.null(factor_scores)) {
    for (k in 1:12) {
      key <- paste0("Factor", k)
      factor_results[[key]] <- append_spearman(
        factor_results[[key]],
        comparison_label = "chosen_Factor",
        age_vec = per_tissue_scaled$AGE,
        score_vec = factor_scores[, k]
      )
    }
  }

  # ----- Each single gene vs AGE -----
  for (g in c("ATF4", "ATF5", "DDIT3", "GDF15")) {
    if (!g %in% gene_cols) next
    if (all(is.na(per_tissue_scaled[[g]]))) next
    gene_results[[g]] <- append_spearman(
      gene_results[[g]],
      comparison_label = g,
      age_vec = per_tissue_scaled$AGE,
      score_vec = per_tissue_scaled[[g]]
    )
  }

  # ----- Proliferation avg vs AGE -----
  if (all(prolif_genes %in% gene_cols)) {
    prolif_mat <- as.matrix(per_tissue_scaled[, prolif_genes, drop = FALSE])
    if (!all(is.na(prolif_mat))) {
      prolif_results <- append_spearman(
        prolif_results,
        comparison_label = "Proliferation_avg",
        age_vec = per_tissue_scaled$AGE,
        score_vec = rowMeans(prolif_mat, na.rm = TRUE)
      )
    }
  }

  # ----- Senescence avg vs AGE -----
  if (all(senes_genes %in% gene_cols)) {
    senes_mat <- as.matrix(per_tissue_scaled[, senes_genes, drop = FALSE])
    if (!all(is.na(senes_mat))) {
      senes_results <- append_spearman(
        senes_results,
        comparison_label = "Senescence_avg",
        age_vec = per_tissue_scaled$AGE,
        score_vec = rowMeans(senes_mat, na.rm = TRUE)
      )
    }
  }
}

# ============================================================================
# Write outputs in dual-comparison format
# ============================================================================
#
# The historical output format (from the original gene-vs-AGE scripts) was:
#   Tissue, Comparison<left_suffix>, Spearman_Rho_chosen_Factor,
#   P_Value<left_suffix>, Asterisks<left_suffix>,
#   Comparison<right_suffix>, Spearman_Rho<right_suffix>,
#   P_Value<right_suffix>, Asterisks<right_suffix>
#
# This is produced by merging two long-format dataframes by Tissue with
# `suffixes = c(<left>, <right>)`, then manually renaming the Spearman_Rho
# column on the left side to "Spearman_Rho_chosen_Factor".
# ============================================================================

merge_dual <- function(left_df, left_suffix, right_df, right_suffix) {
  combined <- merge(left_df, right_df, by = "Tissue",
                    suffixes = c(left_suffix, right_suffix), all = TRUE)
  # Rename Spearman_Rho<left_suffix> -> Spearman_Rho_chosen_Factor (matches
  # the original script's manual rename of column 3).
  spearman_left <- paste0("Spearman_Rho", left_suffix)
  if (spearman_left %in% colnames(combined)) {
    colnames(combined)[colnames(combined) == spearman_left] <- "Spearman_Rho_chosen_Factor"
  }
  combined
}

message("\nWriting output CSVs...")

# ---- BH/Factor1..Factor12 (Factor + GDF15 dual format) ----
for (k in 1:12) {
  out_dir <- here("Data", "gtex", "Age_comparisons",
                  "Comparing_AGE_Other_Factors_BH",
                  paste0("2024_10_28_All_Factors_AGE_Factor", k))
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  combined <- merge_dual(
    left_df = factor_results[[paste0("Factor", k)]],
    left_suffix = paste0("Factor", k),
    right_df = gene_results$GDF15,
    right_suffix = "_GDF15"
  )

  out_file <- file.path(out_dir, paste0("Factor", k, "_vs_AGE_SpearmanRhos.csv"))
  write.csv(combined, out_file)  # row.names = TRUE by default; produces leading row index column
  message("  Wrote: ", out_file)
}

# ---- Single gene files (Factor1 + Gene dual format) ----
gene_paths <- list(
  ATF4  = c("Comparing_AGE_Other_Factors_BH_ATF4",  "2024_10_30_All_Factors_AGE_ATF4"),
  ATF5  = c("Comparing_AGE_Other_Factors_BH_ATF5",  "2024_10_30_All_Factors_AGE_ATF5"),
  DDIT3 = c("Comparing_AGE_Other_Factors_BH_DDIT3", "2024_10_30_All_Factors_AGE_DDIT3"),
  GDF15 = c("Comparing_AGE_Other_Factors_BH_GDF15", "2024_10_30_All_Factors_AGE_GDF15")  # NEW standalone folder
)

for (g in names(gene_paths)) {
  out_dir <- here("Data", "gtex", "Age_comparisons", gene_paths[[g]][1], gene_paths[[g]][2])
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  combined <- merge_dual(
    left_df = factor_results$Factor1,
    left_suffix = "Factor1",
    right_df = gene_results[[g]],
    right_suffix = paste0("_", g)
  )

  out_file <- file.path(out_dir, "Factor1_vs_AGE_SpearmanRhos.csv")
  write.csv(combined, out_file)
  message("  Wrote: ", out_file)
}

# ---- Proliferation and Senescence (Factor1 + avg dual format; avg side uses
# UNSUFFIXED Spearman_Rho/P_Value/Asterisks columns because Figure_3D_Base.R
# reads "Spearman_Rho" directly from these files) ----
for (panel in list(
  list(results = prolif_results,
       out = c("2024_12_02_All_Factors_AGE_AvgOf3ProliferationGenes")),
  list(results = senes_results,
       out = c("2024_12_02_All_Factors_AGE_AvgOf3SenescentGenes"))
)) {
  out_dir <- here("Data", "gtex", "Age_comparisons", panel$out[1])
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  combined <- merge_dual(
    left_df = factor_results$Factor1,
    left_suffix = "Factor1",
    right_df = panel$results,
    right_suffix = ""           # empty suffix → avg-side columns stay as Spearman_Rho/P_Value/Asterisks
  )

  out_file <- file.path(out_dir, "Factor1_vs_AGE_SpearmanRhos.csv")
  write.csv(combined, out_file)
  message("  Wrote: ", out_file)
}

message("\nAll 17 output CSVs written.")
