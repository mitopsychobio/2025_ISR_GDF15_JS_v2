# ============================================================================
# Figure 3D: Journal-Compliant Plots with Individual Data Points
# Journal requirement: show individual data points behind summary values
# ============================================================================
# GTEx version:
#   Summary point = mean Spearman rho across ~54 tissues ± SEM
#   Individual points = per-tissue Spearman rho (one dot per tissue)
#
# Fibroblast version:
#   Summary point = overall Spearman rho (all samples pooled) ± CI
#   Individual points = per-condition (intx) Spearman rho for ALL metrics:
#     Factor1–12  (from FactorX_ISR_Scores_per_Sample_Age_DaysGrown.csv)
#     GDF15, ATF4, ATF5, DDIT3  (from raw expression data)
#     sen_mean (mean of CDKN1A, CDKN2A, CCND2 per sample)
#     prolif_mean (mean of MKI67, RRM2, TOP2A per sample)
#
# Outputs: Figures/Graphs Adhering to Journal Guidelines/Figure_3D/
# ============================================================================

library(here)
library(tidyverse)
library(ggplot2)
library(purrr)
library(stringr)



# ============================================================================
# SOURCE DATA-PREP SCRIPT
# Figure_3D_Base.R is now strictly a data-prep step: it builds All_Spearmans
# (per-tissue Spearman rhos for every metric) and writes the source-data and
# per-tissue CSVs. All plotting is done below in this script.
# ============================================================================
message("\n=== Sourcing Figure 3D data-prep script ===")

suppressMessages(
  source(here("Main Figure Scripts", "Helper Scripts", "Figure_3D_Base.R"), local = TRUE)
)

# Re-define output dir AFTER sourcing (source clears the environment)
out_dir <- here("Results", "Figures", "Figure_3D")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ============================================================================
# SHARED COLOUR + THEME DEFINITIONS (match original exactly)
# ============================================================================
color_mapping <- c(
  "Spearman_Rho_GDF15"           = "pink",
  "Factor1_vs_AGE_SpearmanRhos"  = "red",
  "Factor2_vs_AGE_SpearmanRhos"  = "#808080",
  "Factor3_vs_AGE_SpearmanRhos"  = "#808080",
  "Factor4_vs_AGE_SpearmanRhos"  = "#CCCCCC",
  "Factor5_vs_AGE_SpearmanRhos"  = "#CCCCCC",
  "Factor6_vs_AGE_SpearmanRhos"  = "#CCCCCC",
  "Factor7_vs_AGE_SpearmanRhos"  = "#808080",
  "Factor8_vs_AGE_SpearmanRhos"  = "#CCCCCC",
  "Factor9_vs_AGE_SpearmanRhos"  = "#808080",
  "Factor10_vs_AGE_SpearmanRhos" = "#CCCCCC",
  "Factor11_vs_AGE_SpearmanRhos" = "#808080",
  "Factor12_vs_AGE_SpearmanRhos" = "#808080",
  "Spearman_Rho_ATF4"            = "orange",
  "Spearman_Rho_DDIT3"           = "lightblue",
  "Spearman_Rho_ATF5"            = "brown",
  "Spearman_Rho_Senescence"      = "blue",
  "Spearman_Rho_Proliferation"   = "purple"
)

fb_color_mapping <- c(
  "GDF15"      = "pink",
  "Factor1"    = "red",
  "Factor2"    = "#808080",
  "Factor3"    = "#808080",
  "Factor4"    = "#CCCCCC",
  "Factor5"    = "#CCCCCC",
  "Factor6"    = "#CCCCCC",
  "Factor7"    = "#808080",
  "Factor8"    = "#CCCCCC",
  "Factor9"    = "#808080",
  "Factor10"   = "#CCCCCC",
  "Factor11"   = "#808080",
  "Factor12"   = "#808080",
  "ATF4"       = "orange",
  "DDIT3"      = "lightblue",
  "ATF5"       = "brown",
  "sen_mean"   = "blue",
  "prolif_mean"= "purple"
)

fig3d_theme <- theme(
  plot.title         = element_text(size = 16, face = "bold", hjust = 0.5),
  axis.title         = element_text(size = 14, color = "black"),
  axis.text          = element_text(size = 12, color = "black"),
  axis.ticks         = element_line(linewidth = 0.25, color = "black"),
  panel.grid.major.y = element_blank(),
  panel.grid.minor.y = element_blank(),
  panel.grid.major.x = element_line(linewidth = 0.1, color = "black"),
  panel.grid.minor.x = element_line(linewidth = 0.08, color = "black"),
  panel.background   = element_blank(),
  panel.border       = element_blank(),
  legend.position    = "none"
)

# PLOT ELEMENTS
chosen_width <- 6
chosen_alpha <- .8

# ============================================================================
# GTEX — INDIVIDUAL PER-TISSUE SPEARMAN RHOS
# ============================================================================
message("\n=== GTEx: building per-tissue individual points ===")

plot_columns <- c(
  "Spearman_Rho_GDF15", "Spearman_Rho_Senescence", "Spearman_Rho_Proliferation",
  "Spearman_Rho_ATF4", "Spearman_Rho_ATF5", "Spearman_Rho_DDIT3",
  "Factor1_vs_AGE_SpearmanRhos", "Factor12_vs_AGE_SpearmanRhos",
  "Factor6_vs_AGE_SpearmanRhos",  "Factor2_vs_AGE_SpearmanRhos",
  "Factor4_vs_AGE_SpearmanRhos",  "Factor9_vs_AGE_SpearmanRhos",
  "Factor8_vs_AGE_SpearmanRhos",  "Factor3_vs_AGE_SpearmanRhos",
  "Factor10_vs_AGE_SpearmanRhos", "Factor7_vs_AGE_SpearmanRhos",
  "Factor5_vs_AGE_SpearmanRhos",  "Factor11_vs_AGE_SpearmanRhos"
)
factor_levels_gtex <- rev(plot_columns)

# Per-tissue long format
tissue_long <- All_Spearmans %>%
  select(Tissue, all_of(plot_columns)) %>%
  pivot_longer(cols = all_of(plot_columns),
               names_to  = "Factor",
               values_to = "Spearman_Rho_tissue") %>%
  filter(!is.na(Spearman_Rho_tissue)) %>%
  mutate(Factor = factor(Factor, levels = factor_levels_gtex))

# Summary (mean ± 95% CI)
means <- sapply(All_Spearmans[plot_columns], mean, na.rm = TRUE)
sems  <- sapply(All_Spearmans[plot_columns],
                function(x) sd(x, na.rm = TRUE) / sqrt(length(na.omit(x))))
plot_data_gtex <- data.frame(Factor = names(means), Mean = means, SEM = sems) %>%
  mutate(CI_Lower = Mean - 1.96 * SEM,
         CI_Upper = Mean + 1.96 * SEM,
         Factor = factor(Factor, levels = factor_levels_gtex))

message("  ", nrow(tissue_long), " individual data points (",
        n_distinct(tissue_long$Tissue), " tissues)")

# ---- GTEx plot: per-tissue squares with mean circle (black-outlined) ----
p_gtex <- ggplot() +
  geom_point(
    data  = tissue_long,
    aes(x = Spearman_Rho_tissue, y = Factor, fill = Factor),
    shape = 22, color = NA, size = 2.5, alpha = 0.25,
    position = position_jitter(height = 0.15, width = 0, seed = 42),
    show.legend = FALSE
  ) +
  geom_point(
    data  = tissue_long,
    aes(x = Spearman_Rho_tissue, y = Factor, color = Factor),
    shape = 22, fill = NA, size = 2.5, alpha = 1,
    position = position_jitter(height = 0.15, width = 0, seed = 42),
    show.legend = FALSE
  ) +
  geom_vline(xintercept = 0, linetype = "solid", color = "black") +
  geom_errorbar(
    data  = plot_data_gtex,
    aes(x = Mean, xmin = CI_Lower, xmax = CI_Upper, y = Factor),
    color = "black",
    width = 0.2, show.legend = FALSE
  ) +
  # Mean summary — filled circle with BLACK OUTLINE
  geom_point(
    data  = plot_data_gtex,
    aes(x = Mean, y = Factor, fill = Factor),
    shape = 21, size = 5, alpha = 0.8,
    color = "black", stroke = 0.8,
    show.legend = FALSE
  ) +
  scale_color_manual(values = color_mapping) +
  scale_fill_manual(values  = color_mapping) +
  scale_y_discrete(labels = c(
    "Spearman_Rho_GDF15"           = "GDF15",
    "Spearman_Rho_Senescence"      = "sen_mean",
    "Spearman_Rho_Proliferation"   = "prolif_mean",
    "Spearman_Rho_ATF4"            = "ATF4",
    "Spearman_Rho_ATF5"            = "ATF5",
    "Spearman_Rho_DDIT3"           = "DDIT3",
    "Factor1_vs_AGE_SpearmanRhos"  = "Factor1",
    "Factor2_vs_AGE_SpearmanRhos"  = "Factor2",
    "Factor3_vs_AGE_SpearmanRhos"  = "Factor3",
    "Factor4_vs_AGE_SpearmanRhos"  = "Factor4",
    "Factor5_vs_AGE_SpearmanRhos"  = "Factor5",
    "Factor6_vs_AGE_SpearmanRhos"  = "Factor6",
    "Factor7_vs_AGE_SpearmanRhos"  = "Factor7",
    "Factor8_vs_AGE_SpearmanRhos"  = "Factor8",
    "Factor9_vs_AGE_SpearmanRhos"  = "Factor9",
    "Factor10_vs_AGE_SpearmanRhos" = "Factor10",
    "Factor11_vs_AGE_SpearmanRhos" = "Factor11",
    "Factor12_vs_AGE_SpearmanRhos" = "Factor12"
  )) +
  scale_x_continuous(breaks = seq(-1, 1, by = 0.2),
                     minor_breaks = seq(-1, 1, by = 0.1)) +
  labs(x = "Spearman Rho Mean \u00b1 95% CI",
       y = "Factors and Genes",
       title = "Average Spearman Rho Values with SEM") +
  theme_minimal() + fig3d_theme

print(p_gtex)
ggsave(file.path(out_dir, "Figure_3D_GTEx_with_individual_points_BlackOutline.png"),
       plot = p_gtex, width = chosen_width, height = 8, dpi = 300)
message("  Saved: Figure_3D_GTEx_with_individual_points_BlackOutline.png")

# ============================================================================
# FIBROBLAST — PER-CONDITION SPEARMAN RHOS FOR ALL METRICS
# ============================================================================
message("\n=== Fibroblast: computing per-condition Spearman rhos ===")

fb_order <- c("GDF15", "sen_mean", "prolif_mean",
              "ATF4", "ATF5", "DDIT3",
              "Factor1", "Factor12", "Factor6",  "Factor2",
              "Factor4",  "Factor9",  "Factor8",
              "Factor3",  "Factor10", "Factor7",  "Factor5",
              "Factor11")
factor_levels_fb <- rev(fb_order)

fb_data_dir <- here("Results", "Fibroblast_lifespan", "Factor_Analysis",
                    "Generate_All_Figures", "Spearman_Rhos_Age_vs_Genes")

# ---- Factors 1–12 (from pre-saved per-sample CSVs) ----
per_cond_rhos <- map_dfr(1:12, function(fnum) {
  f <- file.path(fb_data_dir,
                 paste0("Factor", fnum, "_ISR_Scores_per_Sample_Age_DaysGrown.csv"))
  if (!file.exists(f)) { message("  NOT FOUND: ", basename(f)); return(NULL) }
  read.csv(f) %>%
    filter(!is.na(DaysGrown), !is.na(fa_Scores)) %>%
    group_by(intx) %>%
    summarise(rho = tryCatch(
                cor.test(DaysGrown, fa_Scores, method = "spearman")$estimate,
                error = function(e) NA_real_),
              n = n(), .groups = "drop") %>%
    mutate(Factor = paste0("Factor", fnum))
})
message("  Factors 1-12: ", nrow(per_cond_rhos), " per-condition rows")

# ---- Gene metrics (ATF4, ATF5, DDIT3, GDF15, sen_mean, prolif_mean) ----
# Load from raw expression + manifest, matching the Rmd approach exactly
message("\n=== Loading raw fibroblast expression for gene metrics ===")

manifest_path <- here("Data", "Fibroblast_lifespan", "Lifespan_Study_selected_data.csv")
exprs_path    <- here("Data", "Fibroblast_lifespan",
                      "GSE179848_processed_cell_lifespan_RNAseq_data.csv")

gene_rhos <- NULL

if (!file.exists(manifest_path) || !file.exists(exprs_path)) {
  message("  WARNING: manifest or expression file not found — skipping gene metrics")
  message("    manifest:   ", manifest_path)
  message("    expression: ", exprs_path)
} else {
  message("  Reading manifest...")
  manifest_fb <- read.csv(manifest_path, stringsAsFactors = FALSE) %>%
    filter(!is.na(RNAseq_sampleID)) %>%
    mutate(
      SampleID = paste0("Sample_", stringr::str_trim(as.character(RNAseq_sampleID)))
    ) %>%
    rename(
      Cell_line_HCvsSURF = Cell_line,
      Cell_Line = Cell_line_inhouse,
      DaysGrown  = Days_grown_Udays
    ) %>%
    mutate(
      Cell_Line = as.character(Cell_Line),
      Cell_Line = stringr::str_squish(Cell_Line),
      Cell_Line = na_if(Cell_Line, "")
    ) %>%
    filter(!is.na(Cell_Line)) %>%
    mutate(
      Group = case_when(
        Cell_Line %in% c("hFB12","hFB13","hFB14","hFB11") ~ "Control",
        Cell_Line %in% c("hFB6","hFB7","hFB8")            ~ "SURF1",
        TRUE                                               ~ "Ctrl_tech_rep"
      ),
      Experiment = case_when(
        grepl("Mutation_Control",        Cell_line_group) ~ "No_Tx",
        grepl("Normal_Control_ox21",     Cell_line_group) ~ "No_Tx",
        grepl("Mutation_DEX",            Cell_line_group) ~ "DEX",
        grepl("Normal_DEX",              Cell_line_group) ~ "DEX",
        grepl("Normal_Oligomycin_",      Cell_line_group) ~ "Oligo",
        grepl("Normal_Oligomycin\\+",    Cell_line_group) ~ "DEX_Oligo",
        grepl("Normal_mitoNUITs_",       Cell_line_group) ~ "mitoNUITs",
        grepl("Normal_mitoNUITs\\+",     Cell_line_group) ~ "DEX_mitoNUITs",
        grepl("Galactose",               Cell_line_group) ~ "Galactose",
        grepl("2DG",                     Cell_line_group) ~ "2DG",
        grepl("betahydroxybutyrate",     Cell_line_group) ~ "betahydroxybutyrate",
        grepl("Control_ox3",             Cell_line_group) ~ "ox3",
        grepl("Inhibition_ox21",         Cell_line_group) ~ "Contact_Inhibition",
        grepl("Inhibition_ox3",          Cell_line_group) ~ "Contact_Inhibition_ox3",
        TRUE ~ NA_character_
      ),
      intx = paste(Group, Experiment, sep = "_")
    ) %>%
    filter(Group %in% c("Control", "SURF1"))

  message("  Manifest: ", nrow(manifest_fb), " samples (Control + SURF1)")

  message("  Reading expression data (this may take a moment)...")
  exprs_raw <- read.csv(exprs_path, stringsAsFactors = FALSE, check.names = FALSE)
  colnames(exprs_raw)[1] <- "gene"

  # Genes needed for ATF4, ATF5, DDIT3, GDF15, sen_mean, prolif_mean
  target_genes <- c("GDF15", "ATF4", "ATF5", "DDIT3",
                    "CDKN1A", "CDKN2A", "CCND2",   # sen_mean
                    "MKI67",  "RRM2",   "TOP2A")    # prolif_mean

  available_genes <- intersect(target_genes, exprs_raw$gene)
  missing_genes   <- setdiff(target_genes, exprs_raw$gene)
  if (length(missing_genes) > 0)
    message("  WARNING: genes not found in expression matrix: ",
            paste(missing_genes, collapse = ", "))
  message("  Available target genes: ", paste(available_genes, collapse = ", "))

  # Subset to only target genes, pivot to per-sample format
  exprs_sub <- exprs_raw %>%
    filter(gene %in% available_genes) %>%
    column_to_rownames("gene") %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("SampleID")

  # Merge with manifest
  data_gene <- manifest_fb %>%
    select(SampleID, DaysGrown, intx) %>%
    inner_join(exprs_sub, by = "SampleID") %>%
    filter(!is.na(DaysGrown))

  message("  Merged data: ", nrow(data_gene), " samples, ",
          ncol(data_gene) - 3, " gene columns")

  # Scale each gene column (matching Rmd: scale across all retained samples)
  gene_cols <- intersect(available_genes, colnames(data_gene))
  data_gene[, gene_cols] <- scale(data_gene[, gene_cols])

  # Compute per-sample composite scores
  if (all(c("CDKN1A","CDKN2A","CCND2") %in% colnames(data_gene))) {
    data_gene <- data_gene %>%
      rowwise() %>%
      mutate(
        sen_mean = mean(c_across(all_of(c("CDKN1A","CDKN2A","CCND2"))), na.rm = TRUE)
      ) %>%
      ungroup()
    message("  Computed sen_mean (CDKN1A, CDKN2A, CCND2)")
  }

  if (all(c("MKI67","RRM2","TOP2A") %in% colnames(data_gene))) {
    data_gene <- data_gene %>%
      rowwise() %>%
      mutate(
        prolif_mean = mean(c_across(all_of(c("MKI67","RRM2","TOP2A"))), na.rm = TRUE)
      ) %>%
      ungroup()
    message("  Computed prolif_mean (MKI67, RRM2, TOP2A)")
  }

  # Per-condition Spearman rhos for each gene metric
  gene_metrics <- c("GDF15", "ATF4", "ATF5", "DDIT3")
  if ("sen_mean"    %in% colnames(data_gene)) gene_metrics <- c(gene_metrics, "sen_mean")
  if ("prolif_mean" %in% colnames(data_gene)) gene_metrics <- c(gene_metrics, "prolif_mean")

  gene_rhos <- map_dfr(gene_metrics, function(metric) {
    data_gene %>%
      filter(!is.na(.data[[metric]]), !is.na(DaysGrown)) %>%
      group_by(intx) %>%
      summarise(
        rho = tryCatch(
          cor.test(DaysGrown, .data[[metric]], method = "spearman")$estimate,
          error = function(e) NA_real_),
        n = n(),
        .groups = "drop"
      ) %>%
      mutate(Factor = metric)
  })

  message("  Gene metrics: ", nrow(gene_rhos), " per-condition rows (",
          n_distinct(gene_rhos$Factor), " metrics)")
}

# Combine factors + gene metrics
if (!is.null(gene_rhos)) {
  per_cond_rhos <- bind_rows(per_cond_rhos, gene_rhos)
}

per_cond_rhos <- per_cond_rhos %>%
  filter(Factor %in% fb_order, !is.na(rho)) %>%
  mutate(Factor = factor(Factor, levels = factor_levels_fb))

message("\n  Per-condition rows for plot: ", nrow(per_cond_rhos),
        " (", n_distinct(per_cond_rhos$Factor), " factors/metrics, ",
        n_distinct(per_cond_rhos$intx), " conditions)")

# Overall summary (CI) — from Spearman_Rho_Summary_With_CI.csv
fb_summary <- read.csv(
  file.path(fb_data_dir, "Spearman_Rho_Summary_With_CI.csv"),
  stringsAsFactors = FALSE
) %>%
  filter(Column %in% fb_order) %>%
  mutate(Column = factor(Column, levels = factor_levels_fb))

# ---- Fibroblast plot: per-condition squares with mean circle (black-outlined) ----
p_fb <- ggplot() +
  geom_point(
    data  = per_cond_rhos,
    aes(x = rho, y = Factor, fill = Factor),
    shape = 22, color = NA, size = 2.5, alpha = 0.25,
    position = position_jitter(height = 0.15, width = 0, seed = 42),
    show.legend = FALSE
  ) +
  geom_point(
    data  = per_cond_rhos,
    aes(x = rho, y = Factor, color = Factor),
    shape = 22, fill = NA, size = 2.5, alpha = 1,
    position = position_jitter(height = 0.15, width = 0, seed = 42),
    show.legend = FALSE
  ) +
  geom_vline(xintercept = 0, linetype = "solid", color = "gray50") +
  geom_errorbar(
    data  = fb_summary,
    aes(x = Spearman_Rho, xmin = CI_Lower, xmax = CI_Upper,
        y = Column),
    color = "black",
    width = 0.2, show.legend = FALSE
  ) +
  # Mean summary — filled circle with BLACK OUTLINE
  geom_point(
    data  = fb_summary,
    aes(x = Spearman_Rho, y = Column, fill = Column),
    shape = 21, size = 5, alpha = 0.8,
    color = "black", stroke = 0.8,
    show.legend = FALSE
  ) +
  scale_color_manual(values = fb_color_mapping) +
  scale_fill_manual(values  = fb_color_mapping) +
  scale_x_continuous(breaks = seq(-1, 1, by = 0.5)) +
  labs(x = "Spearman Rho",
       y = "Columns",
       title = "Spearman Rho Values for Each Column vs DaysGrown") +
  theme_minimal() + fig3d_theme

print(p_fb)
ggsave(file.path(out_dir, "Figure_3D_Fibroblast_with_individual_points_BlackOutline.png"),
       plot = p_fb, width = chosen_width, height = 8, dpi = 300)
message("  Saved: Figure_3D_Fibroblast_with_individual_points_BlackOutline.png")

# ============================================================================
# EXPORT SOURCE DATA CSVs
# ============================================================================
message("\n=== Exporting source data CSVs ===")

# GTEx CSV: per-tissue rhos with tissue-level mean + 95% CI
gtex_csv <- tissue_long %>%
  rename(Factor_name = Factor, Per_tissue_Spearman_Rho = Spearman_Rho_tissue) %>%
  left_join(
    plot_data_gtex %>% rename(Factor_name = Factor,
                              Mean_Spearman_Rho = Mean,
                              SEM_Spearman_Rho  = SEM),
    by = "Factor_name"
  ) %>%
  arrange(Factor_name, Tissue)

write.csv(gtex_csv,
          file.path(out_dir, "Figure_3D_GTEx_source_data_with_individual_points.csv"),
          row.names = FALSE)
message("  Saved: Figure_3D_GTEx_source_data_with_individual_points.csv (",
        nrow(gtex_csv), " rows)")

# Fibroblast CSV: per-condition rhos with overall rho + CI
fb_csv <- per_cond_rhos %>%
  as.data.frame() %>%
  rename(Factor_name = Factor, Condition = intx,
         Per_condition_Spearman_Rho = rho, N_samples = n) %>%
  left_join(
    fb_summary %>% as.data.frame() %>%
      select(Column, Spearman_Rho, CI_Lower, CI_Upper,
             P_Value, Adjusted_P_Value, Significant) %>%
      rename(Factor_name = Column, Overall_Spearman_Rho = Spearman_Rho),
    by = "Factor_name"
  ) %>%
  arrange(Factor_name, Condition)

write.csv(fb_csv,
          file.path(out_dir, "Figure_3D_Fibroblast_source_data_with_individual_points.csv"),
          row.names = FALSE)
message("  Saved: Figure_3D_Fibroblast_source_data_with_individual_points.csv (",
        nrow(fb_csv), " rows)")

# ============================================================================
# SUMMARY
# ============================================================================
message(paste0("\n", strrep("=", 65)))
message("FIGURE 3D — JOURNAL-COMPLIANT PLOTS COMPLETE")
message(strrep("=", 65))
message("\nOutput folder: ", out_dir)
message("\nGTEx plot:  individual dots = per-tissue Spearman rho (",
        n_distinct(tissue_long$Tissue), " tissues)")
message("Fibro plot: individual dots = per-condition rho (",
        n_distinct(per_cond_rhos$Factor), " factors/metrics)")
message(strrep("=", 65))
