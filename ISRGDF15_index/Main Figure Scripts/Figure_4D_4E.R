# ============================================================================
# Figures 4D & 4E: Journal-Compliant Dot Plots with Individual Donor Points
# Journal requirement: show individual data points behind tissue-level summaries
# ============================================================================
# Generates FOUR plots (two tissue sets × two analytes):
#   [ISR]           Figure_4D_4E_ISR_all_tissues.png
#   [ISR sig]       Figure_4D_4E_ISR_significant_tissues.png
#   [Proliferation] Figure_4D_4E_Proliferation_all_tissues.png
#   [Prolif sig]    Figure_4D_4E_Proliferation_significant_tissues.png
#
# ISR (Factor1) scores come from All_Tissue_data_DTHPLCE.csv (per-donor).
# Proliferation scores (sum of MKI67 + RRM2 + TOP2A TPM per donor) are
# computed on first run from the GTEx TPM file, then cached to:
#   Data/gtex/Processed/proliferation_per_donor.csv
#
# X-axis: tissues, ordered by descending mean score (each plot uses its own metric)
# Dots:   individual donors, translucent and jittered
# Line:   crossbar at tissue mean, opaque
# ============================================================================

library(here)
library(tidyverse)
library(ggplot2)
library(purrr)

# ============================================================================
# SOURCE FIGURE 4D SCRIPT TO GET data_sig (significant tissues)
# ============================================================================
message("\n=== Sourcing Figure 4D script to get data_sig ===")
source(here("Main Figure Scripts", "Helper Scripts", "Figure_4D_Base.R"), local = TRUE)

# Re-define output directory AFTER sourcing (source script clears environment)
out_dir <- here("Results", "Figures", "Figure_4D_4E")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# data_sig: tissue-level summary with Tissue, Avg_Factor1, proliferation_score,
#           Significant, deathplace_comparison, Proliferative_Group
# Standardise tissue names to lowercase
data_sig <- data_sig %>% mutate(Tissue = tolower(Tissue))
sig_tissues <- unique(data_sig$Tissue)

message("Significant tissues in data_sig: ", length(sig_tissues))

# ============================================================================
# LOAD PER-DONOR ISR (FACTOR1) SCORES
# ============================================================================
# Build the per-donor / per-tissue long-format dataframe by reading the
# per-tissue Plot_Data_<tissue>.csv files written by
# Figure_3B_4A_Tissue_Scatter.R. Helper applies NO DTHPLCE filtering; this
# script's main donor_isr below uses every donor (no DTHPLCE restriction),
# while the supplemental Figure 4E DTHPLCE breakdown filters internally.
message("\n=== Loading per-donor ISR scores ===")
source(here("Main Figure Scripts", "Helper Scripts",
            "Load_All_Tissue_DTHPLCE.R"))
all_tissue_data <- load_all_tissue_dthplce()

donor_isr <- all_tissue_data %>%
  filter(fa_vs_gdf15 == "Factor1") %>%
  select(SAMPID, Tissue, value) %>%
  rename(score = value) %>%
  mutate(Tissue = tolower(Tissue))

message("Per-donor ISR rows loaded: ", nrow(donor_isr))
message("Unique tissues in ISR data: ", n_distinct(donor_isr$Tissue))

# ============================================================================
# LOAD OR COMPUTE PER-DONOR PROLIFERATION SCORES
# Cached after first run to avoid re-reading the 2.5 GB GTEx TPM file.
# ============================================================================
per_donor_prolif_file <- here("Data", "gtex", "Processed", "proliferation_per_donor.csv")

if (file.exists(per_donor_prolif_file)) {

  message("\n=== Loading cached per-donor proliferation scores ===")
  donor_prolif <- read.csv(per_donor_prolif_file)
  message("Rows loaded: ", nrow(donor_prolif))

  # Info: raw sums of 0 will become -Inf after log2 at plot time.
  n_zero <- sum(donor_prolif$avg_prolif == 0, na.rm = TRUE)
  if (n_zero > 0) {
    message("  Note: ", n_zero, " samples have raw sum = 0 (all 3 genes = 0 TPM)")
  }

} else {

  message("\n=== Computing per-donor proliferation scores from GTEx TPM ===")
  message("    (This only runs once; results are cached for future use)")

  gtex_tpm_file  <- here("Data", "gtex", "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz")
  annotations_file <- here("Data", "gtex", "Insert_GTEx_SampleAttributes_and_SubjectPhenotypes_Here",
                           "GTEx_Analysis_2017-06-05_v8_Annotations_GTEx_Analysis_2017-06-05_v8_Annotations_SampleAttributesDS.tsv")

  if (!file.exists(gtex_tpm_file)) {
    stop(paste0(
      "\n============================================================\n",
      "ERROR: GTEx TPM file not found.\n",
      "Expected: ", gtex_tpm_file, "\n",
      "Per-donor proliferation scores cannot be computed without it.\n",
      "============================================================\n"
    ))
  }
  if (!file.exists(annotations_file)) {
    stop(paste0(
      "\n============================================================\n",
      "ERROR: GTEx sample annotations file not found.\n",
      "Expected: ", annotations_file, "\n",
      "============================================================\n"
    ))
  }

  message("  Reading GTEx TPM file (may take several minutes)...")
  gtex <- read.delim(gtex_tpm_file, skip = 2)
  exprs <- gtex %>%
    dplyr::select(-Name) %>%
    rename(Gene = Description)

  message("  Filtering to proliferation genes: MKI67, RRM2, TOP2A ...")
  exprs_prolif_genes <- exprs %>%
    as.data.frame() %>%
    filter(Gene %in% c("MKI67", "RRM2", "TOP2A"))

  if (nrow(exprs_prolif_genes) < 3) {
    stop("ERROR: Fewer than 3 proliferation genes found in GTEx TPM file.")
  }

  # Transpose: samples as rows, genes as columns
  exprs_prolif_t <- exprs_prolif_genes %>%
    column_to_rownames("Gene") %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("X")

  # Sum raw TPM of 3 genes per sample (NO log2 — raw sums cached; log2 at plot time)
  exprs_prolif_summed <- exprs_prolif_t %>%
    pivot_longer(cols = c("MKI67", "RRM2", "TOP2A"), names_to = "Gene", values_to = "TPM") %>%
    group_by(X) %>%
    summarise(sum_TPM = sum(TPM, na.rm = TRUE), .groups = "drop") %>%
    mutate(avg_prolif = sum_TPM)

  # Load sample annotations
  message("  Reading sample annotations...")
  ann <- read_tsv(annotations_file, show_col_types = FALSE) %>%
    filter(SMAFRZE == "RNASEQ") %>%
    mutate(SUBJID = substring(SAMPID, 1, 10)) %>%
    mutate(SUBJID = case_when(
      substr(SUBJID, nchar(SUBJID), nchar(SUBJID)) == "-" ~ substr(SUBJID, 1, nchar(SUBJID) - 1),
      TRUE ~ SUBJID
    ))

  # Normalise SAMPID format to match TPM column names
  ann <- ann %>%
    mutate(X = gsub("\\-", ".", SAMPID)) %>%
    select(X, SMTSD, SUBJID, SMRIN)

  # Join and apply RIN filter (matching Tissue_proliferation_index.R)
  donor_prolif_long <- exprs_prolif_summed %>%
    inner_join(ann, by = "X") %>%
    filter(SMRIN >= 5.5) %>%
    select(SAMPID = X, SUBJID, Tissue_raw = SMTSD, avg_prolif) %>%
    distinct()

  # Standardise tissue names (same transformations as Tissue_proliferation_index.R)
  donor_prolif_long <- donor_prolif_long %>%
    mutate(Tissue = gsub(" - ", "_", Tissue_raw)) %>%
    mutate(Tissue = gsub(" ", "_",  Tissue)) %>%
    mutate(Tissue = gsub("-", "_",  Tissue)) %>%
    mutate(Tissue = gsub("\\(", "", Tissue)) %>%
    mutate(Tissue = gsub("\\)", "", Tissue)) %>%
    mutate(Tissue = tolower(Tissue)) %>%
    mutate(Tissue = case_when(
      Tissue == "brain_frontal_cortex_ba9"              ~ "brain_frontal_cortex",
      Tissue == "skin_sun_exposed_lower_leg"            ~ "skin_lower_leg",
      Tissue == "skin_not_sun_exposed_suprapubic"       ~ "skin_suprapubic",
      Tissue == "colon_sigmoid"                         ~ "colon_sigmoid",
      Tissue == "colon_transverse"                      ~ "colon_transverse",
      Tissue == "esophagus_gastroesophageal_junction"   ~ "esophagus_gej",
      Tissue == "heart_left_ventricle"                  ~ "heart_left_ventricle",
      Tissue == "heart_atrial_appendage"                ~ "heart_atrial_appendage",
      TRUE ~ Tissue
    ))

  # Cache per-donor proliferation scores (RAW sums, not log2).
  # log2 is applied at plot time.
  donor_prolif <- donor_prolif_long %>%
    select(SAMPID, Tissue, avg_prolif)

  # Cache for future runs
  write.csv(donor_prolif, per_donor_prolif_file, row.names = FALSE)
  message("  Cached per-donor proliferation scores to: ", per_donor_prolif_file)
  message("  Rows cached: ", nrow(donor_prolif))

}

# ============================================================================
# HARMONISE PROLIFERATION TISSUE NAMES TO MATCH data_sig
# The cached proliferation_per_donor.csv uses GTEx annotation names which
# differ from the ISR analysis names for several tissues.
# ============================================================================
donor_prolif <- donor_prolif %>%
  mutate(Tissue = tolower(Tissue)) %>%
  mutate(Tissue = case_when(
    Tissue == "brain_nucleus_accumbens_basal_ganglia"  ~ "brain_nucelus_accumbens_basal_ganglia",
    Tissue == "brain_spinal_cord_cervical_c_1"         ~ "brain_spinal_cord_cervical_c1",
    Tissue == "breast_mammary_tissue"                  ~ "breast_mammary",
    Tissue == "esophagus_gej"                          ~ "esophagus_gastroesophageal_junction",
    Tissue == "esophagus_gastroesophageal_junction"     ~ "esophagus_gastroesophageal_junction",
    Tissue == "minor_salivary_gland"                   ~ "salivary_gland",
    Tissue == "muscle_skeletal"                        ~ "skeletal_muscle",
    TRUE ~ Tissue
  ))

# Verify matching
prolif_tissues_available <- unique(donor_prolif$Tissue)
missing_from_prolif <- setdiff(sig_tissues, prolif_tissues_available)
if (length(missing_from_prolif) > 0) {
  message("WARNING: ", length(missing_from_prolif), " significant tissues still missing from proliferation data:")
  message(paste("  -", missing_from_prolif, collapse = "\n"))
} else {
  message("All ", length(sig_tissues), " significant tissues found in proliferation data.")
}

# ============================================================================
# BUILD TISSUE-LEVEL SUMMARY TABLES (mean per tissue)
# For ordering and crossbar reference
# ============================================================================
tissue_means_isr <- donor_isr %>%
  group_by(Tissue) %>%
  summarise(mean_score = mean(score, na.rm = TRUE), .groups = "drop")

# Tissue mean for proliferation: avg_prolif contains RAW sums.
# Average raw sums across donors, then apply log2.
tissue_means_prolif <- donor_prolif %>%
  group_by(Tissue) %>%
  summarise(mean_score = log2(mean(avg_prolif, na.rm = TRUE)), .groups = "drop")

# Per-donor display scores: apply log2 to raw sums for plotting.
# Filter out zero-sum samples (log2(0) = -Inf) to avoid axis issues.
donor_prolif_plot <- donor_prolif %>%
  filter(avg_prolif > 0) %>%
  mutate(score = log2(avg_prolif))

# ============================================================================
# PLOTTING HELPER
# ============================================================================
make_dot_plot <- function(donor_data,          # data.frame: SAMPID, Tissue, score
                          tissue_means,         # data.frame: Tissue, mean_score
                          title,
                          ylab_text,
                          dot_color,
                          mean_color,
                          tissues_to_include,   # character vector of tissue names (lowercase)
                          order_by = "desc",
                          show_mean = TRUE) {

  # Restrict to requested tissues
  donor_sub  <- donor_data    %>% filter(Tissue %in% tissues_to_include)
  means_sub  <- tissue_means  %>% filter(Tissue %in% tissues_to_include)

  if (nrow(donor_sub) == 0 || nrow(means_sub) == 0) {
    warning("No data found for the requested tissues in: ", title)
    return(NULL)
  }

  # Compute y-axis limits from data with 10% padding
  y_range <- range(donor_sub$score, na.rm = TRUE)
  y_pad   <- diff(y_range) * 0.1
  y_lim   <- c(y_range[1] - y_pad, y_range[2] + y_pad)

  # Order tissues by descending mean score
  tissue_order <- means_sub %>%
    arrange(desc(mean_score)) %>%
    pull(Tissue)

  donor_sub <- donor_sub %>%
    mutate(Tissue = factor(Tissue, levels = tissue_order))

  means_sub <- means_sub %>%
    mutate(Tissue = factor(Tissue, levels = tissue_order))

  p <- ggplot() +
    # Individual donor dots — translucent, jittered
    geom_jitter(
      data  = donor_sub,
      aes(x = Tissue, y = score),
      color = dot_color, alpha = 0.25, size = 1.8,
      width = 0.22, height = 0
    )

  # Optionally add tissue mean crossbar
  if (show_mean) {
    p <- p +
      geom_crossbar(
        data  = means_sub,
        aes(x = Tissue, y = mean_score, ymin = mean_score, ymax = mean_score),
        color = mean_color, width = 0.55, linewidth = 0.7
      )
  }

  p <- p +
    scale_y_continuous(limits = y_lim, expand = c(0, 0)) +
    scale_x_discrete(expand = expansion(add = 0.6)) +
    labs(title = title, x = "Tissue", y = ylab_text) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title    = element_text(face = "bold"),
      axis.text.x   = element_text(angle = 45, hjust = 1, vjust = 1, colour = "black"),
      axis.text.y   = element_text(colour = "black"),
      axis.title    = element_text(colour = "black"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      panel.grid.major.y = element_line(colour = "grey92"),
      plot.margin   = margin(t = 10, r = 20, b = 5.5, l = 15)
    )

  return(p)
}


# ============================================================================
# FIGURE 4E + SUPPLEMENTAL FIGURE 8
# ============================================================================
# Both plots are computed on `data_sig` (significant tissues only — built by
# Figure_4D_Base.R after the median-proliferation split). Each is a
# two-group dot plot of Avg_Factor1 (per-tissue ISR^GDF15 score):
#
#   Figure 4E             — split by deathplace_comparison
#                           (which side of the H/E comparison the tissue fell on)
#   Supplemental Fig 8    — split by Proliferative_Group
#                           (above vs below median proliferation_score)
#
# Two CSVs (Prism wide format, no tissue labels) and two PNGs are written
# to Results/Figures/Figure_4D_4E/. Figure_4E.png is also copied to the
# shared Results/Figures/ folder; Supplemental_Fig8.png is not.
# ============================================================================
message("\n=== Generating Figure 4E + Supplemental Figure 8 ===")

# --- Helper: build a Prism-style wide CSV (one column per group, no labels) ---
prism_wide <- function(df, group_col, value_col, group_levels, col_names) {
  out <- lapply(group_levels, function(g) df[[value_col]][df[[group_col]] == g])
  n   <- max(lengths(out))
  out <- lapply(out, function(x) c(x, rep(NA_real_, n - length(x))))
  result <- as.data.frame(out, stringsAsFactors = FALSE)
  colnames(result) <- col_names
  result
}

# --- Helper: two-group dot plot in the reference-image style ---
two_group_dot_plot <- function(df, group_col, value_col,
                                group_levels, group_labels, group_colors,
                                title_text, y_lab) {
  df_local <- df %>%
    filter(.data[[group_col]] %in% group_levels) %>%
    mutate(.grp = factor(.data[[group_col]],
                         levels = group_levels,
                         labels = group_labels))

  # Wilcoxon p + Hedges g
  # Sign convention: cohen.d returns positive when level-1's mean exceeds
  # level-2's. We display |g| on the plot (direction is visually obvious),
  # and print the full-precision signed value to the console so the
  # caller can sanity-check the magnitude.
  vals <- split(df_local[[value_col]], df_local[[group_col]])
  v1 <- vals[[group_levels[1]]]
  v2 <- vals[[group_levels[2]]]

  wilcox_p <- tryCatch(
    wilcox.test(v1, v2)$p.value,
    error = function(e) NA_real_
  )
  hedges <- tryCatch(
    effsize::cohen.d(df_local[[value_col]],
                     factor(df_local[[group_col]], levels = group_levels),
                     hedges.correction = TRUE)$estimate,
    error = function(e) NA_real_
  )
  hedges_g <- unname(hedges)

  # ---- Console diagnostics: full-precision summary so 1.00 displayed on
  # ---- the plot can be checked against the underlying value.
  message("  -- ", group_col, " --")
  message("    n(", group_levels[1], ") = ", length(v1),
          ", mean = ", signif(mean(v1, na.rm = TRUE), 5),
          ", median = ", signif(median(v1, na.rm = TRUE), 5))
  message("    n(", group_levels[2], ") = ", length(v2),
          ", mean = ", signif(mean(v2, na.rm = TRUE), 5),
          ", median = ", signif(median(v2, na.rm = TRUE), 5))
  message("    Wilcoxon p (raw) = ", signif(wilcox_p, 6))
  message("    Hedges g (signed, level1 - level2) = ",
          signif(hedges_g, 6))
  message("    |Hedges g| (plot label) = ", signif(abs(hedges_g), 6))

  p_label <- ifelse(is.na(wilcox_p), "ns",
              ifelse(wilcox_p < 0.0001, "****",
                ifelse(wilcox_p < 0.001,  "***",
                  ifelse(wilcox_p < 0.01,   "**",
                    ifelse(wilcox_p < 0.05,    "*", "ns")))))
  # Display with 3 decimals so a true 0.987 / 1.013 isn't masked as "1.00".
  g_label <- if (!is.na(hedges_g)) {
    paste0("g = ", format(round(abs(hedges_g), 3), nsmall = 3))
  } else "g = NA"

  y_max <- max(df_local[[value_col]], na.rm = TRUE)
  y_min <- min(df_local[[value_col]], na.rm = TRUE)
  y_pad <- (y_max - y_min) * 0.18

  ggplot(df_local, aes(x = .grp, y = .data[[value_col]], color = .grp)) +
    geom_jitter(width = 0.15, height = 0, size = 4, alpha = 0.85) +
    stat_summary(fun = median, geom = "crossbar",
                 width = 0.55, linewidth = 0.6, fatten = 1.5) +
    # Significance bracket + stars + g annotation, stacked above data
    annotate("segment",
             x = 1, xend = 2,
             y = y_max + y_pad * 0.7, yend = y_max + y_pad * 0.7,
             color = "black", linewidth = 0.5) +
    annotate("text",
             x = 1.5, y = y_max + y_pad * 1.1,
             label = p_label,
             size = 8, fontface = "bold") +
    annotate("text",
             x = 1.5, y = y_max + y_pad * 0.4,
             label = g_label,
             size = 5) +
    scale_color_manual(values = setNames(group_colors, group_labels),
                       guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.18))) +
    labs(title = title_text, x = NULL, y = y_lab) +
    theme_classic(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5,
                                 margin = margin(b = 10)),
      axis.text  = element_text(color = "black"),
      axis.title = element_text(color = "black")
    )
}

# ---------------------------------------------------------------------------
# Figure 4E: proliferation_score split by deathplace_comparison
# (Matches the Prism reference image with g = 1.17.)
# ---------------------------------------------------------------------------
fig4e_levels  <- c("ER>HI", "HI>ER")
fig4e_labels  <- c("ER > HI", "HI > ER")
fig4e_colors  <- c("#A04C7A", "#E89B3C")  # magenta-pink, orange

# CSV (Prism wide)
fig4e_csv <- prism_wide(
  df           = data_sig %>% filter(deathplace_comparison %in% fig4e_levels),
  group_col    = "deathplace_comparison",
  value_col    = "proliferation_score",
  group_levels = fig4e_levels,
  col_names    = c("ER_greater_than_HI", "HI_greater_than_ER")
)
fig4e_csv_path <- file.path(out_dir, "4E_for_prism.csv")
write.csv(fig4e_csv, fig4e_csv_path, row.names = FALSE, na = "")
message("  Saved: ", fig4e_csv_path)

# Plot
fig4e_plot <- two_group_dot_plot(
  df            = data_sig,
  group_col     = "deathplace_comparison",
  value_col     = "proliferation_score",
  group_levels  = fig4e_levels,
  group_labels  = fig4e_labels,
  group_colors  = fig4e_colors,
  title_text    = "Proliferation scores by place of death",
  y_lab         = "Proliferation score (log2)"
)
fig4e_png <- file.path(shared_out_dir, "Figure_4E.png")
ggsave(fig4e_png, plot = fig4e_plot, width = 5, height = 6, dpi = 300)
message("  Saved: ", fig4e_png)


# ---------------------------------------------------------------------------
# Supplemental Figure 8: Avg_Factor1 split by Proliferative_Group
# ---------------------------------------------------------------------------
suppfig8_levels  <- c("Proliferative", "Non-Proliferative")
suppfig8_labels  <- c("Proliferative", "Non-\nproliferative")
suppfig8_colors  <- c("#A878B7", "#4D7BB7")  # purple, blue

suppfig8_csv <- prism_wide(
  df           = data_sig %>% filter(Proliferative_Group %in% suppfig8_levels),
  group_col    = "Proliferative_Group",
  value_col    = "Avg_Factor1",
  group_levels = suppfig8_levels,
  col_names    = c("Proliferative", "Non_proliferative")
)
suppfig8_csv_path <- file.path(out_dir, "Supplemental_Fig8.csv")
write.csv(suppfig8_csv, suppfig8_csv_path, row.names = FALSE, na = "")
message("  Saved: ", suppfig8_csv_path)

suppfig8_plot <- two_group_dot_plot(
  df            = data_sig,
  group_col     = "Proliferative_Group",
  value_col     = "Avg_Factor1",
  group_levels  = suppfig8_levels,
  group_labels  = suppfig8_labels,
  group_colors  = suppfig8_colors,
  title_text    = expression(bold("ISR"^"GDF15" ~
                                   "scores based on proliferation score")),
  y_lab         = expression("ISR"^"GDF15" ~ "Score")
)
suppfig8_png <- file.path(out_dir, "Supplemental_Fig8.png")
ggsave(suppfig8_png, plot = suppfig8_plot, width = 5, height = 6, dpi = 300)
message("  Saved: ", suppfig8_png)
# (Supplemental_Fig8.png is intentionally NOT copied to Results/Figures/.)

# ============================================================================
# SUMMARY
# ============================================================================
message(paste0("\n", strrep("=", 65)))
message("FIGURES 4D & 4E — JOURNAL-COMPLIANT DOT PLOTS COMPLETE")
message(strrep("=", 65))
message("\nOutput folder: ", out_dir)
message("  ", per_donor_prolif_file)
message("\nFigure 4E + Supplemental Figure 8 outputs:")
message("  ", fig4e_csv_path)
message("  ", fig4e_png)
message("  ", suppfig8_csv_path)
message("  ", suppfig8_png)
message(strrep("=", 65))
