# ============================================================================
# Tissue Proliferation Index — Data Processing Script
# ============================================================================
# Calculates a tissue-level proliferation score from GTEx TPM data using
# three proliferation markers: MKI67, RRM2, TOP2A.
#
# Method:
#   1. Per sample: sum the raw TPM of MKI67, RRM2, TOP2A
#   2. Per tissue: average the raw sums across donors
#   3. log2(tissue average)
#
# Output: proliferation_score is log2-scaled.
# Per-donor cache stores RAW sums (not log2) for downstream plotting.
#
# Saves two CSV files to Data/gtex/Processed/:
#   - proliferation_index_Jack_Devine.csv           (original GTEx tissue names)
#   - proliferation_index_Jack_Devine_namesChanged.csv (names matching ISR scripts)
#
# Also caches per-donor scores to:
#   - Data/gtex/Processed/proliferation_per_donor.csv
#
# Prerequisites:
#   - GTEx TPM file: Data/gtex/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz
#   - GTEx sample annotations TSV
# ============================================================================

library(tidyverse)
library(here)

if (.Platform$OS.type == "unix") {
  Sys.setenv("R_MAX_VSIZE" = "32Gb")
}
# ============================================================================
# OUTPUT DIRECTORY
# ============================================================================
out_dir <- here("Data", "gtex", "Processed")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ============================================================================
# INPUT FILES
# ============================================================================
gtex_file <- here("Data", "gtex", "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz")

annotations_file <- here("Data", "gtex", "Insert_GTEx_SampleAttributes_and_SubjectPhenotypes_Here",
                         "GTEx_Analysis_2017-06-05_v8_Annotations_GTEx_Analysis_2017-06-05_v8_Annotations_SampleAttributesDS.tsv")

if (!file.exists(gtex_file)) {
  stop(paste0(
    "\n============================================================\n",
    "ERROR: GTEx TPM file not found!\n\n",
    "Expected location:\n  ", gtex_file, "\n\n",
    "This file must be downloaded from the GTEx Portal:\n",
    "  https://gtexportal.org/home/datasets\n\n",
    "Look for: GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz\n",
    "============================================================\n"
  ))
}

if (!file.exists(annotations_file)) {
  stop(paste0(
    "\n============================================================\n",
    "ERROR: Sample annotations file not found!\n\n",
    "Expected location:\n  ", annotations_file, "\n",
    "============================================================\n"
  ))
}

# ============================================================================
# READ GTEx TPM DATA
# ============================================================================
message("Reading GTEx TPM file (this may take a few minutes)...")
gtex <- read.delim(gtex_file, skip = 2)

exprs <- gtex %>%
  dplyr::select(-Name) %>%
  rename("Gene" = Description)

# Filter out genes with 0 for every value
exprs <- exprs %>%
  filter(rowSums(across(everything(), ~ . == 0)) != (ncol(exprs) - 1))

message("  Loaded ", nrow(exprs), " genes across ", ncol(exprs) - 1, " samples")

# ============================================================================
# READ SAMPLE ANNOTATIONS
# ============================================================================
message("Reading sample annotations...")
annotations_sample_attributes <- read_tsv(annotations_file, show_col_types = FALSE)

annotations_sample_attributes <- annotations_sample_attributes %>%
  filter(SMAFRZE == "RNASEQ") %>%
  mutate(SUBJID = substring(SAMPID, 1, 10))

annotations_sample_attributes <- annotations_sample_attributes %>%
  mutate(SUBJID_tmp = case_when(
    substr(SUBJID, nchar(SUBJID), nchar(SUBJID)) == "-" ~ substr(SUBJID, 1, nchar(SUBJID) - 1),
    .default = SUBJID
  )) %>%
  dplyr::select(-SUBJID) %>%
  rename(SUBJID = SUBJID_tmp)

annotations_sample_attributes <- annotations_sample_attributes %>%
  rename("X" = SAMPID)
annotations_sample_attributes$X <- gsub("\\-", ".", annotations_sample_attributes$X)

message("  Loaded ", nrow(annotations_sample_attributes), " sample annotations")

# ============================================================================
# CALCULATE PROLIFERATION INDEX
# Using 3 proliferation markers: MKI67, RRM2, TOP2A
# ============================================================================
message("\nCalculating proliferation index from MKI67, RRM2, TOP2A...")
message("  Method: sum(3 raw TPMs) per sample -> average across donors -> log2")

exprs_proliferative <- exprs %>%
  as.data.frame() %>%
  filter(Gene %in% c("MKI67", "RRM2", "TOP2A"))

message("  Found ", nrow(exprs_proliferative), " proliferation genes")

# Transpose: samples as rows, genes as columns
exprs_proliferative <- exprs_proliferative %>%
  as.data.frame() %>%
  column_to_rownames("Gene") %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("X")

# Step 1: sum raw TPM of 3 genes per sample (NO log2 yet)
per_sample_means <- exprs_proliferative %>%
  pivot_longer(cols = c("MKI67", "RRM2", "TOP2A"), names_to = "Gene", values_to = "TPM") %>%
  group_by(X) %>%
  summarise(sum_TPM = sum(TPM, na.rm = TRUE), .groups = "drop") %>%
  mutate(avg_prolif = sum_TPM)   # raw sum stored; log2 applied after tissue averaging

# Step 2: Join with annotations, apply RIN filter
per_donor_data <- per_sample_means %>%
  inner_join(annotations_sample_attributes, by = "X") %>%
  select(X, avg_prolif, SMTSD, SUBJID, SMRIN) %>%
  filter(SMRIN >= 5.5) %>%
  distinct()

message("  Donors after RIN >= 5.5 filter: ", n_distinct(per_donor_data$SUBJID))

# ============================================================================
# SAVE PER-DONOR SCORES (cached for Figure_4D_4E.R)
# ============================================================================
# Standardise tissue names for the per-donor cache
per_donor_cache <- per_donor_data %>%
  mutate(Tissue = gsub(" - ", "_", SMTSD)) %>%
  mutate(Tissue = gsub(" ", "_",  Tissue)) %>%
  mutate(Tissue = gsub("-", "_",  Tissue)) %>%
  mutate(Tissue = gsub("\\(", "", Tissue)) %>%
  mutate(Tissue = gsub("\\)", "", Tissue)) %>%
  mutate(Tissue = tolower(Tissue)) %>%
  mutate(Tissue = case_when(
    Tissue == "brain_frontal_cortex_ba9"              ~ "brain_frontal_cortex",
    Tissue == "skin_sun_exposed_lower_leg"            ~ "skin_lower_leg",
    Tissue == "skin_not_sun_exposed_suprapubic"       ~ "skin_suprapubic",
    TRUE ~ Tissue
  )) %>%
  select(SAMPID = X, Tissue, avg_prolif)

per_donor_file <- file.path(out_dir, "proliferation_per_donor.csv")
write.csv(per_donor_cache, per_donor_file, row.names = FALSE)
message("  Saved per-donor cache: ", per_donor_file)
message("  Rows: ", nrow(per_donor_cache))

# ============================================================================
# TISSUE-LEVEL PROLIFERATION INDEX
# avg_prolif contains RAW sums (sum of 3 TPMs per sample).
# Step 2: Average raw sums across donors per tissue.
# Step 3: Apply log2 to the tissue averages.
# ============================================================================
tissue_level <- per_donor_data %>%
  select(avg_prolif, SMTSD, SUBJID) %>%
  pivot_wider(names_from = SMTSD, values_from = avg_prolif) %>%
  column_to_rownames("SUBJID")

proliferation_index <- tissue_level %>%
  colMeans(na.rm = TRUE) %>%
  as.matrix() %>%
  as.data.frame() %>%
  mutate(V1 = log2(V1)) %>%          # Step 3: log2 AFTER averaging raw sums
  arrange(desc(V1)) %>%
  rename("proliferation_score" = V1) %>%
  rownames_to_column("Tissue")

message("  NOTE: proliferation_score = log2(mean of raw sums across donors)")
message("  Number of tissues: ", nrow(proliferation_index))

# ============================================================================
# SAVE OUTPUT - Original format (GTEx tissue names)
# ============================================================================
output_file <- file.path(out_dir, "proliferation_index_Jack_Devine.csv")
write.csv(proliferation_index, output_file, row.names = FALSE)
message("\nSaved: ", output_file)

# ============================================================================
# SAVE OUTPUT - With tissue names matching ISR analysis scripts
# ============================================================================
proliferation_index_renamed <- proliferation_index %>%
  mutate(Tissue = gsub(" - ", "_", Tissue)) %>%
  mutate(Tissue = gsub(" ", "_", Tissue)) %>%
  mutate(Tissue = gsub("-", "_", Tissue)) %>%
  mutate(Tissue = gsub("\\(", "", Tissue)) %>%
  mutate(Tissue = gsub("\\)", "", Tissue)) %>%
  mutate(Tissue = tolower(Tissue)) %>%
  mutate(Tissue = case_when(
    Tissue == "brain_frontal_cortex_ba9" ~ "brain_frontal_cortex",
    Tissue == "brain_nucleus_accumbens_basal_ganglia" ~ "brain_nucelus_accumbens_basal_ganglia",
    Tissue == "brain_spinal_cord_cervical_c_1" ~ "brain_spinal_cord_cervical_c1",
    Tissue == "breast_mammary_tissue" ~ "breast_mammary",
    Tissue == "minor_salivary_gland" ~ "salivary_gland",
    Tissue == "muscle_skeletal" ~ "skeletal_muscle",
    Tissue == "skin_not_sun_exposed_suprapubic" ~ "skin_suprapubic",
    Tissue == "skin_sun_exposed_lower_leg" ~ "skin_lower_leg",
    TRUE ~ Tissue
  ))

output_file_renamed <- file.path(out_dir, "proliferation_index_Jack_Devine_namesChanged.csv")
write.csv(proliferation_index_renamed, output_file_renamed, row.names = FALSE)
message("Saved: ", output_file_renamed)

message("\n============================================================")
message("TISSUE PROLIFERATION INDEX COMPLETE")
message("  Method: sum(3 raw TPMs) per sample -> average across donors -> log2")
message("  Output is log2-scaled (no further transformation needed)")
message("  Tissues: ", nrow(proliferation_index))
message("  Per-donor cache: ", per_donor_file)
message("============================================================")
