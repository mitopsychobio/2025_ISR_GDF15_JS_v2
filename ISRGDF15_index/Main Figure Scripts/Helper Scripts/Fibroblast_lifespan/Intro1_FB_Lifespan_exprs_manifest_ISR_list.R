


###################################################################################################
# Code with this aesthetic means code that may need to be changed
###################################################################################################

#rm(list = ls())
# install.packages("data.table")

# Load the data.table package
library(data.table)
library(tidyverse)
library(RColorBrewer)
library(circlize)
# NOTE: ComplexHeatmap was previously loaded here but is not used in this
# helper. Loading it has the side effect of masking pheatmap::pheatmap with
# a ComplexHeatmap translation wrapper that drops the `filename`/`width`/
# `height` arguments and renders with different defaults — which broke
# Figure 1E's heatmap. Removed.
# library(ComplexHeatmap)
library(ggplot2)
library(tidyr)
library(purrr)
library(plotly)
library(dplyr)
library(htmlwidgets)
library(broom)
library(limma)
library(ggrepel)
library(gridExtra)  # For arranging multiple plots

#MCP_choice <- "BH"                  #---------------------------------------------------------------------------------------------------
MCP_choice <- "bonferroni"         # CHOOSE type of Multiple Comparisons Procedure to use  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#---------------------------------------------------------------------------------------------------

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# manifest <- read.csv(here("Data", "Fibroblast_lifespan", "Meta_Lifespan_RNAseq_timepoints.csv"))
manifest <- read.csv(here("Data", "Fibroblast_lifespan", "Lifespan_Study_selected_data.csv"))


library(dplyr)
library(stringr)

# --- safety: check column exists
if (!"RNAseq_sampleID" %in% names(manifest)) {
  stop("Column 'RNAseq_sampleID' not found in manifest.")
}

# --- build cleaned table
keep_cols <- c(
  "Cell_line_group",
  "Unique_variable_name",
  "Cell_line",
  "Cell_line_inhouse",
  "Cell_type",
  "Sex",
  "Clinical_condition",
  "Donor_age",
  "Study_part",
  "Replicate_line",
  "Treatments",
  "Treatment_description",
  "Percent_oxygen",
  "Notes",
  "pre_designated_time_point_Udays",
  "Passage",
  "pre_study_passages_Udivisions",
  "Days_grown_Udays"
)

manifest <- manifest %>%
filter(!is.na(RNAseq_sampleID))

manifest <- manifest %>%
  # 1) create SampleID from RNAseq_SampleID
  mutate(
    SampleID = as.character(`RNAseq_sampleID`),
    # remove any existing "Sample_" prefix so we don't double it
    SampleID = str_remove(SampleID, "^Sample_"),
    # trim whitespace and then add the prefix
    SampleID = paste0("Sample_", str_trim(SampleID))
  ) %>%
  # 2) keep only requested columns (with SampleID first)
  select(SampleID, all_of(keep_cols))

# Optional: if you’d like to align names with your `manifest` later:
manifest <- manifest %>%
  rename(DaysGrown = Days_grown_Udays)

# Take a peek
dplyr::glimpse(manifest)


RNAseq_exp_log2 <- read.csv(here("Data", "Fibroblast_lifespan", "GSE179848_processed_cell_lifespan_RNAseq_data.csv"))  # THIS is LOG 2 TRANSFORMED !!!!! <<<<<<<<
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# 
# # Add the column names as new rows at the top of the dataframe
# manifest <- rbind(colnames(manifest), manifest)
# 
# # Rename the columns
# colnames(manifest) <- c("RNAseq_num", "Sample", "SampleID") # only if we use the csv file without the timepoints

manifest <- manifest %>%
  rename(Cell_line_HCvsSURF = Cell_line) %>%
  rename(Cell_Line = Cell_line_inhouse)

manifest <- manifest %>%
  mutate(
    Cell_Line = as.character(Cell_Line),        # in case it's a factor
    Cell_Line = str_squish(Cell_Line),          # trim + collapse internal spaces
    Cell_Line = na_if(Cell_Line, ""),           # "" -> NA
    Cell_Line = case_when(                      # common placeholders -> NA
      !is.na(Cell_Line) & str_to_lower(Cell_Line) %in% c("na","n/a","none","null") ~ NA_character_,
      TRUE ~ Cell_Line
    )
  ) %>%
  filter(!is.na(Cell_Line))


# Create a new column "Group" based on the values in the "Cell_Line" column
manifest$Group <- ifelse(manifest$Cell_Line %in% c("hFB12", "hFB13", "hFB14", "hFB11"), "Control", 
                  ifelse(manifest$Cell_Line %in% c("hFB6", "hFB7", "hFB8"), "SURF1", "Ctrl_tech_rep"))


colnames(RNAseq_exp_log2)[1] <- "gene"

# RNAseq_exp_log2_long <- pivot_longer(RNAseq_exp_log2, cols = -gene, names_to = "SampleID", values_to = "Expression") # Convert from wide to long format
# Exp_long <- merge(RNAseq_exp_log2_long, manifest, by = "SampleID") # Merge with manifest to get additional details for each sample
# 
# #Remove rows with NA values in the "Expression" column
# Exp_long <- Exp_long %>%
#   filter(!is.na(Expression))


manifest <- manifest %>%
  mutate(Experiment = case_when(
    grepl("Mutation_Control", Cell_line_group) ~ "No_Tx",
    grepl("Normal_Control_ox21", Cell_line_group) ~ "No_Tx",
    grepl("Mutation_DEX", Cell_line_group) ~ "DEX",
    grepl("Normal_DEX", Cell_line_group) ~ "DEX",
    grepl("Normal_Oligomycin_", Cell_line_group) ~ "Oligo",
    grepl("Normal_Oligomycin+", Cell_line_group) ~ "DEX_Oligo",
    grepl("Normal_mitoNUITs_", Cell_line_group) ~ "mitoNUITs",
    grepl("Normal_mitoNUITs+", Cell_line_group) ~ "DEX_mitoNUITs",
    grepl("Galactose", Cell_line_group) ~ "Galactose",
    grepl("2DG", Cell_line_group) ~ "2DG",
    grepl("betahydroxybutyrate", Cell_line_group) ~ "betahydroxybutyrate",
    grepl("Control_ox3", Cell_line_group) ~ "ox3",
    grepl("Inhibition_ox21", Cell_line_group) ~ "Contact_Inhibition",
    grepl("Inhibition_ox3", Cell_line_group) ~ "Contact_Inhibition_ox3",
    TRUE ~ NA_character_  # Default value for non-matching cases
  ))

manifest <- manifest %>%
  group_by(Cell_Line, Experiment) %>%
  mutate(Sample_Rep = paste0("Rep", row_number())) %>%
  ungroup()





# using annas code to make the meta data (the explanations of the samples) and the expression data (the expression of the genes)
data_all <- RNAseq_exp_log2 %>%
  pivot_longer(cols = -gene, names_to = "SampleID") %>%
  mutate(SampleID = as.character(SampleID)) %>%
  full_join(manifest %>% mutate(SampleID = as.character(SampleID)), by = "SampleID") %>%
  filter(!is.na(value)) %>%
  group_by(gene) %>%
  pivot_wider(names_from = "gene", values_from = "value")

################################################################################################### 
# data_sub <- data_all %>%                         # Choose only sp2
#   filter(grepl("sp2", Sample))

data_sub <- data_all                               # Choose all experiments
###################################################################################################

meta <- data_sub[,1:15]
exprs_wide <- data_sub[-c(2:15)] %>%
  na.omit(exprs_wide)
original_row_names <- rownames(exprs_wide) # Store the original row names
exprs <- t(exprs_wide)
exprs <- as.data.frame(exprs)
colnames(exprs) <- original_row_names
exprs_wide <- data_sub[-c(2:15)] %>%
  na.omit(exprs_wide)

colnames(exprs) <- exprs[1, ]  # Set the column names to the values in the first row
exprs <- exprs[-1, ]           # Remove the first row from the dataframe

###################################################################################################
# using isr list that includes all of the genes for FB data (igtp is IRGM, NARS to NARS1, WARS to WARS1)
###################################################################################################
ISR <-  read.csv(here("Data", "Fibroblast_lifespan", "ISR_Gene_Lists_updated.csv"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

long_ISR <- ISR %>%
  pivot_longer(cols = everything(), names_to = "Source", values_to = "Gene") %>%
  distinct()  # Remove duplicates if necessary

long_ISR$Gene <- toupper(long_ISR$Gene)

wide_ISR <- long_ISR %>%
  mutate(Presence = "Yes") %>%
  pivot_wider(names_from = Source, values_from = Presence, values_fill = list(Presence = "No"))

Total_List <- wide_ISR %>%
  select(Gene) %>%
  filter(!is.na(Gene) & Gene != "")  # This line removes rows where 'Gene' is NA or an empty string

opt1 <- Total_List$Gene
genes_of_interest <- opt1

# ------------- DO I WANT TO ADD IN GDF15 -----------------------
genes_of_interest <- as.character(genes_of_interest)

# # If we wish to add GDF15
new_gene <- "GDF15"

genes_of_interest_plus_gdf15 <- c(genes_of_interest, new_gene)

genes_of_interest_plus_gdf15 <- as.character(genes_of_interest_plus_gdf15)
genes_of_interest_plus_gdf15 <- genes_of_interest_plus_gdf15[genes_of_interest_plus_gdf15 != ""]
# genes_of_interest_plus_gdf15 <- c(genes_of_interest_plus_gdf15, "GDF15")
# genes_of_interest_plus_gdf15 <- as.character(genes_of_interest_plus_gdf15)



# Function to process datasets


process_datasets <- function(data_sub, genes_of_interest) {
  WT_SURF_No_Tx <- data_sub %>%
    select(SampleID, Experiment, Group, DaysGrown, intx, any_of(genes_of_interest)) %>%
    filter(Experiment == "No_Tx", Group %in% c("Control", "SURF1")) %>%
    na.omit()
  
  Control_No_Tx <- data_sub %>%
    select(SampleID, Experiment, Group, DaysGrown, intx,  any_of(genes_of_interest)) %>%
    filter(Experiment == "No_Tx", Group == "Control") %>%
    na.omit()
  
  SURF1_No_Tx <- data_sub %>%
    select(SampleID, Experiment, Group, DaysGrown, intx,  any_of(genes_of_interest)) %>%
    filter(Experiment == "No_Tx", Group == "SURF1") %>%
    na.omit()
  
  WT_SURF1_Oligo_Dex_No_Tx <- data_sub %>%
    select(SampleID, Experiment, Group, DaysGrown, intx,  any_of(genes_of_interest)) %>%
    filter(Experiment %in% c("No_Tx", "DEX", "Oligo"), Group %in% c("Control", "SURF1")) %>%
    na.omit()
  
  WT_SURF1_all_Txs <- data_sub %>%
    select(SampleID, Group, Experiment, DaysGrown, intx,  any_of(genes_of_interest)) %>%
    filter(Group %in% c("Control", "SURF1")) # %>%
  # na.omit()
  
  OxPhos_deficient_only <- data_sub %>%
    select(SampleID, Experiment, Group, DaysGrown, intx, any_of(genes_of_interest)) %>%
    filter(Experiment %in% c("No_Tx", "Oligo"), Group %in% c("Control", "SURF1")) %>%
    filter(!(Experiment == "No_Tx" & Group == "Control")) %>%
    na.omit()
  
  OxPhos_deficient_and_WT_No_Tx <- data_sub %>%
    select(SampleID, Experiment, Group, DaysGrown, intx,  any_of(genes_of_interest)) %>%
    filter(Experiment %in% c("No_Tx", "Oligo"), Group %in% c("Control", "SURF1")) %>%
    na.omit()
  
  OxPhos_deficient_and_WT_No_Tx_and_Contact_Inhib <- data_sub %>%
    select(SampleID, Experiment, Group, DaysGrown, intx,  any_of(genes_of_interest)) %>%
    filter(Experiment %in% c("No_Tx", "Oligo", "Contact_Inhibition", "ox3", "Contact_Inhibition_ox3"), Group %in% c("Control", "SURF1")) %>%
    na.omit()
  
  datasets <- list(
    WT_SURF_No_Tx = WT_SURF_No_Tx,
    Control_No_Tx = Control_No_Tx,
    SURF1_No_Tx = SURF1_No_Tx,
    WT_SURF1_Oligo_Dex_No_Tx = WT_SURF1_Oligo_Dex_No_Tx,
    WT_SURF1_all_Txs = WT_SURF1_all_Txs,
    OxPhos_deficient_only = OxPhos_deficient_only,
    OxPhos_deficient_and_WT_No_Tx = OxPhos_deficient_and_WT_No_Tx,
    OxPhos_deficient_and_WT_No_Tx_and_Contact_Inhib = OxPhos_deficient_and_WT_No_Tx_and_Contact_Inhib
  )
  
  return(datasets)
}
# How the above function can work:
# data_sub <- your_dataframe_here # Replace with your actual data
# 
# datasets_opt1 <- process_datasets(data_sub, opt1)
# datasets_other_opt <- process_datasets(data_sub, other_opt)



#-----------------------------------------------------------------------------------------------------------------------------------
# What to Keep or remove!
#-----------------------------------------------------------------------------------------------------------------------------------
all_objects <- ls()

# Specify the objects you want to keep
objects_to_keep <- c("df", "exprs", "Total_List", "opt1", "genes_of_interest_plus_gdf15", "manifest", "data_sub", "path", "process_datasets", "folder_path", "Folder_name")

# Create a list of objects to remove
objects_to_remove <- setdiff(all_objects, objects_to_keep)

# Remove all objects except the ones specified
rm(list = objects_to_remove)
rm(all_objects)
rm(objects_to_keep, objects_to_remove)
