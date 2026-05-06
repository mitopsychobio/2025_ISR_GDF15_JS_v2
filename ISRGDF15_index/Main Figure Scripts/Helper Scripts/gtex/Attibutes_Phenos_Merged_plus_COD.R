


#This code reads in the different sample attributes and subject phenotypes, and adds in cause of death column to the merged attributes and phenotpes. 


# BiocManager::install("edgeR")

library(tidyverse)
library(edgeR)
library(corrr)
library(dplyr)
library(tibble)
library(ggpubr)
library(readr)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(beepr)
library(here)   # <—— NEW: project-aware paths

# ------------------------------------------------------------------
# Project-aware paths (assumes your .Rproj is in 2024_08_29_fb_gtex)
# ------------------------------------------------------------------

# Optional, just to confirm in console:
here::here()

# Folder that holds the GTEx input files
gtex_input <- here("Data", "gtex")

# ISR gene list
ISR_read_in <- here("Data", "gtex", "Total_ISR_Gene_List_plus_gdf15.csv")
Total_ISR_List <- read.csv(ISR_read_in) %>%
  select(-X)

# ------------------------------------------------------------------
# Read annotations using here() paths
# ------------------------------------------------------------------

annotations_subject_phenotypes <- read_tsv(
  here("Data", "gtex", "Insert_GTEx_SampleAttributes_and_SubjectPhenotypes_Here",
       "GTEx_Analysis_2017-06-05_v8_Annotations_GTEx_Analysis_2017-06-05_v8_Annotations_SubjectPhenotypesDS.tsv")
)

annotations_sample_attributes <- read_tsv(
  here("Data", "gtex", "Insert_GTEx_SampleAttributes_and_SubjectPhenotypes_Here",
       "GTEx_Analysis_2017-06-05_v8_Annotations_GTEx_Analysis_2017-06-05_v8_Annotations_SampleAttributesDS.tsv")
) %>%
  filter(SMAFRZE == "RNASEQ") %>%
  mutate(SUBJID = substring(SAMPID, 1, 10)) %>%
  filter(SMRIN >= 6) %>%
  mutate(SUBJID_tmp = case_when(
    substr(SUBJID, nchar(SUBJID), nchar(SUBJID)) == "-" ~ substr(SUBJID, 1, nchar(SUBJID) - 1),
    TRUE ~ SUBJID
  )) %>%
  select(-SUBJID) %>%
  rename(SUBJID = SUBJID_tmp)

# Clean up IDs (no more setwd needed)
annotations_subject_phenotypes$SUBJID <- gsub("\\-", ".", annotations_subject_phenotypes$SUBJID)
annotations_sample_attributes$SAMPID <- gsub("\\-", ".", annotations_sample_attributes$SAMPID)

annotations_sample_attributes <- annotations_sample_attributes %>%
  mutate(SUBJID = sub("^([^.]+\\.[^.]+)\\..*$", "\\1", SAMPID))

annotations_merged <- annotations_sample_attributes %>%
  inner_join(annotations_subject_phenotypes, by = "SUBJID")

# Initialize a dataframe to store Spearman's rho and p-values for all tissues
all_results <- data.frame(
  Tissue = character(),
  PC = character(),
  Spearman_Rho = numeric(),
  P_Value = numeric(),
  stringsAsFactors = FALSE
)

# Path to the folder containing individual tissue folders
individual_tissues_folder <- here("2_InputFiles", "gtex", "All_Tissues_Indiv_Folders")

# ------------------------------------------------------------------
# Cause of death recoding
# ------------------------------------------------------------------

Death_causes <- annotations_subject_phenotypes %>%
  select(SUBJID, AGE, WGHT, HGHT, DTHMNNR, DTHCOD, DTHLUCOD, DTHFUCOD)

# Create combined_string column
Death_causes <- Death_causes %>%
  mutate(combined_string = paste(DTHMNNR, DTHCOD, DTHLUCOD, DTHFUCOD, sep = "_"))

# Define the function to determine COD
determine_cod <- function(combined_string) {
  if (grepl("suicide|Suicide", combined_string)) {
    return("suicide")
  } else if (grepl("homicide|Homicide", combined_string)) {
    return("homicide")
  } else if (grepl("accident|Accident", combined_string) & grepl("drug|OD|toxi|overdose", combined_string)) {
    return("OD")
  } else if (grepl("accident|Accident", combined_string) &
             grepl("MVA|mva|motor|vehicle", combined_string) &
             !grepl("non|NON", combined_string)) {
    return("accident")
  } else if (grepl("accident|Accident", combined_string) & grepl("fall|Fall|fell", combined_string)) {
    return("accident")
  } else if (grepl("Undetermined", combined_string) & grepl("overdose", combined_string)) {
    return("undetermined")
  } else if (grepl("Natural|natural", combined_string) & grepl("disease|Disease", combined_string)) {
    return("natural_disease")
  } else if (grepl("Natural|natural", combined_string) & grepl("unknown cause of death", combined_string)) {
    return("natural_unknown")
  } else {
    return("natural") # default
  }
}

# Create COD column using the function
Death_causes <- Death_causes %>%
  mutate(COD = sapply(combined_string, determine_cod))

Death_causes_selected <- Death_causes %>%
  select(SUBJID, COD)

# Join the selected columns from Death_causes to annotations_merged by SUBJID
annotations_merged <- annotations_merged %>%
  inner_join(Death_causes_selected, by = "SUBJID")

# ------------------------------------------------------------------
# Updated July 15th 2024 – DTHCOD recoding
# ------------------------------------------------------------------

#####################################    for DTHCOD    #########################################################
unique_death_causes <- Death_causes %>%
  select(DTHCOD) %>%
  distinct(DTHCOD)

unique_death_causes_new <- unique_death_causes %>%
  mutate(New_DTHCOD = case_when(
    grepl("Cardiac|cardiac|myocardial|Myocardial|arrest|cardiovascular|heart failure|acute mi|Massive Heart Attack|probable mi|heart attack", DTHCOD, ignore.case = TRUE) ~ "cardiac failure",
    grepl("stroke|Stroke|cerebrovascular|CVA|cva|Cva|Cerebral Vascular|cerebral vascular", DTHCOD, ignore.case = TRUE) ~ "cerebrovascular accident",
    grepl("COPD|copd|pulmonary edema", DTHCOD, ignore.case = TRUE) ~ "COPD",
    grepl("alcoholism", DTHCOD, ignore.case = TRUE) ~ "alcoholism",
    grepl("anoxia|Anoxia|ANOXIA|anoxic encephalopathy", DTHCOD, ignore.case = TRUE) ~ "anoxia",
    grepl("als|ALS", DTHCOD, ignore.case = TRUE) ~ "ALS",
    grepl("allergic reaction|allergy|allergic|allergies", DTHCOD, ignore.case = TRUE) ~ "allergic reaction",
    grepl("brain cancer", DTHCOD, ignore.case = TRUE) ~ "brain cancer",
    grepl("cancer", DTHCOD, ignore.case = TRUE) ~ "cancer",
    grepl("head trauma|Head trauma|Head Trauma|blunt injury|blunt force trauma|trauma|trauma due to fall|head trama", DTHCOD, ignore.case = TRUE) ~ "head trauma",
    grepl("heart disease|CAD|congestive failure heart", DTHCOD, ignore.case = TRUE) ~ "heart disease",
    grepl("Dementia", DTHCOD, ignore.case = TRUE) ~ "dementia",
    grepl("end stage liver disease|ESLD|ESRD|renal|liver disease", DTHCOD, ignore.case = TRUE) ~ "liver disease",
    grepl("kidney failure| kidney diseases", DTHCOD, ignore.case = TRUE) ~ "kidney disease",
    grepl("sirs|smoke inhalation-respiratory disease|respiratory diseases|lung disease", DTHCOD, ignore.case = TRUE) ~ "respiratory disease",
    grepl("unknown death|unknown", DTHCOD, ignore.case = TRUE) ~ "unknown",
    grepl("toxic effect of unspecified substance|poisoning|poisoning by overdose of substance|overdose", DTHCOD, ignore.case = TRUE) ~ "poison",
    grepl("suicide by hanging|suicide-hanging|strangulation|asphyiation due to hanging", DTHCOD, ignore.case = TRUE) ~ "suicide by hanging",
    grepl("seizure|epilepsy", DTHCOD, ignore.case = TRUE) ~ "seizures",
    TRUE ~ DTHCOD  # Default case to keep the original value if no match
  ))

unique_death_causes_new <- unique_death_causes_new %>%
  select(New_DTHCOD) %>%
  distinct(New_DTHCOD)


# ADDED IN NOVEMBER 26th! ################## ********************** +++++++++++++ |||||||||||||||||
# Create COD column using the function
# Death_causes <- Death_causes %>%
#   mutate(COD = sapply(combined_string, determine_cod))
# 
# COD_counts <- table(Death_causes$COD)
# print(COD_counts)
# 
# Death_causes_selected <- Death_causes %>%
#   select(SUBJID, COD)
# 
# # Join the selected columns from Death_causes to annotations_merged by SUBJID
# annotations_merged <- annotations_merged %>%
#   inner_join(Death_causes_selected, by = "SUBJID")
# 
# 
# 
# unique_death_causes_DTHFUCOD <- Death_causes %>%
#   distinct(DTHCOD, DTHFUCOD)
# 

annotations_merged <- annotations_merged %>%
  mutate(cardiac = case_when(
    grepl("cardiac|acute mi|myocardial|myocardiac|arrest; cardiac|bradycardia|heart|cardio|chf|aortic dissection|probable mi", 
          DTHCOD, ignore.case = TRUE) |
      grepl("cardiac|acute mi|myocardial|myocardiac|arrest; cardiac|bradycardia|heart|cardio|chf|aortic dissection|probable mi", 
            DTHFUCOD, ignore.case = TRUE) ~ "cardiac",
    TRUE ~ "non cardiac"
  ))


annotations_merged <- annotations_merged %>%
  mutate("heart_disease" = case_when(
    grepl("heart disease", 
          DTHCOD, ignore.case = TRUE) |
      grepl("heart disease", 
            DTHFUCOD, ignore.case = TRUE) ~ "cardiac",
    TRUE ~ "non cardiac"
  ))









#####################################    for DTHCOD    #########################################################


#--------------------------------------
# What to Keep or remove!
#--------------------------------------
all_objects <- ls()

# Specify the objects you want to remove (unchanged logic)
objects_to_remove <- c(
  "Death_causes_selected",
  "Death_causes",
  "annotations_sample_attributes",
  "annotations_subject_phenotypes",
  "determine_cod"
)

rm(list = objects_to_remove)
rm(all_objects)
rm(objects_to_remove)
# Alias for compatibility with scripts that expect "Atts_merged"
Atts_merged <- annotations_merged

