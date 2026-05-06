


###################################################################################################
ISR <- read.csv(here("Data", "Fibroblast_lifespan", "ISR_Gene_Lists_updated.csv")) # using isr list that includes all of the genes for FB data (igtp is IRGM, NARS to NARS1, WARS to WARS1)
###################################################################################################


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

long_ISR <- ISR %>%
  pivot_longer(cols = everything(), names_to = "Source", values_to = "Gene") %>%
  distinct()  # Remove duplicates if necessary

long_ISR$Gene <- toupper(long_ISR$Gene)

wide_ISR <- long_ISR %>%
  mutate(Presence = "Yes") %>%
  pivot_wider(names_from = Source, values_from = Presence, values_fill = list(Presence = "No")) %>%
  filter(!is.na(Gene) & Gene != "")  # This line removes rows where 'Gene' is NA or an empty string

Total_List <- wide_ISR %>%
  select(Gene) %>%
  filter(!is.na(Gene) & Gene != "")  # This line removes rows where 'Gene' is NA or an empty string



only_GO <- wide_ISR %>%
  select(Gene, GO_consortium) 

# Filter the dataframe to keep rows where GO_consortium is "Yes"
only_GO <- only_GO[only_GO$GO_consortium == "Yes", ] %>%
  select(Gene)


only_AnyGenes <- wide_ISR %>%
  select(Gene, AnyGenes) 

# Filter the dataframe to keep rows where GO_consortium is "Yes"
only_AnyGenes <- only_AnyGenes[only_AnyGenes$AnyGenes == "Yes", ] %>%
  select(Gene)


only_JacksonLabs <- wide_ISR %>%
  select(Gene, JacksonLabs) 

# Filter the dataframe to keep rows where GO_consortium is "Yes"
only_JacksonLabs <- only_JacksonLabs[only_JacksonLabs$JacksonLabs == "Yes", ] %>%
  select(Gene)



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


only_GO_plus_gdf15 <- c(only_GO$Gene, new_gene)
only_GO_plus_gdf15 <- as.character(only_GO_plus_gdf15)


only_AnyGenes_plus_gdf15 <- c(only_AnyGenes$Gene, new_gene)
only_AnyGenes_plus_gdf15 <- as.character(only_AnyGenes_plus_gdf15)


only_JacksonLabs_plus_gdf15 <- c(only_JacksonLabs$Gene, new_gene)
only_JacksonLabs_plus_gdf15 <- as.character(only_JacksonLabs_plus_gdf15)
