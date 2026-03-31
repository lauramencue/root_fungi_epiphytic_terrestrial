# Project name: root_fungi_epiphytic_terrestrial
# Author: Laura Méndez
# Date: 15.01.2025


# Goals:
# This script prepares the annotated fungal dataset for downstream analyses by:
# 1. Loading the fully annotated ASV data from script 04
# 2. Filtering and subsetting data into analysis-specific datasets
# 3. Creating summary statistics on taxonomic diversity

# Input files:
# 1. guild_output.Rdata - Complete annotated dataset from script 04

# Output files:
# 1. analyses.Rdata - Workspace containing all prepared data objects

# Main data objects created:
# - fungi: Complete fungal dataset (all guilds)
# - metadata: Sample metadata
# - asv_tax_fungi: Taxonomic summary of all fungal ASVs

# Setup ------------------------------------------------------------------------

# Clear workspace
rm(list = ls())

# Set working directory (adjust as needed)
setwd("C:/Users/mendezcu_lu/Dropbox/Work/03_Marburg/02_epiphites_terrestrial_root_funghi/R/scripts")
# setwd("C:/Users/ploti/Dropbox/Work/03_Marburg/02_epiphites_terrestrial_root_funghi/R")

# Load required libraries
library(tidyverse)
library(openxlsx)

# Load data --------------------------------------------------------------------

load("guild_output.Rdata")

fungi <- guild_output$fungi
metadata <- guild_output$metadata

# Revise complete annotated dataset from script 03
length(unique(fungi$ASV))
length(unique(fungi$SampleID))
colnames(fungi)

# Rename columns
colnames(fungi)[19] <- "host_family"
colnames(fungi)[20] <- "host_genus"
colnames(fungi)[22] <- "host_wcvp_name"
colnames(fungi)[24] <- "host_substrate"
colnames(fungi)[25] <- "host_species_identity"

# Inspect data structure and distributions -------------------------------------

# Check FunctionalGuilds distribution

count(fungi, FunctionalGuilds, sort = TRUE) 

count(fungi, all_functional, sort = TRUE) 

count(fungi, Phylum, sort = TRUE)
count(fungi, Class, sort = TRUE)

# Taxonomic resolution summary -------------------------------------------------

# Get unique ASV-taxonomy combinations for all fungi
asv_tax_fungi <- fungi %>%
  filter(Kingdom == "k__Fungi" | Kingdom == "Fungi") %>%
  distinct(ASV, Phylum, Class, Order, Family, Genus_FungalTraits, Species, FunctionalGuilds,
           BroadSymbiontGroup, all_functional, is_symbiotic)


count(asv_tax_fungi, FunctionalGuilds, sort=T)

# Taxonomic resolution statistics
# Calculate proportion of ASVs assigned at each rank full fungi
tax_assignment_fungi <- tibble(
  Rank = c("Phylum", "Class", "Order", "Family", "Species"),
  Assigned_ASVs = c(
    sum(!is.na(asv_tax_fungi$Phylum)),
    sum(!is.na(asv_tax_fungi$Class)),
    sum(!is.na(asv_tax_fungi$Order)),
    sum(!is.na(asv_tax_fungi$Family)),
    sum(!is.na(asv_tax_fungi$Species))
  )
) %>%
  mutate(
    Total_ASVs = nrow(asv_tax_fungi),
    Percent = round(Assigned_ASVs / nrow(asv_tax_fungi) * 100, 1)
  )

tax_assignment_fungi

# Guild resolution summary -------------------------------------------------

# Calculate proportion of ASVs assigned at each functional guild, full ASVs

guild_coverage <- tibble(
  Metric = c("ASVs with guild", "ASVs without guild"),
  n = c(
    sum(!is.na(asv_tax_fungi$FunctionalGuilds)),
    sum(is.na(asv_tax_fungi$FunctionalGuilds))
  )
) %>%
  mutate(
    Total_ASVs = nrow(asv_tax_fungi),
    Percent = round(n / Total_ASVs * 100, 1)
  )

guild_coverage

symb_summary_all <- asv_tax_fungi %>%
  filter(!is.na(FunctionalGuilds)) %>%
  count(is_symbiotic, name = "ASVs") %>%
  mutate(
    Total_ASVs = sum(ASVs),
    Percent = round(ASVs / Total_ASVs * 100, 1)
  )

symb_summary_all

broad_group_summary <- asv_tax_fungi %>%
  filter(is_symbiotic) %>%
  count(BroadSymbiontGroup, name = "ASVs") %>%
  mutate(
    Total_Symbiotic_ASVs = sum(ASVs),
    Percent = round(ASVs / Total_Symbiotic_ASVs * 100, 1)
  ) %>%
  arrange(desc(ASVs))

broad_group_summary

guild_assignment_noNA <- asv_tax_fungi %>%
  filter(!is.na(FunctionalGuilds)) %>%
  count(FunctionalGuilds, name = "Assigned_ASVs") %>%
  mutate(
    Total_Assigned_ASVs = sum(Assigned_ASVs),
    Percent = round(Assigned_ASVs / Total_Assigned_ASVs * 100, 1)
  ) %>%
  arrange(desc(Assigned_ASVs))

guild_assignment_noNA

# Calculate ASV richness per sample --------------------------------------------

# ASV richness (alpha diversity) for all fungal ASVs
asv_per_sample <- fungi %>%
  distinct(host_species_identity, ASV) %>%
  count(host_species_identity, name = "ASV_richness") %>%
  arrange(desc(ASV_richness))


asv_richness_habitat <- fungi %>%
  group_by(host_substrate) %>%
  summarise(
    total_ASV_richness = n_distinct(ASV),
    n_samples = n_distinct(SampleID),
    .groups = "drop"
  )

asv_richness_habitat

str(fungi)
summary(asv_per_sample$ASV_richness)
count(asv_tax_fungi, Species, sort = TRUE)

# Create metadata table -------------------------------------------------------

colnames(metadata)

# Extract sample metadata for all fungi


colnames(metadata)[1] <- "host_family"
colnames(metadata)[2] <- "host_genus"
colnames(metadata)[4] <- "host_wcvp_name"
colnames(metadata)[6] <- "host_substrate"
colnames(metadata)[7] <- "host_species_identity"
colnames(metadata)[8] <- "SampleID"


metadata <- metadata %>%
  mutate(host_substrate = recode(
    host_substrate,
      "E" = "Epiphytic",
      "T" = "Terrestrial")) 

fungi <- fungi %>%
  mutate(host_substrate = recode(
    host_substrate,
    "E" = "Epiphytic",
    "T" = "Terrestrial")) 

cat(paste("Epiphytic samples:", sum(metadata$host_substrate == "Epiphytic"), "\n"))
cat(paste("Terrestrial samples:", sum(metadata$host_substrate == "Terrestrial"), "\n"))

# Sample distribution by family
count(metadata, host_family) 

# Sample distribution by growth form within families
count(metadata, host_family, host_substrate) 

# Species with both epiphytic and terrestrial samples
metadata %>%
  group_by(host_species_identity, host_family) %>%
  summarise(
    n_habitats = n_distinct(host_substrate),
    habitats = paste(sort(unique(host_substrate)), collapse = ", "),
    n_samples = n()
  ) %>%
  filter(n_habitats > 1)

################################################################################
# Save workspace ---------------------------------------------------------------
################################################################################

analyses <- list(
  asv_tax_fungi = asv_tax_fungi,
  fungi = fungi,
  metadata = metadata)


# Save all prepared data objects for use in subsequent analysis scripts
save(analyses, file="analyses.Rdata")


