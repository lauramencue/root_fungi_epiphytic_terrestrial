# Project name: root_fungi_epiphytic_terrestrial
# Author: Laura Méndez
# Date: 15/12/2025

# Goals:
# This script assigns functional guilds to taxonomically identified ASVs using using 
# two complementary databases:
# 1. FUNGuild (Nguyen et al., 2016) - Guild assignments based on literature
# 2. FungalTraits (Põlme et al., 2020) - Trait-based ecological classifications
# The script then merges these assignments with the ASV table and sample metadata,
# categorizes guilds into broader functional groups, and creates the final annotated dataset.

# Input files:
# 1. taxonomy_output.Rdata from script 03 containing taxonomic assignments
# 2. FungalTraits 1.2_v...xlsx - FungalTraits database

# Output files:
# 1. guild_output.Rdata containing annotated ASV data ready for analysis

# References:
# - Nguyen et al. (2016) Fungal Ecology 20:241-248 (FUNGuild)
# - Põlme et al. (2020) Scientific Data 7:51 (FungalTraits)

################## --------------------------------------------------------------------
# Clear workspace
rm(list = ls())

# Load required libraries

library(FUNGuildR) # packageVersion("FUNGuildR") ‘0.3.0’
library(tidyverse)
library(openxlsx)
library(stringr)

# Load data

load("tax_output.Rdata")

taxa.full <- tax_output$taxa.full
seqtab.full <- tax_output$seqtab.full
metadata <- tax_output$metadata


# Prepare taxonomy string for FUNGuild -----------------------------------------

# FUNGuild requires a semicolon-separated taxonomy string
# Format: Kingdom;Phylum;Class;Order;Family;Genus;Species
fungi_fg <- taxa.full %>%
  as.data.frame() %>%
  mutate(
    # Concatenate all taxonomic ranks
    Taxonomy = paste(Kingdom, Phylum, Class, Order, Family, Genus, Species, sep = ";")
  ) %>%
  mutate(
    # Clean up the taxonomy string:
    # Remove ";NA" patterns (missing ranks)
    Taxonomy = str_replace_all(Taxonomy, ";NA", ""),
    # Remove "NA;" patterns at the beginning
    Taxonomy = str_replace_all(Taxonomy, "NA;", ""),
    # Remove multiple consecutive semicolons
    Taxonomy = str_replace_all(Taxonomy, ";;+", ";"),
    # Remove trailing semicolons
    Taxonomy = str_replace_all(Taxonomy, ";$", "")
  )


# Inspect data structure
colnames(fungi_fg)
nrow(taxa.full)
length(unique(fungi_fg$Phylum))
length(unique(fungi_fg$Genus))
length(unique(fungi_fg$Species))

# Assign functional guilds using FUNGuild --------------------------------------

# Run FUNGuild assignment
# tax_col specifies which column contains the taxonomy string
fungi_guilded <- funguild_assign(fungi_fg, tax_col = "Taxonomy")


# Inspect FUNGuild results
colnames(fungi_guilded)
count(fungi_guilded, confidenceRanking) 
count(fungi_guilded, trophicMode, sort = TRUE) 
count(fungi_guilded, guild, sort = TRUE)

fungi_guilded$Sequence <- rownames(fungi_fg)

# Load and prepare FungalTraits database ---------------------------------------

# Read FungalTraits database
traits <- read.xlsx(
  "../FungalTraits 1.2_vhttps___docs.google.com_spreadsheets_u_0__authuser=0&usp=sheets_weber_16Dec_2020.xlsx"
)

# Inspect FungalTraits database

count(traits, primary_lifestyle, sort = TRUE) 
count(traits, Endophytic_interaction_capability_template, sort = TRUE) 

# Clean taxonomic prefixes from DADA2 output -----------------------------------

# Remove UNITE-style prefixes (k__, p__, c__, etc.) from taxonomy
# These prefixes indicate the taxonomic rank but aren't needed for matching
fungi_clean <- fungi_fg %>%
  rownames_to_column("Sequence")%>%
  mutate(
    Phylum = str_remove(Phylum, "^p__"),
    Class = str_remove(Class, "^c__"),
    Order = str_remove(Order, "^o__"),
    Family = str_remove(Family, "^f__"),
    Genus = str_remove(Genus, "^g__"),
    Species = str_remove(Species, "^s__")
  )

# Join FungalTraits data at genus level ----------------------------------------

# Match fungal ASVs to FungalTraits database by genus
# This provides additional ecological trait information
fungi_traits_genus <- fungi_clean %>%
  left_join(traits, by = c("Genus" = "GENUS"))

# Inspect joined data
count(fungi_traits_genus, primary_lifestyle, sort = TRUE) 

count(fungi_traits_genus, Endophytic_interaction_capability_template, sort = TRUE)
 

################################################################################
# Merge FUNGuild and FungalTraits data -----------------------------------------
################################################################################

# Combine FUNGuild assignments with FungalTraits data
# Join by the row identifier (Sequence column)
fungi_all <- fungi_guilded %>%
  left_join(fungi_traits_genus, by = "Sequence")

################################################################################
# Create clean combined dataset ------------------------------------------------
################################################################################
colnames(fungi_all)
# Select and rename relevant columns for final dataset
# Keep one set of taxonomic columns (from FungalTraits join, .y suffix)
fungi_all_clean <- fungi_all %>%
  transmute(
    Sequence,
    Kingdom = Kingdom.x,                   # Taxonomic ranks
    Phylum = Phylum.x,
    Class = Class.x,
    Order = Order.x,
    Family = Family.x,
    Genus_FungalTraits = Genus.y,          # Genus from FungalTraits
    Genus_FUNGuild = taxon,                # Genus from FUNGuild
    Species = Species.y,
    confidenceRanking_FunGuild = confidenceRanking,
    trophicMode,                           # FUNGuild trophic mode
    guild,                                 # FUNGuild guild assignment
    primary_lifestyle,                      # FungalTraits primary lifestyle
    Secondary_lifestyle,                    # FungalTraits secondary lifestyle
    Endophytic_interaction_capability_template  # FungalTraits endophyte capability
  )

# Verify unique sequences

length(unique(fungi_all_clean$Sequence))

################################################################################
# Load ASV abundance table and merge with annotations --------------------------
################################################################################

### ASV table
seqtab_nochim_long <- seqtab.full %>%
  as.data.frame() %>%
  tibble::rownames_to_column("SampleID") %>%
  pivot_longer(
    cols = -SampleID,
    names_to = "Sequence",
    values_to = "Reads"
  ) %>%
  filter(Reads > 0)


asv_lookup <- seqtab_nochim_long %>%
  distinct(Sequence) %>%
  arrange(Sequence) %>%  
  mutate(ASV = paste0("ASV", row_number()))

seqtab_nochim_long2 <- seqtab_nochim_long %>%
  left_join(asv_lookup, by = "Sequence") 

# Join ASV abundances with taxonomic and functional annotations
ASV_annotated <- seqtab_nochim_long2 %>%
  left_join(fungi_all_clean, by = "Sequence")

# Check for missing data
sum(is.na(ASV_annotated$Kingdom))

# Verify no duplicate sample-ASV combinations
any(duplicated(ASV_annotated[, c("SampleID", "ASV")]))


################################################################################
# Load and merge sample metadata -----------------------------------------------
################################################################################

# Read sample metadata
metadata <- decontam_output$metadata
metadata <- metadata[, 1:10]  # Keep only first 10 columns with relevant data

# Merge ASV data with sample metadata
ASV_annotated <- ASV_annotated %>%
  left_join(metadata, by = c("SampleID" = "name_sequencing"))


################################################################################
# Assign functional guild categories -------------------------------------------
################################################################################

# Create hierarchical functional guild assignments
# Priority: Specific symbiotic guilds > Endophytes > Broad guilds
# This section consolidates multiple database assignments into clear categories


fungi_all <- ASV_annotated %>%
  mutate(FunctionalGuilds = case_when(
      
      # MYCORRHIZAL FUNGI
      # Arbuscular Mycorrhizal (AMF) 
      Phylum == "Glomeromycota" ~ "arbuscular_mycorrhizal",
      
      # Ericoid and Orchid Mycorrhizal (combined)
      str_detect(guild, "Ericoid Mycorrhizal-Orchid Mycorrhizal") ~
        "ericoid_orchid_mycorrhizal",
      
      # Orchid Mycorrhizal (OMF) 
      Family %in% c("Tulasnellaceae", "Ceratobasidiaceae", "Serendipitaceae", "Atractiellales_fam_Incertae_sedis") |
        str_detect(guild, "Orchid Mycorrhizal") ~
        "orchid_mycorrhizal",
      
      # Ericoid Mycorrhizal (ErMF) 
      str_detect(guild, "Ericoid Mycorrhizal") |
        primary_lifestyle == "ericoid_mycorrhizal" ~
        "ericoid_mycorrhizal",
      
      # Ectomycorrhizal (EcMF) 
      str_detect(guild, "Ectomycorrhizal") |
        primary_lifestyle == "ectomycorrhizal" ~
        "ectomycorrhizal",
      
      # ENDOPHYTIC FUNGI
      # Dark septate endophytes (DSE) 
      Endophytic_interaction_capability_template == "root_endophyte_dark_septate" &
        str_detect(guild, "Endophyte") &
        str_detect(trophicMode, "Symbiotroph") ~
        "dark_septate_endophyte",
      
      # MFRE / Endogonales 
      Order == "Endogonales" ~
        "fine_root_endophyte",
      
      # Other root endophytes 
      (Endophytic_interaction_capability_template %in% c("root_endophyte", "root-associated") &
          primary_lifestyle == "root_endophyte" &
          str_detect(coalesce(trophicMode, ""), regex("Symbiotroph", ignore_case = TRUE))) |
        (is.na(Endophytic_interaction_capability_template) &
            str_detect(coalesce(guild, ""), regex("Endophyte", ignore_case = TRUE)) &
            str_detect(coalesce(trophicMode, ""), regex("Symbiotroph", ignore_case = TRUE))) 
      ~ "root_endophyte",
      # Foliar endophytes 
      str_detect(coalesce(guild, ""), "Endophyte") &
        (primary_lifestyle == "foliar_endophyte" |
           Endophytic_interaction_capability_template == "foliar_endophyte") ~
        "foliar_endophyte",
      
      # If no specific category matches, keep as NA
      TRUE ~ NA_character_
    )
  ) %>%
  
  # Fill remaining NAs from guild (FUNGuild) and then from primary_lifestyle
  mutate(FunctionalGuilds = coalesce(FunctionalGuilds, primary_lifestyle, guild))

# Check functional guild distribution

count(fungi_all, FunctionalGuilds, sort = TRUE)  %>%
  print(n = Inf)
count(fungi_all, primary_lifestyle, sort = TRUE) 

################################################################################
# Consolidate functional guilds into broader categories -----------------------
################################################################################

# Group similar guilds into ecologically meaningful categories
# This simplifies analysis while retaining biological meaning
fungi <- fungi_all %>%
  mutate(
    FunctionalGuilds = case_when(
      
      # ANIMAL-ASSOCIATED
      str_detect(FunctionalGuilds, "Animal") |
      FunctionalGuilds %in% c("Endophyte-Insect Pathogen", "animal_parasite",
                              "animal_endosymbiont") ~ 
        "animal-associated",
      
      # PLANT PATHOGENS
      str_detect(FunctionalGuilds, "Plant Pathogen") |
      FunctionalGuilds %in% c(
        "Endophyte-Lichen Parasite-Plant Pathogen-Undefined Saprotroph",
        "Lichen Parasite-Plant Pathogen-Wood Saprotroph",
        "Fungal Parasite-Undefined Saprotroph", 
        "Fungal Parasite",
        "mycoparasite",
        "lichen_parasite",
        "Algal Parasite-Bryophyte Parasite-Fungal Parasite-Undefined Saprotroph",
        "Endophyte-Epiphyte-Fungal Parasite-Insect Parasite",
        "Plant Parasite-Wood Saprotroph",
        "unspecified_pathotroph",
        "algal_parasite",
        "Algal Parasite--Leaf Saprotroph-Wood Saprotroph"
      ) ~ "plant_pathogen",
      
      # SAPROTROPHS (wood, litter, soil decomposers)
      str_detect(FunctionalGuilds, "Dung Saprotroph") |
      FunctionalGuilds %in% c(
        "Plant Saprotroph-Wood Saprotroph",
        "Wood Saprotroph",
        "Leaf Saprotroph-Wood Saprotroph",
        "Litter Saprotroph-Wood Saprotroph",
        "Endophyte-Litter Saprotroph-Soil Saprotroph-Undefined Saprotroph",
        "litter_saprotroph",
        "unspecified_saprotroph",
        "soil_saprotroph",
        "nectar/tap_saprotroph",
        "dung_saprotroph",
        "pollen_saprotroph",
        "wood_saprotroph",
        "Undefined Saprotroph",
        "Epiphyte-Undefined Saprotroph",
        "Endophyte-Undefined Saprotroph",
        "Endophyte-Soil Saprotroph-Undefined Saprotroph",
        "Lichenized-Wood Saprotroph",
        "Plant Saprotroph",
        "Undefined Saprotroph-Wood Saprotroph"
      ) ~ "wood_litter_soil_saprotroph",
      
      # OTHER (miscellaneous groups)
      FunctionalGuilds %in% c(
        "Lichenized",
        "Lichenized-Undefined Saprotroph",
        "lichenized",
        "Epiphyte",
        "epiphyte",
        "sooty_mold",
        "unspecified",
        "protistan_parasite",
        "Lichen Parasite-Lichenized",
        "moss_symbiont"
      ) ~ "other",
      
      # Keep all other assignments as is
      TRUE ~ FunctionalGuilds
    )
  )


count(fungi, FunctionalGuilds, sort = TRUE)  %>%
  print(n = Inf)



#### Assign broader functional guild groups and symbiotic TRUE FALSE:

fungi <- fungi %>%
  mutate(
    BroadSymbiontGroup = case_when(
      FunctionalGuilds == "arbuscular_mycorrhizal" ~ "arbuscular_mycorrhizal",
      
      FunctionalGuilds %in% c(
        "ectomycorrhizal",
        "orchid_mycorrhizal",
        "ericoid_mycorrhizal",
        "ericoid_orchid_mycorrhizal"
      ) ~ "other_mycorrhizal_types",
      
      FunctionalGuilds == "dark_septate_endophyte" ~ "dark_septate_endophytes",
      
      FunctionalGuilds %in% c(
        "root_endophyte",
        "fine_root_endophyte"
      ) ~ "other_root_endophytes",
      
      TRUE ~ NA_character_
    ),
    is_symbiotic = !is.na(BroadSymbiontGroup),
    all_functional = coalesce(BroadSymbiontGroup, FunctionalGuilds)
  )


count(fungi, BroadSymbiontGroup, sort = TRUE)
count(fungi, is_symbiotic)
count(fungi, all_functional, sort = T)

# Final functional guild distribution

fungi %>%
  distinct(ASV, FunctionalGuilds) %>%
  count(FunctionalGuilds, sort = TRUE) %>%
  print()

fungi %>%
  distinct(ASV, BroadSymbiontGroup) %>%
  count(BroadSymbiontGroup, sort = TRUE) %>%
  print()

fungi %>%
  distinct(ASV, is_symbiotic) %>%
  count(is_symbiotic, sort = TRUE) %>%
  print()


length(unique(fungi$ASV))

################################################################################
# Save results
################################################################################

guild_output <- list(
  fungi = fungi,
  metadata = metadata)

save(guild_output, file = "guild_output.Rdata")
