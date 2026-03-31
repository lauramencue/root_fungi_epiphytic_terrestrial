# Project name: root_fungi_epiphytic_terrestrial
# Author: Laura Méndez
# Date: 15/12/2025

# Goals:
# Assign taxonomy to ASVs using the UNITE database for fungal ITS2 sequences
# This creates the taxonomic backbone for all downstream analyses

# Input files:
# 1. decontam_output.Rdata from script 02 containing:
#    - seqtab.full: Full ASV table (cleaned from contaminants)
#    - seqtab.sensitivity: Filtered ASV table (no singletons, for sensitivity analysis)
#    - metadata: Sample metadata
# 2. sh_general_release_dynamic_19.02.2025.fasta - UNITE reference database

# Output files:
# 1. taxonomy_output.Rdata containing:
#    - taxa.full: Taxonomic assignments for all ASVs
#    - boot.full: Bootstrap confidence values for full dataset
#    - seqtab.full: ASV abundance table (passed through)
#    - metadata: Sample metadata (passed through)

# References:
# - UNITE database: Abarenkov et al. (2024) Nucleic Acids Research
#   https://doi.org/10.1093/nar/gkad1039
# - DADA2 assignTaxonomy: Callahan et al. (2016) Nature Methods

################## --------------------------------------------------------------------
# Clear workspace
rm(list = ls())

# Load required libraries

library(dada2) # packageVersion("dada2") ‘1.38.0’
library(tidyverse)
library(openxlsx)
library(stringr)

# Load data

load("decontam_output.Rdata")

seqtab.full <- decontam_output$seqtab.full  # All ASVs (main analysis)
metadata <- decontam_output$metadata

# Path to UNITE reference database

unite.ref <- "input_data/sh_general_release_dynamic_19.02.2025.fasta"  

# The minimum bootstrap confidence for assigning a taxonomic level defaults to 50 (minBoot)
taxa_output <- assignTaxonomy(seqtab.full, unite.ref, multithread = TRUE, 
                              tryRC = TRUE, # try reverse-complement if no forward match
                              outputBootstraps = TRUE)

# Extract taxonomy and bootstrap values
taxa <- taxa_output[[1]]
boot <- taxa_output[[2]]

################################################################################
# Save results
################################################################################

tax_output <- list(
  taxa.full = taxa,
  boot.full = boot,
  seqtab.full = seqtab.full,
  metadata = metadata)

save(tax_output, file = "tax_output.Rdata")
