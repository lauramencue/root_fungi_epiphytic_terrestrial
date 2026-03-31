# Project name: root_fungi_epiphytic_terrestrial
# Author: Laura Méndez
# Date: 06/03/2026

# Goals:
# SENSITIVITY ANALYSIS: Create phyloseq object from ASV data and cluster into OTUs
# This combines ASV data, taxonomy, and metadata into a phyloseq object, then
# clusters ASVs into OTUs at 97% similarity. Taxonomy automatically transfers to OTUs.

# Input files:
# 1. decontam_output.Rdata from script 02 (cleaned ASV table)
# 2. tax_output.Rdata from script 03 (taxonomic assignments)

# Output files:
# 1. OTU_clustering_output.Rdata containing:
#    - otu_table: OTU abundance matrix
#    - cluster_assignments: Mapping of ASVs to OTUs
#    - otu_sequences: Representative sequences for each OTU

################## --------------------------------------------------------------------

# Clear workspace
rm(list = ls())

# Load required libraries

 
library(dplyr)
library(DECIPHER) # packageVersion("DECIPHER") ‘3.6.0’
library(phyloseq)
library(MiscMetabar) # packageVersion("MiscMetabar") ‘0.14.4’
library(Biostrings)
library("magrittr")

# Load data
load("decontam_output.Rdata")
load("tax_output.Rdata")

# Extract data -----------------------------------------------------------------

metadata <- decontam_output$metadata
seqtab.full <- decontam_output$seqtab.full
taxa.full <- tax_output$taxa.full

# Create ASV IDs
asv_ids <- paste0("ASV", seq_along(seqtab.full))

# Rename tables
colnames(seqtab.full) <- asv_ids
rownames(taxa.full)   <- asv_ids

# Build phyloseq object --------------------------------------------------------

otu <- otu_table(seqtab.full, taxa_are_rows = FALSE)
tax <- tax_table(as.matrix(taxa.full))
sam <- sample_data(metadata)

phylo_asv <- phyloseq(otu, tax, sam)

# Add refseq
dna <- DNAStringSet(seq_strings)
names(dna) <- asv_ids

phylo_asv <- merge_phyloseq(phylo_asv, dna)

# Quick visualization
summary_plot_pq(phylo_asv)

# Cluster ASVs into OTUs (97% similarity) --------------------------------------

set.seed(1)

# Use MiscMetabar::asv2otu function
# This automatically:
# - Clusters sequences at 97% similarity
# - Transfers taxonomy to OTUs
# - Aggregates abundances
# - Keeps all metadata
phylo_otu <- asv2otu(
  physeq = physeq_asv,
  cutoff = 0.03,      # 97% similarity (3% divergence)
  method = "complete", # Complete linkage clustering
  nproc = 4           # Number of processors
)

# Visualization of OTU dataset
summary_plot_pq(phylo_otu)

# Optional: Qucik Hill diversity comparison
phylo_otu@sam_data@names[6] <- "host_substrate"

q_otu <- MiscMetabar::hill_pq(phylo_otu, fact = "host_substrate",
                          color_fac = "host_substrate",
                          add_points = T,
                          one_plot = T)
q_otu

q_asv <- MiscMetabar::hill_pq(phylo_asv, fact = "host_substrate",
                           color_fac = "host_substrate",
                           add_points = T,
                           one_plot = T)
q_asv


# Transform back to a normal data frame

otu.table <- as.data.frame(otu_table(phylo_otu))
taxa.table <- as.data.frame(tax_table(phylo_otu))

rownames(taxa.table)

comm_mat_otu <- as(otu_table(otu.table), "matrix")

otu_long <- otu.table %>%
  rownames_to_column("SampleID") %>%
  pivot_longer(
    cols = -SampleID,
    names_to = "ASV",
    values_to = "Reads"
  ) %>%
  filter(Reads > 0) %>%
  left_join(
    taxa.table1 %>%
      rownames_to_column("ASV"),
    by = "ASV"
  )


otu_long <- otu_long %>%
  dplyr::rename(OTU = ASV) %>%
  mutate(OTU = sub("^ASV", "OTU", OTU))


length(unique(otu_long$OTU))

# Save results -----------------------------------------------------------------

otu_clustering_output <- list(
  comm_mat_otu = comm_mat_otu,
  otu_long = otu_long)

save(otu_clustering_output, file = "otu_clustering_output.Rdata")








