# Project name: root_fungi_epiphytic_terrestrial
# Author: Laura Méndez
# Date: 06.02.26

# Goals:
# This script processes ASV table output from dada2, identifies potential contaminants
# using a negative control and the DNA concentration per sample before sequencing and 
# deletes the potential ASV contaminants from the data sets

# Input files:
# 1. dada_ouput.Rdata from 02_DADA2.R that includes:
### - seqtab_full - Initial ASV (all ASVs after general filtering with dada2)
### - seqtab_sensitivity - Filtered ASVs without singletons for potential sensitivity analysis
### - track_reads_ASV_filtering - Read counts through pipeline steps
# 2. samples_DNA.xlsx - all collected samples with metadata

# Output files as decontam_ouput.Rdata:
# 1. seqtab.full = seqtab.clean - cleaned from contaminants full ASV table
# 2. seqtab.sensitivity = seqtab.filt.clean - cleaned for contaminants filtered for singletons ASV table
# 3. metadata - all collected samples with metadata


# References:
# - decontam general vignette: https://benjjneb.github.io/decontam/vignettes/decontam_intro.html#putting-it-all-together
# - Davis, N.M., Proctor, D.M., Holmes, S.P. et al. Simple statistical 
#   identification and removal of contaminant sequences in marker-gene and metagenomics data. 
#   Microbiome 6, 226 (2018). https://doi.org/10.1186/s40168-018-0605-2

################## --------------------------------------------------------------------

rm(list = ls())

# Load required libraries
library(openxlsx)
library(decontam) # packageVersion("decontam") ‘1.30.0’

## Read sample metadata

metadata <- read.xlsx("input_data/samples_DNA.xlsx")

########### Rename samples using the sample_map so they match the names in metadata

sample_map <- dada_output$sample_map

seqtab0 <- dada_output$seqtab_full

seqtab_filt <- dada_output$seqtab_sensitivity

c(n_ASVs_full = ncol(seqtab0), n_ASVs_no_singletons = ncol(seqtab_filt))

rownames(seqtab0) <- sub("_1\\.fq\\.gz$", "", rownames(seqtab0))

rownames(seqtab_filt) <- sub("_1\\.fq\\.gz$", "", rownames(seqtab_filt))

# make sure mapping is aligned
stopifnot(all(rownames(seqtab0) %in% sample_map$sampleName))
stopifnot(all(rownames(seqtab_filt) %in% sample_map$sampleName))

# create a named vector: sampleName -> clientId
name_map <- setNames(sample_map$clientId, sample_map$sampleName)

# rename seqtab rows
rownames(seqtab0) <- name_map[rownames(seqtab0)]
rownames(seqtab_filt) <- name_map[rownames(seqtab_filt)]

str(metadata$DNA_conc)
rownames(metadata) <- metadata$name_sequencing

ctrl_id <- setdiff(rownames(seqtab0), metadata$name_sequencing)

# How many reads are in the negative control?

sum(seqtab0[ctrl_id, ])
# 113155
sum(seqtab_filt[ctrl_id, ])
# 88543

## Create a row for the control in the metadata
ctrl_row <- metadata[1, , drop = FALSE]
ctrl_row[,] <- NA

ctrl_row$name_sequencing <- ctrl_id
ctrl_row$DNA_conc <- 0.001 ## dummy DNA concentration
ctrl_row$is_neg <- TRUE

metadata$is_neg <- FALSE

metadata2 <- rbind(metadata, ctrl_row)
rownames(metadata2) <- metadata2$name_sequencing


## Run decontam for full ASV table

# method = "frequency" - Contaminants are identified by frequency that varies inversely with sample DNA concentration
contam_freq <- isContaminant(
  seqtab0,
  method = "frequency",
  conc = metadata2$DNA_conc)

# method = "prevalence" - Contaminants are identified by increased prevalence in negative controls
contam_prev <- isContaminant(
  seqtab0,
  method = "prevalence",
  neg = metadata2$is_neg,
  threshold=0.55)

# method = "combined" - The frequency and prevalence probabilities are combined with Fisher's method and used to identify contaminants
contam <- isContaminant(
  seqtab0,
  method = "combined",
  conc = metadata2$DNA_conc,
  neg = metadata2$is_neg)



table(contam_freq$contaminant)
table(contam_prev$contaminant)
table(contam$contaminant)

## Run decontam for filtered ASV table

# method = "frequency" - Contaminants are identified by frequency that varies inversely with sample DNA concentration
contam_freq_filt <- isContaminant(
  seqtab_filt,
  method = "frequency",
  conc = metadata2$DNA_conc)

# method = "prevalence" - Contaminants are identified by increased prevalence in negative controls
contam_prev_filt <- isContaminant(
  seqtab_filt,
  method = "prevalence",
  neg = metadata2$is_neg,
  threshold=0.5)

# method = "combined" - The frequency and prevalence probabilities are combined with Fisher's method and used to identify contaminants
contam_filt <- isContaminant(
  seqtab_filt,
  method = "combined",
  conc = metadata2$DNA_conc,
  neg = metadata2$is_neg)


table(contam_freq_filt$contaminant)
table(contam_prev_filt$contaminant)
table(contam_filt$contaminant)

## Contaminant ASVs identity
##Step 1 — Extract contaminant ASV IDs
# full dataset
c_freq_full  <- rownames(contam_freq)[contam_freq$contaminant]
c_prev_full  <- rownames(contam_prev)[contam_prev$contaminant]
c_comb_full  <- rownames(contam)[contam$contaminant]

# no-singleton dataset
c_freq_filt  <- rownames(contam_freq_filt)[contam_freq_filt$contaminant]
c_prev_filt  <- rownames(contam_prev_filt)[contam_prev_filt$contaminant]
c_comb_filt  <- rownames(contam_filt)[contam_filt$contaminant]

# Overlap within each dataset (method comparison)
# Full dataset
length(intersect(c_freq_full, c_prev_full))
# 9
length(intersect(c_freq_full, c_comb_full))
# 70
length(intersect(c_prev_full, c_comb_full))
# 4

# Filtered dataset
length(intersect(c_freq_filt, c_prev_filt))
# 4
length(intersect(c_freq_filt, c_comb_filt))
# 91
length(intersect(c_prev_filt, c_comb_filt))
# 5

## Overlap between datasets (robustness to singleton removal)

# Frequency method robustness
length(intersect(c_freq_full, c_freq_filt))
# 119
length(intersect(c_prev_full, c_prev_filt))
# 32
length(intersect(c_comb_full, c_comb_filt))
# 26

## REMOVE CONTAMINANT samples with the results from the combined method for both datasets separately

seqtab.clean <- seqtab0[, !colnames(seqtab0) %in% c_comb_full, drop = FALSE]

seqtab.filt.clean <- seqtab_filt[, !colnames(seqtab_filt) %in% c_comb_filt, drop = FALSE]

# remove the control sample after filtering
seqtab.clean <- seqtab.clean[rownames(seqtab.clean) != ctrl_id, , drop = FALSE]
seqtab.filt.clean <- seqtab.filt.clean[rownames(seqtab.filt.clean) != ctrl_id, , drop = FALSE]

################################################################################
# Save results
################################################################################

decontam_output <- list(
  seqtab.full = seqtab.clean,
  seqtab.sensitivity = seqtab.filt.clean,
  metadata = metadata)

save(decontam_output, file = "decontam_output.Rdata")
