# Project name: root_fungi_epiphytic_terrestrial
# Author: Laura Méndez
# Date: 15.12.26

# test
#test2

# Goals:
# This script processes raw ITS2 amplicon sequencing data using the DADA2 pipeline
# following the ITS-specific workflow. The pipeline includes:
# 1. Quality filtering and primer removal
# 2. Error learning and denoising
# 3. Merging paired-end reads
# 4. Chimera removal
# 5. Post-processing filters (length, prevalence, abundance)

# Input files:
# 1. Raw sequencing files (*_1.fq.gz and *_2.fq.gz) in ../Sequencing/BMK_data/Data/rawdata
# 2. sampleName_clientId.txt - Sample name mapping


# Output files as dada_ouput.Rdata:
# 1. seqtab_full = seqtab_len - Initial ASV (all ASVs after general filtering with dada2)
# 2. seqtab_sensitivity = seqtab_prev - Filtered ASVs without singletons for sensitivity analysis
# 3. track_reads_ASV_filtering = track2 - Read counts through pipeline steps
# 4. sample_name - map of sample names


# References:
# - DADA2 ITS workflow: https://benjjneb.github.io/dada2/ITS_workflow.html
# - Callahan et al. (2016) Nature Methods 13:581-583

################## --------------------------------------------------------------------

# Load required libraries
library(dada2)       # packageVersion("dada2")   # v1.38.0
library(ShortRead)   # packageVersion("ShortRead") # v1.68.0
library(Biostrings)  # packageVersion("Biostrings") # v2.78.0
library(tidyverse)

# Define file paths to the raw data and primers --------------------------------

load("scripts_analyses/DADA2.Rdata")

# Set path to raw sequencing data
path <- "input_data/rawdata"
list.files(path)  # Verify files are present

# Get forward and reverse read files
fnFs <- sort(list.files(path, pattern = "_1.fq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_2.fq.gz", full.names = TRUE))

# Define ITS2 primers
FWD <- "GTGAATCATCGAATCTTTGAA"  # ITS86F (Turenne et al., 1999)
REV <- "TCCTCCGCTTATTGATATGC"   # ITS4R (White et al., 1990)

# Create primer orientation sequences -----------------------------------------

# Function to generate all possible orientations of a primer sequence
# Needed because primers can appear in any orientation in ITS amplicons
allOrients <- function(primer) {
  require(Biostrings)
  dna <- DNAString(primer)
  orients <- c(
    Forward = dna, 
    Complement = Biostrings::complement(dna), 
    Reverse = Biostrings::reverse(dna),
    RevComp = Biostrings::reverseComplement(dna)
  )
  return(sapply(orients, toString))
}

# Generate all orientations for both primers
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)

# Pre-filter: Remove reads with ambiguous bases (Ns) --------------------------

# Create directory for N-filtered files
fnFs.filtN <- file.path(path, "filtN", basename(fnFs))
fnRs.filtN <- file.path(path, "filtN", basename(fnRs))

# Filter out reads containing N bases
# maxN = 0 means no ambiguous bases allowed
filterAndTrim(
  fnFs, fnFs.filtN, 
  fnRs, fnRs.filtN, 
  maxN = 0, 
  multithread = TRUE
)

# Check for primer presence ----------------------------------------------------

# Function to count primer hits in a FASTQ file
primerHits <- function(primer, fn) {
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}

# Check primer presence in first sample (before removal)
rbind(
  FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), 
  FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]), 
  REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]), 
  REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]])
)

# Remove primers with Cutadapt -------------------------------------------------

# Set path to cutadapt executable (adjust for your system)
cutadapt <- "YOUR/OWN/PATH/Scripts/cutadapt.exe"
system2(cutadapt, args = "--version")  # Verify cutadapt is accessible

# Create output directory for primer-trimmed files
path.cut <- file.path(path, "cutadapt")
if (!dir.exists(path.cut)) dir.create(path.cut)

fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

# Get reverse-complement of primers for trimming
FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)

# Define cutadapt flags
# Trim FWD and reverse-complement of REV from R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC)
# Trim REV and reverse-complement of FWD from R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC)

# Run Cutadapt on all samples
# -n 2 means remove primers from both ends if present
for (i in seq_along(fnFs)) {
  system2(cutadapt, args = c(
    R1.flags, R2.flags, 
    "-n", 2,  # -n 2 required to remove primers from both ends
    "-o", fnFs.cut[i], "-p", fnRs.cut[i],  # Output files
    fnFs.filtN[i], fnRs.filtN[i]  # Input files
  ))
}

# Verify primers have been removed
rbind(
  FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[2]]), 
  FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[2]]), 
  REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[2]]), 
  REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[2]])
)

# Extract sample names ---------------------------------------------------------

# Get cutadapt-processed files
cutFs <- sort(list.files(path.cut, pattern = "_1.fq.gz", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "_2.fq.gz", full.names = TRUE))

# Read sample name mapping file
sample_map <- read.table(
  "input_data/sampleName_clientId.txt", 
  header = TRUE, 
  sep = "\t", 
  stringsAsFactors = FALSE, 
  comment.char = ""
)
colnames(sample_map)[1] <- "sampleName"

# Extract sample names from filenames and map to client IDs
sample_names_full <- basename(fnFs) %>%
  sub("(_good)?_[12]\\.fq\\.gz$", "", .)

sample_names <- data.frame(sampleName = sample_names_full) %>%
  left_join(sample_map, by = "sampleName") %>%
  mutate(clientId = coalesce(clientId, sampleName))  # Keep original if no match

# Store as sample.names for compatibility with downstream code
sample.names <- sample_names$clientId

# Inspect read quality profiles ------------------------------------------------

# Visualize quality profiles of first two samples
# This helps determine appropriate truncation parameters
plotQualityProfile(cutFs[1:2])
plotQualityProfile(cutRs[1:2])

# Filter and trim reads --------------------------------------------------------

# Define output file paths for filtered reads
filtFs <- file.path(path.cut, "filtered", basename(cutFs))
filtRs <- file.path(path.cut, "filtered", basename(cutRs))

# Apply quality filtering
# Parameters chosen based on ITS amplicon best practices
out <- filterAndTrim(
  cutFs, filtFs, 
  cutRs, filtRs, 
  maxN = 0,        # No ambiguous bases
  maxEE = c(2, 2), # Max 2 expected errors for forward and reverse
  truncQ = 2,      # Truncate reads at first base with quality ≤ 2
  minLen = 150,     # Minimum length after filtering
  rm.phix = TRUE,  # Remove PhiX control sequences
  compress = TRUE, 
  multithread = FALSE
)

head(out)  # Review filtering statistics

# Learn error rates ------------------------------------------------------------

# Learn error rates from filtered data
# DADA2 uses machine learning to build error models
errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)

# Visualize learned error rates
plotErrors(errF, nominalQ = TRUE)

# Dereplicate reads ------------------------------------------------------------

# Combine identical sequences to reduce computational burden
# Maintains count information for each unique sequence
derep.Fs <- derepFastq(filtFs, verbose = TRUE)
derep.Rs <- derepFastq(filtRs, verbose = TRUE)

# Sample inference (core DADA2 denoising algorithm) ---------------------------

# Infer true sequence variants using the DADA2 algorithm
# This is the core denoising step that distinguishes ASVs from errors
dadaFs <- dada(derep.Fs, err = errF, multithread = TRUE)
dadaRs <- dada(derep.Rs, err = errR, multithread = TRUE)

# Merge paired-end reads -------------------------------------------------------

# Merge forward and reverse reads
# Default: requires minimum 12 bp overlap without mismatches
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose = TRUE)

# Construct ASV table ----------------------------------------------------------

# Create sequence table (ASV × sample matrix)
seqtab <- makeSequenceTable(mergers)
dim(seqtab)  # Check dimensions

# Inspect distribution of sequence lengths
# ITS2 is naturally variable in length
table(nchar(getSequences(seqtab)))

# Remove chimeric sequences ----------------------------------------------------

# Identify and remove chimeras (PCR artifacts)
# method = "consensus" identifies sequences present in multiple samples
seqtab.nochim <- removeBimeraDenovo(
  seqtab, 
  method = "consensus", 
  multithread = TRUE, 
  verbose = TRUE
)

dim(seqtab.nochim)
sum(seqtab.nochim) / sum(seqtab)  # Proportion of reads retained

# Check sequence length distribution after chimera removal
table(nchar(getSequences(seqtab.nochim)))

# Track reads through pipeline -------------------------------------------------

# Function to extract read counts
getN <- function(x) sum(getUniques(x))

# Create tracking table showing read retention at each step
track <- cbind(
  out,  # Filter and trim output
  sapply(dadaFs, getN),  # After denoising forward
  sapply(dadaRs, getN),  # After denoising reverse
  sapply(mergers, getN),  # After merging
  rowSums(seqtab.nochim)  # After chimera removal
)

colnames(track) <- c(
  "input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim"
)
rownames(track) <- sample.names
head(track)

# OPTIONAL Post-processing filters ------------------------------------------------------

# we inspect the final data set by implementing different filters 

# Store initial unfiltered ASV table
seqtab0 <- seqtab.nochim

# Prevalence filter (remove singletons)
# Remove ASVs present in only 1 sample (likely sequencing noise)
prev <- colSums(seqtab_len > 0)
seqtab_prev <- seqtab_len[, prev >= 2, drop = FALSE]


c(n_ASVs_before = ncol(seqtab_len), n_ASVs_after = ncol(seqtab_prev))

# Abundance filter (remove extremely rare ASVs)
# Remove ASVs with fewer than 5 or 100 total reads across all samples
tot <- colSums(seqtab_len)
seqtab_filt <- seqtab_len[, tot >= 5, drop = FALSE]

seqtab_filt1 <- seqtab_len[, tot >= 100, drop = FALSE]

c(n_ASVs_before = ncol(seqtab_len), n_ASVs_after = ncol(seqtab_filt))

c(n_ASVs_before = ncol(seqtab_len), n_ASVs_after = ncol(seqtab_filt1))


# Per-sample relative abundance filter (0.01%)
threshold <- 0.0001  # 0.01%

# Convert to relative abundance per sample
seqtab_rel <- sweep(seqtab_len, 1, rowSums(seqtab_len), "/")

# Set low-abundance ASVs to zero per sample
seqtab_filt_001 <- seqtab_len
seqtab_filt_001[seqtab_rel < threshold] <- 0

# Remove ASVs that are zero across all samples
seqtab_filt_001 <- seqtab_filt_001[, colSums(seqtab_filt_001) > 0, drop = FALSE]

c(n_ASVs_before = ncol(seqtab_len), n_ASVs_after = ncol(seqtab_filt_001))

# Track filtering steps --------------------------------------------------------

# Helper functions for tracking
asv_richness <- function(tab) rowSums(tab > 0)  # ASVs per sample
reads_total <- function(tab) rowSums(tab)       # Total reads per sample

# Track ASV richness after each filter
track_filters_asv <- cbind(
  ASVs_nonchim = asv_richness(seqtab_len),
  ASVs_noSingletons = asv_richness(seqtab_prev),
  ASVs_rel_abund_001 = asv_richness(seqtab_filt_001))

# Track read counts after each filter
track_filters_reads <- cbind(
  Reads_noSingletons = reads_total(seqtab_prev),
  Reads_rel_abund_001 = reads_total(seqtab_filt_001)
)

# Combine all tracking information
track2 <- cbind(track, track_filters_asv, track_filters_reads)
head(track2)

# Save complete tracking table
#write.csv(track2, "reads_tracking.csv", row.names = TRUE)

################################################################################
# Save results
################################################################################
getwd()
# Save dada outputs
dada_output <- list(
  sample_map = sample_map,
  track_reads_ASV_filtering = track2,
  seqtab_full = seqtab_len,
  seqtab_sensitivity = seqtab_prev)

save(dada_output, file = "dada_output.Rdata")
