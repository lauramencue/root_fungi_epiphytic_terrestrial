# Project name: root_fungi_epiphytic_terrestrial
# Author: Laura Méndez
# Date: 15.02.2025

# Goals:
# Analyze beta diversity (community composition) to test:
# H2: Root-associated fungal community composition differs between rooting substrates, 
# with epiphytic ones being more distinct than terrestrial


# Analyses:
# - NMDS ordination
# - PERMANOVA (community differences)
# - Betadisper (dispersion homogeneity)
# - Venn Diagram (unique and shared ASVs per substrate and family

# Input: analyses.Rdata

# Setup ------------------------------------------------------------------------

rm(list = ls())

library(tidyverse)
library(vegan) # packageVersion("vegan") ‘2.7.3’
library(ggplot2)
library(OmicFlow) # packageVersion("OmicFlow") ‘1.5.0’
library(VennDiagram)

# Load data --------------------------------------------------------------------

load("analyses3.Rdata")
load("beta_diversity_results3.Rdata")

fungi <- analyses$fungi

metadata <- analyses$metadata

asv_tax_fungi <- analyses$asv_tax_fungi

# Calculate beta diversity ----------------------------------------------------

## Sample x ASV matrix (counts) ----
comm <- fungi %>%
  select(SampleID, ASV, Reads) %>%
  group_by(SampleID, ASV) %>%                  # just in case there are duplicates
  summarise(Reads = sum(Reads), .groups = "drop") %>%
  pivot_wider(names_from = ASV, values_from = Reads, values_fill = 0)
str(comm)
# Set rownames and keep numeric matrix
comm_mat <- comm %>%
  as.data.frame()
rownames(comm_mat) <- comm_mat$SampleID
comm_mat$SampleID <- NULL
comm_mat <- as.matrix(comm_mat)


# Make sure order matches comm_mat
meta1 <- metadata %>%
  filter(SampleID %in% rownames(comm_mat)) %>%
  arrange(match(SampleID, rownames(comm_mat)))

stopifnot(all(meta1$SampleID == rownames(comm_mat)))

str(comm_mat)
colnames(metadata)

## Hellinger transform + Bray–Curtis ----
comm_hel <- decostand(comm_mat, method = "hellinger")
bc_hel <- vegdist(comm_hel, method = "bray")


# NMDS ordination --------------------------------------------------------------

set.seed(123)

nmds <- metaMDS(comm_hel, distance = "bray", k = 2, trymax = 100)

# Extract scores and merge with metadata
nmds_df <- scores(nmds, display = "sites") %>%
  as.data.frame() %>%
  rownames_to_column("SampleID") %>%
  left_join(metadata, by = "SampleID")

head(nmds_df)

# PERMANOVA --------------------------------------------------------------------

# Test effects of species, growth form, and family
adonis2(bc_hel ~ host_species_identity, data = metadata, permutations = 999) 

adonis2(bc_hel ~ host_substrate, data = metadata, permutations = 999) 

adonis2(bc_hel ~ host_family, data = metadata, permutations = 999) 

adonis2(bc_hel ~ latitude + longitude, data = metadata, permutations = 999) 

adonis2(bc_hel ~ latitude + longitude + host_substrate + host_family + host_species_identity, by = "margin",
        data = metadata, permutations = 999) 


metadata$family_habitat <- interaction(metadata$host_family, metadata$host_substrate, sep = "_", drop = TRUE)
table(metadata$family_habitat)

pairwise_adonis(bc_hel, 
                groups = as.character(metadata$host_substrate),
                p.adjust.method = "bonferroni",
                perm = 999)

pairwise_adonis(bc_hel, 
                groups = as.character(metadata$host_family),
                p.adjust.method = "bonferroni",
                perm = 999) %>%
  arrange(desc(R2))

family_habitat_pairwise <- pairwise_adonis(bc_hel, 
                                           groups = as.character(metadata$family_habitat),
                                           p.adjust.method = "bonferroni",
                                           perm = 999)

family_habitat_pairwise %>%
  arrange(desc(R2))


# Test homogeneity of dispersions (betadisper) ---------------------------------

# Test by growth form
disp <- betadisper(bc_hel, metadata$host_substrate)

permutest(disp) 

# Test by family
disp1 <- betadisper(bc_hel, metadata$host_family)

permutest(disp1)

# Test family x habitat

disp2 <- betadisper(bc_hel, meta1$family_habitat)

permutest(disp2)

# Test by species
disp3 <- betadisper(bc_hel, metadata$host_species_identity)

permutest(disp3) 


# Venn Diagrams ----------------------------------------------------------------

asvs_by_substrate <- fungi %>%
  distinct(ASV, host_substrate) %>%
  split(.$host_substrate) %>%
  lapply(function(x) x$ASV)

asvs_by_family <- fungi %>% 
  distinct(ASV, host_family) %>% 
  split(.$host_family) %>% 
  lapply(function(x) x$ASV)

venn.diagram(
  x = asvs_by_substrate,
  fill = NA,
  cex = 1.2,
  cat.cex = 1.2,
  filename = NULL
) |> grid::grid.draw()

grid::grid.newpage()

venn.diagram(
  x = asvs_by_family,
  fill = NA,
  cex = 1.2,
  cat.cex = 1.2,
  filename = NULL
) |> grid::grid.draw()

# Save results------------------------------------------------------------------

# Save beta diversity results
beta_diversity_results <- list(
  nmds_df = nmds_df,
  family_habitat_pairwise = family_habitat_pairwise,
  asvs_by_substrate = asvs_by_substrate,
  asvs_by_family = asvs_by_family)

save(beta_diversity_results, "beta_diversity_results.Rdata")
