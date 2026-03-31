# Project name: root_fungi_epiphytic_terrestrial
# Author: Laura Méndez
# Date: [Date]

# Goals:
# SENSITIVITY ANALYSIS: Repeat key diversity analyses using OTUs instead of ASVs
# This tests whether results are robust to taxonomic unit definition
#
# Analyses performed:
# 1. Alpha diversity (OTU richness, Hill numbers, GLMM models)
# 2. Beta diversity (NMDS, PERMANOVA, betadisper, Venn diagrams)
# 3. Functional guild composition (proportions, PERMANOVA)
#


# Input files:
# 1. OTU_annotated.Rdata from script S2

# Output files:
# 1. OTU_diversity_results.Rdata - All sensitivity analysis results

################################################################################
# Setup ------------------------------------------------------------------------
################################################################################

rm(list = ls())

library(tidyverse)
library(vegan)
library(hillR)
library(glmmTMB)
library(performance)
library(car)
library(VennDiagram)

load("OTU_annotated.Rdata")

################################################################################
# Extract OTU data -------------------------------------------------------------
################################################################################

fungi_otu <- OTU_annotated$otu_data
metadata_otu <- fungi_otu %>% distinct(SampleID, host_species_identity, host_substrate, host_family, host_wcvp_name)

################################################################################
# PART 1: ALPHA DIVERSITY (OTU RICHNESS) --------------------------------------
################################################################################

# OTU richness per sample
alpha_otu <- fungi_otu %>%
  group_by(SampleID, host_species_identity, host_substrate, host_family, host_wcvp_name) %>%
  summarise(OTU_richness = n_distinct(OTU_ID), .groups = "drop")

# Summary by growth form
alpha_otu %>%
  group_by(host_substrate) %>%
  summarise(
    mean = mean(OTU_richness),
    sd = sd(OTU_richness),
    n = n(),
    .groups = "drop"
  )

# Summary by family
alpha_otu %>%
  group_by(host_family) %>%
  summarise(
    mean = mean(OTU_richness),
    sd = sd(OTU_richness),
    n = n(),
    .groups = "drop"
  ) %>%
  arrange(desc(mean))

# Summary by family × substrate
alpha_otu %>%
  group_by(host_family, host_substrate) %>%
  summarise(
    mean = mean(OTU_richness),
    sd = sd(OTU_richness),
    n = n(),
    .groups = "drop"
  ) %>%
  arrange(desc(mean))

################################################################################
# Hill diversity (q = 0, 1, 2) -------------------------------------------------
################################################################################

# Create OTU × sample matrix for Hill numbers
comm_otu <- fungi_otu %>%
  select(SampleID, OTU_ID, Reads) %>%
  pivot_wider(names_from = OTU_ID, values_from = Reads, values_fill = 0)

comm_mat_otu <- comm_otu %>%
  as.data.frame() %>%
  column_to_rownames("SampleID") %>%
  as.matrix()

# Calculate Hill numbers (q=0: richness, q=1: Shannon, q=2: Simpson)
hill_otu <- hill_taxa(comm_mat_otu, q = c(0, 1, 2))

hill_otu_df <- hill_otu %>%
  as.data.frame() %>%
  rownames_to_column("SampleID") %>%
  pivot_longer(cols = -SampleID, names_to = "q_order", values_to = "Hill_number") %>%
  left_join(metadata_otu, by = "SampleID")

# Summary by substrate and Hill order
hill_otu_df %>%
  group_by(host_substrate, q_order) %>%
  summarise(
    mean = mean(Hill_number),
    sd = sd(Hill_number),
    .groups = "drop"
  )

################################################################################
# Statistical models - GLMM ----------------------------------------------------
################################################################################

# Check for overdispersion
alpha_otu %>%
  summarise(
    mean_rich = mean(OTU_richness),
    var_rich = var(OTU_richness),
    overdispersion = var_rich / mean_rich
  )

# Negative binomial GLMM with host species as random effect
m_nb_otu <- glmmTMB(
  OTU_richness ~ host_family + host_substrate + host_family:host_substrate + 
    (1 | host_wcvp_name),
  family = nbinom2,
  data = alpha_otu
)

summary(m_nb_otu)

# Model without random effect (for comparison)
m_nb_otu_norand <- glmmTMB(
  OTU_richness ~ host_family + host_substrate + host_family:host_substrate,
  family = nbinom2,
  data = alpha_otu
)

# Likelihood ratio test
anova(m_nb_otu_norand, m_nb_otu)

# Model performance
r2(m_nb_otu)
icc(m_nb_otu)

# Type II Wald chi-square tests
Anova(m_nb_otu, type = "II")

################################################################################
# PART 2: BETA DIVERSITY (COMMUNITY COMPOSITION) ------------------------------
################################################################################

# Ensure metadata matches community matrix
meta_otu <- metadata_otu %>%
  filter(SampleID %in% rownames(comm_mat_otu)) %>%
  arrange(match(SampleID, rownames(comm_mat_otu)))

stopifnot(all(meta_otu$SampleID == rownames(comm_mat_otu)))

# Hellinger transformation + Bray-Curtis dissimilarity
comm_hel_otu <- decostand(comm_mat_otu, method = "hellinger")
bc_hel_otu <- vegdist(comm_hel_otu, method = "bray")

################################################################################
# NMDS ordination --------------------------------------------------------------
################################################################################

set.seed(123)
nmds_otu <- metaMDS(comm_hel_otu, distance = "bray", k = 2, trymax = 100)

# Extract scores
nmds_otu_df <- scores(nmds_otu, display = "sites") %>%
  as.data.frame() %>%
  rownames_to_column("SampleID") %>%
  left_join(metadata_otu, by = "SampleID")

################################################################################
# PERMANOVA --------------------------------------------------------------------
################################################################################

# Test effects on OTU composition
perm_species_otu <- adonis2(bc_hel_otu ~ host_species_identity, 
                            data = meta_otu, permutations = 999)

perm_substrate_otu <- adonis2(bc_hel_otu ~ host_substrate, 
                              data = meta_otu, permutations = 999)

perm_family_otu <- adonis2(bc_hel_otu ~ host_family, 
                           data = meta_otu, permutations = 999)

perm_full_otu <- adonis2(bc_hel_otu ~ host_substrate + host_family + host_species_identity, 
                         by = "margin", data = meta_otu, permutations = 999)

################################################################################
# Betadisper (homogeneity of dispersion) ---------------------------------------
################################################################################

# Test by substrate
disp_substrate_otu <- betadisper(bc_hel_otu, meta_otu$host_substrate)
permutest(disp_substrate_otu)

# Test by family
disp_family_otu <- betadisper(bc_hel_otu, meta_otu$host_family)
permutest(disp_family_otu)

# Test by species
disp_species_otu <- betadisper(bc_hel_otu, meta_otu$host_species_identity)
permutest(disp_species_otu)

################################################################################
# Venn diagrams (OTU overlap) --------------------------------------------------
################################################################################

# OTUs by substrate
otus_by_substrate <- fungi_otu %>%
  distinct(OTU_ID, host_substrate) %>%
  split(.$host_substrate) %>%
  lapply(function(x) x$OTU_ID)

# OTUs by family
otus_by_family <- fungi_otu %>%
  distinct(OTU_ID, host_family) %>%
  split(.$host_family) %>%
  lapply(function(x) x$OTU_ID)

# Calculate overlaps
overlap_substrate <- length(intersect(otus_by_substrate$E, otus_by_substrate$T))
unique_epi <- length(setdiff(otus_by_substrate$E, otus_by_substrate$T))
unique_ter <- length(setdiff(otus_by_substrate$T, otus_by_substrate$E))

################################################################################
# PART 3: FUNCTIONAL GUILD COMPOSITION ----------------------------------------
################################################################################

# Filter for root-associated guilds
fungi_otu_roots <- fungi_otu %>%
  filter(!FunctionalGuilds %in% c(
    "wood_litter_soil_saprotroph",
    "plant_pathogen",
    "foliar_endophyte",
    "animal-associated",
    "other",
    NA
  ))

# Guild proportions by substrate
guild_otu_substrate <- fungi_otu_roots %>%
  distinct(SampleID, OTU_ID, FunctionalGuilds, host_substrate) %>%
  count(host_substrate, FunctionalGuilds, name = "OTUs") %>%
  group_by(host_substrate) %>%
  mutate(prop = OTUs / sum(OTUs)) %>%
  ungroup()

# Guild proportions by family
guild_otu_family <- fungi_otu_roots %>%
  distinct(SampleID, OTU_ID, FunctionalGuilds, host_family) %>%
  count(host_family, FunctionalGuilds, name = "OTUs") %>%
  group_by(host_family) %>%
  mutate(prop = OTUs / sum(OTUs)) %>%
  ungroup()

# Guild proportions by family × substrate
guild_otu_family_substrate <- fungi_otu_roots %>%
  distinct(OTU_ID, host_family, host_substrate, FunctionalGuilds) %>%
  count(host_family, host_substrate, FunctionalGuilds, name = "n_OTUs") %>%
  group_by(host_family, host_substrate) %>%
  mutate(prop = n_OTUs / sum(n_OTUs)) %>%
  ungroup()

################################################################################
# PERMANOVA on guild composition -----------------------------------------------
################################################################################

# Create OTU × sample matrix for root-associated fungi
comm_roots_otu <- fungi_otu_roots %>%
  select(SampleID, OTU_ID, Reads) %>%
  pivot_wider(names_from = OTU_ID, values_from = Reads, values_fill = 0) %>%
  as.data.frame() %>%
  column_to_rownames("SampleID") %>%
  as.matrix()

# Match metadata
meta_roots_otu <- metadata_otu %>%
  filter(SampleID %in% rownames(comm_roots_otu)) %>%
  arrange(match(SampleID, rownames(comm_roots_otu)))

stopifnot(all(meta_roots_otu$SampleID == rownames(comm_roots_otu)))

# Create guild-aggregated matrix
otu_guilds <- fungi_otu_roots %>%
  select(OTU_ID, FunctionalGuilds) %>%
  distinct() %>%
  filter(!is.na(FunctionalGuilds))

guild_mat_otu <- comm_roots_otu %>%
  as.data.frame() %>%
  rownames_to_column("SampleID") %>%
  pivot_longer(-SampleID, names_to = "OTU_ID", values_to = "reads") %>%
  filter(reads > 0) %>%
  left_join(otu_guilds, by = "OTU_ID") %>%
  filter(!is.na(FunctionalGuilds)) %>%
  group_by(SampleID, FunctionalGuilds) %>%
  summarise(richness = n_distinct(OTU_ID), .groups = "drop") %>%
  pivot_wider(names_from = FunctionalGuilds, values_from = richness, values_fill = 0) %>%
  column_to_rownames("SampleID") %>%
  as.matrix()

# Distance matrix
guild_dist_otu <- vegdist(guild_mat_otu, method = "jaccard")

# PERMANOVA on guild composition
perm_guild_substrate_otu <- adonis2(guild_dist_otu ~ host_substrate, 
                                    data = meta_roots_otu, permutations = 999)

perm_guild_family_otu <- adonis2(guild_dist_otu ~ host_family, 
                                 data = meta_roots_otu, permutations = 999)

perm_guild_full_otu <- adonis2(guild_dist_otu ~ host_substrate * host_family, 
                               data = meta_roots_otu, permutations = 999)

################################################################################
# Save all OTU sensitivity results ---------------------------------------------
################################################################################

OTU_diversity_results <- list(
  # Alpha diversity
  alpha_otu = alpha_otu,
  hill_otu = hill_otu_df,
  model_nb_otu = m_nb_otu,
  model_nb_otu_norand = m_nb_otu_norand,
  model_comparison = anova(m_nb_otu_norand, m_nb_otu),
  model_performance = list(r2 = r2(m_nb_otu), icc = icc(m_nb_otu)),
  model_anova = Anova(m_nb_otu, type = "II"),
  
  # Beta diversity
  nmds_otu = nmds_otu,
  nmds_otu_scores = nmds_otu_df,
  distance_matrix = bc_hel_otu,
  permanova_species = perm_species_otu,
  permanova_substrate = perm_substrate_otu,
  permanova_family = perm_family_otu,
  permanova_full = perm_full_otu,
  betadisper_substrate = disp_substrate_otu,
  betadisper_family = disp_family_otu,
  betadisper_species = disp_species_otu,
  
  # Venn diagram data
  otus_by_substrate = otus_by_substrate,
  otus_by_family = otus_by_family,
  overlap_summary = list(
    shared = overlap_substrate,
    epiphytic_unique = unique_epi,
    terrestrial_unique = unique_ter
  ),
  
  # Guild composition
  guild_substrate = guild_otu_substrate,
  guild_family = guild_otu_family,
  guild_family_substrate = guild_otu_family_substrate,
  permanova_guild_substrate = perm_guild_substrate_otu,
  permanova_guild_family = perm_guild_family_otu,
  permanova_guild_full = perm_guild_full_otu
)

save(OTU_diversity_results, file = "OTU_diversity_results.Rdata")

