# Project name: root_fungi_epiphytic_terrestrial
# Author: Laura Méndez
# Date: [Date]

# Goals:
# Analyze functional guild composition across rooting substrates and host families
# to test: H3 - Terrestrial individuals will host higher proportions of symbiotic 
# guilds, specifically AMF and other mycorrhizal types, while epiphytes will host 
# higher proportions of dark septate endophytes and other root endophytes

# Analyses:
# - Guild proportions by growth form and family
# - Compositional differences using PERMANOVA on guild-aggregated data

# Input: analyses.Rdata
# Dependencies: tidyverse, vegan

# Setup ------------------------------------------------------------------------

rm(list = ls())

library(tidyverse)
library(vegan)
library(glmmTMB) # packageVersion("glmmTMB") ‘1.1.13’
library(performance)  # Model diagnostics and R² packageVersion("performance") ‘0.15.3’
library(car)

# Load data --------------------------------------------------------------------

load("analyses.Rdata")


fungi <- analyses$fungi

metadata <- analyses$metadata

asv_tax_fungi <- analyses$asv_tax_fungi

# Calculate guild proportions --------------------------------------------------

func_df <- fungi %>%
  filter(!is.na(BroadSymbiontGroup)) %>%
  distinct(host_species_identity, ASV, BroadSymbiontGroup, host_substrate) %>%
  count(host_substrate, BroadSymbiontGroup, name = "ASVs") %>%
  group_by(host_substrate) %>%
  mutate(prop = ASVs / sum(ASVs)) %>%
  ungroup()

func_df %>%
  select(host_substrate, BroadSymbiontGroup, prop) %>%
  arrange(host_substrate, desc(prop)) %>%
  print(n = 20)

### Of only symbionts

symb_prop_habitat <- fungi %>%
  filter(!is.na(FunctionalGuilds)) %>%
  distinct(host_substrate, ASV, FunctionalGuilds, BroadSymbiontGroup, is_symbiotic) %>%
  group_by(host_substrate) %>%
  summarise(
    total_asv = n(),
    symb_asv  = sum(is_symbiotic),
    prop_symb = symb_asv / total_asv,
    .groups = "drop"
  )

symb_prop_habitat

symb_family_hab_prop <- fungi %>%
  filter(!is.na(FunctionalGuilds)) %>%
  distinct(host_family, host_substrate, ASV, BroadSymbiontGroup, is_symbiotic) %>%
  group_by(host_family, host_substrate) %>%
  summarise(
    total_asv = n(),
    symb_asv  = sum(is_symbiotic),
    prop_symb = symb_asv / total_asv,
    .groups = "drop"
  )

symb_family_hab_prop
str(fungi)
symb_sample_prop <- fungi %>%
  filter(!is.na(FunctionalGuilds)) %>%
  distinct(SampleID, ASV, BroadSymbiontGroup, .keep_all = T) %>%
  group_by(SampleID, host_substrate, host_family, host_species_identity) %>%
  summarise(
    total_asv = n(),
    symb_asv  = sum(is_symbiotic),
    prop_symb = symb_asv / total_asv,
    .groups = "drop"
  )

symb_sample_prop
colnames(symb_sample_prop)

# Binomial GLMMs, proportion of symbiotc guilds---------------------------------

## Model without interaction
?glmmTMB

m_symb <- glmmTMB(prop_symb ~ 
                  host_family + host_substrate + (1 | host_species_identity),
                family = binomial,
                weights = total_asv,
                data = symb_sample_prop)


summary(m_symb)
check_model(m_symb)
check_overdispersion(m_symb)

r2(m_symb); icc(m_symb)

# Model without random effect (for comparison)
m_symb_nospp <- glmmTMB(prop_symb ~ 
                      host_family + host_substrate,
                      family = binomial,
                      weights = total_asv,
                      data = symb_sample_prop)

summary(m_symb_nospp)

# Compare models (test significance of random effect)

anova(m_symb_nospp, m_symb)

# Functional guild composition (proportion of each symbiotic guild) ------------

broad_symbionts<- rev(c("dark_septate_endophytes",
                        "other_root_endophytes",
                        "arbuscular_mycorrhizal",
                        "other_mycorrhizal_types"))

symb_ind_guild <- fungi %>%
  filter(!is.na(BroadSymbiontGroup)) %>%
  count(SampleID, host_substrate, host_family, host_species_identity, BroadSymbiontGroup, name = "n_asv") %>%
  group_by(SampleID, host_substrate, host_family, host_species_identity) %>%
  mutate(total_symb_asv = sum(n_asv)) %>%
  ungroup() %>%
  mutate(
    host_substrate = factor(host_substrate),
    host_family = factor(host_family),
    host_species_identity = factor(host_species_identity),
    BroadSymbiontGroup = factor(BroadSymbiontGroup, levels = broad_symbionts)
  )

fit_nb_add <- function(g) {
  
  # base: all individuals (with symbiotic fungi) + metadata + total_symb_asv
  base_dat <- symb_ind_guild %>%
    distinct(SampleID, host_substrate, host_family, host_species_identity, total_symb_asv)
  
  # guild-specific counts (only individuals where the guild is present)
  g_dat <- symb_ind_guild %>%
    filter(BroadSymbiontGroup == g) %>%
    distinct(SampleID, n_asv)
  
  # join and fill absences with 0
  dat <- base_dat %>%
    left_join(g_dat, by = "SampleID") %>%
    mutate(
      n_asv = ifelse(is.na(n_asv), 0, n_asv),
      log_total = log(total_symb_asv),
      host_substrate = factor(host_substrate),
      host_family = factor(host_family),
      host_species_identity = factor(host_species_identity)
    )
  
  glmmTMB(
    n_asv ~ host_substrate + host_family + offset(log_total) + (1|host_species_identity),
    family = nbinom2,
    data = dat
  )
}


fit_nb_nospp <- function(g) {
  
  # base: all individuals (with symbiotic fungi) + metadata + total_symb_asv
  base_dat <- symb_ind_guild %>%
    distinct(SampleID, host_substrate, host_family, host_species_identity, total_symb_asv)
  
  # guild-specific counts (only individuals where the guild is present)
  g_dat <- symb_ind_guild %>%
    filter(BroadSymbiontGroup == g) %>%
    distinct(SampleID, n_asv)
  
  # join and fill absences with 0
  dat <- base_dat %>%
    left_join(g_dat, by = "SampleID") %>%
    mutate(
      n_asv = ifelse(is.na(n_asv), 0, n_asv),
      log_total = log(total_symb_asv),
      host_substrate = factor(host_substrate),
      host_family = factor(host_family),
      host_species_identity = factor(host_species_identity)
    )
  
  glmmTMB(
    n_asv ~ host_substrate + host_family + offset(log_total),
    family = nbinom2,
    data = dat
  )
}




## Fitting models

fits_nb <- setNames(lapply(broad_symbionts, fit_nb_add), broad_symbionts)

fits_nb_nospp <- setNames(lapply(broad_symbionts, fit_nb_nospp), broad_symbionts)

anova(fits_nb[["arbuscular_mycorrhizal"]], fits_nb_nospp[["arbuscular_mycorrhizal"]])
anova(fits_nb[["dark_septate_endophytes"]], fits_nb_nospp[["dark_septate_endophytes"]])
anova(fits_nb[["other_root_endophytes"]], fits_nb_nospp[["other_root_endophytes"]])
anova(fits_nb[["other_mycorrhizal_types"]], fits_nb_nospp[["other_mycorrhizal_types"]])

summary(fits_nb_nospp[["arbuscular_mycorrhizal"]])


summary(fits_nb[["arbuscular_mycorrhizal"]])
summary(fits_nb[["dark_septate_endophytes"]])
summary(fits_nb[["other_root_endophytes"]])
summary(fits_nb[["other_mycorrhizal_types"]])



r2(fits_nb[["arbuscular_mycorrhizal"]]); icc(fits_nb[["arbuscular_mycorrhizal"]])
r2(fits_nb[["dark_septate_endophytes"]]); icc(fits_nb[["dark_septate_endophytes"]])
r2(fits_nb[["other_root_endophytes"]]); icc(fits_nb[["other_root_endophytes"]])
r2(fits_nb[["other_mycorrhizal_types"]]); icc(fits_nb[["other_mycorrhizal_types"]])

# Save results------------------------------------------------------------------

# Save beta diversity results
guild_composition_results <- list(
  m_symb = m_symb,
  fits_nb = fits_nb)

save(guild_composition_results, "guild_composition_results.Rdata")

