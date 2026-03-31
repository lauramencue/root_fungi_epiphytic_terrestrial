# Project name: root_fungi_epiphytic_terrestrial
# Author: Laura Méndez
# Date: 15.02.2025
# Last modified: [Date]

# Goals:
# Perform alpha diversity analyses (ASV richness) to test the hypothesis:
# H1: Terrestrial plants have higher ASV richness than epiphytic plants and host identity modulates this effect

# Input files:
# 1. analyses.Rdata - Prepared data from script 05

# Output:
# alpha_diversity_results.Rdata

# Key analyses:
# - Summary statistics by growth form, family, species
# - Hill diversity calculation
# - Negative binomial generalized linear mixed models


# Setup ------------------------------------------------------------------------

rm(list = ls())
setwd("C:/Users/mendezcu_lu/Dropbox/Work/03_Marburg/02_epiphites_terrestrial_root_funghi/R/scripts")

library(tidyverse)
library(hillR) # packageVersion("hillR") ‘0.5.2’
library(glmmTMB) # packageVersion("glmmTMB") ‘1.1.13’
library(performance)  # Model diagnostics and R² packageVersion("performance") ‘0.15.3’
library(car)

# Load data --------------------------------------------------------------------

load("analyses.Rdata")


fungi <- analyses$fungi

metadata <- analyses$metadata

asv_tax_fungi <- analyses$asv_tax_fungi

# Calculate alpha diversity ----------------------------------------------------

# ASV richness per sample
alpha_df <- fungi %>%
  group_by(SampleID, host_species_identity, host_substrate, host_family, host_wcvp_name) %>%
  summarise(ASV_richness = n_distinct(ASV), .groups = "drop")

# Summary statistics by growth form

alpha_df %>%
  group_by(host_substrate) %>%
  summarise(
    mean = mean(ASV_richness),
    sd = sd(ASV_richness),
    n = n(),
    .groups = "drop"
  ) %>% print()

# Summary by host family

alpha_df %>%
  group_by(host_family) %>%
  summarise(
    mean = mean(ASV_richness),
    sd = sd(ASV_richness),
    n = n(),
    .groups = "drop"
  ) %>% 
  arrange(desc(mean)) %>% print()

# Summary by family × growth form

alpha_df %>%
  group_by(host_family, host_substrate) %>%
  summarise(
    mean = mean(ASV_richness),
    sd = sd(ASV_richness),
    n = n(),
    .groups = "drop"
  ) %>% 
  arrange(desc(mean))

alpha_df_species <- alpha_df %>%
  group_by(SampleID, host_species_identity, host_family, host_substrate) %>%
  summarise(
    mean = mean(ASV_richness),
    sd = sd(ASV_richness),
    n = n(),
    .groups = "drop"
  ) %>% 
  arrange(desc(mean)) 

alpha_df_species


# Hill numbers -----------------------------------------------------------------

# Build Sample x ASV matrix (counts)
comm_mat <- fungi %>%
  group_by(SampleID, ASV) %>%
  summarise(Reads = sum(Reads), .groups = "drop") %>%
  pivot_wider(names_from = ASV, values_from = Reads, values_fill = 0) %>%
  arrange(SampleID)

# keep sample IDs as rownames and drop the SampleID column for hillR
comm <- as.data.frame(comm_mat)
rownames(comm) <- comm$SampleID
comm$SampleID <- NULL


# Hill numbers (alpha diversity) using hillR
# q = 0 (richness), q = 1 (exp Shannon), q = 2 (inv Simpson)
hill_q0 <- hill_taxa(comm, q = 0)
hill_q1 <- hill_taxa(comm, q = 1)
hill_q2 <- hill_taxa(comm, q = 2)

# hill_taxa returns a numeric vector named by sample (rownames(comm))
hill_df <- tibble(
  SampleID = names(hill_q0),
  hill_q0  = as.numeric(hill_q0),
  hill_q1  = as.numeric(hill_q1),
  hill_q2  = as.numeric(hill_q2)
)

# Add metadata back (one row per sample)
meta_df <- fungi %>%
  distinct(SampleID, host_species_identity, host_substrate, host_family, host_wcvp_name)

hill_df <- hill_df %>%
  left_join(meta_df, by = "SampleID")

## Mean per substrate type
hill_df %>%
  group_by(host_substrate) %>%
  summarise(
    n = n(),
    across(
      c(hill_q0, hill_q1, hill_q2),
      list(mean = ~mean(.x, na.rm = TRUE),
           sd   = ~sd(.x, na.rm = TRUE)),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  )

## Mean per family
hill_df %>%
  group_by(host_family) %>%
  summarise(
    n = n(),
    across(
      c(hill_q0, hill_q1, hill_q2),
      list(mean = ~mean(.x, na.rm = TRUE),
           sd   = ~sd(.x, na.rm = TRUE)),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  )

# quick sanity checks
hill_df %>% summarise(across(starts_with("hill_"), list(min=min, mean=mean, max=max)))

hist(hill_df$hill_q0)
boxplot(hill_df$hill_q0)
hist(hill_df$hill_q2)
boxplot(alpha_df$ASV_richness)

# round down
hill_df <- hill_df %>%
  mutate(
    f_hill_q1 = floor(hill_q1), ## round down - floor to have natural numbers 
    f_hill_q2 = floor(hill_q2)
  )

summary(hill_df)

# Negative binomial GLMMs ------------------------------------------------------

# Check for overdispersion

alpha_df %>%
  summarise(
    mean_rich = mean(ASV_richness),
    var_rich = var(ASV_richness),
    overdispersion_ratio = var_rich / mean_rich
  ) %>% print()
# Overdispersion ratio >> 1 indicates negative binomial is appropriate


# richness (q=0) 

m_q0 <- glmmTMB(
  hill_q0 ~ host_family + host_substrate + (1 | host_species_identity),
  family = nbinom2,
  data = hill_df
)

m_q0_nospp <- glmmTMB(
  hill_q0 ~ host_family + host_substrate,
  family = nbinom2,
  data = hill_df
)

anova(m_q0_nospp, m_q0)
# model with random effect is better

summary(m_q0)

# Check model
check_model(m_q0)
check_overdispersion(m_q0)
Anova(m_q0)


# q1

m_q1 <- glmmTMB(
  f_hill_q1 ~ host_family + host_substrate + (1 | host_species_identity),
  family = nbinom2,
  data = hill_df
)

m_q1_nospp <- glmmTMB(
  f_hill_q1 ~ host_family + host_substrate,
  family = nbinom2,
  data = hill_df
)

anova(m_q1_nospp, m_q1)
# model with random effect is better

summary(m_q1)

# Check model
check_model(m_q1)
check_overdispersion(m_q1)
Anova(m_q1)

# q2

m_q2 <- glmmTMB(
  f_hill_q2 ~ host_family + host_substrate + (1 | host_species_identity),
  family = nbinom2,
  data = hill_df
)

m_q2_nospp <- glmmTMB(
  f_hill_q2 ~ host_family + host_substrate,
  family = nbinom2,
  data = hill_df
)

anova(m_q2_nospp, m_q2)
# model with random effect is better

summary(m_q2)

# Check model
check_model(m_q2)
check_overdispersion(m_q2)
Anova(m_q2)


r2(m_q0); icc(m_q0)
r2(m_q1); icc(m_q1)
r2(m_q2); icc(m_q2)

# Save results -----------------------------------------------------------------


# Save alpha diversity data and models
alpha_diversity_results <- list(
  alpha_df = alpha_df,
  comm = comm,
  m_q0 = m_q0,
  m_q1 = m_q1,
  m_q2 = m_q2)

save(alpha_diversity_results, file = "alpha_diversity_results.Rdata")

