# Project name: root_fungi_epiphytic_terrestrial
# Author: Laura Méndez
# Date: 15.02.26

# Goals: Figure SXX. Symbiotic guild composition by rooting host_substrate within 
# host families

# Input: analyses.Rdata

# Setup ------------------------------------------------------------------------

rm(list = ls())

library(tidyverse)
library(ggplot2)
library(ggdist)
library(forcats)
library(patchwork)

# Load data --------------------------------------------------------------------

load("analyses.Rdata")

fungi <- analyses$fungi
asv_tax_fungi <- analyses$asv_tax_fungi

## Font size across panels

base_axis_size   <- 15
base_legend_size <- 15
base_title_size  <- 15

## Prepare plotting data -------------------------------------------------------
guild_prop_fam_hab <- fungi %>%
  filter(!is.na(BroadSymbiontGroup)) %>%
  mutate(host_substrate = factor(host_substrate, 
    levels = c("Epiphytic", "Terrestrial")),
    BroadSymbiontGroup = factor(BroadSymbiontGroup, levels = broad_symbionts),
    host_family = factor(host_family, 
    levels = c("Araceae","Bromeliaceae","Orchidaceae","Piperaceae"))) %>%
  distinct(ASV, host_family, host_substrate, BroadSymbiontGroup) %>%  
  count(host_family, host_substrate, BroadSymbiontGroup, name = "n") %>%
  group_by(host_family, host_substrate) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

p_guild_family_host_substrate <- ggplot(guild_prop_fam_hab,
  aes(x = host_substrate, y = prop, fill = BroadSymbiontGroup)) +
  geom_col(width = 0.7) +
  facet_wrap(~ host_family, ncol = 2) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = pal_fg, drop = FALSE) +
  labs(y = "Proportion of symbiotic ASVs",
    x = NULL,
    fill = NULL) +
  theme_classic(base_size = 15) +
  theme(
    axis.text.x = element_text(size = base_axis_size),
    axis.text.y = element_text(size = base_axis_size),
    axis.title.x = element_text(size = base_title_size),
    legend.position = "right",
    strip.text = element_text(size = base_axis_size, face = "bold"),
    plot.title = element_text(size = base_title_size, face = "bold"),
    plot.title.position = "plot")

p_guild_family_host_substrate

# Save Figure ------------------------------------------------------------------

ggsave(
  filename = "../Figures/FigureS4.barplots_guild_family.png",
  plot = p_guild_family_habitat,
  width = 8,
  height = 8,
  dpi = 600,
  units = "in",
  bg = "white"
)