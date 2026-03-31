# Project name: root_fungi_epiphytic_terrestrial
# Author: Laura Méndez
# Date: 15.02.26

# Goals: Figure 4. Functional guild availability and composition of 
# root-associated fungal communities

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

# FIRST PANEL: With functional guild vs without ----------------------------

## Prepare plotting data

# Total ASVs (all rows = all ASVs in your ASV taxonomy table)
total_asvs <- nrow(asv_tax_fungi)

# How many have FunctionalGuilds assigned?
n_ident_fg <- asv_tax_fungi %>%
  summarise(n = sum(!is.na(FunctionalGuilds) & FunctionalGuilds != "")) %>%
  pull(n)

n_unident_fg <- total_asvs - n_ident_fg

df_id_fg <- tibble(Status = factor(c("Functional guild available", 
                                     "No functional guild available"),
    levels = c("Functional guild available", "No functional guild available")),
  n = c(n_ident_fg, n_unident_fg))

# Plotting

p1_fg <- ggplot(df_id_fg, aes(x = 1, y = n, fill = Status)) +
  geom_col(width = 0.75) +
  coord_flip() +
  scale_fill_manual(values = c(
    "Functional guild available" = "grey30",
    "No functional guild available" = "grey80")) +
  scale_y_continuous(labels = scales::label_comma(),
  expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL,
    y = "Number of ASVs",
    fill = NULL) +
  theme_classic(base_size = 15) +
  theme(axis.text.x = element_text(size = base_axis_size),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.x = element_text(size = base_title_size,
    margin = margin(t = 10)),
    legend.text = element_text(size = base_legend_size),
    legend.title = element_blank(),
    legend.position = "top",
    plot.title = element_text(size = base_title_size, face = "bold"),
    plot.title.position = "plot")

p1_fg

# SECOND PANEL: All different functional guilds ----------------------------

# Prepare plotting data

guild_levels <- rev(c(
  "wood_litter_soil_saprotroph",
  "plant_pathogen",
  "foliar_endophyte",
  "animal-associated",
  "other",
  "other_mycorrhizal_types",
  "dark_septate_endophytes",
  "other_root_endophytes",
  "arbuscular_mycorrhizal"))

pal_fg <- c("wood_litter_soil_saprotroph" = "#8DD3C7",
  "plant_pathogen" = "#FFFFB3",
  "foliar_endophyte" = "#BEBADA",
  "animal-associated" = "#80B1D3",
  "other" = "#FCCDE5",
  "dark_septate_endophytes" = "#FB8072",
  "other_root_endophytes" = "#B3DE69",
  "other_mycorrhizal_types" = "#FDB462",
  "arbuscular_mycorrhizal" = "#BC80BD")

guild <- asv_tax_fungi %>%
  filter(!is.na(all_functional) & all_functional != "") %>%  # only identified
  count(all_functional, name = "n") %>%
  #mutate(
  # FunctionalGuilds = if_else(n < 100, "other", FunctionalGuilds)
  #) %>%
  group_by(all_functional) %>%
  summarise(n = sum(n), .groups = "drop") %>%
  mutate(all_functional = fct_reorder(all_functional, n, .desc = TRUE))

guild$all_functional <- factor(
  guild$all_functional,
  levels = guild_levels
)

# Plotting

p2_fg <- ggplot(guild, aes(x = 1, y = n, fill = all_functional)) +
  geom_col(width = 0.75) +
  coord_flip() +
  scale_fill_manual(values = pal_fg, drop = FALSE) +
  scale_y_continuous(labels = scales::label_comma(),
  expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "Number of ASVs", fill = NULL) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x = element_text(size = base_axis_size),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.x = element_text(size = base_title_size,
    margin = margin(t = 10)),
    legend.position = "none",
    plot.title = element_text(size = base_title_size, face = "bold"),
    plot.title.position = "plot"
  )

p2_fg


# THIRD PANEL: Symbiotic functional guilds per rooting substrate ---------------

## Prepare plotting data 

sym_level <- rev(c(
  "other_mycorrhizal_types",
  "dark_septate_endophytes",
  "other_root_endophytes",
  "arbuscular_mycorrhizal"))

guild_prop <- fungi %>%
  filter(!is.na(BroadSymbiontGroup)) %>%
  mutate(BroadSymbiontGroup = factor(BroadSymbiontGroup, 
    levels = sym_level)) %>%
  distinct(ASV, host_substrate, BroadSymbiontGroup) %>%
  count(host_substrate, BroadSymbiontGroup) %>%
  group_by(host_substrate) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

# Plotting

p_guild <- ggplot(guild_prop,
  aes(x = host_substrate, y = prop, fill = BroadSymbiontGroup)) +
  coord_flip() +
  geom_col(width = 0.7) +
  scale_y_continuous(
    labels = scales::percent_format(),
    expand = expansion(mult = c(0, 0.02))) +
  scale_fill_manual(values = pal_fg, drop = FALSE) +
  labs(
    y = "Proportion of symbiotic ASVs",
    x = NULL,
    fill = NULL) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x = element_text(size = base_axis_size),
    axis.text.y = element_text(size = base_axis_size),
    axis.title.x = element_text(size = base_title_size,
                                margin = margin(t = 10)),
    legend.position = "none",
    plot.title = element_text(size = base_title_size, face = "bold"),
    plot.title.position = "plot")


p_guild


guild_combined <- p1_fg / p2_fg / p_guild  + plot_layout(heights = c(1, 1, 3))

guild_combined

# Save Figure 4-----------------------------------------------------------------

ggsave(
  filename = "../Figures/Figure4.barplots_guild1.png",
  plot = guild_combined,
  width = 13,
  height = 7,
  dpi = 600,
  units = "in",
  bg = "white"
)

## -> This Figure was further edited using power point to make the final Figure 4
