# Project name: root_fungi_epiphytic_terrestrial
# Author: Laura Méndez
# Date: 15.03.26

# Goals: Figure SXX.Diversity profile of root-associated fungal α-diversity 
# across Hill orders

# Input: alpha_diversity_results.Rdata

# Setup ------------------------------------------------------------------------

rm(list = ls())

library(tidyverse)
library(ggplot2)
library(purrr)
library(hillR)
library(patchwork)

# Load data --------------------------------------------------------------------

load("analyses3.Rdata")

fungi <- analyses$fungi
fungi$FunctionalGuilds

metadata <- analyses$metadata

asv_tax_fungi <- analyses$asv_tax_fungi

#------------------------------------------------------------
# 1) Stacked bars: Identified vs Unidentified at Phylum
#------------------------------------------------------------

total_asvs <- nrow(asv_tax_fungi)

n_ident <- asv_tax_fungi %>%
  summarise(n = sum(!is.na(Phylum) & Phylum != "")) %>%
  pull(n)

n_unident <- total_asvs - n_ident

df_id <- tibble(
  Status = factor(c("Identified", "Unidentified"),
                  levels = c("Identified", "Unidentified")),
  n = c(n_ident, n_unident)
)

########### For the font size
base_axis_size   <- 15
base_legend_size <- 15
base_title_size  <- 15


p1 <- ggplot(df_id, aes(x = 1, y = n, fill = Status)) +
  geom_col(width = 0.75) +
  coord_flip() +
  scale_fill_manual(values = c("Identified" = "grey30", "Unidentified" = "grey80")) +
  scale_y_continuous(labels = scales::label_comma(),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "Number of ASVs", fill = NULL) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = base_axis_size),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.x = element_text(size = base_title_size),
    
    legend.text = element_text(size = base_legend_size),
    legend.title = element_blank(),
    legend.position = "top",
    
    plot.title = element_text(size = base_title_size, face = "bold"),
    plot.title.position = "plot"
  )


p1

#------------------------------------------------------------
# 2) Stacked bar: top 11 Phylum + Other (12 colors, Set3)
#------------------------------------------------------------

phy <- asv_tax_fungi %>%
  filter(!is.na(Phylum) & Phylum != "") %>%          # only identified
  count(Phylum, name = "n") %>%
  mutate(
    Phylum = if_else(n < 100, "Other", Phylum),
    Phylum = fct_reorder(Phylum, n, .desc = TRUE)
  ) %>%
  group_by(Phylum) %>%
  summarise(n = sum(n), .groups = "drop") %>%
  mutate(Phylum = fct_reorder(Phylum, n, .desc = TRUE))

# make sure palette has enough colors (Set3 max = 12)
n_cat <- nlevels(phy$Phylum)
pal <- RColorBrewer::brewer.pal(8, "Set1")



p2 <- ggplot(phy, aes(x = 1, y = n, fill = Phylum)) +
  geom_col(width = 0.75) +
  coord_flip() +
  scale_fill_manual(values = pal) +
  scale_y_continuous(labels = scales::label_comma(),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "Number of ASVs", fill = NULL) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x = element_text(size = base_axis_size),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.x = element_text(size = base_title_size),
    
    legend.text = element_text(size = base_legend_size),
    legend.title = element_blank(),
    legend.position = "bottom",
    
    plot.title = element_text(size = base_title_size, face = "bold"),
    plot.title.position = "plot"
  )


p2

p1 / p2 + plot_layout(heights = c(1.2, 1))

