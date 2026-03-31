# Project name: root_fungi_epiphytic_terrestrial
# Author: Laura Méndez
# Date: 15.02.26

# Goals: Figure 3. Patterns of similarity and differentiation in 
# root-associated fungal communities.

# Input: analyses.Rdata, alpha_diversity_results.Rdata

# Setup ------------------------------------------------------------------------

rm(list = ls())

library(tidyverse)
library(ggplot2)
library(viridis)
library(ggrepel)

# Load data --------------------------------------------------------------------

load("beta_diversity_results.Rdata")

family_habitat_pairwise <- beta_diversity_results$family_habitat_pairwise
asvs_by_substrate <- beta_diversity_results$asvs_by_substrate
asvs_by_family <- beta_diversity_results$asvs_by_family

# Prepare data for plotting ----------------------------------------------------

# Parse the pair names into two group columns 
family_habitat_pairwise <- family_habitat_pairwise %>%
  mutate(group1 = trimws(sub(" vs .*", "", pairs)),
    group2 = trimws(sub(".* vs ", "", pairs)))

#Define a consistent level order

group_levels <- c(
  "Piperaceae_Terrestrial",
  "Piperaceae_Epiphytic",
  "Orchidaceae_Terrestrial",
  "Orchidaceae_Epiphytic",
  "Bromeliaceae_Terrestrial",
  "Bromeliaceae_Epiphytic",
  "Araceae_Epiphytic",
  "Araceae_Terrestrial"
)

# Keep only the upper triangle  
#    Upper triangle = group1 comes BEFORE group2 in group_levels

min(family_habitat_pairwise$R2)
max(family_habitat_pairwise$R2)

pw_upper <- family_habitat_pairwise %>%
  mutate(g1_idx = match(group1, group_levels),
    g2_idx = match(group2, group_levels)) %>%
  mutate(grp1_final = ifelse(g1_idx < g2_idx, group1, group2),
    grp2_final = ifelse(g1_idx < g2_idx, group2, group1)) %>%
  select(group1 = grp1_final, group2 = grp2_final, R2, p.value, p.adj) %>%
  distinct()

# Apply factor levels 

pw_upper <- pw_upper %>%
  mutate(group1 = factor(group1, levels = group_levels),
    group2 = factor(group2, levels = group_levels))

# Plot heatmap of pairwise.adonis results --------------------------------------

ggplot(pw_upper, aes(x = group1, y = group2, fill = R2)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = ifelse(p.adj >= 0.05, "ns", "")),
    size = 8,
    color = "black") +
  scale_fill_viridis(option = "viridis",
    direction = -1,
    name = expression(R^2),
    limits = c(0.05, 0.15),  
    breaks = seq(0.05, 0.15, by = 0.03)) +
  coord_equal() +
  scale_x_discrete(position = "top", limits = group_levels) +
  scale_y_discrete(limits = group_levels) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 14) +
  theme(panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 0))


# Save part of Figure 3 --------------------------------------------------------

ggsave(
  filename = "../Figures/Fig_beta_diversity_heatmap2.png",
  plot = last_plot(),
  width = 13,
  height = 13,
  units = "in",
  dpi = 600,
  bg = "white"
)

# Plot Venn diagrams -----------------------------------------------------------

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

## -> This Figure was further edited using power point to make the final Figure 3 
