# Project name: root_fungi_epiphytic_terrestrial
# Author: Laura Méndez
# Date: 15.02.26

# Goals: Figure SXX. Non-metric multidimensional scaling (NMDS) ordination 
# of root-associated fungal communities 

# Input: beta_diversity_results.Rdata

# Setup ------------------------------------------------------------------------

rm(list = ls())

library(tidyverse)
library(ggplot2)
library(viridis)
library(ggrepel)

# Load data --------------------------------------------------------------------

load("beta_diversity_results.Rdata")

nmds_df <- beta_diversity_results$nmds_df

# Define color and shape schemes -----------------------------------------------

# Color palette for host families
family_cols <- c(
  "Araceae" = "#00CD6C",
  "Bromeliaceae" = "#009ADE",
  "Orchidaceae" = "#AF58BA",
  "Piperaceae" = "#FFC61E"
)

# Shape palette for habitat (shapes that support fill)
# Triangle (24) for epiphytes, circle (21) for terrestrial
form_shapes <- c("Terrestrial" = 16, "Epiphytic" = 17)

# Plot NMDS colored by family and different shapes for epiphytic vs terrestrial

p1 <- ggplot(nmds_df, aes(NMDS1, NMDS2)) +
  geom_point(
    aes(colour = host_family, shape = host_substrate),
    size = 3,
    alpha = 0.9) +
  stat_ellipse(
    aes(colour = host_family, group = host_family),
    linetype = 2,
    level = 0.95,
    linewidth = 0.7,
    show.legend = FALSE) +
  scale_colour_manual(values = family_cols) +
  scale_shape_manual(values = form_shapes) +
  theme_classic() +
  labs(
    title = "NMDS of fungal communities",
    colour = "Host family",
    shape = "Rooting substrate")

p1

# Save Figure ------------------------------------------------------------------

ggsave(
  filename = "../Figures/FigureSx.NMDS.png",
  plot = p1,
  width = 6,
  height = 4,
  dpi = 600,
  units = "in",
  bg = "white"
)
