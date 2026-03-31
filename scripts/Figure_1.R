# Project name: root_fungi_epiphytic_terrestrial
# Author: Laura Méndez
# Date: 15.02.26

# Goals: Figure 1

# Input: analyses.Rdata, alpha_diversity_results.Rdata

# Setup ------------------------------------------------------------------------

rm(list = ls())

library(tidyverse)
library(ggplot2)
library(ggdist)
library(forcats)
library(patchwork)

# Load data --------------------------------------------------------------------

load("alpha_diversity_results.Rdata")

alpha_df <- alpha_diversity_results$alpha_df

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
form_shapes <- c("Epiphytic" = 24, "Terrestrial" = 21)

# Position dodge for non-overlapping error bars
pd <- position_dodge(width = 0.35)

# Recode habitat labels for clarity
# E → Epiphytes, T → Terrestrial
alpha_df <- alpha_df %>%
  mutate(host_substrate = dplyr::recode(host_substrate,
                               "E" = "Epiphytic", "T" = "Terrestrial")) 

# TOP PANEL: Overall habitat comparison ------------------------------------

p_top <- ggplot(
  alpha_df,
  aes(x = ASV_richness, y = host_substrate)) +
  # Individual points with jitter to avoid overplotting
  geom_jitter(aes(shape = host_substrate),
    height = 0.08,
    size = 1.8,
    alpha = 0.5,
    fill = "black",
    color = "black") +
  # Boxplot showing median and quartiles
  geom_boxplot(width = 0.3,
    outlier.shape = NA,  # Don't show outliers (already shown as points)
    alpha = 0) +
  # Position x-axis on top for visual alignment with bottom panel
  scale_x_continuous(position = "top") +
  scale_shape_manual(values = form_shapes) +
  theme_classic() +
  labs(x = NULL, y = NULL) +
  theme(legend.position = "none",
    axis.text = element_text(size = 13),
    axis.title = element_text(size = 15),
    plot.margin = margin(5.5, 5.5, 0, 5.5)  # Remove bottom margin
  )

p_top

# BOTTOM PANEL: Detailed species-level view ------------------------------------

# Order species within families by mean richness
alpha_roots_ordered <- alpha_df %>%
  mutate(host_family = factor(host_family,
      levels = c("Araceae", "Bromeliaceae", "Orchidaceae", "Piperaceae"))) %>%
  group_by(host_family, host_species_identity) %>%
  mutate(mean_richness = mean(ASV_richness)) %>%
  ungroup() %>%
  mutate(
    # Order species by family and mean richness
    host_species_ord = fct_reorder2(host_species_identity, host_family, mean_richness),
    # Add dummy level for summary statistics (will appear at bottom after coord_flip)
    host_species_ord = fct_expand(host_species_ord, "Overall mean per habitat"),
    host_species_ord = fct_relevel(host_species_ord, "Overall mean per habitat", after = 0))

# Calculate growth form means ± SD per family (for summary row)
gf_summary <- alpha_roots_ordered %>%
  group_by(host_family, host_substrate) %>%
  summarise(
    mean_gf = mean(ASV_richness, na.rm = TRUE),
    sd_gf = sd(ASV_richness, na.rm = TRUE),
    .groups = "drop") %>%
  mutate(host_species_ord = factor("Overall mean per habitat",
      levels = levels(alpha_roots_ordered$host_species_ord)))

# Calculate family mean across both growth forms (for reference line)
family_mean <- alpha_roots_ordered %>%
  group_by(host_family) %>%
  summarise(mean_family = mean(ASV_richness, na.rm = TRUE), .groups = "drop") %>%
  mutate(host_species_ord = factor("Overall mean per habitat",
      levels = levels(alpha_roots_ordered$host_species_ord)))

# Create bottom panel
p_bottom <- ggplot(
  alpha_roots_ordered,
  aes(x = host_species_ord, y = ASV_richness, fill = host_family)) +
  # Family mean as dashed horizontal reference line
  geom_hline(
    data = family_mean,
    aes(yintercept = mean_family),
    inherit.aes = FALSE,
    linetype = "dashed",
    linewidth = 0.9,
    color = "black") +
  # Connect samples within species to show variation
  geom_line(aes(group = host_species_identity, color = host_family),
    linewidth = 0.9) +
  # Individual sample points (shape indicates growth form)
  geom_point(aes(shape = host_substrate),
    size = 3,
    alpha = 0.7,
    fill = "black",
    color = "black",
    stroke = 0.3) +
  # Growth form summary: error bars (mean ± SD)
  geom_errorbar(data = gf_summary,
    aes(x = host_species_ord,
      ymin = mean_gf - sd_gf,
      ymax = mean_gf + sd_gf,
      group = host_substrate),
    inherit.aes = FALSE,
    width = 0.12,
    linewidth = 0.6,
    color = "black",
    position = pd) +
  # Growth form summary: mean points
  geom_point(data = gf_summary,
    aes(x = host_species_ord,
      y = mean_gf,
      shape = host_substrate,
      fill = host_family,
      group = host_substrate),
    inherit.aes = FALSE,
    size = 4,
    color = "black",
    stroke = 0.8,
    position = pd) +
  # Aesthetics
  scale_shape_manual(name = "Host rooting substrate", values = form_shapes) +
  scale_fill_manual(name = "Host family", values = family_cols) +
  scale_color_manual(name = "Host family", values = family_cols) +
  coord_flip() +
  theme_classic() +
  labs(y = "ASV richness (α-diversity) Hill q = 0", x = NULL) +
  facet_wrap(~host_family, ncol = 1, scales = "free_y") +
  theme(
    # Facet styling
    strip.background = element_blank(),
    strip.text = element_blank(),  # This removes the family name labels
    # Text sizes
    text = element_text(size = 14),
    axis.title.x = element_text(size = 15),
    axis.text = element_text(size = 13),
    plot.title = element_text(size = 14, face = "bold"),
    # Hide legend (visual elements are self-explanatory with facets)
    legend.position = "none")


p_bottom

# Combine panels ---------------------------------------------------------------

# Stack top and bottom panels with 80:300 height ratio
# Top panel shows overall pattern, bottom shows detailed species view
p_combined <- p_top / p_bottom +
  plot_layout(heights = c(0.8, 3))

# Display combined figure
p_combined

# Save Figure 1 ----------------------------------------------------------------

ggsave(
  filename = "../Figures/Figure2.alpha_diversity_combined.png",
  plot = p_combined,
  width = 10,
  height = 15,
  dpi = 600,
  units = "in",
  bg = "white"
)

