# Project name: root_fungi_epiphytic_terrestrial
# Author: Laura Méndez
# Date: [Date]

# Goals:
# Create publication-quality figures for beta diversity analyses

# Figures created:
# - Fig 1: NMDS ordination colored by family
# - Fig 2: NMDS ordination colored by growth form
# - Fig 3: db-RDA ordination with biplot arrows
# - Fig 4: Betadisper - distance to centroid by species

# Input: analyses.Rdata, beta_diversity_results.Rdata
# Dependencies: ggplot2, tidyverse, ggrepel

################################################################################
# Setup
################################################################################

rm(list = ls())
library(tidyverse)
library(ggplot2)
library(viridis)
library(ggrepel)

load("analyses.Rdata")
load("beta_diversity_results2.Rdata")

# Colors
family_cols <- c(
  "Araceae" = "#00CD6C",
  "Bromeliaceae" = "#009ADE",
  "Orchidaceae" = "#AF58BA",
  "Piperaceae" = "#FFC61E"
)

#form_cols <- c("E" = "#40B0A6", "T" = "#E1BE6A")

################################################################################
# Figure 1: NMDS colored by family
################################################################################


p1 <- ggplot(nmds_df,
             aes(NMDS1, NMDS2, color = host_family, shape = host_habitat)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = family_cols) +
  stat_ellipse(aes(group = host_family), linetype = 2, level = 0.95) +
  theme_classic() +
  labs(
    title = "NMDS of fungal communities",
    color = "Host family",
    shape = "Growth form"
  )

p1

ggsave(
  filename = "../Figures/FigureSx.NMDS.png",
  plot = p1,
  width = 6,
  height = 4,
  dpi = 600,
  units = "in",
  bg = "white"
)

#ggsave("Fig_NMDS_by_family.pdf", p1, width = 8, height = 6)

################################################################################
# Figure 2: heatmap of pairwise.adonis results
################################################################################

# ── 1. Parse the pair names into two group columns ──────────────────────────
family_habitat_pairwise <- family_habitat_pairwise %>%
  mutate(
    group1 = trimws(sub(" vs .*", "", pairs)),
    group2 = trimws(sub(".* vs ", "", pairs))
  )

# ── 2. Define a consistent level order ──────────────────────────────────────
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

# ── 3. Keep only the upper triangle  ────────────────────────────────────────
#    Upper triangle = group1 comes BEFORE group2 in group_levels

min(family_habitat_pairwise$R2)
max(family_habitat_pairwise$R2)

pw_upper <- family_habitat_pairwise %>%
  mutate(
    g1_idx = match(group1, group_levels),
    g2_idx = match(group2, group_levels)
  ) %>%
  mutate(
    grp1_final = ifelse(g1_idx < g2_idx, group1, group2),
    grp2_final = ifelse(g1_idx < g2_idx, group2, group1)
  ) %>%
  select(group1 = grp1_final, group2 = grp2_final, R2, p.value, p.adj) %>%
  distinct()
# ── 4. Apply factor levels ───────────────────────────────────────────────────
pw_upper <- pw_upper %>%
  mutate(
    group1 = factor(group1, levels = group_levels),
    group2 = factor(group2, levels = group_levels)
  )

# ── 5. Plot ──────────────────────────────────────────────────────────────────
ggplot(pw_upper, aes(x = group1, y = group2, fill = R2)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(
    aes(label = ifelse(p.adj >= 0.05, "ns", "")),
    size = 5,
    color = "black"
  ) +
  scale_fill_viridis(
    option = "viridis",
    direction = -1,
    name = expression(R^2),
    limits = c(0.05, 0.15),  # adjust to your data range
    breaks = seq(0.05, 0.15, by = 0.03)  # or list them manually
  ) +
  coord_equal() +
  scale_x_discrete(position = "top", limits = group_levels) +
  scale_y_discrete(limits = group_levels) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 0)
  )




ggsave(
  filename = "../Figures/Fig_beta_diversity_heatmap1.png",
  plot = last_plot(),
  width = 10,
  height = 11,
  units = "in",
  dpi = 600,
  bg = "white"
)



################################################################################
# Figure 3: db-RDA ordination
################################################################################

# Extract components for biplot
dbrda_rmod <- beta_diversity_results$dbrda_model
cap_variance <- beta_diversity_results$cap_variance

# Site scores
sites <- beta_diversity_results$dbrda_sites

# Biplot arrows (constraints)
bp <- as.data.frame(scores(dbrda_rmod, display = "bp", choices = 1:2))
bp$term <- rownames(bp)

# Arrow scaling
arrow_mult <- 1.2

p3 <- ggplot(sites, aes(CAP1, CAP2)) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_vline(xintercept = 0, linewidth = 0.3) +
  stat_ellipse(
    aes(color = host_family, group = host_family),
    type = "t",
    level = 0.95,
    linewidth = 0.8,
    linetype = 2
  ) +
  geom_point(
    aes(color = host_family, shape = host_habitat),
    size = 3, alpha = 0.5
  ) +
  geom_segment(
    data = bp,
    aes(x = 0, y = 0, xend = CAP1 * arrow_mult, yend = CAP2 * arrow_mult),
    linewidth = 0.5,
    arrow = arrow(length = unit(0.22, "cm"))
  ) +
  geom_text_repel(
    data = bp,
    aes(x = CAP1 * arrow_mult, y = CAP2 * arrow_mult, label = term),
    size = 3,
    fontface = "bold",
    segment.color = NA
  ) +
  coord_equal() +
  scale_color_manual(values = family_cols) +
  theme_classic(base_size = 11) +
  labs(
    x = paste0("CAP1 (", cap_variance["CAP1"], "%)"),
    y = paste0("CAP2 (", cap_variance["CAP2"], "%)"),
    color = "Host family",
    shape = "Growth form"
  )

ggsave("Fig_dbRDA.pdf", p3, width = 9, height = 7)

################################################################################
# Figure 4: Betadisper - distance to centroid by species
################################################################################

disp2 <- beta_diversity_results$betadisper_species

disp2_df <- data.frame(
  distance = disp2$distances,
  host_wcvp_name = meta_roots$host_wcvp_name,
  host_family = meta_roots$host_family
) %>%
  mutate(
    host_family = factor(
      host_family,
      levels = c("Araceae", "Bromeliaceae", "Orchidaceae", "Piperaceae")
    ),
    host_wcvp_name = fct_reorder(host_wcvp_name, distance, .fun = median)
  )

p4 <- ggplot(disp2_df, aes(x = host_wcvp_name, y = distance, fill = host_family)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.1, size = 1.8, alpha = 0.7) +
  coord_flip() +
  scale_fill_manual(values = family_cols) +
  theme_classic() +
  labs(
    x = "Host species",
    y = "Distance to species centroid (Jaccard)",
    fill = "Host family"
  ) +
  facet_grid(
    host_family ~ .,
    scales = "free_y",
    space = "free_y",
    switch = "y"
  ) +
  theme(
    strip.background = element_blank(),
    strip.text.y = element_text(face = "bold", angle = 0)
  )

ggsave("Fig_betadisper_by_species.pdf", p4, width = 8, height = 10)

cat("\n✓ Beta diversity figures created!\n")
cat("Figures saved:\n")
cat("  - Fig_NMDS_by_family.pdf\n")
cat("  - Fig_NMDS_by_form.pdf\n")
cat("  - Fig_dbRDA.pdf\n")
cat("  - Fig_betadisper_by_species.pdf\n")
