# Project name: root_fungi_epiphytic_terrestrial
# Author: Laura Méndez
# Date: [Date]

# Goals:
# Create publication-quality figures for alpha diversity analyses

# Figures created:
# - Fig 1: ASV richness by growth form (all fungi)
# - Fig 2: ASV richness by host family (all fungi)
# - Fig 3: ASV richness by growth form × family (all fungi)
# - Fig 4: ASV richness by growth form (root-associated)
# - Fig 5: ASV richness by host family (root-associated)
# - Fig 6: ASV richness by species (root-associated, faceted by family)

# Input: analyses.Rdata, alpha_diversity_results.Rdata
# Dependencies: tidyverse, ggplot2, ggforce, ggdist, forcats

################################################################################
# Setup
################################################################################
#setwd("C:/Users/ploti/Dropbox/Work/03_Marburg/02_epiphites_terrestrial_root_funghi/R/scripts")

rm(list = ls())

library(tidyverse)
library(ggplot2)
library(ggdist)
library(forcats)
library(patchwork)

load("alpha_diversity_results2.Rdata")

################################################################################
# Define color and shape schemes -----------------------------------------------
################################################################################

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

################################################################################
# Prepare data with consistent habitat labels -----------------------------
################################################################################

# Recode habitat labels for clarity
# E → Epiphytes, T → Terrestrial
alpha_all_clean <- alpha_diversity_results$alpha_all %>%
  mutate(host_habitat = recode(host_habitat,
  "E" = "Epiphytic", "T" = "Terrestrial")) 

################################################################################
# TOP PANEL: Overall habitat comparison ------------------------------------
################################################################################

p_top <- ggplot(
  alpha_all_clean,
  aes(x = ASV_richness, y = host_habitat)
) +
  # Individual points with jitter to avoid overplotting
  geom_jitter(
    aes(shape = host_habitat),
    height = 0.08,
    size = 1.8,
    alpha = 0.5,
    fill = "black",
    color = "black"
  ) +
  # Boxplot showing median and quartiles
  geom_boxplot(
    width = 0.3,
    outlier.shape = NA,  # Don't show outliers (already shown as points)
    alpha = 0
  ) +
  # Position x-axis on top for visual alignment with bottom panel
  scale_x_continuous(position = "top") +
  scale_shape_manual(values = form_shapes) +
  theme_classic() +
  labs(x = NULL, y = NULL) +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 13),
    axis.title = element_text(size = 15),
    plot.margin = margin(5.5, 5.5, 0, 5.5)  # Remove bottom margin
  )

p_top

################################################################################
# BOTTOM PANEL: Detailed species-level view ------------------------------------
################################################################################

# Step 1: Order species within families by mean richness
alpha_roots_ordered <- alpha_all_clean %>%
  mutate(
    host_family = factor(
      host_family,
      levels = c("Araceae", "Bromeliaceae", "Orchidaceae", "Piperaceae")
    )
  ) %>%
  group_by(host_family, name_analyses) %>%
  mutate(mean_richness = mean(ASV_richness)) %>%
  ungroup() %>%
  mutate(
    # Order species by family and mean richness
    host_species_ord = fct_reorder2(name_analyses, host_family, mean_richness),
    # Add dummy level for summary statistics (will appear at bottom after coord_flip)
    host_species_ord = fct_expand(host_species_ord, "Overall mean per habitat"),
    host_species_ord = fct_relevel(host_species_ord, "Overall mean per habitat", after = 0)
  )

# Step 2: Calculate growth form means ± SD per family (for summary row)
gf_summary <- alpha_roots_ordered %>%
  group_by(host_family, host_habitat) %>%
  summarise(
    mean_gf = mean(ASV_richness, na.rm = TRUE),
    sd_gf = sd(ASV_richness, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    host_species_ord = factor(
      "Overall mean per habitat",
      levels = levels(alpha_roots_ordered$host_species_ord)
    )
  )

# Step 3: Calculate family mean across both growth forms (for reference line)
family_mean <- alpha_roots_ordered %>%
  group_by(host_family) %>%
  summarise(mean_family = mean(ASV_richness, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    host_species_ord = factor(
      "Overall mean per habitat",
      levels = levels(alpha_roots_ordered$host_species_ord)
    )
  )

# Step 4: Create bottom panel
p_bottom <- ggplot(
  alpha_roots_ordered,
  aes(x = host_species_ord, y = ASV_richness, fill = host_family)
) +
  # Family mean as dashed horizontal reference line
  geom_hline(
    data = family_mean,
    aes(yintercept = mean_family),
    inherit.aes = FALSE,
    linetype = "dashed",
    linewidth = 0.9,
    color = "black"
  ) +
  # Connect samples within species to show variation
  geom_line(
    aes(group = name_analyses, color = host_family),
    linewidth = 0.9
  ) +
  # Individual sample points (shape indicates growth form)
  geom_point(
    aes(shape = host_habitat),
    size = 3,
    alpha = 0.7,
    fill = "black",
    color = "black",
    stroke = 0.3
  ) +
  # Growth form summary: error bars (mean ± SD)
  geom_errorbar(
    data = gf_summary,
    aes(
      x = host_species_ord,
      ymin = mean_gf - sd_gf,
      ymax = mean_gf + sd_gf,
      group = host_habitat
    ),
    inherit.aes = FALSE,
    width = 0.12,
    linewidth = 0.6,
    color = "black",
    position = pd
  ) +
  # Growth form summary: mean points
  geom_point(
    data = gf_summary,
    aes(
      x = host_species_ord,
      y = mean_gf,
      shape = host_habitat,
      fill = host_family,
      group = host_habitat
    ),
    inherit.aes = FALSE,
    size = 4,
    color = "black",
    stroke = 0.8,
    position = pd
  ) +
  # Aesthetics
  scale_shape_manual(name = "Host habitat", values = form_shapes) +
  scale_fill_manual(name = "Host family", values = family_cols) +
  scale_color_manual(name = "Host family", values = family_cols) +
  coord_flip() +
  theme_classic() +
  labs(y = "ASV richness (α-diversity)", x = NULL) +
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
    legend.position = "none"
  )


p_bottom

################################################################################
# Combine panels ---------------------------------------------------------------
################################################################################

# Stack top and bottom panels with 80:300 height ratio
# Top panel shows overall pattern, bottom shows detailed species view
p_combined <- p_top / p_bottom +
  plot_layout(heights = c(0.8, 3))

# Display combined figure
p_combined

getwd()
ggsave(
  filename = "../Figures/Figure2.alpha_diversity_combined.png",
  plot = p_combined,
  width = 10,
  height = 15,
  dpi = 600,
  units = "in",
  bg = "white"
)





### Forest plots for hill numbers models



library(dplyr)
library(stringr)
library(ggplot2)
library(broom.mixed)

get_forest_df <- function(model) {
  broom.mixed::tidy(model, effects = "fixed") %>%
    dplyr::filter(term != "(Intercept)") %>%
    dplyr::mutate(
      IRR = exp(estimate),
      CI_low = exp(estimate - 1.96 * std.error),
      CI_high = exp(estimate + 1.96 * std.error),
      sig = p.value < 0.05,
      alpha_sig = ifelse(sig, 1, 0.35),
      term_clean = case_when(
        term == "host_substrateT" ~ "Substrate: Terrestrial",
        term == "host_familyBromeliaceae" ~ "Family: Bromeliaceae",
        term == "host_familyOrchidaceae" ~ "Family: Orchidaceae",
        term == "host_familyPiperaceae" ~ "Family: Piperaceae",
        TRUE ~ term
      ))
}

  
df_q0 <- get_forest_df(m_q0) %>% mutate(model = "q = 0 (ASV richness)")
df_q1 <- get_forest_df(m_q1) %>% mutate(model = "q = 1 (common ASVs)")
df_q2 <- get_forest_df(m_q2) %>% mutate(model = "q = 2 (dominant ASVs)")

lvl <- c("q = 0 (ASV richness)", "q = 1 (common ASVs)", "q = 2 (dominant ASVs)")

df_all <- bind_rows(df_q0, df_q1, df_q2) %>%
  mutate(model = factor(model, levels = rev(lvl))) 

pd <- position_dodge(width = 0.6)
cols <- c(
  "q = 0 (ASV richness)" = "#C40F5B",
  "q = 1 (common ASVs)" = "#FD8D3C",
  "q = 2 (dominant ASVs)" = "#089099"
)

p_all <- ggplot(df_all, aes(x = IRR, y = term_clean, color = model)) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.7) +
  geom_errorbarh(
    aes(xmin = CI_low, xmax = CI_high, alpha = alpha_sig),
    position = pd, height = 0.15, linewidth = 0.9
  ) +
  geom_point(aes(alpha = alpha_sig), position = pd, size = 2.6) +
  scale_alpha_identity() +
  scale_x_log10() +
  scale_color_manual(values = cols,breaks = lvl,   # legend order
                     limits = lvl    # ensures consistent ordering
  ) +
  guides(color = guide_legend(reverse = FALSE)) +
  labs(
    x = "Incidence rate ratio",
    y = NULL,
    color = "Hill number"
  ) +
  theme_classic(base_size = 15) +
  theme(legend.position = "top")

p_all

ggsave(
  filename = "../Figures/Figure2.forest_plot_ALL_hill.png",
  plot = p_all,
  width = 6,
  height = 7,
  dpi = 600,
  units = "in",
  bg = "white"
)

### One plot per Hill number

fp_q0 <- ggplot(df_q0, aes(x = IRR, y = term_clean)) +
    geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.7) +
    geom_errorbarh(aes(xmin = CI_low, xmax = CI_high, alpha = alpha_sig),
                   height = 0.1, linewidth = 0.9) +
    geom_pointrange(aes(xmin = CI_low, xmax = CI_high, alpha = alpha_sig),
                    linewidth = 0.7) +
    scale_alpha_identity() +
    scale_x_log10() +
    coord_cartesian(clip = "off") +
    labs(
      x = "Incidence rate ratio (IRR)",
      y = NULL,
      title = NULL
    ) +
    theme_classic(base_size = 14) +
    theme(axis.title.x = element_text(margin = margin(t = 8)))




fp_q0

fp_q1 <- ggplot(df_q1, aes(x = IRR, y = term_clean)) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.7) +
  geom_errorbarh(aes(xmin = CI_low, xmax = CI_high, alpha = alpha_sig),
                 height = 0.1, linewidth = 0.9) +
  geom_pointrange(aes(xmin = CI_low, xmax = CI_high, alpha = alpha_sig),
                  linewidth = 0.7) +
  scale_alpha_identity() +
  scale_x_log10() +
  coord_cartesian(clip = "off") +
  labs(
    x = "Incidence rate ratio (IRR)",
    y = NULL,
    title = NULL
  ) +
  theme_classic(base_size = 14) +
  theme(axis.title.x = element_text(margin = margin(t = 8)))




fp_q1

fp_q2 <- ggplot(df_q2, aes(x = IRR, y = term_clean)) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.7) +
  geom_errorbarh(aes(xmin = CI_low, xmax = CI_high, alpha = alpha_sig),
                 height = 0.1, linewidth = 0.9) +
  geom_pointrange(aes(xmin = CI_low, xmax = CI_high, alpha = alpha_sig),
                  linewidth = 0.7) +
  scale_alpha_identity() +
  scale_x_log10() +
  coord_cartesian(clip = "off") +
  labs(
    x = "Incidence rate ratio (IRR)",
    y = NULL,
    title = NULL
  ) +
  theme_classic(base_size = 14) +
  theme(axis.title.x = element_text(margin = margin(t = 8)))




fp_q2



library(patchwork)

plot_forest_hill <- fp_q0 | fp_q1 | fp_q2

plot_forest_hill


ggsave(
  filename = "../Figures/Figure2.forest_plot_hill.png",
  plot = plot_forest_hill,
  width = 15,
  height = 6,
  dpi = 600,
  units = "in",
  bg = "white"
)
############################## diversity profile

library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(hillR)

# q range
qs <- seq(0, 3, by = 0.1)

prof_all <- map_dfr(qs, function(qi) {
  qD <- hill_taxa(comm, q = qi) # vector over samples
  tibble(
    q = qi,
    mean_qD = mean(qD, na.rm = TRUE),
    se_qD   = sd(qD, na.rm = TRUE) / sqrt(sum(!is.na(qD)))
  )
})

pts <- prof_all %>% filter(q %in% c(0, 1, 2))

diversity_profile <- ggplot(prof_all, aes(q, mean_qD)) +
  geom_line(linewidth = 1) +
  geom_point(data = pts, size = 2.5) +
  geom_ribbon(aes(ymin = mean_qD - 1.96*se_qD, ymax = mean_qD + 1.96*se_qD),
              alpha = 0.15, colour = NA) +
  geom_text(
    data = pts,
    aes(label = paste0("q", q, " = ", round(mean_qD, 1))),
    nudge_y = max(prof_all$mean_qD) * 0.04,
    hjust = 0
  ) +
  labs(x = "Order q", y = "Hill number") +
  theme_classic()

diversity_profile

ggsave(
  filename = "../Figures/FigureSx.diversity_profile.png",
  plot = diversity_profile,
  width = 6,
  height = 4,
  dpi = 600,
  units = "in",
  bg = "white"
)


################################################################################
# Figure 1: Alpha diversity by growth form (all fungi)
################################################################################

p1 <- ggplot(
  alpha_diversity_results$alpha_all,
  aes(`epiphitic/terrestrial`, ASV_richness, fill = `epiphitic/terrestrial`)
) +
  stat_halfeye(adjust = 0.5, width = 0.6, .width = 0) +
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.6) +
  geom_jitter(width = 0.08, size = 2, alpha = 0.6) +
  scale_fill_brewer(palette = "Dark2") +
  theme_classic() +
  labs(
    y = "ASV richness",
    x = "Growth form",
    title = "Alpha diversity of root-associated fungi",
    fill = "Growth form"
  )

p1

#ggsave("Fig_alpha_by_form_all.pdf", p1, width = 6, height = 5)

################################################################################
# Figure 2: Alpha diversity by family (all fungi)
################################################################################

p2 <- ggplot(
  alpha_diversity_results$alpha_all,
  aes(host_family, ASV_richness, fill = host_family)
) +
  stat_halfeye(adjust = 0.5, width = 0.6, .width = 0) +
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.6) +
  geom_jitter(width = 0.08, size = 2, alpha = 0.6) +
  theme_classic() +
  labs(
    y = "ASV richness",
    x = "Host Family",
    title = "Alpha diversity of root-associated fungi",
    fill = "Family"
  )

p2

#ggsave("Fig_alpha_by_family_all.pdf", p2, width = 7, height = 5)

################################################################################
# Figure 3: Alpha diversity by growth form within families (all fungi)
################################################################################

p3 <- ggplot(
  alpha_diversity_results$alpha_all,
  aes(`epiphitic/terrestrial`, ASV_richness, fill = `epiphitic/terrestrial`)
) +
  stat_halfeye(adjust = 0.5, width = 0.6, .width = 0) +
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.6) +
  geom_jitter(width = 0.08, size = 2, alpha = 0.6) +
  facet_wrap(~host_family, scales = "free_y") +
  scale_fill_brewer(palette = "Dark2") +
  theme_classic() +
  labs(
    x = "Growth form",
    y = "ASV richness",
    title = "Alpha diversity across growth forms within families",
    fill = "Growth form"
  )


p3


#ggsave("Fig_alpha_by_form_family_all.pdf", p3, width = 10, height = 6)

str(alpha_df)

# Prepare data with ordered species
alpha_roots_plot <- alpha_diversity_results$alpha_all %>%
  mutate(
    host_family = factor(
      host_family,
      levels = c("Araceae", "Bromeliaceae", "Orchidaceae", "Piperaceae")
    )
  ) %>%
  group_by(host_family, name_analyses) %>%
  mutate(mean_richness = mean(ASV_richness)) %>%
  ungroup() %>%
  mutate(
    host_species_ord = fct_reorder2(
      name_analyses,
      host_family,
      mean_richness
    )
  )

p6 <- ggplot(
  alpha_roots_plot,
  aes(host_species_ord, ASV_richness, fill = host_family)
) +
  stat_halfeye(adjust = 0.5, width = 0.6, .width = 0) +
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.6) +
  geom_jitter(
    aes(shape = `epiphitic/terrestrial`),
    width = 0.08, size = 2, alpha = 0.7
  ) +
  scale_fill_manual(values = family_cols) +
  scale_shape_manual(values = c("E" = 17, "T" = 16)) +
  coord_flip() +
  theme_classic() +
  labs(
    y = "ASV richness",
    x = "Host plant species/genus",
    title = "Alpha diversity by species",
    fill = "Family",
    shape = "Growth form"
  ) +
  facet_grid(
    host_family ~ .,
    scales = "free_y",
    space = "free_y",
    switch = "y"
  ) +
  theme(
    strip.background = element_blank(),
    strip.text.y = element_text(angle = 0, face = "bold")
  )


p6




################################################################################
# Figure 4: Alpha diversity by growth form (root-associated, with colors)
################################################################################

p4 <- ggplot(
  alpha_diversity_results$alpha_roots,
  aes(`epiphitic/terrestrial`, ASV_richness, fill = `epiphitic/terrestrial`)
) +
  stat_halfeye(adjust = 0.5, width = 0.6, .width = 0) +
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.6) +
  geom_jitter(width = 0.08, size = 2, alpha = 0.6) +
  scale_fill_manual(values = form_cols) +
  theme_classic() +
  labs(
    y = "ASV richness",
    x = "Growth form",
    title = "Alpha diversity of root-associated fungi",
    fill = "Growth form"
  )

p4


ggsave("Fig_alpha_by_form_roots.pdf", p4, width = 6, height = 5)

################################################################################
# Figure 5: Alpha diversity by family (root-associated, with colors)
################################################################################

p5 <- ggplot(
  alpha_diversity_results$alpha_roots,
  aes(host_family, ASV_richness, fill = host_family)
) +
  stat_halfeye(adjust = 0.5, width = 0.6, .width = 0) +
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.6) +
  geom_jitter(width = 0.08, size = 2, alpha = 0.6) +
  scale_fill_manual(values = family_cols) +
  theme_classic() +
  labs(
    y = "ASV richness",
    x = "Host Family",
    title = "Alpha diversity of root-associated fungi",
    fill = "Family"
  )

p5

ggsave("Fig_alpha_by_family_roots.pdf", p5, width = 7, height = 5)

################################################################################
# Figure 6: Alpha diversity by species (root-associated, faceted by family)
################################################################################

# Prepare data with ordered species
alpha_roots_plot <- alpha_diversity_results$alpha_roots %>%
  mutate(
    host_family = factor(
      host_family,
      levels = c("Araceae", "Bromeliaceae", "Orchidaceae", "Piperaceae")
    )
  ) %>%
  group_by(host_family, name_analyses) %>%
  mutate(mean_richness = mean(ASV_richness)) %>%
  ungroup() %>%
  mutate(
    host_species_ord = fct_reorder2(
      name_analyses,
      host_family,
      mean_richness
    )
  )

p6 <- ggplot(
  alpha_roots_plot,
  aes(host_species_ord, ASV_richness, fill = host_family)
) +
  stat_halfeye(adjust = 0.5, width = 0.6, .width = 0) +
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.6) +
  geom_jitter(
    aes(shape = `epiphitic/terrestrial`),
    width = 0.08, size = 2, alpha = 0.7
  ) +
  scale_fill_manual(values = family_cols) +
  scale_shape_manual(values = c("E" = 17, "T" = 16)) +
  coord_flip() +
  theme_classic() +
  labs(
    y = "ASV richness",
    x = "Host plant species/genus",
    title = "Alpha diversity by species",
    fill = "Family",
    shape = "Growth form"
  ) +
  facet_grid(
    host_family ~ .,
    scales = "free_y",
    space = "free_y",
    switch = "y"
  ) +
  theme(
    strip.background = element_blank(),
    strip.text.y = element_text(angle = 0, face = "bold")
  )


p6

ggsave("Fig_alpha_by_species.pdf", p6, width = 8, height = 12)

cat("\n✓ Alpha diversity figures created!\n")
cat("Figures saved:\n")
cat("  - Fig_alpha_by_form_all.pdf\n")
cat("  - Fig_alpha_by_family_all.pdf\n")
cat("  - Fig_alpha_by_form_family_all.pdf\n")
cat("  - Fig_alpha_by_form_roots.pdf\n")
cat("  - Fig_alpha_by_family_roots.pdf\n")
cat("  - Fig_alpha_by_species.pdf\n")
