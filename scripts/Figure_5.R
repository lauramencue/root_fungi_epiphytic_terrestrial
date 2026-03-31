# Project name: root_fungi_epiphytic_terrestrial
# Author: Laura Méndez
# Date: 15.02.26

# Goals: Figure 5. Rooting substrate and host-family effects on the 
# proportion of symbiotic root-associated fungi and their relative abundance

# Input: guild_composition_results.Rdata

# Setup ------------------------------------------------------------------------

rm(list = ls())

library(tidyverse)
library(ggplot2)
library(broom.mixed)
library(forcats)
library(patchwork)

# Load data --------------------------------------------------------------------

load("guild_composition_results.Rdata")

m_symb <- guild_composition_results$m_symb
fits_nb <- guild_composition_results$fits_nb

# FIRST PANEL: Forest plot from GLMM model output proportion of symbionts ------

# tidy fixed effects from the NO-interaction model
coefs <- broom.mixed::tidy(m_symb, effects = "fixed") %>%
  filter(term != "(Intercept)") %>%
  mutate(
    OR = exp(estimate),
    CI_low  = exp(estimate - 1.96 * std.error),
    CI_high = exp(estimate + 1.96 * std.error),
    significance = ifelse(p.value < 0.05, "Significant", "Not significant"),
    alpha_val = ifelse(significance == "Significant", 1, 0.35),
    term_clean = case_when(
      term == "host_substrateTerrestrial" ~ "Substrate: Terrestrial",
      term == "host_familyBromeliaceae" ~ "Family: Bromeliaceae",
      term == "host_familyOrchidaceae" ~ "Family: Orchidaceae",
      term == "host_familyPiperaceae" ~ "Family: Piperaceae",
      TRUE ~ term),
    term_clean = factor(term_clean, levels = (unique(term_clean))))

forest_plot <- ggplot(coefs, aes(x = OR, y = term_clean)) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.6) +
  geom_errorbarh(aes(xmin = CI_low, xmax = CI_high, alpha = alpha_val),
                 height = 0.05, linewidth = 0.9) +
  geom_point(aes(alpha = alpha_val), size = 3) +
  scale_x_log10(
    breaks = c(0.5, 0.75, 1, 1.4, 2, 3),
    labels = c("0.5","0.75","1","1.4","2","3")) +
  scale_alpha_identity() +
  guides(alpha = "none") +
  coord_cartesian(clip = "off") +
  labs(x = "Odds ratio",
       y = NULL) +
  theme_classic(base_size = 15) +
  theme(axis.title.x = element_text(margin = margin(t = 8)))

forest_plot

# SECOND PANEL: Only rooting substrate results ---------------------------------
# from individual GLMMs of guild proportions -----------------------------------

# Extract model outputs

forest_dat <- imap_dfr(fits_nb, function(model, guild) {
  
  broom.mixed::tidy(model, effects = "fixed", conf.int = TRUE) %>%
    filter(term == "host_substrateTerrestrial") %>%
    mutate(
      guild = guild,
      IRR  = exp(estimate),
      CI_low  = exp(conf.low),
      CI_high = exp(conf.high),
      significance = ifelse(p.value < 0.05, "Significant", "Not significant"))
})

forest_dat

# Prepare plotting data

forest_dat <- forest_dat %>%
  mutate(
    # enforce desired order
    guild = factor(guild, levels = broad_symbionts),
    
    # pretty labels in the same order
    guild_clean = levels(guild) %>%
      str_replace_all("_", " ") %>%
      str_to_sentence()) %>%
  mutate(direction = 
    if_else(IRR > 1, "Terrestrial > Epiphytic", "Epiphytic > Terrestrial"))

# Plotting

forest_plot_guilds <- ggplot(
  forest_dat,
  aes(x = IRR, y = guild, color = direction, alpha = significance)) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.6) +
  geom_errorbarh(aes(xmin = CI_low, xmax = CI_high), 
  height = 0.05, linewidth = 0.9) +
  geom_point(size = 3.5) +
  scale_x_log10(
    breaks = c(0.7, 1, 2.5, 5),
    labels = c("0.7","1","2.5","5")) +
  scale_y_discrete(labels = forest_dat$guild_clean) +
  scale_color_manual(values = c(
    "Terrestrial > Epiphytic" = "#D95F02",
    "Epiphytic > Terrestrial" = "#1B9E77")) +
  scale_alpha_manual(values = c("Significant" = 1, "Not significant" = 0.35)) +
  guides(alpha = "none") +
  labs(x = "Incidence rate ratio",
    y = NULL,
    color = NULL) +
  theme_classic(base_size = 15) +
  theme(legend.position = "top")

forest_plot_guilds


# THIRD PANEL: Only family level results ---------------------------------------
# from individual GLMMs of guild proportions -----------------------------------

# Consistent colors
family_cols <- c(
  "Araceae" = "#00CD6C",
  "Bromeliaceae" = "#009ADE",
  "Orchidaceae" = "#AF58BA",
  "Piperaceae" = "#FFC61E"
)

# Extract model outputs and prepare plotting data

forest_fam <- imap_dfr(fits_nb, function(model, guild) {
  broom.mixed::tidy(model, effects = "fixed", conf.int = TRUE) %>%
    filter(str_detect(term, "^host_family")) %>%
    mutate(
      guild  = guild,
      family = str_remove(term, "^host_family"),
      IRR    = exp(estimate),
      CI_low = exp(conf.low),
      CI_high= exp(conf.high))
}) %>%
  mutate(
    guild = factor(guild, levels = broad_symbionts),
    family = factor(family, levels = names(family_cols)),
    significance = if_else(p.value < 0.05, "Significant", "Not significant"),
    alpha_val    = if_else(p.value < 0.05, 1, 0.35))

label_map <- setNames(
  str_to_sentence(str_replace_all(broad_symbionts, "_", " ")),
  broad_symbionts
)

# Plotting

forest_plot_family_guilds <- ggplot( forest_fam, aes(
  x = IRR, 
  y = guild, 
  color = family) ) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.6) + 
  geom_errorbarh(aes(xmin = CI_low, 
  xmax = CI_high, alpha = alpha_val), 
  height = 0.18, linewidth = 0.9, 
  position = position_dodge(width = 0.8)) +
  geom_point( aes(alpha = alpha_val), size = 3.2, 
  position = position_dodge(width = 0.8)) +
  scale_x_log10(breaks = c(0.05,0.3, 1, 2.0),
        labels = c("0.05", "0.3","1","2.0")) +
  scale_y_discrete(labels = label_map) +
  scale_color_manual(values = family_cols, drop = FALSE) +
  scale_alpha_identity() + guides(alpha = "none") +
  labs( x = "Incidence rate ratio", y = NULL, color = NULL ) +
  theme_classic(base_size = 15) +
  theme(legend.position = "top") 

forest_plot_family_guilds


## Combine the three panels

forest_combined <- forest_plot + forest_plot_guilds + forest_plot_family_guilds 

forest_combined

# Save Figure 5 ----------------------------------------------------------------

ggsave(
  filename = "../Figures/Figure5.all_forest_plots1.png",
  plot = forest_combined,
  width = 14,
  height = 6,
  dpi = 600,
  units = "in",
  bg = "white"
)


