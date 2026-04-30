# Project name: root_fungi_epiphytic_terrestrial
# Author: Laura Méndez
# Date: 15.02.26

# Goals: Figure 2

# Input: alpha_diversity_results.Rdata

# Setup ------------------------------------------------------------------------

rm(list = ls())

library(tidyverse)
library(ggplot2)
library(stringr)
library(broom.mixed)
library(patchwork)

# Load data --------------------------------------------------------------------

load("alpha_diversity_results.Rdata")

m_q0 <- alpha_diversity_results$m_q0
m_q1 <- alpha_diversity_results$m_q1
m_q2 <- alpha_diversity_results$m_q2

# Define color and other schemes -----------------------------------------------

# Color palette for host families

pd <- position_dodge(width = 0.6)

cols <- c("q = 0 (ASV richness)" = "#C40F5B",
  "q = 1 (common ASVs)" = "#FD8D3C",
  "q = 2 (dominant ASVs)" = "#089099")

# Function to extract clean model outputs consistently -------------------------

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
        term == "host_substrateTerrestrial" ~ "Substrate - Terrestrial",
        term == "host_familyBromeliaceae" ~ "Family - Bromeliaceae",
        term == "host_familyOrchidaceae" ~ "Family - Orchidaceae",
        term == "host_familyPiperaceae" ~ "Family - Piperaceae",
        TRUE ~ term))
  }

# Extract model outputs

df_q0 <- get_forest_df(m_q0) %>% mutate(model = "q = 0 (ASV richness)")
df_q1 <- get_forest_df(m_q1) %>% mutate(model = "q = 1 (common ASVs)")
df_q2 <- get_forest_df(m_q2) %>% mutate(model = "q = 2 (dominant ASVs)")

# Levels for plotting

lvl <- c("q = 0 (ASV richness)", "q = 1 (common ASVs)", "q = 2 (dominant ASVs)")

# Combine model outputs ordered by lvl

df_all <- bind_rows(df_q0, df_q1, df_q2) %>%
  mutate(model = factor(model, levels = rev(lvl))) 

# Plot forest plot with all 3 Hill numbers -------------------------------------

p_all <- ggplot(df_all, aes(x = IRR, y = term_clean, color = model)) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.7) +
  geom_errorbarh(aes(xmin = CI_low, xmax = CI_high, alpha = alpha_sig),
                 position = pd, height = 0.15, linewidth = 0.9) +
  geom_point(aes(alpha = alpha_sig), position = pd, size = 2.6) +
  scale_alpha_identity() +
  scale_x_log10(
    breaks = c(0.3, 0.6, 1, 1.5, 2),
    labels = c("0.3", "0.6", "1", "1.5", "2")) +
  scale_color_manual(values = cols, breaks = lvl, limits = lvl) +
  guides(color = guide_legend(reverse = FALSE)) +
  labs(
    x = "Incidence rate ratio",
    y = NULL,
    color = "Hill number") +
  theme_classic(base_size = 15) +
  theme(legend.position = "right",
        legend.text = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.text.x  = element_text(size = 15),
        axis.text.y  = element_text(size = 15))

p_all

# Save Figure 2 ----------------------------------------------------------------

ggsave(
  filename = "../Figures/Figure2.forest_plot_ALL_hill.png",
  plot = p_all,
  width = 6,
  height = 7,
  dpi = 600,
  units = "in",
  bg = "white"
)
