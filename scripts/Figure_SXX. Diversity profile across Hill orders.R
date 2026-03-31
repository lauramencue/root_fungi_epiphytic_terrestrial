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

load("alpha_diversity_results.Rdata")

comm <- alpha_diversity_results$comm

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

# Save Figure ------------------------------------------------------------------

ggsave(
  filename = "../Figures/FigureSx.diversity_profile.png",
  plot = diversity_profile,
  width = 6,
  height = 4,
  dpi = 600,
  units = "in",
  bg = "white"
)
