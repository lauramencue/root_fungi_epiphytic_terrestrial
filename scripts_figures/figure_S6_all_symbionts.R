# Project name: root_fungi_epiphytic_terrestrial
# Author: Laura Méndez
# Date: 15.03.26

# Goals: Figure SXX.Diversity profile of root-associated fungal α-diversity 
# across Hill orders

# Input: alpha_diversity_results.Rdata

# Setup ------------------------------------------------------------------------

rm(list = ls())

library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(scales)

# Load data --------------------------------------------------------------------

load("analyses3.Rdata")

fungi <- analyses$fungi
fungi$FunctionalGuilds

metadata <- analyses$metadata

asv_tax_fungi <- analyses$asv_tax_fungi

symbiotic_guilds <- c(
  "arbuscular_mycorrhizal",
  "ectomycorrhizal",
  "ericoid_mycorrhizal",
  "orchid_mycorrhizal",
  "ericoid_orchid_mycorrhizal",
  "root_endophyte",
  "dark_septate_endophyte",
  "fine_root_endophyte"
)

# keep only symbiotic guilds that are actually present
guild_levels <- fungi %>%
  filter(!is.na(FunctionalGuilds),
         FunctionalGuilds %in% symbiotic_guilds) %>%
  distinct(FunctionalGuilds) %>%
  pull(FunctionalGuilds)

# Dark2 palette
pal_fg <- setNames(
  brewer.pal(length(guild_levels), "Dark2"),
  guild_levels
)

guild_prop_fam_hab <- fungi %>%
  filter(!is.na(FunctionalGuilds),
         FunctionalGuilds %in% symbiotic_guilds) %>%
  mutate(
    host_substrate = factor(host_substrate,
                            levels = c("Epiphytic", "Terrestrial")),
    FunctionalGuilds = factor(FunctionalGuilds, levels = guild_levels),
    host_family = factor(host_family,
                         levels = c("Araceae", "Bromeliaceae", "Orchidaceae", "Piperaceae"))
  ) %>%
  distinct(ASV, host_family, host_substrate, FunctionalGuilds) %>%
  count(host_family, host_substrate, FunctionalGuilds, name = "n") %>%
  group_by(host_family, host_substrate) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

p_guild_family_host_substrate <- ggplot(
  guild_prop_fam_hab,
  aes(x = host_substrate, y = prop, fill = FunctionalGuilds)
) +
  geom_col(width = 0.7) +
  facet_wrap(~ host_family, ncol = 2) +
  scale_y_continuous(labels = percent_format()) +
  scale_fill_manual(values = pal_fg, drop = FALSE) +
  labs(
    y = "Proportion of symbiotic ASVs",
    x = NULL,
    fill = NULL
  ) +
  theme_classic(base_size = 15) +
  theme(
    axis.text.x = element_text(size = base_axis_size),
    axis.text.y = element_text(size = base_axis_size),
    axis.title.x = element_text(size = base_title_size),
    legend.position = "right",
    strip.text = element_text(size = base_axis_size, face = "bold"),
    plot.title = element_text(size = base_title_size, face = "bold"),
    plot.title.position = "plot"
  )

p_guild_family_host_substrate

ggsave(
  filename = "../Figures/FigureS6.all_symbionts.png",
  plot = last_plot(),
  width = 8,
  height = 8,
  dpi = 600,
  units = "in",
  bg = "white"
)
