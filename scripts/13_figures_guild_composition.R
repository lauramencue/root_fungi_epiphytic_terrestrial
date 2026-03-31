# Project name: root_fungi_epiphytic_terrestrial
# Author: Laura Méndez
# Date: [Date]

# Goals:
# Create publication-quality figures for functional guild composition

# Figures created:
# - Fig 1: Guild composition by growth form (stacked bar)
# - Fig 2: Guild composition by family (stacked bar)
# - Fig 3: Guild composition by family × growth form (faceted)
# - Fig 4: Guild composition by species (detailed)
# - Fig 5: Indicator ASV guild composition by family×form (heatmap)

# Input: analyses.Rdata, guild_composition_results.Rdata, indicator_species_results.Rdata
# Dependencies: ggplot2, tidyverse, viridis

################################################################################
# Setup
################################################################################

rm(list = ls())
library(tidyverse)
library(ggplot2)
library(patchwork)
library(ggdist)
library(broom.mixed)

load("analyses2.Rdata")
load("guild_composition_results2.Rdata")


########### For the font size
base_axis_size   <- 15
base_legend_size <- 15
base_title_size  <- 15


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
      term == "host_habitatT" ~ "Substrate: Terrestrial",
      term == "host_familyBromeliaceae" ~ "Family: Bromeliaceae",
      term == "host_familyOrchidaceae" ~ "Family: Orchidaceae",
      term == "host_familyPiperaceae" ~ "Family: Piperaceae",
      TRUE ~ term
    ),
    term_clean = factor(term_clean, levels = rev(unique(term_clean)))
  )

forest_plot <- ggplot(coefs, aes(x = OR, y = term_clean)) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.6) +
  geom_errorbarh(aes(xmin = CI_low, xmax = CI_high, alpha = alpha_val),
                 height = 0.1, linewidth = 0.9) +
  geom_point(aes(alpha = alpha_val), size = 3) +
  scale_x_log10(
    breaks = c(0.5, 0.75, 1, 1.4, 2, 3),
    labels = c("0.5","0.75","1","1.4","2","3")
  ) +
  scale_alpha_identity() +
  guides(alpha = "none") +
  coord_cartesian(clip = "off") +
  labs(
    x = "Odds ratio (log scale)",
    y = NULL
  ) +
  theme_classic(base_size = 13) +
  theme(
    axis.title.x = element_text(margin = margin(t = 8))
  )

forest_plot

str(m_symb)
coefs <- broom.mixed::tidy(m_symb_int, effects = "fixed") %>%
  filter(term != "(Intercept)") %>%
  mutate(
    OR = exp(estimate),
    CI_low = exp(estimate - 1.96 * std.error),
    CI_high = exp(estimate + 1.96 * std.error),
    sig = p.value < 0.05,
    alpha_sig = ifelse(sig, 1, 0.35),
    # nicer labels
    term_clean = case_when(
      term == "host_habitatT" ~ "Habitat: Terrestrial",
      term == "host_familyBromeliaceae" ~ "Family: Bromeliaceae",
      term == "host_familyOrchidaceae" ~ "Family: Orchidaceae",
      term == "host_familyPiperaceae" ~ "Family: Piperaceae",
      term == "host_habitatT:host_familyBromeliaceae" ~ "Interaction: Terrestrial × Bromeliaceae",
      term == "host_habitatT:host_familyOrchidaceae" ~ "Interaction: Terrestrial × Orchidaceae",
      term == "host_habitatT:host_familyPiperaceae" ~ "Interaction: Terrestrial × Piperaceae",
      TRUE ~ term
    )
  ) %>%
  # order: main effects first, then interactions (or keep as you like)
  mutate(
    group = ifelse(str_detect(term, ":"), "Interaction", "Main"),
    term_clean = factor(term_clean, levels = rev(unique(term_clean)))
  )

forest_plot <- ggplot(coefs, aes(x = OR, y = term_clean)) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.7) +
  geom_pointrange(
    aes(xmin = CI_low, xmax = CI_high, alpha = alpha_sig),
    linewidth = 0.7
  ) +
    scale_x_log10(
    breaks = c(0.5, 0.75, 1, 1.4, 2, 3),
    labels = c("0.5","0.75","1","1.4","2","3")
  ) +
  scale_alpha_identity() +
  coord_cartesian(clip = "off") +
  labs(
    x = "Odds ratio (log scale)",
    y = NULL) +
  theme_classic(base_size = 15) +
  theme(
    axis.title.x = element_text(margin = margin(t = 8)))

forest_plot


ggsave(
  filename = "../Figures/Figure3.prop_symbionts_model.png",
  plot = forest_plot,
  width = 5,
  height = 5,
  dpi = 600,
  units = "in",
  bg = "white"
)



############## Dot plots with mean of proportions of symbionts for habitat x family groups


symb_sample_prop <- symb_sample_prop %>%
  mutate(host_habitat = recode(host_habitat,
                               "E" = "Epiphytic", "T" = "Terrestrial")) 


p_symb_raw <- ggplot(symb_sample_prop,
  aes(x = host_family, y = prop_symb)) +
  # Raw points
  geom_jitter(width = 0.12, size = 2, alpha = 0.6, color = "grey40",) +
  
  # Mean ± 95% CI
  stat_summary(
    fun.data = mean_cl_normal,
    geom = "pointrange",
    color = "black",
    size = 0.8
  ) +
  
  #facet_wrap(~host_family) +
  
  scale_y_continuous(labels = scales::percent_format(),
                     limits = c(0,1)) +
  
  labs(
    x = NULL,
    y = "Proportion of symbiotic ASVs"
  ) +
  
  theme_classic(base_size = 15) +
  theme(
    legend.position = "none")


p_symb_raw


#ggsave(
  filename = "../Figures/Figure3.prop_symbionts_raw_dot_plot.png",
  plot = p_symb_raw,
  width = 5,
  height = 8,
  dpi = 600,
  units = "in",
  bg = "white"
)


plot_symb_prop <- p_symb_raw | forest_plot + plot_layout(ncol = 2, widths = c(3, 2))

plot_symb_prop


#ggsave(
  filename = "../Figures/Figure3.dots_forest.svg",
  plot = plot_symb_prop,
  width = 13,
  height = 7,
  dpi = 600,
  units = "in",
  bg = "white"
)


############################### Functional guilds

# Total ASVs (all rows = all ASVs in your ASV taxonomy table)
total_asvs <- nrow(asv_tax_fungi)

# How many have FunctionalGuilds assigned?
n_ident_fg <- asv_tax_fungi %>%
  summarise(n = sum(!is.na(FunctionalGuilds) & FunctionalGuilds != "")) %>%
  pull(n)

n_unident_fg <- total_asvs - n_ident_fg

df_id_fg <- tibble(
  Status = factor(
    c("Functional guild available", "No functional guild available"),
    levels = c("Functional guild available", "No functional guild available")
  ),
  n = c(n_ident_fg, n_unident_fg)
)

p1_fg <- ggplot(df_id_fg, aes(x = 1, y = n, fill = Status)) +
  geom_col(width = 0.75) +
  coord_flip() +
  scale_fill_manual(values = c(
    "Functional guild available" = "grey30",
    "No functional guild available" = "grey80"
  )) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  labs(
    x = NULL,
    y = "Number of ASVs",
    fill = NULL,
    title = "Functional guild assignment"
  ) +
  theme_classic(base_size = 15) +
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

p1_fg






# Set3 has max 12 colors

guild_levels <- rev(c(
  "wood_litter_soil_saprotroph",
  "plant_pathogen",
  "foliar_endophyte",
  "animal-associated",
  "other",
  "dark_septate_endophyte",
  "ectomycorrhizal",
  "root_endophyte",
  "arbuscular_mycorrhizal",
  "ericoid_mycorrhizal",
  "ericoid_orchid_mycorrhizal",
  "orchid_mycorrhizal",
  "fine_root_endophyte"
))

pal_fg <- c(
  "wood_litter_soil_saprotroph" = "#8DD3C7",
  "plant_pathogen" = "#FFFFB3",
  "foliar_endophyte" = "#BEBADA",
  "animal-associated" = "#80B1D3",
  "other" = "#FCCDE5",
  "dark_septate_endophyte" = "#FB8072",
  "ectomycorrhizal" = "#FDB462",
  "root_endophyte" = "#B3DE69",
  "fine_root_endophyte" = "lightpink3",
  "arbuscular_mycorrhizal" = "#BC80BD",
  "ericoid_mycorrhizal" = "#D9D9D9",
  "ericoid_orchid_mycorrhizal" = "#FFED6F",
  "orchid_mycorrhizal" = "#CCEBC5"
)

guild <- asv_tax_fungi %>%
  filter(!is.na(FunctionalGuilds) & FunctionalGuilds != "") %>%  # only identified
  count(FunctionalGuilds, name = "n") %>%
  #mutate(
  # FunctionalGuilds = if_else(n < 100, "other", FunctionalGuilds)
  #) %>%
  group_by(FunctionalGuilds) %>%
  summarise(n = sum(n), .groups = "drop") %>%
  mutate(FunctionalGuilds = fct_reorder(FunctionalGuilds, n, .desc = TRUE))

guild$FunctionalGuilds <- factor(
  guild$FunctionalGuilds,
  levels = guild_levels
)


########### For the font size
base_axis_size   <- 15
base_legend_size <- 15
base_title_size  <- 15

p2_fg <- ggplot(guild, aes(x = 1, y = n, fill = FunctionalGuilds)) +
  geom_col(width = 0.75) +
  coord_flip() +
  scale_fill_manual(values = pal_fg, drop = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "Number of ASVs", fill = NULL,
       title = "Most common identified functional guilds") +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x = element_text(size = base_axis_size),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.x = element_text(size = base_title_size),
    legend.position = "none",
    plot.title = element_text(size = base_title_size, face = "bold"),
    plot.title.position = "plot"
  )

p2_fg

guild_prop <- fungi_symb %>%
  mutate(
    habitat = recode(host_habitat, "E" = "Epiphytic", "T" = "Terrestrial"),
    FunctionalGuilds = factor(FunctionalGuilds, levels = guild_levels)
  ) %>%
  distinct(ASV, habitat, FunctionalGuilds) %>%
  count(habitat, FunctionalGuilds) %>%
  group_by(habitat) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

# 2) Plot with EXACT same palette as p2_fg
p_guild <- ggplot(guild_prop,
                  aes(x = habitat, y = prop, fill = FunctionalGuilds)) +
  coord_flip() +
  geom_col(width = 0.7) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = pal_fg, drop = FALSE) +
  #theme_classic(base_size = 11) +
  labs(
    y = "Proportion of symbiotic ASVs",
    x = NULL,
    fill = NULL,
    title = "Symbiotic functional guild composition by habitat"
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x = element_text(size = base_axis_size),
    axis.text.y = element_text(size = base_axis_size),
    axis.title.x = element_text(size = base_title_size),
    legend.position = "none",
    plot.title = element_text(size = base_title_size, face = "bold"),
    plot.title.position = "plot"
  )

p_guild



guild_combined <- p1_fg / p2_fg / p_guild + plot_layout(heights = c(1, 1, 3))

guild_combined


ggsave(
  filename = "../Figures/Figure4.barplots_guild.svg",
  plot = guild_combined,
  width = 13,
  height = 7,
  dpi = 600,
  units = "in",
  bg = "white"
)


str(fungi_symb)

library(scales)
#library(forcats)

# 1) Summarise to proportions per group
guild_prop_group <- fungi_symb %>%
  filter(!is.na(FunctionalGuilds) & FunctionalGuilds != "") %>%
  count(group, FunctionalGuilds, name = "n_asv") %>%
  group_by(group) %>%
  mutate(prop = n_asv / sum(n_asv)) %>%
  ungroup() %>%
mutate(FunctionalGuilds = factor(FunctionalGuilds, levels = guild_levels))


# OPTIONAL: order the 8 groups nicely (E then T within each family)
guild_prop_group <- guild_prop_group %>%
  mutate(
    group = factor(group, levels = c(
      "Araceae_E","Araceae_T",
      "Bromeliaceae_E","Bromeliaceae_T",
      "Orchidaceae_E","Orchidaceae_T",
      "Piperaceae_E","Piperaceae_T"
    ))
  )

# 2) Plot: vertical stacked bars
ggplot(guild_prop_group, aes(x = group, y = prop, fill = FunctionalGuilds)) +
  geom_col(width = 0.7) +
  scale_y_continuous(labels = percent_format()) +
  scale_fill_manual(values = pal_fg, drop = FALSE) +
  labs(
    y = "Proportion of symbiotic ASVs",
    x = NULL,
    fill = NULL,
    title = "Symbiotic functional guild composition by group"
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.position = "top",
    plot.title = element_text(face = "bold"),
    plot.title.position = "plot"
  )


# Forest plot of guild models

forest_dat <- imap_dfr(fits_nb, function(model, guild) {
  
  broom.mixed::tidy(model, effects = "fixed", conf.int = TRUE) %>%
    filter(term == "host_habitatT") %>%
    mutate(
      guild = guild,
      OR  = exp(estimate),
      CI_low  = exp(conf.low),
      CI_high = exp(conf.high)
    )
})

forest_dat

guild_order <- c(
  "Ericoid/orchid mycorrhizal",
  "Orchid mycorrhizal",
  "Dark septate endophyte",
  "Root endophyte",
  "Ectomycorrhizal",
  "Ericoid mycorrhizal",
  "Arbuscular mycorrhizal"
)


forest_dat <- forest_dat %>%
  mutate(
    guild_clean = str_replace_all(guild, "_", " ") %>%
      str_to_sentence() %>%
      str_replace("^Ericoid orchid mycorrhizal$", "Ericoid/orchid mycorrhizal") %>%
      str_replace("^Root endophyte$", "Root endophyte"),
    direction = ifelse(OR > 1, "Terrestrial > Epiphytic", "Epiphytic > Terrestrial"),
    guild_clean = factor(guild_clean, levels = guild_order),
    significance = ifelse(p.value < 0.05, "Significant", "Not significant")
  )



forest_plot_guilds <- ggplot(forest_dat, aes(x = IRR,
                                             y = guild_clean,
                                             color = direction,
                                             alpha = significance)) +
  
  geom_vline(xintercept = 1, linetype = "dashed", size = 0.6) +
  
  geom_errorbarh(aes(xmin = CI_low, xmax = CI_high),
                 height = 0.2, size = 0.9) +
  
  geom_point(size = 3.5) +
  
  scale_x_log10() +
  
  scale_color_manual(values = c(
    "Terrestrial > Epiphytic" = "#D95F02",
    "Epiphytic > Terrestrial" = "#1B9E77"
  )) +
  
  scale_alpha_manual(values = c(
    "Significant" = 1,
    "Not significant" = 0.35
  )) +
  guides(alpha = "none") +
  labs(x = "Incidence rate ratio (Terrestrial vs Epiphytic)",
       y = NULL,
       color = NULL,
       alpha = NULL) +
  
  theme_classic(base_size = 13) +
  theme(
    legend.position = "top",
    axis.text.y = element_text(size = 11)
  )

forest_plot_guilds

ggsave(
  filename = "../Figures/Figure5.forest_plot_guild.png",
  plot = forest_plot_guilds,
  width = 6,
  height = 6,
  dpi = 600,
  units = "in",
  bg = "white"
)



# ----------------------------
# 2) FAMILY forest plot
# ----------------------------
family_cols <- c(
  "Araceae" = "#00CD6C",
  "Bromeliaceae" = "#009ADE",
  "Orchidaceae" = "#AF58BA",
  "Piperaceae" = "#FFC61E"
)

# Use your model list name here (you used fits_nb before)
# fits_nb should be a named list: names are guilds, values are glmmTMB models

forest_fam <- imap_dfr(fits_nb, function(model, guild) {
  broom.mixed::tidy(model, effects = "fixed", conf.int = TRUE) %>%
    filter(str_detect(term, "^host_family")) %>%
    mutate(
      guild = guild,
      family = str_remove(term, "^host_family"),
      IRR = exp(estimate),
      CI_low  = exp(conf.low),
      CI_high = exp(conf.high)
    )
})

forest_fam <- forest_fam %>%
  mutate(
    # clean guild names + fix the two special cases
    guild_clean = str_replace_all(guild, "_", " ") %>%
      str_to_sentence() %>%
      str_replace("^Ericoid orchid mycorrhizal$", "Ericoid/orchid mycorrhizal") %>%
      str_replace("^Root endophyte$", "Root endophyte"),
    guild_clean = factor(guild_clean, levels = guild_order),
    significance = ifelse(p.value < 0.05, "Significant", "Not significant"),
    alpha_val = ifelse(significance == "Significant", 1, 0.35)
  )


# Grouped (dodged) forest plot: 3 points per guild (one per family coefficient)
forest_plot_family_guilds <- ggplot(
  forest_fam,
  aes(x = IRR, y = guild_clean, color = family)
) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.6) +
  geom_errorbarh(
    aes(xmin = CI_low, xmax = CI_high, alpha = alpha_val),
    height = 0.18,
    linewidth = 0.9,
    position = position_dodge(width = 0.65)
  ) +
  geom_point(
    aes(alpha = alpha_val),
    size = 3.2,
    position = position_dodge(width = 0.65)
  ) +
  scale_x_log10() +
  scale_color_manual(values = family_cols, drop = FALSE) +
  scale_alpha_identity() +
  guides(alpha = "none") +
  labs(
    x = "Incidence rate ratio (Families vs Araceae)",
    y = NULL,
    color = NULL
  ) +
  theme_classic(base_size = 13) +
  theme(legend.position = "top")

forest_plot_family_guilds


ggsave(
  filename = "../Figures/Figure5.forest_plot_family_guild.png",
  plot = forest_plot_family_guilds,
  width = 5,
  height = 6,
  dpi = 600,
  units = "in",
  bg = "white"
)



# Extract habitat term from each model















#

p3 <- ggplot(
  symb_family_hab_prop,
  aes(host_habitat, prop, fill = FunctionalGuilds)
) +
  geom_col(width = 0.7) +
  scale_fill_brewer(palette = "Set3") +
  facet_wrap(~host_family) +
  scale_y_continuous(labels = scales::percent) +
  theme_classic() +
  labs(
    y = "Proportion of ASVs",
    x = "Growth form",
    fill = "Functional guild"
  )

p3




ggplot(plot_dat,
       aes(x = host_habitat,
           y = proportion,
           color = name_analyses)) +
  geom_jitter(width = 0.1, size = 2) +
  facet_wrap(~host_family) +
  theme_classic()


ggplot(symb_sample_prop,
       aes(x = prop_symb,
           y = total_asv)) +
  geom_point() +
  facet_wrap(~host_family) +
  theme_classic()

library(diptest)
dip.test(plot_dat$proportion)

plot_dat %>%
  group_by(host_family, host_habitat) %>%
  summarise(p_dip = dip.test(proportion)$p.value)


pecies_means <- plot_dat %>%
  group_by(name_analyses) %>%
  summarise(mean_prop = mean(proportion))

library(mclust)

species_means <- plot_dat %>%
  group_by(name_analyses) %>%
  summarise(mean_prop = mean(proportion))
Mclust(species_means$mean_prop, G = 1:3)


colnames(fungi)


coords <- fungi %>%
  select(SampleID, latitude, longitude) %>%
  distinct()

plot_dat <- plot_dat %>%
  left_join(coords, by = "SampleID")

library(ggplot2)

ggplot(plot_dat,
       aes(x = longitude,
           y = latitude,
           color = proportion)) +
  geom_point(size = 3, alpha = 0.6) +
  scale_color_viridis_c() +
  theme_classic()


plot_dat <- plot_dat %>%
  mutate(cluster = ifelse(proportion < 0.45, "Low", "High"))

ggplot(plot_dat,
       aes(x = longitude,
           y = latitude,
           color = cluster, alpha = 0.6)) +
  geom_point(size = 3) +
  theme_classic()


library(spdep)

coords <- cbind(plot_dat$longitude, plot_dat$latitude)

# Distance-based neighbours
nb <- dnearneigh(coords, 0, 500)  # adjust max distance (meters)
lw <- nb2listw(nb, style = "W")

moran.test(plot_dat$proportion, lw)



tot0 <- median(symb_ind_guild$total_symb_asv, na.rm = TRUE)
str(tot0)
emm_means_tbl <- imap_dfr(fits_nb, function(m, g) {
  
  # predicted means for all family × habitat combinations
  emm <- emmeans(
    m,
    ~ host_family * host_habitat,
    offset = log(tot0),     # sets the offset used for predictions
    type = "response"       # back-transform to response scale
  )
  
  as.data.frame(emm) %>%
    mutate(
      guild = g,
      guild_lab = str_replace_all(g, "_", " "),
      total_symb_asv_ref = tot0
    ) %>%
    rename(
      mean = response,
      SE = SE,
      CI_low = asymp.LCL,
      CI_high = asymp.UCL
    ) %>%
    select(guild, guild_lab, host_family, host_habitat, mean, SE, CI_low, CI_high, total_symb_asv_ref)
})

emm_means_tbl

contrast_tbl <- imap_dfr(fits_nb, function(m, g) {
  
  emm <- emmeans(
    m,
    ~ host_habitat | host_family,
    offset = log(tot0),
    type = "link"   # force link scale so columns are consistent
  )
  
  con <- contrast(emm, method = "revpairwise")  # T vs E within family on link scale
  
  as.data.frame(con) %>%
    mutate(
      guild = g,
      guild_lab = str_replace_all(g, "_", " "),
      total_symb_asv_ref = tot0,
      # link-scale estimate is log(RR) for log-link models
      RR = exp(estimate),
      RR_low = exp(estimate - 1.96 * SE),
      RR_high = exp(estimate + 1.96 * SE)
    ) %>%
    select(guild, guild_lab, host_family, contrast, RR, RR_low, RR_high, estimate, SE, df, p.value, total_symb_asv_ref)
})

contrast_tbl

plot_contrast <- contrast_tbl %>%
  mutate(
    sig = p.value < 0.05,
    alpha_v = ifelse(sig, 1, 0.25)
  )

ggplot(plot_contrast,
       aes(x = guild_lab, y = RR, ymin = RR_low, ymax = RR_high, alpha = alpha_v)) +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.6) +
  geom_pointrange(linewidth = 0.7, fatten = 1.6) +
  scale_y_log10() +
  scale_alpha_identity() +
  coord_flip() +
  facet_wrap(~ host_family, ncol = 2) +
  labs(
    x = NULL,
    y = "Rate ratio (Terrestrial / Epiphytic)",
    title = "Habitat effect on symbiotic guild richness",
    subtitle = paste0("Additive NB GLMMs; non-significant contrasts shown with higher transparency. Offset set to total_symb_asv = ", tot0, ".")
  ) +
  theme_classic(base_size = 13)

str(asv_tax_fungi)
length(unique(asv_tax_fungi$Phylum))




symbiotic_guilds <- c(
  "arbuscular_mycorrhizal",
  "ectomycorrhizal",
  "ericoid_mycorrhizal",
  "orchid_mycorrhizal",
  "ericoid_orchid_mycorrhizal",
  "root_endophyte_endomycorrhizal",
  "dark_septate_endophyte"
)

symb_prop <- fungi %>%
  filter(!is.na(FunctionalGuilds)) %>%
  distinct(clientId, ASV, FunctionalGuilds,
           host_family, name_analyses, `epiphitic/terrestrial`) %>%
  mutate(
    habitat = recode(`epiphitic/terrestrial`,
                     "E" = "Epiphytic",
                     "T" = "Terrestrial"),
    is_symbiotic = FunctionalGuilds %in% symbiotic_guilds
  ) %>%
  group_by(clientId, habitat, host_family, name_analyses) %>%
  summarise(
    total_asv = n_distinct(ASV),
    symbiotic_asv = n_distinct(ASV[is_symbiotic]),
    prop_symbiotic = symbiotic_asv / total_asv,
    .groups = "drop"
  )

plot_dat <- symb_prop %>%
  mutate(
    non_symbiotic = 1 - prop_symbiotic
  ) %>%
  select(clientId, host_family, habitat, prop_symbiotic, non_symbiotic) %>%
  pivot_longer(
    cols = c(prop_symbiotic, non_symbiotic),
    names_to = "guild_type",
    values_to = "proportion"
  ) %>%
  mutate(
    guild_type = recode(
      guild_type,
      prop_symbiotic = "Symbiotic",
      non_symbiotic  = "Non-symbiotic"
    )
  )

library(ggplot2)

ggplot(plot_dat,
       aes(x = habitat, y = proportion, fill = guild_type)) +
  stat_summary(
    fun = mean,
    geom = "bar",
    position = "fill",
    width = 0.6
  ) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(
    values = c(
      "Symbiotic" = "#1b9e77",
      "Non-symbiotic" = "grey70"
    )
  ) +
  facet_wrap(~host_family) +
  labs(
    x = "Rooting habitat",
    y = "Mean proportion of fungal ASVs",
    fill = "Guild type"
  ) +
  theme_classic(base_size = 13)


################################################################################
# Figure 1: Guild composition by growth form (root-associated)
################################################################################

p1 <- ggplot(
  guild_composition_results$guild_roots_form,
  aes(`epiphitic/terrestrial`, prop, fill = FunctionalGuilds)
) +
  geom_col(width = 0.7) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_brewer(palette = "Set3") +
  theme_classic() +
  labs(
    y = "Proportion of ASVs",
    x = "Growth form",
    fill = "Functional guild",
    title = "Functional composition by growth form"
  )


p1


ggsave("Fig_guild_by_form.pdf", p1, width = 7, height = 5)

################################################################################
# Figure 2: Guild composition by family
################################################################################

p2 <- ggplot(
  guild_composition_results$guild_roots_family,
  aes(host_family, prop, fill = FunctionalGuilds)
) +
  geom_col(width = 0.7) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_brewer(palette = "Set3") +
  theme_classic() +
  labs(
    y = "Proportion of ASVs",
    x = "Host family",
    fill = "Functional guild",
    title = "Functional composition by host family"
  )

p2

ggsave("Fig_guild_by_family.pdf", p2, width = 8, height = 5)

################################################################################
# Figure 3: Guild composition by family × growth form (faceted)
################################################################################

p3 <- ggplot(
  guild_composition_results$guild_roots_family_form,
  aes(`epiphitic/terrestrial`, prop, fill = FunctionalGuilds)
) +
  geom_col(width = 0.7) +
  scale_fill_brewer(palette = "Set3") +
  facet_wrap(~host_family) +
  scale_y_continuous(labels = scales::percent) +
  theme_classic() +
  labs(
    y = "Proportion of ASVs",
    x = "Growth form",
    fill = "Functional guild"
  )

p3

ggsave("Fig_guild_by_family_form.pdf", p3, width = 10, height = 6)

################################################################################
# Figure 4: Indicator ASV guild composition (heatmap)
################################################################################

# Prepare heatmap data
hm_df <- indicator_species_results$indicators_family_form_annotated %>%
  filter(!is.na(indicator_for), !is.na(FunctionalGuilds)) %>%
  count(indicator_for, FunctionalGuilds, name = "n")

p4 <- ggplot(hm_df, aes(x = indicator_for, y = FunctionalGuilds, fill = n)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = n), size = 3) +
  scale_fill_viridis_c(
    option = "inferno",
    direction = 1,
    begin = 0.35,
    end = 1,
    name = "Indicator ASVs"
  ) +
  scale_x_discrete(drop = FALSE) +
  scale_y_discrete(drop = FALSE) +
  labs(x = NULL, y = NULL, fill = "Indicator ASVs") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.ticks = element_blank()
  )


p4

ggsave("Fig_indicator_guild_heatmap.pdf", p4, width = 10, height = 6)

cat("\n✓ Guild composition figures created!\n")
cat("Figures saved:\n")
cat("  - Fig_guild_by_form.pdf\n")
cat("  - Fig_guild_by_family.pdf\n")
cat("  - Fig_guild_by_family_form.pdf\n")
cat("  - Fig_indicator_guild_heatmap.pdf\n")
