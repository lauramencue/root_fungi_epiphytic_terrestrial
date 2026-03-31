


setwd("C:/Users/mendezcu_lu/Dropbox/Work/03_Marburg/02_epiphites_terrestrial_root_funghi/R/scripts")
rm(list = ls())

load("analyses.Rdata")
load("guild_composition_results.Rdata")

library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(RColorBrewer)
library(patchwork)


str(asv_tax_fungi)
#------------------------------------------------------------
# 1) Stacked bars: Identified vs Unidentified at Phylum
#------------------------------------------------------------

total_asvs <- nrow(asv_tax_fungi)

n_ident <- asv_tax_fungi %>%
  summarise(n = sum(!is.na(Phylum) & Phylum != "")) %>%
  pull(n)

n_unident <- total_asvs - n_ident

df_id <- tibble(
  Status = factor(c("Identified", "Unidentified"),
                  levels = c("Identified", "Unidentified")),
  n = c(n_ident, n_unident)
)

########### For the font size
base_axis_size   <- 13
base_legend_size <- 13
base_title_size  <- 13


p1 <- ggplot(df_id, aes(x = 1, y = n, fill = Status)) +
  geom_col(width = 0.75) +
  coord_flip() +
  scale_fill_manual(values = c("Identified" = "grey30", "Unidentified" = "grey80")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "Number of ASVs", fill = NULL,
       title = "Taxonomic assignment to phylum level") +
  theme_classic(base_size = 11) +
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


p1

#------------------------------------------------------------
# 2) Stacked bar: top 11 Phylum + Other (12 colors, Set3)
#------------------------------------------------------------

phy <- asv_tax_fungi %>%
  filter(!is.na(Phylum) & Phylum != "") %>%          # only identified
  count(Phylum, name = "n") %>%
  mutate(
    Phylum = if_else(n < 100, "Other", Phylum),
    Phylum = fct_reorder(Phylum, n, .desc = TRUE)
  ) %>%
  group_by(Phylum) %>%
  summarise(n = sum(n), .groups = "drop") %>%
  mutate(Phylum = fct_reorder(Phylum, n, .desc = TRUE))

# make sure palette has enough colors (Set3 max = 12)
n_cat <- nlevels(phy$Phylum)
pal <- RColorBrewer::brewer.pal(8, "Set1")



p2 <- ggplot(phy, aes(x = 1, y = n, fill = Phylum)) +
  geom_col(width = 0.75) +
  coord_flip() +
  scale_fill_manual(values = pal) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "Number of ASVs", fill = NULL,
       title = "Most common identified phyla (> 100 ASVs)") +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x = element_text(size = base_axis_size),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.x = element_text(size = base_title_size),
    
    legend.text = element_text(size = base_legend_size),
    legend.title = element_blank(),
    legend.position = "bottom",
    
    plot.title = element_text(size = base_title_size, face = "bold"),
    plot.title.position = "plot"
  )


p2

p1 / p2 + plot_layout(heights = c(1.2, 1))

############################### Functional guilds

# Only ASVs with Phylum identified
df_phylum_ident <- asv_tax_fungi %>%
  filter(!is.na(Phylum) & Phylum != "")

total_phylum_ident <- nrow(df_phylum_ident)

# Of those, how many have FunctionalGuilds identified?
n_ident_fg <- df_phylum_ident %>%
  summarise(n = sum(!is.na(FunctionalGuilds) & FunctionalGuilds != "")) %>%
  pull(n)

n_unident_fg <- total_phylum_ident - n_ident_fg

#----------------------------
# Plot 1: Identified vs Unidentified Functional Guilds
#----------------------------
df_id_fg <- tibble(
  Status = factor(c("Functional guild available", "No functional guild available"),
                  levels = c("Functional guild available", "No functional guild available")),
  n = c(n_ident_fg, n_unident_fg)
)

p1_fg <- ggplot(df_id_fg, aes(x = 1, y = n, fill = Status)) +
  geom_col(width = 0.75) +
  coord_flip() +
  scale_fill_manual(values = c("Functional guild available" = "grey30", "No functional guild available" = "grey80")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  labs(
    x = NULL, y = "Number of ASVs", fill = NULL,
    title = "Functional guild assignment"
  ) +
  theme_classic(base_size = 11) +
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

#----------------------------
# Plot 2: Identified guilds only, >=100 else Other
#----------------------------
guild <- asv_tax_fungi %>%
  filter(!is.na(FunctionalGuilds) & FunctionalGuilds != "") %>%  # only identified
  count(FunctionalGuilds, name = "n") %>%
  #mutate(
   # FunctionalGuilds = if_else(n < 100, "other", FunctionalGuilds)
  #) %>%
  group_by(FunctionalGuilds) %>%
  summarise(n = sum(n), .groups = "drop") %>%
  mutate(FunctionalGuilds = fct_reorder(FunctionalGuilds, n, .desc = TRUE))

unique(guild$FunctionalGuilds)

# Set3 has max 12 colors
n_cat <- nlevels(guild$FunctionalGuilds)
guild_levels <- c(
  "wood_litter_soil_saprotroph",
  "plant_pathogen",
  "foliar_endophyte",
  "dark_septate_endophyte",
  "animal-associated",
  "ectomycorrhizal",
  "root_endophyte_endomycorrhizal",
  "other",
  "arbuscular_mycorrhizal",
  "ericoid_mycorrhizal",
  "ericoid_orchid_mycorrhizal",
  "orchid_mycorrhizal")

guild$FunctionalGuilds <- factor(
  guild$FunctionalGuilds,
  levels = guild_levels
)

pal_fg <- setNames(
  RColorBrewer::brewer.pal(12, "Set3"),
  guild_levels
)


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
    
    legend.text = element_text(size = base_legend_size),
    legend.title = element_blank(),
    legend.position = "bottom",
    
    plot.title = element_text(size = base_title_size, face = "bold"),
    plot.title.position = "plot"
  )

p2_fg
#----------------------------
# Combine
#----------------------------
p1_fg / p2_fg + plot_layout(heights = c(1, 1.2))



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
  distinct(SampleID, ASV, FunctionalGuilds,
           host_family, name_analyses, host_habitat) %>%
  mutate(
    habitat = recode(host_habitat,
                     "E" = "Epiphytic",
                     "T" = "Terrestrial"),
    is_symbiotic = FunctionalGuilds %in% symbiotic_guilds
  ) %>%
  group_by(SampleID, habitat, host_family, name_analyses) %>%
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
  select(SampleID, host_family, habitat, prop_symbiotic, non_symbiotic) %>%
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


p_symb <- ggplot(plot_dat,
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
      "Symbiotic" = "grey30",
      "Non-symbiotic" = "grey80"
    )
  ) +
  facet_wrap(~host_family) +
  labs(
    x = " ",
    y = "Mean proportion of fungal ASVs",
    fill = "Guild type"
  ) +
  theme(legend.position = "bottom") +
  theme_classic(base_size = 13)

p_symb



guild_composition_results$guild_roots_family_form$FunctionalGuilds <- factor(
  guild_composition_results$guild_roots_family_form$FunctionalGuilds,
  levels = guild_levels
)

unique(guild_composition_results$guild_roots_family_form$FunctionalGuilds)

p_symb_fg <- ggplot(
  guild_composition_results$guild_roots_family_form,
  aes(`epiphitic/terrestrial`, prop, fill = FunctionalGuilds)
) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = pal_fg, drop = FALSE, guide = "none") +
  scale_x_discrete(
    labels = c("E" = "Epiphytic", "T" = "Terrestrial")) +
  facet_wrap(~host_family) +
  scale_y_continuous(labels = scales::percent) +
  theme_classic() +
  labs(
    y = " ",
    x = " ",
    fill = "Functional guild"
  ) +
  theme_classic(base_size = 13)

p_symb_fg



top_block <- p1 / p2 / p1_fg / p2_fg
bottom_block <- p_symb | p_symb_fg

p_2_2 <- wrap_elements(full = top_block) /
  wrap_elements(full = bottom_block) +
  plot_layout(heights = c(4, 3))

p_2_2


sym_p <- p_symb | p_symb_fg

sym_p

#------------------------------------------------------------
# Define the "starting set" for this step:
# only ASVs with a functional guild assigned
# (OPTIONAL: also require Phylum identified for full consistency)
#------------------------------------------------------------
df_fg <- asv_tax_fungi %>%
  filter(!is.na(FunctionalGuilds) & FunctionalGuilds != "")
# %>% filter(!is.na(Phylum) & Phylum != "")   # <- uncomment if you want that too

total_fg <- nrow(df_fg)

#------------------------------------------------------------
# Define which guild labels count as "symbiotic"
# (edit this vector to match your exact FungalTraits labels)
#------------------------------------------------------------
symbiotic_levels <- c(
  "arbuscular_mycorrhizal",
  "ectomycorrhizal",
  "ericoid_mycorrhizal",
  "orchid_mycorrhizal",
  "root_endophyte",
  "dark_septate_endophyte",
  "ericoid_orchid_mycorrhizal"
)

#------------------------------------------------------------
# Panel 1: Symbiotic vs Non-symbiotic (within guild-assigned ASVs)
#------------------------------------------------------------
df_symbio <- df_fg %>%
  mutate(
    Symbiosis = if_else(FunctionalGuilds %in% symbiotic_levels,
                        "Symbiotic", "Non-symbiotic")
  ) %>%
  count(Symbiosis, name = "n") %>%
  mutate(Symbiosis = factor(Symbiosis, levels = c("Symbiotic", "Non-symbiotic")))

p_symbio_1 <- ggplot(df_symbio, aes(x = 1, y = n, fill = Symbiosis)) +
  geom_col(width = 0.75) +
  coord_flip() +
  scale_fill_manual(values = c("Symbiotic" = "grey30", "Non-symbiotic" = "grey80")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  labs(
    x = NULL, y = "Number of ASVs", fill = NULL,
    title = "Symbiotic vs non-symbiotic functional guilds"
  ) +
  theme_classic(base_size = 11) +
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


p_symbio_1
#------------------------------------------------------------
# Panel 2: Composition of symbiotic guilds (≥ 100 ASVs else Other)
#------------------------------------------------------------
df_sym_guilds <- df_fg %>%
  filter(FunctionalGuilds %in% symbiotic_levels) %>%
  count(FunctionalGuilds, name = "n") %>%
  group_by(FunctionalGuilds) %>%
  summarise(n = sum(n), .groups = "drop") %>%
  mutate(FunctionalGuilds = fct_reorder(FunctionalGuilds, n, .desc = TRUE))

# Set3 max 12 colors
n_cat <- nlevels(df_sym_guilds$FunctionalGuilds)
pal <- RColorBrewer::brewer.pal(7, "Dark2")

p_symbio_2 <- ggplot(df_sym_guilds, aes(x = 1, y = n, fill = FunctionalGuilds)) +
  geom_col(width = 0.75) +
  coord_flip() +
  scale_fill_manual(values = pal) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  labs(
    x = NULL, y = "Number of ASVs", fill = NULL,
    title = "Symbiotic guild composition"
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x = element_text(size = base_axis_size),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.x = element_text(size = base_title_size),
    
    legend.text = element_text(size = base_legend_size),
    legend.title = element_blank(),
    legend.position = "bottom",
    
    plot.title = element_text(size = base_title_size, face = "bold"),
    plot.title.position = "plot"
  )

p_symbio_2

p_all_vertical <- 
  p1 /
  p2 /
  p1_fg /
  p2_fg /
  p_symbio_1 /
  p_symbio_2 +
  plot_layout(heights = c(1, 1.2, 1, 1.2, 1, 1.2))

p_all_vertical

################################################################################
# Figure 2: Guild composition by family
################################################################################

levels_fg <- levels(df_sym_guilds$FunctionalGuilds)

pal <- setNames(
  RColorBrewer::brewer.pal(
    max(3, min(12, length(levels_fg))),
    "Dark2"
  ),
  levels_fg
)

guild_composition_results$guild_roots_family <-
  guild_composition_results$guild_roots_family %>%
  mutate(
    FunctionalGuilds = factor(
      FunctionalGuilds,
      levels = levels_fg
    )
  )

ggplot(
  guild_composition_results$guild_roots_family,
  aes(host_family, prop, fill = FunctionalGuilds)
) +
  geom_col(width = 0.7) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = pal, drop = FALSE) +
  theme_classic(base_size = 11) +
  labs(
    y = "Proportion of ASVs",
    x = "Host family",
    fill = "Functional guild",
    title = "Functional composition by host family"
  )

guild_composition_results$guild_roots_family_form <-
  guild_composition_results$guild_roots_family_form %>%
  mutate(
    FunctionalGuilds = factor(FunctionalGuilds, levels = levels_fg),
    `epiphitic/terrestrial` = factor(`epiphitic/terrestrial`, levels = c("T", "E"))
  )

ggplot(
  guild_composition_results$guild_roots_family_form,
  aes(`epiphitic/terrestrial`, prop, fill = FunctionalGuilds)
) +
  geom_col(width = 0.7) +
  facet_wrap(~host_family) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = pal, drop = FALSE) +
  theme_classic(base_size = 18) +
  labs(
    y = "Proportion of ASVs",
    x = "Growth form",
    fill = "Functional guild"
  )
