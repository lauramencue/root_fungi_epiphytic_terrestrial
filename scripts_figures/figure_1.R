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

# Shape palette for habitat
form_shapes <- c("Epiphytic" = 24, "Terrestrial" = 21)

# Recode habitat labels
alpha_df <- alpha_df %>%
  mutate(host_substrate = recode(host_substrate,
                                 "E" = "Epiphytic",
                                 "T" = "Terrestrial"))

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
  scale_x_continuous(position = "top", labels = scales::comma) +
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

# Label formatter for italic species names ----------------------------------------

format_sci_label <- function(x) {
  x <- gsub("_", " ", x)
  
  sapply(x, function(z) {
    if (z == "Overall mean per habitat") return(z)
    
    parts <- strsplit(z, " +")[[1]]
    
    if (length(parts) >= 2) {
      if (grepl("^sp", parts[2])) {
        return(paste0("<i>", parts[1], "</i> ", paste(parts[-1], collapse = " ")))
      } else {
        return(paste0("<i>", parts[1], " ", parts[2], "</i>"))
      }
    }
    
    paste0("<i>", parts[1], "</i>")
  }, USE.NAMES = FALSE)
}

# Build manual Y axis ----------------------------------------------------------

family_levels <- c("Araceae", "Bromeliaceae", "Orchidaceae", "Piperaceae")

species_order <- alpha_df %>%
  mutate(host_family = factor(host_family, levels = family_levels)) %>%
  group_by(host_family, host_species_identity) %>%
  summarise(mean_val = mean(ASV_richness), .groups = "drop") %>%
  arrange(host_family, mean_val)

gap <- 1.2
current_y <- 0
y_map_list <- list()

for (fam in family_levels) {
  spp <- species_order %>% filter(host_family == fam) %>% pull(host_species_identity)
  
  rows <- c("Overall mean per habitat", spp)
  n <- length(rows)
  
  y_vals <- seq(current_y + n - 1, current_y, by = -1)
  
  y_map_list[[fam]] <- tibble(
    host_family = fam,
    row_id = rows,
    y = y_vals
  )
  
  current_y <- current_y + n + gap
}

y_map <- bind_rows(y_map_list)

# Merge data -------------------------------------------------------------------

alpha_plot <- alpha_df %>%
  left_join(y_map %>% filter(row_id != "Overall mean per habitat"),
            by = c("host_family", "host_species_identity" = "row_id"))

# Growth form summary with manual dodge offset
dodge_offset <- 0.18

gf_summary <- alpha_df %>%
  group_by(host_family, host_substrate) %>%
  summarise(mean = mean(ASV_richness), sd = sd(ASV_richness), .groups = "drop") %>%
  left_join(y_map %>% filter(row_id == "Overall mean per habitat"),
            by = "host_family") %>%
  mutate(y_dodged = ifelse(host_substrate == "Epiphytic", y + dodge_offset, y - dodge_offset))

# Family mean with y range for scoped dashed segments
x_range <- alpha_df %>%
  group_by(host_family) %>%
  summarise(xmin = min(ASV_richness, na.rm = TRUE),
            xmax = max(ASV_richness, na.rm = TRUE),
            .groups = "drop")

family_mean <- alpha_df %>%
  group_by(host_family) %>%
  summarise(mean_family = mean(ASV_richness, na.rm = TRUE), .groups = "drop") %>%
  left_join(y_map %>% filter(row_id == "Overall mean per habitat"), by = "host_family") %>%
  left_join(x_range, by = "host_family") %>%
  left_join(
    y_map %>% group_by(host_family) %>%
      summarise(ymin = min(y), ymax = max(y), .groups = "drop"),
    by = "host_family"
  )

# Axis labels
y_labels <- format_sci_label(y_map$row_id)

# BOTTOM PANEL -----------------------------------------------------------------

p_bottom <- ggplot(alpha_plot, aes(x = ASV_richness, y = y)) +
  
  # Dashed family mean — scoped per family
  geom_segment(data = family_mean,
               aes(x = mean_family, xend = mean_family,
                   y = ymin, yend = ymax),
               inherit.aes = FALSE,
               linetype = "dashed", linewidth = 0.9, color = "black") +
  
  # Lines connecting samples within species
  geom_line(aes(group = host_species_identity, color = host_family),
            linewidth = 0.9) +
  
  # Individual points
  geom_point(aes(shape = host_substrate),
             size = 3, alpha = 0.7,
             fill = "black", color = "black", stroke = 0.3) +
  
  # Error bars for growth form summary
  geom_segment(data = gf_summary,
               aes(x = mean - sd, xend = mean + sd,
                   y = y_dodged, yend = y_dodged),
               inherit.aes = FALSE,
               linewidth = 0.6, color = "black") +
  
  # Summary mean points
  geom_point(data = gf_summary,
             aes(x = mean, y = y_dodged,
                 shape = host_substrate,
                 fill = host_family),
             inherit.aes = FALSE,
             size = 4, color = "black", stroke = 0.8) +
  
  scale_shape_manual(name = "Host rooting substrate", values = form_shapes) +
  scale_fill_manual(name = "Host family", values = family_cols) +
  scale_color_manual(name = "Host family", values = family_cols) +
  scale_y_continuous(
    breaks = y_map$y,
    labels = y_labels
  ) +
  
  theme_classic() +
  labs(x = "ASV richness (α-diversity) Hill q = 0", y = NULL) +
  theme(
    axis.text.y = ggtext::element_markdown(size = 13),
    axis.text.x = element_text(size = 13),
    axis.title.x = element_text(size = 15),
    legend.position = "none"
  ) +
  scale_x_continuous(labels = scales::label_comma())

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