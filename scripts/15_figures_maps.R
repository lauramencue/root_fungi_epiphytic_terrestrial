# Project name: root_fungi_epiphytic_terrestrial
# Author: Laura Méndez
# Date: [Date]

# Goals:
# Create publication-quality maps showing sampling localities

# Figures created:
# - Fig 1: Peru overview map with sampling region
# - Fig 2: Detailed sampling localities map (OpenStreetMap)

# Input: analyses.Rdata, metadata with coordinates
# Dependencies: sf, rnaturalearth, ggplot2, ggspatial, tidyterra, maptiles

################################################################################
# Setup
################################################################################

rm(list = ls())
library(tidyverse)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggplot2)
library(ggspatial)

load("analyses.Rdata")

# Get Peru boundary
world <- ne_countries(scale = "medium", returnclass = "sf")
peru <- world %>% filter(admin == "Peru")

################################################################################
# Figure 1: Overview map with sampling region
################################################################################

# Get unique sampling coordinates
coords_unique <- metadata %>%
  distinct(longitude, latitude, .keep_all = TRUE)

# Convert to sf object
pts <- st_as_sf(coords_unique, coords = c("longitude", "latitude"), crs = 4326)

# Calculate bounding box with margin
bb <- st_bbox(pts)
xpad <- (bb["xmax"] - bb["xmin"]) * 0.15
ypad <- (bb["ymax"] - bb["ymin"]) * 0.15

local_xlim <- c(bb["xmin"] - xpad, bb["xmax"] + xpad)
local_ylim <- c(bb["ymin"] - ypad, bb["ymax"] + ypad)

p1 <- ggplot() +
  geom_sf(data = peru, fill = "grey95", color = "grey70", linewidth = 0.2) +
  geom_sf(data = pts, size = 2.2, shape = 21, fill = "black", color = "white", stroke = 0.3) +
  coord_sf(xlim = local_xlim, ylim = local_ylim, expand = FALSE) +
  annotation_scale(location = "bl", width_hint = 0.35, text_cex = 0.7) +
  annotation_north_arrow(
    location = "bl", 
    which_north = "true",
    style = north_arrow_fancy_orienteering(text_size = 6)
  ) +
  theme_minimal(base_size = 9) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank()
  ) +
  labs(title = "Sampling localities (lowland rainforest, Peru)")

ggsave("Fig_map_overview.pdf", p1, width = 6, height = 7)

cat("\n✓ Map figures created!\n")
cat("Figures saved:\n")
cat("  - Fig_map_overview.pdf\n")
cat("\nNote: For detailed OpenStreetMap figure, uncomment and adapt\n")
cat("the maptiles code from the original script 04 (lines 1260-1280)\n")
