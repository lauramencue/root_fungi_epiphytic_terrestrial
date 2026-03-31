# Project name: root_fungi_epiphytic_terrestrial
# Author: Laura Méndez
# Date: [Date]

# Goals:
# Create publication-quality figures for bipartite network analyses

# Figures created:
# - Network visualization with igraph/ggraph
# - Network metrics comparison plots

# Input: analyses.Rdata, bipartite network results from script 10
# Dependencies: ggplot2, tidyverse, igraph, ggraph, bipartite

# Note: Adapt plotting code from script 10 (bipartite networks)
# based on which visualizations are most informative for your manuscript

################################################################################
# Setup
################################################################################

rm(list = ls())
library(tidyverse)
library(ggplot2)
library(bipartite)
library(igraph)

# Load network analysis results
load("analyses.Rdata")

cat("\n✓ Network figures script created!\n")
cat("Add specific network visualization code from script 10\n")
