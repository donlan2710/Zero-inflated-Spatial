# Script 05: Spatial Weights Matrix Construction
# Purpose: Build spatial weights matrices from Missouri county shapefile
#          Test sensitivity to different weights specifications
#          Prepare spatial structure for regression analysis
# Author: Lan T. Tran
# Date: June 2026

library(sf)
library(spdep)
library(dplyr)

# ── LOAD DATA ─────────────────────────────────────────────────────────────────

counties_mo_proj <- st_read("data/processed/counties_mo.shp")
panel <- read.csv("data/processed/panel_final.csv")

# Confirm county count matches
nrow(counties_mo_proj)        # 115
length(unique(panel$GEOID))   # 115

# ── SORT SHAPEFILE BY GEOID ───────────────────────────────────────────────────
# Critical: weights matrix rows must align with data rows in regression
# Sort both by GEOID to ensure consistent ordering

counties_mo_proj <- counties_mo_proj |>
  arrange(GEOID)

# ── QUEEN CONTIGUITY WEIGHTS ──────────────────────────────────────────────────
# Queen: counties sharing any boundary point are neighbors
# Standard specification for county-level spatial analysis

nb_queen <- poly2nb(counties_mo_proj, queen = TRUE)
summary(nb_queen)

# Check for islands - counties with zero neighbors
islands <- which(card(nb_queen) == 0)
if (length(islands) == 0) {
  message("No islands - all counties have at least one neighbor")
} else {
  message("Islands found: ", counties_mo_proj$NAME[islands])
}

# Distribution of neighbor counts
table(card(nb_queen))

# Convert to row-standardized weights list
# Row standardization: spatial lag = weighted average of neighbor values
w_queen <- nb2listw(nb_queen, style = "W")

# ── ROOK CONTIGUITY COMPARISON ────────────────────────────────────────────────
# Rook: shared edge only (no corner touches)
# Result: identical to queen for Missouri counties
# No county pairs share corner-only boundaries in Missouri
# Robustness check will use distance-based weights instead

nb_rook <- poly2nb(counties_mo_proj, queen = FALSE)

cat("Queen neighbors:", sum(card(nb_queen)), "\n")
cat("Rook neighbors: ", sum(card(nb_rook)), "\n")
cat("Difference:     ", sum(card(nb_queen)) - sum(card(nb_rook)), "\n")

# Queen = Rook for Missouri: no corner-only neighbor relationships exist

# ── DISTANCE-BASED WEIGHTS ────────────────────────────────────────────────────
# Alternative specification: neighbors defined by centroid distance
# Used as robustness check since queen = rook for Missouri
# Threshold: 60km — minimum distance ensuring full connectivity
# (45km threshold leaves 2 disjoint subgraphs)

coords <- st_coordinates(st_centroid(counties_mo_proj))

nb_dist <- dnearneigh(coords, d1 = 0, d2 = 60000)
summary(nb_dist)

# Confirm fully connected
cat("Number of subgraphs:", n.comp.nb(nb_dist)$nc, "\n")

# Convert to row-standardized weights list
w_dist <- nb2listw(nb_dist, style = "W")

# ── WEIGHTS SPECIFICATION COMPARISON ─────────────────────────────────────────
# Queen contiguity: 588 links, avg 5.1 neighbors, 4.4% nonzero
# Distance 60km:   616 links, avg 5.4 neighbors, 4.7% nonzero
# Similar structure — distance adds ~28 additional neighbor pairs
# Both specifications will be used in regression robustness checks

# ── SAVE WEIGHTS OBJECTS ──────────────────────────────────────────────────────
# spdep weights cannot be saved as CSV — use saveRDS

saveRDS(w_queen, "data/processed/weights_queen.rds")
saveRDS(w_dist,  "data/processed/weights_dist60km.rds")
saveRDS(nb_queen, "data/processed/nb_queen.rds")

message("Spatial weights saved successfully")

# ── VISUALIZE WEIGHTS NETWORKS ────────────────────────────────────────────────

png("outputs/map_05_spatial_weights.png",
    width = 12, height = 6, units = "in", res = 300)
par(mfrow = c(1, 2))

plot(st_geometry(counties_mo_proj), border = "grey",
     main = "Queen Contiguity")
plot(nb_queen, coords, add = TRUE, col = "blue", lwd = 0.5)

plot(st_geometry(counties_mo_proj), border = "grey",
     main = "Distance 60km")
plot(nb_dist, coords, add = TRUE, col = "red", lwd = 0.5)

dev.off()

message("Script 05 complete")