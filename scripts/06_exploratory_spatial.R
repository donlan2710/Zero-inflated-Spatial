# Script 06: Exploratory Spatial Analysis
# Purpose: Test for spatial autocorrelation in cropland and transition variables
#          Global Moran's I and Local Moran's I (LISA) cluster analysis
#          Results motivate spatial regression models in Script 07
# Author: Lan T. Tran
# Date: June 2026

library(sf)
library(spdep)
library(dplyr)
library(ggplot2)

# ── LOAD DATA ─────────────────────────────────────────────────────────────────

counties_mo_proj <- st_read("data/processed/counties_mo.shp") |>
  arrange(GEOID)

panel <- read.csv("data/processed/panel_final.csv")
w_queen <- readRDS("data/processed/weights_queen.rds")
nb_queen <- readRDS("data/processed/nb_queen.rds")

# ── PREPARE CROSS-SECTION FOR SPATIAL ANALYSIS ───────────────────────────────
# Moran's I requires a single cross-section not a panel
# Use 2021 as the baseline year for exploratory analysis
# Repeat for 2011 as robustness check

panel_2021 <- panel |>
  filter(year == 2021) |>
  arrange(GEOID)

# Confirm row count matches shapefile
nrow(panel_2021)   # should be 115

# ── GLOBAL MORAN'S I - CROPLAND LEVEL ────────────────────────────────────────
# Tests whether cropland proportion clusters spatially
# H0: no spatial autocorrelation
# Positive I: similar values cluster together
# Negative I: dissimilar values cluster together

moran_cropland <- moran.test(panel_2021$cropland, w_queen)
print(moran_cropland)

# ── GLOBAL MORAN'S I - CROPLAND CHANGE ───────────────────────────────────────
# Tests whether cropland transitions cluster spatially
# Use 2021 period (2011-2021 transitions)

panel_2021_trans <- panel |>
  filter(year == 2021) |>
  arrange(GEOID)

moran_change <- moran.test(panel_2021_trans$cropland_change, 
                           w_queen,
                           na.action = na.omit)
print(moran_change)

# ── GLOBAL MORAN'S I - TRANSITION INDICATOR ──────────────────────────────────

moran_transition <- moran.test(panel_2021_trans$transition,
                               w_queen,
                               na.action = na.omit)
print(moran_transition)

# ── LOCAL MORAN'S I (LISA) ────────────────────────────────────────────────────
# Identifies specific cluster locations
# High-High: high cropland counties surrounded by high cropland neighbors
# Low-Low: low cropland counties surrounded by low cropland neighbors
# High-Low: spatial outliers

lisa <- localmoran(panel_2021$cropland, w_queen)

# Add LISA results to county shapefile
counties_lisa <- counties_mo_proj
counties_lisa$lisa_i    <- lisa[, 1]   # local Moran statistic
counties_lisa$lisa_p    <- lisa[, 5]   # p-value
counties_lisa$cropland  <- panel_2021$cropland

# Classify clusters
mean_crop <- mean(panel_2021$cropland)

counties_lisa <- counties_lisa |>
  mutate(
    lag_cropland = lag.listw(w_queen, cropland),
    cluster = case_when(
      cropland > mean_crop & lag_cropland > mean_crop & lisa_p < 0.05 ~ "High-High",
      cropland < mean_crop & lag_cropland < mean_crop & lisa_p < 0.05 ~ "Low-Low",
      cropland > mean_crop & lag_cropland < mean_crop & lisa_p < 0.05 ~ "High-Low",
      cropland < mean_crop & lag_cropland > mean_crop & lisa_p < 0.05 ~ "Low-High",
      TRUE ~ "Not significant"
    )
  )

# Count cluster types
table(counties_lisa$cluster)

# ── MAP LISA CLUSTERS ─────────────────────────────────────────────────────────

ggplot(counties_lisa) +
  geom_sf(aes(fill = cluster), color = "white", linewidth = 0.3) +
  scale_fill_manual(
    values = c(
      "High-High"       = "#d7191c",
      "Low-Low"         = "#2c7bb6",
      "High-Low"        = "#fdae61",
      "Low-High"        = "#abd9e9",
      "Not significant" = "#f0f0f0"
    ),
    name = "Cluster Type"
  ) +
  labs(
    title = "LISA Clusters — Cropland Share",
    subtitle = "Missouri Counties, 2021",
    caption = "Queen contiguity weights, p < 0.05"
  ) +
  theme_minimal()

ggsave("outputs/map_06_lisa_clusters.png",
       width = 10, height = 7, dpi = 300)

# ── LISA FOR TRANSITION INDICATOR ────────────────────────────────────────────

lisa_trans <- localmoran(panel_2021_trans$transition, w_queen)

counties_lisa$lisa_trans_i <- lisa_trans[, 1]
counties_lisa$lisa_trans_p <- lisa_trans[, 5]
counties_lisa$transition   <- panel_2021_trans$transition

mean_trans <- mean(panel_2021_trans$transition, na.rm = TRUE)

counties_lisa <- counties_lisa |>
  mutate(
    lag_transition = lag.listw(w_queen, transition),
    cluster_trans = case_when(
      transition > mean_trans & lag_transition > mean_trans & lisa_trans_p < 0.05 ~ "High-High",
      transition < mean_trans & lag_transition < mean_trans & lisa_trans_p < 0.05 ~ "Low-Low",
      transition > mean_trans & lag_transition < mean_trans & lisa_trans_p < 0.05 ~ "High-Low",
      transition < mean_trans & lag_transition > mean_trans & lisa_trans_p < 0.05 ~ "Low-High",
      TRUE ~ "Not significant"
    )
  )

table(counties_lisa$cluster_trans)

ggplot(counties_lisa) +
  geom_sf(aes(fill = cluster_trans), color = "white", linewidth = 0.3) +
  scale_fill_manual(
    values = c(
      "High-High"       = "#d7191c",
      "Low-Low"         = "#2c7bb6",
      "High-Low"        = "#fdae61",
      "Low-High"        = "#abd9e9",
      "Not significant" = "#f0f0f0"
    ),
    name = "Cluster Type"
  ) +
  labs(
    title = "LISA Clusters — Cropland Transition",
    subtitle = "Missouri Counties, 2011–2021",
    caption = "Queen contiguity weights, p < 0.05"
  ) +
  theme_minimal()

ggsave("outputs/map_07_lisa_transition_clusters.png",
       width = 10, height = 7, dpi = 300)