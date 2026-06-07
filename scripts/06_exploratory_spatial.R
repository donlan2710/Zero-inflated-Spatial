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

panel        <- read.csv("data/processed/panel_final.csv")
w_queen      <- readRDS("data/processed/weights_queen.rds")
nb_queen     <- readRDS("data/processed/nb_queen.rds")

# ── PREPARE CROSS-SECTIONS ────────────────────────────────────────────────────
# Moran's I requires a single cross-section not a panel
# 2021 is the baseline year for exploratory analysis

panel_2021 <- panel |>
  filter(year == 2021) |>
  arrange(GEOID)

nrow(panel_2021)   # confirm 115 rows match shapefile

# ── GLOBAL MORAN'S I ──────────────────────────────────────────────────────────
# Tests whether values cluster spatially
# H0: no spatial autocorrelation
# Positive I: similar values cluster together

# Cropland level
moran_cropland <- moran.test(panel_2021$cropland, w_queen)
print(moran_cropland)
# Result: Moran I = 0.783, p < 2.2e-16 — strong spatial clustering in cropland levels

# Cropland change
moran_change <- moran.test(panel_2021$cropland_change,
                           w_queen,
                           na.action = na.omit)
print(moran_change)

# Transition indicator
moran_transition <- moran.test(panel_2021$transition,
                               w_queen,
                               na.action = na.omit)
print(moran_transition)
# Result: Moran I = 0.509, p < 2.2e-16 — significant spatial clustering in transitions
# Both results justify spatial regression over standard OLS

# ── LOCAL MORAN'S I (LISA) — CROPLAND LEVEL ───────────────────────────────────
# Identifies specific cluster locations
# High-High: high cropland surrounded by high cropland neighbors
# Low-Low:   low cropland surrounded by low cropland neighbors
# Spatial outliers: High-Low, Low-High

lisa_cropland <- localmoran(panel_2021$cropland, w_queen)

counties_lisa <- counties_mo_proj
counties_lisa$lisa_i   <- lisa_cropland[, 1]
counties_lisa$lisa_p   <- lisa_cropland[, 5]
counties_lisa$cropland <- panel_2021$cropland

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

# Cluster counts
table(counties_lisa$cluster)
# High-High: 14, Low-Low: 23, Not significant: 78
# Low-Low cluster = Ozark highlands = primary structural zero region

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

# ── LOCAL MORAN'S I (LISA) — TRANSITION INDICATOR ────────────────────────────

lisa_trans <- localmoran(panel_2021$transition, w_queen)

counties_lisa$lisa_trans_i <- lisa_trans[, 1]
counties_lisa$lisa_trans_p <- lisa_trans[, 5]
counties_lisa$transition   <- panel_2021$transition

mean_trans <- mean(panel_2021$transition, na.rm = TRUE)

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

# Cluster counts
table(counties_lisa$cluster_trans)
# High-High: 15 (northern Missouri grain belt — transitions clustered)
# Low-High: 4 (spatial outliers — non-transitioning within transitioning cluster)
# Not significant: 96 (structural zero region — south and central Missouri)

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

message("Script 06 complete")