# Script 08: Zero Structure Diagnosis
# Purpose: Formally document where and why standard spatial models fail
#          on sparse land use transition data
#          Characterize the structural zero distribution
#          Motivate the zero-inflated spatial model
# Author: Lan T. Tran
# Date: June 2026

library(sf)
library(spdep)
library(dplyr)
library(ggplot2)

# ── LOAD DATA ─────────────────────────────────────────────────────────────────

counties_mo_proj <- st_read("data/processed/counties_mo.shp") |>
  arrange(GEOID)

panel   <- read.csv("data/processed/panel_final.csv")
w_queen <- readRDS("data/processed/weights_queen.rds")
slag    <- readRDS("data/processed/model_slag.rds")

# ── PREPARE ANALYTICAL SAMPLE ─────────────────────────────────────────────────

df <- panel |>
  filter(year == 2021) |>
  filter(!is.na(n_farms)) |>
  arrange(GEOID) |>
  mutate(
    land_value_z = scale(land_value_acre),
    net_income_z = scale(net_income),
    n_farms_z    = scale(n_farms),
    cropland_z   = scale(cropland),
    acres_z      = scale(acres_operated)
  )

geoids_114 <- df$GEOID
w_114 <- subset(w_queen,
                subset = counties_mo_proj$GEOID %in% geoids_114)

# ── 1. CHARACTERIZE ZERO STRUCTURE ───────────────────────────────────────────
# Zero rate, spatial clustering of zeros

cat("=== ZERO STRUCTURE ANALYSIS ===\n\n")

cat("Transition distribution (2011-2021):\n")
table(df$transition)

cat("\nZero rate:", round(100 * mean(df$transition == 0), 1), "%\n")
cat("Transition rate:", round(100 * mean(df$transition == 1), 1), "%\n")

# Moran's I on zero indicator
# Tests whether zeros cluster spatially rather than distributing randomly
zero_indicator <- as.integer(df$transition == 0)
moran_zeros <- moran.test(zero_indicator, w_114)
cat("\nMoran's I on zero indicator:\n")
print(moran_zeros)
# Result: Moran I = 0.507, p < 2.2e-16
# Zeros are spatially clustered — concentrated in Ozark highlands

# ── 2. SEPARATE STRUCTURAL FROM SAMPLING ZEROS ───────────────────────────────
# Structural zeros: counties in Low-Low LISA cluster for cropland level
# These counties have low cropland AND low-cropland neighbors
# Structurally constrained by terrain — cannot expand cropland regardless
# of economic conditions

lisa_cropland <- localmoran(df$cropland, w_114)
mean_crop     <- mean(df$cropland)
lag_crop      <- lag.listw(w_114, df$cropland)

df <- df |>
  mutate(
    lisa_p          = lisa_cropland[, 5],
    lag_cropland    = lag_crop,
    structural_zero = as.integer(
      cropland < mean_crop &
        lag_cropland < mean_crop &
        lisa_p < 0.05
    )
  )

cat("\nThree-way decomposition:\n")
cat("Structural zeros (Low-Low LISA):", sum(df$structural_zero), "\n")
cat("Sampling zeros (non-transition, not structural):",
    sum(df$transition == 0) - sum(df$structural_zero), "\n")
cat("Transition counties:", sum(df$transition == 1), "\n")

# Cross-tabulation confirms perfect separation
cat("\nCross-tabulation (structural zero x transition):\n")
table(Structural = df$structural_zero,
      Transition = df$transition)
# Key result: structural=1/transition=1 cell = 0
# No structural zero county ever transitioned — confirms two-process structure

# ── 3. STANDARD MODEL FAILURE DIAGNOSTICS ────────────────────────────────────
# Standard spatial lag residuals systematically biased by county type
# Under-predicts transitions, over-predicts for zeros

cat("\n=== STANDARD MODEL FAILURE DIAGNOSTICS ===\n")

df$slag_resid <- residuals(slag)

cat("\nResidual distribution summary:\n")
print(summary(df$slag_resid))

cat("\nMean residual by county type:\n")
df |>
  mutate(county_type = case_when(
    structural_zero == 1 ~ "Structural zero",
    transition == 0      ~ "Sampling zero",
    transition == 1      ~ "Transition"
  )) |>
  group_by(county_type) |>
  summarise(
    n          = n(),
    mean_resid = round(mean(slag_resid), 4),
    sd_resid   = round(sd(slag_resid), 4)
  ) |>
  print()
# Structural zero: mean = -0.114, SD = 0.080 — homogeneous, over-predicted
# Sampling zero:   mean = -0.167, SD = 0.191 — heterogeneous, over-predicted
# Transition:      mean = +0.447, SD = 0.196 — systematically under-predicted
# Classic signature of zero-inflated misspecification

# ── 4. VISUALIZE THREE-WAY DECOMPOSITION ─────────────────────────────────────

df <- df |>
  mutate(county_type = case_when(
    structural_zero == 1 ~ "Structural zero",
    transition == 0      ~ "Sampling zero",
    transition == 1      ~ "Transition"
  ))

counties_diag <- counties_mo_proj |>
  filter(GEOID %in% df$GEOID) |>
  left_join(
    df[, c("GEOID", "county_type")] |>
      mutate(GEOID = as.character(GEOID)),
    by = "GEOID"
  )

ggplot(counties_diag) +
  geom_sf(aes(fill = county_type), color = "white", linewidth = 0.3) +
  scale_fill_manual(
    values = c(
      "Structural zero" = "#2c7bb6",
      "Sampling zero"   = "#f0f0f0",
      "Transition"      = "#d7191c"
    ),
    name = "County Type"
  ) +
  labs(
    title = "Three-Way Decomposition of Missouri Counties",
    subtitle = "Structural zeros, sampling zeros, and transitions — 2011–2021",
    caption = "Structural zeros defined as Low-Low LISA cluster (p < 0.05)"
  ) +
  theme_minimal()

ggsave("outputs/map_08_zero_structure.png",
       width = 10, height = 7, dpi = 300)

# ── 5. WRITE DIAGNOSTIC MEMO ──────────────────────────────────────────────────

sink("docs/diagnostic_memo.txt")
cat("=== ZERO STRUCTURE DIAGNOSTIC MEMO ===\n")
cat("Date: June 2026\n\n")

cat("1. ZERO PREVALENCE\n")
cat("Zero rate (2011-2021): 74.6%\n")
cat("Structural zeros (Low-Low LISA): 23 counties\n")
cat("Sampling zeros: 62 counties\n")
cat("Transitions: 29 counties\n\n")

cat("2. SPATIAL CLUSTERING OF ZEROS\n")
cat("Moran's I on zero indicator: 0.507, p < 2.2e-16\n")
cat("Zeros are not randomly distributed — concentrated in Ozark highlands\n\n")

cat("3. STRUCTURAL VS SAMPLING ZEROS\n")
cat("Perfect separation: zero structural zero counties transitioned\n")
cat("Cross-tabulation: structural=1/transition=1 cell = 0\n")
cat("Confirms two distinct processes in the data\n\n")

cat("4. STANDARD MODEL FAILURE\n")
cat("Spatial lag residuals by county type:\n")
cat("Structural zero: mean = -0.114, SD = 0.080\n")
cat("Sampling zero:   mean = -0.167, SD = 0.191\n")
cat("Transition:      mean = +0.447, SD = 0.196\n")
cat("Systematic bias confirms standard model misspecification\n\n")

cat("5. CONCLUSION\n")
cat("Standard spatial lag model addresses spatial dependence\n")
cat("but cannot separate structural from sampling zeros.\n")
cat("Zero-inflated spatial model required to:\n")
cat("  - Model structural non-transition process separately\n")
cat("  - Recover unbiased estimates of transition drivers\n")
cat("  - Correctly attribute spatial dependence\n")
sink()

message("Script 08 complete - zero structure diagnosis documented")