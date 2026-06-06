# Script 04: Panel Dataset Construction
# Purpose: Merge NLCD land cover panel with USDA NASS economic panel
#          Define cropland transition outcome variable
#          Validate and export final analysis-ready panel
# Author: Lan T. Tran
# Date: June 2026

library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)

# ── LOAD DATA ─────────────────────────────────────────────────────────────────

lc_panel <- read.csv("data/processed/nlcd_panel_2001_2011_2021.csv")
econ_panel <- read.csv("data/processed/econ_panel_mo.csv")

# lc_panel:   345 rows (115 counties x 3 years) — NLCD land cover
# econ_panel: 342 rows (114 counties x 3 years) — USDA NASS economic data
# St. Louis City present in lc_panel but absent in econ_panel (no farms)

# ── MERGE PANELS ──────────────────────────────────────────────────────────────
# Left join retains all 115 counties from lc_panel
# St. Louis City will have NA for all economic variables — expected

panel <- lc_panel |>
  left_join(econ_panel, by = c("GEOID", "year" = "nlcd_year"))

# Confirm: 345 rows, 3 NAs for n_farms (St. Louis City x 3 years)
nrow(panel)
sum(is.na(panel$n_farms))

# Remove duplicate year column and validation total from lc_panel
panel <- panel |>
  select(-year.y, -total)

# ── DEFINE CROPLAND TRANSITION OUTCOME ───────────────────────────────────────
# Calculate cropland change between consecutive NLCD periods
# Sorted by county and year before lagging

panel <- panel |>
  arrange(GEOID, year) |>
  group_by(GEOID) |>
  mutate(
    cropland_lag    = lag(cropland),
    cropland_change = cropland - cropland_lag
  ) |>
  ungroup()

# ── THRESHOLD SENSITIVITY NOTE ────────────────────────────────────────────────
# Transition threshold = 0.02 (2 percentage points) is the baseline definition
# Alternative thresholds tested:
#   0.01 → 100 transitions (43.5% of non-NA observations) — may capture noise
#   0.02 → 50 transitions (21.7%) — baseline choice
#   0.03 → 27 transitions (11.7%) — more conservative
#   0.05 → 2 transitions (0.9%)  — too restrictive
#
# Robustness checks in Phase 4 will re-estimate the zero-inflated spatial model
# under 0.01 and 0.03 thresholds to confirm findings are not threshold-dependent
# If results change substantially across thresholds, report all three
#
# Directional finding: all 50 transitions at 0.02 threshold are cropland GAINS
# Only 1 county shows cropland loss above 0.01 threshold
# Missouri cropland was stable or expanding during 2001-2021

panel <- panel |>
  mutate(
    transition = as.integer(abs(cropland_change) > 0.02),
    transition_dir = case_when(
      cropland_change > 0.02  ~  1L,
      cropland_change < -0.02 ~ -1L,
      is.na(cropland_change)  ~  NA_integer_,
      TRUE                    ~  0L
    )
  )

# ── VALIDATE ──────────────────────────────────────────────────────────────────

# Transition counts
table(panel$transition, useNA = "ifany")

# Summary by period
panel |>
  filter(!is.na(transition)) |>
  group_by(year) |>
  summarise(
    n_counties    = n(),
    n_transitions = sum(transition),
    pct_zero      = round(100 * mean(transition == 0), 1),
    mean_cropland = round(mean(cropland), 3),
    mean_change   = round(mean(cropland_change, na.rm = TRUE), 4)
  )

# Zero rate justification for zero-inflated model:
# 2001-2011: 81.7% zeros — exceeds 70% threshold for ZI modeling
# 2011-2021: 74.8% zeros — exceeds 70% threshold for ZI modeling

# ── EXPORT ────────────────────────────────────────────────────────────────────

write.csv(panel,
          "data/processed/panel_final.csv",
          row.names = FALSE)

message("Final panel saved: ", nrow(panel), " rows, ", ncol(panel), " columns")