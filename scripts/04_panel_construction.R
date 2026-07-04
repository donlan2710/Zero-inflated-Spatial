# Script 04: Panel Dataset Construction
# Purpose: Merge NLCD land cover panel with USDA NASS economic panel
#          Define cropland transition outcome variable
#          Validate and export final analysis-ready panel
#
# Scope: This script does one job. Land cover plus economic data in,
#        clean panel with a transition outcome out. Pasture and CRP
#        are not part of this script. See Variable_Selection_Note.md
#        for why those were suspended.
#
# Author: Lan T. Tran
# Date: June 2026

library(dplyr)
library(tidyr)
library(stringr)

# ── LOAD DATA ─────────────────────────────────────────────────────────────────

lc_panel   <- read.csv("data/processed/nlcd_panel_2001_2011_2021.csv")
econ_panel <- read.csv("data/processed/econ_panel_mo.csv")

# Row counts confirmed on the June 2026 data pull: 345 and 342.
# Read the values below yourself, do not assume they match this note
# on a future rerun with updated raw data.
nrow(lc_panel)
nrow(econ_panel)

# ── MERGE PANELS ──────────────────────────────────────────────────────────────
# Join key: GEOID and lc_panel's year matched against econ_panel's nlcd_year.
# nlcd_year is only ever 2011, 2021, or NA, set in Script 03. It is never 2001.
#
# Confirmed effect on the June 2026 pull:
#   year 2001: 115 of 115 rows missing economic data, expected, no match exists
#   year 2011: 1 of 115 rows missing, St. Louis City, no farm operations
#   year 2021: 1 of 115 rows missing, same reason
# If a future rerun shows a different pattern, check Script 03's nlcd_year
# construction and NASS withheld values first, before assuming the join broke.

panel <- lc_panel |>
  left_join(econ_panel, by = c("GEOID", "year" = "nlcd_year"))

nrow(panel)   # must equal nrow(lc_panel), the join must not add rows

panel |>
  group_by(year) |>
  summarise(n_missing_n_farms = sum(is.na(n_farms)), n_total = n())

# Check for duplicate GEOID and year combinations, a join should never
# produce these here, but this check costs nothing and catches a bad
# upstream file before it silently corrupts everything downstream.
panel |>
  count(GEOID, year) |>
  filter(n > 1)

panel <- panel |>
  select(-year.y, -total)

# ── DEFINE CROPLAND TRANSITION OUTCOME ───────────────────────────────────────
# cropland_lag = cropland share in the previous NLCD period.
#   year 2011 row: cropland_lag equals the 2001 cropland level
#   year 2021 row: cropland_lag equals the 2011 cropland level
#   year 2001 row: cropland_lag is NA, there is no earlier period
# This is the beginning-of-period path dependence control for Part 2
# in Script 09. Confirmed on the June 2026 pull, no missing values in
# cropland_lag for the year 2021 rows, the set Script 09 actually uses.

panel <- panel |>
  arrange(GEOID, year) |>
  group_by(GEOID) |>
  mutate(
    cropland_lag    = lag(cropland),
    cropland_change = cropland - cropland_lag
  ) |>
  ungroup()

# ── DEFINE TRANSITION INDICATOR ───────────────────────────────────────────────
# Baseline threshold is 0.02. Confirmed threshold table from the June 2026
# pull below, read as a record of what was tested, not as a fixed truth
# to assume on future data.
#
#   0.01 -> 100 of 230 valid rows, 43.5 percent
#   0.02 ->  50 of 230 valid rows, 21.7 percent
#   0.03 ->  27 of 230 valid rows, 11.7 percent
#   0.05 ->   2 of 230 valid rows,  0.9 percent
#
# Rerun the loop below every time, do not copy these numbers into a
# later script without checking them again on the current panel.

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

for (t in c(0.01, 0.02, 0.03, 0.05)) {
  n_valid <- sum(!is.na(panel$cropland_change))
  n_trans <- sum(abs(panel$cropland_change) > t, na.rm = TRUE)
  message("Threshold ", t, ": ", n_trans, " transitions out of ",
          n_valid, " valid rows (",
          round(100 * n_trans / n_valid, 1), "%)")
}

# ── VALIDATE ──────────────────────────────────────────────────────────────────

table(panel$transition, useNA = "ifany")

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

# The 2021 row set here still includes St. Louis City, 115 counties.
# St. Louis City is dropped later, in the spatial regression scripts,
# because it has no farm data and no meaningful land use decision to
# model. Do not expect n_counties to equal 114 here, that drop has not
# happened yet at this stage of the pipeline.

panel |>
  filter(year == 2021) |>
  select(GEOID, cropland, cropland_lag, cropland_change) |>
  summary()

# ── EXPORT ────────────────────────────────────────────────────────────────────

write.csv(panel,
          "data/processed/panel_final.csv",
          row.names = FALSE)

message("Final panel saved: ", nrow(panel), " rows, ", ncol(panel), " columns")