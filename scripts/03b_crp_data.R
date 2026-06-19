# Script 03b: CRP Enrollment Data
# Purpose: Load USDA FSA county-level CRP enrollment history
#          Extract 2011 enrollment for Missouri counties
#          Convert to share of county land area for comparability
#          with pasture_2011 (NLCD-derived) in the active transition model
# Data: USDA FSA CRP Enrollment and Rental Payments by County, 1986-2025
# Source: https://www.fsa.usda.gov/tools/informational/reports/conservation-statistics/crp
# File: CRPHistoryCounty86-25.xlsx
# Author: Lan T. Tran
# Date: June 2026

library(sf)
library(readxl)
library(dplyr)

# ── LOAD MISSOURI COUNTY SHAPEFILE FOR LAND AREA ─────────────────────────────
# Land area (ALAND, square meters) is needed to convert CRP acres to a
# proportion of county area, matching the scale of pasture_2011

counties_mo_proj <- st_read("data/processed/counties_mo.shp")

area_lookup <- counties_mo_proj |>
  st_drop_geometry() |>
  select(GEOID, ALAND) |>
  mutate(
    GEOID      = as.integer(GEOID),
    land_acres = ALAND * 0.000247105   # sq meters to acres
  ) |>
  select(GEOID, land_acres)

# ── LOAD CRP COUNTY HISTORY ───────────────────────────────────────────────────
# File header occupies rows 1-3 (title, blank, column names)
# Row 4 onward is data; skip = 3 aligns row 4 as the first data row

crp_county <- read_excel("data/raw/CRPHistoryCounty86-25.xlsx", skip = 3)

# ── FILTER TO MISSOURI, EXTRACT 2011 ENROLLMENT ──────────────────────────────

crp_2011_mo <- crp_county |>
  filter(STATE == "MISSOURI") |>
  select(GEOID = FIPS, crp_acres_2011 = `2011`)

# Confirm full Missouri coverage
nrow(crp_2011_mo)                          # should be 115
sum(is.na(crp_2011_mo$crp_acres_2011))     # should be 0

# ── CONVERT TO SHARE OF COUNTY LAND AREA ─────────────────────────────────────

crp_2011_mo <- crp_2011_mo |>
  left_join(area_lookup, by = "GEOID") |>
  mutate(crp_share_2011 = crp_acres_2011 / land_acres) |>
  select(GEOID, crp_acres_2011, crp_share_2011)

# Validate
summary(crp_2011_mo$crp_share_2011)
sum(is.na(crp_2011_mo$crp_share_2011))     # should be 0

# ── EXPORT ────────────────────────────────────────────────────────────────────

write.csv(crp_2011_mo,
          "data/processed/crp_2011_missouri.csv",
          row.names = FALSE)

message("CRP data saved: ", nrow(crp_2011_mo), " counties")