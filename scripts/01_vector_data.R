# Script 01: Vector Data Operations
# Purpose: Load, inspect, and process US county shapefiles
# Data: Census TIGER/Line Missouri county shapefile
# Author: Lan T. Tran
# Date: June 2026
install.packages("sf")
library(sf)
# Load US county shapefile
counties_us <- st_read("data/raw/tl_2025_us_county/tl_2025_us_county.shp")

# Inspect the data
print(counties_us)
nrow(counties_us)      # number of counties
names(counties_us)     # column names
st_crs(counties_us)    # coordinate reference system

# Filter to Missouri only (FIPS state code = 29)
counties_mo <- counties_us[counties_us$STATEFP == "29", ]

# Confirm
nrow(counties_mo)      
print(counties_mo)

# Result: 115 Missouri counties loaded successfully
# CRS: NAD83 (standard for US data, no reprojection needed yet)
# Key join variable: GEOID (5-digit FIPS: state + county)

# ── INSPECT ──────────────────────────────────────────────────────────────────

# Geometry type
st_geometry_type(counties_mo) |> unique()

# CRS detail
st_crs(counties_mo)

# Attribute fields
names(counties_mo)

# ── SELECT relevant columns only ─────────────────────────────────────────────

counties_mo_clean <- counties_mo[, c("GEOID", "NAME", "ALAND", "AWATER", "geometry")]

# ── REPROJECT to a projected CRS suitable for Missouri ───────────────────────
# NAD83 is geographic (degrees), not projected (meters)
# EPSG 26915 = UTM Zone 15N, standard for Missouri

counties_mo_proj <- st_transform(counties_mo_clean, crs = 26915)

# Confirm new CRS
st_crs(counties_mo_proj)$epsg

# ── CALCULATE AREA ────────────────────────────────────────────────────────────
# Now in projected CRS so area is in square meters

counties_mo_proj$area_sqkm <- as.numeric(st_area(counties_mo_proj)) / 1e6

# Check results
summary(counties_mo_proj$area_sqkm)

# ── EXPORT to data/processed ──────────────────────────────────────────────────

st_write(counties_mo_proj, 
         "data/processed/counties_mo.shp",
         delete_if_exists = TRUE)

# Verify the exported file reads back correctly
counties_check <- st_read("data/processed/counties_mo.shp")
print(counties_check)
summary(counties_check$area_sqkm)