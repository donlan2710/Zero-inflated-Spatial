# Script 02: Raster Data Extraction
# Purpose: Load NLCD land cover rasters, crop to Missouri,
#          extract land cover proportions by county for all years
# Data: NLCD 2001, 2011, 2021 Land Cover
# Author: Lan T. Tran
# Date: June 2026

library(terra)
library(sf)
library(exactextractr)
library(dplyr)
library(ggplot2)

# ── LOAD MISSOURI COUNTIES ────────────────────────────────────────────────────

counties_mo_proj <- st_read("data/processed/counties_mo.shp")

# ── CROP AND SAVE EACH YEAR ───────────────────────────────────────────────────

years <- c(2001, 2011, 2021)

for (yr in years) {
  
  infile <- paste0("data/raw/Annual_NLCD_LndCov_", yr, 
                   "_CU_C1V1/Annual_NLCD_LndCov_", yr, "_CU_C1V1.tif")
  outfile <- paste0("data/raw/nlcd_", yr, "_missouri.tif")
  
  message("Processing year: ", yr)
  
  nlcd_raw <- rast(infile)
  
  # Reproject counties to match NLCD CRS
  counties_mo_nlcd <- st_transform(counties_mo_proj, crs(nlcd_raw))
  
  nlcd_cropped <- crop(nlcd_raw, vect(counties_mo_nlcd))
  nlcd_masked <- mask(nlcd_cropped, vect(counties_mo_nlcd))
  
  writeRaster(nlcd_masked, outfile, overwrite = TRUE)
  
  message("Saved: ", outfile)
}

# ── EXTRACT LAND COVER PROPORTIONS BY COUNTY FOR ALL YEARS ───────────────────

# Store results for each year in a list
extract_list <- list()

for (yr in years) {
  
  message("Extracting year: ", yr)
  
  # Load cropped Missouri raster
  nlcd_mo <- rast(paste0("data/raw/nlcd_", yr, "_missouri.tif"))
  
  # Reproject counties to match raster CRS
  counties_mo_nlcd <- st_transform(counties_mo_proj, crs(nlcd_mo))
  
  # Extract proportions
  lc_extract <- exact_extract(nlcd_mo, counties_mo_nlcd,
                              "frac",
                              append_cols = "GEOID")
  
  # Add year column
  lc_extract$year <- yr
  
  # Store in list
  extract_list[[as.character(yr)]] <- lc_extract
}

# Combine all years into one dataframe
lc_panel <- bind_rows(extract_list)

# Check
nrow(lc_panel)   # should be 115 x 3 = 345
head(lc_panel)

# ── RENAME KEY COLUMNS ────────────────────────────────────────────────────────

lc_panel <- lc_panel |>
  rename(cropland   = frac_82,
         pasture    = frac_81,
         forest     = frac_41,
         developed  = frac_21,
         water      = frac_11,
         wetland    = frac_90)

# ── VALIDATE ──────────────────────────────────────────────────────────────────

# All rows should sum to 1
frac_cols <- names(lc_panel)[grepl("frac_", names(lc_panel))]
named_cols <- c("cropland", "pasture", "forest", "developed", "water", "wetland")
all_cols <- c(frac_cols, named_cols)

lc_panel$total <- rowSums(lc_panel[, intersect(all_cols, names(lc_panel))], 
                          na.rm = TRUE)
summary(lc_panel$total)

# ── EXPORT ────────────────────────────────────────────────────────────────────

write.csv(lc_panel,
          "data/processed/nlcd_panel_2001_2011_2021.csv",
          row.names = FALSE)

message("Panel dataset saved: 115 counties x 3 years")

# ── MAP: CROPLAND 2021 ────────────────────────────────────────────────────────

nlcd_mo_2021 <- rast("data/raw/nlcd_2021_missouri.tif")
counties_mo_nlcd <- st_transform(counties_mo_proj, crs(nlcd_mo_2021))

# Plot raw land cover
png("outputs/map_02_nlcd_2021_missouri.png",
    width = 14, height = 8, units = "in", res = 300)
plot(nlcd_mo_2021, 
     main = "NLCD 2021 Land Cover - Missouri",
     mar = c(3, 3, 3, 12))
dev.off()

# Cropland proportion map
lc_2021 <- lc_panel[lc_panel$year == 2021, c("GEOID", "cropland")]
counties_map <- merge(counties_mo_nlcd, lc_2021, by = "GEOID")

ggplot(counties_map) +
  geom_sf(aes(fill = cropland), color = "white", linewidth = 0.3) +
  scale_fill_viridis_c(name = "Cropland\nProportion",
                       option = "viridis",
                       labels = scales::percent) +
  labs(
    title = "Cultivated Cropland Share by County",
    subtitle = "Missouri, 2021 — NLCD Land Cover Class 82",
    caption = "Source: NLCD 2021"
  ) +
  theme_minimal()

ggsave("outputs/map_03_cropland_2021.png",
       width = 10, height = 7, dpi = 300)

# ── MAP: CROPLAND CHANGE ACROSS YEARS ────────────────────────────────────────

library(patchwork)  # for combining plots side by side

map_year <- function(yr) {
  lc_yr <- lc_panel[lc_panel$year == yr, c("GEOID", "cropland")]
  counties_yr <- merge(counties_mo_nlcd, lc_yr, by = "GEOID")
  
  ggplot(counties_yr) +
    geom_sf(aes(fill = cropland), color = "white", linewidth = 0.3) +
    scale_fill_viridis_c(name = "Cropland",
                         option = "viridis",
                         labels = scales::percent,
                         limits = c(0, 1)) +
    labs(title = paste("Cropland Share", yr)) +
    theme_minimal() +
    theme(legend.position = "none")
}

p2001 <- map_year(2001)
p2011 <- map_year(2011)
p2021 <- map_year(2021)

# Combine
library(patchwork)

combined <- p2001 + p2011 + p2021 +
  plot_annotation(
    title = "Cultivated Cropland Share by County — Missouri",
    subtitle = "NLCD Land Cover Class 82",
    caption = "Source: NLCD 2001, 2011, 2021"
  )

ggsave("outputs/map_04_cropland_change_2001_2011_2021.png",
       combined,
       width = 18, height = 7, dpi = 300)