library(terra)
library(sf)

# Load processed Missouri counties
counties_mo_proj <- st_read("data/processed/counties_mo.shp")

# Load full NLCD 2021 raster
nlcd_conus <- rast("data/raw/Annual_NLCD_LndCov_2021_CU_C1V1/Annual_NLCD_LndCov_2021_CU_C1V1.tif")

# Check CRS of NLCD
crs(nlcd_conus, describe = TRUE)

# Reproject Missouri counties to match NLCD CRS
counties_mo_nlcd <- st_transform(counties_mo_proj, crs(nlcd_conus))

# Crop and mask to Missouri boundary
nlcd_mo <- crop(nlcd_conus, vect(counties_mo_nlcd))
nlcd_mo <- mask(nlcd_mo, vect(counties_mo_nlcd))

# Check result
nlcd_mo

writeRaster(nlcd_mo, 
            "data/raw/nlcd_2021_missouri.tif", 
            overwrite = TRUE)

# Note: full CONUS raster (1.3 GB) stored outside repository
# Only nlcd_2021_missouri.tif (66 MB) needed for analysis
# Future years: repeat crop and save workflow for each year

# Reload from saved file to confirm it works correctly
nlcd_mo <- rast("data/raw/nlcd_2021_missouri.tif")
nlcd_mo

# Quick plot to visually inspect land cover classes
plot(nlcd_mo, main = "NLCD 2021 Land Cover - Missouri")

png("outputs/map_02_nlcd_2021_missouri.png", 
    width = 10, height = 8, units = "in", res = 300)
plot(nlcd_mo, main = "NLCD 2021 Land Cover - Missouri")
dev.off()

library(exactextractr)
library(dplyr)
# Extract land cover proportions by county
# exact_extract computes the fraction of each pixel covered by each polygon
# This gives us the proportion of each land cover class per county
lc_extract <- exact_extract(nlcd_mo, counties_mo_nlcd,
                            "frac",
                            append_cols = "GEOID")

# Check result
head(lc_extract)
ncol(lc_extract)
nrow(lc_extract)

# Rename for clarity
lc_extract <- lc_extract |>
  rename(cropland = frac_82,
         pasture = frac_81,
         forest = frac_41,
         developed_low = frac_21,
         developed_high = frac_24,
         water = frac_11,
         wetland = frac_90)

# Validate: all rows should sum to approximately 1
lc_extract$total <- rowSums(lc_extract[, -1])
summary(lc_extract$total)

# Export to data/processed
write.csv(lc_extract, 
          "data/processed/nlcd_2021_county_extract.csv",
          row.names = FALSE)