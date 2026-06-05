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