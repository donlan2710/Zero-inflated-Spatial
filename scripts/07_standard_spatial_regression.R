# Script 07: Standard Spatial Econometric Regression
# Purpose: Estimate OLS and spatial regression models as baseline
#          Select between spatial lag and spatial error via LM tests
#          Establish benchmark for comparison with zero-inflated spatial model
# Covariates: 2012 Census of Agriculture (beginning-of-period)
#             pasture_bop: 2011 NLCD pasture share (beginning-of-period)
#             Consistent beginning-of-period framework throughout
# Author: Lan T. Tran
# Date: June 2026

library(sf)
library(spdep)
library(spatialreg)
library(dplyr)

# ── LOAD DATA ─────────────────────────────────────────────────────────────────

counties_mo_proj <- st_read("data/processed/counties_mo.shp") |>
  arrange(GEOID)

panel   <- read.csv("data/processed/panel_final.csv")
w_queen <- readRDS("data/processed/weights_queen.rds")

# ── PREPARE ANALYTICAL SAMPLE ─────────────────────────────────────────────────
# Cross-section: 2021 (most recent NLCD year)
# Exclude St. Louis City — missing economic covariates
# year==2021 rows carry 2012 Census economic data (beginning-of-period)
# following Script 03 nlcd_year remapping

df <- panel |>
  filter(year == 2021) |>
  filter(!is.na(n_farms)) |>
  arrange(GEOID)

# Confirm: 114 rows (115 counties minus St. Louis City)
nrow(df)

# ── MERGE BEGINNING-OF-PERIOD PASTURE ─────────────────────────────────────────
# pasture_bop: 2011 NLCD pasture share — beginning of 2011-2021 transition period
# Consistent with 2012 Census economic covariates
# Transitioning counties had nearly twice the pasture share in 2011 (0.459 vs 0.250)
# Lark et al. (2015): grasslands constituted 77% of new cropland nationally 2008-2012
# Lark et al. (2020): Missouri identified as cropland expansion hotspot
# Lubowski, Plantinga, and Stavins (2008): pasture-to-cropland as competing returns decision

pasture_2011 <- panel |>
  filter(year == 2011, !is.na(n_farms)) |>
  select(GEOID, pasture) |>
  rename(pasture_bop = pasture) |>
  mutate(GEOID = as.character(GEOID))

df <- df |>
  mutate(GEOID = as.character(GEOID)) |>
  left_join(pasture_2011, by = "GEOID")

# Confirm no missing values after merge
sum(is.na(df$pasture_bop))

# Confirm no missing values in analytical variables
df |>
  select(transition, land_value_acre, net_income,
         n_farms, cropland, acres_operated, pasture_bop) |>
  summarise(across(everything(), ~sum(is.na(.))))

# ── STANDARDIZE COVARIATES ────────────────────────────────────────────────────
# Necessary to avoid numerical singularity in spatial regression
# Variables on very different scales (net_income in millions, transition 0/1)
# Standardization does not affect Rho, AIC, or residual diagnostics
# Coefficients interpreted as effect of one SD increase in each covariate
# pasture_z: standardized 2011 pasture share (beginning-of-period)

df <- df |>
  mutate(
    land_value_z = as.numeric(scale(land_value_acre)),
    net_income_z = as.numeric(scale(net_income)),
    n_farms_z    = as.numeric(scale(n_farms)),
    cropland_z   = as.numeric(scale(cropland)),
    acres_z      = as.numeric(scale(acres_operated)),
    pasture_z    = as.numeric(scale(pasture_bop))
  )

# ── SUBSET WEIGHTS TO MATCH ANALYTICAL SAMPLE ─────────────────────────────────
# Remove St. Louis City from weights matrix to match 114-row dataset

geoids_114 <- df$GEOID
w_114 <- subset(w_queen,
                subset = counties_mo_proj$GEOID %in% geoids_114)

# ── OLS BASELINE ──────────────────────────────────────────────────────────────
# Linear probability model — ignores spatial dependence
# Benchmark only: OLS residuals confirm severe spatial autocorrelation
# Covariate set includes pasture_z for consistency with ZI Part 2 specification

ols_z <- lm(transition ~ land_value_z + net_income_z +
              n_farms_z + cropland_z + acres_z + pasture_z,
            data = df)

summary(ols_z)

# ── LAGRANGE MULTIPLIER TESTS ─────────────────────────────────────────────────
# Selects between spatial lag and spatial error specification
# Decision rule: if both RSlag and RSerr significant, use robust versions
# adjRSlag significant, adjRSerr not → spatial lag preferred

lm_tests <- lm.RStests(ols_z, w_114,
                       test = c("RSlag", "RSerr", "adjRSlag", "adjRSerr"))
print(lm_tests)

# ── SPATIAL LAG MODEL ─────────────────────────────────────────────────────────
# Primary baseline model: spatial lag selected by robust LM tests
# Rho captures apparent spatial contagion in full 114-county sample
# Includes structural zeros — benchmark for Rho reduction in ZI model
# Note: Rho here reflects both genuine contagion AND spatial clustering
# of structural zeros and pasture availability — will be decomposed in Script 09

slag <- lagsarlm(transition ~ land_value_z + net_income_z +
                   n_farms_z + cropland_z + acres_z + pasture_z,
                 data  = df,
                 listw = w_114)

summary(slag)

# ── SPATIAL ERROR MODEL ───────────────────────────────────────────────────────
# Estimated for robustness comparison
# Lambda captures spatially correlated unobservables

serr <- errorsarlm(transition ~ land_value_z + net_income_z +
                     n_farms_z + cropland_z + acres_z + pasture_z,
                   data  = df,
                   listw = w_114)

summary(serr)

# ── MODEL COMPARISON ──────────────────────────────────────────────────────────

cat("\nAIC comparison:\n")
cat("OLS:          ", AIC(ols_z), "\n")
cat("Spatial lag:  ", AIC(slag), "\n")
cat("Spatial error:", AIC(serr), "\n")

# ── RESIDUAL DIAGNOSTICS ──────────────────────────────────────────────────────

moran_ols_resid  <- moran.test(residuals(ols_z), w_114)
moran_slag_resid <- moran.test(residuals(slag),  w_114)
moran_serr_resid <- moran.test(residuals(serr),  w_114)

cat("\nMoran's I on residuals:\n")
cat("OLS:          ", round(moran_ols_resid$estimate[1], 3),
    " p =", round(moran_ols_resid$p.value, 6), "\n")
cat("Spatial lag:  ", round(moran_slag_resid$estimate[1], 3),
    " p =", round(moran_slag_resid$p.value, 6), "\n")
cat("Spatial error:", round(moran_serr_resid$estimate[1], 3),
    " p =", round(moran_serr_resid$p.value, 6), "\n")

# ── SAVE MODELS AND RESULTS ───────────────────────────────────────────────────

saveRDS(ols_z, "data/processed/model_ols.rds")
saveRDS(slag,  "data/processed/model_slag.rds")
saveRDS(serr,  "data/processed/model_serr.rds")

sink("docs/standard_model_results.txt")
cat("=== OLS BASELINE ===\n")
cat("Covariates: 2012 Census economic data + 2011 NLCD pasture (beginning-of-period)\n\n")
print(summary(ols_z))
cat("\n=== SPATIAL LAG MODEL ===\n")
print(summary(slag))
cat("\n=== SPATIAL ERROR MODEL ===\n")
print(summary(serr))
cat("\n=== MODEL COMPARISON ===\n")
cat("AIC - OLS:          ", AIC(ols_z), "\n")
cat("AIC - Spatial lag:  ", AIC(slag), "\n")
cat("AIC - Spatial error:", AIC(serr), "\n")
cat("\nMoran's I on residuals:\n")
cat("OLS:          ", round(moran_ols_resid$estimate[1], 3),
    " p =", round(moran_ols_resid$p.value, 6), "\n")
cat("Spatial lag:  ", round(moran_slag_resid$estimate[1], 3),
    " p =", round(moran_slag_resid$p.value, 6), "\n")
cat("Spatial error:", round(moran_serr_resid$estimate[1], 3),
    " p =", round(moran_serr_resid$p.value, 6), "\n")
sink()

message("Script 07 complete - standard spatial models with pasture_bop estimated and saved")