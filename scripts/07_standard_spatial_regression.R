# Script 07: Standard Spatial Econometric Regression
# Purpose: Estimate OLS and spatial regression models as baseline
#          Select between spatial lag and spatial error via LM tests
#          Establish benchmark for comparison with zero-inflated spatial model
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
# Exclude 2001 — no transition variable (first period has no lag)

df <- panel |>
  filter(year == 2021) |>
  filter(!is.na(n_farms)) |>
  arrange(GEOID)

# Confirm: 114 rows (115 counties minus St. Louis City)
nrow(df)

# Confirm no missing values in analytical variables
df |>
  select(transition, land_value_acre, net_income,
         n_farms, cropland, acres_operated) |>
  summarise(across(everything(), ~sum(is.na(.))))

# ── STANDARDIZE COVARIATES ────────────────────────────────────────────────────
# Necessary to avoid numerical singularity in spatial regression
# Variables on very different scales (net_income in millions, transition 0/1)
# Standardization does not affect Rho, AIC, or residual diagnostics
# Coefficients interpreted as effect of one SD increase in each covariate

df <- df |>
  mutate(
    land_value_z = scale(land_value_acre),
    net_income_z = scale(net_income),
    n_farms_z    = scale(n_farms),
    cropland_z   = scale(cropland),
    acres_z      = scale(acres_operated)
  )

# ── SUBSET WEIGHTS TO MATCH ANALYTICAL SAMPLE ─────────────────────────────────
# Remove St. Louis City from weights matrix to match 114-row dataset

geoids_114 <- df$GEOID
w_114 <- subset(w_queen,
                subset = counties_mo_proj$GEOID %in% geoids_114)

# ── OLS BASELINE ──────────────────────────────────────────────────────────────
# Linear probability model — ignores spatial dependence
# Benchmark only: OLS residuals show Moran's I = 0.471 (p < 0.001)
# confirming severe spatial autocorrelation

ols_z <- lm(transition ~ land_value_z + net_income_z +
              n_farms_z + cropland_z + acres_z,
            data = df)

summary(ols_z)

# ── LAGRANGE MULTIPLIER TESTS ─────────────────────────────────────────────────
# Selects between spatial lag and spatial error specification
# Decision rule: if both LMlag and LMerr significant, use robust versions
# adjRSlag p = 0.086, adjRSerr p = 0.995 → spatial lag preferred

lm_tests <- lm.LMtests(ols_z, w_114,
                       test = c("LMlag", "LMerr", "RLMlag", "RLMerr"))
print(lm_tests)

# ── SPATIAL LAG MODEL ─────────────────────────────────────────────────────────
# Primary model: spatial lag selected by robust LM tests
# Rho captures spatial contagion — transition probability influenced
# by whether neighboring counties transitioned

slag <- lagsarlm(transition ~ land_value_z + net_income_z +
                   n_farms_z + cropland_z + acres_z,
                 data = df,
                 listw = w_114)

summary(slag)
# Key result: Rho = 0.673, p < 2.2e-16 — strong spatial contagion
# AIC = 89.25 vs OLS 135.76 — large improvement

# ── SPATIAL ERROR MODEL ───────────────────────────────────────────────────────
# Estimated for robustness comparison
# Lambda captures spatially correlated unobservables
# AIC = 88.14 — marginally better than spatial lag but difference < 2
# Both models tell the same substantive story

serr <- errorsarlm(transition ~ land_value_z + net_income_z +
                     n_farms_z + cropland_z + acres_z,
                   data = df,
                   listw = w_114)

summary(serr)

# ── MODEL COMPARISON ──────────────────────────────────────────────────────────

cat("\nAIC comparison:\n")
cat("OLS:          ", AIC(ols_z), "\n")
cat("Spatial lag:  ", AIC(slag), "\n")
cat("Spatial error:", AIC(serr), "\n")

# ── RESIDUAL DIAGNOSTICS ──────────────────────────────────────────────────────
# Spatial models should remove autocorrelation from residuals
# OLS Moran's I = 0.471 (p = 0) — severe autocorrelation
# Spatial lag Moran's I = 0.012 (p = 0.362) — fully resolved
# Spatial error Moran's I = 0.004 (p = 0.415) — fully resolved

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

message("Script 07 complete - standard spatial models estimated and saved")