# Script 07: Standard Spatial Econometric Regression
# Purpose: Estimate OLS, spatial lag, and spatial error models as the
#          naive baseline. This is standard practice, the specification
#          an applied researcher would run without knowing about the
#          structural zero problem. It establishes the Rho value that
#          Script 09's zero inflated model is compared against.
#
# Covariates: land value, net income, number of farms, cropland level,
# acres operated. All from the 2012 Census of Agriculture, beginning of
# period relative to the 2011 to 2021 transition window.
#
# Pasture and CRP are not used here. Both are suspended per the variable
# selection note, pending further validation. This script does not
# require either column to exist in panel_final.csv.
#
# Author: Lan T. Tran
# Date: June 2026

library(sf)
library(spdep)
library(spatialreg)
library(dplyr)

# ── LOAD DATA ─────────────────────────────────────────────────────────────────

counties_mo_proj <- st_read("data/processed/counties_mo.shp") |>
  arrange(GEOID)

panel <- read.csv("data/processed/panel_final.csv")

# weights_queen.rds is expected from an earlier script that builds the
# spatial weights matrix. That script has not been reviewed in this
# conversation, so its contents are not confirmed. Check the object type
# before trusting anything built from it.
w_queen <- readRDS("data/processed/weights_queen.rds")

class(w_queen)
if (!inherits(w_queen, "listw")) {
  stop("w_queen is not a listw object. Check how weights_queen.rds was built.")
}

# ── PREPARE ANALYTICAL SAMPLE ─────────────────────────────────────────────────
# 2021 cross section, excludes St. Louis City, no farm data.

df <- panel |>
  filter(year == 2021, !is.na(n_farms)) |>
  arrange(GEOID)

n_sample <- nrow(df)
message("Analytical sample size: ", n_sample, " counties")
if (n_sample != 114) {
  warning("Sample size is not 114. Check the year and n_farms filter.")
}

dup_check <- df |> count(GEOID) |> filter(n > 1)
if (nrow(dup_check) > 0) {
  warning("Duplicate GEOID in analytical sample. Stop and investigate.")
}

missing_check <- df |>
  select(transition, land_value_acre, net_income, n_farms,
         cropland, acres_operated) |>
  summarise(across(everything(), ~sum(is.na(.))))

cat("\nMissing value check on analytical variables:\n")
print(missing_check)

if (any(missing_check > 0)) {
  warning("At least one analytical variable has missing values. ",
          "lm and lagsarlm will silently drop these rows.")
}

# ── STANDARDIZE COVARIATES ────────────────────────────────────────────────────
# cropland_z here is built from the raw 2021 cropland level, not from
# cropland_lag. This is a specific choice, not an oversight. This script
# represents naive standard practice, the version that does not yet
# know about the timing correction. The corrected version, using
# cropland_lag, belongs in Script 09's Part 2.
#
# This choice has not been confirmed with you directly. Flagging it here
# so it is a decision you make, not one I made quietly on your behalf.

df <- df |>
  mutate(
    land_value_z = as.numeric(scale(land_value_acre)),
    net_income_z = as.numeric(scale(net_income)),
    n_farms_z    = as.numeric(scale(n_farms)),
    cropland_z   = as.numeric(scale(cropland)),
    acres_z      = as.numeric(scale(acres_operated))
  )

# ── SUBSET WEIGHTS TO MATCH ANALYTICAL SAMPLE ─────────────────────────────────

geoids_114 <- df$GEOID
w_114 <- subset(w_queen, subset = counties_mo_proj$GEOID %in% geoids_114)

n_weights <- length(w_114$neighbours)
message("Weights object covers ", n_weights, " counties")
if (n_weights != n_sample) {
  stop("Weights object and analytical sample do not match in size.")
}

# ── MODEL FORMULA ─────────────────────────────────────────────────────────────

model_formula <- transition ~ land_value_z + net_income_z +
  n_farms_z + cropland_z + acres_z

# ── OLS BASELINE ───────────────────────────────────────────────────────────────

ols_baseline <- lm(model_formula, data = df)
summary(ols_baseline)

# ── LAGRANGE MULTIPLIER TESTS ──────────────────────────────────────────────────
# Used to decide between spatial lag and spatial error. Read the output,
# this script does not pick one for you.

lm_tests <- lm.RStests(ols_baseline, w_114,
                       test = c("RSlag", "RSerr", "adjRSlag", "adjRSerr"))
summary(lm_tests)

# ── SPATIAL LAG MODEL ──────────────────────────────────────────────────────────

model_slag <- lagsarlm(model_formula, data = df, listw = w_114)
summary(model_slag)

# ── SPATIAL ERROR MODEL ────────────────────────────────────────────────────────

model_serr <- errorsarlm(model_formula, data = df, listw = w_114)
summary(model_serr)

# ── MODEL COMPARISON ───────────────────────────────────────────────────────────

cat("\nModel fit comparison:\n")
cat("OLS AIC:           ", round(AIC(ols_baseline), 2), "\n")
cat("Spatial lag AIC:   ", round(AIC(model_slag), 2), "\n")
cat("Spatial error AIC: ", round(AIC(model_serr), 2), "\n")

cat("\nSpatial lag Rho:   ", round(model_slag$rho, 3),
    " p value:", round(summary(model_slag)$Wald1$p.value, 4), "\n")

# ── RESIDUAL DIAGNOSTICS ──────────────────────────────────────────────────────

moran_ols_resid  <- moran.test(residuals(ols_baseline), w_114)
moran_slag_resid <- moran.test(residuals(model_slag), w_114)
moran_serr_resid <- moran.test(residuals(model_serr), w_114)

cat("\nMoran's I on residuals:\n")
cat("OLS:           ", round(moran_ols_resid$estimate[1], 3),
    " p =", round(moran_ols_resid$p.value, 4), "\n")
cat("Spatial lag:   ", round(moran_slag_resid$estimate[1], 3),
    " p =", round(moran_slag_resid$p.value, 4), "\n")
cat("Spatial error: ", round(moran_serr_resid$estimate[1], 3),
    " p =", round(moran_serr_resid$p.value, 4), "\n")

# ── SAVE MODELS AND RESULTS ────────────────────────────────────────────────────
# Saved as model_slag.rds. Script 09 loads this exact file name for its
# Rho and AIC comparison. Confirm this matches before running Script 09.

saveRDS(ols_baseline, "data/processed/model_ols.rds")
saveRDS(model_slag,   "data/processed/model_slag.rds")
saveRDS(model_serr,   "data/processed/model_serr.rds")

sink("docs/standard_model_results.txt")
cat("=== STANDARD SPATIAL REGRESSION, NAIVE BASELINE ===\n")
cat("N = ", n_sample, " counties, structural zeros included\n")
cat("No pasture or CRP, both suspended\n\n")

cat("=== OLS ===\n")
print(summary(ols_baseline))

cat("\n=== LAGRANGE MULTIPLIER TESTS ===\n")
print(summary(lm_tests))

cat("\n=== SPATIAL LAG ===\n")
print(summary(model_slag))

cat("\n=== SPATIAL ERROR ===\n")
print(summary(model_serr))

cat("\n=== MODEL FIT COMPARISON ===\n")
cat("OLS AIC:           ", round(AIC(ols_baseline), 2), "\n")
cat("Spatial lag AIC:   ", round(AIC(model_slag), 2), "\n")
cat("Spatial error AIC: ", round(AIC(model_serr), 2), "\n")

cat("\n=== RESIDUAL DIAGNOSTICS ===\n")
cat("OLS Moran's I:           ", round(moran_ols_resid$estimate[1], 3),
    " p =", round(moran_ols_resid$p.value, 4), "\n")
cat("Spatial lag Moran's I:   ", round(moran_slag_resid$estimate[1], 3),
    " p =", round(moran_slag_resid$p.value, 4), "\n")
cat("Spatial error Moran's I: ", round(moran_serr_resid$estimate[1], 3),
    " p =", round(moran_serr_resid$p.value, 4), "\n")
sink()

message("Script 07 complete. Naive baseline model estimated and saved.")