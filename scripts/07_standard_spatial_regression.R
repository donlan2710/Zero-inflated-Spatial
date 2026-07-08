# Script 07: Standard Spatial Baseline
# Purpose: Estimate the standard spatial model on all 114 counties, the
#          version that pools structural zeros with active counties. This is
#          the baseline the two-part model is compared against. The comparison
#          is clean: the baseline uses the SAME two variables, the SAME timing,
#          and the SAME estimator as the Part 2 model. The only difference
#          between this and Part 2 is that here the 23 structural zeros are
#          left in, not removed.
#
# Variables: cropland_lag (2011 initial level) and pasture_2011. Standardized.
#            cropland_z is built from cropland_lag, NOT from the 2021 level.
#            The old version of this script used the 2021 level, which was an
#            error, and it is corrected here.
#
# Estimator: tries the maximum likelihood spatial probit first (ProbitSpatial,
#            SAR). On the full 114 sample the 23 structural zeros are all zeros
#            with no transitions, which may make the probit struggle. If the
#            probit does not converge cleanly, the script runs the linear
#            spatial model (lagsarlm) as the baseline and reports which one to
#            use, so the baseline-to-active comparison stays on one estimator.
#
# Author: Lan T. Tran
# Date: July 2026

library(sf)
library(spdep)
library(spatialreg)
library(dplyr)
library(Matrix)
library(ProbitSpatial)

# ── LOAD DATA ─────────────────────────────────────────────────────────────────

counties_mo_proj <- st_read("data/processed/counties_mo.shp") |>
  arrange(GEOID)

panel   <- read.csv("data/processed/panel_final.csv",
                    colClasses = c(GEOID = "character"))
pasture <- read.csv("data/processed/pasture_by_year.csv",
                    colClasses = c(GEOID = "character"))
w_queen <- readRDS("data/processed/weights_queen.rds")

if (!inherits(w_queen, "listw")) {
  stop("w_queen is not a listw object. Check how weights_queen.rds was built.")
}

# ── ANALYTICAL SAMPLE, 114 COUNTIES ──────────────────────────────────────────
# 2021 cross section, drop St. Louis City on missing farm data.

df <- panel |>
  filter(year == 2021, !is.na(n_farms)) |>
  arrange(GEOID) |>
  left_join(pasture |> select(GEOID, pasture_2011), by = "GEOID")

n_sample <- nrow(df)
message("Analytical sample size: ", n_sample, " counties")
if (n_sample != 114) {
  warning("Sample size is not 114. Check the year and n_farms filter.")
}

# Confirm no missing values on the two model variables.
miss_crop <- sum(is.na(df$cropland_lag))
miss_past <- sum(is.na(df$pasture_2011))
cat("Missing cropland_lag:", miss_crop, "  Missing pasture_2011:", miss_past, "\n")
if (miss_crop > 0 || miss_past > 0) {
  stop("A model variable has missing values. Stop and investigate.")
}

# ── STANDARDIZE THE TWO VARIABLES ────────────────────────────────────────────
# cropland_z is built from cropland_lag (2011 level), the initial input.

df <- df |>
  mutate(
    cropland_z = as.numeric(scale(cropland_lag)),
    pasture_z  = as.numeric(scale(pasture_2011))
  )

cat("cropland_z is built from cropland_lag (2011 level), confirmed by design.\n")

# ── SUBSET WEIGHTS TO THE 114 SAMPLE ─────────────────────────────────────────

geoids_114 <- df$GEOID
w_114 <- subset(w_queen, subset = counties_mo_proj$GEOID %in% geoids_114)

n_weights <- length(w_114$neighbours)
message("Weights object covers ", n_weights, " counties")
if (n_weights != n_sample) {
  stop("Weights object and analytical sample do not match in size.")
}

model_formula <- transition ~ cropland_z + pasture_z

# ── ATTEMPT 1: SPATIAL PROBIT ON 114 ─────────────────────────────────────────
# Build a row normalized weight matrix as dgCMatrix. Check for islands first,
# because ProbitSpatialFit rejects a no-neighbour row.

W_dense <- listw2mat(w_114)
W_dense[is.na(W_dense)] <- 0
rs <- rowSums(W_dense)
n_islands <- sum(abs(rs) < 1e-12)
cat("\nIslands (zero rows) in the 114 weight matrix:", n_islands, "\n")

probit_ok <- FALSE
model_probit <- NULL

if (n_islands == 0) {
  W_dgc <- as(W_dense, "CsparseMatrix")
  model_probit <- tryCatch(
    ProbitSpatialFit(model_formula, data = df, W = W_dgc,
                     DGP = "SAR", method = "conditional", varcov = "varcov"),
    error   = function(e) { cat("PROBIT ERROR:", conditionMessage(e), "\n"); NULL },
    warning = function(w) { cat("PROBIT WARNING:", conditionMessage(w), "\n")
      suppressWarnings(
        ProbitSpatialFit(model_formula, data = df, W = W_dgc,
                         DGP = "SAR", method = "conditional", varcov = "varcov")) }
  )
  if (!is.null(model_probit)) {
    cat("\n=== SPATIAL PROBIT ON 114 ===\n")
    print(summary(model_probit))
    rho_probit <- tryCatch(model_probit@rho, error = function(e) NA_real_)
    cat("\nProbit spatial parameter (rho):", round(as.numeric(rho_probit), 3), "\n")
    if (!is.na(rho_probit) && abs(rho_probit) < 0.999) probit_ok <- TRUE
  }
} else {
  cat("The 114 matrix has an island. The probit cannot take it here.\n")
}

if (probit_ok) {
  cat("\nProbit converged on 114. Use it as the baseline, matching the Part 2\n")
  cat("estimator. The baseline-to-active comparison is then probit vs probit.\n")
} else {
  cat("\nProbit did not give a clean baseline on 114. Falling back to the\n")
  cat("linear spatial model below. The baseline-to-active comparison then runs\n")
  cat("in the linear model, which handles both samples.\n")
}

# ── ATTEMPT 2 / FALLBACK: LINEAR SPATIAL MODEL ON 114 ────────────────────────
# lagsarlm handles the full sample and the island via zero.policy. This is the
# baseline if the probit did not converge, and it is a useful support number
# either way.

model_slag <- lagsarlm(model_formula, data = df, listw = w_114,
                       zero.policy = TRUE)

cat("\n=== LINEAR SPATIAL MODEL ON 114 (lagsarlm) ===\n")
print(summary(model_slag))

cat("\nLinear spatial parameter (rho):", round(model_slag$rho, 3), "\n")
fit <- fitted(model_slag)
cat("Fitted range:", round(min(fit), 3), "to", round(max(fit), 3),
    "  outside [0,1]:", sum(fit < 0 | fit > 1), "of", length(fit), "\n")

moran_resid <- moran.test(residuals(model_slag), w_114, zero.policy = TRUE)
cat("Moran's I on residuals:", round(moran_resid$estimate[1], 3),
    " p =", format.pval(moran_resid$p.value, digits = 3), "\n")

# ── SAVE ─────────────────────────────────────────────────────────────────────
# Save whichever baseline is valid. Part 2 will load the matching estimator.

if (probit_ok) saveRDS(model_probit, "data/processed/model_baseline_probit_114.rds")
saveRDS(model_slag, "data/processed/model_baseline_linear_114.rds")

sink("docs/baseline_114_results.txt")
cat("=== STANDARD SPATIAL BASELINE, 114 COUNTIES ===\n")
cat("Variables: cropland_lag (2011 level) + pasture_2011, standardized\n")
cat("Structural zeros are pooled in. This is the baseline for the two-part model.\n\n")
if (probit_ok) {
  cat("=== SPATIAL PROBIT (baseline estimator) ===\n")
  print(summary(model_probit))
  cat("\n")
} else {
  cat("Probit did not converge cleanly on 114. Linear model is the baseline.\n\n")
}
cat("=== LINEAR SPATIAL MODEL ===\n")
print(summary(model_slag))
cat("\nLinear rho:", round(model_slag$rho, 3), "\n")
cat("Moran's I residuals:", round(moran_resid$estimate[1], 3),
    " p =", format.pval(moran_resid$p.value, digits = 3), "\n")
sink()

message("Script 07 complete. Baseline estimated on 114 counties.")
