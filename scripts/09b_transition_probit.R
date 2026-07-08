# Script 09b: Part 2, the transition probit
# Purpose: model which active counties expanded cropland between 2011 and 2021,
#          and show whether apparent spatial dependence is really spatially
#          clustered convertible land (pasture).
#
# Sample: 90 active counties. The active set is 91 counties (114 minus 23
#         structural zeros). The island, Ste. Genevieve, is dropped because the
#         spatial probit rejects a no-neighbour row. State this one county
#         difference in the paper.
#
# Estimator: maximum likelihood spatial lag probit (ProbitSpatial, DGP SAR,
#            method conditional as primary, full-lik as a check).
#
# Specifications:
#   M1 primary:   transition ~ cropland_z + pasture_z
#   M0 contrast:  transition ~ cropland_z   (no pasture, for the spatial
#                 parameter contrast only; its coefficients are not a
#                 standalone model because removing pasture removes the signal)
#   M2 landvalue: transition ~ cropland_z + pasture_z + landvalue_z
#                 (robustness, justified on theory, land value correlates
#                  0.669 with cropland, read cropland stability when it enters)
#
# Author: Lan T. Tran
# Date: July 2026

library(sf)
library(spdep)
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

# ── BUILD THE ACTIVE SAMPLE (SAME LOGIC AS 07, 08, 09a) ──────────────────────

df <- panel |>
  filter(year == 2021, !is.na(n_farms)) |>
  arrange(GEOID) |>
  left_join(pasture |> select(GEOID, pasture_2011), by = "GEOID")

geoids_114 <- df$GEOID
w_114 <- subset(w_queen, subset = counties_mo_proj$GEOID %in% geoids_114)

lisa_cropland <- localmoran(df$cropland, w_114)
mean_crop     <- mean(df$cropland)
lag_crop      <- lag.listw(w_114, df$cropland)
df$structural_zero <- as.integer(
  df$cropland < mean_crop & lag_crop < mean_crop & lisa_cropland[, 5] < 0.05
)

active_idx <- df$structural_zero == 0
active     <- df[active_idx, ]
w_active   <- subset(w_114, subset = active_idx)

cat("Active counties (with island):", nrow(active), "\n")

# ── CONFIRM NO MISSING ON MODEL VARIABLES ────────────────────────────────────

for (v in c("cropland_lag", "pasture_2011", "land_value_acre")) {
  if (!v %in% names(active)) stop("Missing column: ", v)
  cat("Missing", v, ":", sum(is.na(active[[v]])), "\n")
}

# ── STANDARDIZE ──────────────────────────────────────────────────────────────
# cropland_z is built from cropland_lag (2011 level), the initial input.

active <- active |>
  mutate(
    cropland_z  = as.numeric(scale(cropland_lag)),
    pasture_z   = as.numeric(scale(pasture_2011)),
    landvalue_z = as.numeric(scale(land_value_acre))
  )

# ── BUILD THE WEIGHT MATRIX AND DROP THE ISLAND ──────────────────────────────
# ProbitSpatialFit needs a row normalized dgCMatrix and rejects a zero row.

W_dense <- listw2mat(w_active)
W_dense[is.na(W_dense)] <- 0
rs <- rowSums(W_dense)
n_islands <- sum(abs(rs) < 1e-12)
cat("\nIslands in the active matrix:", n_islands, "\n")

keep <- which(rs > 1e-12)
dropped <- nrow(active) - length(keep)
active_90 <- active[keep, ]
W90_dense <- W_dense[keep, keep]
new_rs <- rowSums(W90_dense)
if (any(abs(new_rs) < 1e-12)) {
  stop("Dropping the island stranded another county. Investigate.")
}
W90_dense <- W90_dense / new_rs
W90 <- as(W90_dense, "CsparseMatrix")

cat("Dropped", dropped, "island county. Probit sample n =", nrow(active_90), "\n")
if (!is(W90, "dgCMatrix")) stop("W90 is not dgCMatrix.")

# ── FIT HELPER ───────────────────────────────────────────────────────────────

fit_probit <- function(label, formula, method = "conditional") {
  cat("\n--------------------------------------------------------------\n")
  cat(label, " [method =", method, "]\n")
  cat("--------------------------------------------------------------\n")
  fit <- tryCatch(
    ProbitSpatialFit(formula, data = active_90, W = W90,
                     DGP = "SAR", method = method, varcov = "varcov"),
    error   = function(e) { cat("FIT ERROR:", conditionMessage(e), "\n"); NULL },
    warning = function(w) { cat("FIT WARNING:", conditionMessage(w), "\n")
      suppressWarnings(ProbitSpatialFit(formula, data = active_90, W = W90,
                     DGP = "SAR", method = method, varcov = "varcov")) }
  )
  if (!is.null(fit)) print(summary(fit))
  fit
}

get_rho <- function(fit) {
  if (is.null(fit)) return(NA_real_)
  as.numeric(tryCatch(fit@rho, error = function(e) NA_real_))
}

# ── M1 PRIMARY: cropland + pasture ───────────────────────────────────────────

m1 <- fit_probit("M1 PRIMARY  transition ~ cropland + pasture",
                 transition ~ cropland_z + pasture_z)

# ── M0 CONTRAST: cropland only (spatial parameter contrast only) ─────────────

m0 <- fit_probit("M0 CONTRAST  transition ~ cropland (no pasture)",
                 transition ~ cropland_z)
cat("\nNote: M0 coefficients are not a standalone model. Removing pasture\n")
cat("removes the signal. M0 exists only to show the spatial parameter contrast.\n")

# ── M2 ROBUSTNESS: add land value ────────────────────────────────────────────

m2 <- fit_probit("M2 ROBUSTNESS  transition ~ cropland + pasture + land value",
                 transition ~ cropland_z + pasture_z + landvalue_z)

# ── FULL-LIK CROSS CHECK ON THE PRIMARY ──────────────────────────────────────

m1_full <- fit_probit("M1 PRIMARY, full-lik cross check",
                      transition ~ cropland_z + pasture_z, method = "full-lik")

# ── SPATIAL PARAMETER SUMMARY ────────────────────────────────────────────────

cat("\n==============================================================\n")
cat("SPATIAL PARAMETER (rho) ACROSS SPECIFICATIONS\n")
cat("==============================================================\n")
cat("M0 no pasture           :", round(get_rho(m0), 3), "\n")
cat("M1 primary, conditional :", round(get_rho(m1), 3), "\n")
cat("M1 primary, full-lik    :", round(get_rho(m1_full), 3), "\n")
cat("M2 with land value      :", round(get_rho(m2), 3), "\n")
cat("\nPasture contrast (M0 minus M1):",
    round(get_rho(m0) - get_rho(m1), 3), "\n")
cat("This movement is the finding. Watch whether it holds in M2.\n")
cat("If auto extraction shows NA, read rho from the summary tables above.\n")

cat("\nWatch also: does the cropland coefficient stay stable from M1 to M2.\n")
cat("A large shift signals the land value collinearity (0.669) biting.\n")

# ── SAVE ─────────────────────────────────────────────────────────────────────

saveRDS(m1, "data/processed/model_part2_primary.rds")
saveRDS(m2, "data/processed/model_part2_landvalue.rds")

sink("docs/part2_results.txt")
cat("=== PART 2 TRANSITION PROBIT ===\n")
cat("Sample: 90 active counties (island dropped). Estimator: ProbitSpatial SAR.\n\n")
cat("M1 PRIMARY (cropland + pasture):\n");            print(summary(m1))
cat("\nM2 ROBUSTNESS (+ land value):\n");              print(summary(m2))
cat("\nSpatial parameter across specifications:\n")
cat("M0 no pasture:", round(get_rho(m0), 3),
    " M1:", round(get_rho(m1), 3),
    " M1 full-lik:", round(get_rho(m1_full), 3),
    " M2:", round(get_rho(m2), 3), "\n")
sink()

message("Script 09b complete. Part 2 transition probit estimated and saved.")
