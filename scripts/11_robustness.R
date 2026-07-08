# Script 11: Robustness Checks — Zero-Inflated Spatial Model
# Purpose: Test sensitivity of Option A findings ordered from least to most invasive
#
#   Sample sensitivity:
#     Check 1 — Ste. Genevieve excluded (GEOID 29186, island county)
#
#   Outcome definition:
#     Check 2 — Transition threshold 0.03
#
#   Spatial structure:
#     Check 3 — Distance 60km weights (active subsample geometry)
#     Check 4 — Distance 80km weights (active subsample geometry)
#
#   Identification assumption — structural zero classification:
#     Check 5 — LISA p < 0.10 (looser threshold, more structural zeros)
#     Check 6 — LISA p < 0.01 (stricter threshold, fewer structural zeros)
#
#   Covariate specification (Part 2 only):
#     Check 7 — Drop net_income_z (marginal predictor, p = 0.094)
#     Check 8 — Drop acres_z (not significant, p = 0.411)
#     Check 9 — Add pasture (theoretically motivated, already in data)
#
# All checks use 2021 cross-section with 2012 Census economic covariates
# Baseline: Rho = 0.551, Moran's I p = 0.293, n=91 active counties
# Author: Lan T. Tran
# Date: June 2026

library(sf)
library(spdep)
library(spatialreg)
library(dplyr)

# ── LOAD DATA ─────────────────────────────────────────────────────────────────

counties_mo_proj <- st_read("data/processed/counties_mo.shp") |>
  arrange(GEOID)

panel    <- read.csv("data/processed/panel_final.csv")
nb_queen <- readRDS("data/processed/nb_queen.rds")
w_queen  <- readRDS("data/processed/weights_queen.rds")

# ── ANALYTICAL SAMPLE ─────────────────────────────────────────────────────────

df_2021 <- panel |>
  filter(year == 2021, !is.na(n_farms)) |>
  arrange(GEOID)

# ── DISTANCE WEIGHTS — ACTIVE SUBSAMPLE GEOMETRY ─────────────────────────────
# Built on 91-county active subsample, not full 115-county sample
# Full connectivity never achieved — Ste. Genevieve permanently isolated
# Connectivity on active subsample:
#   Distance 45km: 250 links, 5 subgraphs
#   Distance 50km: 346 links, 3 subgraphs
#   Distance 55km: 396 links, 2 subgraphs
#   Distance 60km: 444 links, 2 subgraphs
#   Distance 70km: 576 links, 2 subgraphs
#   Distance 80km: 770 links, 2 subgraphs
# 60km and 80km bracket the range — zero.policy = TRUE throughout

coords_active <- st_coordinates(
  st_centroid(
    counties_mo_proj |>
      filter(GEOID %in% df_2021$GEOID[{
        w_sub     <- subset(w_queen,
                            subset = counties_mo_proj$GEOID %in% df_2021$GEOID)
        lisa_temp <- localmoran(df_2021$cropland, w_sub)
        mean_crop <- mean(df_2021$cropland)
        lag_crop  <- lag.listw(w_sub, df_2021$cropland)
        df_2021$cropland < mean_crop &
          lag_crop < mean_crop &
          lisa_temp[, 5] < 0.05
      } == 0])
  )
)

nb_dist60 <- dnearneigh(coords_active, d1 = 0, d2 = 60000)
w_dist60  <- nb2listw(nb_dist60, style = "W", zero.policy = TRUE)

nb_dist80 <- dnearneigh(coords_active, d1 = 0, d2 = 80000)
w_dist80  <- nb2listw(nb_dist80, style = "W", zero.policy = TRUE)

# ── HELPER FUNCTION ───────────────────────────────────────────────────────────
# Arguments:
#   data_in            — panel filtered to relevant year and sample
#   threshold          — transition threshold (0.02, 0.03)
#   nb_full / w_full   — neighbor list and weights for Part 2
#   label              — row label in output table
#   use_prebuilt_active— TRUE for distance weights (already on active geometry)
#   lisa_p_threshold   — LISA p-value cutoff for structural zero classification
#   part2_formula      — formula for Part 2 spatial lag (allows covariate changes)

run_zi_spatial <- function(data_in, threshold, nb_full, w_full, label,
                           use_prebuilt_active = FALSE,
                           lisa_p_threshold    = 0.05,
                           part2_formula       = transition ~ land_value_z +
                             net_income_z + n_farms_z +
                             cropland_z + acres_z) {
  
  message("Running: ", label)
  
  # Redefine transition at specified threshold
  data_in <- data_in |>
    mutate(transition = as.integer(abs(cropland_change) > threshold))
  
  # LISA always uses queen weights for consistency across all checks
  geoids_in <- data_in$GEOID
  w_sub     <- subset(w_queen, subset = counties_mo_proj$GEOID %in% geoids_in)
  
  lisa      <- localmoran(data_in$cropland, w_sub)
  mean_crop <- mean(data_in$cropland)
  lag_crop  <- lag.listw(w_sub, data_in$cropland)
  
  data_in <- data_in |>
    mutate(
      forest_total    = forest + frac_42 + frac_43,
      structural_zero = as.integer(
        cropland < mean_crop &
          lag_crop < mean_crop &
          lisa[, 5] < lisa_p_threshold
      )
    )
  
  n_structural <- sum(data_in$structural_zero)
  n_transition <- sum(data_in$transition)
  zero_rate    <- round(100 * mean(data_in$transition == 0), 1)
  overlap      <- sum(data_in$structural_zero == 1 & data_in$transition == 1)
  
  # Active subsample — standardized within subsample
  data_active <- data_in |>
    filter(structural_zero == 0) |>
    mutate(
      land_value_z = as.numeric(scale(land_value_acre)),
      net_income_z = as.numeric(scale(net_income)),
      n_farms_z    = as.numeric(scale(n_farms)),
      cropland_z   = as.numeric(scale(cropland)),
      acres_z      = as.numeric(scale(acres_operated)),
      pasture_z    = as.numeric(scale(pasture))
    )
  
  # Weights for Part 2
  if (use_prebuilt_active) {
    nb_active <- nb_full
    w_active  <- w_full
  } else {
    nb_active <- subset(nb_full,
                        subset = counties_mo_proj$GEOID %in%
                          data_in$GEOID[data_in$structural_zero == 0])
    w_active  <- nb2listw(nb_active, style = "W", zero.policy = TRUE)
  }
  
  # Part 2 spatial lag
  model <- tryCatch(
    lagsarlm(part2_formula,
             data        = data_active,
             listw       = w_active,
             zero.policy = TRUE),
    error = function(e) {
      message("Model failed for ", label, ": ", e$message)
      return(NULL)
    }
  )
  
  if (is.null(model)) return(NULL)
  
  moran_resid <- moran.test(residuals(model), w_active, zero.policy = TRUE)
  coef_tab    <- summary(model)$Coef
  
  # Extract land value — primary significant covariate — for summary table
  lv_est <- round(coef_tab["land_value_z", 1], 4)
  lv_p   <- round(coef_tab["land_value_z", 4], 4)
  
  data.frame(
    Specification  = label,
    N_total        = nrow(data_in),
    N_structural   = n_structural,
    N_active       = nrow(data_active),
    N_transition   = n_transition,
    Zero_rate      = zero_rate,
    Overlap        = overlap,
    Rho            = round(model$rho, 4),
    LandValue_est  = lv_est,
    LandValue_p    = lv_p,
    MoranI_resid   = round(moran_resid$estimate[1], 4),
    MoranI_p       = round(moran_resid$p.value, 4),
    AIC            = round(AIC(model), 2),
    stringsAsFactors = FALSE
  )
}

# ── BASELINE ──────────────────────────────────────────────────────────────────

r_baseline <- run_zi_spatial(df_2021, 0.02, nb_queen, w_queen,
                             "Baseline — Queen p<0.05 threshold 0.02")

# ── CHECK 1: SAMPLE SENSITIVITY — STE. GENEVIEVE EXCLUDED ────────────────────
# GEOID 29186 — island county, no active neighbors after structural zero removal
# Confirms zero.policy = TRUE treatment does not drive Rho finding

df_no_stegen <- df_2021 |> filter(GEOID != "29186")

r_no_stegen <- run_zi_spatial(df_no_stegen, 0.02, nb_queen, w_queen,
                              "Check 1 — Ste. Genevieve excluded")

# ── CHECK 2: OUTCOME DEFINITION — THRESHOLD 0.03 ─────────────────────────────
# More conservative transition definition — 19 transition events
# Tests whether Rho reduction and covariate signs hold with fewer transitions

r_threshold_03 <- run_zi_spatial(df_2021, 0.03, nb_queen, w_queen,
                                 "Check 2 — Threshold 0.03")

# ── CHECKS 3-4: SPATIAL STRUCTURE — DISTANCE WEIGHTS ─────────────────────────
# Distance weights built on active subsample geometry
# use_prebuilt_active = TRUE — weights passed directly, not subsetted

r_dist60 <- run_zi_spatial(df_2021, 0.02, nb_dist60, w_dist60,
                           "Check 3 — Distance 60km (active)",
                           use_prebuilt_active = TRUE)

r_dist80 <- run_zi_spatial(df_2021, 0.02, nb_dist80, w_dist80,
                           "Check 4 — Distance 80km (active)",
                           use_prebuilt_active = TRUE)

# ── CHECKS 5-6: IDENTIFICATION — LISA THRESHOLD ───────────────────────────────
# Tests sensitivity of structural zero classification to LISA p-value cutoff
# p < 0.10: looser — more counties classified as structural zeros
# p < 0.01: stricter — fewer counties classified as structural zeros
# If overlap > 0 at any threshold, perfect separation breaks — flag it

r_lisa_10 <- run_zi_spatial(df_2021, 0.02, nb_queen, w_queen,
                            "Check 5 — LISA p<0.10",
                            lisa_p_threshold = 0.10)

r_lisa_01 <- run_zi_spatial(df_2021, 0.02, nb_queen, w_queen,
                            "Check 6 — LISA p<0.01",
                            lisa_p_threshold = 0.01)

# ── CHECKS 7-9: COVARIATE SPECIFICATION (PART 2 ONLY) ────────────────────────
# Structural zero classification fixed at baseline (LISA p < 0.05)
# Only Part 2 formula changes

# Check 7: Drop net_income_z — marginal in baseline (p = 0.094)
r_drop_netincome <- run_zi_spatial(
  df_2021, 0.02, nb_queen, w_queen,
  "Check 7 — Drop net_income_z",
  part2_formula = transition ~ land_value_z + n_farms_z + cropland_z + acres_z
)

# Check 8: Drop acres_z — not significant in baseline (p = 0.411)
r_drop_acres <- run_zi_spatial(
  df_2021, 0.02, nb_queen, w_queen,
  "Check 8 — Drop acres_z",
  part2_formula = transition ~ land_value_z + net_income_z +
    n_farms_z + cropland_z
)

# Check 9: Add pasture_z — theoretically motivated (pasture-to-cropland
# conversion is a known transition pathway in Missouri)
# pasture_z standardized within active subsample in helper function
r_add_pasture <- run_zi_spatial(
  df_2021, 0.02, nb_queen, w_queen,
  "Check 9 — Add pasture_z",
  part2_formula = transition ~ land_value_z + net_income_z + n_farms_z +
    cropland_z + acres_z + pasture_z
)

# ── COMBINE ───────────────────────────────────────────────────────────────────

robustness_table <- bind_rows(
  r_baseline,
  r_no_stegen,
  r_threshold_03,
  r_dist60,
  r_dist80,
  r_lisa_10,
  r_lisa_01,
  r_drop_netincome,
  r_drop_acres,
  r_add_pasture
)

print(robustness_table)

# ── SAVE RESULTS ──────────────────────────────────────────────────────────────

sink("docs/robustness_results.txt")

cat("=== ROBUSTNESS CHECKS — ZERO-INFLATED SPATIAL MODEL ===\n")
cat("Date: June 2026\n")
cat("Baseline: 2021 cross-section, threshold 0.02, queen contiguity, LISA p<0.05\n")
cat("Baseline: Rho = 0.551, land_value_z p = 0.0002, Moran's I p = 0.293\n")
cat("Economic covariates: 2012 Census of Agriculture (beginning-of-period)\n\n")

fmt_hdr <- "%-42s %5s %5s %5s %6s %8s %8s %8s %6s\n"
fmt_row <- "%-42s %5d %5d %5d %6.3f %8.4f %8.4f %8.4f %6d\n"
div      <- strrep("-", 100)

print_rows <- function(indices) {
  cat(sprintf(fmt_hdr, "Specification", "N", "Str0", "Trans",
              "Rho", "LandVal", "LandVal_p", "MoranI_p", "Ovlap"))
  cat(div, "\n")
  for (i in indices) {
    r <- robustness_table[i, ]
    cat(sprintf(fmt_row, r$Specification, r$N_total, r$N_structural,
                r$N_transition, r$Rho,
                r$LandValue_est, r$LandValue_p, r$MoranI_p, r$Overlap))
  }
  cat("\n")
}

cat("=== CHECK 1: SAMPLE SENSITIVITY ===\n\n")
print_rows(c(1, 2))

cat("=== CHECK 2: OUTCOME DEFINITION ===\n\n")
print_rows(c(1, 3))

cat("=== CHECKS 3-4: SPATIAL STRUCTURE ===\n\n")
cat("Note: distance weights on active subsample (n=91) geometry\n")
cat("Full connectivity not achieved — Ste. Genevieve permanently isolated\n\n")
print_rows(c(1, 4, 5))

cat("=== CHECKS 5-6: IDENTIFICATION ASSUMPTION — LISA THRESHOLD ===\n\n")
cat("Flag: Overlap > 0 means perfect separation fails at that threshold\n\n")
print_rows(c(1, 6, 7))

cat("=== CHECKS 7-9: COVARIATE SPECIFICATION (PART 2 ONLY) ===\n\n")
cat("Structural zero classification fixed at baseline (LISA p < 0.05)\n\n")
print_rows(c(1, 8, 9, 10))

cat("=== ROBUSTNESS CRITERIA ===\n\n")
cat("Findings are robust if across all specifications:\n")
cat("  (1) Rho remains below standard spatial lag baseline (0.645)\n")
cat("  (2) Land value coefficient sign consistently negative\n")
cat("  (3) Moran's I on residuals insignificant (p > 0.05)\n")
cat("  (4) Overlap = 0 (perfect separation holds)\n\n")
cat("Flag any specification where Overlap > 0 or MoranI_p < 0.05\n")

sink()
# ── VUONG TEST: BASELINE vs PASTURE MODEL ────────────────────────────────────
# Compares two non-nested models on the same outcome and sample
# Positive z-score favors model 1, negative favors model 2
# p < 0.05 indicates one model significantly better than the other
#
# Note: Vuong test is designed for non-nested models estimated by MLE
# Both spatial lag models use MLE — test is appropriate here
# Null hypothesis: both models are equally close to the true DGP
# ── RELOAD DEPENDENCIES ───────────────────────────────────────────────────────

library(sf)
library(spdep)
library(spatialreg)
library(dplyr)
library(pscl)

# Load data
counties_mo_proj <- st_read("data/processed/counties_mo.shp") |>
  arrange(GEOID)

panel   <- read.csv("data/processed/panel_final.csv")
w_queen <- readRDS("data/processed/weights_queen.rds")
nb_queen <- readRDS("data/processed/nb_queen.rds")

# Reconstruct analytical sample
df <- panel |>
  filter(year == 2021, !is.na(n_farms)) |>
  arrange(GEOID)

# Subset weights to 114 counties
geoids_114 <- df$GEOID
w_114 <- subset(w_queen,
                subset = counties_mo_proj$GEOID %in% geoids_114)

# Reconstruct structural zero classification
lisa_cropland <- localmoran(df$cropland, w_114)
mean_crop     <- mean(df$cropland)
lag_crop      <- lag.listw(w_114, df$cropland)

df <- df |>
  mutate(
    structural_zero = as.integer(
      cropland < mean_crop &
        lag_crop < mean_crop &
        lisa_cropland[, 5] < 0.05
    ),
    forest_total = forest + frac_42 + frac_43
  )
# ── ADD PASTURE_2011 FROM NLCD PANEL ─────────────────────────────────────────
# panel_final.csv contains pasture for all years
# Extract 2011 pasture and merge to 2021 cross-section

lc_panel <- read.csv("data/processed/nlcd_panel_2001_2011_2021.csv")

pasture_2011 <- lc_panel |>
  filter(year == 2011) |>
  select(GEOID, pasture_2011 = pasture)

df <- df |>
  left_join(pasture_2011, by = "GEOID")

# Confirm no missing values
sum(is.na(df$pasture_2011))

# Check it looks right
summary(df$pasture_2011)
# Reconstruct active subsample
df_active <- df |>
  filter(structural_zero == 0) |>
  mutate(
    land_value_z = as.numeric(scale(land_value_acre)),
    net_income_z = as.numeric(scale(net_income)),
    n_farms_z    = as.numeric(scale(n_farms)),
    cropland_z   = as.numeric(scale(cropland)),
    acres_z      = as.numeric(scale(acres_operated)),
    pasture_z    = as.numeric(scale(pasture_2011))
  )

# Rebuild active weights
nb_active <- subset(nb_queen,
                    subset = counties_mo_proj$GEOID %in% 
                      df$GEOID[df$structural_zero == 0])
w_active  <- nb2listw(nb_active, style = "W", zero.policy = TRUE)

cat("Reload complete\n")
cat("Active subsample:", nrow(df_active), "counties\n")
cat("Structural zeros:", sum(df$structural_zero), "\n")
cat("Transitions:", sum(df_active$transition), "\n")

# ── FIT BOTH MODELS ───────────────────────────────────────────────────────────

model_baseline <- lagsarlm(
  transition ~ land_value_z + net_income_z + n_farms_z + cropland_z + acres_z,
  data        = df_active,
  listw       = w_active,
  zero.policy = TRUE
)

model_pasture <- lagsarlm(
  transition ~ land_value_z + net_income_z + n_farms_z + 
    cropland_z + acres_z + pasture_z,
  data        = df_active,
  listw       = w_active,
  zero.policy = TRUE
)

# ── AIC AND LOG-LIKELIHOOD COMPARISON ────────────────────────────────────────

cat("=== MODEL COMPARISON ===\n\n")
cat("Baseline — AIC:", round(AIC(model_baseline), 3),
    " Log-lik:", round(model_baseline$LL, 3),
    " Rho:", round(model_baseline$rho, 3), "\n")
cat("Pasture  — AIC:", round(AIC(model_pasture), 3),
    " Log-lik:", round(model_pasture$LL, 3),
    " Rho:", round(model_pasture$rho, 3), "\n\n")

# ── LIKELIHOOD RATIO TEST ─────────────────────────────────────────────────────
# Baseline nested within pasture model — LR test valid

lr_stat <- -2 * (model_baseline$LL - model_pasture$LL)
lr_p    <- pchisq(lr_stat, df = 1, lower.tail = FALSE)

cat("=== LIKELIHOOD RATIO TEST ===\n")
cat("LR statistic:", round(lr_stat, 3), "\n")
cat("p-value:     ", round(lr_p, 6), "\n")
cat("AIC reduction:", round(AIC(model_baseline) - AIC(model_pasture), 3), "\n\n")

# ── COEFFICIENT COMPARISON ────────────────────────────────────────────────────

cat("=== COEFFICIENT COMPARISON ===\n\n")
cat(sprintf("%-16s %10s %10s %10s %10s\n",
            "Variable", "Base coef", "Base p", "Past coef", "Past p"))
cat(strrep("-", 58), "\n")

base_sum <- summary(model_baseline)$Coef
past_sum <- summary(model_pasture)$Coef

vars <- c("(Intercept)", "land_value_z", "net_income_z", 
          "n_farms_z", "cropland_z", "acres_z")

for (v in vars) {
  cat(sprintf("%-16s %10.4f %10.4f %10.4f %10.4f\n",
              v,
              base_sum[v, 1], base_sum[v, 4],
              past_sum[v, 1], past_sum[v, 4]))
}

cat(sprintf("%-16s %10s %10s %10.4f %10.4f\n",
            "pasture_z", "—", "—",
            past_sum["pasture_z", 1], past_sum["pasture_z", 4]))

cat(sprintf("%-16s %10.4f %10s %10.4f %10s\n",
            "Rho",
            model_baseline$rho, "",
            model_pasture$rho, ""))

moran_pasture_resid <- moran.test(residuals(model_pasture), 
                                  w_active, 
                                  zero.policy = TRUE)
print(moran_pasture_resid)

message("Script 11 complete - robustness results saved to docs/robustness_results.txt")

library(readxl)
library(dplyr)

# Skip first 3 rows, use row 3 as header
crp_county <- read_excel("data/raw/CRPHistoryCounty86-25.xlsx", 
                         skip = 3)

# Check column names now
names(crp_county)
head(crp_county)
# Filter to Missouri, extract FIPS and 2011 enrollment
crp_2011_mo <- crp_county |>
  filter(STATE == "MISSOURI") |>
  select(COUNTY, FIPS, crp_acres_2011 = `2011`)

# Check
nrow(crp_2011_mo)   # should be close to 115 (Missouri counties)
head(crp_2011_mo)
summary(crp_2011_mo$crp_acres_2011)

# Check for missing values
sum(is.na(crp_2011_mo$crp_acres_2011))

# ── RELOAD DEPENDENCIES ───────────────────────────────────────────────────────

library(sf)
library(spdep)
library(spatialreg)
library(dplyr)
library(readxl)

# Load data
counties_mo_proj <- st_read("data/processed/counties_mo.shp") |>
  arrange(GEOID)

panel    <- read.csv("data/processed/panel_final.csv")
w_queen  <- readRDS("data/processed/weights_queen.rds")
nb_queen <- readRDS("data/processed/nb_queen.rds")

# Reconstruct analytical sample
df <- panel |>
  filter(year == 2021, !is.na(n_farms)) |>
  arrange(GEOID)

# Subset weights to 114 counties
geoids_114 <- df$GEOID
w_114 <- subset(w_queen,
                subset = counties_mo_proj$GEOID %in% geoids_114)

# Reconstruct structural zero classification
lisa_cropland <- localmoran(df$cropland, w_114)
mean_crop     <- mean(df$cropland)
lag_crop      <- lag.listw(w_114, df$cropland)

df <- df |>
  mutate(
    structural_zero = as.integer(
      cropland < mean_crop &
        lag_crop < mean_crop &
        lisa_cropland[, 5] < 0.05
    ),
    forest_total = forest + frac_42 + frac_43
  )

# Add pasture_2011 from NLCD panel
lc_panel <- read.csv("data/processed/nlcd_panel_2001_2011_2021.csv")
pasture_2011 <- lc_panel |>
  filter(year == 2011) |>
  select(GEOID, pasture_2011 = pasture)

df <- df |>
  left_join(pasture_2011, by = "GEOID")

# ── ADD CRP DATA ──────────────────────────────────────────────────────────────

crp_county <- read_excel("data/raw/CRPHistoryCounty86-25.xlsx", skip = 3)

crp_2011_mo <- crp_county |>
  filter(STATE == "MISSOURI") |>
  select(GEOID = FIPS, crp_acres_2011 = `2011`)

df <- df |>
  left_join(crp_2011_mo, by = "GEOID")

cat("Reload complete\n")
cat("N counties:", nrow(df), "\n")
cat("Missing CRP:", sum(is.na(df$crp_acres_2011)), "\n")
cat("Missing pasture_2011:", sum(is.na(df$pasture_2011)), "\n")

# Check land area variable for converting CRP acres to proportion
names(df)[grepl("ALAND|area", names(df), ignore.case = TRUE)]

# Get land area from shapefile (in square meters, UTM projection)
# Convert to acres: 1 sq meter = 0.000247105 acres

area_lookup <- counties_mo_proj |>
  st_drop_geometry() |>
  select(GEOID, ALAND) |>
  mutate(GEOID = as.integer(GEOID),
         land_acres = ALAND * 0.000247105)

df <- df |>
  left_join(area_lookup, by = "GEOID")

# Confirm merge
sum(is.na(df$land_acres))
summary(df$land_acres)

# Create CRP as proportion of county land area
df <- df |>
  mutate(crp_share_2011 = crp_acres_2011 / land_acres)

summary(df$crp_share_2011)

# Compare CRP share to pasture share
cor(df$crp_share_2011, df$pasture_2011)

# Rebuild active subsample with CRP added
df_active <- df |>
  filter(structural_zero == 0) |>
  mutate(
    land_value_z = as.numeric(scale(land_value_acre)),
    net_income_z = as.numeric(scale(net_income)),
    n_farms_z    = as.numeric(scale(n_farms)),
    cropland_z   = as.numeric(scale(cropland)),
    acres_z      = as.numeric(scale(acres_operated)),
    pasture_z    = as.numeric(scale(pasture_2011)),
    crp_z        = as.numeric(scale(crp_share_2011))
  )

# Rebuild active weights
nb_active <- subset(nb_queen,
                    subset = counties_mo_proj$GEOID %in% 
                      df$GEOID[df$structural_zero == 0])
w_active <- nb2listw(nb_active, style = "W", zero.policy = TRUE)

# Check correlation between pasture and CRP in active subsample specifically
cor(df_active$pasture_z, df_active$crp_z)

# ── MODEL WITH CRP ADDED ──────────────────────────────────────────────────────

model_crp <- lagsarlm(
  transition ~ land_value_z + net_income_z + n_farms_z +
    cropland_z + acres_z + pasture_z + crp_z,
  data        = df_active,
  listw       = w_active,
  zero.policy = TRUE
)

summary(model_crp)

# Compare to pasture-only model
model_pasture <- lagsarlm(
  transition ~ land_value_z + net_income_z + n_farms_z +
    cropland_z + acres_z + pasture_z,
  data        = df_active,
  listw       = w_active,
  zero.policy = TRUE
)

cat("\nAIC pasture only:", round(AIC(model_pasture), 3), "\n")
cat("AIC pasture + CRP:", round(AIC(model_crp), 3), "\n")

lr_stat <- -2 * (model_pasture$LL - model_crp$LL)
lr_p <- pchisq(lr_stat, df = 1, lower.tail = FALSE)
cat("LR test (CRP addition): stat =", round(lr_stat, 3), 
    " p =", round(lr_p, 6), "\n")

moran_crp_resid <- moran.test(residuals(model_crp), w_active, zero.policy = TRUE)
print(moran_crp_resid)