# Script 09: Zero-Inflated Spatial Model — Option A
# Purpose: Two-stage zero-inflated spatial model separating structural
#          non-transition process from active transition process
#          Part 1: probit for structural zero classification (biophysical)
#          Part 2: spatial lag for transition drivers (economic)
#          Compare Rho and fit against standard spatial lag from Script 07
# Author: Lan T. Tran
# Date: June 2026

library(sf)
library(spdep)
library(spatialreg)
library(dplyr)
library(ggplot2)

# ── LOAD DATA ─────────────────────────────────────────────────────────────────

counties_mo_proj <- st_read("data/processed/counties_mo.shp") |>
  arrange(GEOID)

panel    <- read.csv("data/processed/panel_final.csv")
w_queen  <- readRDS("data/processed/weights_queen.rds")
nb_queen <- readRDS("data/processed/nb_queen.rds")

# ── ANALYTICAL SAMPLE ─────────────────────────────────────────────────────────
# 2021 cross-section, exclude St. Louis City (no agricultural operations)
# Same filter as Scripts 07 and 08

df <- panel |>
  filter(year == 2021, !is.na(n_farms)) |>
  arrange(GEOID)

# ── SUBSET WEIGHTS TO 114-COUNTY SAMPLE ───────────────────────────────────────

geoids_114 <- df$GEOID
w_114 <- subset(w_queen, subset = counties_mo_proj$GEOID %in% geoids_114)

# ── RECONSTRUCT STRUCTURAL ZERO CLASSIFICATION ────────────────────────────────
# Low-Low LISA cluster on cropland level (p < 0.05)
# Same logic as Script 08 — 23 Ozark highland counties
# These counties are structurally incapable of cropland transition
# regardless of economic conditions

lisa_cropland <- localmoran(df$cropland, w_114)
mean_crop     <- mean(df$cropland)
lag_crop      <- lag.listw(w_114, df$cropland)

df <- df |>
  mutate(
    structural_zero = as.integer(
      cropland < mean_crop &
        lag_crop < mean_crop &
        lisa_cropland[, 5] < 0.05
    )
  )

# Three-way decomposition
# Structural zeros (Low-Low LISA): 23 counties — Ozark highlands
# Sampling zeros (non-transition, not structural): 62 counties
# Transitions: 29 counties — northern grain belt and southeast bootheel
# Perfect separation: zero structural zero counties transitioned
table(Structural = df$structural_zero, Transition = df$transition)

# ── COMBINE FOREST COVER CLASSES ──────────────────────────────────────────────
# frac_41 (deciduous), frac_42 (evergreen), frac_43 (mixed) are correlated
# frac_42 and frac_43 correlation = 0.868 — combined into forest_total
# forest is already renamed frac_41 in Script 02

df <- df |>
  mutate(forest_total = forest + frac_42 + frac_43)

# ── PART 1: PROBIT — STRUCTURAL ZERO PROCESS ──────────────────────────────────
# Outcome: structural_zero (binary)
# Predictor: forest_total — total forest cover share
# Biophysical predictor only — structural incapacity is terrain-driven
# Wetland dropped: structural zero counties have less wetland than active
# counties (mean 0.001 vs 0.021) — no theoretical or empirical support
# Three separate forest classes dropped: frac_42/frac_43 correlation = 0.868

model_part1 <- glm(structural_zero ~ forest_total,
                   data   = df,
                   family = binomial(link = "probit"))

summary(model_part1)

# Predicted probabilities by county type
# Confirms model assigns high probability to structural zeros
# and near-zero probability to transitioning counties
df$p_structural <- predict(model_part1, type = "response")

df |>
  mutate(county_type = case_when(
    structural_zero == 1 ~ "Structural zero",
    transition == 1      ~ "Transition",
    TRUE                 ~ "Sampling zero"
  )) |>
  group_by(county_type) |>
  summarise(
    n      = n(),
    mean_p = round(mean(p_structural), 3),
    sd_p   = round(sd(p_structural), 3)
  )

# ── PART 2: SPATIAL LAG — ACTIVE COUNTY TRANSITION PROCESS ───────────────────
# Active subsample: 91 counties (structural zeros removed)
# Covariates re-standardized within active subsample
# as.numeric() strips matrix attribute from scale() for lagsarlm compatibility

df_active <- df |>
  filter(structural_zero == 0) |>
  mutate(
    land_value_z = as.numeric(scale(land_value_acre)),
    net_income_z = as.numeric(scale(net_income)),
    n_farms_z    = as.numeric(scale(n_farms)),
    cropland_z   = as.numeric(scale(cropland)),
    acres_z      = as.numeric(scale(acres_operated))
  ) |>
  select(GEOID, transition, land_value_z, net_income_z,
         n_farms_z, cropland_z, acres_z)

# Rebuild neighbor list for 91-county active subsample
# subset() requires logical vector aligned to full 115-county shapefile
nb_active <- subset(nb_queen,
                    subset = counties_mo_proj$GEOID %in% df$GEOID[df$structural_zero == 0])
w_active  <- nb2listw(nb_active, style = "W", zero.policy = TRUE)

# Note: subsetting causes one island — Ste. Genevieve County (region 77)
# Its only neighbors were structural zero counties in the Ozark highlands
# zero.policy = TRUE assigns spatial lag of zero for this county
# Reported as a limitation in paper

part2_formula <- transition ~ land_value_z + net_income_z +
  n_farms_z + cropland_z + acres_z

model_part2 <- lagsarlm(part2_formula,
                        data        = df_active,
                        listw       = w_active,
                        zero.policy = TRUE)

summary(model_part2)

# ── RESIDUAL DIAGNOSTICS ──────────────────────────────────────────────────────

resid_part2       <- residuals(model_part2)
moran_part2_resid <- moran.test(resid_part2, w_active, zero.policy = TRUE)

cat("\nMoran's I on Part 2 residuals:\n")
print(moran_part2_resid)

df_active$part2_resid <- resid_part2

df_active |>
  mutate(county_type = ifelse(transition == 1, "Transition", "Sampling zero")) |>
  group_by(county_type) |>
  summarise(
    n          = n(),
    mean_resid = round(mean(part2_resid), 4),
    sd_resid   = round(sd(part2_resid), 4)
  )

# ── MODEL COMPARISON ──────────────────────────────────────────────────────────
# Compare Rho and fit against standard spatial lag from Script 07

model_slag_full <- readRDS("data/processed/model_slag.rds")

cat("\nRho comparison:\n")
cat("Standard spatial lag (n=114): ", round(model_slag_full$rho, 3), "\n")
cat("ZI Part 2 spatial lag (n=91): ", round(model_part2$rho, 3), "\n")
cat("Rho reduction:                ", round(model_slag_full$rho - model_part2$rho, 3), "\n")

cat("\nAIC comparison:\n")
cat("Standard spatial lag (n=114): ", round(AIC(model_slag_full), 2), "\n")
cat("ZI Part 2 spatial lag (n=91): ", round(AIC(model_part2), 2), "\n")

# ── SAVE MODELS AND RESULTS ───────────────────────────────────────────────────

saveRDS(model_part1, "data/processed/model_part1.rds")
saveRDS(model_part2, "data/processed/model_part2.rds")

sink("docs/zi_model_results.txt")
cat("=== ZERO-INFLATED SPATIAL MODEL — OPTION A ===\n")
cat("Date: June 2026\n\n")

cat("=== PART 1: PROBIT — STRUCTURAL ZERO PROCESS ===\n")
cat("Predictor: forest_total (frac_41 + frac_42 + frac_43)\n")
cat("N = 114 counties\n\n")
print(summary(model_part1))

cat("\n=== PART 2: SPATIAL LAG — ACTIVE COUNTIES ===\n")
cat("N = 91 counties (23 structural zeros removed)\n")
cat("Island: Ste. Genevieve County — neighbors all structural zeros\n\n")
print(summary(model_part2))

cat("\n=== RHO COMPARISON ===\n")
cat("Standard spatial lag (n=114): ", round(model_slag_full$rho, 3), "\n")
cat("ZI Part 2 spatial lag (n=91): ", round(model_part2$rho, 3), "\n")
cat("Rho reduction:                ", round(model_slag_full$rho - model_part2$rho, 3), "\n")

cat("\n=== RESIDUAL DIAGNOSTICS ===\n")
cat("Moran's I Part 2 residuals:   ", round(moran_part2_resid$estimate[1], 3),
    " p =", round(moran_part2_resid$p.value, 3), "\n")
sink()

# ── MAP: STRUCTURAL ZERO PROBABILITIES ───────────────────────────────────────

df <- df |>
  mutate(county_type = case_when(
    structural_zero == 1 ~ "Structural zero",
    transition == 1      ~ "Transition",
    TRUE                 ~ "Sampling zero"
  ))

counties_map <- counties_mo_proj |>
  filter(GEOID %in% df$GEOID) |>
  left_join(
    df |> mutate(GEOID = as.character(GEOID)) |>
      select(GEOID, p_structural, county_type),
    by = "GEOID"
  )

ggplot(counties_map) +
  geom_sf(aes(fill = p_structural), color = "white", linewidth = 0.3) +
  scale_fill_viridis_c(name = "P(Structural\nzero)",
                       option = "plasma",
                       limits = c(0, 1)) +
  labs(
    title    = "Predicted Probability of Structural Zero — Part 1 Probit",
    subtitle = "Missouri Counties, 2021",
    caption  = "Predictor: total forest cover share (NLCD classes 41, 42, 43)"
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text  = element_blank(),
    axis.ticks = element_blank()
  )

ggsave("outputs/map_09_structural_zero_probabilities.png",
       width = 10, height = 7, dpi = 300)

message("Script 09 complete - Option A ZI spatial model estimated and saved")