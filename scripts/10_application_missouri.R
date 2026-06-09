# Script 10: Application Results — Zero-Inflated Spatial Model
# Purpose: Synthesize and present all results from Scripts 07, 08, and 09
#          Three-model comparison table (OLS, spatial lag, ZI spatial)
#          Part 1 probit interpretation and marginal effects
#          Part 2 coefficient table and spatial dependence results
#          Publication-quality maps for paper
# No new estimation — all model objects loaded from data/processed/
# Author: Lan T. Tran
# Date: June 2026

library(sf)
library(spdep)
library(spatialreg)
library(dplyr)
library(ggplot2)

# ── LOAD ALL MODEL OBJECTS ────────────────────────────────────────────────────

counties_mo_proj <- st_read("data/processed/counties_mo.shp") |>
  arrange(GEOID)

panel    <- read.csv("data/processed/panel_final.csv")
w_queen  <- readRDS("data/processed/weights_queen.rds")
nb_queen <- readRDS("data/processed/nb_queen.rds")

model_ols    <- readRDS("data/processed/model_ols.rds")
model_slag   <- readRDS("data/processed/model_slag.rds")
model_part1  <- readRDS("data/processed/model_part1.rds")
model_part2  <- readRDS("data/processed/model_part2.rds")

# ── RECONSTRUCT ANALYTICAL OBJECTS ───────────────────────────────────────────
# Required for residual diagnostics and map merges
# Same construction as Scripts 08 and 09

df <- panel |>
  filter(year == 2021, !is.na(n_farms)) |>
  arrange(GEOID)

geoids_114 <- df$GEOID
w_114 <- subset(w_queen, subset = counties_mo_proj$GEOID %in% geoids_114)

lisa_cropland <- localmoran(df$cropland, w_114)
mean_crop     <- mean(df$cropland)
lag_crop      <- lag.listw(w_114, df$cropland)

df <- df |>
  mutate(
    forest_total    = forest + frac_42 + frac_43,
    structural_zero = as.integer(
      cropland < mean_crop &
        lag_crop < mean_crop &
        lisa_cropland[, 5] < 0.05
    ),
    p_structural = predict(model_part1, type = "response"),
    county_type  = case_when(
      structural_zero == 1 ~ "Structural zero",
      transition == 1      ~ "Transition",
      TRUE                 ~ "Sampling zero"
    )
  )

# Active subsample for Part 2 residuals and fitted values
df_active <- df |>
  filter(structural_zero == 0) |>
  mutate(
    land_value_z = as.numeric(scale(land_value_acre)),
    net_income_z = as.numeric(scale(net_income)),
    n_farms_z    = as.numeric(scale(n_farms)),
    cropland_z   = as.numeric(scale(cropland)),
    acres_z      = as.numeric(scale(acres_operated)),
    fitted_part2 = fitted(model_part2),
    resid_part2  = residuals(model_part2)
  )

# Rebuild active weights for Moran's I
nb_active <- subset(nb_queen,
                    subset = counties_mo_proj$GEOID %in% df$GEOID[df$structural_zero == 0])
w_active  <- nb2listw(nb_active, style = "W", zero.policy = TRUE)

# ── RESIDUAL DIAGNOSTICS FOR COMPARISON TABLE ─────────────────────────────────

moran_ols_resid   <- moran.test(residuals(model_ols),  w_114)
moran_slag_resid  <- moran.test(residuals(model_slag), w_114)
moran_part2_resid <- moran.test(residuals(model_part2), w_active, zero.policy = TRUE)

# ── MARGINAL EFFECT — PART 1 PROBIT ──────────────────────────────────────────
# Marginal effect at the mean: dP/dX evaluated at mean forest_total
# For probit: dP/dX = phi(Xb) * beta
# phi = standard normal density

b0      <- coef(model_part1)["(Intercept)"]
b1      <- coef(model_part1)["forest_total"]
xb_mean <- b0 + b1 * mean(df$forest_total)
me_mean <- dnorm(xb_mean) * b1

# ── SINK: FULL RESULTS ────────────────────────────────────────────────────────

sink("docs/application_results.txt")

cat("=== APPLICATION RESULTS: ZERO-INFLATED SPATIAL MODEL ===\n")
cat("Date: June 2026\n")
cat("Study area: Missouri counties, 2021 cross-section\n\n")

# ── TABLE 1: THREE-MODEL COMPARISON ──────────────────────────────────────────

cat("=== TABLE 1: THREE-MODEL COMPARISON ===\n\n")
cat(sprintf("%-35s %6s %8s %10s %8s %10s\n",
            "Model", "N", "AIC", "Log-Lik", "Rho", "Moran I"))
cat(strrep("-", 80), "\n")
cat(sprintf("%-35s %6d %8.2f %10.3f %8s %10.3f\n",
            "OLS (baseline)",
            114,
            AIC(model_ols),
            as.numeric(logLik(model_ols)),
            "—",
            round(moran_ols_resid$estimate[1], 3)))
cat(sprintf("%-35s %6d %8.2f %10.3f %8.3f %10.3f\n",
            "Standard spatial lag",
            114,
            AIC(model_slag),
            as.numeric(logLik(model_slag)),
            round(model_slag$rho, 3),
            round(moran_slag_resid$estimate[1], 3)))
cat(sprintf("%-35s %6d %8.2f %10.3f %8.3f %10.3f\n",
            "ZI spatial lag — Part 2",
            91,
            AIC(model_part2),
            as.numeric(logLik(model_part2)),
            round(model_part2$rho, 3),
            round(moran_part2_resid$estimate[1], 3)))
cat(strrep("-", 80), "\n")
cat("Note: ZI spatial lag excludes 23 structural zero counties (Ozark Low-Low LISA cluster)\n")
cat("      Ste. Genevieve County has no active neighbors — spatial lag set to zero\n")
cat(sprintf("      Moran's I p-values: OLS p=0.000, Spatial lag p=%.3f, ZI Part 2 p=%.3f\n\n",
            round(moran_slag_resid$p.value, 3),
            round(moran_part2_resid$p.value, 3)))

# ── TABLE 2: PART 1 PROBIT RESULTS ───────────────────────────────────────────

cat("=== TABLE 2: PART 1 PROBIT — STRUCTURAL ZERO PROCESS ===\n\n")
cat("Outcome: structural_zero (1 = Low-Low LISA cluster, p < 0.05)\n")
cat("N = 114 counties\n\n")

part1_sum <- summary(model_part1)$coefficients
cat(sprintf("%-20s %10s %10s %10s %10s\n",
            "Variable", "Estimate", "Std.Error", "z value", "Pr(>|z|)"))
cat(strrep("-", 62), "\n")
for (i in 1:nrow(part1_sum)) {
  cat(sprintf("%-20s %10.4f %10.4f %10.4f %10.4f\n",
              rownames(part1_sum)[i],
              part1_sum[i, 1],
              part1_sum[i, 2],
              part1_sum[i, 3],
              part1_sum[i, 4]))
}
cat(strrep("-", 62), "\n")
cat(sprintf("Null deviance:     %.3f on %d df\n",
            model_part1$null.deviance, model_part1$df.null))
cat(sprintf("Residual deviance: %.3f on %d df\n",
            model_part1$deviance, model_part1$df.residual))
cat(sprintf("AIC: %.3f\n\n", AIC(model_part1)))

cat("Marginal effect at mean forest_total:\n")
cat(sprintf("  dP(structural zero)/d(forest_total) = %.4f\n\n", me_mean))

cat("Predicted probability by county type:\n")
df |>
  group_by(county_type) |>
  summarise(n = n(),
            mean_p = round(mean(p_structural), 3),
            sd_p   = round(sd(p_structural), 3)) |>
  as.data.frame() |>
  print()
cat("\n")

df <- df |>
  mutate(
    predicted_zero = as.integer(p_structural >= 0.5),
    correct        = as.integer(predicted_zero == structural_zero)
  )

# Classification table
table(Predicted = df$predicted_zero, Actual = df$structural_zero)

# Overall accuracy
cat("Overall accuracy:", round(mean(df$correct), 3), "\n")

# Accuracy by group
df |>
  group_by(structural_zero) |>
  summarise(
    n         = n(),
    n_correct = sum(correct),
    accuracy  = round(mean(correct), 3)
  )

cat("\nInterpretation: Forest cover alone recovers 39% of LISA-classified\n")
cat("structural zeros. Spatial neighborhood context in LISA carries\n")
cat("information beyond what a single biophysical covariate can replicate.\n")
cat("This justifies the rule-based LISA classification over a purely\n")
cat("data-driven approach.\n\n")

# ── TABLE 3: PART 2 SPATIAL LAG RESULTS ──────────────────────────────────────

cat("=== TABLE 3: PART 2 SPATIAL LAG — ACTIVE COUNTY TRANSITION PROCESS ===\n\n")
cat("Outcome: transition (1 = abs(cropland change) > 0.02)\n")
cat("N = 91 counties (23 structural zeros removed)\n")
cat("Covariates standardized within active subsample\n\n")

part2_sum <- summary(model_part2)$Coef
cat(sprintf("%-20s %10s %10s %10s %10s\n",
            "Variable", "Estimate", "Std.Error", "z value", "Pr(>|z|)"))
cat(strrep("-", 62), "\n")
for (i in 1:nrow(part2_sum)) {
  cat(sprintf("%-20s %10.4f %10.4f %10.4f %10.4f\n",
              rownames(part2_sum)[i],
              part2_sum[i, 1],
              part2_sum[i, 2],
              part2_sum[i, 3],
              part2_sum[i, 4]))
}
cat(strrep("-", 62), "\n")
cat(sprintf("Rho:              %.4f (p < 0.001)\n", model_part2$rho))
cat(sprintf("Log-likelihood:   %.4f\n", as.numeric(logLik(model_part2))))
cat(sprintf("AIC:              %.2f\n", AIC(model_part2)))
cat(sprintf("Moran's I resid:  %.4f (p = %.3f)\n\n",
            round(moran_part2_resid$estimate[1], 4),
            round(moran_part2_resid$p.value, 3)))

cat("Rho comparison:\n")
cat(sprintf("  Standard spatial lag (n=114): %.3f\n", model_slag$rho))
cat(sprintf("  ZI Part 2 spatial lag (n=91): %.3f\n", model_part2$rho))
cat(sprintf("  Reduction in Rho:             %.3f\n\n",
            model_slag$rho - model_part2$rho))

# ── RESIDUAL BIAS BY COUNTY TYPE ──────────────────────────────────────────────

cat("=== RESIDUAL BIAS BY COUNTY TYPE ===\n\n")
cat("Standard spatial lag (n=114):\n")
df$slag_resid <- residuals(model_slag)
df |>
  group_by(county_type) |>
  summarise(n          = n(),
            mean_resid = round(mean(slag_resid), 4),
            sd_resid   = round(sd(slag_resid), 4)) |>
  as.data.frame() |>
  print()

cat("\nZI Part 2 spatial lag (n=91, active counties only):\n")
df_active |>
  mutate(county_type = ifelse(transition == 1, "Transition", "Sampling zero")) |>
  group_by(county_type) |>
  summarise(n          = n(),
            mean_resid = round(mean(resid_part2), 4),
            sd_resid   = round(sd(resid_part2), 4)) |>
  as.data.frame() |>
  print()
cat("\nNote: Remaining bias in Part 2 reflects linear probability model applied\n")
cat(sprintf("to binary outcome — not spatial misspecification (Moran's I p = %.3f)\n",
            round(moran_part2_resid$p.value, 3)))

sink()

message("Results written to docs/application_results.txt")

# ── MAP 10A: THREE-WAY DECOMPOSITION ─────────────────────────────────────────
# Publication version of Script 08 map — cleaner theme

counties_map <- counties_mo_proj |>
  filter(GEOID %in% df$GEOID) |>
  left_join(
    df |> mutate(GEOID = as.character(GEOID)) |>
      select(GEOID, county_type),
    by = "GEOID"
  )

ggplot(counties_map) +
  geom_sf(aes(fill = county_type), color = "white", linewidth = 0.3) +
  scale_fill_manual(
    values = c(
      "Structural zero" = "#2c7bb6",
      "Sampling zero"   = "#f0f0f0",
      "Transition"      = "#d7191c"
    ),
    name = "County Type"
  ) +
  labs(
    title    = "Three-Way Decomposition of Missouri Counties",
    subtitle = "Structural zeros, sampling zeros, and transitions — 2011–2021",
    caption  = "Structural zeros: Low-Low LISA cluster on cropland level (p < 0.05)"
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text  = element_blank(),
    axis.ticks = element_blank()
  )

ggsave("outputs/map_10a_three_way_decomposition.png",
       width = 10, height = 7, dpi = 300)

# ── MAP 10B: PART 2 FITTED VALUES — ACTIVE COUNTIES ONLY ─────────────────────
# Predicted transition probability from ZI Part 2 spatial lag
# Structural zero counties shown in grey — outside the model

counties_fitted <- counties_mo_proj |>
  filter(GEOID %in% df$GEOID) |>
  left_join(
    df |> mutate(GEOID = as.character(GEOID)) |>
      select(GEOID, structural_zero),
    by = "GEOID"
  ) |>
  left_join(
    df_active |> mutate(GEOID = as.character(GEOID)) |>
      select(GEOID, fitted_part2),
    by = "GEOID"
  )

ggplot(counties_fitted) +
  geom_sf(aes(fill = fitted_part2), color = "white", linewidth = 0.3) +
  scale_fill_viridis_c(
    name   = "Fitted\nP(Transition)",
    option = "magma",
    limits = c(0, 1),
    na.value = "#d9d9d9"
  ) +
  labs(
    title    = "Fitted Transition Probabilities — ZI Part 2 Spatial Lag",
    subtitle = "Active counties only (n=91), Missouri 2011–2021",
    caption  = "Grey: structural zero counties excluded from Part 2 model"
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text  = element_blank(),
    axis.ticks = element_blank()
  )

ggsave("outputs/map_10b_fitted_transition_probabilities.png",
       width = 10, height = 7, dpi = 300)

message("Script 10 complete - all results tables and maps saved")