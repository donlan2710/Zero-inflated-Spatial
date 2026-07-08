# Script 08: Zero Structure Diagnosis
# Purpose: Motivate the two-part model. Show two things. First, the zeros are
#          spatially clustered, so a spatial model is needed. Second, a
#          standard spatial model cannot separate structural zeros from
#          sampling zeros, so it mispredicts in a patterned way, which is why
#          a two-part model is needed.
#
# This script diagnoses the LINEAR spatial baseline from Script 07, the model
# a standard researcher would run first. It reads model_baseline_linear_114.rds.
# If Script 07 is rerun, rerun this script, because the diagnosis reflects
# whatever model is saved in that file.
#
# Note: this runs on the two-variable baseline (cropland_lag + pasture_2011).
# Any residual numbers in older documents came from the old five-variable
# model and do not apply here. Every number below is printed fresh.
#
# Author: Lan T. Tran
# Date: July 2026

library(sf)
library(spdep)
library(dplyr)
library(ggplot2)

# ── LOAD DATA ─────────────────────────────────────────────────────────────────

counties_mo_proj <- st_read("data/processed/counties_mo.shp") |>
  arrange(GEOID)

panel   <- read.csv("data/processed/panel_final.csv",
                    colClasses = c(GEOID = "character"))
pasture <- read.csv("data/processed/pasture_by_year.csv",
                    colClasses = c(GEOID = "character"))
w_queen <- readRDS("data/processed/weights_queen.rds")
slag    <- readRDS("data/processed/model_baseline_linear_114.rds")

# ── ANALYTICAL SAMPLE, 114 COUNTIES ──────────────────────────────────────────
# Same construction as Script 07, so the residuals line up with the model.

df <- panel |>
  filter(year == 2021, !is.na(n_farms)) |>
  arrange(GEOID) |>
  left_join(pasture |> select(GEOID, pasture_2011), by = "GEOID") |>
  mutate(
    cropland_z = as.numeric(scale(cropland_lag)),
    pasture_z  = as.numeric(scale(pasture_2011))
  )

n_sample <- nrow(df)
message("Analytical sample size: ", n_sample, " counties")
if (n_sample != 114) {
  warning("Sample size is not 114. Check the year and n_farms filter.")
}

geoids_114 <- df$GEOID
w_114 <- subset(w_queen, subset = counties_mo_proj$GEOID %in% geoids_114)

# ── 1. ARE THE ZEROS SPATIALLY CLUSTERED ─────────────────────────────────────
# If zeros are scattered at random, a spatial model is not needed. If they
# cluster, it is. Moran's I on the zero indicator answers this.

cat("=== ZERO STRUCTURE ANALYSIS ===\n\n")

cat("Transition distribution 2011 to 2021:\n")
print(table(df$transition))

zero_rate  <- round(100 * mean(df$transition == 0), 1)
trans_rate <- round(100 * mean(df$transition == 1), 1)
cat("\nZero rate:", zero_rate, "percent\n")
cat("Transition rate:", trans_rate, "percent\n")

zero_indicator <- as.integer(df$transition == 0)
moran_zeros <- moran.test(zero_indicator, w_114, zero.policy = TRUE)
cat("\nMoran's I on the zero indicator:\n")
cat("  Moran's I:", round(moran_zeros$estimate[1], 3),
    "  p =", format.pval(moran_zeros$p.value, digits = 3), "\n")
cat("  A positive, significant value means zeros cluster in space.\n")

# ── 2. SEPARATE STRUCTURAL FROM SAMPLING ZEROS ───────────────────────────────
# Structural zero: Low-Low LISA cluster on cropland level, p below 0.05.
# Same classification as Script 07 and Part 2.

lisa_cropland <- localmoran(df$cropland, w_114)
mean_crop     <- mean(df$cropland)
lag_crop      <- lag.listw(w_114, df$cropland)

df <- df |>
  mutate(
    lisa_p          = lisa_cropland[, 5],
    lag_cropland    = lag_crop,
    structural_zero = as.integer(
      cropland < mean_crop & lag_cropland < mean_crop & lisa_p < 0.05
    )
  )

n_structural <- sum(df$structural_zero)
n_sampling   <- sum(df$transition == 0) - n_structural
n_transition <- sum(df$transition == 1)

cat("\nThree way decomposition:\n")
cat("  Structural zeros:", n_structural, "\n")
cat("  Sampling zeros:  ", n_sampling, "\n")
cat("  Transitions:     ", n_transition, "\n")

cat("\nCross tabulation, structural zero by transition:\n")
cross_tab <- table(Structural = df$structural_zero, Transition = df$transition)
print(cross_tab)

if (cross_tab["1", "1"] != 0) {
  warning("A structural zero county shows a transition. Perfect separation ",
          "no longer holds. Check the classification.")
} else {
  cat("\nPerfect separation holds: no structural zero county transitioned.\n")
  cat("This is the sharpest evidence that two processes exist in the data.\n")
}

# ── 3. STANDARD MODEL FAILURE, RESIDUALS BY COUNTY TYPE ───────────────────────
# The linear baseline has no way to separate the two zero types. If it
# mispredicts in a patterned way across county types, that patterned failure
# is what motivates the two-part model.

cat("\n=== STANDARD MODEL FAILURE DIAGNOSTICS ===\n")

df$slag_resid <- residuals(slag)

cat("\nResidual summary, linear baseline on 114:\n")
print(summary(df$slag_resid))

df <- df |>
  mutate(county_type = case_when(
    structural_zero == 1 ~ "Structural zero",
    transition == 0      ~ "Sampling zero",
    transition == 1      ~ "Transition"
  ))

resid_summary <- df |>
  group_by(county_type) |>
  summarise(
    n          = n(),
    mean_resid = round(mean(slag_resid), 4),
    sd_resid   = round(sd(slag_resid), 4),
    .groups = "drop"
  )

cat("\nMean residual by county type (printed fresh, not from old documents):\n")
print(as.data.frame(resid_summary))
cat("\nA patterned bias, for example zeros over-predicted and transitions\n")
cat("under-predicted, shows the standard model cannot handle the zero\n")
cat("structure and motivates the two-part model.\n")

# ── 4. MAP THE THREE WAY DECOMPOSITION ───────────────────────────────────────

counties_diag <- counties_mo_proj |>
  filter(GEOID %in% df$GEOID) |>
  left_join(df |> select(GEOID, county_type), by = "GEOID")

ggplot(counties_diag) +
  geom_sf(aes(fill = county_type), color = "white", linewidth = 0.3) +
  scale_fill_manual(
    values = c(
      "Structural zero" = "#2c7bb6",
      "Sampling zero"   = "#f0f0f0",
      "Transition"      = "#d7191c"
    ),
    name = "County type"
  ) +
  labs(
    title    = "Three way decomposition of Missouri counties",
    subtitle = "Structural zeros, sampling zeros, and transitions, 2011 to 2021",
    caption  = "Structural zeros: Low-Low LISA cluster on cropland level, p < 0.05"
  ) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.text  = element_blank(),
        axis.ticks = element_blank())

ggsave("outputs/map_08_zero_structure.png", width = 10, height = 7, dpi = 300)

# ── 5. WRITE DIAGNOSTIC MEMO ──────────────────────────────────────────────────
# Every number below is pulled from an object computed above. Nothing is
# typed as a fixed value. Rerun after any change to Script 07 or the panel.

sink("docs/diagnostic_memo.txt")
cat("=== ZERO STRUCTURE DIAGNOSTIC MEMO ===\n")
cat("Date: July 2026\n")
cat("Model diagnosed: linear spatial baseline, two variables, 114 counties\n\n")

cat("1. ZERO PREVALENCE\n")
cat("Sample size:", n_sample, "counties\n")
cat("Zero rate:", zero_rate, "percent\n")
cat("Structural zeros:", n_structural, "\n")
cat("Sampling zeros:", n_sampling, "\n")
cat("Transitions:", n_transition, "\n\n")

cat("2. SPATIAL CLUSTERING OF ZEROS\n")
cat("Moran's I on zero indicator:", round(moran_zeros$estimate[1], 3),
    " p =", format.pval(moran_zeros$p.value, digits = 3), "\n\n")

cat("3. STRUCTURAL VERSUS SAMPLING ZEROS\n")
print(cross_tab)
cat("Structural zero and transition cell:", cross_tab["1", "1"], "\n\n")

cat("4. STANDARD MODEL FAILURE\n")
print(as.data.frame(resid_summary))
cat("\n")

cat("5. CONCLUSION\n")
cat("The standard spatial model addresses spatial dependence but cannot\n")
cat("separate structural from sampling zeros. A two-part model is needed to\n")
cat("model the structural non-transition process separately from the active\n")
cat("transition process.\n")
sink()

message("Script 08 complete. Zero structure diagnosis written.")
