# Script 09a: Part 1, the structural zero probit
# Purpose: classify which counties are structurally unable to expand cropland.
#          Predictor is forest cover. Forest reflects the terrain constraint
#          of the Ozark region, which is why it belongs here and pasture does
#          not. Pasture is a Part 2 variable, it explains transition, not
#          structural incapacity.
#
# This script runs on all 114 counties, because Part 1 classifies every county.
# Part 2 (the transition probit) is a separate script and runs on the active
# counties only.
#
# A diagnostic block adds pasture to the Part 1 probit as a ONE TIME CHECK, to
# confirm forest carries the structural signal and pasture does not add to it.
# The check does not change the model. Forest alone stays the specification.
# The decision to keep forest alone is theoretical, not statistical: structural
# incapacity is a terrain matter, and pasture belongs to the transition frame.
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

# ── ANALYTICAL SAMPLE, 114 COUNTIES ──────────────────────────────────────────

df <- panel |>
  filter(year == 2021, !is.na(n_farms)) |>
  arrange(GEOID) |>
  left_join(pasture |> select(GEOID, pasture_2011), by = "GEOID")

n_sample <- nrow(df)
message("Analytical sample size: ", n_sample, " counties")
if (n_sample != 114) {
  warning("Sample size is not 114. Check the year and n_farms filter.")
}

geoids_114 <- df$GEOID
w_114 <- subset(w_queen, subset = counties_mo_proj$GEOID %in% geoids_114)

# ── BUILD FOREST TOTAL, WITH A COLUMN CHECK ──────────────────────────────────
# The old script built forest_total from forest + frac_42 + frac_43. Confirm
# those columns exist before using them, so the script stops clearly rather
# than failing halfway.

need_cols <- c("forest", "frac_42", "frac_43")
have_cols <- need_cols %in% names(df)
if (!all(have_cols)) {
  stop("Missing forest columns: ",
       paste(need_cols[!have_cols], collapse = ", "),
       ". Check the panel before building forest_total.")
}

df <- df |> mutate(forest_total = forest + frac_42 + frac_43)

cat("forest_total built from forest + frac_42 + frac_43.\n")
cat("Summary of forest_total:\n")
print(summary(df$forest_total))

# ── STRUCTURAL ZERO CLASSIFICATION ───────────────────────────────────────────
# Low-Low LISA cluster on cropland level, p below 0.05. Same rule as 07 and 08.

lisa_cropland <- localmoran(df$cropland, w_114)
mean_crop     <- mean(df$cropland)
lag_crop      <- lag.listw(w_114, df$cropland)

df <- df |>
  mutate(
    structural_zero = as.integer(
      cropland < mean_crop & lag_crop < mean_crop & lisa_cropland[, 5] < 0.05
    )
  )

cat("\nStructural zeros:", sum(df$structural_zero), "\n")
cat("Cross tabulation, structural zero by transition:\n")
print(table(Structural = df$structural_zero, Transition = df$transition))

# ── PART 1 MODEL: PROBIT, FOREST ONLY ────────────────────────────────────────

model_part1 <- glm(structural_zero ~ forest_total,
                   data   = df,
                   family = binomial(link = "probit"))

cat("\n=== PART 1 PROBIT, FOREST ONLY (the model) ===\n")
print(summary(model_part1))

# Predicted probability by county type
df$p_structural <- predict(model_part1, type = "response")

df <- df |>
  mutate(county_type = case_when(
    structural_zero == 1 ~ "Structural zero",
    transition == 1      ~ "Transition",
    TRUE                 ~ "Sampling zero"
  ))

cat("\nPredicted P(structural zero) by county type:\n")
df |>
  group_by(county_type) |>
  summarise(n = n(),
            mean_p = round(mean(p_structural), 3),
            sd_p   = round(sd(p_structural), 3),
            .groups = "drop") |>
  as.data.frame() |>
  print()

# Marginal effect at the mean of forest_total
b0 <- coef(model_part1)["(Intercept)"]
b1 <- coef(model_part1)["forest_total"]
xb_mean <- b0 + b1 * mean(df$forest_total)
me_mean <- dnorm(xb_mean) * b1
cat("\nMarginal effect at mean forest_total:", round(me_mean, 4), "\n")

# Classification accuracy
df <- df |>
  mutate(pred_zero = as.integer(p_structural >= 0.5),
         correct   = as.integer(pred_zero == structural_zero))
cat("\nClassification, forest only:\n")
print(table(Predicted = df$pred_zero, Actual = df$structural_zero))
cat("Overall accuracy:", round(mean(df$correct), 3), "\n")
df |>
  group_by(structural_zero) |>
  summarise(n = n(), accuracy = round(mean(correct), 3), .groups = "drop") |>
  as.data.frame() |>
  print()

# ── DIAGNOSTIC: ADD PASTURE, ONE TIME CHECK ONLY ─────────────────────────────
# This tests whether pasture adds structural zero information beyond forest.
# It is a check, not the model. The decision to keep forest alone rests on
# theory: structural incapacity is terrain, pasture belongs to the transition
# frame. Read whether pasture is significant here, then set it aside.

cat("\n=== DIAGNOSTIC ONLY: forest + pasture (NOT the model) ===\n")
df <- df |> mutate(pasture_z = as.numeric(scale(pasture_2011)))
check_p1 <- glm(structural_zero ~ forest_total + pasture_z,
                data = df, family = binomial(link = "probit"))
print(summary(check_p1))
cat("\nReading: if pasture is not significant here, this confirms forest alone\n")
cat("carries the structural signal. Even if it were significant, the model\n")
cat("keeps forest alone, because pasture belongs to the transition frame, not\n")
cat("to structural incapacity. This block is a check, not a specification.\n")

# ── SAVE ─────────────────────────────────────────────────────────────────────

saveRDS(model_part1, "data/processed/model_part1.rds")

sink("docs/part1_results.txt")
cat("=== PART 1 PROBIT, STRUCTURAL ZERO PROCESS ===\n")
cat("Predictor: forest_total. N = 114 counties.\n\n")
print(summary(model_part1))
cat("\nMarginal effect at mean forest_total:", round(me_mean, 4), "\n")
cat("Overall classification accuracy:", round(mean(df$correct), 3), "\n\n")
cat("Diagnostic (not the model): forest + pasture\n")
print(summary(check_p1))
sink()

# ── MAP: PREDICTED STRUCTURAL ZERO PROBABILITY ───────────────────────────────

counties_map <- counties_mo_proj |>
  filter(GEOID %in% df$GEOID) |>
  left_join(df |> select(GEOID, p_structural), by = "GEOID")

ggplot(counties_map) +
  geom_sf(aes(fill = p_structural), color = "white", linewidth = 0.3) +
  scale_fill_viridis_c(name = "P(structural\nzero)", option = "plasma",
                       limits = c(0, 1)) +
  labs(
    title    = "Predicted probability of structural zero, Part 1 probit",
    subtitle = "Missouri counties, forest cover predictor",
    caption  = "Predictor: total forest cover share, NLCD classes 41, 42, 43"
  ) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.text  = element_blank(),
        axis.ticks = element_blank())

ggsave("outputs/map_09a_structural_zero_probability.png",
       width = 10, height = 7, dpi = 300)

message("Script 09a complete. Part 1 probit estimated and saved.")
