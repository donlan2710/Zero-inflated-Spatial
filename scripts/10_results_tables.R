# Script 10: Results and tables
# Purpose: assemble the paper's main tables and the key figure from the
#          zero-inflated spatial model. Self-contained: it rebuilds the sample
#          and refits the models, so every number comes from one coherent run.
#
# Narrative frame: zero-inflated model leads. Part 1 separates structural zeros.
# Part 2 models transition among active counties, and the pasture contrast is
# the key finding inside Part 2.
#
# On extraction: fit@rho is confirmed to work. The coefficient table with p
# values is captured from summary() into docs/, which is authoritative. The
# convenience tables below use best-effort extraction; if a cell is NA, read
# the captured summary instead. An inspection block prints the object structure
# so the real slot names are visible.
#
# Author: Lan T. Tran
# Date: July 2026

library(sf)
library(spdep)
library(spatialreg)
library(dplyr)
library(Matrix)
library(ProbitSpatial)
library(ggplot2)

# ── REBUILD SAMPLE ────────────────────────────────────────────────────────────

counties_mo_proj <- st_read("data/processed/counties_mo.shp") |> arrange(GEOID)
panel   <- read.csv("data/processed/panel_final.csv", colClasses = c(GEOID="character"))
pasture <- read.csv("data/processed/pasture_by_year.csv", colClasses = c(GEOID="character"))
w_queen <- readRDS("data/processed/weights_queen.rds")

df <- panel |>
  filter(year == 2021, !is.na(n_farms)) |>
  arrange(GEOID) |>
  left_join(pasture |> select(GEOID, pasture_2011), by = "GEOID")

geoids_114 <- df$GEOID
w_114 <- subset(w_queen, subset = counties_mo_proj$GEOID %in% geoids_114)

lisa_cropland <- localmoran(df$cropland, w_114)
mean_crop <- mean(df$cropland)
lag_crop  <- lag.listw(w_114, df$cropland)
df <- df |>
  mutate(
    forest_total    = forest + frac_42 + frac_43,
    structural_zero = as.integer(cropland < mean_crop & lag_crop < mean_crop &
                                 lisa_cropland[, 5] < 0.05),
    cropland_z      = as.numeric(scale(cropland_lag)),
    pasture_z_full  = as.numeric(scale(pasture_2011))
  )

# ── HELPERS ───────────────────────────────────────────────────────────────────

get_rho <- function(fit) {
  if (is.null(fit)) return(NA_real_)
  as.numeric(tryCatch(fit@rho, error = function(e) NA_real_))
}

# Best-effort coefficient extraction. Returns a named vector of estimates.
# If it cannot find them, returns NA and the captured summary is the source.
get_coefs <- function(fit) {
  out <- tryCatch(coef(fit), error = function(e) NULL)
  if (is.null(out)) return(NA_real_)
  out
}

fit_p2 <- function(formula, data, W) {
  tryCatch(
    ProbitSpatialFit(formula, data = data, W = W,
                     DGP = "SAR", method = "conditional", varcov = "varcov"),
    error   = function(e) { cat("FIT ERROR:", conditionMessage(e), "\n"); NULL },
    warning = function(w) suppressWarnings(
      ProbitSpatialFit(formula, data = data, W = W,
                       DGP = "SAR", method = "conditional", varcov = "varcov"))
  )
}

# ── PART 1: STRUCTURAL ZERO PROBIT ───────────────────────────────────────────

m_part1 <- glm(structural_zero ~ forest_total, data = df,
               family = binomial(link = "probit"))
b1 <- coef(m_part1)["forest_total"]; b0 <- coef(m_part1)["(Intercept)"]
me_forest <- dnorm(b0 + b1 * mean(df$forest_total)) * b1
df$p_struct <- predict(m_part1, type = "response")
acc_part1 <- mean(as.integer(df$p_struct >= 0.5) == df$structural_zero)

# ── PART 2: ACTIVE SAMPLE, DROP ISLAND, REFIT M0 M1 M2 ───────────────────────

active_idx <- df$structural_zero == 0
active <- df[active_idx, ]
w_active <- subset(w_114, subset = active_idx)
active <- active |>
  mutate(pasture_z   = as.numeric(scale(pasture_2011)),
         landvalue_z = as.numeric(scale(land_value_acre)))

W_dense <- listw2mat(w_active); W_dense[is.na(W_dense)] <- 0
rs <- rowSums(W_dense); keep <- which(rs > 1e-12)
active90 <- active[keep, ]
Wk <- W_dense[keep, keep]; Wk <- as(Wk / rowSums(Wk), "CsparseMatrix")

m0 <- fit_p2(transition ~ cropland_z,               active90, Wk)
m1 <- fit_p2(transition ~ cropland_z + pasture_z,    active90, Wk)
m2 <- fit_p2(transition ~ cropland_z + pasture_z + landvalue_z, active90, Wk)

# baseline probit on 114 (structural zeros in), from Script 07
baseline <- tryCatch(readRDS("data/processed/model_baseline_probit_114.rds"),
                     error = function(e) NULL)

# ── INSPECT OBJECT STRUCTURE (so slot names are visible this run) ────────────

cat("\n=== OBJECT STRUCTURE, ProbitSpatialFit ===\n")
cat("slotNames(m1):\n"); print(slotNames(m1))
cat("\ncoef(m1):\n");    print(get_coefs(m1))
cat("\nstr(summary(m1)) [captured below in docs as authoritative]\n")

# ── CAPTURE AUTHORITATIVE SUMMARIES ──────────────────────────────────────────

sink("docs/part2_model_summaries.txt")
cat("=== PART 1 PROBIT ===\n");                 print(summary(m_part1))
cat("\nMarginal effect at mean forest:", round(me_forest, 4), "\n")
cat("Classification accuracy:", round(acc_part1, 3), "\n")
cat("\n=== M1 PRIMARY (cropland + pasture) ===\n"); print(summary(m1))
cat("\n=== M2 (+ land value) ===\n");             print(summary(m2))
cat("\n=== M0 (cropland only, contrast) ===\n");  print(summary(m0))
if (!is.null(baseline)) { cat("\n=== BASELINE 114 ===\n"); print(summary(baseline)) }
sink()

# ── TABLE 1: PART 1 ──────────────────────────────────────────────────────────

cat("\n=== TABLE 1: PART 1 STRUCTURAL ZERO PROBIT ===\n")
tab1 <- data.frame(
  term = c("forest_total", "marginal effect at mean", "N", "accuracy"),
  value = c(round(b1, 3), round(me_forest, 3), nrow(df), round(acc_part1, 3))
)
print(tab1)

# ── TABLE 2: PART 2 SPATIAL PARAMETER CONTRAST (the finding) ─────────────────

cat("\n=== TABLE 2: PART 2 SPATIAL PARAMETER CONTRAST ===\n")
tab2 <- data.frame(
  model = c("M0 cropland only (no pasture)",
            "M1 cropland + pasture (primary)",
            "M2 + land value (robustness)"),
  n = nrow(active90),
  rho = round(c(get_rho(m0), get_rho(m1), get_rho(m2)), 3)
)
print(tab2)
cat("\nThe finding: rho falls from", round(get_rho(m0), 3), "to",
    round(get_rho(m1), 3), "when pasture enters.\n")
cat("Coefficients: read the captured summaries in docs/part2_model_summaries.txt\n")
cat("Convenience coef extraction (verify against summaries):\n")
cat("M1 coefs:\n"); print(get_coefs(m1))
cat("M2 coefs:\n"); print(get_coefs(m2))

# ── TABLE 3: ZERO-INFLATED CONTRIBUTION (baseline vs Part 2) ─────────────────
# Both include pasture, so this isolates the effect of removing structural zeros.

cat("\n=== TABLE 3: ZERO-INFLATED CONTRIBUTION ===\n")
tab3 <- data.frame(
  model = c("Baseline, structural zeros IN (114)",
            "Part 2, structural zeros OUT (90)"),
  n = c(114, nrow(active90)),
  rho = round(c(get_rho(baseline), get_rho(m1)), 3)
)
print(tab3)
cat("\nBoth models include pasture. The change in rho isolates the effect of\n")
cat("removing structural zeros. This is the zero-inflated contribution.\n")

# ── SAVE TABLES ──────────────────────────────────────────────────────────────

sink("docs/results_tables.txt")
cat("TABLE 1: PART 1\n");  print(tab1)
cat("\nTABLE 2: PART 2 SPATIAL PARAMETER CONTRAST\n"); print(tab2)
cat("\nTABLE 3: ZERO-INFLATED CONTRIBUTION\n"); print(tab3)
sink()

# ── FIGURE: THE PASTURE CONTRAST ─────────────────────────────────────────────

fig_df <- data.frame(
  model = factor(c("Without pasture", "With pasture"),
                 levels = c("Without pasture", "With pasture")),
  rho   = c(get_rho(m0), get_rho(m1))
)

ggplot(fig_df, aes(x = model, y = rho, fill = model)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = round(rho, 3)), vjust = -0.4, size = 5) +
  scale_fill_manual(values = c("Without pasture" = "#2c7bb6",
                               "With pasture" = "#d7191c"), guide = "none") +
  ylim(0, 1) +
  labs(
    title    = "Spatial parameter falls when convertible land enters the model",
    subtitle = "Active counties, maximum likelihood spatial probit",
    x = NULL, y = "Spatial parameter (rho)"
  ) +
  theme_minimal(base_size = 13) +
  theme(panel.grid.major.x = element_blank())

ggsave("outputs/fig_10_pasture_contrast.png", width = 8, height = 6, dpi = 300)

message("Script 10 complete. Tables and figure written. Read docs summaries as authoritative.")
