# ── RHO PROGRESSION CHART FOR README ─────────────────────────────────────────
# Visualizes how spatial dependence (Rho) declines as the model correctly
# accounts for structural zeros and land cover opportunity variables

library(ggplot2)

rho_progression <- data.frame(
  step = factor(
    c("Standard\nspatial lag\n(n=114)",
      "ZI baseline\n(structural zeros\nremoved, n=91)",
      "+ Pasture\navailability",
      "+ CRP\nenrollment"),
    levels = c("Standard\nspatial lag\n(n=114)",
               "ZI baseline\n(structural zeros\nremoved, n=91)",
               "+ Pasture\navailability",
               "+ CRP\nenrollment")
  ),
  rho = c(0.645, 0.551, 0.173, -0.043),
  significant = c("Significant", "Significant", "Not significant", "Not significant")
)

ggplot(rho_progression, aes(x = step, y = rho, fill = significant)) +
  geom_col(width = 0.6) +
  geom_hline(yintercept = 0, color = "grey40", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.3f", rho),
                y = ifelse(rho >= 0, rho + 0.025, rho - 0.025)),
            size = 5, fontface = "bold") +
  scale_fill_manual(
    values = c("Significant" = "#d7191c", "Not significant" = "#2c7bb6"),
    name = NULL
  ) +
  scale_y_continuous(limits = c(-0.12, 0.75), 
                     breaks = seq(0, 0.6, 0.2)) +
  labs(
    title = "Apparent Spatial Contagion Fully Explained by Land Cover Opportunity",
    subtitle = "Spatial dependence parameter (Rho) across successive model refinements",
    x = NULL,
    y = "Rho (spatial dependence parameter)",
    caption = "Missouri cropland transitions, 2011–2021"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "top",
    panel.grid.major.x = element_blank(),
    plot.title = element_text(face = "bold", size = 14),
    plot.margin = margin(t = 20, r = 15, b = 25, l = 15)
  )

ggsave("outputs/chart_rho_progression.png",
       width = 9, height = 6.5, dpi = 300)