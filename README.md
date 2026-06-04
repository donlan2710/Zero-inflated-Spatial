# Zero-inflated-Spatial
Extending zero-inflated modeling into spatial econometrics: separating the structural non-transition process from the transition process in land use data. Companion to gRDA-Optimizer in addressing sparsity across different methodological contexts.

## Overview
Standard spatial econometric models assume well-behaved outcome distributions.
Land use transition data and many other spatial outcomes in economics are
structurally sparse: most observations report no change in any given period.
This project develops and applies a zero-inflated spatial econometric framework
that models the zero-generating process separately from the transition process.

Missouri county-level cropland transitions from the National Land Cover Database
(2001–2021) serve as the primary application. The methodological contribution
travels to any sparse spatial outcome: flood adaptation, environmental compliance,
agricultural technology adoption.

## Repository Structure
- `/data/raw` — original downloaded files (not version controlled if large)
- `/data/processed` — cleaned analytical datasets
- `/scripts` — numbered R scripts in execution order
- `/outputs` — tables, figures, maps
- `/docs` — data construction notes, model specification, diagnostic memos

## Data Sources
- National Land Cover Database (NLCD): https://www.mrlc.gov
- USDA ERS county economic data: https://www.ers.usda.gov
- Census TIGER/Line county shapefiles: https://www.census.gov/geographies/mapping-files.html

## Software
R version 4.x. Key packages: sf, terra, exactextractr, spdep, spatialreg.
Full session information: see `docs/session_info.txt`

## Related Work
- [gRDA-Optimizer](https://github.com/YOUR_USERNAME/gRDA-Optimizer):
  sparsity-promoting methods in deep learning — the same sparsity problem
  in a different methodological context
- Tran, L (2026). Measuring what matters: How Zero-Inflated Latent Class Analysis reveals hidden heterogeneity in U.S.
  Household Climate Vulnerability. AAEA Annual Meeting, Kansas City, MO. July 2026 (forthcoming).
- Tran, L. et al. (2023). Measuring pesticide overuse: Evidence from Vietnamese
  rice and fruit farms. AJARE, 67(4), 417–437.

## Status
Active development. Started June 2026.
