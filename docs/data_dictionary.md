# Data Dictionary: panel_final.csv
**File location:** `data/processed/panel_final.csv`  
**Dimensions:** 345 rows × 28 columns (115 Missouri counties × 3 years)  
**Unit of observation:** County-year  
**Years:** 2001, 2011, 2021  
**Last updated:** June 2026  

---

## Identifiers

| Column | Type | Description |
|---|---|---|
| `GEOID` | integer | 5-digit county FIPS code (state + county). Primary merge key. |
| `county_name` | character | County name in uppercase (from USDA NASS) |
| `state_fips_code` | integer | State FIPS code (29 = Missouri throughout) |
| `county_code` | integer | 3-digit county FIPS code |
| `year` | integer | NLCD snapshot year: 2001, 2011, or 2021 |
| `year.x` | integer | Census of Agriculture year aligned to NLCD year (2002, 2012, 2022) — redundant with `year`, may be dropped in analysis |

---

## Land Cover Variables
**Source:** National Land Cover Database (NLCD), years 2001, 2011, 2021  
**Method:** Proportions extracted to county polygons using `exactextractr`  
**Unit:** Proportion of county area (0 to 1); all land cover columns sum to 1.0 per county-year  
**NLCD class codes:** See [mrlc.gov](https://www.mrlc.gov) for full classification

| Column | NLCD Class | Description |
|---|---|---|
| `water` | 11 | Open water |
| `developed` | 21 | Developed, open space (low impervious surface) |
| `frac_22` | 22 | Developed, low intensity |
| `frac_23` | 23 | Developed, medium intensity |
| `frac_24` | 24 | Developed, high intensity |
| `frac_31` | 31 | Barren land (rock, sand, clay) |
| `forest` | 41 | Deciduous forest |
| `frac_42` | 42 | Evergreen forest |
| `frac_43` | 43 | Mixed forest |
| `frac_52` | 52 | Shrub/scrub |
| `frac_71` | 71 | Grassland/herbaceous |
| `pasture` | 81 | Pasture/hay |
| `cropland` | 82 | Cultivated crops — **primary outcome variable** |
| `wetland` | 90 | Woody wetlands |
| `frac_95` | 95 | Emergent herbaceous wetlands |

---

## Economic Variables
**Source:** USDA NASS Census of Agriculture  
**Census years aligned to NLCD years:**
- Census 2002 → NLCD 2001
- Census 2012 → NLCD 2011
- Census 2022 → NLCD 2021

**Note:** St. Louis City (GEOID 29510) has NA for all economic variables — it is an independent city with no agricultural operations recorded in the Census of Agriculture.

| Column | Unit | Description |
|---|---|---|
| `n_farms` | count | Number of farm operations in the county |
| `acres_operated` | acres | Total acres operated by farms in the county |
| `net_income` | USD | Net cash farm income of operations, total for county |
| `land_value_acre` | USD/acre | Average agricultural land and buildings asset value per acre |

---

## Transition Variables
**Constructed in:** `04_panel_construction.R`  
**Based on:** Change in `cropland` proportion between consecutive NLCD periods

| Column | Type | Description |
|---|---|---|
| `cropland_lag` | numeric | Cropland proportion in the previous NLCD period (NA for 2001 — no prior period) |
| `cropland_change` | numeric | Change in cropland proportion: `cropland - cropland_lag` (NA for 2001) |
| `transition` | integer | Binary indicator: 1 if `abs(cropland_change) > 0.02`, 0 otherwise, NA for 2001 |
| `transition_dir` | integer | Direction: 1 = cropland gain, -1 = cropland loss, 0 = no transition, NA for 2001 |

---

## Threshold Decision
The 0.02 (2 percentage point) threshold for `transition` is the baseline definition.  
See `docs/threshold_decision.md` for full sensitivity analysis and justification.

**Zero rates at baseline threshold:**
- 2001–2011 period: 81.7% zeros
- 2011–2021 period: 74.8% zeros

Both exceed the 70% benchmark commonly cited for zero-inflated modeling.

---

## Notes
- All land cover proportions sum to 1.0 per county-year (validated in `04_panel_construction.R`)
- All transitions at the 0.02 threshold are cropland gains — no county lost more than 2 percentage points of cropland over any 10-year period
- The `total` validation column from `02_raster_extraction.R` was dropped during panel construction
