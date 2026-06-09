# Script 03: Economic Data Construction
# Purpose: Download county-level agricultural economic data from USDA NASS
#          Census of Agriculture via API
# Census years: 2002, 2012, 2022
# Alignment: Census year used as BEGINNING-OF-PERIOD economic conditions
# Census 2002 → explains 2001-2011 transition (measured at period start)
# Census 2012 → explains 2011-2021 transition (measured at period start)
# Census 2022 → not used as predictor (end-of-period, reverse causality risk)
# Data: USDA NASS Quick Stats API
# Author: Lan T. Tran
# Date: June 2026

library(dplyr)
library(tidyr)
library(httr)
library(jsonlite)
library(stringr)

# ── API KEY ───────────────────────────────────────────────────────────────────
# API key stored in scripts/api_key.R (excluded from Git via .gitignore)
source("scripts/api_key.R")

# ── DOWNLOAD FUNCTION ─────────────────────────────────────────────────────────
# Queries USDA NASS Quick Stats API for a single variable and year
# Filters to TOTAL domain only to avoid size class breakdowns
# Returns one row per Missouri county

get_nass <- function(short_desc, year) {
  
  url <- paste0(
    "https://quickstats.nass.usda.gov/api/api_GET/?key=", api_key,
    "&source_desc=CENSUS",
    "&sector_desc=ECONOMICS",
    "&agg_level_desc=COUNTY",
    "&state_alpha=MO",
    "&year=", year,
    "&short_desc=", URLencode(short_desc),
    "&format=JSON"
  )
  
  response <- GET(url)
  
  if (status_code(response) != 200) {
    message("Failed: ", short_desc, " ", year)
    return(NULL)
  }
  
  content_raw <- content(response, as = "text", encoding = "UTF-8")
  parsed <- fromJSON(content_raw)
  
  if (is.null(parsed$data) || nrow(parsed$data) == 0) {
    message("No data: ", short_desc, " ", year)
    return(NULL)
  }
  
  parsed$data |>
    filter(domain_desc == "TOTAL") |>
    group_by(county_name, state_fips_code, county_code, year) |>
    slice(1) |>
    ungroup() |>
    select(county_name, state_fips_code, county_code, year, Value) |>
    mutate(variable = short_desc)
}

# ── VARIABLES AND YEARS ───────────────────────────────────────────────────────

years <- c(2002, 2012, 2022)

variables <- c(
  "FARM OPERATIONS - NUMBER OF OPERATIONS",
  "FARM OPERATIONS - ACRES OPERATED",
  "INCOME, NET CASH FARM, OF OPERATIONS - NET INCOME, MEASURED IN $",
  "AG LAND, INCL BUILDINGS - ASSET VALUE, MEASURED IN $ / ACRE"
)

# ── DOWNLOAD ──────────────────────────────────────────────────────────────────

all_data <- list()

for (yr in years) {
  for (var in variables) {
    message("Downloading: ", var, " - ", yr)
    result <- get_nass(var, yr)
    if (!is.null(result)) {
      all_data[[paste(yr, var)]] <- result
    }
    Sys.sleep(0.5)
  }
}

raw_panel <- bind_rows(all_data)

# Confirm: 114 counties x 4 variables x 3 years = 1,368 rows
# 114 counties (St. Louis City excluded - no agricultural operations)
nrow(raw_panel)

# ── CLEAN VALUE COLUMN ────────────────────────────────────────────────────────
# NASS returns numeric values as strings with commas
# "(D)" indicates withheld data - converted to NA

raw_panel <- raw_panel |>
  mutate(Value = as.numeric(gsub(",", "", Value)))

# Check missing values
raw_panel |>
  group_by(variable) |>
  summarise(n_missing = sum(is.na(Value)), n_total = n())

# ── RESHAPE TO WIDE FORMAT ────────────────────────────────────────────────────

econ_wide <- raw_panel |>
  mutate(variable = case_when(
    grepl("NUMBER OF OPERATIONS", variable) ~ "n_farms",
    grepl("ACRES OPERATED",       variable) ~ "acres_operated",
    grepl("NET INCOME",           variable) ~ "net_income",
    grepl("ASSET VALUE",          variable) ~ "land_value_acre",
    TRUE ~ variable
  )) |>
  pivot_wider(
    id_cols = c(county_name, state_fips_code, county_code, year),
    names_from = variable,
    values_from = Value
  )

# ── CREATE GEOID AND ALIGN YEARS ──────────────────────────────────────────────
# GEOID = 5-digit FIPS code (state + county) for merging with spatial data
# nlcd_year maps Census year to corresponding NLCD snapshot year

econ_wide <- econ_wide |>
  mutate(
    GEOID = paste0(state_fips_code, str_pad(county_code, 3, pad = "0")),
    nlcd_year = case_when(
      year == 2002 ~ 2011,   # Census 2002 predicts 2001→2011 transition outcome
      year == 2012 ~ 2021,   # Census 2012 predicts 2011→2021 transition outcome
      year == 2022 ~ NA_real_ # not used as predictor — end of period, reverse causality risk
    )
  )

# ── EXPORT ────────────────────────────────────────────────────────────────────

write.csv(econ_wide,
          "data/processed/econ_panel_mo.csv",
          row.names = FALSE)

message("Economic panel saved: ", nrow(econ_wide), " rows")
# Expected: 342 rows (114 counties x 3 years)