# ── THRESHOLD SENSITIVITY NOTE ────────────────────────────────────────────────
# Transition threshold = 0.02 (2 percentage points) is the baseline definition
# Alternative thresholds tested:
#   0.01 → 100 transitions (43.5% of non-NA observations) — may capture noise
#   0.02 → 50 transitions (21.7%) — baseline choice
#   0.03 → 27 transitions (11.7%) — more conservative
#   0.05 → 2 transitions (0.9%)  — too restrictive
#
# Robustness checks in Phase 4 will re-estimate the zero-inflated spatial model
# under 0.01 and 0.03 thresholds to confirm findings are not threshold-dependent
# If results change substantially across thresholds, report all three

## Zero Rate by Period
- 2001-2011: 81.7% structural zeros
- 2011-2021: 74.8% structural zeros
- Both periods exceed the 70% threshold commonly cited as 
  justifying zero-inflated modeling approaches
- Transition rate increased between periods (18.3% to 25.2%)
  suggesting acceleration of cropland expansion in recent decade