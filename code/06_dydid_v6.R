# run_analysis.R
# Execution script for the portable Sun-Abraham / TWFE DiD pipeline
# Claude-generated update
# -------------------------------------------------------
# This script is the only file that should need to change when adapting the
# pipeline to a new dataset or new set of research questions. All estimation
# and output logic lives in portable_sunab.R.
#
# Workflow:
#   1. Environment setup
#   2. Dataset spec          — what do the columns mean in this dataset?
#   3. Analysis subset specs — which data slices do we analyse?
#   4. Outcome specs         — which outcome columns?
#   5. Treatment group specs — how do we partition treated units?
#   6. Model specs           — what formula / estimator?
#   7. Vcov specs            — which standard error estimators?
#   8. Run experiment
#   9. Rebuild merged tables
# -------------------------------------------------------


# ══════════════════════════════════════════════════════════════════════════════
# 1. Environment setup
# ══════════════════════════════════════════════════════════════════════════════

rm(list = ls())

if (!requireNamespace("here", quietly = TRUE)) install.packages("here")
library(here)

# Install any missing packages before sourcing the functions file
required_script_pkgs <- c(
  "tidyverse",   # dplyr, tidyr, purrr, readr, stringr, etc.
  "fixest",
  "arrow",
  "glue",
  "here"
)
missing_script_pkgs <- required_script_pkgs[
  !vapply(required_script_pkgs, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_script_pkgs) > 0) {
  install.packages(missing_script_pkgs)
}

library(tidyverse)

source(here::here("code", "portable_sunab3.R"))

seed <- 1234
set.seed(seed)


# ── Directory layout ──────────────────────────────────────────────────────────
# One shared dir_out per project. Multiple run_experiment() calls over time all
# write into the same directory; rebuild_estimation_tables() merges them later.

run_name  <- "GEE_resilience_v6_operational_ss500_ts50000"
results_version <- "v6"

# Adapt these paths to your project layout / CyVerse toggle as needed.
cyverse <- FALSE

if (cyverse) {
  dir_base    <- file.path(
    "~/data-store/data/iplant/home/shared/earthlab/macrosystems/tlm", run_name
  )
  dir_data    <- file.path(dir_base, "data", "derived")
  dir_raw     <- file.path(dir_base, "data", "raw")
  dir_manual  <- file.path(dir_base, "data", "manual")
  dir_results <- file.path(dir_base, "results", results_version)
  dir_figs    <- file.path(dir_base, "figs")
} else {
  dir_data    <- here::here("data", "derived", run_name)
  dir_raw     <- here::here("data", "raw")
  dir_manual  <- here::here("data", "manual")
  dir_results <- here::here("results", results_version)
  dir_figs    <- here::here("figs")
}

# The parquet data lives in a single directory; all filters are pushed to Arrow.
dir_parquet <- file.path(dir_data, "parquet")

dir_ensure_local(c(dir_data, dir_parquet, dir_raw, dir_manual, dir_results, dir_figs))




# ══════════════════════════════════════════════════════════════════════════════
# 2. Dataset spec
# ══════════════════════════════════════════════════════════════════════════════
# Captures all column-name mappings that are properties of the DATA, not of
# any particular model. Change this block when pointing at a new dataset.

dataset_spec <- make_dataset_spec(
  unit_id    = "pt_id",        # unit (pixel / plot) identifier
  time_var   = "year",         # calendar time
  trt_col    = "fire",         # binary treatment indicator (0/1)
  cohort_var = "FirstTreat",   # treatment cohort fed to sunab(); NA if plain TWFE
  event_id   = "fireid"        # event identifier; NA if not applicable
)


# ══════════════════════════════════════════════════════════════════════════════
# 3. Analysis subset specs
# ══════════════════════════════════════════════════════════════════════════════
# Each row = one data slice. Filters are pushed to Arrow before collect().
#
# Three approaches are shown; mix and match as needed.
#
#   A) Manually specify each subset with a tidy filter expression.
#   B) Auto-expand by unique values of a categorical column (e.g. ecoregion).
#   C) A single "all data" subset with only a time filter.
#
# Combine approaches with dplyr::bind_rows().

# ── A) Manual subsets — e.g. specific burn-year windows ──────────────────────

manual_subset_specs <- dplyr::bind_rows(
  make_analysis_subset_spec(
    subset_id   = "burnyear_2000_2009",
    data_source = dir_parquet,
    data_filter = ~ year >= 1997 & burn_year >= 2000 & burn_year < 2010
  ),
  make_analysis_subset_spec(
    subset_id   = "burnyear_2010_2019",
    data_source = dir_parquet,
    data_filter = ~ year >= 1997 & burn_year >= 2010 & burn_year < 2020
  )
)

# ── B) Auto-expand by ecoregion (or any other categorical column) ─────────────
# expand_analysis_subset_specs_by_col() reads the unique values of the split
# column directly from Arrow (no full collect), then builds one spec row per
# value. Use `values` to restrict to a known list; use `base_filter` to add a
# time or quality filter on top of each per-value filter.
#
# The ecoregion codes we want to include (mirrors the original script's
# forested_ecoregions tibble, minus those excluded for low USFS coverage):

ecoregion_code_names = c(
  "bluemtns", "cascades", "coastrange", "eastcascades",
  "klamathmtns", "northcascades", "pugetlowland",
  "willamettevalley", "centralcaliforniamtns", "sierranevada",
  "southerncaliforniamtns", "canadianrockies", "idahobatholith",
  "middlerockies", "northernrockies", "southernrockies",
  "wasatchuintamtns", "aznmmtns", "coloradoplateaus"
)

ecoregion_subset_specs <- expand_analysis_subset_specs_by_col(
  data_source = dir_parquet,
  split_col   = "ecoregion_code_name",
  id_prefix   = "ecor",
  base_filter = ~ year >= 1997,
  values      = ecoregion_code_names
)

# ── C) Single "all ecoregions pooled" subset ──────────────────────────────────

all_eco_subset_spec <- make_analysis_subset_spec(
  subset_id   = "all_ecoregions",
  data_source = dir_parquet,
  data_filter = ~ year >= 1997
)

# ── Combine into the full analysis_subset_specs table ─────────────────────────
# Comment out or add rows here to control which subsets are actually run.

analysis_subset_specs <- dplyr::bind_rows(
  ecoregion_subset_specs[1,]#,
  #all_eco_subset_spec
  # manual_subset_specs   # uncomment to also run burn-year window subsets
)

# Quick sanity check — print the run grid size before committing:
message(glue::glue("Analysis subsets defined: {nrow(analysis_subset_specs)}"))


# ══════════════════════════════════════════════════════════════════════════════
# 4. Outcome specs
# ══════════════════════════════════════════════════════════════════════════════

outcome_specs <- tibble::tibble(
  outcome = c("rap_tree",
              "vcf_tree"
              )
)


# ══════════════════════════════════════════════════════════════════════════════
# 5. Treatment group specs
# ══════════════════════════════════════════════════════════════════════════════
# Each row defines one way of partitioning treated units into subgroups.
#
# group_fun contract:
#   - Receives: df (data frame), group_col (string), plus everything in group_args
#   - Must return: df with a new column named exactly `group_col`
#   - The reference level ("f") is set inside set_cd_groups via relevel()
#
# group_col is specified at the top level (not inside group_args) so the
# pipeline can read it without knowing the internals of group_fun.

set_cd_groups <- function(df,
                          group_col,   # injected by the pipeline
                          b_nm,
                          d_nm,
                          b_threshold,
                          d_threshold) {
  
  if (d_threshold < 0) {
    df_new <- df |>
      dplyr::mutate(
        "{group_col}" := dplyr::case_when(
          .data[[b_nm]] <  b_threshold & fire == 1 & .data[[d_nm]] >  d_threshold ~ "f",
          .data[[b_nm]] >= b_threshold & fire == 1 & .data[[d_nm]] >  d_threshold ~ "bf",
          .data[[b_nm]] <  b_threshold & fire == 1 & .data[[d_nm]] <= d_threshold ~ "df",
          .data[[b_nm]] >= b_threshold & fire == 1 & .data[[d_nm]] <= d_threshold ~ "bdf",
          TRUE ~ NA_character_
        )
      )
  } else {
    df_new <- df |>
      dplyr::mutate(
        "{group_col}" := dplyr::case_when(
          .data[[b_nm]] <  b_threshold & fire == 1 & .data[[d_nm]] <  d_threshold ~ "f",
          .data[[b_nm]] >= b_threshold & fire == 1 & .data[[d_nm]] <  d_threshold ~ "bf",
          .data[[b_nm]] <  b_threshold & fire == 1 & .data[[d_nm]] >= d_threshold ~ "df",
          .data[[b_nm]] >= b_threshold & fire == 1 & .data[[d_nm]] >= d_threshold ~ "bdf",
          TRUE ~ NA_character_
        )
      )
  }
  
  # Set "f" (fire-only, low biotic stress, low drought) as the reference level
  df_new <- df_new |>
    dplyr::mutate(
      "{group_col}" := relevel(factor(.data[[group_col]]), ref = "f")
    )
  
  df_new
}


treatment_group_specs <- dplyr::bind_rows(
  make_treatment_group_spec(
    group_id  = "b10_pdsin4t1",
    group_col = "b10_pdsin4t1",
    group_fun = set_cd_groups,
    group_args = list(
      b_nm        = "biotic_relaxedforestnorm_5_yrs_prior_sum_yot",
      d_nm        = "pdsi_annual_5_yrs_prior_threshold_n4_yot",
      b_threshold = 10,
      d_threshold = 1
    )
  ),
  make_treatment_group_spec(
    group_id  = "b10_pdsin3t1",
    group_col = "b10_pdsin3t1",
    group_fun = set_cd_groups,
    group_args = list(
      b_nm        = "biotic_relaxedforestnorm_5_yrs_prior_sum_yot",
      d_nm        = "pdsi_annual_5_yrs_prior_threshold_n3_yot",
      b_threshold = 10,
      d_threshold = 1
    )
  ),
  make_treatment_group_spec(
    group_id  = "b10_pdsisumn10",
    group_col = "b10_pdsisumn10",
    group_fun = set_cd_groups,
    group_args = list(
      b_nm        = "biotic_relaxedforestnorm_5_yrs_prior_sum_yot",
      d_nm        = "pdsi_annual_5_yrs_prior_sum_yot",
      b_threshold = 10,
      d_threshold = -10
    )
  ),
  make_treatment_group_spec(
    group_id  = "b25_pdsisumn10",
    group_col = "b25_pdsisumn10",
    group_fun = set_cd_groups,
    group_args = list(
      b_nm        = "biotic_relaxedforestnorm_5_yrs_prior_sum_yot",
      d_nm        = "pdsi_annual_5_yrs_prior_sum_yot",
      b_threshold = 25,
      d_threshold = -10
    )
  )
)


treatment_group_specs <- treatment_group_specs[3,]

# ══════════════════════════════════════════════════════════════════════════════
# 6. Model specs
# ══════════════════════════════════════════════════════════════════════════════
# formula_template uses {outcome} as the only glue placeholder.
# estimator_type controls post-processing: "sunab" extracts event-time support;
# "twfe" skips that step.
#
# Both a Sun-Abraham spec and a plain TWFE spec are defined here. Comment out
# whichever you don't need, or pass just one to run_experiment().

model_specs <- dplyr::bind_rows(
  
  # ── Sun-Abraham (preferred for staggered treatment) ────────────────────────
  make_model_spec(
    model_id         = "b_pdsi_sunab_twfe",
    formula_template = paste0(
      "{outcome} ~ biotic_relaxedforestnorm + pdsi_annual",
      " + sunab(FirstTreat, year, ref.p = -6)",
      " | pt_id + year"
    ),
    estimator_type   = "sunab",
    term_pattern     = "^year::"   # regex flagging event-time coefficients
  ),
  
  # ── Plain TWFE (for comparison / robustness checks) ────────────────────────
  make_model_spec(
    model_id         = "b_pdsi_twfe",
    formula_template = paste0(
      "{outcome} ~ biotic_relaxedforestnorm + pdsi_annual",
      " + i(rel_year, ref = -6)",   # assumes a pre-built rel_year column in data
      " | pt_id + year"
    ),
    estimator_type   = "twfe",
    term_pattern     = "^rel_year::"
  )
  
)

model_specs <- model_specs[1,]


# ══════════════════════════════════════════════════════════════════════════════
# 7. Vcov specs
# ══════════════════════════════════════════════════════════════════════════════
# vcov_vars lists every column the SE estimator needs (used to ensure those
# columns are retained before the model is fitted).
# lat/lon column names stay here (not in dataset_spec) for flexibility.

vcov_specs <- tibble::tibble(
  vcov_id = c(
    "iid",
    "cluster_pt",
    "cluster_eco4",
    "conley_125km_20km"
  ),
  vcov_label = c(
    "IID Classical",
    "Clustered SEs: pt_id",
    "Clustered SEs: us_l4code",
    "Conley SEs: 125 km cutoff, 20 km pixel"
  ),
  vcov = list(
    "iid",
    stats::as.formula("~ pt_id"),
    stats::as.formula("~ us_l4code"),
    fixest::vcov_conley(
      lat      = "lat",
      lon      = "long",
      cutoff   = 125,
      pixel    = 20,
      distance = "triangular"
    )
  ),
  vcov_vars = list(
    character(0),
    "pt_id",
    "us_l4code",
    c("lat", "long")
  )
)

vcov_specs <- vcov_specs[3,]

# ══════════════════════════════════════════════════════════════════════════════
# 8. Optional: group palette
# ══════════════════════════════════════════════════════════════════════════════
# Named character vector mapping subgroup labels to hex colours used in plots.
# NULL is fine if you don't use colours downstream.

group_palette <- c(
  "f"   = "#E69F00",   # fire only
  "bf"  = "#56B4E9",   # high biotic stress + fire
  "df"  = "#009E73",   # drought + fire
  "bdf" = "#CC79A7"    # high biotic stress + drought + fire
)


# ══════════════════════════════════════════════════════════════════════════════
# 9. Preview the run grid
# ══════════════════════════════════════════════════════════════════════════════
# Inspect what will be run before committing. Nothing is executed here.

preview_grid <- tidyr::crossing(
  analysis_subset_specs |> dplyr::select(subset_id),
  outcome_specs,
  treatment_group_specs |> dplyr::select(group_id),
  model_specs |> dplyr::select(model_id, estimator_type)
) |>
  dplyr::mutate(
    run_id = glue::glue("{subset_id}__{outcome}__{group_id}__{model_id}")
  )

message(glue::glue("Total runs planned: {nrow(preview_grid)}"))
print(preview_grid |> dplyr::select(subset_id, outcome, group_id, model_id, estimator_type))


# ══════════════════════════════════════════════════════════════════════════════
# 10. Run experiment — ecoregion subsets, Sun-Abraham only
# ══════════════════════════════════════════════════════════════════════════════
# Split into multiple run_experiment() calls if needed (e.g. different subsets
# on different days, or a subset of models for a quick robustness check).
# All calls write into the same dir_results; rebuild() merges them.

results_sunab <- run_experiment(
  dataset_spec           = dataset_spec,
  analysis_subset_specs  = analysis_subset_specs,
  outcome_specs          = outcome_specs,
  treatment_group_specs  = treatment_group_specs,
  model_specs            = model_specs |> dplyr::filter(estimator_type == "sunab"),
  vcov_specs             = vcov_specs,
  dir_out                = dir_results,
  group_palette          = group_palette,
  ci_level               = 0.95,
  run_estimation         = TRUE,
  run_descriptive        = TRUE,
  descriptive_args       = list(
    treated_year_var = "burn_year",
    control_year_var = "mock_burn_year"
  ),
  skip_existing          = TRUE,
  write_vcov             = FALSE,
  verbose_timing         = TRUE,
  .progress              = TRUE
)

# Inspect which runs failed (if any):
failed_sunab <- purrr::keep(results_sunab$run_results, \(r) !is.null(r$error))
if (length(failed_sunab) > 0) {
  message("Failed runs:")
  purrr::walk(failed_sunab, \(r) message("  ", r$run_id, ": ", r$error))
}


# ══════════════════════════════════════════════════════════════════════════════
# 11. Run experiment — burn-year window subsets (manual), both estimators
# ══════════════════════════════════════════════════════════════════════════════
# Example of a second experiment call writing into the same dir_results.
# Uncomment when ready to run.

# results_burnyear <- run_experiment(
#   dataset_spec           = dataset_spec,
#   analysis_subset_specs  = manual_subset_specs,
#   outcome_specs          = outcome_specs,
#   treatment_group_specs  = treatment_group_specs,
#   model_specs            = model_specs,   # runs both sunab and twfe
#   vcov_specs             = vcov_specs,
#   dir_out                = dir_results,
#   group_palette          = group_palette,
#   ci_level               = 0.95,
#   run_estimation         = TRUE,
#   run_descriptive        = FALSE,
#   skip_existing          = TRUE,
#   write_vcov             = FALSE,
#   verbose_timing         = TRUE,
#   .progress              = TRUE
# )


# ══════════════════════════════════════════════════════════════════════════════
# 12. Rebuild merged output tables
# ══════════════════════════════════════════════════════════════════════════════
# Call rebuild at any point after one or more run_experiment() calls.
# It scans all by_run folders under dir_results, deduplicates on run_id
# (keeping the most recently written file), and writes merged tables to
# dir_results/tables/all/ and dir_results/descriptive/all/.

all_estimation_tables <- rebuild_estimation_tables(
  dir_out   = dir_results,
  write_csv = TRUE
)

all_descriptive_tables <- rebuild_descriptive_tables(
  dir_out   = dir_results,
  write_csv = TRUE
)

# Quick summary of what was rebuilt:
message(glue::glue(
  "Estimation table rows: {nrow(all_estimation_tables$coef_expanded_vcov)}\n",
  "Unique run_ids:        {dplyr::n_distinct(all_estimation_tables$run_registry$run_id)}\n",
  "Descriptive rows:      {nrow(all_descriptive_tables$event_time_trajectory)}"
))
