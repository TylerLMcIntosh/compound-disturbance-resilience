# run_analysis.R
# Execution script for the portable Sun-Abraham / TWFE DiD pipeline
# -------------------------------------------------------
# This script is the only file that should need to change when adapting the
# pipeline to a new dataset or new set of research questions. All estimation
# and output logic lives in portable_sunab.R.
#
# Workflow:
#   1.  Environment setup
#   2.  Dataset spec           — what do the columns mean in this dataset?
#   3.  Analysis subset specs  — which data slices do we analyse?
#   4.  Outcome specs          — which outcome columns?
#   5.  Treatment group specs  — how do we partition treated units?
#   6.  Weighting specs        — which propensity-score weighting approaches?
#   7.  Model specs            — what formula / estimator / weight column?
#   8.  Vcov specs             — which standard error estimators?
#   9.  Preview run grid
#   10. Run weighting experiment
#   11. Run estimation experiment
#   12. Rebuild merged tables
# -------------------------------------------------------


# ══════════════════════════════════════════════════════════════════════════════
# 1. Environment setup
# ══════════════════════════════════════════════════════════════════════════════

rm(list = ls())

if (!requireNamespace("here", quietly = TRUE)) install.packages("here")
library(here)

required_script_pkgs <- c(
  "tidyverse", "fixest", "arrow", "glue", "here"
)
missing_script_pkgs <- required_script_pkgs[
  !vapply(required_script_pkgs, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_script_pkgs) > 0) install.packages(missing_script_pkgs)

library(tidyverse)

source(here::here("code", "portable_sunab3.R"))

seed <- 1234
set.seed(seed)


# ── Directory layout ──────────────────────────────────────────────────────────
# One shared dir_results per project. Multiple run_experiment() calls over time
# all write into the same directory; rebuild_*() merges them later.

run_name <- "GEE_resilience_v6_operational_ss500_ts50000"
version <- "v6"

cyverse <- FALSE

if (cyverse) {
  dir_base    <- file.path(
    "~/data-store/data/iplant/home/shared/earthlab/macrosystems/tlm", run_name
  )
  dir_data    <- file.path(dir_base, "data", "derived")
  dir_raw     <- file.path(dir_base, "data", "raw")
  dir_manual  <- file.path(dir_base, "data", "manual")
  dir_results <- file.path(dir_base, "results", version)
  dir_figs    <- file.path(dir_base, "figs", version)
} else {
  dir_data    <- here::here("data", "derived", run_name)
  dir_raw     <- here::here("data", "raw")
  dir_manual  <- here::here("data", "manual")
  dir_results <- here::here("results", version)
  dir_figs    <- here::here("figs", version)
}

# Long (panel) data — one row per unit × year, used for estimation
dir_parquet_long <- file.path(dir_data, "parquet_long_filtered")

# Short (cross-sectional) data — one row per unit, used for weighting.
# Set to NULL and omit short_data_source in subset specs when not weighting.
dir_parquet_short <- file.path(dir_data, "parquet_short_filtered")

dir_ensure_local(c(dir_data, dir_parquet_long, dir_raw, dir_manual, dir_results, dir_figs))


# ══════════════════════════════════════════════════════════════════════════════
# 2. Dataset spec
# ══════════════════════════════════════════════════════════════════════════════

dataset_spec <- make_dataset_spec(
  unit_id    = "pt_id",       # unit (pixel / plot) identifier
  time_var   = "year",        # calendar time
  trt_col    = "fire",        # binary treatment indicator (0/1)
  cohort_var = "FirstTreat",  # treatment cohort for sunab(); NA = plain TWFE
  event_id   = "fireid"       # event identifier; NA if not applicable
)


# ══════════════════════════════════════════════════════════════════════════════
# 3. Analysis subset specs
# ══════════════════════════════════════════════════════════════════════════════
# Each row = one data slice. data_filter is pushed to Arrow before collect()
# for both long and short sources.
#
# long_data_source  — required; the panel data used for estimation.
# short_data_source — optional; the cross-sectional data used for weighting.
#                     Set to NULL (or omit) when not weighting.
# data_filter       — applied to BOTH sources.
#
# Three approaches shown; combine with dplyr::bind_rows().
#
#   A) Auto-expand by unique values of a categorical column (ecoregion).
#   B) Single "all data" subset.
#   C) Manual subsets (e.g. burn-year windows).

# Ecoregion codes to include (forested ecoregions with sufficient USFS coverage)
ecoregion_code_names = c(
  "bluemtns", "cascades", "coastrange", "eastcascades",
  "klamathmtns", "northcascades", "pugetlowland",
  "willamettevalley", "centralcaliforniamtns", "sierranevada",
  "southerncaliforniamtns", "canadianrockies", "idahobatholith",
  "middlerockies", "northernrockies", "southernrockies",
  "wasatchuintamtns", "aznmmtns", "coloradoplateaus"
)


# ── A) Auto-expand by ecoregion ───────────────────────────────────────────────
# expand_analysis_subset_specs_by_col() reads unique values from the long data
# via Arrow (no full collect), then builds one row per value.
# short_data_source is the same directory for every ecoregion subset.

ecoregion_subset_specs <- expand_analysis_subset_specs_by_col(
  long_data_source  = dir_parquet_long,
  split_col         = "ecoregion_code_name",
  id_prefix         = "ecor",
  #base_filter       = ~ year >= 1997,
  values            = ecoregion_code_names,
  short_data_source = dir_parquet_short   # NULL if no weighting planned
)

# ── B) Auto-expand by nfg_group ───────────────────────────────────────────────

forestgroup_subset_specs <- expand_analysis_subset_specs_by_col(
  long_data_source  = dir_parquet_long,
  split_col         = "nfg_factor",
  id_prefix         = "nfg",
  #base_filter       = ~ year >= 1997,
  #values            = nfg_values,
  short_data_source = dir_parquet_short,   # NULL if no weighting planned
  check_all_files   = TRUE
)

# ── C) All ecoregions pooled ──────────────────────────────────────────────────

all_eco_subset_spec <- make_analysis_subset_spec(
  subset_id         = "all_ecoregions",
  long_data_source  = dir_parquet_long,
  #data_filter       = ~ year >= 1997,
  short_data_source = dir_parquet_short   # NULL if no weighting planned
)

# ── D) Manual subsets — e.g. specific burn-year windows ──────────────────────
# short_data_source can be NULL here if weighting is not needed for these
# subsets; they will simply be skipped by run_weighting_experiment().

manual_subset_specs <- dplyr::bind_rows(
  make_analysis_subset_spec(
    subset_id        = "burnyear_2000_2009",
    long_data_source = dir_parquet_long,
    data_filter      = ~ year >= 1997 & burn_year >= 2000 & burn_year < 2010
    # short_data_source = NULL  (default — no weighting for these subsets)
  ),
  make_analysis_subset_spec(
    subset_id        = "burnyear_2010_2019",
    long_data_source = dir_parquet_long,
    data_filter      = ~ year >= 1997 & burn_year >= 2010 & burn_year < 2020
  )
)

# ── Combine ───────────────────────────────────────────────────────────────────

analysis_subset_specs <- dplyr::bind_rows(
  ecoregion_subset_specs,
  #all_eco_subset_spec
  # manual_subset_specs   # uncomment to also run burn-year window subsets
)

message(glue::glue("Analysis subsets defined: {nrow(analysis_subset_specs)}"))

analysis_subset_specs <- analysis_subset_specs[11,]

# ══════════════════════════════════════════════════════════════════════════════
# 4. Outcome specs
# ══════════════════════════════════════════════════════════════════════════════

outcome_specs <- tibble::tibble(
  outcome = c("rap_tree"#,
              #"vcf_tree"
              )
)


# ══════════════════════════════════════════════════════════════════════════════
# 5. Treatment group specs
# ══════════════════════════════════════════════════════════════════════════════
# group_fun contract:
#   Receives: df, group_col (injected), plus everything in group_args.
#             When include_control = TRUE (injected by weighting pipeline),
#             must assign "control" to rows where fire == 0.
#   Returns:  df with a new column named exactly group_col.
#   Reference level "f" is set inside set_cd_groups via relevel().

set_cd_groups <- function(df,
                          group_col,        # injected by pipeline
                          b_nm,
                          d_nm,
                          b_threshold,
                          d_threshold,
                          include_control = FALSE) {
  
  if (d_threshold < 0) {
    df_new <- df |>
      dplyr::mutate(
        "{group_col}" := dplyr::case_when(
          fire == 1 & .data[[b_nm]] <  b_threshold & .data[[d_nm]] >  d_threshold ~ "f",
          fire == 1 & .data[[b_nm]] >= b_threshold & .data[[d_nm]] >  d_threshold ~ "bf",
          fire == 1 & .data[[b_nm]] <  b_threshold & .data[[d_nm]] <= d_threshold ~ "df",
          fire == 1 & .data[[b_nm]] >= b_threshold & .data[[d_nm]] <= d_threshold ~ "bdf",
          include_control & fire == 0                                               ~ "control",
          TRUE ~ NA_character_
        )
      )
  } else {
    df_new <- df |>
      dplyr::mutate(
        "{group_col}" := dplyr::case_when(
          fire == 1 & .data[[b_nm]] <  b_threshold & .data[[d_nm]] <  d_threshold ~ "f",
          fire == 1 & .data[[b_nm]] >= b_threshold & .data[[d_nm]] <  d_threshold ~ "bf",
          fire == 1 & .data[[b_nm]] <  b_threshold & .data[[d_nm]] >= d_threshold ~ "df",
          fire == 1 & .data[[b_nm]] >= b_threshold & .data[[d_nm]] >= d_threshold ~ "bdf",
          include_control & fire == 0                                               ~ "control",
          TRUE ~ NA_character_
        )
      )
  }
  
  # Set "f" (fire-only, low biotic stress, low drought) as the reference level.
  # "control" is not in the reference-level logic because WeightIt uses
  # group_col as a multi-valued treatment indicator, not for feols().
  df_new <- df_new |>
    dplyr::mutate(
      "{group_col}" := relevel(factor(.data[[group_col]]), ref = "f")
    )
  
  df_new
}


treatment_group_specs <- dplyr::bind_rows(
  
  make_treatment_group_spec(
    group_id   = "b10_pdsin4t1",
    group_col  = "b10_pdsin4t1",
    group_fun  = set_cd_groups,
    group_args = list(
      b_nm        = "biotic_relaxedforestnorm_5_yrs_prior_sum_yot",
      d_nm        = "pdsi_annual_5_yrs_prior_threshold_n4_yot",
      b_threshold = 10,
      d_threshold = 1
    )
  ),
  
  make_treatment_group_spec(
    group_id   = "b10_pdsin3t1",
    group_col  = "b10_pdsin3t1",
    group_fun  = set_cd_groups,
    group_args = list(
      b_nm        = "biotic_relaxedforestnorm_5_yrs_prior_sum_yot",
      d_nm        = "pdsi_annual_5_yrs_prior_threshold_n3_yot",
      b_threshold = 10,
      d_threshold = 1
    )
  ),
  
  make_treatment_group_spec(
    group_id   = "b10_pdsisumn10",
    group_col  = "b10_pdsisumn10",
    group_fun  = set_cd_groups,
    group_args = list(
      b_nm        = "biotic_relaxedforestnorm_5_yrs_prior_sum_yot",
      d_nm        = "pdsi_annual_5_yrs_prior_sum_yot",
      b_threshold = 10,
      d_threshold = -10
    )
  ),
  
  make_treatment_group_spec(
    group_id   = "b25_pdsisumn10",
    group_col  = "b25_pdsisumn10",
    group_fun  = set_cd_groups,
    group_args = list(
      b_nm        = "biotic_relaxedforestnorm_5_yrs_prior_sum_yot",
      d_nm        = "pdsi_annual_5_yrs_prior_sum_yot",
      b_threshold = 25,
      d_threshold = -10
    )
  )
  
)

treatment_group_specs <- treatment_group_specs[2,]


# ══════════════════════════════════════════════════════════════════════════════
# 6. Weighting specs
# ══════════════════════════════════════════════════════════════════════════════
# Defines propensity-score weighting approaches run BEFORE estimation.
# The weighting pipeline crosses analysis_subset_specs x treatment_group_specs
# x weighting_specs, producing one weights parquet per combination.
#
# Weight column names are constructed automatically:
#   "{weighting_id}_weights"               when weighting_name is NA (default)
#   "{weighting_id}_{weighting_name}_weights"  when weighting_name is set
#
# Skip this section and leave weighting_specs undefined (or set to NULL) if you
# are running an unweighted analysis only. In that case also set
# weights_col = NA in all model_specs (the default).

weighting_specs <- dplyr::bind_rows(
  
  # GLM logistic regression, ATO estimand (average treatment in the overlap)
  make_weighting_spec(
    weighting_id   = "glm_ato",
    weight_formula = ~ aet + srtm + tpi + def + chili + nfg_factor,
    method         = "glm",
    estimand       = "ATO",
    weighting_name = "topoclimnfg"
  )
  
  # Add further weighting specs here as needed, e.g.:
  # make_weighting_spec(
  #   weighting_id   = "gbm_ate",
  #   weight_formula = ~ biotic_relaxedforestnorm_5_yrs_prior_sum_yot
  #                      + pdsi_annual_5_yrs_prior_sum_yot
  #                      + aet + srtm + tpi,
  #   method         = "gbm",
  #   estimand       = "ATE",
  #   weighting_name = "v2"    # → column "gbm_ate_v2_weights"
  # )
  
)


# ══════════════════════════════════════════════════════════════════════════════
# 7. Model specs
# ══════════════════════════════════════════════════════════════════════════════
# formula_template uses {outcome} as the only glue placeholder.
# weights_col: name of the weight column to pass to feols(). Must match a
#   column produced by run_weighting_experiment() for the same subset x group
#   combination. NA_character_ (default) = unweighted.

model_specs <- dplyr::bind_rows(
  
  # ── Unweighted Sun-Abraham ─────────────────────────────────────────────────
  make_model_spec(
    model_id         = "b_pdsi_sunab_twfe",
    formula_template = paste0(
      "{outcome} ~ biotic_relaxedforestnorm + pdsi_annual",
      " + sunab(FirstTreat, year, ref.p = -6)",
      " | pt_id + year"
    ),
    estimator_type   = "sunab",
    term_pattern     = "^year::",
    weights_col      = NA_character_   # unweighted
  ),
  
  # ── Weighted Sun-Abraham (GLM ATO) ────────────────────────────────────────
  make_model_spec(
    model_id         = "b_pdsi_sunab_twfe_glmatotopoclimnfg",
    formula_template = paste0(
      "{outcome} ~ biotic_relaxedforestnorm + pdsi_annual",
      " + sunab(FirstTreat, year, ref.p = -6)",
      " | pt_id + year"
    ),
    estimator_type   = "sunab",
    term_pattern     = "^year::",
    weights_col      = "glm_ato_topoclimnfg_weights"
  ),
  
  # ── Plain TWFE (unweighted, for robustness) ────────────────────────────────
  # make_model_spec(
  #   model_id         = "b_pdsi_twfe",
  #   formula_template = paste0(
  #     "{outcome} ~ biotic_relaxedforestnorm + pdsi_annual",
  #     " + i(rel_year, ref = -6)",
  #     " | pt_id + year"
  #   ),
  #   estimator_type   = "twfe",
  #   term_pattern     = "^rel_year::",
  #   weights_col      = NA_character_
  # )
  
)


# ══════════════════════════════════════════════════════════════════════════════
# 8. Vcov specs
# ══════════════════════════════════════════════════════════════════════════════

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
# Optional: group palette
# ══════════════════════════════════════════════════════════════════════════════

group_palette <- c(
  "f"   = "#E69F00",   # fire only
  "bf"  = "#56B4E9",   # high biotic stress + fire
  "df"  = "#009E73",   # drought + fire
  "bdf" = "#CC79A7"    # high biotic stress + drought + fire
)


# ══════════════════════════════════════════════════════════════════════════════
# 9. Preview run grid
# ══════════════════════════════════════════════════════════════════════════════

preview_grid <- tidyr::crossing(
  analysis_subset_specs |> dplyr::select(subset_id),
  outcome_specs,
  treatment_group_specs |> dplyr::select(group_id),
  model_specs |> dplyr::select(model_id, estimator_type, weights_col)
) |>
  dplyr::mutate(
    run_id = glue::glue("{subset_id}__{outcome}__{group_id}__{model_id}")
  )

message(glue::glue("Total runs planned: {nrow(preview_grid)}"))
print(preview_grid |> dplyr::select(subset_id, outcome, group_id, model_id, weights_col))


# ══════════════════════════════════════════════════════════════════════════════
# 10. Run weighting experiment
# ══════════════════════════════════════════════════════════════════════════════
# Must run BEFORE run_experiment() for any model_specs rows with weights_col
# set to a non-NA value.
#
# Only subsets with a non-NULL short_data_source are included here. Manual
# subsets without short_data_source are not passed to this call.
#
# Comment out this entire block if running an unweighted-only analysis.

# weighting_subsets <- dplyr::bind_rows(
#   ecoregion_subset_specs,
#   all_eco_subset_spec
# )


weighting_subsets <- analysis_subset_specs

results_weighting <- run_weighting_experiment(
  dataset_spec          = dataset_spec,
  analysis_subset_specs = weighting_subsets,
  treatment_group_specs = treatment_group_specs,
  weighting_specs       = weighting_specs,
  dir_out               = dir_results,
  skip_existing         = TRUE,
  verbose_timing        = TRUE,
  .progress             = TRUE
)

failed_weighting <- purrr::keep(results_weighting$run_results, \(r) !is.null(r$error))
if (length(failed_weighting) > 0) {
  message("Failed weighting runs:")
  purrr::walk(failed_weighting, \(r) message("  ", r$weight_run_id, ": ", r$error))
}

# Rebuild merged weight registry after weighting is complete:
rebuild_weighting_tables(dir_out = dir_results, write_csv = TRUE)


# ══════════════════════════════════════════════════════════════════════════════
# 11. Run estimation experiment
# ══════════════════════════════════════════════════════════════════════════════
# For weighted model specs, resolve_weights_parquet() (called internally)
# locates the correct weights parquet from weights/by_run using subset_id x
# group_id x weights_col. It errors informatively if the weighting step was
# not run first for a given combination.
#
# Split into multiple run_experiment() calls if needed (e.g. different subsets
# on different days). All calls write into the same dir_results.

results_sunab <- run_experiment(
  dataset_spec          = dataset_spec,
  analysis_subset_specs = analysis_subset_specs,
  outcome_specs         = outcome_specs,
  treatment_group_specs = treatment_group_specs,
  model_specs           = model_specs |> dplyr::filter(estimator_type == "sunab"),
  vcov_specs            = vcov_specs,
  dir_out               = dir_results,
  group_palette         = group_palette,
  ci_level              = 0.95,
  run_estimation        = TRUE,
  run_descriptive       = FALSE,
  descriptive_args      = list(
    treated_year_var = "burn_year",
    control_year_var = "mock_burn_year"
  ),
  skip_existing         = TRUE,
  write_vcov            = FALSE,
  verbose_timing        = TRUE,
  .progress             = TRUE
)

failed_sunab <- purrr::keep(results_sunab$run_results, \(r) !is.null(r$error))
if (length(failed_sunab) > 0) {
  message("Failed estimation runs:")
  purrr::walk(failed_sunab, \(r) message("  ", r$run_id, ": ", r$error))
}

# ── Optional second call — burn-year window subsets, both estimators ──────────
# (unweighted because manual_subset_specs have no short_data_source)

# results_burnyear <- run_experiment(
#   dataset_spec          = dataset_spec,
#   analysis_subset_specs = manual_subset_specs,
#   outcome_specs         = outcome_specs,
#   treatment_group_specs = treatment_group_specs,
#   model_specs           = model_specs,   # runs all models in model_specs
#   vcov_specs            = vcov_specs,
#   dir_out               = dir_results,
#   group_palette         = group_palette,
#   ci_level              = 0.95,
#   run_estimation        = TRUE,
#   run_descriptive       = FALSE,
#   skip_existing         = TRUE,
#   write_vcov            = FALSE,
#   verbose_timing        = TRUE,
#   .progress             = TRUE
# )


# ══════════════════════════════════════════════════════════════════════════════
# 12. Rebuild merged output tables
# ══════════════════════════════════════════════════════════════════════════════
# Call at any point after one or more run_experiment() calls. Scans all
# by_run folders under dir_results, deduplicates on run_id (most recent file),
# writes merged tables to tables/all/ and descriptive/all/.

all_estimation_tables <- rebuild_estimation_tables(
  dir_out   = dir_results,
  write_csv = TRUE
)

all_descriptive_tables <- rebuild_descriptive_tables(
  dir_out   = dir_results,
  write_csv = TRUE
)

message(glue::glue(
  "Estimation rows:   {nrow(all_estimation_tables$coef_expanded_vcov)}\n",
  "Unique run_ids:    {dplyr::n_distinct(all_estimation_tables$run_registry$run_id)}\n",
  "Descriptive rows:  {nrow(all_descriptive_tables$event_time_trajectory)}"
))


