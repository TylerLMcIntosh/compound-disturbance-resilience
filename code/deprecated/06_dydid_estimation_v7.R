# 06_dydid_estimation_v7.R
# Execution script: weighting and DiD coefficient estimation
# -------------------------------------------------------
# Workflow:
#   1.  Environment setup
#   2.  Dataset spec
#   3.  Analysis subset specs
#   4.  Outcome specs
#   5.  Treatment group specs  (group_fun produces group_col + all dummy_cols)
#   6.  Weighting specs
#   7.  Model specs            (formula_template contains full dummy interaction)
#   8.  Vcov specs
#   9.  Preview run grid
#   10. Run weighting experiment
#   11. Run estimation experiment
#   12. Rebuild merged tables
# -------------------------------------------------------

Sys.setenv(LD_LIBRARY_PATH = paste("/opt/conda/lib", Sys.getenv("LD_LIBRARY_PATH"), sep = ":"))
Sys.setenv(PATH = paste("/usr/bin:/bin:/usr/local/bin", Sys.getenv("PATH"), sep = ":"))
Sys.setenv(PKG_CONFIG_PATH = "/usr/lib/x86_64-linux-gnu/pkgconfig")


# ══════════════════════════════════════════════════════════════════════════════
# 1. Environment setup ----
# ══════════════════════════════════════════════════════════════════════════════

rm(list = ls())

if (!requireNamespace("here", quietly = TRUE)) install.packages("here")
library(here)
here::i_am("code/06_dydid_estimation_v7.R")

required_pkgs <- c(
  "dplyr", "ggplot2", "tidyr", "readr", "purrr", "tibble", "stringr",
  "forcats", "fixest", "arrow", "glue", "here", "WeightIt"
)
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) install.packages(missing_pkgs)

library(dplyr); library(ggplot2); library(tidyr); library(readr)
library(tibble); library(purrr);  library(stringr); library(forcats)
library(fixest); library(arrow);  library(glue);    library(WeightIt)

source(here::here("code", "weight_dydid_pipeline_v7.R"))

seed <- 1234
set.seed(seed)


# ── Directory layout ──────────────────────────────────────────────────────────

run_name <- "GEE_resilience_v7_operational_ss500_ts50000"
version  <- "v7"
cyverse  <- FALSE

if (cyverse) {
  dir_base    <- file.path("~/data-store/data/iplant/home/shared/earthlab/macrosystems/tlm", run_name)
  dir_data    <- file.path(dir_base, "data", "derived")
  dir_raw     <- file.path(dir_base, "data", "raw")
  dir_manual  <- file.path(dir_base, "data", "manual")
  dir_results <- file.path(dir_base, "results", version)
  dir_figs    <- file.path(dir_base, "figs",    version)
} else {
  dir_data    <- here::here("data", "derived", run_name)
  dir_raw     <- here::here("data", "raw")
  dir_manual  <- here::here("data", "manual")
  dir_results <- here::here("results", version)
  dir_figs    <- here::here("figs",    version)
}

dir_parquet_long  <- file.path(dir_data, "parquet_long_filtered")
dir_parquet_short <- file.path(dir_data, "parquet_short_filtered")

dir_ensure_local(c(dir_data, dir_parquet_long, dir_raw, dir_manual, dir_results, dir_figs))


# ══════════════════════════════════════════════════════════════════════════════
# 2. Dataset spec ----
# ══════════════════════════════════════════════════════════════════════════════

dataset_spec <- make_dataset_spec(
  unit_id    = "pt_id",
  time_var   = "year",
  trt_col    = "fire",
  cohort_var = "FirstTreat",
  event_id   = "fireid"
)


# ══════════════════════════════════════════════════════════════════════════════
# 3. Analysis subset specs ----
# ══════════════════════════════════════════════════════════════════════════════

ecoregion_code_names <- c(
  "bluemtns", "cascades", "coastrange", "eastcascades",
  "klamathmtns", "northcascades", "pugetlowland",
  "willamettevalley", "centralcaliforniamtns", "sierranevada",
  "southerncaliforniamtns", "canadianrockies", "idahobatholith",
  "middlerockies", "northernrockies", "southernrockies",
  "wasatchuintamtns", "aznmmtns", "coloradoplateaus"
)

ecoregion_subset_specs <- expand_analysis_subset_specs_by_col(
  long_data_source  = dir_parquet_long,
  split_col         = "ecoregion_code_name",
  id_prefix         = "ecor",
  values            = ecoregion_code_names,
  short_data_source = dir_parquet_short
)

forestgroup_subset_specs <- expand_analysis_subset_specs_by_col(
  long_data_source  = dir_parquet_long,
  split_col         = "nfg_factor",
  id_prefix         = "nfg",
  short_data_source = dir_parquet_short,
  check_all_files   = TRUE
)

all_data_subset_spec <- make_analysis_subset_spec(
  subset_id         = "all_ecoregions",
  long_data_source  = dir_parquet_long,
  short_data_source = dir_parquet_short
)

all_data_temporalsplit_subset_specs <- dplyr::bind_rows(
  make_analysis_subset_spec(
    subset_id         = "burnyear_2000_2009",
    long_data_source  = dir_parquet_long,
    data_filter       = ~ burn_year >= 2000 & burn_year < 2010,
    short_data_source = dir_parquet_short
  ),
  make_analysis_subset_spec(
    subset_id         = "burnyear_2010_2019",
    long_data_source  = dir_parquet_long,
    data_filter       = ~ burn_year >= 2010 & burn_year < 2020,
    short_data_source = dir_parquet_short
  )
)


# ══════════════════════════════════════════════════════════════════════════════
# 4. Outcome specs ----
# ══════════════════════════════════════════════════════════════════════════════

outcome_specs <- tibble::tibble(outcome = c("rap_tree", "vcf_tree"))


# ══════════════════════════════════════════════════════════════════════════════
# 5. Treatment group specs ----
# ══════════════════════════════════════════════════════════════════════════════
#
# group_fun (set_cd_groups) must produce:
#   - group_col: multi-valued factor (f/bf/df/bdf/control) used by WeightIt
#   - dummy_cols: binary integer columns (cd_f, cd_bf, cd_df, cd_bdf) for feols
#
# The four categories are designed to be exhaustive for all treated units
# (fire == 1 will always be assigned to exactly one of f/bf/df/bdf).
# Units where fire == 1 and all dummies == 0 after replacement indicate a data
# issue and will be dropped with a warning by run_one_estimation().

set_cd_groups <- function(df,
                          group_col,
                          b_nm,
                          d_nm,
                          b_threshold,
                          d_threshold,
                          dummy_cols      = c("cd_f", "cd_bf", "cd_df", "cd_bdf"),
                          include_control = FALSE) {
  
  # multi-valued group assignment for WeightIt
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
  
  # "f" as reference level for WeightIt factor
  df_new <- df_new |>
    dplyr::mutate("{group_col}" := relevel(factor(.data[[group_col]]), ref = "f"))
  
  # binary dummies for feols; NAs become 0 (control units get 0 for all dummies)
  level_names <- c("f", "bf", "df", "bdf")
  for (i in seq_along(dummy_cols)) {
    lv <- level_names[i]
    df_new[[dummy_cols[i]]] <- tidyr::replace_na(
      as.integer(!is.na(df_new[[group_col]]) & as.character(df_new[[group_col]]) == lv),
      0L
    )
  }
  
  df_new
}


initial_treatment_group_specs <- dplyr::bind_rows(
  make_treatment_group_spec(
    group_id   = "b10_pdsisumn10",
    group_col  = "b10_pdsisumn10",
    dummy_cols = c("cd_f", "cd_bf", "cd_df", "cd_bdf"),
    group_fun  = set_cd_groups,
    group_args = list(
      b_nm        = "biotic_relaxedforestnorm_5_yrs_prior_sum_yot",
      d_nm        = "pdsi_annual_5_yrs_prior_sum_yot",
      b_threshold = 10,
      d_threshold = -10
    )
  )
)

expanded_treatment_group_specs <- dplyr::bind_rows(
  make_treatment_group_spec(
    group_id   = "b10_pdsin4t1",
    group_col  = "b10_pdsin4t1",
    dummy_cols = c("cd_f", "cd_bf", "cd_df", "cd_bdf"),
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
    dummy_cols = c("cd_f", "cd_bf", "cd_df", "cd_bdf"),
    group_fun  = set_cd_groups,
    group_args = list(
      b_nm        = "biotic_relaxedforestnorm_5_yrs_prior_sum_yot",
      d_nm        = "pdsi_annual_5_yrs_prior_threshold_n3_yot",
      b_threshold = 10,
      d_threshold = 1
    )
  ),
  make_treatment_group_spec(
    group_id   = "b25_pdsisumn10",
    group_col  = "b25_pdsisumn10",
    dummy_cols = c("cd_f", "cd_bf", "cd_df", "cd_bdf"),
    group_fun  = set_cd_groups,
    group_args = list(
      b_nm        = "biotic_relaxedforestnorm_5_yrs_prior_sum_yot",
      d_nm        = "pdsi_annual_5_yrs_prior_sum_yot",
      b_threshold = 25,
      d_threshold = -10
    )
  )
)


# ══════════════════════════════════════════════════════════════════════════════
# 6. Weighting specs ----
# ══════════════════════════════════════════════════════════════════════════════

weighting_specs <- dplyr::bind_rows(
  make_weighting_spec(
    weighting_id   = "glm_ato_topoclimnfg",
    weight_formula = ~ aet + srtm + tpi + def + chili + nfg_factor,
    method         = "glm",
    estimand       = "ATO"
  ),
  make_weighting_spec(
    weighting_id   = "glm_ato_topoclimnfgrap",
    weight_formula = ~ aet + srtm + tpi + def + chili + nfg_factor + gam_rap_tree_pre6_fit,
    method         = "glm",
    estimand       = "ATO"
  )
)


# ══════════════════════════════════════════════════════════════════════════════
# 7. Model specs ----
# ══════════════════════════════════════════════════════════════════════════════
#
# formula_template contains the full unified dummy interaction structure.
# no_agg = TRUE gives cohort-specific coefficients as the primary output,
# which are needed for cohort-group aggregation in Script 3.
#
# primary_vcov_id identifies which vcov spec is stored in the lean model
# and used for cohort-group aggregation in Script 3.
#
# mem.clean = TRUE is strongly recommended for models of this size.

# ── Unweighted Sun-Abraham (primary model) ─────────────────────────────────

sunab_formula_b10_pdsisumn10 <- paste0(
  "{outcome} ~ ",
  "sunab(FirstTreat, year, ref.p = -6, no_agg = TRUE):cd_f + ",
  "sunab(FirstTreat, year, ref.p = -6, no_agg = TRUE):cd_bf + ",
  "sunab(FirstTreat, year, ref.p = -6, no_agg = TRUE):cd_df + ",
  "sunab(FirstTreat, year, ref.p = -6, no_agg = TRUE):cd_bdf",
  " | pt_id + year"
)

initial_model_specs <- dplyr::bind_rows(
  
  make_model_spec(
    model_id         = "sunab_twfe_unweighted",
    formula_template = sunab_formula_b10_pdsisumn10,
    estimator_type   = "sunab",
    term_pattern     = "^year::",
    weights_col      = NA_character_,
    primary_vcov_id  = "cluster_pt",
    feols_args       = list(mem.clean = TRUE)
  ),
  
  make_model_spec(
    model_id         = "sunab_twfe_glmatotopoclimnfg",
    formula_template = sunab_formula_b10_pdsisumn10,
    estimator_type   = "sunab",
    term_pattern     = "^year::",
    weights_col      = "glm_ato_topoclimnfg_weights",
    primary_vcov_id  = "cluster_pt",
    feols_args       = list(mem.clean = TRUE)
  )
  
)


# ══════════════════════════════════════════════════════════════════════════════
# 8. Vcov specs ----
# ══════════════════════════════════════════════════════════════════════════════

ecor_vcov_specs <- tibble::tibble(
  vcov_id = c(
    "cluster_pt",
    "cluster_eco4",
    "conley_125km_20km"
  ),
  vcov_label = c(
    "Clustered SEs: pt_id",
    "Clustered SEs: us_l4code",
    "Conley SEs: 125 km cutoff, 20 km pixel"
  ),
  vcov = list(
    stats::as.formula("~ pt_id"),
    stats::as.formula("~ us_l4code"),
    fixest::vcov_conley(lat = "lat", lon = "long", cutoff = 125, pixel = 20, distance = "triangular")
  ),
  vcov_vars = list("pt_id", "us_l4code", c("lat", "long"))
)

regionwide_vcov_specs <- tibble::tibble(
  vcov_id = c(
    "cluster_pt",
    "cluster_eco3",
    "conley_125km_20km"
  ),
  vcov_label = c(
    "Clustered SEs: pt_id",
    "Clustered SEs: us_l3code",
    "Conley SEs: 125 km cutoff, 20 km pixel"
  ),
  vcov = list(
    stats::as.formula("~ pt_id"),
    stats::as.formula("~ us_l3code"),
    fixest::vcov_conley(lat = "lat", lon = "long", cutoff = 125, pixel = 20, distance = "triangular")
  ),
  vcov_vars = list("pt_id", "us_l3code", c("lat", "long"))
)


# ══════════════════════════════════════════════════════════════════════════════
# Optional: group palette (keyed on dummy column names)
# ══════════════════════════════════════════════════════════════════════════════

group_palette <- c(
  "cd_f"   = "#E69F00",
  "cd_bf"  = "#56B4E9",
  "cd_df"  = "#009E73",
  "cd_bdf" = "#CC79A7"
)


# ══════════════════════════════════════════════════════════════════════════════
# 9. Preview run grid ----
# ══════════════════════════════════════════════════════════════════════════════

preview_run_grid <- function(subset_specs, outcome_specs, treatment_group_specs,
                             model_specs, vcov_specs) {
  grid <- tidyr::crossing(
    subset_specs |> dplyr::select(subset_id),
    outcome_specs,
    treatment_group_specs |> dplyr::select(group_id, dummy_cols),
    model_specs |> dplyr::select(model_id, estimator_type, weights_col, primary_vcov_id)
  ) |>
    dplyr::mutate(run_id = glue::glue("{subset_id}__{outcome}__{group_id}__{model_id}"),
                  dummy_cols = purrr::map_chr(dummy_cols, \(dc) paste(dc, collapse = ",")))
  
  message(glue::glue("Total runs planned: {nrow(grid)}"))
  message(glue::glue("Vcov specs per run: {paste(vcov_specs$vcov_id, collapse = ', ')}"))
  print(grid |> dplyr::select(subset_id, outcome, group_id, model_id, weights_col,
                              primary_vcov_id, dummy_cols))
}

preview_run_grid(ecoregion_subset_specs, outcome_specs, initial_treatment_group_specs,
                 initial_model_specs, ecor_vcov_specs)

preview_run_grid(all_data_subset_spec, outcome_specs, initial_treatment_group_specs,
                 initial_model_specs, regionwide_vcov_specs)


# ══════════════════════════════════════════════════════════════════════════════
# 10. Run weighting experiment ----
# ══════════════════════════════════════════════════════════════════════════════

weighting_subsets <- dplyr::bind_rows(
  ecoregion_subset_specs,
  all_data_temporalsplit_subset_specs
)

results_weighting <- run_weighting_experiment(
  dataset_spec          = dataset_spec,
  analysis_subset_specs = weighting_subsets,
  treatment_group_specs = initial_treatment_group_specs,
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

rebuild_weighting_tables(dir_out = dir_results, write_csv = TRUE)


# ══════════════════════════════════════════════════════════════════════════════
# 11. Run estimation experiment ----
# ══════════════════════════════════════════════════════════════════════════════

results_sunab_ecor <- run_experiment(
  dataset_spec          = dataset_spec,
  analysis_subset_specs = ecoregion_subset_specs,
  outcome_specs         = outcome_specs,
  treatment_group_specs = initial_treatment_group_specs,
  model_specs           = initial_model_specs,
  vcov_specs            = ecor_vcov_specs,
  dir_out               = dir_results,
  group_palette         = group_palette,
  ci_level              = 0.95,
  run_estimation        = TRUE,
  run_descriptive       = FALSE,
  descriptive_args      = list(treated_year_var = "burn_year", control_year_var = "mock_burn_year"),
  skip_existing         = TRUE,
  verbose_timing        = TRUE,
  .progress             = TRUE
)

failed_ecor <- purrr::keep(results_sunab_ecor$run_results, \(r) !is.null(r$error))
if (length(failed_ecor) > 0) {
  message("Failed estimation runs:")
  purrr::walk(failed_ecor, \(r) message("  ", r$run_id, ": ", r$error))
}


# ── Optional: pooled all-ecoregion run ────────────────────────────────────────

# results_sunab_all <- run_experiment(
#   dataset_spec          = dataset_spec,
#   analysis_subset_specs = all_data_subset_spec,
#   outcome_specs         = outcome_specs,
#   treatment_group_specs = initial_treatment_group_specs,
#   model_specs           = initial_model_specs,
#   vcov_specs            = regionwide_vcov_specs,
#   dir_out               = dir_results,
#   group_palette         = group_palette,
#   ci_level              = 0.95,
#   run_estimation        = TRUE,
#   run_descriptive       = FALSE,
#   skip_existing         = TRUE,
#   verbose_timing        = TRUE,
#   .progress             = TRUE
# )


# ══════════════════════════════════════════════════════════════════════════════
# 12. Rebuild merged output tables ----
# ══════════════════════════════════════════════════════════════════════════════

all_estimation_tables <- rebuild_estimation_tables(dir_out = dir_results, write_csv = TRUE)
all_descriptive_tables <- rebuild_descriptive_tables(dir_out = dir_results, write_csv = TRUE)

message(glue::glue(
  "Coef rows:            {nrow(all_estimation_tables$coef_tbl)}\n",
  "Agg event-study rows: {nrow(all_estimation_tables$agg_eventstudy)}\n",
  "Unique run_ids:       {dplyr::n_distinct(all_estimation_tables$run_registry$run_id)}"
))