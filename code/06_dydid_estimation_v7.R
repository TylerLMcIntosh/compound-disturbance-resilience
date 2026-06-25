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
#   7.  Model specs
#   8.  Vcov specs
#   9.  Aggregation specs      (agg_spec x vcov_spec run inside estimation)
#   10. Preview run grid
#   11. Run weighting experiment
#   12. Run estimation experiment
#   13. Rebuild merged tables
# -------------------------------------------------------


# Monitoring resources on Exosphere
# 1TB RAM instance:
# Highest RAM usage for all runs was ~25% (250GB)
# EXCEPT for all-data run - which spiked to 1TB and crashed...
# appeared to cruise for a while at 850GB
# Full dataset: 14.27 million pixel-time units, 425,743 pt_ids; 1993-2025 (~30 years)


# Sys.setenv(LD_LIBRARY_PATH = paste("/opt/conda/lib", Sys.getenv("LD_LIBRARY_PATH"), sep = ":"))
# Sys.setenv(PATH = paste("/usr/bin:/bin:/usr/local/bin", Sys.getenv("PATH"), sep = ":"))
# Sys.setenv(PKG_CONFIG_PATH = "/usr/lib/x86_64-linux-gnu/pkgconfig")
# 

# ══════════════════════════════════════════════════════════════════════════════
# 1. Environment setup ----
# ══════════════════════════════════════════════════════════════════════════════

rm(list = ls())


if (!requireNamespace("here", quietly = TRUE)) install.packages("here")
library(here)
here::i_am("code/06_dydid_estimation_v7.R")

required_pkgs <- c(
  "dplyr", "ggplot2", "tidyr", "readr", "purrr", "tibble", "stringr",
  "forcats", "fixest", "arrow", "glue", "here", "WeightIt", "tictoc"
)
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) install.packages(missing_pkgs)

library(dplyr); library(ggplot2); library(tidyr);  library(readr)
library(tibble); library(purrr);  library(stringr); library(forcats)
library(fixest); library(arrow);  library(glue);    library(WeightIt)
library(tictoc)

# sunab_aggregate_vcov must be sourced before the pipeline file
source(here::here("code", "sunab_aggregate_vcov.R"))
source(here::here("code", "weight_dydid_pipeline_v7.R"))

seed <- 1234
set.seed(seed)

# Set number of cores to use in FEOLS call
#fixest::setFixest_nthreads(48)
fixest::getFixest_nthreads()

# log any unhandled errors to the status file before R exits
options(error = function() {
  err_msg <- geterrmessage()
  line    <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                    " | FATAL | ", err_msg)
  write(line, file = here::here("logs/status_file.txt"), append = TRUE)
})

# ── Directory layout ──────────────────────────────────────────────────────────

run_name <- "GEE_resilience_v6_operational_ss500_ts50000"
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
  dir_results <- here::here("results-exo", version)
  dir_figs    <- here::here("figs",    version)
}

dir_parquet_long  <- file.path(dir_data, "parquet_long_filtered")
dir_parquet_short <- file.path(dir_data, "parquet_short_filtered")

dir_ensure_local(c(dir_data, dir_parquet_long, dir_raw, dir_manual, dir_results, dir_figs))

x <- arrow::open_dataset(dir_parquet_long) |> collect()
summary <- x |> group_by(nfg_factor) |> summarize(n = n())
length(unique(x$pt_id))

x_s <- arrow::open_dataset(dir_parquet_short) |> collect()
summary_s <- x_s |> group_by(nfg_factor, fire) |> summarize(n = n())

analysis_portion <- c(1:15)

# ══════════════════════════════════════════════════════════════════════════════
# 1.5. Test datasets ----
# ══════════════════════════════════════════════════════════════════════════════
# 
# 
# create_test_set <- function(long, short) {
# 
#   dir_parquet_long_test  <- file.path(dir_data, "parquet_long_filtered_test")
#   dir_parquet_short_test <- file.path(dir_data, "parquet_short_filtered_test")
#   dir_ensure_local(c(dir_parquet_long_test, dir_parquet_short_test))
# 
#   short_file <- file.path(dir_parquet_short_test, basename(short))
#   long_file <- file.path(dir_parquet_long_test, basename(long))
# 
#   if(!file.exists(short_file) | !file.exists(long_file)) {
# 
#     dl <- arrow::read_parquet(long)
#     ds <- arrow::read_parquet(short)
# 
#     set.seed(123)
# 
#     sample_pt_ids <- sample(unique(ds$pt_id), round(dplyr::n_distinct(ds$pt_id) * 0.1))
# 
#     dl_t <- dl |>
#       dplyr::mutate(event_time = year - FirstTreat) |>
#       dplyr::filter(
#         dplyr::between(event_time, -15, 20) | FirstTreat == 1000    # Filtering to -15-20 drops about 1.3 million points
#       ) |>
#       dplyr::filter(pt_id %in% sample_pt_ids)
# 
#     ds_t <- ds |>
#       dplyr::filter(pt_id %in% sample_pt_ids)
# 
# 
#     arrow::write_parquet(ds_t, short_file)
#     arrow::write_parquet(dl_t, long_file)
#   }
# 
# }
# 
# 
# # run to create the test dataset
# 
# long_files <- list.files(dir_parquet_long, full.names = TRUE)
# short_files <- list.files(dir_parquet_short, full.names = TRUE)
# 
# stopifnot(length(long_files) == length(short_files))
# 
# purrr::pmap(
#   .l = list(
#     long = long_files,
#     short = short_files
#   ),
#   .f = create_test_set
# )
# 
# dir_parquet_long  <- file.path(dir_data, "parquet_long_filtered_test")
# dir_parquet_short <- file.path(dir_data, "parquet_short_filtered_test")
# 
# dir_results <- file.path(dir_results, "test_10perc")
# dir_figs <- file.path(dir_results, "test_10perc")

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


conifer_nfg <- c(
  "Pinyon/Juniper Group",
  "Douglas-fir Group",
  "Ponderosa Pine Group",
  "Fir/Spruce/Mountain Hemlock Group",
  "Other Western Softwood Group",
  "Lodgepole Pine Group",
  "Western Larch Group",
  "Hemlock/Sitka Spruce Group",
  "California Mixed Conifer Group",
  "Western White Pine Group",
  "Redwood Group"
)

broadleaf_nfg <- c(
  "Western Oak Group",
  "Other Western Hardwoods Group",
  "Aspen/Birch Group",
  "Elm/Ash/Cottonwood Group",
  "Alder/Maple Group",
  "Tanoak/Laurel Group",
  "Oak/Hickory Group",
  "Tanoak/Laurel Group"
)



ecoregion_code_names <- c(
  "bluemtns", "cascades", "coastrange", "eastcascades",
  "klamathmtns", "northcascades", "pugetlowland",
  "willamettevalley", "centralcaliforniamtns", "sierranevada",
  "southerncaliforniamtns", "canadianrockies", "idahobatholith",
  "middlerockies", "northernrockies", "southernrockies",
  "wasatchuintamtns", "aznmmtns", "coloradoplateaus"
)

# ecoregion_subset_specs <- expand_analysis_subset_specs_by_col(
#   long_data_source  = dir_parquet_long,
#   split_col         = "ecoregion_code_name",
#   id_prefix         = "ecor",
#   values            = ecoregion_code_names,
#   short_data_source = dir_parquet_short
# )

# Only analyze forest groups with over 10,000 unique pt_ids
nfg_code_names <- arrow::open_dataset(dir_parquet_short) |>
  collect() |>
  group_by(nfg_factor) |>
  summarize(n = n()) |>
  arrange(desc(n)) |>
  filter(n >= 10000) |>
  pull(nfg_factor)

extended_nfg_code_names <- arrow::open_dataset(dir_parquet_short) |>
  collect() |>
  group_by(nfg_factor) |>
  summarize(n = n()) |>
  arrange(desc(n)) |>
  pull(nfg_factor)

forestgroup_subset_specs <- expand_analysis_subset_specs_by_col(
  long_data_source  = dir_parquet_long,
  split_col         = "nfg_factor",
  id_prefix         = "nfg",
  values            = nfg_code_names,
  short_data_source = dir_parquet_short,
  check_all_files   = TRUE,
  base_filter       = ~ fire == 0 | 
    (fire == 1 & year_from_fire_index >= -20 & 
       year_from_fire_index <= 20 & 
       analysis_subset %in% analysis_portion)
)

all_data_subset_spec <- make_analysis_subset_spec(
  subset_id   = "all_ecoregions",
  long_data_source  = dir_parquet_long,
  short_data_source = dir_parquet_short
)


# all_data_temporalsplit_subset_specs <- dplyr::bind_rows(
#   make_analysis_subset_spec(
#     subset_id         = "burnyear_2000_2009",
#     long_data_source  = dir_parquet_long,
#     data_filter       = ~ burn_year >= 2000 & burn_year < 2010,
#     short_data_source = dir_parquet_short
#   ),
#   make_analysis_subset_spec(
#     subset_id         = "burnyear_2010_2019",
#     long_data_source  = dir_parquet_long,
#     data_filter       = ~ burn_year >= 2010 & burn_year < 2020,
#     short_data_source = dir_parquet_short
#   )
# )


# nfg_temporalsplit_subset_specs <- dplyr::bind_rows(
#   # early burn cohort (2000-2009) x nfg group
#   expand_analysis_subset_specs_by_col(
#     long_data_source  = dir_parquet_long,
#     split_col         = "nfg_factor",
#     id_prefix         = "nfg_early",
#     values            = nfg_code_names,
#     short_data_source = dir_parquet_short,
#     check_all_files   = TRUE,
#     base_filter       = ~ burn_year >= 2000 & burn_year < 2010
#   ),
#   # late burn cohort (2010-2019) x nfg group
#   expand_analysis_subset_specs_by_col(
#     long_data_source  = dir_parquet_long,
#     split_col         = "nfg_factor",
#     id_prefix         = "nfg_late",
#     values            = nfg_code_names,
#     short_data_source = dir_parquet_short,
#     check_all_files   = TRUE,
#     base_filter       = ~ burn_year >= 2010 & burn_year < 2020
#   )
# )

# ══════════════════════════════════════════════════════════════════════════════
# 4. Outcome specs ----
# ══════════════════════════════════════════════════════════════════════════════

outcome_specs <- tibble::tibble(outcome = c("rap_tree", "vcf_tree"))


# ══════════════════════════════════════════════════════════════════════════════
# 5. Treatment group specs ----
# ══════════════════════════════════════════════════════════════════════════════
#
# set_cd_groups() produces:
#   - group_col: multi-valued factor (f/bf/df/bdf/control) for WeightIt
#   - dummy_cols: binary integer columns (cd_f, cd_bf, cd_df, cd_bdf) for feols
#
# The four treatment categories are exhaustive: every fire==1 unit belongs to
# exactly one. Units where fire==1 and all dummies==0 after grouping indicate
# a data issue and will be dropped with a warning by run_one_estimation().

set_cd_groups <- function(df,
                           group_col,
                           b_nm,
                           d_nm,
                           b_threshold,
                           d_threshold,
                           dummy_cols      = c("cd_f", "cd_bf", "cd_df", "cd_bdf"),
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

  # f as reference level for WeightIt factor
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

forest_group_weighting_specs <- dplyr::bind_rows(
  make_weighting_spec(
    weighting_id   = "glm_ato_topoclimnfg",
    weight_formula = ~ aet + srtm + tpi + def + chili,
    method         = "glm",
    estimand       = "ATO"
  ),
  make_weighting_spec(
    weighting_id   = "glm_ato_topoclimnfgrap",
    weight_formula = ~ aet + srtm + tpi + def + chili + gam_rap_tree_pre6_fit,
    method         = "glm",
    estimand       = "ATO"
  )
)


# ══════════════════════════════════════════════════════════════════════════════
# 7. Model specs ----
# ══════════════════════════════════════════════════════════════════════════════
#
# formula_template contains the full unified dummy interaction structure.
# no_agg = TRUE gives cohort-specific coefficients required by agg_specs.
# mem.clean = TRUE is recommended for models of this size.

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
    feols_args       = list(mem.clean = TRUE)
  ),

  make_model_spec(
    model_id         = "sunab_twfe_glmatotopoclimnfg",
    formula_template = sunab_formula_b10_pdsisumn10,
    estimator_type   = "sunab",
    term_pattern     = "^year::",
    weights_col      = "glm_ato_topoclimnfg_weights",
    feols_args       = list(mem.clean = TRUE)
  )

)


# ══════════════════════════════════════════════════════════════════════════════
# 8. Vcov specs ----
# ══════════════════════════════════════════════════════════════════════════════

vcov_specs <- tibble::tibble(
  vcov_id = c(
    "cluster_pt",
    "cluster_h3",
    "conley_75km_5km"
  ),
  vcov_label = c(
    "Clustered SEs: pt_id",
    "Clustered SEs: h3jsr_5",
    "Conley SEs: 75 km cutoff, 5 km pixel"
  ),
  vcov = list(
    stats::as.formula("~ pt_id"),
    stats::as.formula("~ h3jsr_5"),
    fixest::vcov_conley(lat = "lat", lon = "long", cutoff = 75, pixel = 5, distance = "triangular")
  ),
  vcov_vars = list("pt_id", "h3jsr_5", c("lat", "long"))
)



# ══════════════════════════════════════════════════════════════════════════════
# 9. Aggregation specs ----
# ══════════════════════════════════════════════════════════════════════════════
#
# Each spec runs sunab_aggregate_vcov() for every (run x vcov_spec) combination.
# group_fun must produce event_time (integer) and dummy_group (character) at
# minimum — these column names are required by Script 3 inference functions.
# Additional columns (e.g. cohort_bin) define finer aggregation strata loaded
# in Script 3 via subset_agg_by_cohort_group().
#
# Coefficient names look like: year::2:cohort::2005:cd_bf

agg_specs <- list(

  # Standard event study: cohort-averaged ATT per event_time x dummy_group.
  # Primary aggregation for pre-trend tests, ATT windows, and pairwise comparisons.
  make_agg_spec(
    id    = "event_study",
    agg   = "(year::-?[0-9]+):cohort::[0-9]+:(cd_.*)",
    group_fun = function(x) {
      x |>
        dplyr::mutate(
          event_time  = as.integer(stringr::str_extract(group_1, "-?[0-9]+")),
          dummy_group = group_2
        ) |>
        dplyr::select(term, event_time, dummy_group)
    },
    label = "Cohort-averaged event study by dummy group"
  )#,

  # # Cohort-bin event study: same estimates further stratified by cohort bin.
  # # Script 3 subsets by cohort_bin to run cohort-stratified inference.
  # # Adjust min/max cohort years to match your data.
  # make_agg_spec(
  #   id    = "cohort_early_late",
  #   agg   = "(year::-?[0-9]+):cohort::([0-9]+):(cd_.*)",
  #   group_fun = function(x) {
  #     x |>
  #       dplyr::mutate(
  #         event_time  = as.integer(stringr::str_extract(group_1, "-?[0-9]+")),
  #         cohort      = as.integer(group_2),
  #         dummy_group = group_3,
  #         cohort_bin  = dplyr::case_when(
  #           cohort >= 2000 & cohort <= 2010 ~ "early",
  #           cohort >= 2011 & cohort <= 2020 ~ "late",
  #           TRUE ~ NA_character_
  #         )
  #       ) |>
  #       dplyr::filter(!is.na(cohort_bin)) |>
  #       dplyr::select(term, event_time, cohort_bin, dummy_group)
  #   },
  #   label = "Event study by cohort bin (early 2000-2010 / late 2011-2020) and dummy group"
  # )

)


# ══════════════════════════════════════════════════════════════════════════════
# Group palette (keyed on dummy column names)
# ══════════════════════════════════════════════════════════════════════════════

group_palette <- c(
  "cd_f"   = "#E69F00",
  "cd_bf"  = "#56B4E9",
  "cd_df"  = "#009E73",
  "cd_bdf" = "#CC79A7"
)


# ══════════════════════════════════════════════════════════════════════════════
# 10. Preview run grid ----
# ══════════════════════════════════════════════════════════════════════════════

preview_run_grid <- function(subset_specs, outcome_specs, treatment_group_specs,
                              model_specs, vcov_specs, agg_specs) {
  grid <- tidyr::crossing(
    subset_specs |> dplyr::select(subset_id),
    outcome_specs,
    treatment_group_specs |> dplyr::select(group_id, dummy_cols),
    model_specs |> dplyr::select(model_id, estimator_type, weights_col)
  ) |>
    dplyr::mutate(run_id     = glue::glue("{subset_id}__{outcome}__{group_id}__{model_id}"),
                  dummy_cols = purrr::map_chr(dummy_cols, \(dc) paste(dc, collapse = ",")))

  message(glue::glue("Total runs planned: {nrow(grid)}"))
  message(glue::glue("Vcov specs per run: {paste(vcov_specs$vcov_id, collapse = ', ')}"))
  message(glue::glue("Agg specs per run:  {paste(purrr::map_chr(agg_specs, 'id'), collapse = ', ')}"))
  print(grid |> dplyr::select(subset_id, outcome, group_id, model_id,
                               weights_col, dummy_cols))
}


# SIMPLIFY

# ecoregion_subset_specs <- ecoregion_subset_specs[2,] #11 - sierra nevada; 
# outcome_specs <- outcome_specs[1,]
# vcov_specs <- vcov_specs[1,]

preview_run_grid(ecoregion_subset_specs, outcome_specs, initial_treatment_group_specs,
                 initial_model_specs, vcov_specs, agg_specs)

preview_run_grid(all_data_subset_spec, outcome_specs, initial_treatment_group_specs,
                 initial_model_specs, vcov_specs, agg_specs)


# ══════════════════════════════════════════════════════════════════════════════
# 11. Run weighting experiment ----
# ══════════════════════════════════════════════════════════════════════════════

tic('weighting ecoregions')
ecoregion_results_weighting <- run_weighting_experiment(
  dataset_spec          = dataset_spec,
  analysis_subset_specs = ecoregion_subset_specs,
  treatment_group_specs = initial_treatment_group_specs,
  weighting_specs       = weighting_specs,
  dir_out               = dir_results,
  skip_existing         = TRUE,
  verbose_timing        = TRUE,
  .progress             = TRUE
)
toc()

failed_weighting <- purrr::keep(ecoregion_results_weighting$run_results, \(r) !is.null(r$error))
if (length(failed_weighting) > 0) {
  message("Failed weighting runs:")
  purrr::walk(failed_weighting, \(r) message("  ", r$weight_run_id, ": ", r$error))
}

tic('weighting forest groups')
forestgroup_results_weighting <- run_weighting_experiment(
  dataset_spec          = dataset_spec,
  analysis_subset_specs = forestgroup_subset_specs,
  treatment_group_specs = initial_treatment_group_specs,
  weighting_specs       = forest_group_weighting_specs,
  dir_out               = dir_results,
  skip_existing         = TRUE,
  verbose_timing        = TRUE,
  .progress             = TRUE
)
toc()

failed_weighting <- purrr::keep(forestgroup_results_weighting$run_results, \(r) !is.null(r$error))
if (length(failed_weighting) > 0) {
  message("Failed weighting runs:")
  purrr::walk(failed_weighting, \(r) message("  ", r$weight_run_id, ": ", r$error))
}

tic('weighting all')
all_results_weighting <- run_weighting_experiment(
  dataset_spec          = dataset_spec,
  analysis_subset_specs = all_data_subset_spec,
  treatment_group_specs = initial_treatment_group_specs,
  weighting_specs       = weighting_specs,
  dir_out               = dir_results,
  skip_existing         = TRUE,
  verbose_timing        = TRUE,
  .progress             = TRUE
)
toc()

failed_weighting <- purrr::keep(all_results_weighting$run_results, \(r) !is.null(r$error))
if (length(failed_weighting) > 0) {
  message("Failed weighting runs:")
  purrr::walk(failed_weighting, \(r) message("  ", r$weight_run_id, ": ", r$error))
}

rebuild_weighting_tables(dir_out = dir_results, write_csv = TRUE)


# ══════════════════════════════════════════════════════════════════════════════
# 12. Run estimation experiment ----
# ══════════════════════════════════════════════════════════════════════════════

tic('sunab estimation ecoregions')
results_sunab_ecor <- run_experiment(
  dataset_spec          = dataset_spec,
  analysis_subset_specs = ecoregion_subset_specs,
  outcome_specs         = outcome_specs,
  treatment_group_specs = initial_treatment_group_specs,
  model_specs           = initial_model_specs,
  vcov_specs            = vcov_specs,
  agg_specs             = agg_specs,
  dir_out               = dir_results,
  group_palette         = group_palette,
  ci_level              = 0.95,
  run_estimation        = TRUE,
  run_descriptive       = TRUE,
  descriptive_args      = list(treated_year_var = "burn_year", control_year_var = "mock_burn_year"),
  skip_existing         = TRUE,
  verbose_timing        = TRUE,
  .progress             = TRUE
)
toc()

failed_ecor <- purrr::keep(results_sunab_ecor$run_results, \(r) !is.null(r$error))
if (length(failed_ecor) > 0) {
  message("Failed estimation runs:")
  purrr::walk(failed_ecor, \(r) message("  ", r$run_id, ": ", r$error))
}



tic('sunab estimation forest groups')
results_sunab_forestgroup <- run_experiment(
  dataset_spec          = dataset_spec,
  analysis_subset_specs = forestgroup_subset_specs,
  outcome_specs         = outcome_specs,
  treatment_group_specs = initial_treatment_group_specs,
  model_specs           = initial_model_specs,
  vcov_specs            = vcov_specs,
  agg_specs             = agg_specs,
  dir_out               = dir_results,
  group_palette         = group_palette,
  ci_level              = 0.95,
  run_estimation        = TRUE,
  run_descriptive       = TRUE,
  descriptive_args      = list(treated_year_var = "burn_year", control_year_var = "mock_burn_year"),
  skip_existing         = TRUE,
  verbose_timing        = TRUE,
  .progress             = TRUE
)
toc()

failed_forestgroup <- purrr::keep(results_sunab_forestgroup$run_results, \(r) !is.null(r$error))
if (length(failed_forestgroup) > 0) {
  message("Failed estimation runs:")
  purrr::walk(failed_forestgroup, \(r) message("  ", r$run_id, ": ", r$error))
}



tic('sunab estimation all')
results_sunab_all <- run_experiment(
  dataset_spec          = dataset_spec,
  analysis_subset_specs = all_data_subset_spec,
  outcome_specs         = outcome_specs,
  treatment_group_specs = initial_treatment_group_specs,
  model_specs           = initial_model_specs,
  vcov_specs            = vcov_specs,
  agg_specs             = agg_specs,
  dir_out               = dir_results,
  group_palette         = group_palette,
  ci_level              = 0.95,
  run_estimation        = TRUE,
  run_descriptive       = TRUE,
  descriptive_args      = list(treated_year_var = "burn_year", control_year_var = "mock_burn_year"),
  skip_existing         = TRUE,
  verbose_timing        = TRUE,
  .progress             = TRUE
)
toc()

failed_all <- purrr::keep(results_sunab_all$run_results, \(r) !is.null(r$error))
if (length(failed_all) > 0) {
  message("Failed estimation runs:")
  purrr::walk(failed_all, \(r) message("  ", r$run_id, ": ", r$error))
}




# ══════════════════════════════════════════════════════════════════════════════
# 13. Rebuild merged output tables ----
# ══════════════════════════════════════════════════════════════════════════════

all_estimation_tables <- rebuild_estimation_tables(dir_out = dir_results, write_csv = TRUE)
all_descriptive_tables <- rebuild_descriptive_tables(dir_out = dir_results, write_csv = TRUE)

message(glue::glue(
  "Coef rows:       {nrow(all_estimation_tables$coef_tbl)}\n",
  "Agg specs found: {paste(names(all_estimation_tables$agg_tbls), collapse = ', ')}\n",
  "Unique run_ids:  {dplyr::n_distinct(all_estimation_tables$run_registry$run_id)}"
))

relativize_result_paths(dir_results = dir_results, base = here::here())

