
# 06_dydid_estimation_v8.R
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


# v8: started run 11:34pm on 250gb machine 
# had crash next day at 4:22pm on fir-spruce-hemlock 6-yr run. ~ 12 hours for full 3-yr run

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
here::i_am("code/06_dydid_estimation_v8.R")

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
fixest::setFixest_nthreads(48)
#fixest::getFixest_nthreads()
ram_size <- 250 # either 250 or 1000

# log any unhandled errors to the status file before R exits
options(error = function() {
  err_msg <- geterrmessage()
  line    <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                    " | FATAL | ", err_msg)
  write(line, file = here::here("logs/status_file.txt"), append = TRUE)
})

# ── Directory layout ──────────────────────────────────────────────────────────

run_name <- "GEE_resilience_v6_operational_ss500_ts50000"
version  <- "v8"
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

dir_parquet_long  <- file.path(dir_data, "parquet_long")
dir_parquet_short <- file.path(dir_data, "parquet_short")

dir_ensure_local(c(dir_data, dir_parquet_long, dir_raw, dir_manual, dir_results, dir_figs))

#x <- arrow::open_dataset(dir_parquet_long) |> collect()
#xx <- x[1:300,]
# summary <- x |> group_by(nfg_factor_clean) |> summarize(n = n())
# length(unique(x$pt_id))
# 
#x_s <- arrow::open_dataset(dir_parquet_short) |> collect()
# summary_s <- x_s |> group_by(nfg_factor_clean, fire) |> summarize(n = n())

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
  event_id   = "fireid" #pass NA_character_ to get NAs for event_time_support instead of error if needed
)

extended_dataset_spec <- make_dataset_spec(
  unit_id    = "pt_id",
  time_var   = "year",
  trt_col    = "extended_treat",
  cohort_var = "FirstTreat",
  event_id   = NA_character_ #pass NA_character_ to get NAs for event_time_support instead of error if needed
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
  group_by(nfg_factor_clean) |>
  summarize(n = n()) |>
  arrange(desc(n)) |>
  filter(n >= 10000) |>
  pull(nfg_factor_clean)

# extended_nfg_code_names <- arrow::open_dataset(dir_parquet_short) |>
#   collect() |>
#   group_by(nfg_factor_clean) |>
#   summarize(n = n()) |>
#   arrange(desc(n)) |>
#   pull(nfg_factor_clean)


nfg_subset_specs <- expand_analysis_subset_specs_by_col(
  long_data_source  = dir_parquet_long,
  split_col         = "nfg_factor_clean",
  id_prefix         = "nfg",
  values            = nfg_code_names[nfg_code_names != "fir_spruce_mountain_hemlock" & nfg_code_names != "douglas_fir"],
  short_data_source = dir_parquet_short,
  check_all_files   = TRUE,
  base_filter       = ~ fire == 0 | 
    (fire == 1 & #year_from_fire_index >= -20 & 
       #year_from_fire_index <= 20 & 
       analysis_subset %in% analysis_portion)
)


nfg_large_subset_specs <- expand_analysis_subset_specs_by_col(
  long_data_source  = dir_parquet_long,
  split_col         = "nfg_factor_clean",
  id_prefix         = "nfg",
  values            = c("fir_spruce_mountain_hemlock", "douglas_fir"),
  short_data_source = dir_parquet_short,
  check_all_files   = TRUE,
  base_filter       = ~ fire == 0 | 
    (fire == 1 & #year_from_fire_index >= -20 & 
       #year_from_fire_index <= 20 & 
       analysis_subset %in% analysis_portion)
)

broadtype_subset_specs <- expand_analysis_subset_specs_by_col(
  long_data_source  = dir_parquet_long,
  split_col         = "nfg_broad_type",
  id_prefix         = "broad",
  short_data_source = dir_parquet_short,
  check_all_files   = TRUE,
  base_filter       = ~ fire == 0 | 
    (fire == 1 & #year_from_fire_index >= -20 & 
       #year_from_fire_index <= 20 & 
       analysis_subset %in% analysis_portion)
)

# 
# broadtype_temporalsplit_subset_specs <- dplyr::bind_rows(
#   # early burn cohort x nfg group
#   expand_analysis_subset_specs_by_col(
#     long_data_source  = dir_parquet_long,
#     split_col         = "nfg_broad_type",
#     id_prefix         = "broad_early",
#     short_data_source = dir_parquet_short,
#     check_all_files   = TRUE,
#     base_filter = ~ fire == 0 | 
#       (fire == 1 & #year_from_fire_index >= -20 & 
#          #year_from_fire_index <= 20 & 
#          analysis_subset %in% analysis_portion & burn_year >= 2002 & burn_year <= 2011)
#   ),
#   # late burn cohort x nfg group
#   expand_analysis_subset_specs_by_col(
#     long_data_source  = dir_parquet_long,
#     split_col         = "nfg_broad_type",
#     id_prefix         = "broad_late",
#     short_data_source = dir_parquet_short,
#     check_all_files   = TRUE,
#     base_filter = ~ fire == 0 | 
#       (fire == 1 & #year_from_fire_index >= -20 & 
#          #year_from_fire_index <= 20 & 
#          analysis_subset %in% analysis_portion & burn_year >= 2012 & burn_year <= 2021)
#   )
# )
# 
# 
# nfg_temporalsplit_subset_specs <- dplyr::bind_rows(
#   # early burn cohort x nfg group
#   expand_analysis_subset_specs_by_col(
#     long_data_source  = dir_parquet_long,
#     split_col         = "nfg_factor_clean",
#     id_prefix         = "nfg_early",
#     values            = nfg_code_names,
#     short_data_source = dir_parquet_short,
#     check_all_files   = TRUE,
#     base_filter = ~ fire == 0 | 
#       (fire == 1 & #year_from_fire_index >= -20 & 
#          #year_from_fire_index <= 20 & 
#          analysis_subset %in% analysis_portion & burn_year >= 2002 & burn_year <= 2011)
#   ),
#   # late burn cohort x nfg group
#   expand_analysis_subset_specs_by_col(
#     long_data_source  = dir_parquet_long,
#     split_col         = "nfg_factor_clean",
#     id_prefix         = "nfg_late",
#     values            = nfg_code_names,
#     short_data_source = dir_parquet_short,
#     check_all_files   = TRUE,
#     base_filter = ~ fire == 0 | 
#       (fire == 1 & #year_from_fire_index >= -20 & 
#          #year_from_fire_index <= 20 & 
#          analysis_subset %in% analysis_portion & burn_year >= 2012 & burn_year <= 2021)
#   )
# )


# small_machine_specs <- dplyr::bind_rows(
#   nfg_subset_specs
# )
# 
# large_machine_specs <- dplyr::bind_rows(
#   broadtype_subset_specs,
#   nfg_large_subset_specs
# )


# ══════════════════════════════════════════════════════════════════════════════
# 4. Outcome specs ----
# ══════════════════════════════════════════════════════════════════════════════

outcome_specs <- tibble::tibble(outcome = c("rap_tree"#,
                                            #"vcf_tree"
                                            )
                                )


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







set_cd_groups_extended_strict <- function(df,
                                          group_col,
                                          b_nm,
                                          d_nm,
                                          b_rollmax_nm,
                                          d_rollmax_nm,
                                          b_threshold,
                                          d_threshold,
                                          dummy_cols      = c("cd_f", "cd_bf", "cd_df", "cd_bdf", "cd_bd", "cd_b", "cd_d"),
                                          include_control = FALSE) {
  
  if(d_threshold < 0) {
    df_new <- df |>
      dplyr::mutate(
        "{group_col}" := dplyr::case_when(
          .data[[b_rollmax_nm]] <  b_threshold &
            .data[[b_nm]] <  b_threshold &
            fire == 1 &
            .data[[d_nm]] > d_threshold &
            .data[[d_rollmax_nm]] > d_threshold ~ "f",
          .data[[b_nm]] >=  b_threshold &
            fire == 1 &
            .data[[d_nm]] > d_threshold &
            .data[[d_rollmax_nm]] > d_threshold ~ "bf",
          .data[[b_nm]] <  b_threshold &
            .data[[b_rollmax_nm]] < b_threshold &
            fire == 1 &
            .data[[d_nm]] <= d_threshold ~ "df",
          .data[[b_nm]] >= b_threshold & 
            fire == 1 & 
            .data[[d_nm]] <= d_threshold ~ "bdf",
          .data[[b_nm]] >= b_threshold & 
            fire == 0 & 
            .data[[d_nm]] <= d_threshold ~ "bd",
          .data[[b_nm]] >= b_threshold & 
            fire == 0 & 
            .data[[d_rollmax_nm]] > d_threshold &
            .data[[d_nm]] > d_threshold ~ "b",
          .data[[b_rollmax_nm]] < b_threshold & 
            .data[[b_nm]] < b_threshold &
            fire == 0 & 
            .data[[d_nm]] <= d_threshold ~ "d",
          .data[[b_rollmax_nm]] < b_threshold & 
            .data[[b_nm]] < b_threshold &
            fire == 0 & 
            .data[[d_nm]] > d_threshold &
            .data[[d_rollmax_nm]] > d_threshold ~ ifelse(include_control, "control", NA_character_),
          TRUE ~ "DROP"
        )
      )
  } else {
    stop("please use drought threshold < 0")
    
  }
  
  df_new <- df_new |>
    dplyr::filter(is.na(.data[[group_col]]) | .data[[group_col]] != "DROP") |>
    dplyr::mutate("{group_col}" := relevel(factor(.data[[group_col]]), ref = "f"))
  
  
  level_names <- c("f", "bf", "df", "bdf", "bd", "b", "d")
  
  
  # fix the dydid columns to adapt bd, d, b for the long data; these units will as of now
  # not have FirstTreat, etc data that is necessary for the estimator
  if("year" %in% names(df_new)) { #only exists in long data
    df_new <- df_new |>
      dplyr::mutate(
        treated        = dplyr::if_else(year_from_fire_index >= 0 & 
                                          .data[[group_col]] %in% level_names, 1L, 0L),
        extended_treat = dplyr::if_else(.data[[group_col]] %in% level_names, 1L, 0L),
        FirstTreat     = dplyr::if_else(.data[[group_col]] %in% level_names, 
                                        as.integer(mock_burn_year), 1000L)
      )
  }
  
  # binary dummies for feols; NAs become 0 (control units get 0 for all dummies)
  for (i in seq_along(dummy_cols)) {
    lv <- level_names[i]
    df_new[[dummy_cols[i]]] <- tidyr::replace_na(
      as.integer(!is.na(df_new[[group_col]]) & as.character(df_new[[group_col]]) == lv),
      0L
    )
  }
  
  return(df_new)
}




# standard definitions

sixyr_treatment_group_specs <- dplyr::bind_rows(
  make_treatment_group_spec(
    group_id   = "sixyr_b90global_pdsisum90global",
    group_col  = "sixyr_b90global_pdsisum90global",
    dummy_cols = c("cd_f", "cd_bf", "cd_df", "cd_bdf"),
    group_fun  = set_cd_groups,
    group_args = list(
      b_nm        = "biotic_relaxedforestnorm_sum_n5_to_p0",
      d_nm        = "pdsi_annual_sum_n5_to_p0",
      b_threshold = 9.7,
      d_threshold = -11.8
    )
  )
)

threeyr_treatment_group_specs <- dplyr::bind_rows(
  make_treatment_group_spec(
    group_id   = "threeyr_b90global_pdsisum90global",
    group_col  = "threeyr_b90global_pdsisum90global",
    dummy_cols = c("cd_f", "cd_bf", "cd_df", "cd_bdf"),
    group_fun  = set_cd_groups,
    group_args = list(
      b_nm        = "biotic_relaxedforestnorm_sum_n2_to_p0",
      d_nm        = "pdsi_annual_sum_n2_to_p0",
      b_threshold = 4,
      d_threshold = -7.9
    )    
  )
)


# extended definitions
sixyr_treatment_group_specs_extended_strict <- dplyr::bind_rows(
  make_treatment_group_spec(
    group_id   = "sixyr_b90global_pdsisum90global_extended_strict",
    group_col  = "sixyr_b90global_pdsisum90global_extended_strict",
    dummy_cols = c("cd_f", "cd_bf", "cd_df", "cd_bdf", "cd_bd", "cd_b", "cd_d"),
    group_fun  = set_cd_groups_extended_strict ,
    group_args = list(
      b_nm        = "biotic_relaxedforestnorm_sum_n5_to_p0",
      d_nm        = "pdsi_annual_sum_n5_to_p0",
      b_threshold = 9.7,
      d_threshold = -11.8,
      b_rollmax_nm = "biotic_relaxedforestnorm_rollsum6max",
      d_rollmax_nm = "pdsi_annual_rollsum6max"
    )
  )
)

threeyr_treatment_group_specs_extended_strict <- dplyr::bind_rows(
  make_treatment_group_spec(
    group_id   = "threeyr_b90global_pdsisum90global_extended_strict",
    group_col  = "threeyr_b90global_pdsisum90global_extended_strict",
    dummy_cols = c("cd_f", "cd_bf", "cd_df", "cd_bdf", "cd_bd", "cd_b", "cd_d"),
    group_fun  = set_cd_groups_extended_strict ,
    group_args = list(
      b_nm        = "biotic_relaxedforestnorm_sum_n2_to_p0",
      d_nm        = "pdsi_annual_sum_n2_to_p0",
      b_threshold = 4,
      d_threshold = -7.9,
      b_rollmax_nm = "biotic_relaxedforestnorm_rollsum3max",
      d_rollmax_nm = "pdsi_annual_rollsum3max"
    )    
  )
)




# ══════════════════════════════════════════════════════════════════════════════
# 6. Weighting specs ----
# ══════════════════════════════════════════════════════════════════════════════

broad_weighting_specs <- dplyr::bind_rows(
  # make_weighting_spec(
  #   weighting_id   = "glm_ato_topoclimnfg",
  #   weight_formula = ~ tt_normal_aet + srtm + tpi + tt_normal_def + nfg_factor,
  #   method         = "glm",
  #   estimand       = "ATO"
  # ),
  make_weighting_spec(
    weighting_id   = "glm_ato_topoclimnfgrap",
    weight_formula = ~ tt_normal_aet + srtm + tpi + tt_normal_def + nfg_factor + gam_rap_tree_pre6_fit,
    method         = "glm",
    estimand       = "ATO"
  )
)

nfg_weighting_specs <- dplyr::bind_rows(
  # make_weighting_spec(
  #   weighting_id   = "glm_ato_topoclimnfg",
  #   weight_formula = ~ tt_normal_aet + srtm + tpi + tt_normal_def,
  #   method         = "glm",
  #   estimand       = "ATO"
  # ),
  make_weighting_spec(
    weighting_id   = "glm_ato_topoclimnfgrap",
    weight_formula = ~ tt_normal_aet + srtm + tpi + tt_normal_def + gam_rap_tree_pre6_fit,
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

sunab_formula_6yr <- paste0(
  "{outcome} ~ ",
  "sunab(FirstTreat, year, ref.p = -6, no_agg = TRUE):cd_f + ",
  "sunab(FirstTreat, year, ref.p = -6, no_agg = TRUE):cd_bf + ",
  "sunab(FirstTreat, year, ref.p = -6, no_agg = TRUE):cd_df + ",
  "sunab(FirstTreat, year, ref.p = -6, no_agg = TRUE):cd_bdf",
  " | pt_id + year"
)

sunab_formula_3yr <- paste0(
  "{outcome} ~ ",
  "sunab(FirstTreat, year, ref.p = -3, no_agg = TRUE):cd_f + ",
  "sunab(FirstTreat, year, ref.p = -3, no_agg = TRUE):cd_bf + ",
  "sunab(FirstTreat, year, ref.p = -3, no_agg = TRUE):cd_df + ",
  "sunab(FirstTreat, year, ref.p = -3, no_agg = TRUE):cd_bdf",
  " | pt_id + year"
)

sunab_formula_6yr_extended <- paste0(
  "{outcome} ~ ",
  "sunab(FirstTreat, year, ref.p = -6, no_agg = TRUE):cd_f + ",
  "sunab(FirstTreat, year, ref.p = -6, no_agg = TRUE):cd_bf + ",
  "sunab(FirstTreat, year, ref.p = -6, no_agg = TRUE):cd_df + ",
  "sunab(FirstTreat, year, ref.p = -6, no_agg = TRUE):cd_bdf + ",
  "sunab(FirstTreat, year, ref.p = -6, no_agg = TRUE):cd_bd + ",
  "sunab(FirstTreat, year, ref.p = -6, no_agg = TRUE):cd_b + ",
  "sunab(FirstTreat, year, ref.p = -6, no_agg = TRUE):cd_d",
  " | pt_id + year"
)

sunab_formula_3yr_extended <- paste0(
  "{outcome} ~ ",
  "sunab(FirstTreat, year, ref.p = -3, no_agg = TRUE):cd_f + ",
  "sunab(FirstTreat, year, ref.p = -3, no_agg = TRUE):cd_bf + ",
  "sunab(FirstTreat, year, ref.p = -3, no_agg = TRUE):cd_df + ",
  "sunab(FirstTreat, year, ref.p = -3, no_agg = TRUE):cd_bdf + ",
  "sunab(FirstTreat, year, ref.p = -3, no_agg = TRUE):cd_bd + ",
  "sunab(FirstTreat, year, ref.p = -3, no_agg = TRUE):cd_b + ",
  "sunab(FirstTreat, year, ref.p = -3, no_agg = TRUE):cd_d",
  " | pt_id + year"
)



sixyr_model_specs <- dplyr::bind_rows(

  make_model_spec(
    model_id         = "sunab_twfe_unweighted_6yr",
    formula_template = sunab_formula_6yr,
    estimator_type   = "sunab",
    term_pattern     = "^year::",
    weights_col      = NA_character_,
    feols_args       = list(mem.clean = TRUE)
  ),

  make_model_spec(
    model_id         = "sunab_twfe_glm_ato_topoclimnfgrap_6yr",
    formula_template = sunab_formula_6yr,
    estimator_type   = "sunab",
    term_pattern     = "^year::",
    weights_col      = "glm_ato_topoclimnfgrap_weights",
    feols_args       = list(mem.clean = TRUE)
  )

)

threeyr_model_specs <- dplyr::bind_rows(
  
  make_model_spec(
    model_id         = "sunab_twfe_unweighted_3yr",
    formula_template = sunab_formula_3yr,
    estimator_type   = "sunab",
    term_pattern     = "^year::",
    weights_col      = NA_character_,
    feols_args       = list(mem.clean = TRUE)
  ),
  
  make_model_spec(
    model_id         = "sunab_twfe_glm_ato_topoclimnfgrap_3yr",
    formula_template = sunab_formula_3yr,
    estimator_type   = "sunab",
    term_pattern     = "^year::",
    weights_col      = "glm_ato_topoclimnfgrap_weights",
    feols_args       = list(mem.clean = TRUE)
  )
  
)



# extended versions

sixyr_model_specs_extended <- dplyr::bind_rows(
  
  make_model_spec(
    model_id         = "sunab_twfe_unweighted_6yr_extended",
    formula_template = sunab_formula_6yr_extended,
    estimator_type   = "sunab",
    term_pattern     = "^year::",
    weights_col      = NA_character_,
    feols_args       = list(mem.clean = TRUE)
  ),
  
  make_model_spec(
    model_id         = "sunab_twfe_glm_ato_topoclimnfgrap_6yr_extended",
    formula_template = sunab_formula_6yr_extended,
    estimator_type   = "sunab",
    term_pattern     = "^year::",
    weights_col      = "glm_ato_topoclimnfgrap_weights",
    feols_args       = list(mem.clean = TRUE)
  )
  
)

threeyr_model_specs_extended <- dplyr::bind_rows(
  
  make_model_spec(
    model_id         = "sunab_twfe_unweighted_3yr_extended",
    formula_template = sunab_formula_3yr_extended,
    estimator_type   = "sunab",
    term_pattern     = "^year::",
    weights_col      = NA_character_,
    feols_args       = list(mem.clean = TRUE)
  ),
  
  make_model_spec(
    model_id         = "sunab_twfe_glm_ato_topoclimnfgrap_3yr_extended",
    formula_template = sunab_formula_3yr_extended,
    estimator_type   = "sunab",
    term_pattern     = "^year::",
    weights_col      = "glm_ato_topoclimnfgrap_weights",
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
  ),

  # Cohort-bin event study: same estimates further stratified by cohort bin.
  # Script 3 subsets by cohort_bin to run cohort-stratified inference.
  # Adjust min/max cohort years to match your data.
  make_agg_spec(
    id    = "cohort_early_late",
    agg   = "(year::-?[0-9]+):cohort::([0-9]+):(cd_.*)",
    group_fun = function(x) {
      x |>
        dplyr::mutate(
          event_time  = as.integer(stringr::str_extract(group_1, "-?[0-9]+")),
          cohort      = as.integer(group_2),
          dummy_group = group_3,
          cohort_bin  = dplyr::case_when(
            cohort >= 2002 & cohort <= 2011 ~ "early",
            cohort >= 2012 & cohort <= 2021 ~ "late",
            TRUE ~ NA_character_
          )
        ) |>
        dplyr::filter(!is.na(cohort_bin)) |>
        dplyr::select(term, event_time, cohort_bin, dummy_group)
    },
    label = "Event study by cohort bin (early 2002-2011 / late 2012-2021) and dummy group"
  )

)


# ══════════════════════════════════════════════════════════════════════════════
# Group palette (keyed on dummy column names)
# ══════════════════════════════════════════════════════════════════════════════

group_palette <- c(
  "cd_f"   = "#E69F00",
  "cd_bf"  = "#56B4E9",
  "cd_df"  = "#009E73",
  "cd_bdf" = "#CC79A7",
  "cd_bd"  = "#0072B2",
  "cd_b"   = "#D55E00",
  "cd_d"   = "#F0E442"
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

if(ram_size == 250) {
  
  small_3 <- preview_run_grid(nfg_subset_specs, outcome_specs, threeyr_treatment_group_specs,
                              threeyr_model_specs, vcov_specs, agg_specs)
  
  small_6 <- preview_run_grid(nfg_subset_specs, outcome_specs, sixyr_treatment_group_specs,
                              sixyr_model_specs, vcov_specs, agg_specs)
  
  small_3_extended <- preview_run_grid(nfg_subset_specs, outcome_specs, threeyr_treatment_group_specs_extended_strict,
                                       threeyr_model_specs_extended, vcov_specs, agg_specs)
  
  small_6_extended <- preview_run_grid(nfg_subset_specs, outcome_specs, sixyr_treatment_group_specs_extended_strict,
                                       sixyr_model_specs_extended, vcov_specs, agg_specs)
  
}



if(ram_size == 1000) {
  
  large_3 <- preview_run_grid(broadtype_subset_specs, outcome_specs, threeyr_treatment_group_specs,
                              threeyr_model_specs, vcov_specs, agg_specs)
  
  large_6 <- preview_run_grid(broadtype_subset_specs, outcome_specs, sixyr_treatment_group_specs,
                              sixyr_model_specs, vcov_specs, agg_specs)
  
  large_3_extended <- preview_run_grid(broadtype_subset_specs, outcome_specs, threeyr_treatment_group_specs_extended_strict,
                                       threeyr_model_specs_extended, vcov_specs, agg_specs)
  
  large_6_extended <- preview_run_grid(broadtype_subset_specs, outcome_specs, sixyr_treatment_group_specs_extended_strict,
                                       sixyr_model_specs_extended, vcov_specs, agg_specs)
   
}



# ══════════════════════════════════════════════════════════════════════════════
# 11. Run weighting experiment ----
# ══════════════════════════════════════════════════════════════════════════════

# tic('weighting ecoregions')
# ecoregion_results_weighting <- run_weighting_experiment(
#   dataset_spec          = dataset_spec,
#   analysis_subset_specs = ecoregion_subset_specs,
#   treatment_group_specs = initial_treatment_group_specs,
#   weighting_specs       = weighting_specs,
#   dir_out               = dir_results,
#   skip_existing         = TRUE,
#   verbose_timing        = TRUE,
#   .progress             = TRUE
# )
# toc()
# 
# failed_weighting <- purrr::keep(ecoregion_results_weighting$run_results, \(r) !is.null(r$error))
# if (length(failed_weighting) > 0) {
#   message("Failed weighting runs:")
#   purrr::walk(failed_weighting, \(r) message("  ", r$weight_run_id, ": ", r$error))
# }

if(ram_size == 250) {
  
  # normal groups
  tic('weighting 3-yr forest groups')
  nfg_results_weighting <- run_weighting_experiment(
    dataset_spec          = dataset_spec,
    analysis_subset_specs = nfg_subset_specs,
    treatment_group_specs = threeyr_treatment_group_specs,
    weighting_specs       = nfg_weighting_specs,
    dir_out               = dir_results,
    skip_existing         = TRUE,
    verbose_timing        = TRUE,
    .progress             = TRUE
  )
  toc()
  
  failed_weighting <- purrr::keep(nfg_results_weighting$run_results, \(r) !is.null(r$error))
  if (length(failed_weighting) > 0) {
    message("Failed weighting runs:")
    purrr::walk(failed_weighting, \(r) message("  ", r$weight_run_id, ": ", r$error))
  }
  
  
  tic('weighting 6-yr forest groups')
  nfg_results_weighting <- run_weighting_experiment(
    dataset_spec          = dataset_spec,
    analysis_subset_specs = nfg_subset_specs,
    treatment_group_specs = sixyr_treatment_group_specs,
    weighting_specs       = nfg_weighting_specs,
    dir_out               = dir_results,
    skip_existing         = TRUE,
    verbose_timing        = TRUE,
    .progress             = TRUE
  )
  toc()
  
  failed_weighting <- purrr::keep(nfg_results_weighting$run_results, \(r) !is.null(r$error))
  if (length(failed_weighting) > 0) {
    message("Failed weighting runs:")
    purrr::walk(failed_weighting, \(r) message("  ", r$weight_run_id, ": ", r$error))
  }
  
  
  # extended groups
  tic('weighting 3-yr forest groups')
  nfg_results_weighting <- run_weighting_experiment(
    dataset_spec          = extended_dataset_spec,
    analysis_subset_specs = nfg_subset_specs,
    treatment_group_specs = threeyr_treatment_group_specs_extended_strict,
    weighting_specs       = nfg_weighting_specs,
    dir_out               = dir_results,
    skip_existing         = TRUE,
    verbose_timing        = TRUE,
    .progress             = TRUE
  )
  toc()
  
  failed_weighting <- purrr::keep(nfg_results_weighting$run_results, \(r) !is.null(r$error))
  if (length(failed_weighting) > 0) {
    message("Failed weighting runs:")
    purrr::walk(failed_weighting, \(r) message("  ", r$weight_run_id, ": ", r$error))
  }
  
  
  tic('weighting 6-yr forest groups')
  nfg_results_weighting <- run_weighting_experiment(
    dataset_spec          = extended_dataset_spec,
    analysis_subset_specs = nfg_subset_specs,
    treatment_group_specs = sixyr_treatment_group_specs_extended_strict,
    weighting_specs       = nfg_weighting_specs,
    dir_out               = dir_results,
    skip_existing         = TRUE,
    verbose_timing        = TRUE,
    .progress             = TRUE
  )
  toc()
  
  failed_weighting <- purrr::keep(nfg_results_weighting$run_results, \(r) !is.null(r$error))
  if (length(failed_weighting) > 0) {
    message("Failed weighting runs:")
    purrr::walk(failed_weighting, \(r) message("  ", r$weight_run_id, ": ", r$error))
  }
  
}





if(ram_size == 1000) {
  
  # BROADTYPE
  
  # normal groups
  tic('weighting 3-yr broad groups')
  broad_results_weighting <- run_weighting_experiment(
    dataset_spec          = dataset_spec,
    analysis_subset_specs = broadtype_subset_specs,
    treatment_group_specs = threeyr_treatment_group_specs,
    weighting_specs       = broad_weighting_specs,
    dir_out               = dir_results,
    skip_existing         = TRUE,
    verbose_timing        = TRUE,
    .progress             = TRUE
  )
  toc()
  
  failed_weighting <- purrr::keep(broad_results_weighting$run_results, \(r) !is.null(r$error))
  if (length(failed_weighting) > 0) {
    message("Failed weighting runs:")
    purrr::walk(failed_weighting, \(r) message("  ", r$weight_run_id, ": ", r$error))
  }
  
  
  tic('weighting 6-yr broad groups')
  broad_results_weighting <- run_weighting_experiment(
    dataset_spec          = dataset_spec,
    analysis_subset_specs = broadtype_subset_specs,
    treatment_group_specs = sixyr_treatment_group_specs,
    weighting_specs       = broad_weighting_specs,
    dir_out               = dir_results,
    skip_existing         = TRUE,
    verbose_timing        = TRUE,
    .progress             = TRUE
  )
  toc()
  
  failed_weighting <- purrr::keep(broad_results_weighting$run_results, \(r) !is.null(r$error))
  if (length(failed_weighting) > 0) {
    message("Failed weighting runs:")
    purrr::walk(failed_weighting, \(r) message("  ", r$weight_run_id, ": ", r$error))
  }
  
  
  # extended groups
  tic('weighting 3-yr broad groups')
  broad_results_weighting <- run_weighting_experiment(
    dataset_spec          = extended_dataset_spec,
    analysis_subset_specs = broadtype_subset_specs,
    treatment_group_specs = threeyr_treatment_group_specs_extended_strict,
    weighting_specs       = broad_weighting_specs,
    dir_out               = dir_results,
    skip_existing         = TRUE,
    verbose_timing        = TRUE,
    .progress             = TRUE
  )
  toc()
  
  failed_weighting <- purrr::keep(broad_results_weighting$run_results, \(r) !is.null(r$error))
  if (length(failed_weighting) > 0) {
    message("Failed weighting runs:")
    purrr::walk(failed_weighting, \(r) message("  ", r$weight_run_id, ": ", r$error))
  }
  
  
  tic('weighting 6-yr broad groups')
  broad_results_weighting <- run_weighting_experiment(
    dataset_spec          = extended_dataset_spec,
    analysis_subset_specs = broadtype_subset_specs,
    treatment_group_specs = sixyr_treatment_group_specs_extended_strict,
    weighting_specs       = broad_weighting_specs,
    dir_out               = dir_results,
    skip_existing         = TRUE,
    verbose_timing        = TRUE,
    .progress             = TRUE
  )
  toc()
  
  failed_weighting <- purrr::keep(broad_results_weighting$run_results, \(r) !is.null(r$error))
  if (length(failed_weighting) > 0) {
    message("Failed weighting runs:")
    purrr::walk(failed_weighting, \(r) message("  ", r$weight_run_id, ": ", r$error))
  }
  
  
  
  #### LARGE NFG
  
  # normal groups
  tic('weighting 3-yr forest groups')
  nfg_results_weighting <- run_weighting_experiment(
    dataset_spec          = dataset_spec,
    analysis_subset_specs = nfg_large_subset_specs,
    treatment_group_specs = threeyr_treatment_group_specs,
    weighting_specs       = nfg_weighting_specs,
    dir_out               = dir_results,
    skip_existing         = TRUE,
    verbose_timing        = TRUE,
    .progress             = TRUE
  )
  toc()
  
  failed_weighting <- purrr::keep(nfg_results_weighting$run_results, \(r) !is.null(r$error))
  if (length(failed_weighting) > 0) {
    message("Failed weighting runs:")
    purrr::walk(failed_weighting, \(r) message("  ", r$weight_run_id, ": ", r$error))
  }
  
  
  tic('weighting 6-yr forest groups')
  nfg_results_weighting <- run_weighting_experiment(
    dataset_spec          = dataset_spec,
    analysis_subset_specs = nfg_large_subset_specs,
    treatment_group_specs = sixyr_treatment_group_specs,
    weighting_specs       = nfg_weighting_specs,
    dir_out               = dir_results,
    skip_existing         = TRUE,
    verbose_timing        = TRUE,
    .progress             = TRUE
  )
  toc()
  
  failed_weighting <- purrr::keep(nfg_results_weighting$run_results, \(r) !is.null(r$error))
  if (length(failed_weighting) > 0) {
    message("Failed weighting runs:")
    purrr::walk(failed_weighting, \(r) message("  ", r$weight_run_id, ": ", r$error))
  }
  
  
  # extended groups
  tic('weighting 3-yr forest groups')
  nfg_results_weighting <- run_weighting_experiment(
    dataset_spec          = extended_dataset_spec,
    analysis_subset_specs = nfg_large_subset_specs,
    treatment_group_specs = threeyr_treatment_group_specs_extended_strict,
    weighting_specs       = nfg_weighting_specs,
    dir_out               = dir_results,
    skip_existing         = TRUE,
    verbose_timing        = TRUE,
    .progress             = TRUE
  )
  toc()
  
  failed_weighting <- purrr::keep(nfg_results_weighting$run_results, \(r) !is.null(r$error))
  if (length(failed_weighting) > 0) {
    message("Failed weighting runs:")
    purrr::walk(failed_weighting, \(r) message("  ", r$weight_run_id, ": ", r$error))
  }
  
  
  tic('weighting 6-yr forest groups')
  nfg_results_weighting <- run_weighting_experiment(
    dataset_spec          = extended_dataset_spec,
    analysis_subset_specs = nfg_large_subset_specs,
    treatment_group_specs = sixyr_treatment_group_specs_extended_strict,
    weighting_specs       = nfg_weighting_specs,
    dir_out               = dir_results,
    skip_existing         = TRUE,
    verbose_timing        = TRUE,
    .progress             = TRUE
  )
  toc()
  
  failed_weighting <- purrr::keep(nfg_results_weighting$run_results, \(r) !is.null(r$error))
  if (length(failed_weighting) > 0) {
    message("Failed weighting runs:")
    purrr::walk(failed_weighting, \(r) message("  ", r$weight_run_id, ": ", r$error))
  }
  
}


rebuild_weighting_tables(dir_out = dir_results, write_csv = TRUE)


# ══════════════════════════════════════════════════════════════════════════════
# 12. Run estimation experiment ----
# ══════════════════════════════════════════════════════════════════════════════


if(ram_size == 250) {
  
  tic('sunab estimation 3 yr nfg')
  results_sunab_ecor <- run_experiment(
    dataset_spec          = dataset_spec,
    analysis_subset_specs = nfg_subset_specs,
    outcome_specs         = outcome_specs,
    treatment_group_specs = threeyr_treatment_group_specs,
    model_specs           = threeyr_model_specs,
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
  
  
  tic('sunab estimation 6 yr nfg')
  results_sunab_ecor <- run_experiment(
    dataset_spec          = dataset_spec,
    analysis_subset_specs = nfg_subset_specs,
    outcome_specs         = outcome_specs,
    treatment_group_specs = sixyr_treatment_group_specs,
    model_specs           = sixyr_model_specs,
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
  
  
  tic('sunab estimation 3 yr nfg EXTENDED')
  results_sunab_ecor <- run_experiment(
    dataset_spec          = extended_dataset_spec,
    analysis_subset_specs = nfg_subset_specs,
    outcome_specs         = outcome_specs,
    treatment_group_specs = threeyr_treatment_group_specs_extended_strict,
    model_specs           = threeyr_model_specs_extended,
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
  
  
  tic('sunab estimation 6 yr nfg EXTENDED')
  results_sunab_ecor <- run_experiment(
    dataset_spec          = extended_dataset_spec,
    analysis_subset_specs = nfg_subset_specs,
    outcome_specs         = outcome_specs,
    treatment_group_specs = sixyr_treatment_group_specs_extended_strict,
    model_specs           = sixyr_model_specs_extended,
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
  
  
}





if(ram_size == 1000) {
  
  
  ########## LARGE NFG
  
  
  tic('sunab estimation 3 yr nfg')
  results_sunab_ecor <- run_experiment(
    dataset_spec          = dataset_spec,
    analysis_subset_specs = nfg_large_subset_specs,
    outcome_specs         = outcome_specs,
    treatment_group_specs = threeyr_treatment_group_specs,
    model_specs           = threeyr_model_specs,
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
  
  
  tic('sunab estimation 6 yr nfg')
  results_sunab_ecor <- run_experiment(
    dataset_spec          = dataset_spec,
    analysis_subset_specs = nfg_large_subset_specs,
    outcome_specs         = outcome_specs,
    treatment_group_specs = sixyr_treatment_group_specs,
    model_specs           = sixyr_model_specs,
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
  
  
  tic('sunab estimation 3 yr nfg EXTENDED')
  results_sunab_ecor <- run_experiment(
    dataset_spec          = extended_dataset_spec,
    analysis_subset_specs = nfg_large_subset_specs,
    outcome_specs         = outcome_specs,
    treatment_group_specs = threeyr_treatment_group_specs_extended_strict,
    model_specs           = threeyr_model_specs_extended,
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
  
  
  tic('sunab estimation 6 yr nfg EXTENDED')
  results_sunab_ecor <- run_experiment(
    dataset_spec          = extended_dataset_spec,
    analysis_subset_specs = nfg_large_subset_specs,
    outcome_specs         = outcome_specs,
    treatment_group_specs = sixyr_treatment_group_specs_extended_strict,
    model_specs           = sixyr_model_specs_extended,
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
  
  
  ############ BROAD
  
  tic('sunab estimation 3 yr broad')
  results_sunab_ecor <- run_experiment(
    dataset_spec          = dataset_spec,
    analysis_subset_specs = broadtype_subset_specs,
    outcome_specs         = outcome_specs,
    treatment_group_specs = threeyr_treatment_group_specs,
    model_specs           = threeyr_model_specs,
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
  
  
  tic('sunab estimation 6 yr broad')
  results_sunab_ecor <- run_experiment(
    dataset_spec          = dataset_spec,
    analysis_subset_specs = broadtype_subset_specs,
    outcome_specs         = outcome_specs,
    treatment_group_specs = sixyr_treatment_group_specs,
    model_specs           = sixyr_model_specs,
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
  
  
  tic('sunab estimation 3 yr broad EXTENDED')
  results_sunab_ecor <- run_experiment(
    dataset_spec          = extended_dataset_spec,
    analysis_subset_specs = broadtype_subset_specs,
    outcome_specs         = outcome_specs,
    treatment_group_specs = threeyr_treatment_group_specs_extended_strict,
    model_specs           = threeyr_model_specs_extended,
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
  
  
  tic('sunab estimation 6 yr broad EXTENDED')
  results_sunab_ecor <- run_experiment(
    dataset_spec          = extended_dataset_spec,
    analysis_subset_specs = broadtype_subset_specs,
    outcome_specs         = outcome_specs,
    treatment_group_specs = sixyr_treatment_group_specs_extended_strict,
    model_specs           = sixyr_model_specs_extended,
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
  
}


# ══════════════════════════════════════════════════════════════════════════════
# Rebuild merged output tables ----
# ══════════════════════════════════════════════════════════════════════════════

all_estimation_tables <- rebuild_estimation_tables(dir_out = dir_results, write_csv = TRUE)
all_descriptive_tables <- rebuild_descriptive_tables(dir_out = dir_results, write_csv = TRUE)

message(glue::glue(
  "Coef rows:       {nrow(all_estimation_tables$coef_tbl)}\n",
  "Agg specs found: {paste(names(all_estimation_tables$agg_tbls), collapse = ', ')}\n",
  "Unique run_ids:  {dplyr::n_distinct(all_estimation_tables$run_registry$run_id)}"
))

relativize_result_paths(dir_results = dir_results, base = here::here())

