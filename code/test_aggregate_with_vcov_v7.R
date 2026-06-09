# This script tests sunab_beta_cvc, a function for aggregating sun-abraham coefficients
# in a way that exposes the vcov matrix for post-hoc hypothesis testing and formal comparisons

# The function reproduces the estimates and standard errors returned by
# aggregate() when applied to a fixest object, i.e. the aggregate.fixest()
# method used for Sun-Abraham aggregations, while also returning the full transformed
# covariance matrix A V A'. The full covariance matrix is not directly exposed
# by aggregate(), so validation is based on:
#   1. matching aggregate() estimates,
#   2. matching aggregate() standard errors,
#   3. direct verification of pairwise covariance entries.

# The function passes all assigned tests, including event period aggregation,
# multi-year aggregation, and use of custom group_fun.


rm(list = ls())

cyverse = FALSE

if(!requireNamespace("here", quietly = TRUE)) {
  install.packages("here")
}
library(here)

source(here::here("utils", "functions.R"))


install_load_packages(c(
  "tidyverse",
  "tictoc",
  "glue",
  "fixest",
  "arrow",
  "did",
  "didimputation"
))

source(here::here("code", "portable_sunab2.R"))

run <- "GEE_resilience_v6_operational_ss500_ts50000"


if(cyverse) {
  dir_figs <- file.path("~/data-store/data/iplant/home/shared/earthlab/macrosystems/tlm", run, "figs")
  dir_derived <- file.path("~/data-store/data/iplant/home/shared/earthlab/macrosystems/tlm", run, "data/derived")
  dir_parquet <- here::here(dir_derived, 'parquet')
  dir_manual <- file.path("~/data-store/data/iplant/home/shared/earthlab/macrosystems/tlm", run, "data/manual")
  dir_raw <- file.path("~/data-store/data/iplant/home/shared/earthlab/macrosystems/tlm", run, "data/raw")
  dir_results <- file.path("~/data-store/data/iplant/home/shared/earthlab/macrosystems/tlm", run, "results")
} else {
  dir_figs <- here::here('figs', run)
  dir_derived <- here::here('data', 'derived', run)
  dir_parquet <- here::here(dir_derived, 'parquet')
  dir_raw <- here::here('data', 'raw')
  dir_manual <- here::here('data', 'manual')
  dir_results <- here::here('results')
}

dir_ensure(c(dir_figs,
             dir_derived,
             dir_parquet,
             dir_manual,
             dir_raw,
             dir_results))


seed = 1234
set.seed(seed)



set_cd_groups <- function(df,
                          b_nm,
                          d_nm,
                          cd_nm,
                          b_threshold,
                          d_threshold) {
  
  if(d_threshold < 0) {
    df_new <- df |>
      dplyr::mutate(
        "{cd_nm}" := dplyr::case_when(
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
        "{cd_nm}" := dplyr::case_when(
          .data[[b_nm]] <  b_threshold & fire == 1 & .data[[d_nm]] <  d_threshold ~ "f",
          .data[[b_nm]] >= b_threshold & fire == 1 & .data[[d_nm]] <  d_threshold ~ "bf",
          .data[[b_nm]] <  b_threshold & fire == 1 & .data[[d_nm]] >= d_threshold ~ "df",
          .data[[b_nm]] >= b_threshold & fire == 1 & .data[[d_nm]] >= d_threshold ~ "bdf",
          TRUE ~ NA_character_
        )
      )
    
  }
  
  df_new <- df_new |>
    dplyr::mutate(
      "{cd_nm}" := forcats::fct_relevel(.data[[cd_nm]], c("f", "bf", "df", "bdf"))
    )
  
  return(df_new)
}


set_cd_groups_with_control <- function(df,
                                       b_nm,
                                       d_nm,
                                       cd_nm,
                                       b_threshold,
                                       d_threshold) {
  
  if(d_threshold < 0) {
    df_new <- df |>
      dplyr::mutate(
        "{cd_nm}" := dplyr::case_when(
          .data[[b_nm]] <  b_threshold & fire == 1 & .data[[d_nm]] >  d_threshold ~ "f",
          .data[[b_nm]] >= b_threshold & fire == 1 & .data[[d_nm]] >  d_threshold ~ "bf",
          .data[[b_nm]] <  b_threshold & fire == 1 & .data[[d_nm]] <= d_threshold ~ "df",
          .data[[b_nm]] >= b_threshold & fire == 1 & .data[[d_nm]] <= d_threshold ~ "bdf",
          fire == 0 ~ "control",
          TRUE ~ NA_character_
        )
      )
  } else {
    df_new <- df |>
      dplyr::mutate(
        "{cd_nm}" := dplyr::case_when(
          .data[[b_nm]] <  b_threshold & fire == 1 & .data[[d_nm]] <  d_threshold ~ "f",
          .data[[b_nm]] >= b_threshold & fire == 1 & .data[[d_nm]] <  d_threshold ~ "bf",
          .data[[b_nm]] <  b_threshold & fire == 1 & .data[[d_nm]] >= d_threshold ~ "df",
          .data[[b_nm]] >= b_threshold & fire == 1 & .data[[d_nm]] >= d_threshold ~ "bdf",
          fire == 0 ~ "control",
          TRUE ~ NA_character_
        )
      )
    
  }
  
  df_new <- df_new |>
    dplyr::mutate(
      "{cd_nm}" := forcats::fct_relevel(.data[[cd_nm]], c("control", "f", "bf", "df", "bdf"))
    )
  
  return(df_new)
}




test_dats_sn <- arrow::read_parquet(here(dir_derived, "parquet_long", "dats_long_sierranevada.parquet"))
test_dats_sn_short <- arrow::read_parquet(here(dir_derived, "parquet_short", "dats_short_sierranevada.parquet"))
test_dats_nc <- arrow::read_parquet(here(dir_derived, "parquet_long", "dats_long_northcascades.parquet"))

#test_dats_sr <- arrow::read_parquet(here(dir_derived, "parquet", "dats_long_southernrockies.parquet"))



test_dats_sn <- test_dats_sn |>
  set_cd_groups(
    b_nm = "biotic_relaxedforestnorm_5_yrs_prior_sum_yot",
    d_nm = "pdsi_annual_5_yrs_prior_threshold_n3_yot",
    cd_nm = "b10_pdsin3t1",
    b_threshold = 10,
    d_threshold = 1
  ) |>
  filter(year >= 1997) #|>
#mutate(use_tv_covars = ifelse(Time_to_Treatment <= -5 & Time_to_Treatment <= 0, 0, 1))


test_dats_nc <- test_dats_nc |>
  set_cd_groups(
    b_nm = "biotic_relaxedforestnorm_5_yrs_prior_sum_yot",
    d_nm = "pdsi_annual_5_yrs_prior_threshold_n3_yot",
    cd_nm = "b10_pdsin3t1",
    b_threshold = 10,
    d_threshold = 1
  ) |>
  filter(year >= 1997) #|>
#mutate(use_tv_covars = ifelse(Time_to_Treatment <= -5 & Time_to_Treatment <= 0, 0, 1))





test_dats_sn <- test_dats_sn |>
  dplyr::mutate(
    b10_pdsin3t1 = factor(b10_pdsin3t1)
  )

table(test_dats_sn$b10_pdsin3t1, useNA = "ifany")

with(
  test_dats_sn,
  table(b10_pdsin3t1, is.na(FirstTreat), useNA = "ifany")
)


# DUMMY VERSION

test_dats_sn2 <- test_dats_sn |>
  dplyr::mutate(
    cd_f   = as.integer(b10_pdsin3t1 == "f"),
    cd_bf  = as.integer(b10_pdsin3t1 == "bf"),
    cd_df  = as.integer(b10_pdsin3t1 == "df"),
    cd_bdf = as.integer(b10_pdsin3t1 == "bdf"),
    dplyr::across(
      c(cd_f, cd_bf, cd_df, cd_bdf),
      ~ tidyr::replace_na(.x, 0L)
    )
  )

test_dats_sn2_small <- test_dats_sn2 |>
  dplyr::mutate(
    event_time = year - FirstTreat
  ) |>
  dplyr::filter(
    dplyr::between(event_time, -7, 10) |
      FirstTreat == 1000
  ) |>
  dplyr::filter(
    pt_id %in% sample(unique(test_dats_sn2$pt_id), 10000)
  )


# Function to test

# This function reproduces the estimates and standard errors returned by
# aggregate() when applied to a fixest object, i.e. the aggregate.fixest()
# method used for Sun-Abraham aggregations, while also returning the full
# transformed covariance matrix A V A'. The full covariance matrix is not
# directly exposed by aggregate(), so validation is based on:
#   1. matching aggregate() estimates,
#   2. matching aggregate() standard errors,
#   3. direct verification of pairwise covariance entries.
#
# Function base/motivation:
# This function builds on the approach discussed in fixest GitHub Issue #295:
# https://github.com/lrberge/fixest/issues/295

#' Extract aggregated Sun-Abraham coefficients and their covariance matrix
#'
#' @description
#' `sunab_aggregate_vcov()` aggregates non-aggregated Sun-Abraham coefficients
#' from a `fixest::feols()` model and returns both the aggregated coefficient
#' estimates and the full variance-covariance matrix of the aggregated estimates.
#'
#' The function is intended for models estimated with `fixest::sunab(...,
#' no_agg = TRUE)`, especially when the Sun-Abraham terms are interacted with
#' subgroup indicators and the user needs post-estimation inference on aggregated
#' effects. It reproduces the estimates and standard errors returned by
#' `aggregate()` for regex-defined aggregations, while also returning the full
#' transformed covariance matrix:
#'
#' `Sigma_agg = A %*% Sigma %*% t(A)`
#'
#' where `A` is the aggregation matrix and `Sigma` is the covariance matrix of
#' the non-aggregated model coefficients.
#'
#' This is useful for Wald tests, custom contrasts, subgroup comparisons,
#' event-window averages, cohort-bin comparisons, and other post-estimation tests
#' that require the covariance between aggregated coefficients, not just their
#' individual standard errors.
#'
#' @details
#' The aggregation is defined by a regular expression supplied to `agg`.
#' Coefficients whose names match `agg` are selected for aggregation. Capture
#' groups in the regex define the aggregation groups.
#'
#' A key rule is:
#'
#' - To preserve a dimension in the output, capture it with parentheses.
#' - To aggregate over a dimension, match it but do not capture it.
#'
#' For example, with coefficient names like:
#'
#' `year::2:cohort::2005:cd_bf`
#'
#' the regex:
#'
#' `"(year::-?[0-9]+):cohort::[0-9]+:(cd_.*)"`
#'
#' captures:
#'
#' - `group_1`: the event-time term, e.g. `"year::2"`
#' - `group_2`: the subgroup term, e.g. `"cd_bf"`
#'
#' and aggregates over cohorts within each event-time-by-subgroup cell.
#'
#' A regex such as:
#'
#' `"year::[2-6]:cohort::[0-9]+:(cd_.*)"`
#'
#' captures only the subgroup term and therefore aggregates over both cohorts
#' and event times 2 through 6, producing one estimate per subgroup.
#'
#' For wider event-time windows, non-capturing groups can be useful. For example:
#'
#' `"year::(?:[2-9]|1[0-5]):cohort::[0-9]+:(cd_.*)"`
#'
#' matches event times 2 through 15 but captures only the subgroup, so the output
#' is aggregated over event time and cohort, with one row per subgroup.
#'
#' @section Aggregation weights:
#' By default, `weight_method = "model_matrix"` reconstructs aggregation weights
#' using the same model-matrix logic used by `aggregate()` for `fixest` objects:
#'
#' `colSums(sign(model.matrix(...)))`
#'
#' or, for weighted models:
#'
#' `colSums(weights * sign(model.matrix(...)))`.
#'
#' The resulting coefficient-level weights are normalized within each aggregation
#' group. This means that the function produces model-matrix-weighted averages,
#' not simple equal-weighted averages across coefficient names.
#'
#' The `"model_matrix"` path is the safest reference path, but it can be
#' memory-intensive for large models because it may materialize an `N x K` model
#' matrix, where `N` is the number of observations and `K` is the number of
#' selected non-aggregated coefficients.
#'
#' If `weight_method = "data_count"`, the function avoids calling
#' `model.matrix()` and instead computes aggregation weights from `df_est` by
#' counting, or weighted-counting, observations in each coefficient cell. This is
#' much more memory-efficient, but it assumes coefficient names follow the usual
#' `sunab(..., no_agg = TRUE):dummy` pattern, e.g.:
#'
#' `year::2:cohort::2005:cd_bf`
#'
#' In this case, the function counts observations where:
#'
#' - `period_var - cohort_var == 2`;
#' - `cohort_var == 2005`;
#' - `cd_bf != 0`.
#'
#' If the interaction column is signed rather than strictly 0/1, the data-count
#' path uses `sign(interaction_value)`, matching the model-matrix weighting logic.
#'
#' @section Data-count weighting assumptions:
#' The `data_count` path is a fast path for dummy- or signed-dummy-interacted
#' Sun-Abraham designs where each selected coefficient corresponds to a cell
#' defined by:
#'
#' - event time, parsed from the coefficient name;
#' - cohort, parsed from the coefficient name;
#' - an interaction column in `df_est`, parsed from the coefficient name.
#'
#' It assumes that `df_est` is the exact estimation sample used by `feols()`.
#' If `feols()` dropped rows due to missingness, singleton fixed effects,
#' collinearity handling, weights, or other preprocessing, `df_est` must already
#' reflect those dropped rows.
#'
#' If the original `feols()` model was estimated with weights and
#' `weight_method = "data_count"`, `weight_var` should be supplied and should
#' point to the same weight variable used in `feols()`. Otherwise the aggregation
#' weights will be unweighted and may not match `aggregate()`.
#'
#' The `data_count` path may not reproduce the model-matrix weights for
#' continuous interactions, transformed variables, or specifications where the
#' selected coefficient columns are not simple cohort-by-event-time-by-interaction
#' cells.
#'
#' Always inspect `names(sunab_fixest$coefficients)` before using
#' `weight_method = "data_count"`. If the coefficient names differ from the
#' defaults, adjust `event_time_regex`, `cohort_regex`, and `interaction_regex`.
#'
#' The `data_count` path should be validated against `weight_method =
#' "model_matrix"` on a smaller dataset before production use.
#'
#' @section Post-hoc covariance matrices:
#' If `vcov_mat` is supplied, the function uses that covariance matrix instead
#' of `sunab_fixest$cov.scaled`. This allows the same aggregation matrix to be
#' applied to post-hoc covariance matrices, such as alternative clustered,
#' Conley, heteroskedastic, or user-supplied covariance matrices. The covariance
#' matrix must correspond to the same non-aggregated coefficient vector and must
#' have row and column names matching the model coefficient names.
#'
#' @section `group_fun`:
#' `group_fun` can be used to modify the captured grouping variables before
#' aggregation. This is useful for custom aggregations that cannot be expressed
#' with regex capture groups alone, such as recoding cohort years into bins
#' before aggregation.
#'
#' The function passed to `group_fun` receives a data frame containing the
#' captured groups plus a `term` column. It must return a data frame that
#' includes `term` and one or more grouping columns. Rows may be filtered to drop
#' terms from the aggregation, but `group_fun` may not add terms that were not
#' selected by the original `agg` regex.
#'
#' The grouping columns returned by `group_fun`, excluding `term`, define the
#' rows of the aggregated output. For example, returning columns `term`,
#' `event_time`, `cohort_bin`, and `cd` produces one output row per
#' event-time-by-cohort-bin-by-subgroup cell.
#'
#' @section Important:
#' This function should be used with models estimated using
#' `sunab(..., no_agg = TRUE)`. It is not intended for models where `coef()`
#' returns already-aggregated `sunab()` coefficients. In particular,
#' `sunab(..., no_agg = FALSE)` can return an aggregated coefficient view while
#' still retaining an underlying non-aggregated covariance structure, which is
#' not the target use case for this function.
#'
#' @section Dependencies:
#' This function uses `dplyr` internally. The `data_count` path additionally
#' uses `stringr`, `tibble`, and `tidyr`. These packages should be installed and
#' available.
#'
#' The examples use the native R pipe `|>` and additional `dplyr` verbs.
#'
#' @section Limitations:
#' This function reconstructs the model-matrix-weighted aggregation used by
#' `aggregate()` for `fixest` Sun-Abraham coefficients. It is not a general-purpose
#' coefficient-averaging function unless model-matrix weights are the desired
#' estimand.
#'
#' For ordinary TWFE event-study models, equal-weighted aggregation across event
#' periods may sometimes be more appropriate than model-matrix-weighted
#' aggregation.
#'
#' The `model_matrix` path assumes that `model.matrix(sunab_fixest)` can be
#' reconstructed. It may fail for lean model objects or objects where the model
#' matrix/data needed by `model.matrix()` have been removed.
#'
#' The supplied `vcov_mat`, if used, must correspond to the same non-aggregated
#' coefficient vector and must use coefficient names matching
#' `names(sunab_fixest$coefficients)`.
#'
#' `aggregate()` does not directly expose the full aggregated covariance matrix,
#' so the full `sigma` matrix returned here cannot be compared to a built-in
#' full aggregated covariance matrix. Validation should instead check that:
#'
#' - aggregated estimates match `aggregate()`;
#' - aggregated standard errors match `aggregate()`;
#' - pairwise covariance entries satisfy
#'   `sigma[i, j] = A[i, ] %*% V %*% A[j, ]`.
#'
#' @section Attribution:
#' This function builds on the approach discussed in `fixest` GitHub Issue #295,
#' "Post-estimation linear hypothesis testing using sunab()", which raised the
#' need to obtain the variance-covariance matrix for aggregated `sunab()`
#' coefficients in order to conduct post-estimation linear hypothesis tests:
#' <https://github.com/lrberge/fixest/issues/295>.
#'
#' This implementation generalizes that idea by allowing arbitrary regex-defined
#' aggregation groups, optional post-hoc covariance matrices, optional
#' user-defined recoding of aggregation groups through `group_fun`, and an
#' optional data-count weighting path to avoid materializing the model matrix.
#'
#' @param sunab_fixest A `fixest` model object, typically returned by
#'   `fixest::feols()`, estimated with one or more `fixest::sunab(...,
#'   no_agg = TRUE)` terms.
#'
#' @param agg Character string. A Perl-compatible regular expression used to
#'   select and group non-aggregated Sun-Abraham coefficient names. Coefficients
#'   matching `agg` are selected. Capture groups in `agg` define the aggregation
#'   groups. At least one capture group is required.
#'
#' @param vcov_mat Optional numeric covariance matrix. If `NULL`, the function
#'   uses `sunab_fixest$cov.scaled`. If supplied, `vcov_mat` should be a
#'   covariance matrix for the non-aggregated model coefficients, with row and
#'   column names matching `names(sunab_fixest$coefficients)`.
#'
#' @param group_fun Optional function used to modify or recode the captured
#'   grouping variables before aggregation. The function receives a data frame
#'   with columns `group_1`, `group_2`, ..., and `term`. It must return a data
#'   frame containing `term` and at least one grouping column. Rows may be
#'   filtered to drop terms from the aggregation, but `group_fun` may not add
#'   terms that were not selected by the original `agg` regex.
#'
#' @param weight_method Character. Either `"model_matrix"` or `"data_count"`.
#'   `"model_matrix"` is the default and reproduces `aggregate()` weighting
#'   using `model.matrix()`. `"data_count"` computes weights directly from
#'   `df_est` and avoids materializing the model matrix.
#'
#' @param df_est Data frame used when `weight_method = "data_count"`. This
#'   should be the exact estimation sample used by the `fixest` model.
#'
#' @param cohort_var Character. Name of the cohort/treatment-timing variable in
#'   `df_est`, e.g. `"FirstTreat"`. Required for `weight_method = "data_count"`.
#'
#' @param period_var Character. Name of the time-period variable in `df_est`,
#'   e.g. `"year"`. Required for `weight_method = "data_count"`.
#'
#' @param weight_var Optional character. Name of a column in `df_est` containing
#'   estimation weights. If `NULL`, unweighted counts are used in the
#'   `data_count` path. Use this only if the model was estimated with the same
#'   weights.
#'
#' @param event_time_regex Character regex used by the `data_count` path to
#'   parse event time from coefficient names. The first capture group must be
#'   the event-time value. The default assumes coefficient names contain strings
#'   like `"year::2"`.
#'
#' @param cohort_regex Character regex used by the `data_count` path to parse
#'   cohort from coefficient names. The first capture group must be the cohort
#'   value. The default assumes coefficient names contain strings like
#'   `"cohort::2005"`.
#'
#' @param interaction_regex Character regex used by the `data_count` path to
#'   parse the interacted variable from coefficient names. The first capture
#'   group must be the interaction variable name, e.g. `"cd_bf"`. The default
#'   assumes the interaction variable is the final colon-delimited component of
#'   the coefficient name.
#'
#' @returns
#' A list with the following elements:
#'
#' \describe{
#'   \item{beta}{A matrix of aggregated coefficient estimates. Rows correspond
#'   to aggregation groups and the single column is named `"estimate"`.}
#'
#'   \item{sigma}{The full variance-covariance matrix of the aggregated
#'   coefficient estimates, computed as `A %*% V %*% t(A)`.}
#'
#'   \item{transform}{The aggregation matrix `A`. Multiplying `A` by the
#'   non-aggregated coefficient vector produces the aggregated coefficients.}
#'
#'   \item{groups}{A data frame describing the aggregation groups, including
#'   the grouping variables, a `key` column, the aggregated estimate, and its
#'   standard error.}
#'
#'   \item{coef_names}{The original non-aggregated coefficient names used in
#'   the aggregation.}
#'
#'   \item{parsed_terms}{A data frame mapping each selected coefficient term to
#'   its parsed and, if applicable, recoded aggregation group.}
#'
#'   \item{coef_weights}{The unnormalized coefficient-level weights used to
#'   build the rows of `A`. These are useful for debugging and for validating
#'   whether `model_matrix` and `data_count` produce the same aggregation
#'   weights.}
#' }
#'
#' @examples
#' \dontrun{
#' # Disaggregated Sun-Abraham model with subgroup-specific effects
#' est_sunab_dummy <- fixest::feols(
#'   rap_tree ~
#'     sunab(FirstTreat, year, ref.p = -6, no_agg = TRUE):cd_f +
#'     sunab(FirstTreat, year, ref.p = -6, no_agg = TRUE):cd_bf +
#'     sunab(FirstTreat, year, ref.p = -6, no_agg = TRUE):cd_df +
#'     sunab(FirstTreat, year, ref.p = -6, no_agg = TRUE):cd_bdf |
#'     pt_id + year,
#'   data = test_dats_sn2_small,
#'   cluster = ~ pt_id
#' )
#'
#' # Aggregate over cohorts within event-time-by-subgroup cells.
#' agg_cd_event <- sunab_aggregate_vcov(
#'   est_sunab_dummy,
#'   agg = "(year::-?[0-9]+):cohort::[0-9]+:(cd_.*)"
#' )
#'
#' # Aggregate over cohorts and event times 2 through 6, separately by subgroup.
#' # The event-time window is matched but not captured.
#' att_2_6_cd <- sunab_aggregate_vcov(
#'   est_sunab_dummy,
#'   agg = "year::[2-6]:cohort::[0-9]+:(cd_.*)"
#' )
#'
#' # Use a post-hoc covariance matrix.
#' V_alt <- vcov(
#'   est_sunab_dummy,
#'   vcov = ~ pt_id + year
#' )
#'
#' att_2_6_cd_alt <- sunab_aggregate_vcov(
#'   est_sunab_dummy,
#'   agg = "year::[2-6]:cohort::[0-9]+:(cd_.*)",
#'   vcov_mat = V_alt
#' )
#'
#' # Use the faster data-count path. `df_est` should be the exact estimation
#' # sample used by feols().
#' agg_cd_event_data <- sunab_aggregate_vcov(
#'   est_sunab_dummy,
#'   agg = "(year::-?[0-9]+):cohort::[0-9]+:(cd_.*)",
#'   weight_method = "data_count",
#'   df_est = test_dats_sn2_small,
#'   cohort_var = "FirstTreat",
#'   period_var = "year"
#' )
#'
#' # Validate the data-count path against the model-matrix path on a smaller
#' # dataset before using it in production.
#' all.equal(
#'   agg_cd_event$transform,
#'   agg_cd_event_data$transform,
#'   tolerance = 1e-12
#' )
#'
#' all.equal(
#'   agg_cd_event$beta,
#'   agg_cd_event_data$beta,
#'   tolerance = 1e-10
#' )
#'
#' all.equal(
#'   agg_cd_event$sigma,
#'   agg_cd_event_data$sigma,
#'   tolerance = 1e-10
#' )
#'
#' # Example: compare cd_bf versus cd_f for the event-time 2--6 average.
#' g <- att_2_6_cd$groups
#'
#' i_bf <- which(g$key == "cd_bf")
#' i_f  <- which(g$key == "cd_f")
#'
#' L <- matrix(0, nrow = 1, ncol = nrow(g))
#' L[1, i_bf] <- 1
#' L[1, i_f]  <- -1
#'
#' est <- as.numeric(L %*% att_2_6_cd$beta)
#' se <- sqrt(as.numeric(L %*% att_2_6_cd$sigma %*% t(L)))
#' z <- est / se
#' p <- 2 * pnorm(abs(z), lower.tail = FALSE)
#'
#' tibble::tibble(
#'   contrast = "cd_bf - cd_f, event times 2--6",
#'   estimate = est,
#'   se = se,
#'   z = z,
#'   p = p
#' )
#'
#' # Example: aggregate by event time, cohort bin, and subgroup.
#' es_by_cohort_bin <- sunab_aggregate_vcov(
#'   est_sunab_dummy,
#'   agg = "(year::-?[0-9]+):cohort::([0-9]+):(cd_.*)",
#'   group_fun = function(x) {
#'     x |>
#'       dplyr::mutate(
#'         event_time = group_1,
#'         cohort = as.integer(group_2),
#'         cd = group_3,
#'         cohort_bin = dplyr::case_when(
#'           cohort >= 2000 & cohort <= 2010 ~ "cohort_2000_2010",
#'           cohort >= 2011 & cohort <= 2020 ~ "cohort_2011_2020",
#'           TRUE ~ NA_character_
#'         )
#'       ) |>
#'       dplyr::filter(!is.na(cohort_bin)) |>
#'       dplyr::select(term, event_time, cohort_bin, cd)
#'   }
#' )
#'
#' # Example: aggregate event times 2 through 6 by cohort bin and subgroup.
#' # The event-time window is matched but not captured, so the output rows are
#' # cohort_bin x subgroup rather than event_time x cohort_bin x subgroup.
#' att_2_6_by_cohort_bin <- sunab_aggregate_vcov(
#'   est_sunab_dummy,
#'   agg = "year::[2-6]:cohort::([0-9]+):(cd_.*)",
#'   group_fun = function(x) {
#'     x |>
#'       dplyr::mutate(
#'         cohort = as.integer(group_1),
#'         cd = group_2,
#'         cohort_bin = dplyr::case_when(
#'           cohort >= 2000 & cohort <= 2010 ~ "cohort_2000_2010",
#'           cohort >= 2011 & cohort <= 2020 ~ "cohort_2011_2020",
#'           TRUE ~ NA_character_
#'         )
#'       ) |>
#'       dplyr::filter(!is.na(cohort_bin)) |>
#'       dplyr::select(term, cohort_bin, cd)
#'   }
#' )
#'
#' # Validate estimates and SEs against aggregate().
#' agg_check <- aggregate(
#'   est_sunab_dummy,
#'   agg = "(year::-?[0-9]+):cohort::[0-9]+:(cd_.*)"
#' )
#'
#' all.equal(
#'   as.numeric(agg_check[, "Estimate"]),
#'   as.numeric(agg_cd_event$beta),
#'   tolerance = 1e-8
#' )
#'
#' all.equal(
#'   as.numeric(agg_check[, "Std. Error"]),
#'   sqrt(diag(agg_cd_event$sigma)),
#'   tolerance = 1e-8
#' )
#'
#' # Validate selected pairwise covariance entries.
#' A <- agg_cd_event$transform
#' V <- est_sunab_dummy$cov.scaled[colnames(A), colnames(A), drop = FALSE]
#'
#' i <- 1
#' j <- 2
#'
#' direct_cov_ij <- as.numeric(
#'   A[i, , drop = FALSE] %*% V %*% t(A[j, , drop = FALSE])
#' )
#'
#' stored_cov_ij <- agg_cd_event$sigma[i, j]
#'
#' all.equal(direct_cov_ij, stored_cov_ij, tolerance = 1e-8)
#' }
#'
#' @seealso
#' [fixest::sunab()], [aggregate()], [fixest::vcov.fixest()]
#'
#' @export
sunab_aggregate_vcov <- function(
    sunab_fixest,
    agg,
    vcov_mat = NULL,
    group_fun = NULL,
    weight_method = c("model_matrix", "data_count"),
    df_est = NULL,
    cohort_var = NULL,
    period_var = NULL,
    weight_var = NULL,
    event_time_regex = "year::(-?[0-9]+)",
    cohort_regex = "cohort::([0-9]+)",
    interaction_regex = ":([^:]+)$"
) {
  
  weight_method <- match.arg(weight_method)
  
  coef_vec <- sunab_fixest$coefficients
  coef_names <- names(coef_vec)
  
  if (is.null(coef_names)) {
    stop("The model coefficient vector must have names.")
  }
  
  is_match <- grepl(agg, coef_names, perl = TRUE)
  agg_coef_names <- coef_names[is_match]
  
  if (length(agg_coef_names) == 0) {
    stop("No coefficients matched the supplied `agg` regex.")
  }
  
  matches <- regexec(agg, agg_coef_names, perl = TRUE)
  captures <- regmatches(agg_coef_names, matches)
  
  n_captures <- length(captures[[1]]) - 1
  
  if (n_captures < 1) {
    stop("`agg` must contain at least one capture group.")
  }
  
  capture_mat <- do.call(
    rbind,
    lapply(captures, function(x) x[-1])
  )
  
  groups_raw <- as.data.frame(capture_mat, stringsAsFactors = FALSE)
  names(groups_raw) <- paste0("group_", seq_len(n_captures))
  groups_raw$term <- agg_coef_names
  
  if (!is.null(group_fun)) {
    groups <- group_fun(groups_raw)
  } else {
    groups <- groups_raw
  }
  
  if (!"term" %in% names(groups)) {
    stop("`group_fun` must return a data frame that includes the `term` column.")
  }
  
  groups <- groups |>
    dplyr::filter(!is.na(term))
  
  if (nrow(groups) == 0) {
    stop("No terms remain after applying `group_fun`.")
  }
  
  if (anyDuplicated(groups$term)) {
    stop("`group_fun` returned duplicated coefficient terms.")
  }
  
  invalid_terms <- setdiff(groups$term, groups_raw$term)
  
  if (length(invalid_terms) > 0) {
    stop(
      "`group_fun` returned terms that were not selected by the original `agg` regex: ",
      paste(utils::head(invalid_terms, 5), collapse = ", "),
      if (length(invalid_terms) > 5) " ..."
    )
  }
  
  missing_terms <- setdiff(groups$term, coef_names)
  
  if (length(missing_terms) > 0) {
    stop(
      "`group_fun` returned terms that are not model coefficients: ",
      paste(utils::head(missing_terms, 5), collapse = ", "),
      if (length(missing_terms) > 5) " ..."
    )
  }
  
  agg_coef_names <- groups$term
  
  group_cols <- setdiff(names(groups), "term")
  
  if (length(group_cols) == 0) {
    stop("No grouping columns remain after applying `group_fun`.")
  }
  
  group_key <- apply(
    groups[, group_cols, drop = FALSE],
    1,
    paste,
    collapse = "::"
  )
  
  unique_groups <- unique(groups[, group_cols, drop = FALSE])
  unique_groups$key <- apply(
    unique_groups,
    1,
    paste,
    collapse = "::"
  )
  
  if (weight_method == "model_matrix") {
    
    mm <- model.matrix(sunab_fixest)[, agg_coef_names, drop = FALSE]
    
    if (!is.null(sunab_fixest$weights)) {
      coef_wgt <- colSums(sunab_fixest$weights * sign(mm))
    } else {
      coef_wgt <- colSums(sign(mm))
    }
    
  } else if (weight_method == "data_count") {
    
    coef_wgt <- compute_sunab_coef_weights_from_data(
      coef_names = agg_coef_names,
      df_est = df_est,
      cohort_var = cohort_var,
      period_var = period_var,
      weight_var = weight_var,
      event_time_regex = event_time_regex,
      cohort_regex = cohort_regex,
      interaction_regex = interaction_regex
    )
  }
  
  A <- matrix(
    0,
    nrow = nrow(unique_groups),
    ncol = length(coef_vec),
    dimnames = list(unique_groups$key, coef_names)
  )
  
  for (i in seq_len(nrow(unique_groups))) {
    idx_local <- which(group_key == unique_groups$key[i])
    idx_global <- match(agg_coef_names[idx_local], coef_names)
    
    w <- coef_wgt[idx_local]
    
    if (sum(w) == 0) {
      stop("Aggregation weights sum to zero for group: ", unique_groups$key[i])
    }
    
    A[i, idx_global] <- w / sum(w)
  }
  
  if (is.null(vcov_mat)) {
    vcov_mat <- sunab_fixest$cov.scaled
  }
  
  if (is.null(vcov_mat)) {
    stop("No covariance matrix found. Provide `vcov_mat` or use a model with `cov.scaled`.")
  }
  
  if (is.null(rownames(vcov_mat)) || is.null(colnames(vcov_mat))) {
    stop("`vcov_mat` must have rownames and colnames matching coefficient names.")
  }
  
  missing_vcov_names <- setdiff(coef_names, rownames(vcov_mat))
  
  if (length(missing_vcov_names) > 0) {
    stop(
      "`vcov_mat` is missing model coefficients: ",
      paste(utils::head(missing_vcov_names, 5), collapse = ", "),
      if (length(missing_vcov_names) > 5) " ..."
    )
  }
  
  V <- vcov_mat[coef_names, coef_names, drop = FALSE]
  
  beta <- A %*% cbind(coef_vec)
  sigma <- A %*% V %*% t(A)
  
  colnames(beta) <- "estimate"
  
  unique_groups$estimate <- as.numeric(beta)
  unique_groups$se <- sqrt(diag(sigma))
  
  list(
    beta = beta,
    sigma = sigma,
    transform = A,
    groups = unique_groups,
    coef_names = agg_coef_names,
    parsed_terms = groups,
    coef_weights = coef_wgt
  )
}


#' Compute Sun-Abraham coefficient aggregation weights from estimation data
#'
#' @description
#' Internal helper used by `sunab_aggregate_vcov()` when
#' `weight_method = "data_count"`. It parses event time, cohort, and interaction
#' variable names from non-aggregated Sun-Abraham coefficient names, then
#' computes coefficient-level weights directly from the estimation data.
#'
#' @details
#' This helper is intended for coefficient names like:
#'
#' `year::2:cohort::2005:cd_bf`
#'
#' and data columns like:
#'
#' - `FirstTreat`, the cohort/treatment-timing variable;
#' - `year`, the period variable;
#' - `cd_bf`, the interacted dummy/signed-dummy variable.
#'
#' The returned weights are designed to match:
#'
#' `colSums(sign(model.matrix(model)[, coef_names]))`
#'
#' or, when `weight_var` is supplied:
#'
#' `colSums(weights * sign(model.matrix(model)[, coef_names]))`.
#'
#' This helper should generally not be called directly by users.
#'
#' @param coef_names Character vector of coefficient names.
#' @param df_est Data frame containing the exact estimation sample.
#' @param cohort_var Character name of the cohort/treatment-timing variable.
#' @param period_var Character name of the time-period variable.
#' @param weight_var Optional character name of the model weight variable.
#' @param event_time_regex Regex with first capture group identifying event time.
#' @param cohort_regex Regex with first capture group identifying cohort.
#' @param interaction_regex Regex with first capture group identifying the
#'   interaction variable.
#'
#' @returns
#' A named numeric vector of unnormalized coefficient-level aggregation weights.
compute_sunab_coef_weights_from_data <- function(
    coef_names,
    df_est,
    cohort_var,
    period_var,
    weight_var = NULL,
    event_time_regex = "year::(-?[0-9]+)",
    cohort_regex = "cohort::([0-9]+)",
    interaction_regex = ":([^:]+)$"
) {
  
  if (is.null(df_est)) {
    stop("`df_est` must be supplied when `weight_method = 'data_count'`.")
  }
  
  if (is.null(cohort_var) || is.null(period_var)) {
    stop(
      "`cohort_var` and `period_var` must be supplied when ",
      "`weight_method = 'data_count'`."
    )
  }
  
  required_vars <- c(cohort_var, period_var)
  
  if (!is.null(weight_var)) {
    required_vars <- c(required_vars, weight_var)
  }
  
  missing_required <- setdiff(required_vars, names(df_est))
  
  if (length(missing_required) > 0) {
    stop(
      "`df_est` is missing required columns: ",
      paste(missing_required, collapse = ", ")
    )
  }
  
  event_time <- stringr::str_match(coef_names, event_time_regex)[, 2]
  cohort <- stringr::str_match(coef_names, cohort_regex)[, 2]
  interaction_var <- stringr::str_match(coef_names, interaction_regex)[, 2]
  
  parsed <- tibble::tibble(
    term = coef_names,
    event_time = as.integer(event_time),
    cohort = as.integer(cohort),
    interaction_var = interaction_var
  )
  
  if (anyNA(parsed$event_time)) {
    bad <- parsed$term[is.na(parsed$event_time)]
    stop(
      "Could not parse event time from some coefficient names. Examples: ",
      paste(utils::head(bad, 5), collapse = ", "),
      if (length(bad) > 5) " ..."
    )
  }
  
  if (anyNA(parsed$cohort)) {
    bad <- parsed$term[is.na(parsed$cohort)]
    stop(
      "Could not parse cohort from some coefficient names. Examples: ",
      paste(utils::head(bad, 5), collapse = ", "),
      if (length(bad) > 5) " ..."
    )
  }
  
  if (anyNA(parsed$interaction_var)) {
    bad <- parsed$term[is.na(parsed$interaction_var)]
    stop(
      "Could not parse interaction variable from some coefficient names. Examples: ",
      paste(utils::head(bad, 5), collapse = ", "),
      if (length(bad) > 5) " ..."
    )
  }
  
  interaction_vars <- unique(parsed$interaction_var)
  missing_interactions <- setdiff(interaction_vars, names(df_est))
  
  if (length(missing_interactions) > 0) {
    stop(
      "`df_est` is missing interaction/dummy columns parsed from coefficient names: ",
      paste(utils::head(missing_interactions, 10), collapse = ", "),
      if (length(missing_interactions) > 10) " ..."
    )
  }
  
  df_counts <- df_est |>
    dplyr::mutate(
      .event_time = .data[[period_var]] - .data[[cohort_var]],
      .cohort = .data[[cohort_var]],
      .row_weight = if (is.null(weight_var)) 1 else .data[[weight_var]]
    ) |>
    tidyr::pivot_longer(
      cols = dplyr::all_of(interaction_vars),
      names_to = "interaction_var",
      values_to = ".interaction_value"
    ) |>
    dplyr::filter(.interaction_value != 0, !is.na(.interaction_value)) |>
    dplyr::group_by(.event_time, .cohort, interaction_var) |>
    dplyr::summarise(
      coef_wgt = sum(.row_weight * sign(.interaction_value), na.rm = TRUE),
      .groups = "drop"
    )
  
  weights <- parsed |>
    dplyr::left_join(
      df_counts,
      by = c(
        "event_time" = ".event_time",
        "cohort" = ".cohort",
        "interaction_var" = "interaction_var"
      )
    ) |>
    dplyr::mutate(
      coef_wgt = tidyr::replace_na(coef_wgt, 0)
    ) |>
    dplyr::pull(coef_wgt)
  
  stats::setNames(weights, coef_names)
}

#Operate

# Run disaggregated, then compute estimates and vcov with new function
# at the period level for each subgroup
est_sunab_dummy <- feols(
  rap_tree ~
    sunab(FirstTreat, year, ref.p = -6, no_agg = TRUE):cd_f +
    sunab(FirstTreat, year, ref.p = -6, no_agg = TRUE):cd_bf +
    sunab(FirstTreat, year, ref.p = -6, no_agg = TRUE):cd_df +
    sunab(FirstTreat, year, ref.p = -6, no_agg = TRUE):cd_bdf |
    pt_id + year,
  data = test_dats_sn2_small,
  cluster = ~ pt_id
)


# Helpers ----

fixest_agg_to_df <- function(x) {
  
  out <- as.data.frame(x)
  out$term <- rownames(out)
  rownames(out) <- NULL
  
  est_col <- intersect(
    c("Estimate", "estimate", "Coef.", "Coefficient"),
    names(out)
  )[1]
  
  se_col <- intersect(
    c("Std. Error", "Std. error", "SE", "se"),
    names(out)
  )[1]
  
  if (is.na(est_col) || is.na(se_col)) {
    stop(
      "Could not identify estimate/SE columns in aggregate() output. ",
      "Column names were: ",
      paste(names(out), collapse = ", ")
    )
  }
  
  out |>
    dplyr::transmute(
      term,
      estimate_fixest = .data[[est_col]],
      se_fixest = .data[[se_col]]
    )
}


custom_agg_to_df <- function(x) {
  
  x$groups |>
    dplyr::transmute(
      term = key,
      estimate_custom = estimate,
      se_custom = se
    )
}


compare_custom_to_fixest <- function(custom, fixest_agg, tolerance = 1e-8) {
  
  custom_df <- custom_agg_to_df(custom)
  fixest_df <- fixest_agg_to_df(fixest_agg)
  
  comp <- dplyr::full_join(
    custom_df,
    fixest_df,
    by = "term"
  ) |>
    dplyr::mutate(
      estimate_diff = estimate_custom - estimate_fixest,
      se_diff = se_custom - se_fixest,
      estimate_abs_diff = abs(estimate_diff),
      se_abs_diff = abs(se_diff),
      estimate_match = estimate_abs_diff < tolerance,
      se_match = se_abs_diff < tolerance
    )
  
  unmatched <- comp |>
    dplyr::filter(
      is.na(estimate_custom) |
        is.na(estimate_fixest) |
        is.na(se_custom) |
        is.na(se_fixest)
    )
  
  matched <- comp |>
    dplyr::filter(
      !is.na(estimate_custom),
      !is.na(estimate_fixest),
      !is.na(se_custom),
      !is.na(se_fixest)
    )
  
  list(
    comparison = comp,
    unmatched = unmatched,
    estimates_match = nrow(unmatched) == 0 &&
      nrow(matched) > 0 &&
      all(matched$estimate_match),
    ses_match = nrow(unmatched) == 0 &&
      nrow(matched) > 0 &&
      all(matched$se_match),
    max_estimate_abs_diff = if (nrow(matched) == 0) NA_real_ else max(matched$estimate_abs_diff),
    max_se_abs_diff = if (nrow(matched) == 0) NA_real_ else max(matched$se_abs_diff)
  )
}


compare_custom_outputs <- function(custom, reference, tolerance = 1e-8) {
  
  custom_names <- rownames(custom$beta)
  reference_names <- rownames(reference$beta)
  
  missing_from_reference <- setdiff(custom_names, reference_names)
  missing_from_custom <- setdiff(reference_names, custom_names)
  
  common_names <- intersect(custom_names, reference_names)
  
  if (length(common_names) == 0) {
    stop("No common coefficient names between custom and reference outputs.")
  }
  
  beta_custom <- custom$beta[common_names, , drop = FALSE]
  beta_reference <- reference$beta[common_names, , drop = FALSE]
  
  sigma_custom <- custom$sigma[common_names, common_names, drop = FALSE]
  sigma_reference <- reference$sigma[common_names, common_names, drop = FALSE]
  
  A_custom <- custom$transform[common_names, , drop = FALSE]
  A_reference <- reference$transform[common_names, , drop = FALSE]
  
  se_custom <- sqrt(diag(sigma_custom))
  se_reference <- sqrt(diag(sigma_reference))
  
  comparison <- tibble::tibble(
    term = common_names,
    estimate_custom = as.numeric(beta_custom),
    estimate_reference = as.numeric(beta_reference),
    se_custom = se_custom,
    se_reference = se_reference,
    estimate_diff = estimate_custom - estimate_reference,
    se_diff = se_custom - se_reference,
    estimate_abs_diff = abs(estimate_diff),
    se_abs_diff = abs(se_diff)
  )
  
  sigma_diff <- sigma_custom - sigma_reference
  A_diff <- A_custom - A_reference
  
  list(
    comparison = comparison,
    missing_from_reference = missing_from_reference,
    missing_from_custom = missing_from_custom,
    estimates_match = length(missing_from_reference) == 0 &&
      length(missing_from_custom) == 0 &&
      all(comparison$estimate_abs_diff < tolerance),
    ses_match = length(missing_from_reference) == 0 &&
      length(missing_from_custom) == 0 &&
      all(comparison$se_abs_diff < tolerance),
    vcov_match = length(missing_from_reference) == 0 &&
      length(missing_from_custom) == 0 &&
      all(abs(sigma_diff) < tolerance),
    transform_match = length(missing_from_reference) == 0 &&
      length(missing_from_custom) == 0 &&
      all(abs(A_diff) < tolerance),
    max_estimate_abs_diff = max(comparison$estimate_abs_diff),
    max_se_abs_diff = max(comparison$se_abs_diff),
    max_vcov_abs_diff = max(abs(sigma_diff)),
    max_transform_abs_diff = max(abs(A_diff))
  )
}


# Test against aggregate() ----

agg_cd_event <- "(year::-?[0-9]+):cohort::[0-9]+:(cd.*)"
agg_2_6_cd <- "year::[2-6]:cohort::[0-9]+:(cd_.*)"

custom_cd_event <- sunab_aggregate_vcov(
  est_sunab_dummy,
  agg = agg_cd_event,
  weight_method = "model_matrix"
)

fixest_cd_event <- aggregate(
  est_sunab_dummy,
  agg = agg_cd_event
)

check_cd_event <- compare_custom_to_fixest(
  custom = custom_cd_event,
  fixest_agg = fixest_cd_event
)

custom_2_6_cd <- sunab_aggregate_vcov(
  est_sunab_dummy,
  agg = agg_2_6_cd,
  weight_method = "model_matrix"
)

fixest_2_6_cd <- aggregate(
  est_sunab_dummy,
  agg = agg_2_6_cd
)

check_2_6_cd <- compare_custom_to_fixest(
  custom = custom_2_6_cd,
  fixest_agg = fixest_2_6_cd
)


# Test post-hoc covariance matrix ----

vcov_alt <- ~ pt_id + year

V_alt <- vcov(
  est_sunab_dummy,
  vcov = vcov_alt
)

est_sunab_dummy_alt <- summary(
  est_sunab_dummy,
  vcov = vcov_alt
)

custom_cd_event_alt <- sunab_aggregate_vcov(
  est_sunab_dummy,
  agg = agg_cd_event,
  vcov_mat = V_alt,
  weight_method = "model_matrix"
)

fixest_cd_event_alt <- aggregate(
  est_sunab_dummy_alt,
  agg = agg_cd_event
)

check_cd_event_alt <- compare_custom_to_fixest(
  custom = custom_cd_event_alt,
  fixest_agg = fixest_cd_event_alt
)

custom_2_6_cd_alt <- sunab_aggregate_vcov(
  est_sunab_dummy,
  agg = agg_2_6_cd,
  vcov_mat = V_alt,
  weight_method = "model_matrix"
)

fixest_2_6_cd_alt <- aggregate(
  est_sunab_dummy_alt,
  agg = agg_2_6_cd
)

check_2_6_cd_alt <- compare_custom_to_fixest(
  custom = custom_2_6_cd_alt,
  fixest_agg = fixest_2_6_cd_alt
)


# Summary of aggregate() checks ----

test_summary <- tibble::tibble(
  test = c(
    "period x cd, original vcov",
    "event times 2-6 x cd, original vcov",
    "period x cd, post-hoc vcov",
    "event times 2-6 x cd, post-hoc vcov"
  ),
  estimates_match = c(
    check_cd_event$estimates_match,
    check_2_6_cd$estimates_match,
    check_cd_event_alt$estimates_match,
    check_2_6_cd_alt$estimates_match
  ),
  ses_match = c(
    check_cd_event$ses_match,
    check_2_6_cd$ses_match,
    check_cd_event_alt$ses_match,
    check_2_6_cd_alt$ses_match
  ),
  max_estimate_abs_diff = c(
    check_cd_event$max_estimate_abs_diff,
    check_2_6_cd$max_estimate_abs_diff,
    check_cd_event_alt$max_estimate_abs_diff,
    check_2_6_cd_alt$max_estimate_abs_diff
  ),
  max_se_abs_diff = c(
    check_cd_event$max_se_abs_diff,
    check_2_6_cd$max_se_abs_diff,
    check_cd_event_alt$max_se_abs_diff,
    check_2_6_cd_alt$max_se_abs_diff
  )
)

print(test_summary)


# Validate off-diagonal vcov ----

validate_agg_vcov <- function(
    custom,
    model,
    vcov_mat = NULL,
    tolerance = 1e-8
) {
  
  A <- custom$transform
  
  b <- coef(model)
  
  if (is.null(vcov_mat)) {
    V <- model$cov.scaled
  } else {
    V <- vcov_mat
  }
  
  if (is.null(V)) {
    stop("No covariance matrix found.")
  }
  
  if (is.null(rownames(V)) || is.null(colnames(V))) {
    stop("`V` must have rownames and colnames.")
  }
  
  if (!all(colnames(A) %in% names(b))) {
    stop("Some columns of `A` are not in the model coefficient vector.")
  }
  
  if (!all(colnames(A) %in% rownames(V))) {
    stop("Some columns of `A` are not in the covariance matrix.")
  }
  
  b <- b[colnames(A)]
  V <- V[colnames(A), colnames(A), drop = FALSE]
  
  beta_recomputed <- A %*% cbind(b)
  sigma_recomputed <- A %*% V %*% t(A)
  
  beta_diff <- beta_recomputed - custom$beta[rownames(A), , drop = FALSE]
  sigma_diff <- sigma_recomputed - custom$sigma[rownames(A), rownames(A), drop = FALSE]
  
  pairwise_cov <- matrix(
    NA_real_,
    nrow = nrow(A),
    ncol = nrow(A),
    dimnames = list(rownames(A), rownames(A))
  )
  
  for (i in seq_len(nrow(A))) {
    for (j in seq_len(nrow(A))) {
      ai <- A[i, , drop = FALSE]
      aj <- A[j, , drop = FALSE]
      pairwise_cov[i, j] <- as.numeric(ai %*% V %*% t(aj))
    }
  }
  
  pairwise_diff <- pairwise_cov - custom$sigma[rownames(A), rownames(A), drop = FALSE]
  
  eig <- eigen(
    (custom$sigma + t(custom$sigma)) / 2,
    symmetric = TRUE,
    only.values = TRUE
  )$values
  
  row_sums_nonzero <- rowSums(A != 0)
  row_weight_sums <- rowSums(A)
  
  list(
    beta_match = max(abs(beta_diff)) < tolerance,
    sigma_match = max(abs(sigma_diff)) < tolerance,
    pairwise_cov_match = max(abs(pairwise_diff)) < tolerance,
    symmetric = max(abs(custom$sigma - t(custom$sigma))) < tolerance,
    min_eigenvalue = min(eig),
    positive_semidefinite_approx = min(eig) > -tolerance,
    max_beta_abs_diff = max(abs(beta_diff)),
    max_sigma_abs_diff = max(abs(sigma_diff)),
    max_pairwise_cov_abs_diff = max(abs(pairwise_diff)),
    row_weight_sums = row_weight_sums,
    row_sums_nonzero = row_sums_nonzero
  )
}


vcov_check_cd_event <- validate_agg_vcov(
  custom = custom_cd_event,
  model = est_sunab_dummy
)

vcov_check_2_6_cd <- validate_agg_vcov(
  custom = custom_2_6_cd,
  model = est_sunab_dummy
)

vcov_check_cd_event_alt <- validate_agg_vcov(
  custom = custom_cd_event_alt,
  model = est_sunab_dummy,
  vcov_mat = V_alt
)

vcov_check_2_6_cd_alt <- validate_agg_vcov(
  custom = custom_2_6_cd_alt,
  model = est_sunab_dummy,
  vcov_mat = V_alt
)

vcov_validation_summary <- tibble::tibble(
  test = c(
    "period x cd, original vcov",
    "event times 2-6 x cd, original vcov",
    "period x cd, post-hoc vcov",
    "event times 2-6 x cd, post-hoc vcov"
  ),
  beta_match = c(
    vcov_check_cd_event$beta_match,
    vcov_check_2_6_cd$beta_match,
    vcov_check_cd_event_alt$beta_match,
    vcov_check_2_6_cd_alt$beta_match
  ),
  sigma_match = c(
    vcov_check_cd_event$sigma_match,
    vcov_check_2_6_cd$sigma_match,
    vcov_check_cd_event_alt$sigma_match,
    vcov_check_2_6_cd_alt$sigma_match
  ),
  pairwise_cov_match = c(
    vcov_check_cd_event$pairwise_cov_match,
    vcov_check_2_6_cd$pairwise_cov_match,
    vcov_check_cd_event_alt$pairwise_cov_match,
    vcov_check_2_6_cd_alt$pairwise_cov_match
  ),
  symmetric = c(
    vcov_check_cd_event$symmetric,
    vcov_check_2_6_cd$symmetric,
    vcov_check_cd_event_alt$symmetric,
    vcov_check_2_6_cd_alt$symmetric
  ),
  positive_semidefinite_approx = c(
    vcov_check_cd_event$positive_semidefinite_approx,
    vcov_check_2_6_cd$positive_semidefinite_approx,
    vcov_check_cd_event_alt$positive_semidefinite_approx,
    vcov_check_2_6_cd_alt$positive_semidefinite_approx
  ),
  max_pairwise_cov_abs_diff = c(
    vcov_check_cd_event$max_pairwise_cov_abs_diff,
    vcov_check_2_6_cd$max_pairwise_cov_abs_diff,
    vcov_check_cd_event_alt$max_pairwise_cov_abs_diff,
    vcov_check_2_6_cd_alt$max_pairwise_cov_abs_diff
  ),
  min_eigenvalue = c(
    vcov_check_cd_event$min_eigenvalue,
    vcov_check_2_6_cd$min_eigenvalue,
    vcov_check_cd_event_alt$min_eigenvalue,
    vcov_check_2_6_cd_alt$min_eigenvalue
  )
)

print(vcov_validation_summary)


# Check that group_fun parameter works against raw binned feols ----

test_dats_sn2_small_2000_2020 <- test_dats_sn2_small |>
  dplyr::filter(
    FirstTreat == 1000 |
      dplyr::between(FirstTreat, 2000, 2020)
  )

est_sunab_dummy_2000_2020 <- feols(
  rap_tree ~
    sunab(FirstTreat, year, ref.p = -6, no_agg = TRUE):cd_f +
    sunab(FirstTreat, year, ref.p = -6, no_agg = TRUE):cd_bf +
    sunab(FirstTreat, year, ref.p = -6, no_agg = TRUE):cd_df +
    sunab(FirstTreat, year, ref.p = -6, no_agg = TRUE):cd_bdf |
    pt_id + year,
  data = test_dats_sn2_small_2000_2020,
  cluster = ~ pt_id
)
test_dats_sn2_small_bins <- test_dats_sn2_small_2000_2020 |>
  dplyr::mutate(
    cohort_bin = dplyr::case_when(
      FirstTreat >= 2000 & FirstTreat <= 2010 ~ "cohort_2000_2010",
      FirstTreat >= 2011 & FirstTreat <= 2020 ~ "cohort_2011_2020",
      TRUE ~ NA_character_
    ),
    
    cd_f_cohort_2000_2010 =
      as.integer(cd_f == 1 & cohort_bin == "cohort_2000_2010"),
    cd_bf_cohort_2000_2010 =
      as.integer(cd_bf == 1 & cohort_bin == "cohort_2000_2010"),
    cd_df_cohort_2000_2010 =
      as.integer(cd_df == 1 & cohort_bin == "cohort_2000_2010"),
    cd_bdf_cohort_2000_2010 =
      as.integer(cd_bdf == 1 & cohort_bin == "cohort_2000_2010"),
    
    cd_f_cohort_2011_2020 =
      as.integer(cd_f == 1 & cohort_bin == "cohort_2011_2020"),
    cd_bf_cohort_2011_2020 =
      as.integer(cd_bf == 1 & cohort_bin == "cohort_2011_2020"),
    cd_df_cohort_2011_2020 =
      as.integer(cd_df == 1 & cohort_bin == "cohort_2011_2020"),
    cd_bdf_cohort_2011_2020 =
      as.integer(cd_bdf == 1 & cohort_bin == "cohort_2011_2020")
  )

bin_dummy_counts <- test_dats_sn2_small_bins |>
  dplyr::summarise(
    dplyr::across(
      dplyr::starts_with("cd_"),
      ~ sum(.x, na.rm = TRUE)
    )
  )

print(bin_dummy_counts)

est_sunab_dummy_bins_raw <- feols(
  rap_tree ~
    sunab(FirstTreat, year, ref.p = -6, no_agg = TRUE):cd_f_cohort_2000_2010 +
    sunab(FirstTreat, year, ref.p = -6, no_agg = TRUE):cd_bf_cohort_2000_2010 +
    sunab(FirstTreat, year, ref.p = -6, no_agg = TRUE):cd_df_cohort_2000_2010 +
    sunab(FirstTreat, year, ref.p = -6, no_agg = TRUE):cd_bdf_cohort_2000_2010 +
    sunab(FirstTreat, year, ref.p = -6, no_agg = TRUE):cd_f_cohort_2011_2020 +
    sunab(FirstTreat, year, ref.p = -6, no_agg = TRUE):cd_bf_cohort_2011_2020 +
    sunab(FirstTreat, year, ref.p = -6, no_agg = TRUE):cd_df_cohort_2011_2020 +
    sunab(FirstTreat, year, ref.p = -6, no_agg = TRUE):cd_bdf_cohort_2011_2020 |
    pt_id + year,
  data = test_dats_sn2_small_bins,
  cluster = ~ pt_id
)

agg_cd_event_cohort <- "(year::-?[0-9]+):cohort::([0-9]+):(cd_.*)"

cohort_bin_fun <- function(x) {
  
  x |>
    dplyr::mutate(
      event_time = group_1,
      cohort = as.integer(group_2),
      cd = group_3,
      cohort_bin = dplyr::case_when(
        cohort >= 2000 & cohort <= 2010 ~ "cohort_2000_2010",
        cohort >= 2011 & cohort <= 2020 ~ "cohort_2011_2020",
        TRUE ~ NA_character_
      )
    ) |>
    dplyr::filter(!is.na(cohort_bin)) |>
    dplyr::select(term, event_time, cohort_bin, cd)
}



custom_cd_event_cohort_bins <- sunab_aggregate_vcov(
  est_sunab_dummy_2000_2020,
  agg = agg_cd_event_cohort,
  group_fun = cohort_bin_fun,
  weight_method = "model_matrix"
)

agg_cd_event_bins_raw <- "(year::-?[0-9]+):cohort::[0-9]+:(cd_.*_cohort_[0-9_]+)"

fixest_cd_event_bins_raw <- aggregate(
  est_sunab_dummy_bins_raw,
  agg = agg_cd_event_bins_raw
)

normalize_raw_bin_terms <- function(x) {
  
  x |>
    stringr::str_replace(
      "^(year::-?[0-9]+)::(cd_[a-z]+)_cohort_([0-9]{4}_[0-9]{4})$",
      "\\1::cohort_\\3::\\2"
    )
}

fixest_cd_event_bins_raw_df <- fixest_agg_to_df(
  fixest_cd_event_bins_raw
) |>
  dplyr::mutate(
    term = normalize_raw_bin_terms(term)
  )

custom_cd_event_cohort_bins_df <- custom_cd_event_cohort_bins$groups |>
  dplyr::transmute(
    term = key,
    estimate_custom = estimate,
    se_custom = se
  )

check_group_fun_against_raw_feols <- custom_cd_event_cohort_bins_df |>
  dplyr::full_join(
    fixest_cd_event_bins_raw_df,
    by = "term"
  ) |>
  dplyr::mutate(
    estimate_diff = estimate_custom - estimate_fixest,
    se_diff = se_custom - se_fixest,
    estimate_abs_diff = abs(estimate_diff),
    se_abs_diff = abs(se_diff),
    estimate_match = estimate_abs_diff < 1e-8,
    se_match = se_abs_diff < 1e-8
  )


# NOTE:
# The failure of this text is expected. Splitting the interaction variables into cohort-bin-specific
# regressors changes the design matrix and can induce a different pattern of
# collinearity drops. Here, the original model and raw binned model dropped
# substantially different numbers of coefficients.

group_fun_raw_feols_summary <- tibble::tibble(
  test = "group_fun cohort bins vs raw binned feols",
  estimates_match = all(check_group_fun_against_raw_feols$estimate_match, na.rm = FALSE),
  ses_match = all(check_group_fun_against_raw_feols$se_match, na.rm = FALSE),
  max_estimate_abs_diff = max(check_group_fun_against_raw_feols$estimate_abs_diff, na.rm = TRUE),
  max_se_abs_diff = max(check_group_fun_against_raw_feols$se_abs_diff, na.rm = TRUE),
  n_unmatched = sum(
    is.na(check_group_fun_against_raw_feols$estimate_custom) |
      is.na(check_group_fun_against_raw_feols$estimate_fixest)
  )
)

print(group_fun_raw_feols_summary)




length(est_sunab_dummy_2000_2020$collin.var)
length(est_sunab_dummy_bins_raw$collin.var)

head(est_sunab_dummy_2000_2020$collin.var, 50)
head(est_sunab_dummy_bins_raw$collin.var, 50)

length(coef(est_sunab_dummy_2000_2020))
length(coef(est_sunab_dummy_bins_raw))






# Matrix-less data_count path checks ----
# These validate that weight_method = "data_count" reproduces the reference
# weight_method = "model_matrix" path on the same estimation data.

custom_cd_event_data <- sunab_aggregate_vcov(
  est_sunab_dummy,
  agg = agg_cd_event,
  weight_method = "data_count",
  df_est = test_dats_sn2_small,
  cohort_var = "FirstTreat",
  period_var = "year"
)

custom_2_6_cd_data <- sunab_aggregate_vcov(
  est_sunab_dummy,
  agg = agg_2_6_cd,
  weight_method = "data_count",
  df_est = test_dats_sn2_small,
  cohort_var = "FirstTreat",
  period_var = "year"
)

custom_cd_event_alt_data <- sunab_aggregate_vcov(
  est_sunab_dummy,
  agg = agg_cd_event,
  vcov_mat = V_alt,
  weight_method = "data_count",
  df_est = test_dats_sn2_small,
  cohort_var = "FirstTreat",
  period_var = "year"
)

custom_2_6_cd_alt_data <- sunab_aggregate_vcov(
  est_sunab_dummy,
  agg = agg_2_6_cd,
  vcov_mat = V_alt,
  weight_method = "data_count",
  df_est = test_dats_sn2_small,
  cohort_var = "FirstTreat",
  period_var = "year"
)

check_data_cd_event <- compare_custom_outputs(
  custom = custom_cd_event_data,
  reference = custom_cd_event
)

check_data_2_6_cd <- compare_custom_outputs(
  custom = custom_2_6_cd_data,
  reference = custom_2_6_cd
)

check_data_cd_event_alt <- compare_custom_outputs(
  custom = custom_cd_event_alt_data,
  reference = custom_cd_event_alt
)

check_data_2_6_cd_alt <- compare_custom_outputs(
  custom = custom_2_6_cd_alt_data,
  reference = custom_2_6_cd_alt
)

data_count_test_summary <- tibble::tibble(
  test = c(
    "period x cd, original vcov, data_count vs model_matrix",
    "event times 2-6 x cd, original vcov, data_count vs model_matrix",
    "period x cd, post-hoc vcov, data_count vs model_matrix",
    "event times 2-6 x cd, post-hoc vcov, data_count vs model_matrix"
  ),
  estimates_match = c(
    check_data_cd_event$estimates_match,
    check_data_2_6_cd$estimates_match,
    check_data_cd_event_alt$estimates_match,
    check_data_2_6_cd_alt$estimates_match
  ),
  ses_match = c(
    check_data_cd_event$ses_match,
    check_data_2_6_cd$ses_match,
    check_data_cd_event_alt$ses_match,
    check_data_2_6_cd_alt$ses_match
  ),
  vcov_match = c(
    check_data_cd_event$vcov_match,
    check_data_2_6_cd$vcov_match,
    check_data_cd_event_alt$vcov_match,
    check_data_2_6_cd_alt$vcov_match
  ),
  transform_match = c(
    check_data_cd_event$transform_match,
    check_data_2_6_cd$transform_match,
    check_data_cd_event_alt$transform_match,
    check_data_2_6_cd_alt$transform_match
  ),
  max_estimate_abs_diff = c(
    check_data_cd_event$max_estimate_abs_diff,
    check_data_2_6_cd$max_estimate_abs_diff,
    check_data_cd_event_alt$max_estimate_abs_diff,
    check_data_2_6_cd_alt$max_estimate_abs_diff
  ),
  max_se_abs_diff = c(
    check_data_cd_event$max_se_abs_diff,
    check_data_2_6_cd$max_se_abs_diff,
    check_data_cd_event_alt$max_se_abs_diff,
    check_data_2_6_cd_alt$max_se_abs_diff
  ),
  max_vcov_abs_diff = c(
    check_data_cd_event$max_vcov_abs_diff,
    check_data_2_6_cd$max_vcov_abs_diff,
    check_data_cd_event_alt$max_vcov_abs_diff,
    check_data_2_6_cd_alt$max_vcov_abs_diff
  ),
  max_transform_abs_diff = c(
    check_data_cd_event$max_transform_abs_diff,
    check_data_2_6_cd$max_transform_abs_diff,
    check_data_cd_event_alt$max_transform_abs_diff,
    check_data_2_6_cd_alt$max_transform_abs_diff
  )
)

print(data_count_test_summary)


# Optional: direct comparison of raw coefficient weights
# This is useful when data_count fails, because it isolates the problem to the
# coefficient-weight construction rather than the covariance propagation.

coef_weight_check_cd_event <- tibble::tibble(
  term = names(custom_cd_event$coef_weights),
  weight_model_matrix = as.numeric(custom_cd_event$coef_weights),
  weight_data_count = as.numeric(custom_cd_event_data$coef_weights[names(custom_cd_event$coef_weights)]),
  weight_diff = weight_data_count - weight_model_matrix,
  weight_abs_diff = abs(weight_diff)
)

coef_weight_check_cd_event |>
  dplyr::arrange(dplyr::desc(weight_abs_diff)) |>
  print(n = 20)

coef_weight_summary_cd_event <- coef_weight_check_cd_event |>
  dplyr::summarise(
    weights_match = all(weight_abs_diff < 1e-8),
    max_weight_abs_diff = max(weight_abs_diff, na.rm = TRUE)
  )

print(coef_weight_summary_cd_event)



