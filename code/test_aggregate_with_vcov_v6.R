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

#' Extract aggregated Sun-Abraham coefficients and their covariance matrix
#'
#' @description
#' `sunab_beta_vcv()` aggregates non-aggregated Sun-Abraham coefficients from a
#' `fixest::feols()` model and returns both the aggregated coefficient estimates
#' and the full variance-covariance matrix of the aggregated estimates.
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
#' Aggregation weights are reconstructed from the model matrix following the
#' same logic used by `aggregate()` for `fixest` objects: the function counts
#' the observations contributing to each selected coefficient using
#' `colSums(sign(model.matrix(...)))`. If the original model was estimated with
#' weights, the function uses weighted counts:
#'
#' `colSums(weights * sign(model.matrix(...)))`.
#'
#' The resulting weights are normalized within each aggregation group. This means
#' that the function produces model-matrix-weighted averages, not simple
#' equal-weighted averages across coefficient names.
#'
#' If `vcov_mat` is supplied, the function uses that covariance matrix instead
#' of `sunab_fixest$cov.scaled`. This allows the same aggregation matrix to be
#' applied to post-hoc covariance matrices, such as alternative clustered,
#' Conley, heteroskedastic, or user-supplied covariance matrices. The covariance
#' matrix must correspond to the same non-aggregated coefficient vector and must
#' have row and column names matching the model coefficient names.
#'
#' `group_fun` can be used to modify the captured grouping variables before
#' aggregation. This is useful for custom aggregations that cannot be expressed
#' with regex capture groups alone, such as recoding cohort years into bins
#' before aggregation. The function passed to `group_fun` receives a data frame
#' containing the captured groups plus a `term` column. It must return a data
#' frame that includes `term` and one or more grouping columns. Rows may be
#' filtered to drop terms from the aggregation.
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
#' This function uses `dplyr::filter()` internally and the examples use several
#' additional `dplyr` verbs. The `dplyr` package should be installed and
#' available.
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
#' The function assumes that `model.matrix(sunab_fixest)` can be reconstructed.
#' It may fail for lean model objects or objects where the model matrix/data
#' needed by `model.matrix()` have been removed.
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
#' aggregation groups, optional post-hoc covariance matrices, and optional
#' user-defined recoding of aggregation groups through `group_fun`.
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
#' # Aggregate over cohorts within event-time-by-subgroup cells
#' agg_cd_event <- sunab_beta_vcv(
#'   est_sunab_dummy,
#'   agg = "(year::-?[0-9]+):cohort::[0-9]+:(cd_.*)"
#' )
#'
#' # Aggregate over cohorts and event times 2 through 6, separately by subgroup
#' att_2_6_cd <- sunab_beta_vcv(
#'   est_sunab_dummy,
#'   agg = "year::[2-6]:cohort::[0-9]+:(cd_.*)"
#' )
#'
#' # Use a post-hoc covariance matrix
#' V_alt <- vcov(
#'   est_sunab_dummy,
#'   vcov = ~ pt_id + year
#' )
#'
#' att_2_6_cd_alt <- sunab_beta_vcv(
#'   est_sunab_dummy,
#'   agg = "year::[2-6]:cohort::[0-9]+:(cd_.*)",
#'   vcov_mat = V_alt
#' )
#'
#' # Example: compare cd_bf versus cd_f for the event-time 2--6 average
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
#' # Example: aggregate by event time, cohort bin, and subgroup
#' es_by_cohort_bin <- sunab_beta_vcv(
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
#' att_2_6_by_cohort_bin <- sunab_beta_vcv(
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
#' # Validate estimates and SEs against aggregate()
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
#' # Validate selected pairwise covariance entries
#' A <- agg_cd_event$transform
#' V <- est_sunab_dummy$cov.scaled[colnames(A), colnames(A), drop = FALSE]
#'
#' i <- 1
#' j <- 2
#'
#' direct_cov_ij <- as.numeric(A[i, , drop = FALSE] %*% V %*% t(A[j, , drop = FALSE]))
#' stored_cov_ij <- agg_cd_event$sigma[i, j]
#'
#' all.equal(direct_cov_ij, stored_cov_ij, tolerance = 1e-8)
#' }
#'
#' @seealso
#' [fixest::sunab()], [aggregate()], [fixest::vcov.fixest()]
#'
#' @export
sunab_beta_vcv <- function(
    sunab_fixest,
    agg,
    vcov_mat = NULL,
    group_fun = NULL
) {
  
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
  
  mm <- model.matrix(sunab_fixest)[, agg_coef_names, drop = FALSE]
  
  if (!is.null(sunab_fixest$weights)) {
    coef_wgt <- colSums(sunab_fixest$weights * sign(mm))
  } else {
    coef_wgt <- colSums(sign(mm))
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
    parsed_terms = groups
  )
}

#Operate

# Run disagregated, then compute estimates and vcov with new function at the period level for each subgroup
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

# Run aggregation in feols() call and compare estimates and resulting vcov

# Now compute a different cluster error (post-hoc), use it to aggregate using new funciton, and then compare it to computing that cluster error in main FEOLS call with built-in aggregation



# Helpers

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
      "Could not identify estimate/SE columns in fixest::aggregate() output. ",
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
  
  list(
    comparison = comp,
    estimates_match = all(comp$estimate_match, na.rm = FALSE),
    ses_match = all(comp$se_match, na.rm = FALSE),
    max_estimate_abs_diff = max(comp$estimate_abs_diff, na.rm = TRUE),
    max_se_abs_diff = max(comp$se_abs_diff, na.rm = TRUE)
  )
}



# Test
# ---- Test sunab_beta_vcv() against fixest::aggregate() ----

agg_cd_event <- "(year::-?[0-9]+):cohort::[0-9]+:(cd.*)"
agg_2_6_cd <- "year::[2-6]:cohort::[0-9]+:(cd_.*)"

custom_cd_event <- sunab_beta_vcv(
  est_sunab_dummy,
  agg = agg_cd_event
)

fixest_cd_event <- aggregate(
  est_sunab_dummy,
  agg = agg_cd_event
)

check_cd_event <- compare_custom_to_fixest(
  custom = custom_cd_event,
  fixest_agg = fixest_cd_event
)

custom_2_6_cd <- sunab_beta_vcv(
  est_sunab_dummy,
  agg = agg_2_6_cd
)

fixest_2_6_cd <- aggregate(
  est_sunab_dummy,
  agg = agg_2_6_cd
)

check_2_6_cd <- compare_custom_to_fixest(
  custom = custom_2_6_cd,
  fixest_agg = fixest_2_6_cd
)

# ---- Test post-hoc covariance matrix ----

vcov_alt <- ~ pt_id + year

V_alt <- vcov(
  est_sunab_dummy,
  vcov = vcov_alt
)

est_sunab_dummy_alt <- summary(
  est_sunab_dummy,
  vcov = vcov_alt
)

custom_cd_event_alt <- sunab_beta_vcv(
  est_sunab_dummy,
  agg = agg_cd_event,
  vcov_mat = V_alt
)

fixest_cd_event_alt <- aggregate(
  est_sunab_dummy_alt,
  agg = agg_cd_event
)

check_cd_event_alt <- compare_custom_to_fixest(
  custom = custom_cd_event_alt,
  fixest_agg = fixest_cd_event_alt
)

custom_2_6_cd_alt <- sunab_beta_vcv(
  est_sunab_dummy,
  agg = agg_2_6_cd,
  vcov_mat = V_alt
)

fixest_2_6_cd_alt <- aggregate(
  est_sunab_dummy_alt,
  agg = agg_2_6_cd
)

check_2_6_cd_alt <- compare_custom_to_fixest(
  custom = custom_2_6_cd_alt,
  fixest_agg = fixest_2_6_cd_alt
)

# ---- Summary of checks ----

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





# validate off-diagonal vcov ----

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
  
  # Align everything
  b <- b[colnames(A)]
  V <- V[colnames(A), colnames(A), drop = FALSE]
  
  # Recompute from scratch
  beta_recomputed <- A %*% cbind(b)
  sigma_recomputed <- A %*% V %*% t(A)
  
  # Compare against stored custom output
  beta_diff <- beta_recomputed - custom$beta[rownames(A), , drop = FALSE]
  sigma_diff <- sigma_recomputed - custom$sigma[rownames(A), rownames(A), drop = FALSE]
  
  # Direct pairwise covariance checks
  # Cov(a_i'b, a_j'b) = a_i' V a_j
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


# Check that custom group_fun parameter works ----


test_dats_sn2_small_bins <- test_dats_sn2_small |>
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

x <- test_dats_sn2_small_bins |>
  dplyr::summarise(
    dplyr::across(
      dplyr::starts_with("cd_"),
      ~ sum(.x, na.rm = TRUE)
    )
  )


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

custom_cd_event_cohort_bins <- sunab_beta_vcv(
  est_sunab_dummy,
  agg = agg_cd_event_cohort,
  group_fun = cohort_bin_fun
)


agg_cd_event_bins_raw <- "(year::-?[0-9]+):cohort::[0-9]+:(cd_.*_cohort_[0-9_]+)"

fixest_cd_event_bins_raw <- aggregate(
  est_sunab_dummy_bins_raw,
  agg = agg_cd_event_bins_raw
)




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
      "Could not identify estimate/SE columns in fixest::aggregate() output. ",
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

