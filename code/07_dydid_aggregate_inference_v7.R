# 07_dydid_aggregate_inference_v7.R
# Post-estimation aggregation and inference
# -------------------------------------------------------
# Inputs (from 06_dydid_estimation_v7.R output, per run):
#   tables/by_run/{run_stub}/
#     agg_{agg_id}__{vcov_id}.rds  list(coef, vcov, groups) for inference
#     run_spec.rds                 run spec including agg_specs and vcov_specs
#
# Outputs (tables/by_run/{run_stub}/inference/):
#   pretrend_flat_{vcov_id}.parquet
#   pretrend_slope_{vcov_id}.parquet
#   att_windows_{vcov_id}.parquet
#   post_slopes_{vcov_id}.parquet
#   att_comparisons_{vcov_id}.parquet     (with BH + Holm corrections)
#   slope_comparisons_{vcov_id}.parquet   (with BH + Holm corrections)
#   cohortgroup_att_windows_{vcov_id}.parquet
#   cohortgroup_post_slopes_{vcov_id}.parquet
#
# Merged into tables/inference/all/ by rebuild_inference_tables().
# -------------------------------------------------------


# Sys.setenv(LD_LIBRARY_PATH = paste("/opt/conda/lib", Sys.getenv("LD_LIBRARY_PATH"), sep = ":"))
# Sys.setenv(PATH = paste("/usr/bin:/bin:/usr/local/bin", Sys.getenv("PATH"), sep = ":"))

rm(list = ls())

if (!requireNamespace("here", quietly = TRUE)) install.packages("here")
library(here)
here::i_am("code/07_dydid_aggregate_inference_v7.R")

required_pkgs <- c("dplyr", "readr", "purrr", "tibble", "stringr", "arrow", "glue", "here")
missing_pkgs  <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) install.packages(missing_pkgs)

library(dplyr); library(readr); library(purrr); library(tibble)
library(stringr); library(arrow); library(glue)

source(here::here("code", "sunab_aggregate_vcov.R"))
source(here::here("code", "weight_dydid_pipeline_v7.R"))


# ══════════════════════════════════════════════════════════════════════════════
# 1. Configuration ----
# ══════════════════════════════════════════════════════════════════════════════

run_name    <- "GEE_resilience_v7_operational_ss500_ts50000"
version     <- "v7"
dir_results <- here::here("results", version, "test_10perc")
dir_infout  <- file.path(dir_results, "tables", "inference", "all")
dir.create(dir_infout, recursive = TRUE, showWarnings = FALSE)


# ══════════════════════════════════════════════════════════════════════════════
# 2. Inference specification ----
# ══════════════════════════════════════════════════════════════════════════════

# Which agg_spec id to use for pre-trend and post-treatment ATT inference.
# Must match an id in the agg_specs list used in Script 06.
main_agg_id <- "event_study"

# Which agg_spec id to use for cohort-stratified inference.
# Set to NULL to skip cohort-group inference entirely.
cohortgroup_agg_id   <- "cohort_early_late"
cohortgroup_agg_col  <- "cohort_bin"          # column in groups tibble identifying strata
cohortgroup_spec_ids <- c("early", "late")    # expected values of cohort_bin

# Pre-trend windows
# NOTE: do not include ref.p (here -6) in pre_years; that estimate is
# mechanically zero and should not be tested.
pretrend_specs <- list(
  list(id = "pre_7_to_15", years = -15:-7)
)

# Post-treatment windows
posttreatment_window_specs <- list(
  list(id = "years_2_6",   years = 2:6),
  list(id = "years_7_11",  years = 7:11),
  list(id = "years_12_16", years = 12:16)
)

# Pairwise comparison pairs (reference dummy col first)
comparison_pairs <- list(
  c("cd_f", "cd_bf"),
  c("cd_f", "cd_df"),
  c("cd_f", "cd_bdf")
)

ci_level <- 0.95


# ══════════════════════════════════════════════════════════════════════════════
# 3. Discover runs ----
# ══════════════════════════════════════════════════════════════════════════════

dir_by_run <- file.path(dir_results, "tables", "by_run")
run_spec_files <- list.files(dir_by_run, pattern = "^run_spec\\.rds$",
                              recursive = TRUE, full.names = TRUE)

if (length(run_spec_files) == 0) stop("No run_spec.rds files found under: ", dir_by_run)
message(glue::glue("Discovered {length(run_spec_files)} estimation runs to process."))


# ══════════════════════════════════════════════════════════════════════════════
# 4. Inference helpers ----
# ══════════════════════════════════════════════════════════════════════════════

#' Attach run/vcov metadata columns to any output tibble
add_inference_meta <- function(tbl, meta, vcov_id) {
  dplyr::mutate(tbl,
    run_id    = meta$run_id,
    subset_id = meta$subset_id,
    outcome   = meta$outcome,
    group_id  = meta$group_id,
    model_id  = meta$model_id,
    vcov_id   = vcov_id
  )
}


#' Apply BH and Holm multiple testing corrections
#' Corrections applied within each (window_id x vcov_id) group.
apply_mtc <- function(tbl, p_col = "p", grouping = "window_id") {
  if (nrow(tbl) == 0 || !p_col %in% names(tbl)) return(tbl)
  tbl |>
    dplyr::group_by(dplyr::across(dplyr::all_of(c(grouping, "vcov_id")))) |>
    dplyr::mutate(
      p_bh   = stats::p.adjust(.data[[p_col]], method = "BH"),
      p_holm = stats::p.adjust(.data[[p_col]], method = "holm")
    ) |>
    dplyr::ungroup()
}


#' Only write non-empty tibbles
write_if_nonempty <- function(tbl, path) {
  if (is.null(tbl) || nrow(tbl) == 0) return(invisible(NULL))
  arrow::write_parquet(tbl, path)
  invisible(path)
}


#' Subset a cohort-group agg_obj to one stratum
#'
#' Extracts the rows for one cohort_group_id value and drops the cohort group
#' column, returning a standard agg_obj with event_time and dummy_group only.
#'
#' @param cg_agg_obj      agg_obj from a cohort-group agg_spec RDS.
#' @param cohort_group_id Value of cohort_group_col to keep.
#' @param cohort_group_col Column name identifying strata (default "cohort_bin").
subset_agg_by_cohort_group <- function(cg_agg_obj,
                                        cohort_group_id,
                                        cohort_group_col = "cohort_bin") {
  if (!cohort_group_col %in% names(cg_agg_obj$groups)) {
    warning(glue::glue("Column '{cohort_group_col}' not found in agg_obj$groups. ",
                       "Available: ", paste(names(cg_agg_obj$groups), collapse = ", ")))
    return(NULL)
  }

  idx <- which(cg_agg_obj$groups[[cohort_group_col]] == cohort_group_id)
  if (length(idx) == 0) {
    warning(glue::glue("No rows for cohort_group_id='{cohort_group_id}' ",
                       "in column '{cohort_group_col}'."))
    return(NULL)
  }

  # drop cohort_group_col (and key/estimate/se bookkeeping cols) so inference
  # functions only see event_time and dummy_group
  drop_cols <- intersect(names(cg_agg_obj$groups),
                         c(cohort_group_col, "key", "estimate", "se"))
  list(
    coef   = cg_agg_obj$coef[idx],
    vcov   = cg_agg_obj$vcov[idx, idx, drop = FALSE],
    groups = cg_agg_obj$groups[idx, setdiff(names(cg_agg_obj$groups), drop_cols),
                                drop = FALSE]
  )
}


# ══════════════════════════════════════════════════════════════════════════════
# 5. Per-run inference worker ----
# ══════════════════════════════════════════════════════════════════════════════

run_one_inference <- function(run_spec_file,
                               main_agg_id,
                               cohortgroup_agg_id,
                               cohortgroup_agg_col,
                               cohortgroup_spec_ids,
                               pretrend_specs,
                               posttreatment_window_specs,
                               comparison_pairs,
                               ci_level,
                               skip_existing = TRUE) {

  run_spec  <- readRDS(run_spec_file)
  run_dir   <- dirname(run_spec_file)
  inf_dir   <- file.path(run_dir, "inference")
  dir.create(inf_dir, recursive = TRUE, showWarnings = FALSE)

  run_id     <- run_spec$run_id
  dummy_cols <- run_spec$dummy_cols
  vcov_specs <- run_spec$vcov_specs

  meta <- list(
    run_id    = run_id,       subset_id = run_spec$subset_id,
    outcome   = run_spec$outcome, group_id  = run_spec$group_id,
    model_id  = run_spec$model_id
  )

  message(glue::glue("  Processing: {run_id}"))

  # verify that the requested agg_ids were computed in this run
  run_agg_ids <- purrr::map_chr(run_spec$agg_specs, "id")
  if (!main_agg_id %in% run_agg_ids) {
    warning(glue::glue("[{run_id}] main_agg_id='{main_agg_id}' not in run agg_ids: ",
                       paste(run_agg_ids, collapse = ", ")))
    return(invisible(NULL))
  }

  results <- list()

  # ── Per-vcov inference on the main agg_spec ─────────────────────────────────
  for (vi in seq_len(nrow(vcov_specs))) {
    v_id         <- vcov_specs$vcov_id[[vi]]
    agg_rds_file <- file.path(run_dir,
      glue::glue("agg_{safe_path_component(main_agg_id)}__{safe_path_component(v_id)}.rds"))

    if (!file.exists(agg_rds_file)) {
      warning(glue::glue("[{run_id}] agg RDS not found: {basename(agg_rds_file)}"))
      next
    }

    agg_obj <- readRDS(agg_rds_file)

    # pre-trend flatness
    pretrend_flat_rows <- purrr::map_dfr(pretrend_specs, function(ps) {
      purrr::map_dfr(dummy_cols, function(dc) {
        wald_pretrend_flat(agg_obj, dummy_group = dc, pre_years = ps$years) |>
          dplyr::mutate(pretrend_id = ps$id)
      })
    }) |> add_inference_meta(meta, v_id)

    # pre-trend slope
    pretrend_slope_rows <- purrr::map_dfr(pretrend_specs, function(ps) {
      purrr::map_dfr(dummy_cols, function(dc) {
        wald_pretrend_slope(agg_obj, dummy_group = dc, pre_years = ps$years) |>
          dplyr::mutate(pretrend_id = ps$id)
      })
    }) |> add_inference_meta(meta, v_id)

    # ATT windows
    att_window_rows <- purrr::map_dfr(posttreatment_window_specs, function(ws) {
      purrr::map_dfr(dummy_cols, function(dc) {
        att_window(agg_obj, dummy_group = dc, years = ws$years) |>
          dplyr::mutate(window_id = ws$id)
      })
    }) |> add_inference_meta(meta, v_id)

    # post-treatment slopes
    post_slope_rows <- purrr::map_dfr(posttreatment_window_specs, function(ws) {
      purrr::map_dfr(dummy_cols, function(dc) {
        gls_slope_one_group(agg_obj, dummy_group = dc, years = ws$years) |>
          dplyr::mutate(window_id = ws$id)
      })
    }) |> add_inference_meta(meta, v_id)

    # pairwise ATT comparisons with MTC within each window
    att_comp_rows <- purrr::map_dfr(posttreatment_window_specs, function(ws) {
      purrr::map_dfr(comparison_pairs, function(pair) {
        wald_compare_att(agg_obj, group_a = pair[1], group_b = pair[2], years = ws$years) |>
          dplyr::mutate(window_id = ws$id)
      })
    }) |>
      add_inference_meta(meta, v_id) |>
      apply_mtc(p_col = "p", grouping = "window_id")

    # pairwise slope comparisons with MTC within each window
    slope_comp_rows <- purrr::map_dfr(posttreatment_window_specs, function(ws) {
      purrr::map_dfr(comparison_pairs, function(pair) {
        compare_gls_slopes(agg_obj, group_a = pair[1], group_b = pair[2], years = ws$years) |>
          dplyr::mutate(window_id = ws$id)
      })
    }) |>
      add_inference_meta(meta, v_id) |>
      apply_mtc(p_col = "p", grouping = "window_id")

    write_if_nonempty(pretrend_flat_rows,  file.path(inf_dir, glue::glue("pretrend_flat_{safe_path_component(v_id)}.parquet")))
    write_if_nonempty(pretrend_slope_rows, file.path(inf_dir, glue::glue("pretrend_slope_{safe_path_component(v_id)}.parquet")))
    write_if_nonempty(att_window_rows,     file.path(inf_dir, glue::glue("att_windows_{safe_path_component(v_id)}.parquet")))
    write_if_nonempty(post_slope_rows,     file.path(inf_dir, glue::glue("post_slopes_{safe_path_component(v_id)}.parquet")))
    write_if_nonempty(att_comp_rows,       file.path(inf_dir, glue::glue("att_comparisons_{safe_path_component(v_id)}.parquet")))
    write_if_nonempty(slope_comp_rows,     file.path(inf_dir, glue::glue("slope_comparisons_{safe_path_component(v_id)}.parquet")))

    results[[v_id]] <- list(
      pretrend_flat     = pretrend_flat_rows,
      pretrend_slope    = pretrend_slope_rows,
      att_windows       = att_window_rows,
      post_slopes       = post_slope_rows,
      att_comparisons   = att_comp_rows,
      slope_comparisons = slope_comp_rows
    )
  }

  # ── Cohort-group inference ────────────────────────────────────────────────────
  # agg_cohort_early_late__{vcov_id}.rds was pre-computed in Script 06.
  # We load it per vcov_spec and subset by cohort_bin to run inference per stratum.

  if (!is.null(cohortgroup_agg_id) && cohortgroup_agg_id %in% run_agg_ids) {
    results$cohortgroup <- list()

    for (vi in seq_len(nrow(vcov_specs))) {
      v_id <- vcov_specs$vcov_id[[vi]]
      cg_rds_file <- file.path(run_dir,
        glue::glue("agg_{safe_path_component(cohortgroup_agg_id)}__{safe_path_component(v_id)}.rds"))

      if (!file.exists(cg_rds_file)) next
      cg_agg_obj <- readRDS(cg_rds_file)

      # ATT windows per cohort stratum
      cg_att_rows <- purrr::map_dfr(cohortgroup_spec_ids, function(cg_id) {
        cg_sub <- subset_agg_by_cohort_group(cg_agg_obj, cg_id, cohortgroup_agg_col)
        if (is.null(cg_sub)) return(tibble::tibble())
        purrr::map_dfr(posttreatment_window_specs, function(ws) {
          purrr::map_dfr(dummy_cols, function(dc) {
            att_window(cg_sub, dummy_group = dc, years = ws$years) |>
              dplyr::mutate(window_id = ws$id, cohort_group = cg_id)
          })
        })
      }) |> add_inference_meta(meta, v_id)

      # slopes per cohort stratum
      cg_slope_rows <- purrr::map_dfr(cohortgroup_spec_ids, function(cg_id) {
        cg_sub <- subset_agg_by_cohort_group(cg_agg_obj, cg_id, cohortgroup_agg_col)
        if (is.null(cg_sub)) return(tibble::tibble())
        purrr::map_dfr(posttreatment_window_specs, function(ws) {
          purrr::map_dfr(dummy_cols, function(dc) {
            gls_slope_one_group(cg_sub, dummy_group = dc, years = ws$years) |>
              dplyr::mutate(window_id = ws$id, cohort_group = cg_id)
          })
        })
      }) |> add_inference_meta(meta, v_id)

      write_if_nonempty(cg_att_rows,   file.path(inf_dir, glue::glue("cohortgroup_att_windows_{safe_path_component(v_id)}.parquet")))
      write_if_nonempty(cg_slope_rows, file.path(inf_dir, glue::glue("cohortgroup_post_slopes_{safe_path_component(v_id)}.parquet")))

      results$cohortgroup[[v_id]] <- list(att_windows = cg_att_rows, post_slopes = cg_slope_rows)
    }
  }

  invisible(results)
}


# ══════════════════════════════════════════════════════════════════════════════
# 6. Process all runs ----
# ══════════════════════════════════════════════════════════════════════════════

message(glue::glue("[{Sys.time()}] Running inference for {length(run_spec_files)} runs..."))

inference_results <- purrr::map(run_spec_files, function(f) {
  tryCatch(
    run_one_inference(
      run_spec_file             = f,
      main_agg_id               = main_agg_id,
      cohortgroup_agg_id        = cohortgroup_agg_id,
      cohortgroup_agg_col       = cohortgroup_agg_col,
      cohortgroup_spec_ids      = cohortgroup_spec_ids,
      pretrend_specs            = pretrend_specs,
      posttreatment_window_specs = posttreatment_window_specs,
      comparison_pairs          = comparison_pairs,
      ci_level                  = ci_level,
      skip_existing             = TRUE
    ),
    error = function(e) {
      message(glue::glue("  [ERROR] {basename(dirname(f))}: {e$message}"))
      NULL
    }
  )
})

message(glue::glue("[{Sys.time()}] Done."))



run_one_inference(
  run_spec_file             = run_spec_files[[1]],
  main_agg_id               = main_agg_id,
  cohortgroup_agg_id        = cohortgroup_agg_id,
  cohortgroup_agg_col       = cohortgroup_agg_col,
  cohortgroup_spec_ids      = cohortgroup_spec_ids,
  pretrend_specs            = pretrend_specs,
  posttreatment_window_specs = posttreatment_window_specs,
  comparison_pairs          = comparison_pairs,
  ci_level                  = ci_level,
  skip_existing             = FALSE
)


# ══════════════════════════════════════════════════════════════════════════════
# 7. Rebuild merged inference tables ----
# ══════════════════════════════════════════════════════════════════════════════

rebuild_inference_tables <- function(dir_results, write_csv = TRUE) {
  dir_by_run <- file.path(dir_results, "tables", "by_run")
  dir_inf    <- file.path(dir_results, "tables", "inference", "all")
  dir.create(dir_inf, recursive = TRUE, showWarnings = FALSE)

  inf_files <- list.files(file.path(dir_by_run), pattern = "\\.parquet$",
                           recursive = TRUE, full.names = TRUE) |>
    grep("inference/", x = _, value = TRUE, fixed = TRUE)

  if (length(inf_files) == 0) {
    warning("No inference parquets found.")
    return(invisible(NULL))
  }

  # group by file type prefix (strip vcov_id suffix to get the type name)
  file_type_prefixes <- c(
    "pretrend_flat",
    "pretrend_slope",
    "att_windows",
    "post_slopes",
    "att_comparisons",
    "slope_comparisons",
    "cohortgroup_att_windows",
    "cohortgroup_post_slopes"
  )

  read_and_merge <- function(pattern_files) {
    if (length(pattern_files) == 0) return(NULL)
    purrr::map_dfr(pattern_files, \(f) {
      tbl <- arrow::read_parquet(f)
      tbl$.mtime <- file.mtime(f)
      tbl
    }) |>
      dplyr::group_by(run_id, vcov_id, dplyr::across(dplyr::any_of(
        c("dummy_group", "pretrend_id", "window_id", "contrast", "cohort_group")
      ))) |>
      dplyr::filter(.mtime == max(.mtime)) |>
      dplyr::ungroup() |>
      dplyr::select(-.mtime)
  }

  write_pair <- function(tbl, stem) {
    if (is.null(tbl) || nrow(tbl) == 0) return(NULL)
    pq  <- file.path(dir_inf, paste0(stem, ".parquet"))
    csv <- if (write_csv) file.path(dir_inf, paste0(stem, ".csv")) else NULL
    arrow::write_parquet(tbl, pq)
    if (write_csv) readr::write_csv(tbl, csv)
    list(parquet = pq, csv = csv)
  }

  merged <- purrr::set_names(
    purrr::map(file_type_prefixes, function(prefix) {
      prefix_files <- inf_files[stringr::str_detect(basename(inf_files),
                                                      paste0("^", prefix, "_"))]
      read_and_merge(prefix_files)
    }),
    file_type_prefixes
  )

  files <- purrr::imap(merged, function(tbl, nm) write_pair(tbl, nm))

  message("Inference tables merged:")
  purrr::iwalk(merged, function(tbl, nm) {
    if (!is.null(tbl)) {
      message(glue::glue("  {nm}: {nrow(tbl)} rows, {dplyr::n_distinct(tbl$run_id)} runs"))
    }
  })

  invisible(list(tables = merged, files = files))
}


all_inference <- rebuild_inference_tables(dir_results = dir_results, write_csv = TRUE)
