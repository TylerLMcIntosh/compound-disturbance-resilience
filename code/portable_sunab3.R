# portable_sunab.R
# Portable Sun-Abraham / TWFE DiD estimation pipeline
# -------------------------------------------------------
# Design principles:
#   - All dataset-specific column names live in a single dataset_spec object
#   - Analysis subsets are defined by tidy filter expressions pushed to Arrow
#   - treatment_group_specs define how treated units are partitioned
#   - Per-run folders are the source of truth; rebuild merges them
#   - One shared dir_out per project; multiple experiment calls write into it
# -------------------------------------------------------


# ── Required packages ────────────────────────────────────────────────────────

.required_pkgs <- c(
  "arrow",
  "dplyr",
  "fixest",
  "glue",
  "purrr",
  "readr",
  "rlang",
  "stringr",
  "tibble",
  "tidyr",
  "jsonlite"
)

.missing_pkgs <- .required_pkgs[
  !vapply(.required_pkgs, requireNamespace, logical(1), quietly = TRUE)
]

if (length(.missing_pkgs) > 0) {
  stop(
    "Missing required packages: ",
    paste(.missing_pkgs, collapse = ", "),
    call. = FALSE
  )
}

rm(.required_pkgs, .missing_pkgs)


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 1 — Dataset spec
# ══════════════════════════════════════════════════════════════════════════════

#' Build a dataset spec
#'
#' A dataset spec captures all column-name mappings that are properties of the
#' *data*, not of any particular model or grouping scheme. Pass one spec per
#' experiment call (all subsets must share the same column layout).
#'
#' @param unit_id      Name of the unit (e.g. pixel / plot) identifier column.
#' @param time_var     Name of the calendar-time column (e.g. "year").
#' @param trt_col      Name of the binary treatment indicator column (0/1).
#' @param cohort_var   Name of the treatment-cohort column fed to sunab()
#'                     (e.g. "FirstTreat"). NA_character_ for plain TWFE.
#' @param event_id     Name of the event identifier column (e.g. "fireid").
#'                     NA_character_ if not applicable.
#'
#' @return A named list validated for required fields.
make_dataset_spec <- function(unit_id,
                              time_var,
                              trt_col,
                              cohort_var   = NA_character_,
                              event_id     = NA_character_) {
  
  stopifnot(
    is.character(unit_id)    && length(unit_id)    == 1,
    is.character(time_var)   && length(time_var)   == 1,
    is.character(trt_col)    && length(trt_col)    == 1
  )
  
  list(
    unit_id    = unit_id,
    time_var   = time_var,
    trt_col    = trt_col,
    cohort_var = normalize_optional_colname(cohort_var),
    event_id   = normalize_optional_colname(event_id)
  )
}


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 2 — Analysis subset specs
# ══════════════════════════════════════════════════════════════════════════════

#' Build analysis subset specs from a data source and optional filter
#'
#' Each row defines one analysis subset. Filters are pushed to Arrow before
#' collect(), so only the needed rows are loaded into memory.
#'
#' @param subset_id     A short human-readable identifier for this subset.
#' @param data_source   Either a path to a directory (passed to
#'                      arrow::open_dataset) or a character vector of parquet
#'                      file paths.
#' @param data_filter   A one-sided formula filter, e.g. ~ year >= 1997, or
#'                      NULL for no filter. Applied via dplyr::filter() on the
#'                      Arrow Dataset before collect().
#'
#' @return A one-row tibble suitable for rbind-ing into an analysis_subset_specs
#'         table, which is then passed to run_experiment().
make_analysis_subset_spec <- function(subset_id,
                                      data_source,
                                      data_filter = NULL) {
  
  if (!is.character(subset_id) || length(subset_id) != 1 || is.na(subset_id)) {
    stop("subset_id must be a length-1 non-NA character string.")
  }
  
  if (!is.character(data_source) || length(data_source) == 0) {
    stop("data_source must be a non-empty character vector of file paths or a single directory path.")
  }
  
  # Validate that source exists
  if (length(data_source) == 1 && dir.exists(data_source)) {
    source_type <- "directory"
  } else {
    missing_files <- data_source[!file.exists(data_source)]
    if (length(missing_files) > 0) {
      warning(
        glue::glue(
          "subset_id '{subset_id}': these source files do not exist yet:\n",
          paste(missing_files, collapse = "\n")
        )
      )
    }
    source_type <- "files"
  }
  
  if (!is.null(data_filter) && !inherits(data_filter, "formula")) {
    stop("data_filter must be NULL or a one-sided formula, e.g. ~ year >= 1997")
  }
  
  tibble::tibble(
    subset_id   = subset_id,
    data_source = list(data_source),
    source_type = source_type,
    data_filter = list(data_filter)
  )
}


#' Auto-expand analysis subset specs by unique values of a column
#'
#' Reads the unique values of `split_col` from `data_source`, then builds one
#' spec row per value (optionally combined with a base filter).
#'
#' @param data_source   Directory or file vector (same as make_analysis_subset_spec).
#' @param split_col     Column name to split on.
#' @param id_prefix     Prefix for auto-generated subset_id values.
#' @param base_filter   Optional one-sided formula applied on top of the
#'                      per-value filter (e.g. ~ year >= 1997).
#' @param values        If not NULL, restrict to this subset of unique values.
#'
#' @return A tibble of spec rows, one per unique value.
expand_analysis_subset_specs_by_col <- function(data_source,
                                                split_col,
                                                id_prefix       = split_col,
                                                base_filter     = NULL,
                                                values          = NULL) {
  
  ds <- open_arrow_source(data_source)
  
  unique_vals <- ds |>
    dplyr::select(dplyr::all_of(split_col)) |>
    dplyr::distinct() |>
    dplyr::collect() |>
    dplyr::pull(1) |>
    stats::na.omit() |>
    sort()
  
  if (!is.null(values)) {
    unique_vals <- intersect(unique_vals, values)
    if (length(unique_vals) == 0) {
      stop("No overlap between provided values and unique values of '", split_col, "' in data.")
    }
  }
  
  purrr::map_dfr(
    unique_vals,
    function(val) {
      val_filter <- build_equality_filter(split_col, val)
      combined_filter <- combine_filters(base_filter, val_filter)
      subset_id <- glue::glue("{id_prefix}_{safe_path_component(val)}")
      make_analysis_subset_spec(
        subset_id   = subset_id,
        data_source = data_source,
        data_filter = combined_filter
      )
    }
  )
}


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 3 — Treatment group specs
# ══════════════════════════════════════════════════════════════════════════════

#' Build a treatment group spec row
#'
#' Each row defines one way of partitioning the treated units into subgroups.
#' The grouping function receives the full data frame plus all group_args, and
#' must return a data frame with a new column named `group_col`.
#'
#' @param group_id    Short identifier for this grouping scheme.
#' @param group_col   Name of the subgroup column that group_fun will create.
#'                    This is the column the model loop iterates over.
#' @param group_fun   A function with signature f(df, ...) that adds `group_col`
#'                    to df. All additional arguments are passed via group_args.
#' @param group_args  A named list of additional arguments forwarded to group_fun.
#'                    Must NOT include a key named "group_col" (that is injected
#'                    automatically from the top-level field).
#'
#' @return A one-row tibble suitable for rbind-ing into treatment_group_specs.
make_treatment_group_spec <- function(group_id,
                                      group_col,
                                      group_fun,
                                      group_args = list()) {
  
  if (!is.character(group_id) || length(group_id) != 1) {
    stop("group_id must be a length-1 character string.")
  }
  if (!is.character(group_col) || length(group_col) != 1) {
    stop("group_col must be a length-1 character string.")
  }
  if (!is.function(group_fun)) {
    stop("group_fun must be a function.")
  }
  if ("group_col" %in% names(group_args)) {
    stop("group_args must not contain 'group_col'; set it via the top-level group_col argument.")
  }
  
  tibble::tibble(
    group_id   = group_id,
    group_col  = group_col,
    group_fun  = list(group_fun),
    group_args = list(group_args)
  )
}


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 4 — Model specs
# ══════════════════════════════════════════════════════════════════════════════

#' Build a model spec row
#'
#' @param model_id         Short identifier for this model.
#' @param formula_template A glue-style string with {outcome} as the only
#'                         placeholder, e.g.
#'                         "{outcome} ~ x1 + sunab(FirstTreat, year) | id + year"
#' @param estimator_type   One of "sunab" or "twfe". Controls which
#'                         post-processing steps run (event-time support
#'                         extraction is only meaningful for "sunab").
#' @param term_pattern     Regex matched against coefficient names to flag which
#'                         terms are event-time coefficients. Default ".*".
#'
#' @return A one-row tibble suitable for rbind-ing into model_specs.
make_model_spec <- function(model_id,
                            formula_template,
                            estimator_type = c("sunab", "twfe"),
                            term_pattern   = ".*") {
  
  estimator_type <- match.arg(estimator_type)
  
  if (!is.character(model_id) || length(model_id) != 1) {
    stop("model_id must be a length-1 character string.")
  }
  if (!is.character(formula_template) || length(formula_template) != 1) {
    stop("formula_template must be a length-1 character string.")
  }
  
  tibble::tibble(
    model_id         = model_id,
    formula_template = formula_template,
    estimator_type   = estimator_type,
    term_pattern     = term_pattern
  )
}


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 5 — Top-level experiment runners
# ══════════════════════════════════════════════════════════════════════════════

#' Run a full DiD experiment
#'
#' Crosses analysis_subset_specs × outcome_specs × treatment_group_specs ×
#' model_specs and runs estimation (and optionally descriptive summaries) for
#' each combination.
#'
#' Failed runs are caught, logged, and skipped; the loop always completes.
#'
#' @param dataset_spec            Output of make_dataset_spec().
#' @param analysis_subset_specs   Tibble of subset specs (make_analysis_subset_spec rows).
#' @param outcome_specs           Tibble with column "outcome".
#' @param treatment_group_specs   Tibble of treatment group specs (make_treatment_group_spec rows).
#' @param model_specs             Tibble of model specs (make_model_spec rows).
#' @param vcov_specs              Tibble with columns vcov_id, vcov, vcov_label,
#'                                and optionally vcov_vars.
#' @param dir_out                 Root output directory for this project.
#' @param group_palette           Optional named character vector mapping subgroup
#'                                values to hex colors.
#' @param ci_level                Confidence level for intervals. Default 0.95.
#' @param run_estimation          Logical. Run estimation pipeline? Default TRUE.
#' @param run_descriptive         Logical. Run descriptive pipeline? Default FALSE.
#' @param descriptive_args        Named list of extra args for the descriptive
#'                                pipeline: treated_year_var, control_year_var.
#'                                Both default to dataset_spec$time_var if omitted.
#' @param skip_existing           If TRUE, skip runs whose output files already exist.
#' @param write_vcov              If TRUE, write full vcov matrices per run.
#'                                Default FALSE.
#' @param verbose_timing          If TRUE, print timing for each analysis subset.
#' @param .progress               Passed to purrr::pmap for a progress bar.
#'
#' @return Invisibly: a list with run_grid, run_results, and metadata_files.
run_experiment <- function(dataset_spec,
                           analysis_subset_specs,
                           outcome_specs,
                           treatment_group_specs,
                           model_specs,
                           vcov_specs,
                           dir_out,
                           group_palette      = NULL,
                           ci_level           = 0.95,
                           run_estimation     = TRUE,
                           run_descriptive    = FALSE,
                           descriptive_args   = list(),
                           skip_existing      = TRUE,
                           write_vcov         = FALSE,
                           verbose_timing     = FALSE,
                           .progress          = FALSE) {
  
  # ── Validate inputs ────────────────────────────────────────────────────────
  
  validate_dataset_spec(dataset_spec)
  validate_spec_table(analysis_subset_specs,  c("subset_id", "data_source", "data_filter"), "analysis_subset_specs")
  validate_spec_table(outcome_specs,           c("outcome"),                                 "outcome_specs")
  validate_spec_table(treatment_group_specs,   c("group_id", "group_col", "group_fun", "group_args"), "treatment_group_specs")
  validate_spec_table(model_specs,             c("model_id", "formula_template"),            "model_specs")
  validate_spec_table(vcov_specs,              c("vcov_id", "vcov", "vcov_label"),           "vcov_specs")
  
  if (!"vcov_vars" %in% names(vcov_specs)) {
    vcov_specs$vcov_vars <- replicate(nrow(vcov_specs), character(0), simplify = FALSE)
  }
  
  if (!run_estimation && !run_descriptive) {
    stop("At least one of run_estimation or run_descriptive must be TRUE.")
  }
  
  # ── Resolve descriptive args defaults ─────────────────────────────────────
  
  desc_args <- list(
    treated_year_var  = dataset_spec$time_var,
    control_year_var  = dataset_spec$time_var
  )
  desc_args[names(descriptive_args)] <- descriptive_args
  
  # ── Build run grid ─────────────────────────────────────────────────────────
  
  dir_ensure_local(c(
    dir_out,
    file.path(dir_out, "tables", "by_run"),
    file.path(dir_out, "descriptive", "by_run"),
    file.path(dir_out, "metadata")
  ))
  
  run_grid <- tidyr::crossing(
    analysis_subset_specs,
    outcome_specs,
    treatment_group_specs,
    model_specs
  ) |>
    dplyr::mutate(
      run_id = glue::glue(
        "{subset_id}__{outcome}__{group_id}__{model_id}"
      ),
      run_stub = purrr::pmap_chr(
        list(subset_id, outcome, group_id, model_id),
        \(s, o, g, m) make_run_stub(s, o, g, m)
      )
    )
  
  ts <- timestamp_now()
  
  run_grid_file <- file.path(dir_out, "metadata",
                             glue::glue("run_grid__{ts}.rds"))
  specs_file    <- file.path(dir_out, "metadata",
                             glue::glue("specs_snapshot__{ts}.rds"))
  info_file     <- file.path(dir_out, "metadata",
                             glue::glue("experiment_info__{ts}.rds"))
  
  saveRDS(run_grid, run_grid_file)
  saveRDS(
    list(
      dataset_spec           = dataset_spec,
      analysis_subset_specs  = analysis_subset_specs,
      outcome_specs          = outcome_specs,
      treatment_group_specs  = treatment_group_specs,
      model_specs            = model_specs,
      vcov_specs             = vcov_specs,
      group_palette          = group_palette,
      ci_level               = ci_level,
      run_estimation         = run_estimation,
      run_descriptive        = run_descriptive,
      descriptive_args       = desc_args,
      skip_existing          = skip_existing,
      write_vcov             = write_vcov
    ),
    specs_file
  )
  
  # ── Execute runs ───────────────────────────────────────────────────────────
  
  if (verbose_timing) {
    message(glue::glue(
      "[{timestamp_now()}] Starting experiment: {nrow(run_grid)} runs planned."
    ))
  }
  
  run_results <- purrr::pmap(
    list(
      data_source  = run_grid$data_source,
      data_filter  = run_grid$data_filter,
      subset_id    = run_grid$subset_id,
      outcome      = run_grid$outcome,
      group_id     = run_grid$group_id,
      group_col    = run_grid$group_col,
      group_fun    = run_grid$group_fun,
      group_args   = run_grid$group_args,
      model_id     = run_grid$model_id,
      formula_template = run_grid$formula_template,
      estimator_type   = run_grid$estimator_type,
      term_pattern     = run_grid$term_pattern,
      run_id       = run_grid$run_id,
      run_stub     = run_grid$run_stub
    ),
    function(data_source, data_filter, subset_id, outcome,
             group_id, group_col, group_fun, group_args,
             model_id, formula_template, estimator_type, term_pattern,
             run_id, run_stub) {
      
      result <- list(
        run_id      = run_id,
        estimation  = NULL,
        descriptive = NULL,
        error       = NULL,
        skipped     = FALSE
      )
      
      t_start <- proc.time()
      
      tryCatch({
        
        if (run_estimation) {
          result$estimation <- run_one_estimation(
            data_source      = data_source,
            data_filter      = data_filter,
            dataset_spec     = dataset_spec,
            subset_id        = subset_id,
            outcome          = outcome,
            group_id         = group_id,
            group_col        = group_col,
            group_fun        = group_fun,
            group_args       = group_args,
            model_id         = model_id,
            formula_template = formula_template,
            estimator_type   = estimator_type,
            term_pattern     = term_pattern,
            vcov_specs       = vcov_specs,
            dir_out          = dir_out,
            run_id           = run_id,
            run_stub         = run_stub,
            group_palette    = group_palette,
            ci_level         = ci_level,
            skip_existing    = skip_existing,
            write_vcov       = write_vcov
          )
          result$skipped <- isTRUE(result$estimation$skipped_existing)
        }
        
        if (run_descriptive) {
          result$descriptive <- run_one_descriptive(
            data_source       = data_source,
            data_filter       = data_filter,
            dataset_spec      = dataset_spec,
            subset_id         = subset_id,
            outcome           = outcome,
            group_id          = group_id,
            group_col         = group_col,
            group_fun         = group_fun,
            group_args        = group_args,
            model_id          = model_id,
            formula_template  = formula_template,
            dir_out           = dir_out,
            run_id            = run_id,
            run_stub          = run_stub,
            group_palette     = group_palette,
            treated_year_var  = desc_args$treated_year_var,
            control_year_var  = desc_args$control_year_var,
            skip_existing     = skip_existing
          )
        }
        
      }, error = function(e) {
        result$error <<- conditionMessage(e)
        message(glue::glue("[ERROR] run_id = {run_id}: {conditionMessage(e)}"))
      })
      
      if (verbose_timing) {
        elapsed <- (proc.time() - t_start)[["elapsed"]]
        status  <- if (!is.null(result$error)) "FAILED" else if (result$skipped) "SKIPPED" else "OK"
        message(glue::glue(
          "[{timestamp_now()}] {run_id} | {status} | {round(elapsed, 1)}s"
        ))
      }
      
      result
    },
    .progress = .progress
  )
  
  # ── Log failures ───────────────────────────────────────────────────────────
  
  failed <- purrr::keep(run_results, \(r) !is.null(r$error))
  if (length(failed) > 0) {
    message(glue::glue(
      "\n{length(failed)} run(s) failed:\n",
      paste(purrr::map_chr(failed, "run_id"), collapse = "\n")
    ))
  }
  
  # ── Save experiment info ───────────────────────────────────────────────────
  
  experiment_info <- list(
    dir_out         = dir_out,
    timestamp       = ts,
    n_runs_planned  = nrow(run_grid),
    n_runs_ok       = sum(purrr::map_lgl(run_results, \(r) is.null(r$error))),
    n_runs_failed   = length(failed),
    n_runs_skipped  = sum(purrr::map_lgl(run_results, \(r) isTRUE(r$skipped))),
    failed_run_ids  = purrr::map_chr(failed, "run_id"),
    metadata_files  = list(
      run_grid_rds      = run_grid_file,
      specs_snapshot    = specs_file
    )
  )
  saveRDS(experiment_info, info_file)
  
  invisible(list(
    run_grid       = run_grid,
    run_results    = run_results,
    metadata_files = list(
      run_grid_rds   = run_grid_file,
      specs_snapshot = specs_file,
      experiment_info = info_file
    )
  ))
}


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 6 — Single-run workers
# ══════════════════════════════════════════════════════════════════════════════

# ── Estimation ────────────────────────────────────────────────────────────────

run_one_estimation <- function(data_source,
                               data_filter,
                               dataset_spec,
                               subset_id,
                               outcome,
                               group_id,
                               group_col,
                               group_fun,
                               group_args,
                               model_id,
                               formula_template,
                               estimator_type,
                               term_pattern,
                               vcov_specs,
                               dir_out,
                               run_id,
                               run_stub,
                               group_palette  = NULL,
                               ci_level       = 0.95,
                               skip_existing  = TRUE,
                               write_vcov     = FALSE) {
  
  run_dir <- file.path(dir_out, "tables", "by_run", run_stub)
  dir_ensure_local(run_dir)
  
  coef_file      <- file.path(run_dir, "coef.parquet")
  subgroup_file  <- file.path(run_dir, "subgroup.parquet")
  support_file   <- file.path(run_dir, "support.parquet")
  vcov_file      <- file.path(run_dir, "vcov.parquet")
  registry_file  <- file.path(run_dir, "registry.parquet")
  run_spec_file  <- file.path(run_dir, "run_spec.rds")
  
  required_files <- c(coef_file, subgroup_file, support_file, registry_file, run_spec_file)
  if (write_vcov) required_files <- c(required_files, vcov_file)
  
  if (skip_existing && all(file.exists(required_files))) {
    return(list(
      coef_file      = coef_file,
      subgroup_file  = subgroup_file,
      support_file   = support_file,
      vcov_file      = if (write_vcov) vcov_file else NULL,
      registry_file  = registry_file,
      run_spec_file  = run_spec_file,
      skipped_existing = TRUE
    ))
  }
  
  fit_started <- Sys.time()
  
  # ── Load data ──────────────────────────────────────────────────────────────
  
  df <- load_arrow_data(data_source, data_filter)
  n_rows_read        <- nrow(df)
  n_rows_after_filter <- nrow(df)   # filter already pushed to Arrow
  
  # ── Apply grouping ─────────────────────────────────────────────────────────
  
  df_grouped <- apply_treatment_grouping(
    df        = df,
    group_fun = group_fun,
    group_col = group_col,
    group_args = group_args
  )
  
  # ── Build formula ──────────────────────────────────────────────────────────
  
  model_formula <- build_model_formula(formula_template, outcome)
  
  # ── Identify columns needed ────────────────────────────────────────────────
  
  needed_cols <- get_needed_columns(
    formula     = model_formula,
    trt_col     = dataset_spec$trt_col,
    group_col   = group_col,
    unit_id     = dataset_spec$unit_id,
    event_id    = dataset_spec$event_id,
    vcov_vars   = vcov_specs$vcov_vars
  )
  
  check_required_columns(df_grouped, needed_cols, context = run_id)
  
  df_grouped <- df_grouped |> dplyr::select(dplyr::all_of(needed_cols))
  
  # ── Fit subgroup models ────────────────────────────────────────────────────
  
  subgroup_models <- fit_subgroup_models(
    df          = df_grouped,
    group_col   = group_col,
    dataset_spec = dataset_spec,
    formula     = model_formula,
    estimator_type = estimator_type
  )
  
  fit_finished <- Sys.time()
  
  # ── Extract tables ─────────────────────────────────────────────────────────
  
  shared_meta <- list(
    subset_id        = subset_id,
    outcome          = outcome,
    group_id         = group_id,
    group_col        = group_col,
    model_id         = model_id,
    formula_template = formula_template,
    estimator_type   = estimator_type,
    group_palette    = group_palette
  )
  
  subgroup_tbl <- make_subgroup_run_summary(
    models = subgroup_models, meta = shared_meta
  ) |> dplyr::mutate(run_id = run_id)
  
  support_tbl <- make_event_time_support_run_summary(
    models = subgroup_models, meta = shared_meta
  ) |> dplyr::mutate(run_id = run_id)
  
  coef_tbl <- extract_all_coef_tables(
    models     = subgroup_models,
    meta       = shared_meta,
    vcov_specs = vcov_specs,
    term_pattern = term_pattern,
    ci_level   = ci_level
  ) |> dplyr::mutate(run_id = run_id)
  
  vcov_tbl <- if (write_vcov) {
    extract_all_vcov_tables(
      models     = subgroup_models,
      meta       = shared_meta,
      vcov_specs = vcov_specs
    ) |> dplyr::mutate(run_id = run_id)
  } else {
    NULL
  }
  
  registry_tbl <- make_estimation_registry(
    run_id           = run_id,
    subset_id        = subset_id,
    outcome          = outcome,
    group_id         = group_id,
    model_id         = model_id,
    group_col        = group_col,
    dataset_spec     = dataset_spec,
    formula_template = formula_template,
    estimator_type   = estimator_type,
    term_pattern     = term_pattern,
    data_filter      = data_filter,
    data_source      = data_source,
    vcov_specs       = vcov_specs,
    group_args       = group_args,
    n_rows_read      = n_rows_read,
    n_rows_after_filter = n_rows_after_filter,
    subgroup_tbl     = subgroup_tbl,
    support_tbl      = support_tbl,
    fit_started      = fit_started,
    fit_finished     = fit_finished,
    coef_file        = coef_file,
    subgroup_file    = subgroup_file,
    support_file     = support_file,
    vcov_file        = if (write_vcov) vcov_file else NA_character_,
    registry_file    = registry_file,
    run_spec_file    = run_spec_file
  )
  
  run_spec <- list(
    run_id           = run_id,
    run_stub         = run_stub,
    subset_id        = subset_id,
    outcome          = outcome,
    group_id         = group_id,
    model_id         = model_id,
    group_col        = group_col,
    dataset_spec     = dataset_spec,
    formula_template = formula_template,
    estimator_type   = estimator_type,
    term_pattern     = term_pattern,
    data_filter      = data_filter,
    data_source      = data_source,
    group_fun        = group_fun,
    group_args       = group_args,
    vcov_specs       = vcov_specs,
    group_palette    = group_palette,
    ci_level         = ci_level,
    write_vcov       = write_vcov
  )
  
  # ── Write outputs ──────────────────────────────────────────────────────────
  
  arrow::write_parquet(subgroup_tbl, subgroup_file)
  arrow::write_parquet(support_tbl,  support_file)
  arrow::write_parquet(coef_tbl,     coef_file)
  arrow::write_parquet(registry_tbl, registry_file)
  if (write_vcov && !is.null(vcov_tbl)) {
    arrow::write_parquet(vcov_tbl, vcov_file)
  }
  saveRDS(run_spec, run_spec_file)
  
  rm(subgroup_models, df, df_grouped,
     subgroup_tbl, support_tbl, coef_tbl, vcov_tbl, registry_tbl, run_spec)
  gc()
  
  list(
    coef_file      = coef_file,
    subgroup_file  = subgroup_file,
    support_file   = support_file,
    vcov_file      = if (write_vcov) vcov_file else NULL,
    registry_file  = registry_file,
    run_spec_file  = run_spec_file,
    skipped_existing = FALSE
  )
}


# ── Descriptive ───────────────────────────────────────────────────────────────

run_one_descriptive <- function(data_source,
                                data_filter,
                                dataset_spec,
                                subset_id,
                                outcome,
                                group_id,
                                group_col,
                                group_fun,
                                group_args,
                                model_id,
                                formula_template,
                                dir_out,
                                run_id,
                                run_stub,
                                group_palette    = NULL,
                                treated_year_var,
                                control_year_var,
                                skip_existing    = TRUE) {
  
  run_dir <- file.path(dir_out, "descriptive", "by_run", run_stub)
  dir_ensure_local(run_dir)
  
  trajectory_file <- file.path(run_dir, "event_time_trajectory.parquet")
  registry_file   <- file.path(run_dir, "registry.parquet")
  run_spec_file   <- file.path(run_dir, "descriptive_spec.rds")
  
  if (skip_existing &&
      file.exists(trajectory_file) &&
      file.exists(registry_file) &&
      file.exists(run_spec_file)) {
    return(list(
      trajectory_file  = trajectory_file,
      registry_file    = registry_file,
      run_spec_file    = run_spec_file,
      skipped_existing = TRUE
    ))
  }
  
  run_started <- Sys.time()
  
  df <- load_arrow_data(data_source, data_filter)
  n_rows_read         <- nrow(df)
  n_rows_after_filter <- nrow(df)
  
  df_grouped <- apply_treatment_grouping(
    df         = df,
    group_fun  = group_fun,
    group_col  = group_col,
    group_args = group_args
  )
  
  shared_meta <- list(
    subset_id        = subset_id,
    outcome          = outcome,
    group_id         = group_id,
    group_col        = group_col,
    model_id         = model_id,
    formula_template = formula_template,
    group_palette    = group_palette
  )
  
  trajectory_tbl <- make_descriptive_trajectory(
    df               = df_grouped,
    meta             = shared_meta,
    dataset_spec     = dataset_spec,
    treated_year_var = treated_year_var,
    control_year_var = control_year_var
  ) |> dplyr::mutate(run_id = run_id)
  
  run_finished <- Sys.time()
  
  registry_tbl <- tibble::tibble(
    run_id              = run_id,
    subset_id           = subset_id,
    outcome             = outcome,
    group_id            = group_id,
    model_id            = model_id,
    group_col           = group_col,
    unit_id             = dataset_spec$unit_id,
    event_id            = dataset_spec$event_id,
    trt_col             = dataset_spec$trt_col,
    formula_template    = formula_template,
    data_filter         = data_filter_to_chr(data_filter),
    input_source        = paste(unlist(data_source), collapse = " | "),
    treated_year_var    = treated_year_var,
    control_year_var    = control_year_var,
    time_var            = dataset_spec$time_var,
    n_rows_read         = n_rows_read,
    n_rows_after_filter = n_rows_after_filter,
    n_subgroups         = dplyr::n_distinct(trajectory_tbl$subgroup),
    subgroup_levels     = paste(sort(unique(trajectory_tbl$subgroup)), collapse = ","),
    run_started         = run_started,
    run_finished        = run_finished,
    trajectory_file     = trajectory_file,
    registry_file       = registry_file,
    run_spec_file       = run_spec_file
  )
  
  run_spec <- list(
    run_id           = run_id,
    run_stub         = run_stub,
    subset_id        = subset_id,
    outcome          = outcome,
    group_id         = group_id,
    model_id         = model_id,
    group_col        = group_col,
    dataset_spec     = dataset_spec,
    formula_template = formula_template,
    data_filter      = data_filter,
    data_source      = data_source,
    group_fun        = group_fun,
    group_args       = group_args,
    group_palette    = group_palette,
    treated_year_var = treated_year_var,
    control_year_var = control_year_var
  )
  
  arrow::write_parquet(trajectory_tbl, trajectory_file)
  arrow::write_parquet(registry_tbl,   registry_file)
  saveRDS(run_spec, run_spec_file)
  
  rm(df, df_grouped, trajectory_tbl, registry_tbl, run_spec)
  gc()
  
  list(
    trajectory_file  = trajectory_file,
    registry_file    = registry_file,
    run_spec_file    = run_spec_file,
    skipped_existing = FALSE
  )
}


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 7 — Rebuild (merge across runs)
# ══════════════════════════════════════════════════════════════════════════════

#' Rebuild estimation tables by merging all by_run outputs
#'
#' Scans dir_out/tables/by_run recursively, stack-binds all per-run parquets,
#' deduplicates on run_id (keeping the most recently written file), and writes
#' merged outputs to dir_out/tables/all/.
#'
#' @param dir_out    Root output directory.
#' @param write_csv  Also write CSV versions? Default TRUE.
#' @param recursive  Scan subdirectories? Default TRUE.
#'
#' @return Invisibly: a list of merged tibbles and output file paths.
rebuild_estimation_tables <- function(dir_out,
                                      write_csv = TRUE,
                                      recursive = TRUE) {
  
  dir_by_run <- file.path(dir_out, "tables", "by_run")
  dir_all    <- file.path(dir_out, "tables", "all")
  
  if (!dir.exists(dir_by_run)) {
    stop("by_run directory does not exist: ", dir_by_run)
  }
  dir.create(dir_all, recursive = TRUE, showWarnings = FALSE)
  
  find_parquets <- function(pattern) {
    list.files(dir_by_run, pattern = pattern,
               recursive = recursive, full.names = TRUE)
  }
  
  coef_files     <- find_parquets("^coef\\.parquet$")
  subgroup_files <- find_parquets("^subgroup\\.parquet$")
  support_files  <- find_parquets("^support\\.parquet$")
  registry_files <- find_parquets("^registry\\.parquet$")
  run_spec_files <- find_parquets("^run_spec\\.rds$") # pattern still uses list.files
  run_spec_files <- list.files(dir_by_run, pattern = "^run_spec\\.rds$",
                               recursive = recursive, full.names = TRUE)
  vcov_files     <- find_parquets("^vcov\\.parquet$")
  
  for (nm in c("coef_files", "subgroup_files", "support_files", "registry_files")) {
    if (length(get(nm)) == 0) {
      stop("No ", sub("_files", ".parquet", nm), " files found under: ", dir_by_run)
    }
  }
  
  read_and_dedup <- function(files, id_col = "run_id") {
    purrr::map_dfr(
      files,
      \(f) arrow::read_parquet(f) |>
        dplyr::mutate(.source_mtime = file.mtime(f))
    ) |>
      dplyr::group_by(dplyr::across(dplyr::all_of(id_col))) |>
      dplyr::filter(.source_mtime == max(.source_mtime)) |>
      dplyr::ungroup() |>
      dplyr::select(-.source_mtime)
  }
  
  # For coef and vcov, the dedup key includes vcov_id so we keep one row per
  # run × term × vcov combination from the most recent run file
  coef_tbl     <- read_and_dedup(coef_files,     id_col = c("run_id", "subgroup", "term", "vcov_id"))
  subgroup_tbl <- read_and_dedup(subgroup_files,  id_col = c("run_id", "subgroup"))
  support_tbl  <- read_and_dedup(support_files,   id_col = c("run_id", "subgroup", "event_time"))
  registry_tbl <- read_and_dedup(registry_files,  id_col = "run_id")
  
  vcov_tbl <- if (length(vcov_files) > 0) {
    read_and_dedup(vcov_files, id_col = c("run_id", "subgroup", "term_i", "term_j", "vcov_id"))
  } else {
    NULL
  }
  
  # ── run_spec index ──────────────────────────────────────────────────────────
  
  spec_index_tbl <- NULL
  if (length(run_spec_files) > 0) {
    spec_index_tbl <- purrr::map_dfr(
      run_spec_files,
      \(f) {
        spec <- readRDS(f)
        tibble::tibble(
          run_id           = spec$run_id %||% NA_character_,
          subset_id        = spec$subset_id %||% NA_character_,
          outcome          = spec$outcome %||% NA_character_,
          group_id         = spec$group_id %||% NA_character_,
          model_id         = spec$model_id %||% NA_character_,
          group_col        = spec$group_col %||% NA_character_,
          estimator_type   = spec$estimator_type %||% NA_character_,
          formula_template = spec$formula_template %||% NA_character_,
          term_pattern     = spec$term_pattern %||% NA_character_,
          data_filter      = data_filter_to_chr(spec$data_filter),
          data_source      = paste(unlist(spec$data_source), collapse = " | "),
          group_args_json  = as.character(serialize_object_json(spec$group_args)),
          vcov_ids         = paste(spec$vcov_specs$vcov_id, collapse = " | "),
          run_spec_file    = f,
          .source_mtime    = file.mtime(f)
        )
      }
    ) |>
      dplyr::group_by(run_id) |>
      dplyr::filter(.source_mtime == max(.source_mtime)) |>
      dplyr::ungroup() |>
      dplyr::select(-.source_mtime)
  }
  
  # ── Write outputs ───────────────────────────────────────────────────────────
  
  write_outputs <- function(tbl, stem) {
    if (is.null(tbl)) return(list(parquet = NULL, csv = NULL))
    pq  <- file.path(dir_all, paste0(stem, ".parquet"))
    csv <- if (write_csv) file.path(dir_all, paste0(stem, ".csv")) else NULL
    arrow::write_parquet(tbl, pq)
    if (write_csv) readr::write_csv(tbl, csv)
    list(parquet = pq, csv = csv)
  }
  
  out_coef     <- write_outputs(coef_tbl,       "coef_expanded_vcov")
  out_subgroup <- write_outputs(subgroup_tbl,    "subgroup_run_summary")
  out_support  <- write_outputs(support_tbl,     "event_time_support")
  out_registry <- write_outputs(registry_tbl,    "run_registry")
  out_vcov     <- write_outputs(vcov_tbl,        "vcov_matrices")
  out_spec_idx <- write_outputs(spec_index_tbl,  "run_spec_index")
  
  invisible(list(
    coef_expanded_vcov   = coef_tbl,
    subgroup_run_summary = subgroup_tbl,
    event_time_support   = support_tbl,
    run_registry         = registry_tbl,
    vcov_matrices        = vcov_tbl,
    run_spec_index       = spec_index_tbl,
    files = list(
      coef_expanded_vcov   = out_coef,
      subgroup_run_summary = out_subgroup,
      event_time_support   = out_support,
      run_registry         = out_registry,
      vcov_matrices        = out_vcov,
      run_spec_index       = out_spec_idx
    )
  ))
}


#' Rebuild descriptive tables by merging all by_run outputs
#'
#' @inheritParams rebuild_estimation_tables
rebuild_descriptive_tables <- function(dir_out,
                                       write_csv = TRUE,
                                       recursive = TRUE) {
  
  dir_by_run <- file.path(dir_out, "descriptive", "by_run")
  dir_all    <- file.path(dir_out, "descriptive", "all")
  
  if (!dir.exists(dir_by_run)) {
    stop("descriptive by_run directory does not exist: ", dir_by_run)
  }
  dir.create(dir_all, recursive = TRUE, showWarnings = FALSE)
  
  traj_files     <- list.files(dir_by_run, "^event_time_trajectory\\.parquet$",
                               recursive = recursive, full.names = TRUE)
  registry_files <- list.files(dir_by_run, "^registry\\.parquet$",
                               recursive = recursive, full.names = TRUE)
  
  if (length(traj_files) == 0)     stop("No event_time_trajectory.parquet files found under: ", dir_by_run)
  if (length(registry_files) == 0) stop("No registry.parquet files found under: ", dir_by_run)
  
  read_and_dedup <- function(files, id_col) {
    purrr::map_dfr(
      files,
      \(f) arrow::read_parquet(f) |>
        dplyr::mutate(.source_mtime = file.mtime(f))
    ) |>
      dplyr::group_by(dplyr::across(dplyr::all_of(id_col))) |>
      dplyr::filter(.source_mtime == max(.source_mtime)) |>
      dplyr::ungroup() |>
      dplyr::select(-.source_mtime)
  }
  
  traj_tbl     <- read_and_dedup(traj_files,     c("run_id", "subgroup", "series", "event_time"))
  registry_tbl <- read_and_dedup(registry_files,  "run_id")
  
  write_pair <- function(tbl, stem) {
    pq  <- file.path(dir_all, paste0(stem, ".parquet"))
    csv <- if (write_csv) file.path(dir_all, paste0(stem, ".csv")) else NULL
    arrow::write_parquet(tbl, pq)
    if (write_csv) readr::write_csv(tbl, csv)
    list(parquet = pq, csv = csv)
  }
  
  out_traj <- write_pair(traj_tbl,     "event_time_trajectory")
  out_reg  <- write_pair(registry_tbl, "run_registry")
  
  invisible(list(
    event_time_trajectory = traj_tbl,
    run_registry          = registry_tbl,
    files = list(
      event_time_trajectory = out_traj,
      run_registry          = out_reg
    )
  ))
}


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 8 — Model fitting internals
# ══════════════════════════════════════════════════════════════════════════════

fit_subgroup_models <- function(df,
                                group_col,
                                dataset_spec,
                                formula,
                                estimator_type,
                                ...) {
  
  trt_col  <- dataset_spec$trt_col
  unit_id  <- dataset_spec$unit_id
  event_id <- dataset_spec$event_id
  
  groups <- df |>
    dplyr::filter(!is.na(.data[[group_col]])) |>
    dplyr::distinct(.data[[group_col]]) |>
    dplyr::pull(.data[[group_col]]) |>
    as.character() |>
    sort()
  
  if (length(groups) == 0) {
    stop(glue::glue("No non-NA values found in group column '{group_col}'."))
  }
  
  models <- purrr::map(
    groups,
    \(g) fit_one_subgroup_model(
      df             = df,
      group_col      = group_col,
      group_value    = g,
      trt_col        = trt_col,
      formula        = formula,
      unit_id        = unit_id,
      event_id       = event_id,
      estimator_type = estimator_type,
      ...
    )
  )
  
  names(models) <- groups
  models
}


fit_one_subgroup_model <- function(df,
                                   group_col,
                                   group_value,
                                   trt_col,
                                   formula,
                                   unit_id,
                                   event_id,
                                   estimator_type,
                                   ...) {
  
  group_sym <- rlang::sym(group_col)
  trt_sym   <- rlang::sym(trt_col)
  
  df_controls <- df |> dplyr::filter(!!trt_sym == 0)
  df_treated  <- df |> dplyr::filter(!!trt_sym == 1, !!group_sym == group_value)
  df_model    <- dplyr::bind_rows(df_controls, df_treated)
  
  model <- fixest::feols(formula, data = df_model, ...)
  
  obs_used <- fixest::obs(model)
  df_est   <- df_model[obs_used, , drop = FALSE]
  
  # Attach metadata attributes
  attr(model, "group_value")    <- group_value
  attr(model, "group_col")      <- group_col
  attr(model, "unit_id")        <- unit_id
  attr(model, "event_id")       <- normalize_optional_colname(event_id)
  attr(model, "estimator_type") <- estimator_type
  
  # Sun-Abraham specific post-processing
  if (estimator_type == "sunab") {
    sunab_vars <- parse_sunab_vars(formula)
    attr(model, "cohort_var") <- sunab_vars$cohort_var
    attr(model, "time_var")   <- sunab_vars$time_var
    
    event_time_support <- make_event_time_support_one_model(
      df_est     = df_est,
      trt_col    = trt_col,
      unit_id    = unit_id,
      event_id   = event_id,
      cohort_var = sunab_vars$cohort_var,
      time_var   = sunab_vars$time_var
    )
  } else {
    attr(model, "cohort_var") <- NA_character_
    attr(model, "time_var")   <- NA_character_
    event_time_support        <- empty_event_time_support()
  }
  
  attr(model, "event_time_support") <- event_time_support
  
  # Sample size attributes
  attr(model, "n_treated_units") <- df_est |>
    dplyr::filter(.data[[trt_col]] == 1) |>
    dplyr::summarise(n = dplyr::n_distinct(.data[[unit_id]])) |>
    dplyr::pull(n)
  
  attr(model, "n_treated_events") <- df_est |>
    dplyr::filter(.data[[trt_col]] == 1) |>
    (\(x) compute_n_distinct_optional(x, event_id))()
  
  attr(model, "n_control_units") <- df_est |>
    dplyr::filter(.data[[trt_col]] == 0) |>
    dplyr::summarise(n = dplyr::n_distinct(.data[[unit_id]])) |>
    dplyr::pull(n)
  
  attr(model, "n_total_units")     <- dplyr::n_distinct(df_est[[unit_id]])
  attr(model, "n_rows_model_data") <- nrow(df_est)
  
  model
}


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 9 — Table extraction internals
# ══════════════════════════════════════════════════════════════════════════════

extract_all_coef_tables <- function(models,
                                    meta,
                                    vcov_specs,
                                    term_pattern,
                                    ci_level) {
  purrr::imap_dfr(
    models,
    \(model, group_nm) {
      purrr::pmap_dfr(
        list(
          vcov       = vcov_specs$vcov,
          vcov_id    = vcov_specs$vcov_id,
          vcov_label = vcov_specs$vcov_label
        ),
        \(vcov, vcov_id, vcov_label) {
          extract_coef_table_one_vcov(
            model        = model,
            vcov         = vcov,
            vcov_id      = vcov_id,
            vcov_label   = vcov_label,
            meta         = meta,
            term_pattern = term_pattern,
            ci_level     = ci_level
          )
        }
      )
    }
  )
}


extract_coef_table_one_vcov <- function(model,
                                        vcov,
                                        vcov_id,
                                        vcov_label,
                                        meta,
                                        term_pattern,
                                        ci_level) {
  
  ct      <- fixest::coeftable(model, vcov = vcov)
  ct_tbl  <- tibble::as_tibble(as.data.frame(ct), rownames = "term")
  
  se_col <- if ("Std. Error" %in% names(ct_tbl)) "Std. Error" else names(ct_tbl)[3]
  t_col  <- names(ct_tbl)[stringr::str_detect(names(ct_tbl), "^t value$|^t-value$")]
  p_col  <- names(ct_tbl)[stringr::str_detect(names(ct_tbl), "^Pr\\(>\\|t\\|\\)$|^p-value$")]
  z_val  <- stats::qnorm(1 - (1 - ci_level) / 2)
  
  group_value      <- attr(model, "group_value")
  event_time_support <- attr(model, "event_time_support")
  
  subgroup_color <- get_group_color(group_value, meta$group_palette)
  
  subgroup_run_id <- make_subgroup_run_id(meta, group_value)
  coef_run_id     <- glue::glue("{subgroup_run_id}__{vcov_id}")
  
  out <- ct_tbl |>
    dplyr::mutate(
      subset_id        = meta$subset_id,
      outcome          = meta$outcome,
      group_id         = meta$group_id,
      group_col        = meta$group_col,
      subgroup         = group_value,
      model_id         = meta$model_id,
      estimator_type   = meta$estimator_type,
      vcov_id          = vcov_id,
      vcov_label       = vcov_label,
      formula_template = meta$formula_template,
      estimate         = .data[["Estimate"]],
      std_error        = .data[[se_col]],
      ci_lower         = estimate - z_val * std_error,
      ci_upper         = estimate + z_val * std_error,
      term_matches_pattern = stringr::str_detect(term, term_pattern),
      term_value       = extract_first_number(term),
      n_treated_units  = attr(model, "n_treated_units"),
      n_treated_events = attr(model, "n_treated_events"),
      n_control_units  = attr(model, "n_control_units"),
      n_total_units    = attr(model, "n_total_units"),
      n_rows_model_data = attr(model, "n_rows_model_data"),
      subgroup_color   = subgroup_color,
      subgroup_run_id  = subgroup_run_id,
      coef_run_id      = coef_run_id
    )
  
  if (length(t_col) == 1) out$t_value <- ct_tbl[[t_col]] else out$t_value <- NA_real_
  if (length(p_col) == 1) out$p_value <- ct_tbl[[p_col]] else out$p_value <- NA_real_
  
  if (!is.null(event_time_support) && nrow(event_time_support) > 0) {
    out <- out |>
      dplyr::left_join(
        event_time_support |> dplyr::rename(term_value = event_time),
        by = "term_value"
      )
  } else {
    out <- out |>
      dplyr::mutate(
        n_ptids        = NA_integer_,
        n_fireids      = NA_integer_,
        n_rows_treated = NA_integer_
      )
  }
  
  out |>
    dplyr::select(
      subset_id, outcome, group_id, group_col, subgroup, model_id,
      estimator_type, vcov_id, vcov_label, formula_template,
      term, estimate, std_error, t_value, p_value, ci_lower, ci_upper,
      term_matches_pattern, term_value,
      n_ptids, n_fireids, n_rows_treated,
      n_treated_units, n_treated_events, n_control_units,
      n_total_units, n_rows_model_data,
      subgroup_color, subgroup_run_id, coef_run_id
    )
}


extract_all_vcov_tables <- function(models, meta, vcov_specs) {
  purrr::imap_dfr(
    models,
    \(model, group_nm) {
      purrr::pmap_dfr(
        list(
          vcov       = vcov_specs$vcov,
          vcov_id    = vcov_specs$vcov_id,
          vcov_label = vcov_specs$vcov_label
        ),
        \(vcov, vcov_id, vcov_label) {
          extract_vcov_table_one_vcov(
            model      = model,
            vcov       = vcov,
            vcov_id    = vcov_id,
            vcov_label = vcov_label,
            meta       = meta
          )
        }
      )
    }
  )
}


extract_vcov_table_one_vcov <- function(model, vcov, vcov_id, vcov_label, meta) {
  vc     <- stats::vcov(model, vcov = vcov)
  vc_mat <- as.matrix(vc)
  
  if (is.null(rownames(vc_mat)) || is.null(colnames(vc_mat))) {
    stop("vcov matrix must have row and column names.")
  }
  
  group_value     <- attr(model, "group_value")
  subgroup_run_id <- make_subgroup_run_id(meta, group_value)
  
  expand.grid(
    term_i = rownames(vc_mat),
    term_j = colnames(vc_mat),
    stringsAsFactors = FALSE
  ) |>
    tibble::as_tibble() |>
    dplyr::mutate(
      subset_id        = meta$subset_id,
      outcome          = meta$outcome,
      group_id         = meta$group_id,
      group_col        = meta$group_col,
      subgroup         = group_value,
      model_id         = meta$model_id,
      vcov_id          = vcov_id,
      vcov_label       = vcov_label,
      formula_template = meta$formula_template,
      vcov_value       = purrr::map2_dbl(term_i, term_j, \(i, j) vc_mat[i, j]),
      subgroup_run_id  = subgroup_run_id,
      vcov_run_id      = glue::glue("{subgroup_run_id}__{vcov_id}")
    ) |>
    dplyr::select(
      subset_id, outcome, group_id, group_col, subgroup, model_id,
      vcov_id, vcov_label, formula_template,
      term_i, term_j, vcov_value, subgroup_run_id, vcov_run_id
    )
}


make_subgroup_run_summary <- function(models, meta) {
  purrr::imap_dfr(
    models,
    \(model, group_nm) {
      tibble::tibble(
        subset_id        = meta$subset_id,
        outcome          = meta$outcome,
        group_id         = meta$group_id,
        group_col        = meta$group_col,
        subgroup         = group_nm,
        model_id         = meta$model_id,
        formula_template = meta$formula_template,
        estimator_type   = meta$estimator_type,
        n_treated_units  = attr(model, "n_treated_units"),
        n_treated_events = attr(model, "n_treated_events"),
        n_control_units  = attr(model, "n_control_units"),
        n_total_units    = attr(model, "n_total_units"),
        n_rows_model_data = attr(model, "n_rows_model_data"),
        subgroup_color   = get_group_color(group_nm, meta$group_palette),
        subgroup_run_id  = make_subgroup_run_id(meta, group_nm)
      )
    }
  )
}


make_event_time_support_run_summary <- function(models, meta) {
  purrr::imap_dfr(
    models,
    \(model, group_nm) {
      support_tbl <- attr(model, "event_time_support")
      cohort_var  <- attr(model, "cohort_var")
      time_var    <- attr(model, "time_var")
      
      base_cols <- list(
        subset_id        = meta$subset_id,
        outcome          = meta$outcome,
        group_id         = meta$group_id,
        group_col        = meta$group_col,
        subgroup         = group_nm,
        model_id         = meta$model_id,
        formula_template = meta$formula_template,
        estimator_type   = meta$estimator_type,
        cohort_var       = normalize_optional_colname(cohort_var),
        time_var         = normalize_optional_colname(time_var),
        subgroup_color   = get_group_color(group_nm, meta$group_palette),
        subgroup_run_id  = make_subgroup_run_id(meta, group_nm)
      )
      
      if (is.null(support_tbl) || nrow(support_tbl) == 0) {
        return(
          tibble::as_tibble(c(
            base_cols,
            list(event_time = numeric(), n_ptids = integer(),
                 n_fireids = integer(), n_rows_treated = integer())
          ))
        )
      }
      
      support_tbl |>
        dplyr::mutate(!!!base_cols) |>
        dplyr::select(
          subset_id, outcome, group_id, group_col, subgroup,
          model_id, formula_template, estimator_type,
          cohort_var, time_var, event_time,
          n_ptids, n_fireids, n_rows_treated,
          subgroup_color, subgroup_run_id
        )
    }
  )
}


make_estimation_registry <- function(run_id, subset_id, outcome, group_id,
                                     model_id, group_col, dataset_spec,
                                     formula_template, estimator_type,
                                     term_pattern, data_filter, data_source,
                                     vcov_specs, group_args,
                                     n_rows_read, n_rows_after_filter,
                                     subgroup_tbl, support_tbl,
                                     fit_started, fit_finished,
                                     coef_file, subgroup_file, support_file,
                                     vcov_file, registry_file, run_spec_file) {
  
  model_formula <- build_model_formula(formula_template, outcome)
  
  tibble::tibble(
    run_id              = run_id,
    subset_id           = subset_id,
    outcome             = outcome,
    group_id            = group_id,
    model_id            = model_id,
    group_col           = group_col,
    unit_id             = dataset_spec$unit_id,
    event_id            = dataset_spec$event_id,
    trt_col             = dataset_spec$trt_col,
    cohort_var          = dataset_spec$cohort_var,
    time_var            = dataset_spec$time_var,
    formula_template    = formula_template,
    formula_resolved    = paste(deparse(model_formula), collapse = " "),
    estimator_type      = estimator_type,
    term_pattern        = term_pattern,
    data_filter         = data_filter_to_chr(data_filter),
    input_source        = paste(unlist(data_source), collapse = " | "),
    group_args_json     = as.character(serialize_object_json(group_args)),
    vcov_ids            = paste(vcov_specs$vcov_id, collapse = " | "),
    vcov_labels         = paste(vcov_specs$vcov_label, collapse = " | "),
    vcov_vars_json      = as.character(serialize_object_json(vcov_specs$vcov_vars)),
    n_rows_read         = n_rows_read,
    n_rows_after_filter = n_rows_after_filter,
    n_subgroups         = dplyr::n_distinct(subgroup_tbl$subgroup),
    subgroup_levels     = paste(sort(unique(subgroup_tbl$subgroup)), collapse = ","),
    n_support_rows      = nrow(support_tbl),
    fit_started         = fit_started,
    fit_finished        = fit_finished,
    coef_file           = coef_file,
    subgroup_file       = subgroup_file,
    support_file        = support_file,
    vcov_file           = vcov_file,
    registry_file       = registry_file,
    run_spec_file       = run_spec_file
  )
}


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 10 — Descriptive internals
# ══════════════════════════════════════════════════════════════════════════════

make_descriptive_trajectory <- function(df,
                                        meta,
                                        dataset_spec,
                                        treated_year_var,
                                        control_year_var) {
  
  trt_col   <- dataset_spec$trt_col
  unit_id   <- dataset_spec$unit_id
  event_id  <- dataset_spec$event_id
  time_var  <- dataset_spec$time_var
  group_col <- meta$group_col
  outcome   <- meta$outcome
  
  subgroup_levels <- df |>
    dplyr::filter(.data[[trt_col]] == 1, !is.na(.data[[group_col]])) |>
    dplyr::distinct(.data[[group_col]]) |>
    dplyr::pull(.data[[group_col]]) |>
    as.character() |>
    sort()
  
  purrr::map_dfr(
    subgroup_levels,
    \(grp) {
      tbl <- make_descriptive_one_subgroup(
        df               = df,
        group_col        = group_col,
        group_value      = grp,
        trt_col          = trt_col,
        outcome          = outcome,
        unit_id          = unit_id,
        event_id         = event_id,
        time_var         = time_var,
        treated_year_var = treated_year_var,
        control_year_var = control_year_var
      )
      tbl |>
        dplyr::mutate(
          subset_id        = meta$subset_id,
          outcome          = meta$outcome,
          group_id         = meta$group_id,
          group_col        = group_col,
          model_id         = meta$model_id,
          formula_template = meta$formula_template,
          treated_year_var = treated_year_var,
          control_year_var = control_year_var,
          time_var         = time_var,
          subgroup_color   = get_group_color(grp, meta$group_palette),
          subgroup_run_id  = make_subgroup_run_id(meta, grp),
          descriptive_run_id = glue::glue(
            "{make_subgroup_run_id(meta, grp)}__{tbl$series}"
          )
        ) |>
        dplyr::select(
          subset_id, outcome, group_id, group_col, subgroup,
          model_id, formula_template,
          treated_year_var, control_year_var, time_var,
          event_time, series,
          mean_outcome, sd_outcome, n_rows, n_ptids, n_fireids,
          subgroup_color, subgroup_run_id, descriptive_run_id
        )
    }
  )
}


make_descriptive_one_subgroup <- function(df,
                                          group_col,
                                          group_value,
                                          trt_col,
                                          outcome,
                                          unit_id,
                                          event_id,
                                          time_var,
                                          treated_year_var,
                                          control_year_var) {
  
  event_id <- normalize_optional_colname(event_id)
  
  needed_cols <- unique(stats::na.omit(c(
    group_col, trt_col, outcome, unit_id, event_id,
    time_var, treated_year_var, control_year_var
  )))
  check_required_columns(df, needed_cols, context = glue::glue("descriptive: {group_value}"))
  
  d_treated <- df |>
    dplyr::filter(.data[[trt_col]] == 1, .data[[group_col]] == group_value) |>
    dplyr::mutate(
      event_time = as.numeric(.data[[time_var]]) - as.numeric(.data[[treated_year_var]]),
      series     = "treated"
    ) |>
    dplyr::filter(!is.na(event_time), is.finite(event_time))
  
  d_control <- df |>
    dplyr::filter(.data[[trt_col]] == 0) |>
    dplyr::mutate(
      event_time = as.numeric(.data[[time_var]]) - as.numeric(.data[[control_year_var]]),
      series     = "control_mock"
    ) |>
    dplyr::filter(!is.na(event_time), is.finite(event_time))
  
  d_plot <- dplyr::bind_rows(d_treated, d_control)
  
  if (nrow(d_plot) == 0) {
    return(tibble::tibble(
      subgroup = character(), event_time = numeric(), series = character(),
      mean_outcome = numeric(), sd_outcome = numeric(),
      n_rows = integer(), n_ptids = integer(), n_fireids = integer()
    ))
  }
  
  grp_vars <- c("series", "event_time")
  
  if (is.na(event_id)) {
    d_plot |>
      dplyr::group_by(dplyr::across(dplyr::all_of(grp_vars))) |>
      dplyr::summarise(
        mean_outcome = mean(.data[[outcome]], na.rm = TRUE),
        sd_outcome   = stats::sd(.data[[outcome]], na.rm = TRUE),
        n_rows       = dplyr::n(),
        n_ptids      = dplyr::n_distinct(.data[[unit_id]]),
        n_fireids    = NA_integer_,
        .groups = "drop"
      ) |>
      dplyr::mutate(subgroup = group_value) |>
      dplyr::arrange(subgroup, series, event_time)
  } else {
    d_plot |>
      dplyr::group_by(dplyr::across(dplyr::all_of(grp_vars))) |>
      dplyr::summarise(
        mean_outcome = mean(.data[[outcome]], na.rm = TRUE),
        sd_outcome   = stats::sd(.data[[outcome]], na.rm = TRUE),
        n_rows       = dplyr::n(),
        n_ptids      = dplyr::n_distinct(.data[[unit_id]]),
        n_fireids    = dplyr::n_distinct(.data[[event_id]]),
        .groups = "drop"
      ) |>
      dplyr::mutate(subgroup = group_value) |>
      dplyr::arrange(subgroup, series, event_time)
  }
}


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 11 — Sun-Abraham support helpers
# ══════════════════════════════════════════════════════════════════════════════

parse_sunab_vars <- function(formula) {
  txt <- paste(deparse(formula), collapse = " ")
  txt <- gsub("\\s+", " ", txt)
  
  m <- stringr::str_match(
    txt,
    "sunab\\s*\\(\\s*([^,]+?)\\s*,\\s*([^,\\)]+?)\\s*(?:,|\\))"
  )
  
  if (length(m) == 0 || all(is.na(m[1, ]))) {
    warning(
      "Could not parse sunab() vars from formula: ", txt,
      ". event_time_support will be empty. ",
      "Check that the formula contains sunab(cohort_var, time_var, ...).",
      call. = FALSE
    )
    return(list(cohort_var = NA_character_, time_var = NA_character_))
  }
  
  list(
    cohort_var = trimws(m[1, 2]),
    time_var   = trimws(m[1, 3])
  )
}


make_event_time_support_one_model <- function(df_est,
                                              trt_col,
                                              unit_id,
                                              event_id   = NA_character_,
                                              cohort_var,
                                              time_var) {
  
  event_id   <- normalize_optional_colname(event_id)
  cohort_var <- normalize_optional_colname(cohort_var)
  time_var   <- normalize_optional_colname(time_var)
  
  if (is.na(cohort_var) || is.na(time_var)) {
    return(empty_event_time_support())
  }
  
  needed <- unique(stats::na.omit(c(trt_col, unit_id, event_id, cohort_var, time_var)))
  check_required_columns(df_est, needed, context = "event_time_support")
  
  d_treated <- df_est |>
    dplyr::filter(.data[[trt_col]] == 1) |>
    dplyr::mutate(
      event_time = suppressWarnings(
        as.numeric(.data[[time_var]]) - as.numeric(.data[[cohort_var]])
      )
    ) |>
    # Guard against never-treated sentinels (Inf, large values, NA)
    dplyr::filter(!is.na(event_time), is.finite(event_time))
  
  if (nrow(d_treated) == 0) {
    return(empty_event_time_support())
  }
  
  if (is.na(event_id)) {
    d_treated |>
      dplyr::group_by(event_time) |>
      dplyr::summarise(
        n_ptids        = dplyr::n_distinct(.data[[unit_id]]),
        n_fireids      = NA_integer_,
        n_rows_treated = dplyr::n(),
        .groups = "drop"
      ) |>
      dplyr::arrange(event_time)
  } else {
    d_treated |>
      dplyr::group_by(event_time) |>
      dplyr::summarise(
        n_ptids        = dplyr::n_distinct(.data[[unit_id]]),
        n_fireids      = dplyr::n_distinct(.data[[event_id]]),
        n_rows_treated = dplyr::n(),
        .groups = "drop"
      ) |>
      dplyr::arrange(event_time)
  }
}


empty_event_time_support <- function() {
  tibble::tibble(
    event_time     = numeric(),
    n_ptids        = integer(),
    n_fireids      = integer(),
    n_rows_treated = integer()
  )
}


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 12 — Arrow data loading helpers
# ══════════════════════════════════════════════════════════════════════════════

#' Open an Arrow Dataset from a directory or vector of file paths
open_arrow_source <- function(data_source) {
  src <- unlist(data_source)
  
  if (length(src) == 1 && dir.exists(src)) {
    arrow::open_dataset(src, format = "parquet")
  } else {
    arrow::open_dataset(src, format = "parquet")
  }
}


#' Load data from an Arrow source, pushing the filter before collect()
#'
#' @param data_source  Character vector of file paths or a single directory.
#' @param data_filter  One-sided formula or NULL.
load_arrow_data <- function(data_source, data_filter) {
  ds <- open_arrow_source(data_source)
  
  if (!is.null(data_filter)) {
    filter_quo <- rlang::new_quosure(
      rlang::f_rhs(data_filter),
      env = rlang::f_env(data_filter)
    )
    ds <- ds |> dplyr::filter(!!filter_quo)
  }
  
  ds |> dplyr::collect()
}


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 13 — Treatment grouping helper
# ══════════════════════════════════════════════════════════════════════════════

#' Apply a grouping function and validate its output
#'
#' Injects group_col into group_args (under the key "group_col") so that
#' group_fun can use it, then checks the returned data frame has the expected
#' column.
apply_treatment_grouping <- function(df, group_fun, group_col, group_args) {
  args <- c(list(df = df, group_col = group_col), group_args)
  df_out <- do.call(group_fun, args)
  
  if (!group_col %in% names(df_out)) {
    stop(glue::glue(
      "group_fun did not produce expected column '{group_col}'. ",
      "Check that your grouping function uses the group_col argument to name its output."
    ))
  }
  
  df_out
}


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 14 — Utility functions
# ══════════════════════════════════════════════════════════════════════════════

# ── Formula helpers ───────────────────────────────────────────────────────────

build_model_formula <- function(formula_template, outcome) {
  if (!is.character(formula_template) || length(formula_template) != 1) {
    stop("formula_template must be a length-1 character string.")
  }
  stats::as.formula(glue::glue(formula_template, outcome = outcome))
}


parse_formula_vars <- function(formula) {
  all.vars(formula)
}


# ── Column helpers ────────────────────────────────────────────────────────────

normalize_optional_colname <- function(x) {
  if (is.null(x) || length(x) == 0) return(NA_character_)
  if (length(x) > 1) {
    warning("normalize_optional_colname received length > 1 input; using first element only.")
    x <- x[1]
  }
  if (is.na(x)) return(NA_character_)
  as.character(x)
}


get_needed_columns <- function(formula,
                               trt_col,
                               group_col,
                               unit_id,
                               event_id  = NA_character_,
                               vcov_vars = NULL) {
  cols <- c(
    all.vars(formula),
    trt_col,
    group_col,
    unit_id,
    event_id,
    unlist(vcov_vars, use.names = FALSE)
  )
  cols |> unique() |> stats::na.omit() |> as.character()
}


check_required_columns <- function(df, required, context = "") {
  missing <- setdiff(required, names(df))
  if (length(missing) > 0) {
    stop(glue::glue(
      "{if (nchar(context) > 0) paste0('[', context, '] ') else ''}",
      "Missing required columns: ",
      paste(missing, collapse = ", ")
    ))
  }
  invisible(TRUE)
}


compute_n_distinct_optional <- function(df, col_nm) {
  col_nm <- normalize_optional_colname(col_nm)
  if (is.na(col_nm)) return(NA_integer_)
  if (!col_nm %in% names(df)) {
    stop(glue::glue("Column '{col_nm}' not found in data."))
  }
  dplyr::n_distinct(stats::na.omit(df[[col_nm]]))
}


# ── Path / ID helpers ─────────────────────────────────────────────────────────

safe_path_component <- function(x) {
  gsub("[^[:alnum:]_\\-]+", "_", as.character(x))
}


make_run_stub <- function(subset_id, outcome, group_id, model_id) {
  paste(
    safe_path_component(subset_id),
    safe_path_component(outcome),
    safe_path_component(group_id),
    safe_path_component(model_id),
    sep = "__"
  )
}


make_subgroup_run_id <- function(meta, group_value) {
  glue::glue(
    "{meta$subset_id}__{meta$outcome}__{meta$group_id}__{meta$model_id}__{group_value}"
  )
}


# ── Misc helpers ──────────────────────────────────────────────────────────────

extract_first_number <- function(x) {
  suppressWarnings(as.numeric(stringr::str_extract(x, "-?\\d+")))
}


dir_ensure_local <- function(paths) {
  purrr::walk(paths, \(p) {
    dir.create(p, recursive = TRUE, showWarnings = FALSE)
    if (!dir.exists(p)) stop("Failed to create directory: ", p)
  })
}


timestamp_now <- function() {
  format(Sys.time(), "%Y%m%d_%H%M%S")
}


serialize_object_json <- function(x) {
  # jsonlite::toJSON returns class "json", not plain character; as.character() strips it
  as.character(
    jsonlite::toJSON(x, auto_unbox = TRUE, null = "null", digits = NA, pretty = FALSE)
  )
}


data_filter_to_chr <- function(data_filter = NULL) {
  if (is.null(data_filter))               return(NA_character_)
  if (inherits(data_filter, "formula"))   return(paste(deparse(data_filter), collapse = " "))
  as.character(data_filter)
}


get_group_color <- function(group_value, group_palette) {
  if (!is.null(group_palette) && group_value %in% names(group_palette)) {
    return(unname(group_palette[[group_value]]))
  }
  NA_character_
}


# ── Spec validators ───────────────────────────────────────────────────────────

validate_dataset_spec <- function(spec) {
  required <- c("unit_id", "time_var", "trt_col", "cohort_var", "event_id")
  missing  <- setdiff(required, names(spec))
  if (length(missing) > 0) {
    stop("dataset_spec is missing fields: ", paste(missing, collapse = ", "))
  }
  invisible(TRUE)
}


validate_spec_table <- function(tbl, required_cols, name) {
  if (!is.data.frame(tbl)) {
    stop(name, " must be a data frame / tibble.")
  }
  missing <- setdiff(required_cols, names(tbl))
  if (length(missing) > 0) {
    stop(name, " is missing columns: ", paste(missing, collapse = ", "))
  }
  invisible(TRUE)
}


# ── Filter construction helpers ───────────────────────────────────────────────

#' Build a one-sided equality filter formula for a single column value.
#' Intended for use inside expand_analysis_subset_specs_by_col.
#' @param col   Column name (string).
#' @param value Value to match (will be compared with ==).
#' @return A one-sided formula, e.g. ~ forest_type == "pine"
build_equality_filter <- function(col, value) {
  expr <- call("==", as.name(col), value)
  rlang::new_formula(lhs = NULL, rhs = expr, env = baseenv())
}


#' Combine two one-sided filter formulas with &.
#' Either argument may be NULL (in which case the other is returned as-is).
combine_filters <- function(filter_a, filter_b) {
  if (is.null(filter_a)) return(filter_b)
  if (is.null(filter_b)) return(filter_a)
  
  expr <- call("&", rlang::f_rhs(filter_a), rlang::f_rhs(filter_b))
  rlang::new_formula(lhs = NULL, rhs = expr, env = rlang::f_env(filter_a))
}


# ── Null-coalescing operator (base-R safe) ────────────────────────────────────

`%||%` <- function(a, b) if (!is.null(a)) a else b