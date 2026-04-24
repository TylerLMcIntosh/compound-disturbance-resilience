# portable_sunab.R
# Portable Sun-Abraham / TWFE DiD estimation pipeline
# -------------------------------------------------------
# Design principles:
#   - All dataset-specific column names live in a single dataset_spec object
#   - Analysis subsets are defined by tidy filter expressions pushed to Arrow
#   - Each subset carries both a long data source (estimation) and an optional
#     short data source (one row per unit, for weighting)
#   - treatment_group_specs define how treated units are partitioned
#   - Weighting is a separate pre-computation stage (subset x group x weighting
#     spec); outputs are weights parquets keyed on unit_id, joined into the
#     long panel at estimation time
#   - Per-run folders are the source of truth; rebuild merges them
#   - One shared dir_out per project; multiple experiment calls write into it
# -------------------------------------------------------


# ── Required packages ────────────────────────────────────────────────────────

.required_pkgs <- c(
  "arrow", "dplyr", "fixest", "glue", "purrr",
  "readr", "rlang", "stringr", "tibble", "tidyr", "jsonlite"
)
.missing_pkgs <- .required_pkgs[
  !vapply(.required_pkgs, requireNamespace, logical(1), quietly = TRUE)
]
if (length(.missing_pkgs) > 0) {
  stop("Missing required packages: ", paste(.missing_pkgs, collapse = ", "), call. = FALSE)
}
rm(.required_pkgs, .missing_pkgs)


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 1 — Dataset spec
# ══════════════════════════════════════════════════════════════════════════════

#' Build a dataset spec
#'
#' Captures all column-name mappings that are properties of the DATA, not of
#' any particular model or grouping scheme. One spec per experiment call.
#'
#' @param unit_id    Unit identifier column (e.g. "pt_id").
#' @param time_var   Calendar-time column (e.g. "year").
#' @param trt_col    Binary treatment indicator column (0/1).
#' @param cohort_var Treatment-cohort column for sunab() (e.g. "FirstTreat").
#'                   NA_character_ for plain TWFE.
#' @param event_id   Event identifier column (e.g. "fireid").
#'                   NA_character_ if not applicable.
make_dataset_spec <- function(unit_id,
                              time_var,
                              trt_col,
                              cohort_var = NA_character_,
                              event_id   = NA_character_) {
  stopifnot(
    is.character(unit_id)  && length(unit_id)  == 1,
    is.character(time_var) && length(time_var) == 1,
    is.character(trt_col)  && length(trt_col)  == 1
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

#' Build one analysis subset spec row
#'
#' Each row defines one data slice. data_filter is pushed to Arrow before
#' collect() for both long and short sources.
#'
#' @param subset_id          Short human-readable identifier.
#' @param long_data_source   Directory path or character vector of parquet
#'                           file paths for the long (panel) data used in
#'                           estimation and descriptive pipelines.
#' @param data_filter        One-sided formula (e.g. ~ year >= 1997) or NULL.
#'                           Applied to BOTH long and short sources.
#' @param short_data_source  Optional. Directory or file vector for the
#'                           cross-sectional (one row per unit) data used in
#'                           the weighting pipeline. NULL when no weighting
#'                           is planned for this subset.
make_analysis_subset_spec <- function(subset_id,
                                      long_data_source,
                                      data_filter       = NULL,
                                      short_data_source = NULL) {
  
  if (!is.character(subset_id) || length(subset_id) != 1 || is.na(subset_id)) {
    stop("subset_id must be a length-1 non-NA character string.")
  }
  validate_data_source(long_data_source,
                       label = glue::glue("subset '{subset_id}' long_data_source"))
  if (!is.null(short_data_source)) {
    validate_data_source(short_data_source,
                         label = glue::glue("subset '{subset_id}' short_data_source"))
  }
  if (!is.null(data_filter) && !inherits(data_filter, "formula")) {
    stop("data_filter must be NULL or a one-sided formula, e.g. ~ year >= 1997")
  }
  
  tibble::tibble(
    subset_id         = subset_id,
    long_data_source  = list(long_data_source),
    short_data_source = list(short_data_source),
    data_filter       = list(data_filter)
  )
}


#' Auto-expand analysis subset specs by unique values of a column
#'
#' Reads unique values of split_col from long_data_source via Arrow (no full
#' collect), builds one spec row per value optionally combined with base_filter.
#'
#' @param long_data_source   Directory or file vector for the long data.
#' @param split_col          Column to split on.
#' @param id_prefix          Prefix for auto-generated subset_id values.
#' @param base_filter        Optional one-sided formula applied on top of the
#'                           per-value filter.
#' @param values             If not NULL, restrict to this subset of values.
#' @param short_data_source  Passed through to make_analysis_subset_spec.
#'                           NULL if no weighting is planned.
expand_analysis_subset_specs_by_col <- function(long_data_source,
                                                split_col,
                                                id_prefix         = split_col,
                                                base_filter       = NULL,
                                                values            = NULL,
                                                short_data_source = NULL) {
  
  ds <- open_arrow_source(long_data_source)
  
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
      stop("No overlap between provided values and unique values of '",
           split_col, "' in data.")
    }
  }
  
  purrr::map_dfr(unique_vals, function(val) {
    combined_filter <- combine_filters(base_filter, build_equality_filter(split_col, val))
    make_analysis_subset_spec(
      subset_id         = as.character(glue::glue("{id_prefix}_{safe_path_component(val)}")),
      long_data_source  = long_data_source,
      data_filter       = combined_filter,
      short_data_source = short_data_source
    )
  })
}


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 3 — Treatment group specs
# ══════════════════════════════════════════════════════════════════════════════

#' Build one treatment group spec row
#'
#' group_fun contract:
#'   Receives: df, group_col (injected by pipeline), plus everything in
#'   group_args. Must return df with a new column named exactly group_col.
#'   When include_control = TRUE (injected during weighting), must also assign
#'   the value "control" to rows where trt_col == 0.
#'
#' @param group_id    Short identifier.
#' @param group_col   Name of the subgroup column group_fun will create.
#' @param group_fun   Function with signature f(df, group_col, ...).
#' @param group_args  Named list forwarded to group_fun. Must not contain
#'                    "group_col".
make_treatment_group_spec <- function(group_id,
                                      group_col,
                                      group_fun,
                                      group_args = list()) {
  
  if (!is.character(group_id)  || length(group_id)  != 1) stop("group_id must be length-1 character.")
  if (!is.character(group_col) || length(group_col) != 1) stop("group_col must be length-1 character.")
  if (!is.function(group_fun))                            stop("group_fun must be a function.")
  if ("group_col" %in% names(group_args)) {
    stop("group_args must not contain 'group_col'; set it via the top-level argument.")
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

#' Build one model spec row
#'
#' @param model_id         Short identifier.
#' @param formula_template Glue string with {outcome} placeholder.
#' @param estimator_type   "sunab" or "twfe". Controls post-processing.
#' @param term_pattern     Regex flagging event-time coefficient names.
#' @param weights_col      Name of the weight column to pass to feols()
#'                         (e.g. "glm_ato_weights"). NA_character_ = unweighted.
make_model_spec <- function(model_id,
                            formula_template,
                            estimator_type = c("sunab", "twfe"),
                            term_pattern   = ".*",
                            weights_col    = NA_character_) {
  
  estimator_type <- match.arg(estimator_type)
  if (!is.character(model_id)         || length(model_id)         != 1) stop("model_id must be length-1 character.")
  if (!is.character(formula_template) || length(formula_template) != 1) stop("formula_template must be length-1 character.")
  
  tibble::tibble(
    model_id         = model_id,
    formula_template = formula_template,
    estimator_type   = estimator_type,
    term_pattern     = term_pattern,
    weights_col      = normalize_optional_colname(weights_col)
  )
}


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 5 — Weighting specs
# ══════════════════════════════════════════════════════════════════════════════

#' Build one weighting spec row
#'
#' Defines how propensity-score weights are estimated for one weighting
#' approach. The weighting pipeline crosses analysis_subset_specs x
#' treatment_group_specs x weighting_specs.
#'
#' The output weight column name is:
#'   "{weighting_id}_{weighting_name}_weights"  when weighting_name is not NA
#'   "{weighting_id}_weights"                    otherwise
#'
#' @param weighting_id      Short identifier (e.g. "glm_ato").
#' @param weight_formula    One-sided formula of covariates, e.g.
#'                          ~ chili + def + aet + srtm + tpi.
#'                          The LHS (treatment variable = group_col) is
#'                          constructed automatically.
#' @param method            WeightIt method string (e.g. "glm", "gbm").
#' @param estimand          WeightIt estimand string (e.g. "ATO", "ATE").
#' @param weighting_name    Optional suffix for the weight column name.
#'                          NA_character_ = no suffix.
make_weighting_spec <- function(weighting_id,
                                weight_formula,
                                method         = "glm",
                                estimand       = "ATO",
                                weighting_name = NA_character_) {
  
  if (!is.character(weighting_id) || length(weighting_id) != 1) {
    stop("weighting_id must be a length-1 character string.")
  }
  if (!inherits(weight_formula, "formula")) {
    stop("weight_formula must be a formula, e.g. ~ chili + def + aet")
  }
  if (length(attr(stats::terms(weight_formula), "term.labels")) == 0) {
    stop("weight_formula must contain at least one covariate on the RHS.")
  }
  
  tibble::tibble(
    weighting_id   = weighting_id,
    weight_formula = list(weight_formula),
    method         = method,
    estimand       = estimand,
    weighting_name = normalize_optional_colname(weighting_name)
  )
}


#' Construct the weight column name for a weighting_id / weighting_name pair
make_weight_col_name <- function(weighting_id, weighting_name) {
  if (!is.na(weighting_name) && nchar(weighting_name) > 0) {
    paste(weighting_id, weighting_name, "weights", sep = "_")
  } else {
    paste(weighting_id, "weights", sep = "_")
  }
}


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 6 — Weighting experiment runner
# ══════════════════════════════════════════════════════════════════════════════

#' Run the weighting pre-computation stage
#'
#' Crosses analysis_subset_specs x treatment_group_specs x weighting_specs.
#' For each combination:
#'   1. Loads the short (cross-sectional) data via Arrow with data_filter.
#'   2. Applies treatment grouping with include_control = TRUE so untreated
#'      units receive the value "control" in group_col.
#'   3. Fits weights with WeightIt using group_col as the treatment variable.
#'   4. Saves a weights parquet (unit_id + weight_col), a cobalt balance RDS,
#'      a love plot PNG, and a registry parquet.
#'
#' @param dataset_spec          Output of make_dataset_spec().
#' @param analysis_subset_specs Subset specs. Each row must have a non-NULL
#'                              short_data_source.
#' @param treatment_group_specs Treatment group specs.
#' @param weighting_specs       Weighting specs.
#' @param dir_out               Root output directory (same as run_experiment).
#' @param skip_existing         Skip combinations whose outputs already exist.
#' @param verbose_timing        Print timing per combination.
#' @param .progress             Progress bar for purrr::pmap.
run_weighting_experiment <- function(dataset_spec,
                                     analysis_subset_specs,
                                     treatment_group_specs,
                                     weighting_specs,
                                     dir_out,
                                     skip_existing  = TRUE,
                                     verbose_timing = FALSE,
                                     .progress      = FALSE) {
  
  for (pkg in c("WeightIt", "cobalt", "ggplot2")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("Package '", pkg, "' is required for the weighting pipeline. ",
           "Install with install.packages('", pkg, "').")
    }
  }
  
  validate_dataset_spec(dataset_spec)
  validate_spec_table(analysis_subset_specs, c("subset_id", "long_data_source", "short_data_source", "data_filter"), "analysis_subset_specs")
  validate_spec_table(treatment_group_specs,  c("group_id", "group_col", "group_fun", "group_args"),                  "treatment_group_specs")
  validate_spec_table(weighting_specs,         c("weighting_id", "weight_formula", "method", "estimand"),              "weighting_specs")
  
  missing_short <- analysis_subset_specs |>
    dplyr::filter(purrr::map_lgl(short_data_source, is.null)) |>
    dplyr::pull(subset_id)
  if (length(missing_short) > 0) {
    stop(
      "These subsets have no short_data_source and cannot be weighted:\n",
      paste(missing_short, collapse = "\n"),
      "\nProvide short_data_source in make_analysis_subset_spec()."
    )
  }
  
  dir_ensure_local(c(
    dir_out,
    file.path(dir_out, "weights", "by_run"),
    file.path(dir_out, "metadata")
  ))
  
  run_grid <- tidyr::crossing(
    analysis_subset_specs,
    treatment_group_specs,
    weighting_specs
  ) |>
    dplyr::mutate(
      weight_run_id = as.character(glue::glue("{subset_id}__{group_id}__{weighting_id}")),
      weight_col    = purrr::map2_chr(weighting_id, weighting_name, make_weight_col_name)
    )
  
  ts <- timestamp_now()
  saveRDS(
    list(
      dataset_spec          = dataset_spec,
      analysis_subset_specs = analysis_subset_specs,
      treatment_group_specs = treatment_group_specs,
      weighting_specs       = weighting_specs,
      skip_existing         = skip_existing
    ),
    file.path(dir_out, "metadata",
              glue::glue("weighting_specs_snapshot__{ts}.rds"))
  )
  
  if (verbose_timing) {
    message(glue::glue(
      "[{timestamp_now()}] Starting weighting: {nrow(run_grid)} combinations."
    ))
  }
  
  results <- purrr::pmap(
    list(
      short_data_source = run_grid$short_data_source,
      data_filter       = run_grid$data_filter,
      subset_id         = run_grid$subset_id,
      group_id          = run_grid$group_id,
      group_col         = run_grid$group_col,
      group_fun         = run_grid$group_fun,
      group_args        = run_grid$group_args,
      weighting_id      = run_grid$weighting_id,
      weight_formula    = run_grid$weight_formula,
      method            = run_grid$method,
      estimand          = run_grid$estimand,
      weight_col        = run_grid$weight_col,
      weight_run_id     = run_grid$weight_run_id
    ),
    function(short_data_source, data_filter, subset_id,
             group_id, group_col, group_fun, group_args,
             weighting_id, weight_formula, method, estimand,
             weight_col, weight_run_id) {
      
      result <- list(weight_run_id = weight_run_id, error = NULL, skipped = FALSE)
      t_start <- proc.time()
      
      tryCatch({
        result <- run_one_weighting(
          short_data_source = short_data_source,
          data_filter       = data_filter,
          dataset_spec      = dataset_spec,
          subset_id         = subset_id,
          group_id          = group_id,
          group_col         = group_col,
          group_fun         = group_fun,
          group_args        = group_args,
          weighting_id      = weighting_id,
          weight_formula    = weight_formula,
          method            = method,
          estimand          = estimand,
          weight_col        = weight_col,
          weight_run_id     = weight_run_id,
          dir_out           = dir_out,
          skip_existing     = skip_existing
        )
      }, error = function(e) {
        result$error <<- conditionMessage(e)
        message(glue::glue("[ERROR] weight_run_id = {weight_run_id}: {conditionMessage(e)}"))
      })
      
      if (verbose_timing) {
        elapsed <- (proc.time() - t_start)[["elapsed"]]
        status  <- if (!is.null(result$error)) "FAILED" else if (isTRUE(result$skipped)) "SKIPPED" else "OK"
        message(glue::glue("[{timestamp_now()}] {weight_run_id} | {status} | {round(elapsed, 1)}s"))
      }
      
      result
    },
    .progress = .progress
  )
  
  failed <- purrr::keep(results, \(r) !is.null(r$error))
  if (length(failed) > 0) {
    message(glue::glue(
      "\n{length(failed)} weighting run(s) failed:\n",
      paste(purrr::map_chr(failed, "weight_run_id"), collapse = "\n")
    ))
  }
  
  invisible(list(run_grid = run_grid, run_results = results))
}


# ── Single weighting run ──────────────────────────────────────────────────────

run_one_weighting <- function(short_data_source,
                              data_filter,
                              dataset_spec,
                              subset_id,
                              group_id,
                              group_col,
                              group_fun,
                              group_args,
                              weighting_id,
                              weight_formula,
                              method,
                              estimand,
                              weight_col,
                              weight_run_id,
                              dir_out,
                              skip_existing = TRUE) {
  
  run_stub <- as.character(glue::glue(
    "{safe_path_component(subset_id)}__{safe_path_component(group_id)}__{safe_path_component(weighting_id)}"
  ))
  run_dir <- file.path(dir_out, "weights", "by_run", run_stub)
  dir_ensure_local(run_dir)
  
  weights_file   <- file.path(run_dir, "weights.parquet")
  registry_file  <- file.path(run_dir, "registry.parquet")
  bal_rds_file   <- file.path(run_dir, "balance_objects.rds")
  love_plot_file <- file.path(run_dir, "love_plot.png")
  run_spec_file  <- file.path(run_dir, "weight_run_spec.rds")
  
  if (skip_existing && all(file.exists(c(weights_file, registry_file,
                                         bal_rds_file, love_plot_file,
                                         run_spec_file)))) {
    return(list(
      weight_run_id  = weight_run_id,
      weights_file   = weights_file,
      registry_file  = registry_file,
      bal_rds_file   = bal_rds_file,
      love_plot_file = love_plot_file,
      run_spec_file  = run_spec_file,
      skipped        = TRUE,
      error          = NULL
    ))
  }
  
  run_started <- Sys.time()
  
  # ── Load short data ────────────────────────────────────────────────────────
  df_short <- load_arrow_data(short_data_source, data_filter)
  
  # ── Apply grouping with include_control = TRUE ─────────────────────────────
  # include_control is injected here; group_fun must accept and use it
  df_grouped <- apply_treatment_grouping(
    df         = df_short,
    group_fun  = group_fun,
    group_col  = group_col,
    group_args = c(group_args, list(include_control = TRUE))
  )
  
  df_grouped <- df_grouped |> dplyr::filter(!is.na(.data[[group_col]]))
  
  if (nrow(df_grouped) == 0) {
    stop("No rows remain after grouping with include_control = TRUE for ",
         weight_run_id)
  }
  
  # ── Build full WeightIt formula (LHS = group_col, RHS from weight_formula) ─
  rhs_chr   <- sub("^~\\s*", "", paste(deparse(weight_formula), collapse = " "))
  w_formula <- stats::as.formula(paste(group_col, "~", rhs_chr))
  
  # ── Fit weights ────────────────────────────────────────────────────────────
  w_out <- WeightIt::weightit(
    formula  = w_formula,
    data     = df_grouped,
    method   = method,
    estimand = estimand
  )
  
  # ── Weights parquet ────────────────────────────────────────────────────────
  unit_id     <- dataset_spec$unit_id
  weights_tbl <- tibble::tibble(
    !!unit_id      := df_grouped[[unit_id]],
    !!weight_col   := w_out$weights,
    subset_id      = subset_id,
    group_id       = group_id,
    weighting_id   = weighting_id,
    weight_run_id  = weight_run_id
  )
  
  # ── Balance diagnostics ────────────────────────────────────────────────────
  bal_obj <- cobalt::bal.tab(
    w_out,
    un         = TRUE,
    pairwise   = TRUE,
    abs        = TRUE,
    stats      = c("mean.diffs", "variance.ratios"),
    thresholds = c(m = 0.1, v = 2)
  )
  
  love_p <- cobalt::love.plot(
    w_out,
    stats      = "mean.diffs",
    abs        = TRUE,
    pairwise   = TRUE,
    thresholds = c(m = 0.1),
    var.order  = "unadjusted"
  ) +
    ggplot2::labs(
      title = as.character(glue::glue(
        "{subset_id} | {group_id} | {weighting_id} ({method}, {estimand})"
      ))
    )
  
  ggplot2::ggsave(
    love_plot_file,
    plot  = love_p,
    width = 2500, height = 2000, units = "px", bg = "white"
  )
  
  # ── Registry ───────────────────────────────────────────────────────────────
  run_finished <- Sys.time()
  
  registry_tbl <- tibble::tibble(
    weight_run_id     = weight_run_id,
    subset_id         = subset_id,
    group_id          = group_id,
    weighting_id      = weighting_id,
    unit_id           = unit_id,
    group_col         = group_col,
    weight_col        = weight_col,
    method            = method,
    estimand          = estimand,
    weight_formula    = paste(deparse(w_formula), collapse = " "),
    short_data_source = paste(unlist(short_data_source), collapse = " | "),
    data_filter       = data_filter_to_chr(data_filter),
    n_units           = nrow(df_grouped),
    n_groups          = dplyr::n_distinct(df_grouped[[group_col]]),
    group_levels      = paste(sort(unique(df_grouped[[group_col]])), collapse = ","),
    run_started       = run_started,
    run_finished      = run_finished,
    weights_file      = weights_file,
    registry_file     = registry_file,
    bal_rds_file      = bal_rds_file,
    love_plot_file    = love_plot_file,
    run_spec_file     = run_spec_file
  )
  
  run_spec <- list(
    weight_run_id     = weight_run_id,
    run_stub          = run_stub,
    subset_id         = subset_id,
    group_id          = group_id,
    weighting_id      = weighting_id,
    group_col         = group_col,
    weight_col        = weight_col,
    dataset_spec      = dataset_spec,
    weight_formula    = w_formula,
    method            = method,
    estimand          = estimand,
    short_data_source = short_data_source,
    data_filter       = data_filter,
    group_fun         = group_fun,
    group_args        = group_args
  )
  
  arrow::write_parquet(weights_tbl,  weights_file)
  arrow::write_parquet(registry_tbl, registry_file)
  saveRDS(bal_obj,  bal_rds_file)
  saveRDS(run_spec, run_spec_file)
  
  rm(df_short, df_grouped, weights_tbl, registry_tbl, run_spec, w_out, bal_obj)
  gc()
  
  list(
    weight_run_id  = weight_run_id,
    weights_file   = weights_file,
    registry_file  = registry_file,
    bal_rds_file   = bal_rds_file,
    love_plot_file = love_plot_file,
    run_spec_file  = run_spec_file,
    skipped        = FALSE,
    error          = NULL
  )
}


#' Rebuild weighting registry by merging all by_run outputs
rebuild_weighting_tables <- function(dir_out, write_csv = TRUE) {
  
  dir_by_run <- file.path(dir_out, "weights", "by_run")
  dir_all    <- file.path(dir_out, "weights", "all")
  
  if (!dir.exists(dir_by_run)) stop("Weights by_run directory does not exist: ", dir_by_run)
  dir.create(dir_all, recursive = TRUE, showWarnings = FALSE)
  
  registry_files <- list.files(dir_by_run, pattern = "^registry\\.parquet$",
                               recursive = TRUE, full.names = TRUE)
  if (length(registry_files) == 0) {
    stop("No weight registry.parquet files found under: ", dir_by_run)
  }
  
  registry_tbl <- purrr::map_dfr(
    registry_files,
    \(f) arrow::read_parquet(f) |> dplyr::mutate(.mtime = file.mtime(f))
  ) |>
    dplyr::group_by(weight_run_id) |>
    dplyr::filter(.mtime == max(.mtime)) |>
    dplyr::ungroup() |>
    dplyr::select(-.mtime)
  
  pq <- file.path(dir_all, "weight_registry.parquet")
  arrow::write_parquet(registry_tbl, pq)
  if (write_csv) readr::write_csv(registry_tbl, file.path(dir_all, "weight_registry.csv"))
  
  invisible(list(weight_registry = registry_tbl, files = list(parquet = pq)))
}


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 7 — Weights path resolver
# ══════════════════════════════════════════════════════════════════════════════

#' Locate the weights parquet for a given subset x group x weight column
#'
#' Returns NULL when weights_col is NA (unweighted model).
#' Errors informatively if a match is expected but not found.
resolve_weights_parquet <- function(weights_col, subset_id, group_id, dir_out) {
  
  if (is.na(weights_col)) return(NULL)
  
  dir_weights_run <- file.path(dir_out, "weights", "by_run")
  
  if (!dir.exists(dir_weights_run)) {
    stop(glue::glue(
      "Model requires weights_col = '{weights_col}' but no weights/by_run ",
      "directory exists under '{dir_out}'. Run run_weighting_experiment() first."
    ))
  }
  
  registry_files <- list.files(dir_weights_run, pattern = "^registry\\.parquet$",
                               recursive = TRUE, full.names = TRUE)
  if (length(registry_files) == 0) {
    stop("No weight registry files found. Run run_weighting_experiment() first.")
  }
  
  registry_tbl <- purrr::map_dfr(registry_files, arrow::read_parquet)
  
  match_row <- registry_tbl |>
    dplyr::filter(
      .data$subset_id  == !!subset_id,
      .data$group_id   == !!group_id,
      .data$weight_col == !!weights_col
    )
  
  if (nrow(match_row) == 0) {
    stop(glue::glue(
      "No weights found for subset_id='{subset_id}', group_id='{group_id}', ",
      "weight_col='{weights_col}'. ",
      "Check that run_weighting_experiment() completed for this combination."
    ))
  }
  
  # Keep most recently written if duplicated
  match_row <- match_row |>
    dplyr::mutate(.mtime = purrr::map_dbl(weights_file, \(f) as.numeric(file.mtime(f)))) |>
    dplyr::slice_max(.mtime, n = 1, with_ties = FALSE)
  
  match_row$weights_file[[1]]
}


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 8 — Top-level estimation experiment runner
# ══════════════════════════════════════════════════════════════════════════════

#' Run a full DiD experiment
#'
#' Crosses analysis_subset_specs x outcome_specs x treatment_group_specs x
#' model_specs. For weighted models (weights_col not NA in model_specs),
#' locates the matching weights parquet from weights/by_run and joins it into
#' the panel data before fitting.
#'
#' @param dataset_spec            Output of make_dataset_spec().
#' @param analysis_subset_specs   Tibble from make_analysis_subset_spec rows.
#' @param outcome_specs           Tibble with column "outcome".
#' @param treatment_group_specs   Tibble from make_treatment_group_spec rows.
#' @param model_specs             Tibble from make_model_spec rows.
#' @param vcov_specs              Tibble with vcov_id, vcov, vcov_label,
#'                                optionally vcov_vars.
#' @param dir_out                 Root output directory.
#' @param group_palette           Optional named colour vector for subgroups.
#' @param ci_level                CI level. Default 0.95.
#' @param run_estimation          Run estimation pipeline? Default TRUE.
#' @param run_descriptive         Run descriptive pipeline? Default FALSE.
#' @param descriptive_args        List with treated_year_var, control_year_var.
#' @param skip_existing           Skip runs whose outputs already exist.
#' @param write_vcov              Write full vcov matrices? Default FALSE.
#' @param verbose_timing          Print per-run timing? Default FALSE.
#' @param .progress               Progress bar for purrr::pmap.
run_experiment <- function(dataset_spec,
                           analysis_subset_specs,
                           outcome_specs,
                           treatment_group_specs,
                           model_specs,
                           vcov_specs,
                           dir_out,
                           group_palette    = NULL,
                           ci_level         = 0.95,
                           run_estimation   = TRUE,
                           run_descriptive  = FALSE,
                           descriptive_args = list(),
                           skip_existing    = TRUE,
                           write_vcov       = FALSE,
                           verbose_timing   = FALSE,
                           .progress        = FALSE) {
  
  validate_dataset_spec(dataset_spec)
  validate_spec_table(analysis_subset_specs, c("subset_id", "long_data_source", "data_filter"),         "analysis_subset_specs")
  validate_spec_table(outcome_specs,          c("outcome"),                                               "outcome_specs")
  validate_spec_table(treatment_group_specs,  c("group_id", "group_col", "group_fun", "group_args"),     "treatment_group_specs")
  validate_spec_table(model_specs,            c("model_id", "formula_template"),                          "model_specs")
  validate_spec_table(vcov_specs,             c("vcov_id", "vcov", "vcov_label"),                         "vcov_specs")
  
  if (!"vcov_vars"    %in% names(vcov_specs))  vcov_specs$vcov_vars  <- replicate(nrow(vcov_specs), character(0), simplify = FALSE)
  if (!"weights_col"  %in% names(model_specs)) model_specs$weights_col <- NA_character_
  
  if (!run_estimation && !run_descriptive) {
    stop("At least one of run_estimation or run_descriptive must be TRUE.")
  }
  
  desc_args <- list(
    treated_year_var = dataset_spec$time_var,
    control_year_var = dataset_spec$time_var
  )
  desc_args[names(descriptive_args)] <- descriptive_args
  
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
      run_id   = as.character(glue::glue("{subset_id}__{outcome}__{group_id}__{model_id}")),
      run_stub = purrr::pmap_chr(
        list(subset_id, outcome, group_id, model_id),
        \(s, o, g, m) make_run_stub(s, o, g, m)
      )
    )
  
  ts            <- timestamp_now()
  run_grid_file <- file.path(dir_out, "metadata", glue::glue("run_grid__{ts}.rds"))
  specs_file    <- file.path(dir_out, "metadata", glue::glue("specs_snapshot__{ts}.rds"))
  info_file     <- file.path(dir_out, "metadata", glue::glue("experiment_info__{ts}.rds"))
  
  saveRDS(run_grid, run_grid_file)
  saveRDS(
    list(
      dataset_spec          = dataset_spec,
      analysis_subset_specs = analysis_subset_specs,
      outcome_specs         = outcome_specs,
      treatment_group_specs = treatment_group_specs,
      model_specs           = model_specs,
      vcov_specs            = vcov_specs,
      group_palette         = group_palette,
      ci_level              = ci_level,
      run_estimation        = run_estimation,
      run_descriptive       = run_descriptive,
      descriptive_args      = desc_args,
      skip_existing         = skip_existing,
      write_vcov            = write_vcov
    ),
    specs_file
  )
  
  if (verbose_timing) {
    message(glue::glue("[{timestamp_now()}] Starting experiment: {nrow(run_grid)} runs planned."))
  }
  
  run_results <- purrr::pmap(
    list(
      long_data_source = run_grid$long_data_source,
      data_filter      = run_grid$data_filter,
      subset_id        = run_grid$subset_id,
      outcome          = run_grid$outcome,
      group_id         = run_grid$group_id,
      group_col        = run_grid$group_col,
      group_fun        = run_grid$group_fun,
      group_args       = run_grid$group_args,
      model_id         = run_grid$model_id,
      formula_template = run_grid$formula_template,
      estimator_type   = run_grid$estimator_type,
      term_pattern     = run_grid$term_pattern,
      weights_col      = run_grid$weights_col,
      run_id           = run_grid$run_id,
      run_stub         = run_grid$run_stub
    ),
    function(long_data_source, data_filter, subset_id, outcome,
             group_id, group_col, group_fun, group_args,
             model_id, formula_template, estimator_type, term_pattern,
             weights_col, run_id, run_stub) {
      
      result <- list(
        run_id = run_id, estimation = NULL, descriptive = NULL,
        error = NULL, skipped = FALSE
      )
      t_start <- proc.time()
      
      tryCatch({
        
        weights_parquet_path <- resolve_weights_parquet(
          weights_col = weights_col,
          subset_id   = subset_id,
          group_id    = group_id,
          dir_out     = dir_out
        )
        
        if (run_estimation) {
          result$estimation <- run_one_estimation(
            data_source          = long_data_source,
            data_filter          = data_filter,
            dataset_spec         = dataset_spec,
            subset_id            = subset_id,
            outcome              = outcome,
            group_id             = group_id,
            group_col            = group_col,
            group_fun            = group_fun,
            group_args           = group_args,
            model_id             = model_id,
            formula_template     = formula_template,
            estimator_type       = estimator_type,
            term_pattern         = term_pattern,
            weights_col          = weights_col,
            weights_parquet_path = weights_parquet_path,
            vcov_specs           = vcov_specs,
            dir_out              = dir_out,
            run_id               = run_id,
            run_stub             = run_stub,
            group_palette        = group_palette,
            ci_level             = ci_level,
            skip_existing        = skip_existing,
            write_vcov           = write_vcov
          )
          result$skipped <- isTRUE(result$estimation$skipped_existing)
        }
        
        if (run_descriptive) {
          result$descriptive <- run_one_descriptive(
            data_source      = long_data_source,
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
            dir_out          = dir_out,
            run_id           = run_id,
            run_stub         = run_stub,
            group_palette    = group_palette,
            treated_year_var = desc_args$treated_year_var,
            control_year_var = desc_args$control_year_var,
            skip_existing    = skip_existing
          )
        }
        
      }, error = function(e) {
        result$error <<- conditionMessage(e)
        message(glue::glue("[ERROR] run_id = {run_id}: {conditionMessage(e)}"))
      })
      
      if (verbose_timing) {
        elapsed <- (proc.time() - t_start)[["elapsed"]]
        status  <- if (!is.null(result$error)) "FAILED" else if (result$skipped) "SKIPPED" else "OK"
        message(glue::glue("[{timestamp_now()}] {run_id} | {status} | {round(elapsed, 1)}s"))
      }
      
      result
    },
    .progress = .progress
  )
  
  failed <- purrr::keep(run_results, \(r) !is.null(r$error))
  if (length(failed) > 0) {
    message(glue::glue(
      "\n{length(failed)} run(s) failed:\n",
      paste(purrr::map_chr(failed, "run_id"), collapse = "\n")
    ))
  }
  
  saveRDS(
    list(
      dir_out        = dir_out,
      timestamp      = ts,
      n_runs_planned = nrow(run_grid),
      n_runs_ok      = sum(purrr::map_lgl(run_results, \(r) is.null(r$error))),
      n_runs_failed  = length(failed),
      n_runs_skipped = sum(purrr::map_lgl(run_results, \(r) isTRUE(r$skipped))),
      failed_run_ids = purrr::map_chr(failed, "run_id")
    ),
    info_file
  )
  
  invisible(list(
    run_grid       = run_grid,
    run_results    = run_results,
    metadata_files = list(
      run_grid_rds    = run_grid_file,
      specs_snapshot  = specs_file,
      experiment_info = info_file
    )
  ))
}


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 9 — Single-run estimation worker
# ══════════════════════════════════════════════════════════════════════════════

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
                               weights_col,
                               weights_parquet_path,
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
  
  coef_file     <- file.path(run_dir, "coef.parquet")
  subgroup_file <- file.path(run_dir, "subgroup.parquet")
  support_file  <- file.path(run_dir, "support.parquet")
  vcov_file     <- file.path(run_dir, "vcov.parquet")
  registry_file <- file.path(run_dir, "registry.parquet")
  run_spec_file <- file.path(run_dir, "run_spec.rds")
  
  required_files <- c(coef_file, subgroup_file, support_file, registry_file, run_spec_file)
  if (write_vcov) required_files <- c(required_files, vcov_file)
  
  if (skip_existing && all(file.exists(required_files))) {
    return(list(
      coef_file        = coef_file,
      subgroup_file    = subgroup_file,
      support_file     = support_file,
      vcov_file        = if (write_vcov) vcov_file else NULL,
      registry_file    = registry_file,
      run_spec_file    = run_spec_file,
      skipped_existing = TRUE
    ))
  }
  
  fit_started <- Sys.time()
  
  # ── Load long panel ────────────────────────────────────────────────────────
  df <- load_arrow_data(data_source, data_filter)
  n_rows_read         <- nrow(df)
  n_rows_after_filter <- nrow(df)
  
  # ── Join weights when applicable ───────────────────────────────────────────
  if (!is.na(weights_col) && !is.null(weights_parquet_path)) {
    weights_tbl <- arrow::read_parquet(weights_parquet_path) |>
      dplyr::select(dplyr::all_of(c(dataset_spec$unit_id, weights_col)))
    df <- df |> dplyr::left_join(weights_tbl, by = dataset_spec$unit_id)
    n_na_weights <- sum(is.na(df[[weights_col]]))
    if (n_na_weights > 0) {
      warning(glue::glue(
        "[{run_id}] {n_na_weights} rows have NA weights after join and will be dropped by feols()."
      ))
    }
  }
  
  # ── Apply grouping ─────────────────────────────────────────────────────────
  df_grouped <- apply_treatment_grouping(
    df = df, group_fun = group_fun, group_col = group_col, group_args = group_args
  )
  
  # ── Build formula and identify needed columns ──────────────────────────────
  model_formula <- build_model_formula(formula_template, outcome)
  
  needed_cols <- get_needed_columns(
    formula     = model_formula,
    trt_col     = dataset_spec$trt_col,
    group_col   = group_col,
    unit_id     = dataset_spec$unit_id,
    event_id    = dataset_spec$event_id,
    vcov_vars   = vcov_specs$vcov_vars,
    weights_col = weights_col
  )
  
  check_required_columns(df_grouped, needed_cols, context = run_id)
  df_grouped <- df_grouped |> dplyr::select(dplyr::all_of(needed_cols))
  
  # ── Fit subgroup models ────────────────────────────────────────────────────
  subgroup_models <- fit_subgroup_models(
    df             = df_grouped,
    group_col      = group_col,
    dataset_spec   = dataset_spec,
    formula        = model_formula,
    estimator_type = estimator_type,
    weights_col    = weights_col
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
    weights_col      = weights_col,
    group_palette    = group_palette
  )
  
  subgroup_tbl <- make_subgroup_run_summary(models = subgroup_models, meta = shared_meta) |>
    dplyr::mutate(run_id = run_id)
  
  support_tbl <- make_event_time_support_run_summary(models = subgroup_models, meta = shared_meta) |>
    dplyr::mutate(run_id = run_id)
  
  coef_tbl <- extract_all_coef_tables(
    models       = subgroup_models,
    meta         = shared_meta,
    vcov_specs   = vcov_specs,
    term_pattern = term_pattern,
    ci_level     = ci_level
  ) |> dplyr::mutate(run_id = run_id)
  
  vcov_tbl <- if (write_vcov) {
    extract_all_vcov_tables(
      models = subgroup_models, meta = shared_meta, vcov_specs = vcov_specs
    ) |> dplyr::mutate(run_id = run_id)
  } else {
    NULL
  }
  
  registry_tbl <- make_estimation_registry(
    run_id               = run_id,
    subset_id            = subset_id,
    outcome              = outcome,
    group_id             = group_id,
    model_id             = model_id,
    group_col            = group_col,
    dataset_spec         = dataset_spec,
    formula_template     = formula_template,
    estimator_type       = estimator_type,
    term_pattern         = term_pattern,
    weights_col          = weights_col,
    weights_parquet_path = weights_parquet_path,
    data_filter          = data_filter,
    data_source          = data_source,
    vcov_specs           = vcov_specs,
    group_args           = group_args,
    n_rows_read          = n_rows_read,
    n_rows_after_filter  = n_rows_after_filter,
    subgroup_tbl         = subgroup_tbl,
    support_tbl          = support_tbl,
    fit_started          = fit_started,
    fit_finished         = fit_finished,
    coef_file            = coef_file,
    subgroup_file        = subgroup_file,
    support_file         = support_file,
    vcov_file            = if (write_vcov) vcov_file else NA_character_,
    registry_file        = registry_file,
    run_spec_file        = run_spec_file
  )
  
  run_spec <- list(
    run_id               = run_id,
    run_stub             = run_stub,
    subset_id            = subset_id,
    outcome              = outcome,
    group_id             = group_id,
    model_id             = model_id,
    group_col            = group_col,
    dataset_spec         = dataset_spec,
    formula_template     = formula_template,
    estimator_type       = estimator_type,
    term_pattern         = term_pattern,
    weights_col          = weights_col,
    weights_parquet_path = weights_parquet_path,
    data_filter          = data_filter,
    data_source          = data_source,
    group_fun            = group_fun,
    group_args           = group_args,
    vcov_specs           = vcov_specs,
    group_palette        = group_palette,
    ci_level             = ci_level,
    write_vcov           = write_vcov
  )
  
  arrow::write_parquet(subgroup_tbl, subgroup_file)
  arrow::write_parquet(support_tbl,  support_file)
  arrow::write_parquet(coef_tbl,     coef_file)
  arrow::write_parquet(registry_tbl, registry_file)
  if (write_vcov && !is.null(vcov_tbl)) arrow::write_parquet(vcov_tbl, vcov_file)
  saveRDS(run_spec, run_spec_file)
  
  rm(subgroup_models, df, df_grouped,
     subgroup_tbl, support_tbl, coef_tbl, vcov_tbl, registry_tbl, run_spec)
  gc()
  
  list(
    coef_file        = coef_file,
    subgroup_file    = subgroup_file,
    support_file     = support_file,
    vcov_file        = if (write_vcov) vcov_file else NULL,
    registry_file    = registry_file,
    run_spec_file    = run_spec_file,
    skipped_existing = FALSE
  )
}


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 10 — Single-run descriptive worker
# ══════════════════════════════════════════════════════════════════════════════

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
    df = df, group_fun = group_fun, group_col = group_col, group_args = group_args
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
# SECTION 11 — Rebuild (merge across runs)
# ══════════════════════════════════════════════════════════════════════════════

rebuild_estimation_tables <- function(dir_out, write_csv = TRUE, recursive = TRUE) {
  
  dir_by_run <- file.path(dir_out, "tables", "by_run")
  dir_all    <- file.path(dir_out, "tables", "all")
  
  if (!dir.exists(dir_by_run)) stop("by_run directory does not exist: ", dir_by_run)
  dir.create(dir_all, recursive = TRUE, showWarnings = FALSE)
  
  find_parquets <- function(pattern) {
    list.files(dir_by_run, pattern = pattern, recursive = recursive, full.names = TRUE)
  }
  
  coef_files     <- find_parquets("^coef\\.parquet$")
  subgroup_files <- find_parquets("^subgroup\\.parquet$")
  support_files  <- find_parquets("^support\\.parquet$")
  registry_files <- find_parquets("^registry\\.parquet$")
  run_spec_files <- list.files(dir_by_run, "^run_spec\\.rds$", recursive = recursive, full.names = TRUE)
  vcov_files     <- find_parquets("^vcov\\.parquet$")
  
  for (nm in c("coef_files", "subgroup_files", "support_files", "registry_files")) {
    if (length(get(nm)) == 0) {
      stop("No ", sub("_files", ".parquet", nm), " files found under: ", dir_by_run)
    }
  }
  
  read_and_dedup <- function(files, id_col) {
    out <- purrr::map(files, \(f) {
      x <- arrow::read_parquet(f)
      
      x[] <- lapply(x, function(col) {
        if (inherits(col, "json") || is.list(col)) {
          as.character(col)
        } else {
          col
        }
      })
      
      x$.mtime <- file.mtime(f)
      x
    })
    
    dplyr::bind_rows(out) |>
      dplyr::group_by(dplyr::across(dplyr::all_of(id_col))) |>
      dplyr::filter(.mtime == max(.mtime)) |>
      dplyr::ungroup() |>
      dplyr::select(-.mtime)
  }
  
  coef_tbl     <- read_and_dedup(coef_files,     c("run_id", "subgroup", "term", "vcov_id"))
  subgroup_tbl <- read_and_dedup(subgroup_files,  c("run_id", "subgroup"))
  support_tbl  <- read_and_dedup(support_files,   c("run_id", "subgroup", "event_time"))
  registry_tbl <- read_and_dedup(registry_files,  "run_id")
  
  vcov_tbl <- if (length(vcov_files) > 0) {
    read_and_dedup(vcov_files, c("run_id", "subgroup", "term_i", "term_j", "vcov_id"))
  } else {
    NULL
  }
  
  spec_index_tbl <- NULL
  if (length(run_spec_files) > 0) {
    spec_index_tbl <- purrr::map_dfr(run_spec_files, \(f) {
      spec <- readRDS(f)
      tibble::tibble(
        run_id               = spec$run_id               %||% NA_character_,
        subset_id            = spec$subset_id             %||% NA_character_,
        outcome              = spec$outcome               %||% NA_character_,
        group_id             = spec$group_id              %||% NA_character_,
        model_id             = spec$model_id              %||% NA_character_,
        group_col            = spec$group_col             %||% NA_character_,
        estimator_type       = spec$estimator_type        %||% NA_character_,
        formula_template     = spec$formula_template      %||% NA_character_,
        term_pattern         = spec$term_pattern          %||% NA_character_,
        weights_col          = spec$weights_col           %||% NA_character_,
        weights_parquet_path = spec$weights_parquet_path  %||% NA_character_,
        data_filter          = data_filter_to_chr(spec$data_filter),
        data_source          = paste(unlist(spec$data_source), collapse = " | "),
        group_args_json      = as.character(serialize_object_json(spec$group_args)),
        vcov_ids             = paste(spec$vcov_specs$vcov_id, collapse = " | "),
        run_spec_file        = f,
        .mtime               = file.mtime(f)
      )
    }) |>
      dplyr::group_by(run_id) |>
      dplyr::filter(.mtime == max(.mtime)) |>
      dplyr::ungroup() |>
      dplyr::select(-.mtime)
  }
  
  write_pair <- function(tbl, stem) {
    if (is.null(tbl)) return(list(parquet = NULL, csv = NULL))
    pq  <- file.path(dir_all, paste0(stem, ".parquet"))
    csv <- if (write_csv) file.path(dir_all, paste0(stem, ".csv")) else NULL
    arrow::write_parquet(tbl, pq)
    if (write_csv) readr::write_csv(tbl, csv)
    list(parquet = pq, csv = csv)
  }
  
  invisible(list(
    coef_expanded_vcov   = coef_tbl,
    subgroup_run_summary = subgroup_tbl,
    event_time_support   = support_tbl,
    run_registry         = registry_tbl,
    vcov_matrices        = vcov_tbl,
    run_spec_index       = spec_index_tbl,
    files = list(
      coef_expanded_vcov   = write_pair(coef_tbl,      "coef_expanded_vcov"),
      subgroup_run_summary = write_pair(subgroup_tbl,   "subgroup_run_summary"),
      event_time_support   = write_pair(support_tbl,    "event_time_support"),
      run_registry         = write_pair(registry_tbl,   "run_registry"),
      vcov_matrices        = write_pair(vcov_tbl,       "vcov_matrices"),
      run_spec_index       = write_pair(spec_index_tbl, "run_spec_index")
    )
  ))
}


rebuild_descriptive_tables <- function(dir_out, write_csv = TRUE, recursive = TRUE) {
  
  dir_by_run <- file.path(dir_out, "descriptive", "by_run")
  dir_all    <- file.path(dir_out, "descriptive", "all")
  
  if (!dir.exists(dir_by_run)) stop("descriptive by_run directory does not exist: ", dir_by_run)
  dir.create(dir_all, recursive = TRUE, showWarnings = FALSE)
  
  traj_files     <- list.files(dir_by_run, "^event_time_trajectory\\.parquet$",
                               recursive = recursive, full.names = TRUE)
  registry_files <- list.files(dir_by_run, "^registry\\.parquet$",
                               recursive = recursive, full.names = TRUE)
  
  if (length(traj_files) == 0)     stop("No event_time_trajectory.parquet found under: ", dir_by_run)
  if (length(registry_files) == 0) stop("No registry.parquet found under: ", dir_by_run)
  
  read_and_dedup <- function(files, id_col) {
    purrr::map_dfr(
      files,
      \(f) arrow::read_parquet(f) |> dplyr::mutate(.mtime = file.mtime(f))
    ) |>
      dplyr::group_by(dplyr::across(dplyr::all_of(id_col))) |>
      dplyr::filter(.mtime == max(.mtime)) |>
      dplyr::ungroup() |>
      dplyr::select(-.mtime)
  }
  
  traj_tbl     <- read_and_dedup(traj_files,    c("run_id", "subgroup", "series", "event_time"))
  registry_tbl <- read_and_dedup(registry_files, "run_id")
  
  write_pair <- function(tbl, stem) {
    pq  <- file.path(dir_all, paste0(stem, ".parquet"))
    csv <- if (write_csv) file.path(dir_all, paste0(stem, ".csv")) else NULL
    arrow::write_parquet(tbl, pq)
    if (write_csv) readr::write_csv(tbl, csv)
    list(parquet = pq, csv = csv)
  }
  
  invisible(list(
    event_time_trajectory = traj_tbl,
    run_registry          = registry_tbl,
    files = list(
      event_time_trajectory = write_pair(traj_tbl,     "event_time_trajectory"),
      run_registry          = write_pair(registry_tbl, "run_registry")
    )
  ))
}


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 12 — Model fitting internals
# ══════════════════════════════════════════════════════════════════════════════

fit_subgroup_models <- function(df, group_col, dataset_spec, formula,
                                estimator_type, weights_col = NA_character_, ...) {
  
  trt_col  <- dataset_spec$trt_col
  unit_id  <- dataset_spec$unit_id
  event_id <- dataset_spec$event_id
  
  groups <- df |>
    dplyr::filter(!is.na(.data[[group_col]])) |>
    dplyr::distinct(.data[[group_col]]) |>
    dplyr::pull(.data[[group_col]]) |>
    as.character() |>
    sort()
  
  if (length(groups) == 0) stop(glue::glue("No non-NA values in group column '{group_col}'."))
  
  models <- purrr::map(groups, \(g) fit_one_subgroup_model(
    df = df, group_col = group_col, group_value = g, trt_col = trt_col,
    formula = formula, unit_id = unit_id, event_id = event_id,
    estimator_type = estimator_type, weights_col = weights_col, ...
  ))
  
  names(models) <- groups
  models
}


fit_one_subgroup_model <- function(df, group_col, group_value, trt_col,
                                   formula, unit_id, event_id,
                                   estimator_type, weights_col = NA_character_, ...) {
  
  group_sym <- rlang::sym(group_col)
  trt_sym   <- rlang::sym(trt_col)
  
  df_model <- dplyr::bind_rows(
    df |> dplyr::filter(!!trt_sym == 0),
    df |> dplyr::filter(!!trt_sym == 1, !!group_sym == group_value)
  )
  
  feols_weights <- if (!is.na(weights_col) && weights_col %in% names(df_model)) {
    stats::as.formula(paste("~", weights_col))
  } else {
    NULL
  }
  
  model  <- fixest::feols(formula, data = df_model, weights = feols_weights, ...)
  df_est <- df_model[fixest::obs(model), , drop = FALSE]
  
  attr(model, "group_value")    <- group_value
  attr(model, "group_col")      <- group_col
  attr(model, "unit_id")        <- unit_id
  attr(model, "event_id")       <- normalize_optional_colname(event_id)
  attr(model, "estimator_type") <- estimator_type
  attr(model, "weights_col")    <- normalize_optional_colname(weights_col)
  
  if (estimator_type == "sunab") {
    sv <- parse_sunab_vars(formula)
    attr(model, "cohort_var")       <- sv$cohort_var
    attr(model, "time_var")         <- sv$time_var
    attr(model, "event_time_support") <- make_event_time_support_one_model(
      df_est = df_est, trt_col = trt_col, unit_id = unit_id,
      event_id = event_id, cohort_var = sv$cohort_var, time_var = sv$time_var
    )
  } else {
    attr(model, "cohort_var")         <- NA_character_
    attr(model, "time_var")           <- NA_character_
    attr(model, "event_time_support") <- empty_event_time_support()
  }
  
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
# SECTION 13 — Table extraction internals
# ══════════════════════════════════════════════════════════════════════════════

extract_all_coef_tables <- function(models, meta, vcov_specs, term_pattern, ci_level) {
  purrr::imap_dfr(models, \(model, group_nm) {
    purrr::pmap_dfr(
      list(vcov = vcov_specs$vcov, vcov_id = vcov_specs$vcov_id, vcov_label = vcov_specs$vcov_label),
      \(vcov, vcov_id, vcov_label) {
        extract_coef_table_one_vcov(
          model = model, vcov = vcov, vcov_id = vcov_id, vcov_label = vcov_label,
          meta = meta, term_pattern = term_pattern, ci_level = ci_level
        )
      }
    )
  })
}


extract_coef_table_one_vcov <- function(model, vcov, vcov_id, vcov_label,
                                        meta, term_pattern, ci_level) {
  
  ct     <- fixest::coeftable(model, vcov = vcov)
  ct_tbl <- tibble::as_tibble(as.data.frame(ct), rownames = "term")
  
  se_col <- if ("Std. Error" %in% names(ct_tbl)) "Std. Error" else names(ct_tbl)[3]
  t_col  <- names(ct_tbl)[stringr::str_detect(names(ct_tbl), "^t value$|^t-value$")]
  p_col  <- names(ct_tbl)[stringr::str_detect(names(ct_tbl), "^Pr\\(>\\|t\\|\\)$|^p-value$")]
  z_val  <- stats::qnorm(1 - (1 - ci_level) / 2)
  
  group_value        <- attr(model, "group_value")
  event_time_support <- attr(model, "event_time_support")
  subgroup_run_id    <- make_subgroup_run_id(meta, group_value)
  
  out <- ct_tbl |>
    dplyr::mutate(
      subset_id            = meta$subset_id,
      outcome              = meta$outcome,
      group_id             = meta$group_id,
      group_col            = meta$group_col,
      subgroup             = group_value,
      model_id             = meta$model_id,
      estimator_type       = meta$estimator_type,
      weights_col          = meta$weights_col %||% NA_character_,
      vcov_id              = vcov_id,
      vcov_label           = vcov_label,
      formula_template     = meta$formula_template,
      estimate             = .data[["Estimate"]],
      std_error            = .data[[se_col]],
      ci_lower             = estimate - z_val * std_error,
      ci_upper             = estimate + z_val * std_error,
      term_matches_pattern = stringr::str_detect(term, term_pattern),
      term_value           = extract_first_number(term),
      n_treated_units      = attr(model, "n_treated_units"),
      n_treated_events     = attr(model, "n_treated_events"),
      n_control_units      = attr(model, "n_control_units"),
      n_total_units        = attr(model, "n_total_units"),
      n_rows_model_data    = attr(model, "n_rows_model_data"),
      subgroup_color       = get_group_color(group_value, meta$group_palette),
      subgroup_run_id      = subgroup_run_id,
      coef_run_id          = as.character(glue::glue("{subgroup_run_id}__{vcov_id}"))
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
      dplyr::mutate(n_ptids = NA_integer_, n_fireids = NA_integer_, n_rows_treated = NA_integer_)
  }
  
  out |>
    dplyr::select(
      subset_id, outcome, group_id, group_col, subgroup, model_id,
      estimator_type, weights_col, vcov_id, vcov_label, formula_template,
      term, estimate, std_error, t_value, p_value, ci_lower, ci_upper,
      term_matches_pattern, term_value,
      n_ptids, n_fireids, n_rows_treated,
      n_treated_units, n_treated_events, n_control_units,
      n_total_units, n_rows_model_data,
      subgroup_color, subgroup_run_id, coef_run_id
    )
}


extract_all_vcov_tables <- function(models, meta, vcov_specs) {
  purrr::imap_dfr(models, \(model, group_nm) {
    purrr::pmap_dfr(
      list(vcov = vcov_specs$vcov, vcov_id = vcov_specs$vcov_id, vcov_label = vcov_specs$vcov_label),
      \(vcov, vcov_id, vcov_label) {
        extract_vcov_table_one_vcov(
          model = model, vcov = vcov, vcov_id = vcov_id,
          vcov_label = vcov_label, meta = meta
        )
      }
    )
  })
}


extract_vcov_table_one_vcov <- function(model, vcov, vcov_id, vcov_label, meta) {
  vc     <- stats::vcov(model, vcov = vcov)
  vc_mat <- as.matrix(vc)
  if (is.null(rownames(vc_mat)) || is.null(colnames(vc_mat))) {
    stop("vcov matrix must have row and column names.")
  }
  group_value     <- attr(model, "group_value")
  subgroup_run_id <- make_subgroup_run_id(meta, group_value)
  
  expand.grid(term_i = rownames(vc_mat), term_j = colnames(vc_mat),
              stringsAsFactors = FALSE) |>
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
      vcov_run_id      = as.character(glue::glue("{subgroup_run_id}__{vcov_id}"))
    ) |>
    dplyr::select(
      subset_id, outcome, group_id, group_col, subgroup, model_id,
      vcov_id, vcov_label, formula_template,
      term_i, term_j, vcov_value, subgroup_run_id, vcov_run_id
    )
}


make_subgroup_run_summary <- function(models, meta) {
  purrr::imap_dfr(models, \(model, group_nm) {
    tibble::tibble(
      subset_id         = meta$subset_id,
      outcome           = meta$outcome,
      group_id          = meta$group_id,
      group_col         = meta$group_col,
      subgroup          = group_nm,
      model_id          = meta$model_id,
      formula_template  = meta$formula_template,
      estimator_type    = meta$estimator_type,
      weights_col       = meta$weights_col %||% NA_character_,
      n_treated_units   = attr(model, "n_treated_units"),
      n_treated_events  = attr(model, "n_treated_events"),
      n_control_units   = attr(model, "n_control_units"),
      n_total_units     = attr(model, "n_total_units"),
      n_rows_model_data = attr(model, "n_rows_model_data"),
      subgroup_color    = get_group_color(group_nm, meta$group_palette),
      subgroup_run_id   = make_subgroup_run_id(meta, group_nm)
    )
  })
}


make_event_time_support_run_summary <- function(models, meta) {
  purrr::imap_dfr(models, \(model, group_nm) {
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
      return(tibble::as_tibble(c(
        base_cols,
        list(event_time = numeric(), n_ptids = integer(),
             n_fireids = integer(), n_rows_treated = integer())
      )))
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
  })
}


make_estimation_registry <- function(run_id, subset_id, outcome, group_id,
                                     model_id, group_col, dataset_spec,
                                     formula_template, estimator_type,
                                     term_pattern, weights_col,
                                     weights_parquet_path,
                                     data_filter, data_source,
                                     vcov_specs, group_args,
                                     n_rows_read, n_rows_after_filter,
                                     subgroup_tbl, support_tbl,
                                     fit_started, fit_finished,
                                     coef_file, subgroup_file, support_file,
                                     vcov_file, registry_file, run_spec_file) {
  
  model_formula <- build_model_formula(formula_template, outcome)
  
  tibble::tibble(
    run_id               = run_id,
    subset_id            = subset_id,
    outcome              = outcome,
    group_id             = group_id,
    model_id             = model_id,
    group_col            = group_col,
    unit_id              = dataset_spec$unit_id,
    event_id             = dataset_spec$event_id,
    trt_col              = dataset_spec$trt_col,
    cohort_var           = dataset_spec$cohort_var,
    time_var             = dataset_spec$time_var,
    formula_template     = formula_template,
    formula_resolved     = paste(deparse(model_formula), collapse = " "),
    estimator_type       = estimator_type,
    term_pattern         = term_pattern,
    weights_col          = weights_col          %||% NA_character_,
    weights_parquet_path = weights_parquet_path %||% NA_character_,
    data_filter          = data_filter_to_chr(data_filter),
    input_source         = paste(unlist(data_source), collapse = " | "),
    group_args_json      = as.character(serialize_object_json(group_args)),
    vcov_ids             = paste(vcov_specs$vcov_id,    collapse = " | "),
    vcov_labels          = paste(vcov_specs$vcov_label, collapse = " | "),
    vcov_vars_json       = as.character(serialize_object_json(vcov_specs$vcov_vars)),
    n_rows_read          = n_rows_read,
    n_rows_after_filter  = n_rows_after_filter,
    n_subgroups          = dplyr::n_distinct(subgroup_tbl$subgroup),
    subgroup_levels      = paste(sort(unique(subgroup_tbl$subgroup)), collapse = ","),
    n_support_rows       = nrow(support_tbl),
    fit_started          = fit_started,
    fit_finished         = fit_finished,
    coef_file            = coef_file,
    subgroup_file        = subgroup_file,
    support_file         = support_file,
    vcov_file            = vcov_file,
    registry_file        = registry_file,
    run_spec_file        = run_spec_file
  )
}


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 14 — Descriptive internals
# ══════════════════════════════════════════════════════════════════════════════

make_descriptive_trajectory <- function(df, meta, dataset_spec,
                                        treated_year_var, control_year_var) {
  
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
  
  purrr::map_dfr(subgroup_levels, \(grp) {
    tbl <- make_descriptive_one_subgroup(
      df = df, group_col = group_col, group_value = grp,
      trt_col = trt_col, outcome = outcome, unit_id = unit_id,
      event_id = event_id, time_var = time_var,
      treated_year_var = treated_year_var, control_year_var = control_year_var
    )
    srp <- make_subgroup_run_id(meta, grp)
    tbl |>
      dplyr::mutate(
        subset_id          = meta$subset_id,
        outcome            = meta$outcome,
        group_id           = meta$group_id,
        group_col          = group_col,
        model_id           = meta$model_id,
        formula_template   = meta$formula_template,
        treated_year_var   = treated_year_var,
        control_year_var   = control_year_var,
        time_var           = time_var,
        subgroup_color     = get_group_color(grp, meta$group_palette),
        subgroup_run_id    = srp,
        descriptive_run_id = as.character(glue::glue("{srp}__{tbl$series}"))
      ) |>
      dplyr::select(
        subset_id, outcome, group_id, group_col, subgroup,
        model_id, formula_template,
        treated_year_var, control_year_var, time_var,
        event_time, series,
        mean_outcome, sd_outcome, n_rows, n_ptids, n_fireids,
        subgroup_color, subgroup_run_id, descriptive_run_id
      )
  })
}


make_descriptive_one_subgroup <- function(df, group_col, group_value, trt_col,
                                          outcome, unit_id, event_id, time_var,
                                          treated_year_var, control_year_var) {
  
  event_id <- normalize_optional_colname(event_id)
  needed   <- unique(stats::na.omit(c(group_col, trt_col, outcome, unit_id,
                                      event_id, time_var, treated_year_var,
                                      control_year_var)))
  check_required_columns(df, needed, context = glue::glue("descriptive: {group_value}"))
  
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
# SECTION 15 — Sun-Abraham support helpers
# ══════════════════════════════════════════════════════════════════════════════

parse_sunab_vars <- function(formula) {
  txt <- gsub("\\s+", " ", paste(deparse(formula), collapse = " "))
  m   <- stringr::str_match(
    txt, "sunab\\s*\\(\\s*([^,]+?)\\s*,\\s*([^,\\)]+?)\\s*(?:,|\\))"
  )
  if (length(m) == 0 || all(is.na(m[1, ]))) {
    warning("Could not parse sunab() vars from formula: ", txt,
            ". event_time_support will be empty.", call. = FALSE)
    return(list(cohort_var = NA_character_, time_var = NA_character_))
  }
  list(cohort_var = trimws(m[1, 2]), time_var = trimws(m[1, 3]))
}


make_event_time_support_one_model <- function(df_est, trt_col, unit_id,
                                              event_id   = NA_character_,
                                              cohort_var, time_var) {
  
  event_id   <- normalize_optional_colname(event_id)
  cohort_var <- normalize_optional_colname(cohort_var)
  time_var   <- normalize_optional_colname(time_var)
  
  if (is.na(cohort_var) || is.na(time_var)) return(empty_event_time_support())
  
  needed <- unique(stats::na.omit(c(trt_col, unit_id, event_id, cohort_var, time_var)))
  check_required_columns(df_est, needed, context = "event_time_support")
  
  d_treated <- df_est |>
    dplyr::filter(.data[[trt_col]] == 1) |>
    dplyr::mutate(
      event_time = suppressWarnings(
        as.numeric(.data[[time_var]]) - as.numeric(.data[[cohort_var]])
      )
    ) |>
    dplyr::filter(!is.na(event_time), is.finite(event_time))
  
  if (nrow(d_treated) == 0) return(empty_event_time_support())
  
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
# SECTION 16 — Arrow data loading helpers
# ══════════════════════════════════════════════════════════════════════════════

open_arrow_source <- function(data_source) {
  arrow::open_dataset(unlist(data_source), format = "parquet")
}

# load_arrow_data <- function(data_source, data_filter) {
#   ds <- open_arrow_source(data_source)
#   if (!is.null(data_filter)) {
#     filter_quo <- rlang::new_quosure(rlang::f_rhs(data_filter), env = rlang::f_env(data_filter))
#     ds <- ds |> dplyr::filter(!!filter_quo)
#   }
#   ds |> dplyr::collect()
# }

load_arrow_data <- function(data_source, data_filter = NULL) {
  ds <- open_arrow_source(data_source)
  
  if (is.null(data_filter)) {
    return(dplyr::collect(ds))
  }
  
  filter_expr <- rlang::f_rhs(data_filter)
  
  tryCatch(
    {
      ds |>
        dplyr::filter(!!filter_expr) |>
        dplyr::collect()
    },
    error = function(e) {
      message(
        "[load_arrow_data] Arrow could not push down filter — collecting full dataset and filtering in R. Filter: ",
        data_filter_to_chr(data_filter)
      )
      
      df <- dplyr::collect(ds)
      
      df |>
        dplyr::filter(!!filter_expr)
    }
  )
}


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 17 — Treatment grouping helper
# ══════════════════════════════════════════════════════════════════════════════

#' Apply a grouping function and validate its output
#'
#' Injects group_col and any extra args (including include_control when called
#' from the weighting pipeline) into group_fun's argument list.
apply_treatment_grouping <- function(df, group_fun, group_col, group_args) {
  df_out <- do.call(group_fun, c(list(df = df, group_col = group_col), group_args))
  if (!group_col %in% names(df_out)) {
    stop(glue::glue(
      "group_fun did not produce expected column '{group_col}'. ",
      "Check that your grouping function uses the group_col argument."
    ))
  }
  df_out
}


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 18 — Utility functions
# ══════════════════════════════════════════════════════════════════════════════

build_model_formula <- function(formula_template, outcome) {
  if (!is.character(formula_template) || length(formula_template) != 1) {
    stop("formula_template must be a length-1 character string.")
  }
  stats::as.formula(glue::glue(formula_template, outcome = outcome))
}

normalize_optional_colname <- function(x) {
  if (is.null(x) || length(x) == 0) return(NA_character_)
  if (length(x) > 1) {
    warning("normalize_optional_colname: length > 1 input; using first element only.")
    x <- x[1]
  }
  if (is.na(x)) return(NA_character_)
  as.character(x)
}

get_needed_columns <- function(formula, trt_col, group_col, unit_id,
                               event_id    = NA_character_,
                               vcov_vars   = NULL,
                               weights_col = NA_character_) {
  cols <- c(
    all.vars(formula), trt_col, group_col, unit_id, event_id,
    unlist(vcov_vars, use.names = FALSE),
    if (!is.na(weights_col)) weights_col else NULL
  )
  cols |> unique() |> stats::na.omit() |> as.character()
}

check_required_columns <- function(df, required, context = "") {
  missing <- setdiff(required, names(df))
  if (length(missing) > 0) {
    stop(
      if (nchar(context) > 0) paste0("[", context, "] ") else "",
      "Missing required columns: ", paste(missing, collapse = ", ")
    )
  }
  invisible(TRUE)
}

compute_n_distinct_optional <- function(df, col_nm) {
  col_nm <- normalize_optional_colname(col_nm)
  if (is.na(col_nm)) return(NA_integer_)
  if (!col_nm %in% names(df)) stop(glue::glue("Column '{col_nm}' not found in data."))
  dplyr::n_distinct(stats::na.omit(df[[col_nm]]))
}

safe_path_component <- function(x) {
  gsub("[^[:alnum:]_\\-]+", "_", as.character(x))
}

make_run_stub <- function(subset_id, outcome, group_id, model_id) {
  paste(safe_path_component(subset_id), safe_path_component(outcome),
        safe_path_component(group_id),  safe_path_component(model_id), sep = "__")
}

make_subgroup_run_id <- function(meta, group_value) {
  as.character(glue::glue(
    "{meta$subset_id}__{meta$outcome}__{meta$group_id}__{meta$model_id}__{group_value}"
  ))
}

extract_first_number <- function(x) {
  suppressWarnings(as.numeric(stringr::str_extract(x, "-?\\d+")))
}

dir_ensure_local <- function(paths) {
  purrr::walk(paths, \(p) {
    dir.create(p, recursive = TRUE, showWarnings = FALSE)
    if (!dir.exists(p)) stop("Failed to create directory: ", p)
  })
}

timestamp_now <- function() format(Sys.time(), "%Y%m%d_%H%M%S")

serialize_object_json <- function(x) {
  as.character(
    jsonlite::toJSON(x, auto_unbox = TRUE, null = "null", digits = NA, pretty = FALSE)
  )
}

data_filter_to_chr <- function(data_filter = NULL) {
  if (is.null(data_filter))             return(NA_character_)
  if (inherits(data_filter, "formula")) return(paste(deparse(data_filter), collapse = " "))
  as.character(data_filter)
}

get_group_color <- function(group_value, group_palette) {
  if (!is.null(group_palette) && group_value %in% names(group_palette)) {
    return(unname(group_palette[[group_value]]))
  }
  NA_character_
}

validate_data_source <- function(src, label = "data_source") {
  if (!is.character(src) || length(src) == 0) {
    stop(label, " must be a non-empty character vector of paths or a directory path.")
  }
  if (length(src) == 1 && dir.exists(src)) return(invisible(TRUE))
  missing_files <- src[!file.exists(src)]
  if (length(missing_files) > 0) {
    warning(label, ": these files do not exist yet:\n", paste(missing_files, collapse = "\n"))
  }
  invisible(TRUE)
}

validate_dataset_spec <- function(spec) {
  required <- c("unit_id", "time_var", "trt_col", "cohort_var", "event_id")
  missing  <- setdiff(required, names(spec))
  if (length(missing) > 0) stop("dataset_spec is missing fields: ", paste(missing, collapse = ", "))
  invisible(TRUE)
}

validate_spec_table <- function(tbl, required_cols, name) {
  if (!is.data.frame(tbl)) stop(name, " must be a data frame / tibble.")
  missing <- setdiff(required_cols, names(tbl))
  if (length(missing) > 0) {
    stop(name, " is missing columns: ", paste(missing, collapse = ", "))
  }
  invisible(TRUE)
}

build_equality_filter <- function(col, value) {
  rlang::new_formula(lhs = NULL, rhs = call("==", as.name(col), value), env = baseenv())
}

combine_filters <- function(filter_a, filter_b) {
  if (is.null(filter_a)) return(filter_b)
  if (is.null(filter_b)) return(filter_a)
  rlang::new_formula(
    lhs = NULL,
    rhs = call("&", rlang::f_rhs(filter_a), rlang::f_rhs(filter_b)),
    env = rlang::f_env(filter_a)
  )
}

`%||%` <- function(a, b) if (!is.null(a)) a else b