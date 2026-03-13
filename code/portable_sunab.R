# A set of portable functions for running large-scale Sun-Abraham analyses

required_pkgs <- c(
  "arrow",
  "dplyr",
  "fixest",
  "glue",
  "purrr",
  "readr",
  "rlang",
  "stringr",
  "tibble",
  "tidyr"
)

missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_pkgs) > 0) {
  stop(
    "Missing required packages: ",
    paste(missing_pkgs, collapse = ", "),
    call. = FALSE
  )
}


build_model_formula <- function(formula_template, outcome) {
  if (!is.character(formula_template) || length(formula_template) != 1) {
    stop("formula_template must be a length-1 character string.")
  }
  
  stats::as.formula(
    glue::glue(formula_template, outcome = outcome)
  )
}


safe_path_component <- function(x) {
  x |>
    as.character() |>
    gsub("[^[:alnum:]_\\-]+", "_", x = _)
}

apply_optional_filter <- function(df, data_filter = NULL) {
  if (is.null(data_filter)) {
    return(df)
  }
  
  if (inherits(data_filter, "formula")) {
    filter_quo <- rlang::new_quosure(
      rlang::f_rhs(data_filter),
      env = rlang::f_env(data_filter)
    )
    
    return(
      df |>
        dplyr::filter(!!filter_quo)
    )
  }
  
  stop("data_filter must be NULL or a one-sided formula, e.g. ~ year >= 1993")
}

get_needed_columns <- function(formula,
                               trt_col,
                               subgroup_col,
                               unit_id,
                               event_id = NA_character_,
                               vcov_vars = NULL) {
  
  cols <- c(
    all.vars(formula),
    trt_col,
    subgroup_col,
    unit_id,
    event_id,
    unlist(vcov_vars, use.names = FALSE)
  )
  
  cols |>
    unique() |>
    stats::na.omit() |>
    as.character()
}

normalize_optional_colname <- function(x) {
  if (is.null(x)) {
    return(NA_character_)
  }
  
  if (length(x) == 0) {
    return(NA_character_)
  }
  
  if (is.na(x[1])) {
    return(NA_character_)
  }
  
  as.character(x[1])
}


compute_n_distinct_optional <- function(df, col_nm) {
  col_nm <- normalize_optional_colname(col_nm)
  
  if (is.na(col_nm)) {
    return(NA_integer_)
  }
  
  if (!col_nm %in% names(df)) {
    stop(glue::glue("Column '{col_nm}' not found in data."))
  }
  
  dplyr::n_distinct(stats::na.omit(df[[col_nm]]))
}


extract_first_number <- function(x) {
  out <- stringr::str_extract(x, "-?\\d+")
  suppressWarnings(as.numeric(out))
}


dir_ensure_local <- function(paths) {
  purrr::walk(
    paths,
    \(p) dir.create(p, recursive = TRUE, showWarnings = FALSE)
  )
}


set_cd_groups <- function(df,
                          b_nm,
                          d_nm,
                          cd_nm,
                          b_threshold,
                          d_threshold) {
  
  if (d_threshold < 0) {
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
  
  df_new
}


fit_one_subgroup_model <- function(df,
                                   subgroup_col,
                                   subgroup_value,
                                   trt_col,
                                   formula,
                                   unit_id,
                                   event_id = NA_character_,
                                   ...) {
  
  subgroup_sym <- rlang::sym(subgroup_col)
  trt_sym <- rlang::sym(trt_col)
  
  df_controls <- df |>
    dplyr::filter(!!trt_sym == 0)
  
  df_treated_subgroup <- df |>
    dplyr::filter(!!trt_sym == 1, !!subgroup_sym == subgroup_value)
  
  df_model <- dplyr::bind_rows(df_controls, df_treated_subgroup)
  
  model <- fixest::feols(
    formula,
    data = df_model,
    ...
  )
  
  obs_used <- fixest::obs(model)
  df_est <- df_model[obs_used, , drop = FALSE]
  
  attr(model, "subgroup") <- subgroup_value
  attr(model, "subgroup_col") <- subgroup_col
  attr(model, "unit_id") <- unit_id
  attr(model, "event_id") <- normalize_optional_colname(event_id)
  
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
  
  attr(model, "n_total_units") <- dplyr::n_distinct(df_est[[unit_id]])
  attr(model, "n_rows_model_data") <- nrow(df_est)
  
  model
}


fit_subgroup_models <- function(df,
                                subgroup_col,
                                trt_col,
                                formula,
                                unit_id,
                                event_id = NA_character_,
                                groups = NULL,
                                ...) {
  
  if (!subgroup_col %in% names(df)) {
    stop(glue::glue("Subgroup column '{subgroup_col}' not found in data."))
  }
  
  if (!trt_col %in% names(df)) {
    stop(glue::glue("Treatment column '{trt_col}' not found in data."))
  }
  
  if (!unit_id %in% names(df)) {
    stop(glue::glue("Unit ID column '{unit_id}' not found in data."))
  }
  
  event_id <- normalize_optional_colname(event_id)
  if (!is.na(event_id) && !event_id %in% names(df)) {
    stop(glue::glue("Event ID column '{event_id}' not found in data."))
  }
  
  if (is.null(groups)) {
    groups <- df |>
      dplyr::filter(!is.na(.data[[subgroup_col]])) |>
      dplyr::distinct(.data[[subgroup_col]]) |>
      dplyr::pull(.data[[subgroup_col]]) |>
      as.character() |>
      sort()
  }
  
  models <- purrr::map(
    groups,
    \(g) {
      fit_one_subgroup_model(
        df = df,
        subgroup_col = subgroup_col,
        subgroup_value = g,
        trt_col = trt_col,
        formula = formula,
        unit_id = unit_id,
        event_id = event_id,
        ...
      )
    }
  )
  
  names(models) <- groups
  models
}


extract_coef_table_one_vcov <- function(model,
                                        vcov,
                                        vcov_id,
                                        vcov_label,
                                        analysis_id,
                                        outcome,
                                        group_id,
                                        subgroup_col,
                                        model_id,
                                        formula_template,
                                        term_pattern = ".*",
                                        ci_level = 0.95,
                                        group_palette = NULL) {
  
  ct <- fixest::coeftable(model, vcov = vcov)
  ct_tbl <- tibble::as_tibble(as.data.frame(ct), rownames = "term")
  
  se_col <- if ("Std. Error" %in% names(ct_tbl)) "Std. Error" else names(ct_tbl)[3]
  t_col <- names(ct_tbl)[stringr::str_detect(names(ct_tbl), "^t value$|^t-value$")]
  p_col <- names(ct_tbl)[stringr::str_detect(names(ct_tbl), "^Pr\\(>\\|t\\|\\)$|^p-value$")]
  
  z_val <- stats::qnorm(1 - (1 - ci_level) / 2)
  
  subgroup <- attr(model, "subgroup")
  n_treated_units <- attr(model, "n_treated_units")
  n_treated_events <- attr(model, "n_treated_events")
  n_control_units <- attr(model, "n_control_units")
  n_total_units <- attr(model, "n_total_units")
  n_rows_model_data <- attr(model, "n_rows_model_data")
  
  subgroup_color <- NA_character_
  if (!is.null(group_palette) && subgroup %in% names(group_palette)) {
    subgroup_color <- unname(group_palette[[subgroup]])
  }
  
  out <- ct_tbl |>
    dplyr::mutate(
      analysis_id = analysis_id,
      outcome = outcome,
      group_id = group_id,
      subgroup_col = subgroup_col,
      subgroup = subgroup,
      model_id = model_id,
      vcov_id = vcov_id,
      vcov_label = vcov_label,
      formula_template = formula_template,
      estimate = .data[["Estimate"]],
      std_error = .data[[se_col]],
      ci_lower = estimate - z_val * std_error,
      ci_upper = estimate + z_val * std_error,
      term_matches_pattern = stringr::str_detect(term, term_pattern),
      term_value = extract_first_number(term),
      n_treated_units = n_treated_units,
      n_treated_events = n_treated_events,
      n_control_units = n_control_units,
      n_total_units = n_total_units,
      n_rows_model_data = n_rows_model_data,
      subgroup_color = subgroup_color,
      subgroup_run_id = glue::glue("{analysis_id}__{outcome}__{group_id}__{model_id}__{subgroup}"),
      coef_run_id = glue::glue("{analysis_id}__{outcome}__{group_id}__{model_id}__{subgroup}__{vcov_id}")
    )
  
  if (length(t_col) == 1) {
    out$t_value <- ct_tbl[[t_col]]
  } else {
    out$t_value <- NA_real_
  }
  
  if (length(p_col) == 1) {
    out$p_value <- ct_tbl[[p_col]]
  } else {
    out$p_value <- NA_real_
  }
  
  out |>
    dplyr::select(
      analysis_id,
      outcome,
      group_id,
      subgroup_col,
      subgroup,
      model_id,
      vcov_id,
      vcov_label,
      formula_template,
      term,
      estimate,
      std_error,
      t_value,
      p_value,
      ci_lower,
      ci_upper,
      term_matches_pattern,
      term_value,
      n_treated_units,
      n_treated_events,
      n_control_units,
      n_total_units,
      n_rows_model_data,
      subgroup_color,
      subgroup_run_id,
      coef_run_id
    )
}


make_subgroup_run_summary <- function(models,
                                      analysis_id,
                                      outcome,
                                      group_id,
                                      subgroup_col,
                                      model_id,
                                      formula_template,
                                      group_palette = NULL) {
  
  purrr::imap_dfr(
    models,
    \(model, subgroup_nm) {
      subgroup_color <- NA_character_
      if (!is.null(group_palette) && subgroup_nm %in% names(group_palette)) {
        subgroup_color <- unname(group_palette[[subgroup_nm]])
      }
      
      tibble::tibble(
        analysis_id = analysis_id,
        outcome = outcome,
        group_id = group_id,
        subgroup_col = subgroup_col,
        subgroup = subgroup_nm,
        model_id = model_id,
        formula_template = formula_template,
        n_treated_units = attr(model, "n_treated_units"),
        n_treated_events = attr(model, "n_treated_events"),
        n_control_units = attr(model, "n_control_units"),
        n_total_units = attr(model, "n_total_units"),
        n_rows_model_data = attr(model, "n_rows_model_data"),
        subgroup_color = subgroup_color,
        subgroup_run_id = glue::glue("{analysis_id}__{outcome}__{group_id}__{model_id}__{subgroup_nm}")
      )
    }
  )
}


run_one_dydid_spec <- function(parquet_files,
                               analysis_id,
                               outcome,
                               group_id,
                               group_fun,
                               group_args,
                               model_id,
                               formula_template,
                               term_pattern,
                               unit_id,
                               event_id,
                               vcov_specs,
                               dir_out,
                               group_palette = NULL,
                               trt_col = "fire",
                               ci_level = 0.95,
                               mem.clean = TRUE,
                               lean = FALSE,
                               data_filter = NULL) {
  
  if (length(parquet_files) == 0) {
    stop(glue::glue("No parquet files supplied for analysis_id = '{analysis_id}'."))
  }
  
  if (!all(file.exists(parquet_files))) {
    missing_files <- parquet_files[!file.exists(parquet_files)]
    stop(
      glue::glue(
        "These parquet files do not exist:\n{paste(missing_files, collapse = '\n')}"
      )
    )
  }
  
  if (is.null(group_args$cd_nm)) {
    stop("group_args must include cd_nm.")
  }
  
  subgroup_col <- group_args$cd_nm
  event_id <- normalize_optional_colname(event_id)
  
  run_id <- glue::glue("{analysis_id}__{outcome}__{group_id}__{model_id}")
  
  print(glue::glue("Running {run_id}"))
  
  safe_analysis <- safe_path_component(analysis_id)
  safe_outcome <- safe_path_component(outcome)
  safe_group <- safe_path_component(group_id)
  safe_model <- safe_path_component(model_id)
  safe_run_id <- safe_path_component(run_id)
  
  run_dir <- file.path(
    dir_out,
    "tables",
    "by_run",
    paste0("analysis_id=", safe_analysis),
    paste0("outcome=", safe_outcome),
    paste0("group_id=", safe_group),
    paste0("model_id=", safe_model)
  )
  
  dir_ensure_local(c(run_dir))
  
  df <- arrow::open_dataset(parquet_files, format = "parquet") |>
    dplyr::collect()
  
  df <- apply_optional_filter(df, data_filter = data_filter)
  
  df_grouped <- do.call(
    group_fun,
    c(list(df = df), group_args)
  )
  
  if (!subgroup_col %in% names(df_grouped)) {
    stop(glue::glue("Grouping function did not create column '{subgroup_col}'."))
  }
  
  model_formula <- build_model_formula(
    formula_template = formula_template,
    outcome = outcome
  )
  
  needed_cols <- get_needed_columns(
    formula = model_formula,
    trt_col = trt_col,
    subgroup_col = subgroup_col,
    unit_id = unit_id,
    event_id = event_id,
    vcov_vars = vcov_specs$vcov_vars
  )
  
  missing_cols <- setdiff(needed_cols, names(df_grouped))
  if (length(missing_cols) > 0) {
    stop(
      "These required columns are missing after grouping: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  df_grouped <- df_grouped |>
    dplyr::select(dplyr::all_of(needed_cols))
  
  fit_started <- Sys.time()
  
  subgroup_models <- fit_subgroup_models(
    df = df_grouped,
    subgroup_col = subgroup_col,
    trt_col = trt_col,
    formula = model_formula,
    unit_id = unit_id,
    event_id = event_id,
    mem.clean = mem.clean,
    lean = lean
  )
  
  fit_finished <- Sys.time()
  
  subgroup_summary_tbl <- make_subgroup_run_summary(
    models = subgroup_models,
    analysis_id = analysis_id,
    outcome = outcome,
    group_id = group_id,
    subgroup_col = subgroup_col,
    model_id = model_id,
    formula_template = formula_template,
    group_palette = group_palette
  ) |>
    dplyr::mutate(run_id = run_id)
  
  coef_tbl <- purrr::imap_dfr(
    subgroup_models,
    \(model, subgroup_nm) {
      purrr::pmap_dfr(
        list(
          vcov = vcov_specs$vcov,
          vcov_id = vcov_specs$vcov_id,
          vcov_label = vcov_specs$vcov_label
        ),
        \(vcov, vcov_id, vcov_label) {
          extract_coef_table_one_vcov(
            model = model,
            vcov = vcov,
            vcov_id = vcov_id,
            vcov_label = vcov_label,
            analysis_id = analysis_id,
            outcome = outcome,
            group_id = group_id,
            subgroup_col = subgroup_col,
            model_id = model_id,
            formula_template = formula_template,
            term_pattern = term_pattern,
            ci_level = ci_level,
            group_palette = group_palette
          )
        }
      )
    }
  ) |>
    dplyr::mutate(run_id = run_id)
  
  registry_tbl <- tibble::tibble(
    run_id = run_id,
    analysis_id = analysis_id,
    outcome = outcome,
    group_id = group_id,
    model_id = model_id,
    subgroup_col = subgroup_col,
    unit_id = unit_id,
    event_id = event_id,
    trt_col = trt_col,
    formula_template = formula_template,
    formula_resolved = paste(deparse(model_formula), collapse = " "),
    term_pattern = term_pattern,
    n_subgroups = dplyr::n_distinct(subgroup_summary_tbl$subgroup),
    subgroup_levels = paste(sort(unique(subgroup_summary_tbl$subgroup)), collapse = ","),
    fit_started = fit_started,
    fit_finished = fit_finished
  )
  
  print(glue::glue("Saving to: {run_dir}"))
  
  subgroup_summary_file <- file.path(
    run_dir,
    glue::glue("subgroup_run_summary__{safe_run_id}.parquet")
  )
  
  coef_file <- file.path(
    run_dir,
    glue::glue("coef_expanded_vcov__{safe_run_id}.parquet")
  )
  
  registry_file <- file.path(
    run_dir,
    glue::glue("run_registry__{safe_run_id}.parquet")
  )
  
  arrow::write_parquet(subgroup_summary_tbl, subgroup_summary_file)
  arrow::write_parquet(coef_tbl, coef_file)
  
  registry_tbl <- registry_tbl |>
    dplyr::mutate(
      subgroup_run_summary_file = subgroup_summary_file,
      coef_expanded_vcov_file = coef_file,
      run_registry_file = registry_file
    )
  
  arrow::write_parquet(registry_tbl, registry_file)
  
  rm(subgroup_models, df, df_grouped)
  gc()
  
  list(
    subgroup_run_summary = subgroup_summary_tbl,
    coef_expanded_vcov = coef_tbl,
    run_registry = registry_tbl
  )
}


run_dydid_experiment <- function(analysis_specs,
                                 outcome_specs,
                                 group_specs,
                                 model_specs,
                                 vcov_specs,
                                 dir_out,
                                 group_palette = NULL,
                                 trt_col = "fire",
                                 ci_level = 0.95,
                                 mem.clean = TRUE,
                                 lean = FALSE,
                                 data_filter = NULL,
                                 .progress = FALSE) {
  
  required_analysis_cols <- c("analysis_id", "parquet_files")
  required_outcome_cols <- c("outcome")
  required_group_cols <- c("group_id", "group_fun", "group_args")
  required_model_cols <- c("model_id", "formula_template")
  required_vcov_cols <- c("vcov_id", "vcov", "vcov_label")
  
  if (!all(required_analysis_cols %in% names(analysis_specs))) {
    stop("analysis_specs must contain: analysis_id, parquet_files")
  }
  
  if (!all(required_outcome_cols %in% names(outcome_specs))) {
    stop("outcome_specs must contain: outcome")
  }
  
  if (!all(required_group_cols %in% names(group_specs))) {
    stop("group_specs must contain: group_id, group_fun, group_args")
  }
  
  if (!all(required_model_cols %in% names(model_specs))) {
    stop("model_specs must contain: model_id, formula_template")
  }
  
  if (!all(required_vcov_cols %in% names(vcov_specs))) {
    stop("vcov_specs must contain: vcov_id, vcov, vcov_label")
  }
  
  if (!"term_pattern" %in% names(model_specs)) {
    model_specs$term_pattern <- ".*"
  }
  
  if (!"unit_id" %in% names(model_specs)) {
    model_specs$unit_id <- "pt_id"
  }
  
  if (!"event_id" %in% names(model_specs)) {
    model_specs$event_id <- NA_character_
  }
  
  if (!"vcov_vars" %in% names(vcov_specs)) {
    vcov_specs$vcov_vars <- replicate(nrow(vcov_specs), character(0), simplify = FALSE)
  }
  
  dir_ensure_local(c(
    dir_out,
    file.path(dir_out, "tables"),
    file.path(dir_out, "tables", "all"),
    file.path(dir_out, "tables", "by_run"),
    file.path(dir_out, "metadata")
  ))
  
  run_grid <- tidyr::crossing(
    analysis_specs,
    outcome_specs,
    group_specs,
    model_specs
  ) |>
    dplyr::mutate(
      run_id = glue::glue("{analysis_id}__{outcome}__{group_id}__{model_id}")
    )
  
  results <- purrr::pmap(
    list(
      parquet_files = run_grid$parquet_files,
      analysis_id = run_grid$analysis_id,
      outcome = run_grid$outcome,
      group_id = run_grid$group_id,
      group_fun = run_grid$group_fun,
      group_args = run_grid$group_args,
      model_id = run_grid$model_id,
      formula_template = run_grid$formula_template,
      term_pattern = run_grid$term_pattern,
      unit_id = run_grid$unit_id,
      event_id = run_grid$event_id
    ),
    \(parquet_files,
      analysis_id,
      outcome,
      group_id,
      group_fun,
      group_args,
      model_id,
      formula_template,
      term_pattern,
      unit_id,
      event_id) {
      run_one_dydid_spec(
        parquet_files = parquet_files,
        analysis_id = analysis_id,
        outcome = outcome,
        group_id = group_id,
        group_fun = group_fun,
        group_args = group_args,
        model_id = model_id,
        formula_template = formula_template,
        term_pattern = term_pattern,
        unit_id = unit_id,
        event_id = event_id,
        vcov_specs = vcov_specs,
        dir_out = dir_out,
        group_palette = group_palette,
        trt_col = trt_col,
        ci_level = ci_level,
        mem.clean = mem.clean,
        lean = lean,
        data_filter = data_filter
      )
    },
    .progress = .progress
  )
  
  subgroup_summary_all <- purrr::map_dfr(results, "subgroup_run_summary")
  coef_all <- purrr::map_dfr(results, "coef_expanded_vcov")
  registry_all <- purrr::map_dfr(results, "run_registry")
  
  subgroup_summary_parquet <- file.path(
    dir_out, "tables", "all", "subgroup_run_summary.parquet"
  )
  subgroup_summary_csv <- file.path(
    dir_out, "tables", "all", "subgroup_run_summary.csv"
  )
  
  coef_all_parquet <- file.path(
    dir_out, "tables", "all", "coef_expanded_vcov.parquet"
  )
  coef_all_csv <- file.path(
    dir_out, "tables", "all", "coef_expanded_vcov.csv"
  )
  
  registry_all_parquet <- file.path(
    dir_out, "tables", "all", "run_registry.parquet"
  )
  registry_all_csv <- file.path(
    dir_out, "tables", "all", "run_registry.csv"
  )
  
  run_grid_csv <- file.path(
    dir_out, "metadata", "run_grid.csv"
  )
  specs_rds <- file.path(
    dir_out, "metadata", "specs_snapshot.rds"
  )
  experiment_info_rds <- file.path(
    dir_out, "metadata", "dydid_experiment_info.rds"
  )
  
  arrow::write_parquet(subgroup_summary_all, subgroup_summary_parquet)
  arrow::write_parquet(coef_all, coef_all_parquet)
  arrow::write_parquet(registry_all, registry_all_parquet)
  
  readr::write_csv(subgroup_summary_all, subgroup_summary_csv)
  readr::write_csv(coef_all, coef_all_csv)
  readr::write_csv(registry_all, registry_all_csv)
  readr::write_csv(run_grid, run_grid_csv)
  
  saveRDS(
    list(
      analysis_specs = analysis_specs,
      outcome_specs = outcome_specs,
      group_specs = group_specs,
      model_specs = model_specs,
      vcov_specs = vcov_specs,
      group_palette = group_palette
    ),
    specs_rds
  )
  
  experiment_info <- list(
    dir_out = dir_out,
    n_runs = nrow(run_grid),
    files = list(
      subgroup_run_summary_parquet = subgroup_summary_parquet,
      subgroup_run_summary_csv = subgroup_summary_csv,
      coef_expanded_vcov_parquet = coef_all_parquet,
      coef_expanded_vcov_csv = coef_all_csv,
      run_registry_parquet = registry_all_parquet,
      run_registry_csv = registry_all_csv,
      run_grid_csv = run_grid_csv,
      specs_snapshot_rds = specs_rds
    )
  )
  
  saveRDS(experiment_info, experiment_info_rds)
  
  invisible(
    list(
      subgroup_run_summary = subgroup_summary_all,
      coef_expanded_vcov = coef_all,
      run_registry = registry_all,
      run_grid = run_grid,
      experiment_info = experiment_info
    )
  )
}