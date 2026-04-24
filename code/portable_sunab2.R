
# ESTIMATION ----

# Portable large-scale Sun-Abraham pipeline
# Per-run folders are the source of truth.
# Top-level metadata is append-only via timestamped files.

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
  "tidyr",
  "jsonlite",
  "tictoc"
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


serialize_object_json <- function(x) {
  jsonlite::toJSON(
    x,
    auto_unbox = TRUE,
    null = "null",
    digits = NA,
    pretty = FALSE
  )
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


data_filter_to_chr <- function(data_filter = NULL) {
  if (is.null(data_filter)) {
    return(NA_character_)
  }
  
  if (inherits(data_filter, "formula")) {
    return(paste(deparse(data_filter), collapse = " "))
  }
  
  as.character(data_filter)
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
    \(p) {
      dir.create(p, recursive = TRUE, showWarnings = FALSE)
      if (!dir.exists(p)) {
        stop("Failed to create directory: ", p)
      }
    }
  )
}


timestamp_now <- function() {
  format(Sys.time(), "%Y%m%d_%H%M%S")
}


make_run_stub <- function(analysis_id, outcome, group_id, model_id) {
  paste(
    safe_path_component(analysis_id),
    safe_path_component(outcome),
    safe_path_component(group_id),
    safe_path_component(model_id),
    sep = "__"
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
  
  sunab_vars <- parse_sunab_vars(formula)
  
  event_time_support <- make_event_time_support_one_model(
    df_est = df_est,
    trt_col = trt_col,
    unit_id = unit_id,
    event_id = event_id,
    cohort_var = sunab_vars$cohort_var,
    time_var = sunab_vars$time_var
  )
  
  attr(model, "subgroup") <- subgroup_value
  attr(model, "subgroup_col") <- subgroup_col
  attr(model, "unit_id") <- unit_id
  attr(model, "event_id") <- normalize_optional_colname(event_id)
  attr(model, "cohort_var") <- sunab_vars$cohort_var
  attr(model, "time_var") <- sunab_vars$time_var
  attr(model, "event_time_support") <- event_time_support
  
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

extract_vcov_table_one_vcov <- function(model,
                                        vcov,
                                        vcov_id,
                                        vcov_label,
                                        analysis_id,
                                        outcome,
                                        group_id,
                                        subgroup_col,
                                        model_id,
                                        formula_template) {
  
  vc <- stats::vcov(model, vcov = vcov)
  vc_mat <- as.matrix(vc)
  
  if (is.null(rownames(vc_mat)) || is.null(colnames(vc_mat))) {
    stop("vcov matrix must have row and column names.")
  }
  
  subgroup <- attr(model, "subgroup")
  
  expand.grid(
    term_i = rownames(vc_mat),
    term_j = colnames(vc_mat),
    stringsAsFactors = FALSE
  ) |>
    tibble::as_tibble() |>
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
      vcov_value = purrr::map2_dbl(
        term_i,
        term_j,
        \(i, j) vc_mat[i, j]
      ),
      subgroup_run_id = glue::glue(
        "{analysis_id}__{outcome}__{group_id}__{model_id}__{subgroup}"
      ),
      vcov_run_id = glue::glue(
        "{analysis_id}__{outcome}__{group_id}__{model_id}__{subgroup}__{vcov_id}"
      )
    ) |>
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
      term_i,
      term_j,
      vcov_value,
      subgroup_run_id,
      vcov_run_id
    )
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
  tic(glue::glue("VCOV: {vcov_label}"))
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
  event_time_support <- attr(model, "event_time_support")
  
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
  
  if (!is.null(event_time_support) && nrow(event_time_support) > 0) {
    out <- out |>
      dplyr::left_join(
        event_time_support |>
          dplyr::rename(term_value = event_time),
        by = "term_value"
      )
  } else {
    out <- out |>
      dplyr::mutate(
        n_ptids = NA_integer_,
        n_fireids = NA_integer_,
        n_rows_treated = NA_integer_
      )
  }
  
  toc()
  
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
      n_ptids,
      n_fireids,
      n_rows_treated,
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
                               data_filter = NULL,
                               skip_existing = TRUE) {
  
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
  
  run_stub <- make_run_stub(
    analysis_id = analysis_id,
    outcome = outcome,
    group_id = group_id,
    model_id = model_id
  )
  
  run_dir <- file.path(
    dir_out,
    "tables",
    "by_run",
    run_stub
  )
  
  dir_ensure_local(c(run_dir))
  
  coef_file <- file.path(run_dir, "coef.parquet")
  subgroup_file <- file.path(run_dir, "subgroup.parquet")
  support_file <- file.path(run_dir, "support.parquet")
  vcov_file <- file.path(run_dir, "vcov.parquet")
  registry_file <- file.path(run_dir, "registry.parquet")
  run_spec_file <- file.path(run_dir, "run_spec.rds")
  
  if (skip_existing &&
      file.exists(coef_file) &&
      file.exists(subgroup_file) &&
      file.exists(support_file) &&
      file.exists(vcov_file) &&
      file.exists(registry_file) &&
      file.exists(run_spec_file)) {
    
    return(
      list(
        subgroup_run_summary_file = subgroup_file,
        event_time_support_file = support_file,
        coef_expanded_vcov_file = coef_file,
        vcov_file = vcov_file,
        run_registry_file = registry_file,
        run_spec_file = run_spec_file,
        skipped_existing = TRUE
      )
    )
  }
  
  fit_started <- Sys.time()
  
  df <- arrow::open_dataset(parquet_files, format = "parquet") |>
    dplyr::collect()
  
  n_rows_read <- nrow(df)
  df <- apply_optional_filter(df, data_filter = data_filter)
  n_rows_after_filter <- nrow(df)
  
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
  
  subgroup_tbl <- make_subgroup_run_summary(
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
  
  support_tbl <- make_event_time_support_run_summary(
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
  
  vcov_tbl <- purrr::imap_dfr(
    subgroup_models,
    \(model, subgroup_nm) {
      purrr::pmap_dfr(
        list(
          vcov = vcov_specs$vcov,
          vcov_id = vcov_specs$vcov_id,
          vcov_label = vcov_specs$vcov_label
        ),
        \(vcov, vcov_id, vcov_label) {
          extract_vcov_table_one_vcov(
            model = model,
            vcov = vcov,
            vcov_id = vcov_id,
            vcov_label = vcov_label,
            analysis_id = analysis_id,
            outcome = outcome,
            group_id = group_id,
            subgroup_col = subgroup_col,
            model_id = model_id,
            formula_template = formula_template
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
    data_filter = data_filter_to_chr(data_filter),
    input_parquet_files = paste(parquet_files, collapse = " | "),
    input_parquet_n = length(parquet_files),
    group_fun_name = deparse(substitute(group_fun)),
    group_args_json = serialize_object_json(group_args),
    vcov_ids = paste(vcov_specs$vcov_id, collapse = " | "),
    vcov_labels = paste(vcov_specs$vcov_label, collapse = " | "),
    vcov_vars_json = serialize_object_json(vcov_specs$vcov_vars),
    n_rows_read = n_rows_read,
    n_rows_after_filter = n_rows_after_filter,
    n_subgroups = dplyr::n_distinct(subgroup_tbl$subgroup),
    subgroup_levels = paste(sort(unique(subgroup_tbl$subgroup)), collapse = ","),
    n_support_rows = nrow(support_tbl),
    n_vcov_rows = nrow(vcov_tbl),
    fit_started = fit_started,
    fit_finished = fit_finished,
    subgroup_run_summary_file = subgroup_file,
    event_time_support_file = support_file,
    coef_expanded_vcov_file = coef_file,
    vcov_file = vcov_file,
    run_registry_file = registry_file,
    run_spec_file = run_spec_file
  )
  
  run_spec <- list(
    run_id = run_id,
    run_stub = run_stub,
    analysis_id = analysis_id,
    outcome = outcome,
    group_id = group_id,
    model_id = model_id,
    subgroup_col = subgroup_col,
    unit_id = unit_id,
    event_id = event_id,
    trt_col = trt_col,
    formula_template = formula_template,
    formula_resolved = model_formula,
    term_pattern = term_pattern,
    data_filter = data_filter,
    parquet_files = parquet_files,
    group_fun = group_fun,
    group_args = group_args,
    vcov_specs = vcov_specs,
    group_palette = group_palette,
    ci_level = ci_level,
    mem.clean = mem.clean,
    lean = lean
  )
  
  tic("Writing files")
  arrow::write_parquet(subgroup_tbl, subgroup_file)
  arrow::write_parquet(support_tbl, support_file)
  arrow::write_parquet(coef_tbl, coef_file)
  arrow::write_parquet(vcov_tbl, vcov_file)
  arrow::write_parquet(registry_tbl, registry_file)
  saveRDS(run_spec, run_spec_file)
  toc()
  
  rm(
    subgroup_models, df, df_grouped,
    subgroup_tbl, support_tbl, coef_tbl, vcov_tbl, registry_tbl, run_spec
  )
  gc()
  
  list(
    subgroup_run_summary_file = subgroup_file,
    event_time_support_file = support_file,
    coef_expanded_vcov_file = coef_file,
    vcov_file = vcov_file,
    run_registry_file = registry_file,
    run_spec_file = run_spec_file,
    skipped_existing = FALSE
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
                                 skip_existing = TRUE,
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
      run_id = glue::glue("{analysis_id}__{outcome}__{group_id}__{model_id}"),
      run_stub = purrr::pmap_chr(
        list(analysis_id, outcome, group_id, model_id),
        \(analysis_id, outcome, group_id, model_id) {
          make_run_stub(analysis_id, outcome, group_id, model_id)
        }
      )
    )
  
  ts <- timestamp_now()
  
  run_grid_file <- file.path(
    dir_out,
    "metadata",
    glue::glue("run_grid__{ts}.csv")
  )
  
  specs_snapshot_file <- file.path(
    dir_out,
    "metadata",
    glue::glue("specs_snapshot__{ts}.rds")
  )
  
  experiment_info_file <- file.path(
    dir_out,
    "metadata",
    glue::glue("dydid_experiment_info__{ts}.rds")
  )
  
  readr::write_csv(run_grid, run_grid_file)
  
  saveRDS(
    list(
      analysis_specs = analysis_specs,
      outcome_specs = outcome_specs,
      group_specs = group_specs,
      model_specs = model_specs,
      vcov_specs = vcov_specs,
      group_palette = group_palette,
      trt_col = trt_col,
      ci_level = ci_level,
      mem.clean = mem.clean,
      lean = lean,
      data_filter = data_filter,
      skip_existing = skip_existing
    ),
    specs_snapshot_file
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
        data_filter = data_filter,
        skip_existing = skip_existing
      )
    },
    .progress = .progress
  )
  
  experiment_info <- list(
    dir_out = dir_out,
    timestamp = ts,
    n_runs_planned = nrow(run_grid),
    n_runs_returned = length(results),
    metadata_files = list(
      run_grid_csv = run_grid_file,
      specs_snapshot_rds = specs_snapshot_file
    ),
    per_run_mode_only = TRUE,
    tables_all_written = FALSE
  )
  
  saveRDS(experiment_info, experiment_info_file)
  
  invisible(
    list(
      run_grid = run_grid,
      run_results = results,
      metadata_files = list(
        run_grid_csv = run_grid_file,
        specs_snapshot_rds = specs_snapshot_file,
        experiment_info_rds = experiment_info_file
      )
    )
  )
}


parse_sunab_vars <- function(formula) {
  txt <- paste(deparse(formula), collapse = " ")
  txt <- gsub("\\s+", " ", txt)
  
  m <- stringr::str_match(
    txt,
    "sunab\\s*\\(\\s*([^,]+?)\\s*,\\s*([^,\\)]+?)\\s*(?:,|\\))"
  )
  
  if (nrow(m) == 0 || all(is.na(m[1, ]))) {
    return(
      list(
        cohort_var = NA_character_,
        time_var = NA_character_
      )
    )
  }
  
  list(
    cohort_var = trimws(m[1, 2]),
    time_var = trimws(m[1, 3])
  )
}

make_event_time_support_one_model <- function(df_est,
                                              trt_col,
                                              unit_id,
                                              event_id = NA_character_,
                                              cohort_var,
                                              time_var) {
  
  event_id <- normalize_optional_colname(event_id)
  cohort_var <- normalize_optional_colname(cohort_var)
  time_var <- normalize_optional_colname(time_var)
  
  if (is.na(cohort_var) || is.na(time_var)) {
    return(
      tibble::tibble(
        event_time = numeric(),
        n_ptids = integer(),
        n_fireids = integer(),
        n_rows_treated = integer()
      )
    )
  }
  
  needed <- c(trt_col, unit_id, event_id, cohort_var, time_var)
  needed <- stats::na.omit(unique(needed))
  
  missing_cols <- setdiff(needed, names(df_est))
  if (length(missing_cols) > 0) {
    stop(
      "These columns needed for event-time support are missing: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  d_treated <- df_est |>
    dplyr::filter(.data[[trt_col]] == 1) |>
    dplyr::mutate(
      event_time = suppressWarnings(
        as.numeric(.data[[time_var]]) - as.numeric(.data[[cohort_var]])
      )
    ) |>
    dplyr::filter(!is.na(event_time))
  
  if (nrow(d_treated) == 0) {
    return(
      tibble::tibble(
        event_time = numeric(),
        n_ptids = integer(),
        n_fireids = integer(),
        n_rows_treated = integer()
      )
    )
  }
  
  if (is.na(event_id)) {
    d_treated |>
      dplyr::group_by(event_time) |>
      dplyr::summarise(
        n_ptids = dplyr::n_distinct(.data[[unit_id]]),
        n_fireids = NA_integer_,
        n_rows_treated = dplyr::n(),
        .groups = "drop"
      ) |>
      dplyr::arrange(event_time)
  } else {
    d_treated |>
      dplyr::group_by(event_time) |>
      dplyr::summarise(
        n_ptids = dplyr::n_distinct(.data[[unit_id]]),
        n_fireids = dplyr::n_distinct(.data[[event_id]]),
        n_rows_treated = dplyr::n(),
        .groups = "drop"
      ) |>
      dplyr::arrange(event_time)
  }
}

make_event_time_support_run_summary <- function(models,
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
      
      support_tbl <- attr(model, "event_time_support")
      cohort_var <- attr(model, "cohort_var")
      time_var <- attr(model, "time_var")
      
      if (is.null(support_tbl) || nrow(support_tbl) == 0) {
        return(
          tibble::tibble(
            analysis_id = analysis_id,
            outcome = outcome,
            group_id = group_id,
            subgroup_col = subgroup_col,
            subgroup = subgroup_nm,
            model_id = model_id,
            formula_template = formula_template,
            cohort_var = normalize_optional_colname(cohort_var),
            time_var = normalize_optional_colname(time_var),
            event_time = numeric(),
            n_ptids = integer(),
            n_fireids = integer(),
            n_rows_treated = integer(),
            subgroup_color = character(),
            subgroup_run_id = character()
          )
        )
      }
      
      support_tbl |>
        dplyr::mutate(
          analysis_id = analysis_id,
          outcome = outcome,
          group_id = group_id,
          subgroup_col = subgroup_col,
          subgroup = subgroup_nm,
          model_id = model_id,
          formula_template = formula_template,
          cohort_var = normalize_optional_colname(cohort_var),
          time_var = normalize_optional_colname(time_var),
          subgroup_color = subgroup_color,
          subgroup_run_id = glue::glue(
            "{analysis_id}__{outcome}__{group_id}__{model_id}__{subgroup_nm}"
          )
        ) |>
        dplyr::select(
          analysis_id,
          outcome,
          group_id,
          subgroup_col,
          subgroup,
          model_id,
          formula_template,
          cohort_var,
          time_var,
          event_time,
          n_ptids,
          n_fireids,
          n_rows_treated,
          subgroup_color,
          subgroup_run_id
        )
    }
  )
}


rebuild_dydid_tables_from_by_run <- function(dir_out,
                                             write_csv = TRUE,
                                             recursive = TRUE) {
  
  dir_by_run <- file.path(dir_out, "tables", "by_run")
  dir_all <- file.path(dir_out, "tables", "all")
  
  if (!dir.exists(dir_by_run)) {
    stop("by_run directory does not exist: ", dir_by_run)
  }
  
  dir.create(dir_all, recursive = TRUE, showWarnings = FALSE)
  
  coef_files <- list.files(
    dir_by_run,
    pattern = "^coef\\.parquet$",
    recursive = recursive,
    full.names = TRUE
  )
  
  subgroup_files <- list.files(
    dir_by_run,
    pattern = "^subgroup\\.parquet$",
    recursive = recursive,
    full.names = TRUE
  )
  
  support_files <- list.files(
    dir_by_run,
    pattern = "^support\\.parquet$",
    recursive = recursive,
    full.names = TRUE
  )
  
  registry_files <- list.files(
    dir_by_run,
    pattern = "^registry\\.parquet$",
    recursive = recursive,
    full.names = TRUE
  )
  
  run_spec_files <- list.files(
    dir_by_run,
    pattern = "^run_spec\\.rds$",
    recursive = recursive,
    full.names = TRUE
  )
  
  if (length(coef_files) == 0) {
    stop("No coef.parquet files found under: ", dir_by_run)
  }
  
  if (length(subgroup_files) == 0) {
    stop("No subgroup.parquet files found under: ", dir_by_run)
  }
  
  if (length(support_files) == 0) {
    stop("No support.parquet files found under: ", dir_by_run)
  }
  
  if (length(registry_files) == 0) {
    stop("No registry.parquet files found under: ", dir_by_run)
  }
  
  coef_tbl <- purrr::map_dfr(
    coef_files,
    arrow::read_parquet
  )
  
  subgroup_tbl <- purrr::map_dfr(
    subgroup_files,
    arrow::read_parquet
  )
  
  support_tbl <- purrr::map_dfr(
    support_files,
    arrow::read_parquet
  )
  
  registry_tbl <- purrr::map_dfr(
    registry_files,
    arrow::read_parquet
  )
  
  coef_parquet <- file.path(dir_all, "coef_expanded_vcov.parquet")
  subgroup_parquet <- file.path(dir_all, "subgroup_run_summary.parquet")
  support_parquet <- file.path(dir_all, "event_time_support.parquet")
  registry_parquet <- file.path(dir_all, "run_registry.parquet")
  
  arrow::write_parquet(coef_tbl, coef_parquet)
  arrow::write_parquet(subgroup_tbl, subgroup_parquet)
  arrow::write_parquet(support_tbl, support_parquet)
  arrow::write_parquet(registry_tbl, registry_parquet)
  
  coef_csv <- NULL
  subgroup_csv <- NULL
  support_csv <- NULL
  registry_csv <- NULL
  
  if (write_csv) {
    coef_csv <- file.path(dir_all, "coef_expanded_vcov.csv")
    subgroup_csv <- file.path(dir_all, "subgroup_run_summary.csv")
    support_csv <- file.path(dir_all, "event_time_support.csv")
    registry_csv <- file.path(dir_all, "run_registry.csv")
    
    readr::write_csv(coef_tbl, coef_csv)
    readr::write_csv(subgroup_tbl, subgroup_csv)
    readr::write_csv(support_tbl, support_csv)
    readr::write_csv(registry_tbl, registry_csv)
  }
  
  spec_index_tbl <- NULL
  spec_index_parquet <- NULL
  spec_index_csv <- NULL
  
  if (length(run_spec_files) > 0) {
    spec_index_tbl <- tibble::tibble(
      run_spec_file = run_spec_files
    ) |>
      dplyr::mutate(
        run_spec = purrr::map(run_spec_file, readRDS),
        run_id = purrr::map_chr(run_spec, "run_id"),
        analysis_id = purrr::map_chr(run_spec, "analysis_id"),
        outcome = purrr::map_chr(run_spec, "outcome"),
        group_id = purrr::map_chr(run_spec, "group_id"),
        model_id = purrr::map_chr(run_spec, "model_id"),
        subgroup_col = purrr::map_chr(run_spec, "subgroup_col"),
        unit_id = purrr::map_chr(run_spec, "unit_id"),
        event_id = purrr::map_chr(
          run_spec,
          \(x) {
            if (is.null(x$event_id) || length(x$event_id) == 0 || is.na(x$event_id)) {
              return(NA_character_)
            }
            as.character(x$event_id)
          }
        ),
        trt_col = purrr::map_chr(run_spec, "trt_col"),
        formula_template = purrr::map_chr(run_spec, "formula_template"),
        term_pattern = purrr::map_chr(run_spec, "term_pattern"),
        data_filter = purrr::map_chr(
          run_spec,
          \(x) {
            if (is.null(x$data_filter)) {
              return(NA_character_)
            }
            paste(deparse(x$data_filter), collapse = " ")
          }
        ),
        parquet_files = purrr::map(
          run_spec,
          "parquet_files"
        ),
        parquet_files_chr = purrr::map_chr(
          parquet_files,
          \(x) paste(x, collapse = " | ")
        ),
        group_args_json = purrr::map_chr(
          run_spec,
          \(x) jsonlite::toJSON(
            x$group_args,
            auto_unbox = TRUE,
            null = "null",
            digits = NA,
            pretty = FALSE
          )
        ),
        vcov_ids = purrr::map_chr(
          run_spec,
          \(x) paste(x$vcov_specs$vcov_id, collapse = " | ")
        )
      ) |>
      dplyr::select(
        -run_spec,
        -parquet_files
      )
    
    spec_index_parquet <- file.path(dir_all, "run_spec_index.parquet")
    arrow::write_parquet(spec_index_tbl, spec_index_parquet)
    
    if (write_csv) {
      spec_index_csv <- file.path(dir_all, "run_spec_index.csv")
      readr::write_csv(spec_index_tbl, spec_index_csv)
    }
  }
  
  invisible(
    list(
      coef_expanded_vcov = coef_tbl,
      subgroup_run_summary = subgroup_tbl,
      event_time_support = support_tbl,
      run_registry = registry_tbl,
      run_spec_index = spec_index_tbl,
      files = list(
        coef_expanded_vcov_parquet = coef_parquet,
        subgroup_run_summary_parquet = subgroup_parquet,
        event_time_support_parquet = support_parquet,
        run_registry_parquet = registry_parquet,
        coef_expanded_vcov_csv = coef_csv,
        subgroup_run_summary_csv = subgroup_csv,
        event_time_support_csv = support_csv,
        run_registry_csv = registry_csv,
        run_spec_index_parquet = spec_index_parquet,
        run_spec_index_csv = spec_index_csv
      )
    )
  )
}






# DESCRIPTIVE ----

make_event_time_descriptive_one_subgroup_simple <- function(df,
                                                            subgroup_col,
                                                            subgroup_value,
                                                            trt_col,
                                                            outcome,
                                                            unit_id,
                                                            event_id = NA_character_,
                                                            time_var = "year",
                                                            treated_year_var = "burn_year",
                                                            control_year_var = "mock_burn_year") {
  
  event_id <- normalize_optional_colname(event_id)
  
  needed_cols <- c(
    subgroup_col,
    trt_col,
    outcome,
    unit_id,
    event_id,
    time_var,
    treated_year_var,
    control_year_var
  ) |>
    unique() |>
    stats::na.omit() |>
    as.character()
  
  missing_cols <- setdiff(needed_cols, names(df))
  if (length(missing_cols) > 0) {
    stop(
      "These required columns are missing for descriptive summaries: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  d_treated <- df |>
    dplyr::filter(
      .data[[trt_col]] == 1,
      .data[[subgroup_col]] == subgroup_value
    ) |>
    dplyr::mutate(
      event_time = as.numeric(.data[[time_var]]) - as.numeric(.data[[treated_year_var]]),
      series = "treated"
    ) |>
    dplyr::filter(!is.na(event_time))
  
  d_control <- df |>
    dplyr::filter(.data[[trt_col]] == 0) |>
    dplyr::mutate(
      event_time = as.numeric(.data[[time_var]]) - as.numeric(.data[[control_year_var]]),
      series = "control_mock"
    ) |>
    dplyr::filter(!is.na(event_time))
  
  d_plot <- dplyr::bind_rows(d_treated, d_control)
  
  if (nrow(d_plot) == 0) {
    return(
      tibble::tibble(
        subgroup = character(),
        event_time = numeric(),
        series = character(),
        mean_outcome = numeric(),
        sd_outcome = numeric(),
        n_rows = integer(),
        n_ptids = integer(),
        n_fireids = integer()
      )
    )
  }
  
  if (is.na(event_id)) {
    d_plot |>
      dplyr::group_by(series, event_time) |>
      dplyr::summarise(
        mean_outcome = mean(.data[[outcome]], na.rm = TRUE),
        sd_outcome = stats::sd(.data[[outcome]], na.rm = TRUE),
        n_rows = dplyr::n(),
        n_ptids = dplyr::n_distinct(.data[[unit_id]]),
        n_fireids = NA_integer_,
        .groups = "drop"
      ) |>
      dplyr::mutate(subgroup = subgroup_value) |>
      dplyr::select(
        subgroup,
        event_time,
        series,
        mean_outcome,
        sd_outcome,
        n_rows,
        n_ptids,
        n_fireids
      ) |>
      dplyr::arrange(subgroup, series, event_time)
  } else {
    d_plot |>
      dplyr::group_by(series, event_time) |>
      dplyr::summarise(
        mean_outcome = mean(.data[[outcome]], na.rm = TRUE),
        sd_outcome = stats::sd(.data[[outcome]], na.rm = TRUE),
        n_rows = dplyr::n(),
        n_ptids = dplyr::n_distinct(.data[[unit_id]]),
        n_fireids = dplyr::n_distinct(.data[[event_id]]),
        .groups = "drop"
      ) |>
      dplyr::mutate(subgroup = subgroup_value) |>
      dplyr::select(
        subgroup,
        event_time,
        series,
        mean_outcome,
        sd_outcome,
        n_rows,
        n_ptids,
        n_fireids
      ) |>
      dplyr::arrange(subgroup, series, event_time)
  }
}


make_event_time_descriptive_run_summary_simple <- function(df,
                                                           analysis_id,
                                                           outcome,
                                                           group_id,
                                                           subgroup_col,
                                                           model_id,
                                                           formula_template,
                                                           trt_col,
                                                           unit_id,
                                                           event_id = NA_character_,
                                                           treated_year_var = "burn_year",
                                                           control_year_var = "mock_burn_year",
                                                           time_var = "year",
                                                           group_palette = NULL) {
  
  subgroup_levels <- df |>
    dplyr::filter(
      .data[[trt_col]] == 1,
      !is.na(.data[[subgroup_col]])
    ) |>
    dplyr::distinct(.data[[subgroup_col]]) |>
    dplyr::pull(.data[[subgroup_col]]) |>
    as.character() |>
    sort()
  
  purrr::map_dfr(
    subgroup_levels,
    \(subgroup_value) {
      subgroup_color <- NA_character_
      if (!is.null(group_palette) && subgroup_value %in% names(group_palette)) {
        subgroup_color <- unname(group_palette[[subgroup_value]])
      }
      
      make_event_time_descriptive_one_subgroup_simple(
        df = df,
        subgroup_col = subgroup_col,
        subgroup_value = subgroup_value,
        trt_col = trt_col,
        outcome = outcome,
        unit_id = unit_id,
        event_id = event_id,
        time_var = time_var,
        treated_year_var = treated_year_var,
        control_year_var = control_year_var
      ) |>
        dplyr::mutate(
          analysis_id = analysis_id,
          outcome = outcome,
          group_id = group_id,
          subgroup_col = subgroup_col,
          model_id = model_id,
          formula_template = formula_template,
          treated_year_var = treated_year_var,
          control_year_var = control_year_var,
          time_var = time_var,
          subgroup_color = subgroup_color,
          subgroup_run_id = glue::glue(
            "{analysis_id}__{outcome}__{group_id}__{model_id}__{subgroup_value}"
          ),
          descriptive_run_id = glue::glue(
            "{analysis_id}__{outcome}__{group_id}__{model_id}__{subgroup_value}__{series}"
          )
        ) |>
        dplyr::select(
          analysis_id,
          outcome,
          group_id,
          subgroup_col,
          subgroup,
          model_id,
          formula_template,
          treated_year_var,
          control_year_var,
          time_var,
          event_time,
          series,
          mean_outcome,
          sd_outcome,
          n_rows,
          n_ptids,
          n_fireids,
          subgroup_color,
          subgroup_run_id,
          descriptive_run_id
        )
    }
  )
}


run_one_descriptive_spec_simple <- function(parquet_files,
                                            analysis_id,
                                            outcome,
                                            group_id,
                                            group_fun,
                                            group_args,
                                            model_id,
                                            formula_template,
                                            unit_id,
                                            event_id,
                                            dir_out,
                                            group_palette = NULL,
                                            trt_col = "fire",
                                            data_filter = NULL,
                                            treated_year_var = "burn_year",
                                            control_year_var = "mock_burn_year",
                                            time_var = "year",
                                            skip_existing = TRUE) {
  
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
  run_stub <- make_run_stub(
    analysis_id = analysis_id,
    outcome = outcome,
    group_id = group_id,
    model_id = model_id
  )
  
  run_dir <- file.path(
    dir_out,
    "descriptive",
    "by_run",
    run_stub
  )
  
  dir_ensure_local(run_dir)
  
  trajectory_file <- file.path(run_dir, "event_time_trajectory.parquet")
  registry_file <- file.path(run_dir, "registry.parquet")
  run_spec_file <- file.path(run_dir, "descriptive_spec.rds")
  
  if (skip_existing &&
      file.exists(trajectory_file) &&
      file.exists(registry_file) &&
      file.exists(run_spec_file)) {
    
    return(
      list(
        event_time_trajectory_file = trajectory_file,
        run_registry_file = registry_file,
        run_spec_file = run_spec_file,
        skipped_existing = TRUE
      )
    )
  }
  
  run_started <- Sys.time()
  
  df <- arrow::open_dataset(parquet_files, format = "parquet") |>
    dplyr::collect()
  
  n_rows_read <- nrow(df)
  
  df <- apply_optional_filter(df, data_filter = data_filter)
  n_rows_after_filter <- nrow(df)
  
  df_grouped <- do.call(
    group_fun,
    c(list(df = df), group_args)
  )
  
  if (!subgroup_col %in% names(df_grouped)) {
    stop(glue::glue("Grouping function did not create column '{subgroup_col}'."))
  }
  
  trajectory_tbl <- make_event_time_descriptive_run_summary_simple(
    df = df_grouped,
    analysis_id = analysis_id,
    outcome = outcome,
    group_id = group_id,
    subgroup_col = subgroup_col,
    model_id = model_id,
    formula_template = formula_template,
    trt_col = trt_col,
    unit_id = unit_id,
    event_id = event_id,
    treated_year_var = treated_year_var,
    control_year_var = control_year_var,
    time_var = time_var,
    group_palette = group_palette
  ) |>
    dplyr::mutate(run_id = run_id)
  
  run_finished <- Sys.time()
  
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
    data_filter = data_filter_to_chr(data_filter),
    input_parquet_files = paste(parquet_files, collapse = " | "),
    input_parquet_n = length(parquet_files),
    group_fun_name = deparse(substitute(group_fun)),
    group_args_json = serialize_object_json(group_args),
    treated_year_var = treated_year_var,
    control_year_var = control_year_var,
    time_var = time_var,
    n_rows_read = n_rows_read,
    n_rows_after_filter = n_rows_after_filter,
    n_subgroups = dplyr::n_distinct(trajectory_tbl$subgroup),
    subgroup_levels = paste(sort(unique(trajectory_tbl$subgroup)), collapse = ","),
    run_started = run_started,
    run_finished = run_finished,
    event_time_trajectory_file = trajectory_file,
    run_registry_file = registry_file,
    run_spec_file = run_spec_file
  )
  
  run_spec <- list(
    run_id = run_id,
    run_stub = run_stub,
    analysis_id = analysis_id,
    outcome = outcome,
    group_id = group_id,
    model_id = model_id,
    subgroup_col = subgroup_col,
    unit_id = unit_id,
    event_id = event_id,
    trt_col = trt_col,
    formula_template = formula_template,
    data_filter = data_filter,
    parquet_files = parquet_files,
    group_fun = group_fun,
    group_args = group_args,
    group_palette = group_palette,
    treated_year_var = treated_year_var,
    control_year_var = control_year_var,
    time_var = time_var
  )
  
  arrow::write_parquet(trajectory_tbl, trajectory_file)
  arrow::write_parquet(registry_tbl, registry_file)
  saveRDS(run_spec, run_spec_file)
  
  rm(df, df_grouped, trajectory_tbl, registry_tbl, run_spec)
  gc()
  
  list(
    event_time_trajectory_file = trajectory_file,
    run_registry_file = registry_file,
    run_spec_file = run_spec_file,
    skipped_existing = FALSE
  )
}

run_descriptive_experiment_simple <- function(analysis_specs,
                                              outcome_specs,
                                              group_specs,
                                              model_specs,
                                              dir_out,
                                              group_palette = NULL,
                                              trt_col = "fire",
                                              data_filter = NULL,
                                              treated_year_var = "burn_year",
                                              control_year_var = "mock_burn_year",
                                              time_var = "year",
                                              skip_existing = TRUE,
                                              .progress = FALSE) {
  
  required_analysis_cols <- c("analysis_id", "parquet_files")
  required_outcome_cols <- c("outcome")
  required_group_cols <- c("group_id", "group_fun", "group_args")
  required_model_cols <- c("model_id", "formula_template")
  
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
  
  if (!"unit_id" %in% names(model_specs)) {
    model_specs$unit_id <- "pt_id"
  }
  
  if (!"event_id" %in% names(model_specs)) {
    model_specs$event_id <- NA_character_
  }
  
  dir_ensure_local(c(
    dir_out,
    file.path(dir_out, "descriptive"),
    file.path(dir_out, "descriptive", "by_run"),
    file.path(dir_out, "descriptive", "metadata")
  ))
  
  run_grid <- tidyr::crossing(
    analysis_specs,
    outcome_specs,
    group_specs,
    model_specs
  ) |>
    dplyr::mutate(
      run_id = glue::glue("{analysis_id}__{outcome}__{group_id}__{model_id}"),
      run_stub = purrr::pmap_chr(
        list(analysis_id, outcome, group_id, model_id),
        \(analysis_id, outcome, group_id, model_id) {
          make_run_stub(analysis_id, outcome, group_id, model_id)
        }
      )
    )
  
  ts <- timestamp_now()
  
  run_grid_file <- file.path(
    dir_out,
    "descriptive",
    "metadata",
    glue::glue("run_grid__{ts}.rds")
  )
  
  experiment_info_file <- file.path(
    dir_out,
    "descriptive",
    "metadata",
    glue::glue("descriptive_experiment_info__{ts}.rds")
  )
  
  saveRDS(run_grid, run_grid_file)
  
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
      unit_id,
      event_id) {
      run_one_descriptive_spec_simple(
        parquet_files = parquet_files,
        analysis_id = analysis_id,
        outcome = outcome,
        group_id = group_id,
        group_fun = group_fun,
        group_args = group_args,
        model_id = model_id,
        formula_template = formula_template,
        unit_id = unit_id,
        event_id = event_id,
        dir_out = dir_out,
        group_palette = group_palette,
        trt_col = trt_col,
        data_filter = data_filter,
        treated_year_var = treated_year_var,
        control_year_var = control_year_var,
        time_var = time_var,
        skip_existing = skip_existing
      )
    },
    .progress = .progress
  )
  
  experiment_info <- list(
    dir_out = dir_out,
    timestamp = ts,
    n_runs_planned = nrow(run_grid),
    n_runs_returned = length(results),
    metadata_files = list(
      run_grid_rds = run_grid_file
    )
  )
  
  saveRDS(experiment_info, experiment_info_file)
  
  invisible(
    list(
      run_grid = run_grid,
      run_results = results,
      metadata_files = list(
        run_grid_rds = run_grid_file,
        experiment_info_rds = experiment_info_file
      )
    )
  )
}

rebuild_descriptive_tables_from_by_run_simple <- function(dir_out,
                                                          write_csv = TRUE,
                                                          recursive = TRUE) {
  
  dir_by_run <- file.path(dir_out, "descriptive", "by_run")
  dir_all <- file.path(dir_out, "descriptive", "all")
  
  if (!dir.exists(dir_by_run)) {
    stop("descriptive by_run directory does not exist: ", dir_by_run)
  }
  
  dir.create(dir_all, recursive = TRUE, showWarnings = FALSE)
  
  trajectory_files <- list.files(
    dir_by_run,
    pattern = "^event_time_trajectory\\.parquet$",
    recursive = recursive,
    full.names = TRUE
  )
  
  registry_files <- list.files(
    dir_by_run,
    pattern = "^registry\\.parquet$",
    recursive = recursive,
    full.names = TRUE
  )
  
  if (length(trajectory_files) == 0) {
    stop("No event_time_trajectory.parquet files found under: ", dir_by_run)
  }
  
  if (length(registry_files) == 0) {
    stop("No registry.parquet files found under: ", dir_by_run)
  }
  
  trajectory_tbl <- purrr::map_dfr(
    trajectory_files,
    arrow::read_parquet
  )
  
  registry_tbl <- purrr::map_dfr(
    registry_files,
    arrow::read_parquet
  )
  
  trajectory_parquet <- file.path(dir_all, "event_time_trajectory.parquet")
  registry_parquet <- file.path(dir_all, "run_registry.parquet")
  
  arrow::write_parquet(trajectory_tbl, trajectory_parquet)
  arrow::write_parquet(registry_tbl, registry_parquet)
  
  trajectory_csv <- NULL
  registry_csv <- NULL
  
  if (write_csv) {
    trajectory_csv <- file.path(dir_all, "event_time_trajectory.csv")
    registry_csv <- file.path(dir_all, "run_registry.csv")
    
    readr::write_csv(trajectory_tbl, trajectory_csv)
    readr::write_csv(registry_tbl, registry_csv)
  }
  
  invisible(
    list(
      event_time_trajectory = trajectory_tbl,
      run_registry = registry_tbl,
      files = list(
        event_time_trajectory_parquet = trajectory_parquet,
        run_registry_parquet = registry_parquet,
        event_time_trajectory_csv = trajectory_csv,
        run_registry_csv = registry_csv
      )
    )
  )
}