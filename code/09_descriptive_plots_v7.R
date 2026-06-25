# 08b_descriptive_plots_v7.R
# Descriptive event-time trajectory plots
# -------------------------------------------------------
# Reads from:
#   descriptive/by_run/{run_stub}/event_time_trajectory.parquet
#   (merged into descriptive/all/ by rebuild_descriptive_tables())
#
# Outcomes plotted here:
#   mean_outcome, normalized to the mean of a pre-period baseline window,
#   for the treated and control_mock series overlaid by dummy_group.
# -------------------------------------------------------

# ══════════════════════════════════════════════════════════════════════════════
# Config ----
# ══════════════════════════════════════════════════════════════════════════════

rm(list = ls())

library(here)
here::i_am("code/09_descriptive_plots_v7.R")

library(dplyr); library(ggplot2); library(readr); library(purrr)
library(tibble); library(stringr); library(arrow); library(glue)

source(here::here("code", "weight_dydid_pipeline_v7.R"))


run_name    <- "GEE_resilience_v6_operational_ss500_ts50000"
version     <- "v7"
dir_results <- here::here("results-exo", version)
dir_figs    <- here::here("figs", version, "exo", "descriptive")

dir.create(dir_figs, recursive = TRUE, showWarnings = FALSE)

# Pre-period window used to compute the normalization baseline.
# Excludes ref.p (-6) to match the estimation reference period — the mean
# across this window becomes 0 after normalization.
baseline_window <- -15:-7


rebuild_descriptive_tables <- function(dir_out, write_csv = TRUE, recursive = TRUE) {
  dir_by_run <- file.path(dir_out, "descriptive", "by_run")
  dir_all    <- file.path(dir_out, "descriptive", "all")
  if (!dir.exists(dir_by_run)) stop("descriptive by_run directory does not exist: ", dir_by_run)
  dir.create(dir_all, recursive = TRUE, showWarnings = FALSE)
  
  find_pq <- function(pattern) {
    list.files(dir_by_run, pattern = pattern, recursive = recursive, full.names = TRUE)
  }
  
  traj_files     <- find_pq("^event_time_trajectory\\.parquet$")
  registry_files <- find_pq("^registry\\.parquet$")
  
  if (length(traj_files) == 0) {
    message("No event_time_trajectory.parquet files found under: ", dir_by_run)
  }
  
  # mirrors the read_and_dedup pattern from rebuild_estimation_tables();
  # keeps the most recently written file when the same run appears twice
  read_and_dedup <- function(files, id_cols) {
    if (length(files) == 0) return(NULL)
    purrr::map(files, \(f) {
      x <- arrow::read_parquet(f)
      x[] <- lapply(x, \(col) if (is.list(col)) as.character(col) else col)
      x$.mtime <- file.mtime(f)
      x
    }) |>
      dplyr::bind_rows() |>
      dplyr::group_by(dplyr::across(dplyr::all_of(id_cols))) |>
      dplyr::filter(.mtime == max(.mtime)) |>
      dplyr::ungroup() |>
      dplyr::select(-.mtime)
  }
  
  write_pair <- function(tbl, stem) {
    if (is.null(tbl) || nrow(tbl) == 0) return(list(parquet = NULL, csv = NULL))
    pq  <- file.path(dir_all, paste0(stem, ".parquet"))
    csv <- if (write_csv) file.path(dir_all, paste0(stem, ".csv")) else NULL
    arrow::write_parquet(tbl, pq)
    if (write_csv) readr::write_csv(tbl, csv)
    list(parquet = pq, csv = csv)
  }
  
  # trajectory dedup key: one row per run x dummy_group x series x event_time
  traj_tbl     <- read_and_dedup(traj_files,     c("run_id", "dummy_group", "series", "event_time"))
  registry_tbl <- read_and_dedup(registry_files, c("run_id"))
  
  files <- list(
    event_time_trajectory = write_pair(traj_tbl,     "event_time_trajectory"),
    run_registry          = write_pair(registry_tbl, "descriptive_run_registry")
  )
  
  if (!is.null(traj_tbl)) {
    message(glue::glue(
      "Descriptive tables merged: {nrow(traj_tbl)} trajectory rows, ",
      "{dplyr::n_distinct(traj_tbl$run_id)} runs"
    ))
  }
  
  invisible(list(
    traj_tbl     = traj_tbl,
    registry_tbl = registry_tbl,
    files        = files
  ))
}
all_descriptive_tables <- rebuild_descriptive_tables(dir_out = dir_results, write_csv = TRUE)


# ══════════════════════════════════════════════════════════════════════════════
# Load merged descriptive table ----
# ══════════════════════════════════════════════════════════════════════════════

# rebuild_descriptive_tables() merges all per-run event_time_trajectory.parquets
# into one file; the call above does this, so by the time we get here the file
# is guaranteed to exist if any runs completed successfully.
traj_file <- file.path(dir_results, "descriptive", "all", "event_time_trajectory.parquet")

if (!file.exists(traj_file)) {
  stop(
    "Merged trajectory file not found: ", traj_file, "\n",
    "Run rebuild_descriptive_tables() in script 06 first."
  )
}

traj_all <- arrow::read_parquet(traj_file)


# ══════════════════════════════════════════════════════════════════════════════
# Normalization helper ----
# ══════════════════════════════════════════════════════════════════════════════

# Normalize mean_outcome within each series x dummy_group x run combination
# to the mean of mean_outcome over baseline_window. The normalization is
# applied separately to treated and control_mock so both series start at 0
# in the baseline period, making pre-trend visual comparison straightforward.
#
# SE is unaffected by the constant baseline subtraction (sd / sqrt(n) is
# computed from the raw sd_outcome and n_rows columns written by the pipeline).
normalize_to_baseline <- function(traj,
                                  baseline_window,
                                  outcome_col = "mean_outcome") {
  
  baseline_means <- traj |>
    dplyr::filter(event_time %in% baseline_window) |>
    dplyr::group_by(run_id, subset_id, outcome, group_id, model_id,
                    dummy_group, series) |>
    dplyr::summarise(
      baseline_mean = mean(.data[[outcome_col]], na.rm = TRUE),
      .groups = "drop"
    )
  
  traj |>
    dplyr::left_join(
      baseline_means,
      by = c("run_id", "subset_id", "outcome", "group_id",
             "model_id", "dummy_group", "series")
    ) |>
    dplyr::mutate(
      mean_outcome_norm = .data[[outcome_col]] - baseline_mean,
      se_outcome        = sd_outcome / sqrt(n_rows)
    )
}


# ══════════════════════════════════════════════════════════════════════════════
# Plot helper ----
# ══════════════════════════════════════════════════════════════════════════════

dummy_group_labels <- c(
  "cd_f"   = "Fire only",
  "cd_bf"  = "Biotic stress + fire",
  "cd_df"  = "Drought + fire",
  "cd_bdf" = "Drought + biotic + fire"
)

dummy_group_palette <- c(
  "cd_f"   = "maroon3",
  "cd_bf"  = "forestgreen",
  "cd_df"  = "goldenrod2",
  "cd_bdf" = "purple2"
)

# Line types distinguish treated from control_mock series. Using a single
# color per dummy_group keeps the two series visually linked.
series_linetypes <- c(
  "treated"      = "solid",
  "control_mock" = "dashed"
)

series_labels <- c(
  "treated"      = "Treated",
  "control_mock" = "Control (mock)"
)


#' Plot normalized descriptive event-time trajectories
#'
#' Overlays treated and control_mock series, normalized to a pre-period
#' baseline, for one combination of subset x outcome x group x model.
#' SE ribbons are drawn as +/- se_outcome (i.e. +/- 1 SE, not a CI),
#' computed upstream in normalize_to_baseline() as sd_outcome / sqrt(n_rows).
#'
#' @param traj              Merged trajectory tibble (already normalized).
#' @param subset_id_filter  Character. subset_id to display.
#' @param outcome_filter    Character. outcome to display.
#' @param group_id_filter   Character. group_id to display.
#' @param model_id_filter   Character. model_id to display.
#' @param event_time_range  Integer vector of length 2.
#' @param baseline_window   Integer vector of event times used for baseline
#'                          shading (normalization is done upstream).
#' @param se_ribbon         Logical. If TRUE, draw +/- 1 SE ribbons.
#' @param ribbon_alpha      Alpha for SE ribbons. Lower = more transparent.
#' @param min_n_points      Minimum n_ptids to include a point. NULL = no filter.
#' @param min_n_events      Minimum n_fireids to include a point. NULL = no filter.
#' @param palette           Named color vector keyed on dummy_group.
#' @param group_labels      Named label vector keyed on dummy_group.
#' @param title             Plot title. NULL = auto-generated from filters.
#' @param facet_by_dummy    Logical. If TRUE, facet by dummy_group.
#' @param free_y            Logical. If TRUE, facets use free y scales.
plot_descriptive_trajectory <- function(traj,
                                        subset_id_filter,
                                        outcome_filter,
                                        group_id_filter,
                                        model_id_filter,
                                        event_time_range  = c(-15, 20),
                                        baseline_window   = -15:-7,
                                        se_ribbon         = TRUE,
                                        ribbon_alpha      = 0.12,
                                        min_n_points      = NULL,
                                        min_n_events      = NULL,
                                        palette           = dummy_group_palette,
                                        group_labels      = dummy_group_labels,
                                        title             = NULL,
                                        facet_by_dummy    = TRUE,
                                        free_y            = FALSE) {
  
  # ── Filter ──────────────────────────────────────────────────────────────────
  d <- traj |>
    dplyr::filter(
      subset_id  == subset_id_filter,
      outcome    == outcome_filter,
      group_id   == group_id_filter,
      model_id   == model_id_filter,
      event_time >= event_time_range[1],
      event_time <= event_time_range[2]
    )
  
  if (nrow(d) == 0) {
    warning("No rows after filtering. Check subset/outcome/group/model filters.")
    return(NULL)
  }
  
  if (!is.null(min_n_points)) d <- dplyr::filter(d, !is.na(n_ptids),   n_ptids   >= min_n_points)
  if (!is.null(min_n_events)) d <- dplyr::filter(d, !is.na(n_fireids), n_fireids >= min_n_events)
  
  if (nrow(d) == 0) {
    warning("No rows remain after support filtering.")
    return(NULL)
  }
  
  # ── Labels ──────────────────────────────────────────────────────────────────
  d <- d |>
    dplyr::mutate(
      dummy_group_label = factor(
        dplyr::recode(dummy_group, !!!group_labels),
        levels = unname(group_labels)
      ),
      series_label = factor(
        dplyr::recode(series, !!!series_labels),
        levels = unname(series_labels)
      )
    )
  
  if (is.null(title)) {
    title <- glue::glue(
      "{subset_id_filter} | {outcome_filter} | {group_id_filter} | {model_id_filter}"
    )
  }
  
  # ── Build plot ───────────────────────────────────────────────────────────────
  p <- ggplot2::ggplot(
    d,
    ggplot2::aes(
      x        = event_time,
      y        = mean_outcome_norm,
      color    = dummy_group,
      fill     = dummy_group,
      linetype = series,
      group    = interaction(dummy_group, series)
    )
  ) +
    # shade the baseline normalization window
    ggplot2::annotate(
      "rect",
      xmin  = min(baseline_window),
      xmax  = max(baseline_window),
      ymin  = -Inf,
      ymax  = Inf,
      alpha = 0.06,
      fill  = "grey50"
    ) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    ggplot2::geom_vline(xintercept = -0.5, linetype = "dotted", color = "grey60")
  
  # SE ribbons drawn before lines/points so they sit underneath
  if (se_ribbon) {
    p <- p + ggplot2::geom_ribbon(
      ggplot2::aes(
        ymin = mean_outcome_norm - se_outcome,
        ymax = mean_outcome_norm + se_outcome
      ),
      alpha = ribbon_alpha,
      color = NA   # no ribbon border — would clash with linetype mapping
    )
  }
  
  p <- p +
    ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::geom_point(size = 1.5) +
    ggplot2::scale_color_manual(
      values = palette,
      labels = group_labels,
      name   = NULL
    ) +
    ggplot2::scale_fill_manual(
      values = palette,
      labels = group_labels,
      name   = NULL
    ) +
    ggplot2::scale_linetype_manual(
      values = series_linetypes,
      labels = series_labels,
      name   = NULL
    ) +
    ggplot2::labs(
      x     = "Event time (years since fire)",
      y     = glue::glue("Mean {outcome_filter} (normalized to baseline)"),
      title = title
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "bottom")
  
  # ── Faceting ─────────────────────────────────────────────────────────────────
  if (facet_by_dummy) {
    p <- p + ggplot2::facet_wrap(
      ggplot2::vars(dummy_group_label),
      scales = if (free_y) "free_y" else "fixed"
    )
  }
  
  p
}


# ══════════════════════════════════════════════════════════════════════════════
# Normalize and plot one run ----
# ══════════════════════════════════════════════════════════════════════════════

# Normalize the full table once upfront; the plot function just filters into it.
traj_norm <- normalize_to_baseline(traj_all, baseline_window = baseline_window)

# Single example run — swap these filters to explore other combinations.
plot_descriptive_trajectory(
  traj             = traj_norm,
  subset_id_filter = "nfg_Lodgepole_Pine_Group",
  outcome_filter   = "rap_tree",
  group_id_filter  = "b10_pdsisumn10",
  model_id_filter  = "sunab_twfe_unweighted",
  event_time_range = c(-15, 20),
  baseline_window  = baseline_window,
  se_ribbon        = TRUE,
  min_n_points     = 20,
  min_n_events     = 3,
  facet_by_dummy   = FALSE
)


