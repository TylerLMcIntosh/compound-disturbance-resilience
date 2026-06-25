# 08_dydid_plots_v7.R
# Plotting — skeleton
# -------------------------------------------------------
# Reads from:
#   tables/all/agg_eventstudy.parquet          — aggregated event studies
#   tables/inference/all/*.parquet             — ATT windows, comparisons, etc.
#   tables/all/dummy_group_summary.parquet     — N by dummy_group
#   tables/all/event_time_support.parquet      — N by event_time x dummy_group
# -------------------------------------------------------

rm(list = ls())

library(here)
here::i_am("code/08_dydid_plots_v7.R")

library(dplyr); library(ggplot2); library(readr); library(purrr)
library(tibble); library(stringr); library(arrow); library(glue)

source(here::here("code", "weight_dydid_pipeline_v7.R"))

run_name    <- "GEE_resilience_v7_operational_ss500_ts50000"
version     <- "v7"
dir_results <- here::here("results-exo", version)
dir_figs    <- here::here("figs",    version, "exo")
dir_figs_ecor <- here(dir_figs, "ecor")
dir_figs_nfg <- here(dir_figs, "nfg")

dir.create(dir_figs, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_figs_ecor, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_figs_nfg, recursive = TRUE, showWarnings = FALSE)


# ══════════════════════════════════════════════════════════════════════════════
# Load merged tables ----
# ══════════════════════════════════════════════════════════════════════════════

dir_all      <- file.path(dir_results, "tables", "all")
dir_inf_all  <- file.path(dir_results, "tables", "inference", "all")

agg_es_tbl  <- arrow::read_parquet(file.path(dir_all, "agg_event_study.parquet"))
registry    <- arrow::read_parquet(file.path(dir_all, "run_registry.parquet"))
dg_summary  <- arrow::read_parquet(file.path(dir_all, "dummy_group_summary.parquet"))
support <- arrow::read_parquet(file.path(dir_all, "event_time_support.parquet"))


att_windows <- arrow::read_parquet(file.path(dir_inf_all, "att_windows.parquet"))
att_comps   <- arrow::read_parquet(file.path(dir_inf_all, "att_comparisons.parquet"))
post_slopes <- arrow::read_parquet(file.path(dir_inf_all, "post_slopes.parquet"))


# ══════════════════════════════════════════════════════════════════════════════
# Plot helpers ----
# ══════════════════════════════════════════════════════════════════════════════

# Dummy group labels for axes
dummy_group_labels <- c(
  "cd_f"   = "Fire only",
  "cd_bf"  = "Biotic stress + fire",
  "cd_df"  = "Drought + fire",
  "cd_bdf" = "Drought + biotic + fire"
)

# Default palette (override via group_palette from Script 6)
dummy_group_palette <- c(
  "cd_f"   = "maroon3",
  "cd_bf"  = "forestgreen",
  "cd_df"  = "goldenrod2",
  "cd_bdf" = "purple2"
)


# ══════════════════════════════════════════════════════════════════════════════
# Event study plot ----
# ══════════════════════════════════════════════════════════════════════════════

plot_event_study <- function(agg_es,
                             subset_id_filter    = NULL,
                             outcome_filter,
                             group_id_filter,
                             model_id_filter,
                             vcov_id_filter      = NULL,
                             event_time_range    = c(-15, 20),
                             ref_period          = -6,
                             palette             = dummy_group_palette,
                             group_labels        = dummy_group_labels,
                             title               = NULL,
                             facet_by_dummy      = TRUE,
                             facet               = NULL,
                             free_y              = FALSE,
                             support             = NULL,
                             min_n_events        = NULL,
                             min_n_points        = NULL) {
  
  # ── Filter ────────────────────────────────────────────────────────────────
  d <- agg_es |>
    dplyr::filter(
      if (is.null(subset_id_filter)) TRUE else subset_id %in% subset_id_filter,
      outcome    %in% outcome_filter,
      group_id   %in% group_id_filter,
      model_id   %in% model_id_filter,
      event_time >= event_time_range[1],
      event_time <= event_time_range[2]
    )
  
  if (!is.null(vcov_id_filter)) {
    d <- dplyr::filter(d, vcov_id %in% vcov_id_filter)
  }
  
  if (nrow(d) == 0) {
    warning("No rows after filtering. Check subset/outcome/group/model filters.")
    return(NULL)
  }
  
  d <- d |>
    dplyr::mutate(
      dummy_group_label = factor(
        dplyr::recode(dummy_group, !!!group_labels),
        levels = unname(group_labels)
      )
    )
  
  # ── Support threshold filtering ───────────────────────────────────────────
  d_plot <- d
  
  if (!is.null(support) && (!is.null(min_n_events) || !is.null(min_n_points))) {
    support_sub <- support |>
      dplyr::filter(
        if (is.null(subset_id_filter)) TRUE else subset_id %in% subset_id_filter,
        outcome    %in% outcome_filter,
        group_id   %in% group_id_filter,
        model_id   %in% model_id_filter,
        event_time >= event_time_range[1],
        event_time <= event_time_range[2]
      ) |>
      dplyr::select(
        subset_id, outcome, group_id, model_id,
        event_time, dummy_group, n_fireids, n_ptids
      )
    
    d_plot <- d |>
      dplyr::left_join(
        support_sub,
        by = c(
          "subset_id", "outcome", "group_id", "model_id",
          "event_time", "dummy_group"
        )
      )
    
    if (!is.null(min_n_events)) {
      d_plot <- dplyr::filter(d_plot, !is.na(n_fireids), n_fireids >= min_n_events)
    }
    
    if (!is.null(min_n_points)) {
      d_plot <- dplyr::filter(d_plot, !is.na(n_ptids), n_ptids >= min_n_points)
    }
  }
  
  if (nrow(d_plot) == 0) {
    warning("No rows remain after support filtering. Try relaxing min_n_events or min_n_points.")
    return(NULL)
  }
  
  # ── Build plot ────────────────────────────────────────────────────────────
  p <- ggplot2::ggplot(
    d_plot,
    ggplot2::aes(
      x = event_time,
      y = estimate,
      color = dummy_group,
      fill = dummy_group,
      group = dummy_group
    )
  ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = ci_lower, ymax = ci_upper),
      alpha = 0.15,
      color = NA
    ) +
    ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::geom_point(size = 1.5) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    ggplot2::geom_vline(xintercept = ref_period, linetype = "dotted", color = "grey60") +
    ggplot2::annotate(
      "rect",
      xmin = min(d_plot$event_time),
      xmax = -0.5,
      ymin = -Inf,
      ymax = Inf,
      alpha = 0.03,
      fill = "grey50"
    ) +
    ggplot2::scale_color_manual(values = palette, labels = group_labels, name = NULL) +
    ggplot2::scale_fill_manual(values = palette, labels = group_labels, name = NULL) +
    ggplot2::labs(
      x     = "Event time (years since fire)",
      y     = "ATT",
      title = title
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "bottom")
  
  # ── Faceting ──────────────────────────────────────────────────────────────
  facet_vars <- c(
    if (facet_by_dummy) "dummy_group_label" else NULL,
    facet
  )
  
  if (length(facet_vars) > 0) {
    p <- p + ggplot2::facet_wrap(
      ggplot2::vars(!!!rlang::syms(facet_vars)),
      scales = if (free_y) "free_y" else "fixed"
    )
  }
  
  p
}

# ══════════════════════════════════════════════════════════════════════════════
# ATT window dot plot ----
# ══════════════════════════════════════════════════════════════════════════════

plot_att_windows <- function(att_windows,
                             outcome_filter,
                             group_id_filter,
                             model_id_filter,
                             vcov_id_filter,
                             subset_id_filter  = NULL,
                             window_id_filter  = NULL,
                             palette           = dummy_group_palette,
                             group_labels      = dummy_group_labels,
                             title             = NULL,
                             facet_by_window   = TRUE,
                             facet             = NULL,
                             free_x            = FALSE,
                             dodge_width       = 0.6,
                             sig_threshold     = NULL) {
  
  # ── Filter ────────────────────────────────────────────────────────────────
  d <- att_windows |>
    dplyr::filter(
      outcome  %in% outcome_filter,
      group_id %in% group_id_filter,
      model_id %in% model_id_filter,
      vcov_id  %in% vcov_id_filter
    )
  
  if (!is.null(subset_id_filter)) d <- dplyr::filter(d, subset_id %in% subset_id_filter)
  if (!is.null(window_id_filter)) d <- dplyr::filter(d, window_id %in% window_id_filter)
  
  if (nrow(d) == 0) {
    warning("No rows after filtering. Check outcome/group/model/vcov filters.")
    return(NULL)
  }
  
  d <- d |>
    dplyr::mutate(
      dummy_group_label = factor(
        dplyr::recode(dummy_group, !!!group_labels),
        levels = unname(group_labels)
      ),
      # clean subset label for y-axis
      subset_label = stringr::str_remove(subset_id, "^ecor_")
    )
  
  # optionally flag significance
  if (!is.null(sig_threshold)) {
    d <- d |> dplyr::mutate(significant = !is.na(p) & p < sig_threshold)
  }
  
  # ── Build plot ────────────────────────────────────────────────────────────
  pos <- ggplot2::position_dodge(width = dodge_width)
  
  p <- ggplot2::ggplot(
    d,
    ggplot2::aes(
      x     = estimate,
      y     = stats::reorder(subset_label, estimate),
      color = dummy_group,
      shape = if (!is.null(sig_threshold)) significant else NULL
    )
  ) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
    ggplot2::geom_pointrange(
      ggplot2::aes(xmin = ci_lower, xmax = ci_upper),
      position = pos,
      linewidth = 0.5,
      size      = 0.4
    ) +
    ggplot2::scale_color_manual(values = palette, labels = group_labels, name = NULL) +
    ggplot2::labs(
      x     = "ATT (95% CI)",
      y     = NULL,
      title = title,
      shape = if (!is.null(sig_threshold)) glue::glue("p < {sig_threshold}") else NULL
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position  = "bottom",
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor   = ggplot2::element_blank()
    )
  
  # ── Faceting ──────────────────────────────────────────────────────────────
  facet_vars <- c(
    if (facet_by_window && dplyr::n_distinct(d$window_id) > 1) "window_id" else NULL,
    facet
  )
  
  if (length(facet_vars) > 0) {
    p <- p + ggplot2::facet_wrap(
      ggplot2::vars(!!!rlang::syms(facet_vars)),
      scales = if (free_x) "free_x" else "fixed"
    )
  }
  
  p
}


# ══════════════════════════════════════════════════════════════════════════════
# Pairwise comparison dot plot ----
# ══════════════════════════════════════════════════════════════════════════════

# plot_att_comparisons <- function(att_comps,
#                                  outcome_filter,
#                                  group_id_filter,
#                                  model_id_filter,
#                                  vcov_id_filter,
#                                  subset_id_filter  = NULL,
#                                  window_id_filter  = NULL,
#                                  sig_col           = "p_bh",
#                                  sig_threshold     = 0.05,
#                                  facet_by_contrast = TRUE,
#                                  facet_by_window   = FALSE,
#                                  facet             = NULL,
#                                  free_x            = FALSE,
#                                  dodge_width       = 0.6,
#                                  title             = NULL) {
#   
#   # ── Filter ────────────────────────────────────────────────────────────────
#   d <- att_comps |>
#     dplyr::filter(
#       outcome  %in% outcome_filter,
#       group_id %in% group_id_filter,
#       model_id %in% model_id_filter,
#       vcov_id  %in% vcov_id_filter
#     )
#   
#   if (!is.null(subset_id_filter)) d <- dplyr::filter(d, subset_id %in% subset_id_filter)
#   if (!is.null(window_id_filter)) d <- dplyr::filter(d, window_id %in% window_id_filter)
#   
#   if (nrow(d) == 0) {
#     warning("No rows after filtering. Check outcome/group/model/vcov filters.")
#     return(NULL)
#   }
#   
#   if (!sig_col %in% names(d)) {
#     stop(glue::glue("sig_col '{sig_col}' not found. Available: ",
#                     paste(names(d), collapse = ", ")))
#   }
#   
#   d <- d |>
#     dplyr::mutate(
#       subset_label = stringr::str_remove(subset_id, "^ecor_"),
#       significant  = !is.na(.data[[sig_col]]) & .data[[sig_col]] < sig_threshold,
#       sig_label    = dplyr::if_else(significant,
#                                     glue::glue("p < {sig_threshold}"),
#                                     glue::glue("p \u2265 {sig_threshold}")) |>
#         factor(levels = c(glue::glue("p < {sig_threshold}"),
#                           glue::glue("p \u2265 {sig_threshold}")))
#     )
#   
#   # ── Build plot ────────────────────────────────────────────────────────────
#   pos <- ggplot2::position_dodge(width = dodge_width)
#   
#   p <- ggplot2::ggplot(
#     d,
#     ggplot2::aes(
#       x     = estimate,
#       y     = stats::reorder(subset_label, estimate),
#       color = sig_label,
#       group = window_id
#     )
#   ) +
#     ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
#     ggplot2::geom_pointrange(
#       ggplot2::aes(xmin = ci_lower, xmax = ci_upper),
#       position  = pos,
#       linewidth = 0.5,
#       size      = 0.4
#     ) +
#     ggplot2::scale_color_manual(
#       values = stats::setNames(
#         c("#D62728", "grey60"),
#         c(glue::glue("p < {sig_threshold}"),
#           glue::glue("p \u2265 {sig_threshold}"))
#       ),
#       drop = FALSE,
#       name = glue::glue("{sig_col}")
#     ) +
#     ggplot2::labs(
#       x     = "Difference in ATT (95% CI)",
#       y     = NULL,
#       title = title
#     ) +
#     ggplot2::theme_bw() +
#     ggplot2::theme(
#       legend.position    = "bottom",
#       panel.grid.major.y = ggplot2::element_blank(),
#       panel.grid.minor   = ggplot2::element_blank()
#     )
#   
#   # ── Faceting ──────────────────────────────────────────────────────────────
#   facet_vars <- c(
#     if (facet_by_contrast && dplyr::n_distinct(d$contrast) > 1) "contrast"  else NULL,
#     if (facet_by_window   && dplyr::n_distinct(d$window_id) > 1) "window_id" else NULL,
#     facet
#   )
#   
#   if (length(facet_vars) > 0) {
#     p <- p + ggplot2::facet_wrap(
#       ggplot2::vars(!!!rlang::syms(facet_vars)),
#       scales = if (free_x) "free_x" else "fixed"
#     )
#   }
#   
#   p
# }


plot_att_comparisons <- function(att_comps,
                                 outcome_filter,
                                 group_id_filter,
                                 model_id_filter,
                                 vcov_id_filter,
                                 subset_id_filter  = NULL,
                                 window_id_filter  = NULL,
                                 y_var             = "subset_id",
                                 color_var         = "sig_label",
                                 shape_var         = NULL,
                                 color_values      = NULL,
                                 sig_col           = "p_bh",
                                 sig_threshold     = 0.05,
                                 facet_by_contrast = TRUE,
                                 facet_by_window   = FALSE,
                                 facet             = NULL,
                                 free_x            = FALSE,
                                 dodge_width       = 0.6,
                                 title             = NULL) {
  
  # ── Filter ────────────────────────────────────────────────────────────────
  d <- att_comps |>
    dplyr::filter(
      outcome  %in% outcome_filter,
      group_id %in% group_id_filter,
      model_id %in% model_id_filter,
      vcov_id  %in% vcov_id_filter
    )
  
  if (!is.null(subset_id_filter)) d <- dplyr::filter(d, subset_id %in% subset_id_filter)
  if (!is.null(window_id_filter)) d <- dplyr::filter(d, window_id %in% window_id_filter)
  
  if (nrow(d) == 0) {
    warning("No rows after filtering. Check outcome/group/model/vcov filters.")
    return(NULL)
  }
  
  for (v in c(y_var, color_var, shape_var)) {
    if (!is.null(v) && !v %in% names(d)) {
      stop(glue::glue("Column '{v}' not found in att_comps."))
    }
  }
  
  if (!sig_col %in% names(d)) {
    stop(glue::glue("sig_col '{sig_col}' not found. Available: ",
                    paste(names(d), collapse = ", ")))
  }
  
  sig_lab_pass <- glue::glue("p < {sig_threshold}")
  sig_lab_fail <- glue::glue("p \u2265 {sig_threshold}")
  
  d <- d |>
    dplyr::mutate(
      subset_label = stringr::str_remove(subset_id, "^ecor_"),
      significant  = !is.na(.data[[sig_col]]) & .data[[sig_col]] < sig_threshold,
      sig_label    = factor(
        dplyr::if_else(significant, sig_lab_pass, sig_lab_fail),
        levels = c(sig_lab_pass, sig_lab_fail)
      )
    )
  
  # ── Y variable ───────────────────────────────────────────────────────────
  # if the column is already a factor, respect its levels; otherwise reorder
  # by estimate only when y_var is a subset identifier
  y_col <- if (y_var == "subset_id") "subset_label" else y_var
  
  if (!is.factor(d[[y_col]])) {
    if (y_var == "subset_id") {
      d <- dplyr::mutate(d, !!y_col := stats::reorder(.data[[y_col]], estimate))
    }
    # for other y_vars leave the natural order intact
  }
  
  # ── Color variable and palette ────────────────────────────────────────────
  # if color_var is sig_label, apply the sig/non-sig coloring scheme;
  # if it's a user column, apply color_values with non-significant entries
  # overridden to grey
  if (color_var == "sig_label" || is.null(color_var)) {
    color_col    <- "sig_label"
    color_levels <- levels(d$sig_label)
    final_colors <- stats::setNames(c("#D62728", "grey60"), color_levels)
    color_name   <- sig_col
    
  } else {
    color_col    <- color_var
    color_levels <- if (is.factor(d[[color_col]])) levels(d[[color_col]]) else sort(unique(d[[color_col]]))
    
    if (is.null(color_values)) {
      # default palette if none supplied
      color_values <- stats::setNames(
        scales::hue_pal()(length(color_levels)),
        color_levels
      )
    } else {
      if (length(color_values) < length(color_levels)) {
        stop("color_values must have at least as many entries as unique values of color_var.")
      }
      color_values <- stats::setNames(color_values[seq_along(color_levels)], color_levels)
    }
    
    # non-significant points: override their color to grey regardless of color_var
    d <- d |>
      dplyr::mutate(
        .color_actual = dplyr::if_else(significant,
                                       as.character(.data[[color_col]]),
                                       "grey60")
      )
    # build expanded palette including grey60
    final_colors <- c(color_values, "grey60" = "grey60")
    color_col    <- ".color_actual"
    color_name   <- color_var
  }
  
  # ── Build plot ────────────────────────────────────────────────────────────
  pos <- ggplot2::position_dodge(width = dodge_width)
  
  aes_args <- list(
    x     = ggplot2::aes(x = estimate),
    y     = ggplot2::aes(y = .data[[y_col]]),
    color = ggplot2::aes(color = .data[[color_col]])
  )
  
  base_aes <- ggplot2::aes(
    x     = estimate,
    y     = .data[[y_col]],
    color = .data[[color_col]],
    group = .data[[color_col]]
  )
  
  if (!is.null(shape_var)) {
    base_aes <- utils::modifyList(
      base_aes,
      ggplot2::aes(shape = .data[[shape_var]])
    )
  }
  
  p <- ggplot2::ggplot(d, base_aes) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
    ggplot2::geom_pointrange(
      ggplot2::aes(xmin = ci_lower, xmax = ci_upper),
      position  = pos,
      linewidth = 0.5,
      size      = 0.4
    ) +
    ggplot2::scale_color_manual(
      values = final_colors,
      drop   = FALSE,
      name   = color_name
    ) +
    ggplot2::labs(
      x     = "Difference in ATT (95% CI)",
      y     = NULL,
      title = title,
      shape = shape_var
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position    = "bottom",
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor   = ggplot2::element_blank()
    )
  
  # ── Faceting ──────────────────────────────────────────────────────────────
  facet_vars <- c(
    if (facet_by_contrast && dplyr::n_distinct(d$contrast) > 1) "contrast"  else NULL,
    if (facet_by_window   && dplyr::n_distinct(d$window_id) > 1) "window_id" else NULL,
    facet
  )
  
  if (length(facet_vars) > 0) {
    p <- p + ggplot2::facet_wrap(
      ggplot2::vars(!!!rlang::syms(facet_vars)),
      scales = if (free_x) "free_x" else "fixed"
    )
  }
  
  p
}




# ══════════════════════════════════════════════════════════════════════════════
# Vcov sensitivity comparison ----
# ══════════════════════════════════════════════════════════════════════════════

# TODO: overlay event study from multiple vcov_ids (one line per vcov_id)
# for a single subset x outcome x group x model, to visualize SE sensitivity.

plot_vcov_sensitivity <- function(agg_es,
                                   subset_id_filter,
                                   outcome_filter,
                                   group_id_filter,
                                   model_id_filter,
                                   dummy_group_filter) {
  stop("plot_vcov_sensitivity() not yet implemented.")
}



# Comparisons to F line ----


# F1 b

plot_att_comparisons_vs_f <- function(att_comps,
                                      outcome_filter,
                                      group_id_filter,
                                      model_id_filter,
                                      vcov_id_filter,
                                      subset_id_filter = NULL,
                                      window_order     = c("years_2_6", "years_7_11", "years_12_16"),
                                      window_labels    = c("years_2_6"   = "Years 2\u20136",
                                                           "years_7_11"  = "Years 7\u201311",
                                                           "years_12_16" = "Years 12\u201316"),
                                      palette          = dummy_group_palette,
                                      group_labels     = dummy_group_labels,
                                      sig_col          = "p_bh",
                                      sig_threshold    = 0.05,
                                      dodge_width      = 0.5,
                                      title            = NULL) {
  
  # ── Filter ────────────────────────────────────────────────────────────────
  d <- att_comps |>
    dplyr::filter(
      outcome  %in% outcome_filter,
      group_id %in% group_id_filter,
      model_id %in% model_id_filter,
      vcov_id  %in% vcov_id_filter,
      group_a  == "cd_f"            # only comparisons against fire-only
    )
  
  if (!is.null(subset_id_filter)) d <- dplyr::filter(d, subset_id %in% subset_id_filter)
  
  if (nrow(d) == 0) {
    warning("No rows after filtering. Check filters and that group_a == 'cd_f' exists.")
    return(NULL)
  }
  
  if (!sig_col %in% names(d)) {
    stop(glue::glue("sig_col '{sig_col}' not found. Available: ",
                    paste(names(d), collapse = ", ")))
  }
  
  # ── Prepare ───────────────────────────────────────────────────────────────
  # window factor: earliest at top means reversed order on y axis
  # ggplot2 y-axis: first factor level = bottom, last = top
  window_levels_btm_to_top <- rev(window_order)
  
  # facet labels: readable name of group_b (the group being compared to F)
  # drop cd_f from group_labels since it's never group_b here
  facet_level_order <- unname(group_labels[names(group_labels) != "cd_f"])
  
  d <- d |>
    dplyr::mutate(
      significant  = !is.na(.data[[sig_col]]) & .data[[sig_col]] < sig_threshold,
      # significant: use group_b's color; non-significant: grey
      point_color  = dplyr::if_else(significant, group_b, "grey60"),
      facet_label  = factor(
        dplyr::recode(group_b, !!!group_labels),
        levels = facet_level_order
      ),
      window_id    = factor(
        dplyr::recode(window_id, !!!window_labels),
        levels = dplyr::recode(window_levels_btm_to_top, !!!window_labels)
      ),
      subset_label = stringr::str_remove(subset_id, "^ecor_")
    )
  
  # palette: one entry per group_b value + grey for non-significant
  sig_colors <- palette[unique(d$group_b[d$significant])]
  all_colors <- c(sig_colors, "grey60" = "grey60")
  
  # ── Build plot ────────────────────────────────────────────────────────────
  pos <- ggplot2::position_dodge(width = dodge_width)
  
  p <- ggplot2::ggplot(
    d,
    ggplot2::aes(
      x     = estimate,
      y     = window_id,
      color = point_color,
      group = subset_label
    )
  ) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
    ggplot2::geom_pointrange(
      ggplot2::aes(xmin = ci_lower, xmax = ci_upper),
      position  = pos,
      linewidth = 0.5,
      size      = 0.4
    ) +
    ggplot2::scale_color_manual(values = all_colors, guide = "none") +
    ggplot2::labs(
      x     = "Difference in ATT vs fire only (95% CI)",
      y     = "Post-treatment window",
      title = title
    ) +
    ggplot2::facet_wrap(~ facet_label, ncol = 3) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position    = "none",
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor   = ggplot2::element_blank()
    )
  
  p
}


plot_att_comparisons_vs_f_windows <- function(att_comps,
                                              outcome_filter,
                                              group_id_filter,
                                              model_id_filter,
                                              vcov_id_filter,
                                              subset_id_filter = NULL,
                                              window_order     = c("years_2_6", "years_7_11", "years_12_16"),
                                              window_labels    = c("years_2_6"   = "Years 2\u20136",
                                                                   "years_7_11"  = "Years 7\u201311",
                                                                   "years_12_16" = "Years 12\u201316"),
                                              window_colors    = c("years_2_6"   = "#FFD700",
                                                                   "years_7_11"  = "#FF8C00",
                                                                   "years_12_16" = "#CC0000"),
                                              window_shapes    = c("years_2_6"   = 15,
                                                                   "years_7_11"  = 17,
                                                                   "years_12_16" = 19),
                                              group_labels     = dummy_group_labels,
                                              sig_col          = "p_bh",
                                              sig_threshold    = 0.05,
                                              dodge_width      = 0.5,
                                              title            = NULL,
                                              sort             = TRUE,
                                              xlim             = NULL) {
  
  # ── Filter ────────────────────────────────────────────────────────────────
  d <- att_comps |>
    dplyr::filter(
      outcome   %in% outcome_filter,
      group_id  %in% group_id_filter,
      model_id  %in% model_id_filter,
      vcov_id   %in% vcov_id_filter,
      group_a   == "cd_f",
      window_id %in% window_order
    )
  
  if (!is.null(subset_id_filter)) {
    d <- dplyr::filter(d, subset_id %in% subset_id_filter)
  }
  
  if (nrow(d) == 0) {
    warning("No rows after filtering.")
    return(NULL)
  }
  
  if (!sig_col %in% names(d)) {
    stop(glue::glue("sig_col '{sig_col}' not found. Available: ",
                    paste(names(d), collapse = ", ")))
  }
  
  if (!is.null(xlim) && length(xlim) != 2) {
    stop("xlim must be NULL or a numeric vector of length 2, e.g. c(-10, 10).")
  }
  
  # ── Prepare ───────────────────────────────────────────────────────────────
  window_label_vec    <- window_labels[window_order]
  window_levels_y     <- rev(unname(window_label_vec))
  window_color_mapped <- stats::setNames(
    unname(window_colors[window_order]),
    unname(window_label_vec)
  )
  window_shape_mapped <- stats::setNames(
    unname(window_shapes[window_order]),
    unname(window_label_vec)
  )
  
  facet_level_order <- unname(group_labels[names(group_labels) != "cd_f"])
  
  d <- d |>
    dplyr::mutate(
      significant   = !is.na(.data[[sig_col]]) & .data[[sig_col]] < sig_threshold,
      window_label  = factor(
        dplyr::recode(window_id, !!!window_label_vec),
        levels = window_levels_y
      ),
      point_color   = dplyr::if_else(
        significant,
        dplyr::recode(window_id, !!!window_colors[window_order]),
        "grey60"
      ),
      facet_label   = factor(
        dplyr::recode(group_b, !!!group_labels),
        levels = facet_level_order
      ),
      subset_label  = stringr::str_remove(subset_id, "^ecor_")
    )
  
  # ── Sort subsets ──────────────────────────────────────────────────────────
  if (sort) {
    subset_order <- d |>
      dplyr::group_by(subset_label) |>
      dplyr::summarise(
        n_nonsig    = sum(!significant),
        mean_effect = mean(estimate, na.rm = TRUE),
        .groups     = "drop"
      ) |>
      dplyr::arrange(dplyr::desc(n_nonsig), mean_effect) |>
      dplyr::pull(subset_label)
  } else {
    subset_order <- sort(unique(d$subset_label))
  }
  
  d <- d |>
    dplyr::mutate(
      subset_label = factor(subset_label, levels = subset_order)
    )
  
  all_colors <- c(window_color_mapped, "grey60" = "grey60")
  all_shapes <- window_shape_mapped
  
  # ── Build plot ────────────────────────────────────────────────────────────
  pos <- ggplot2::position_dodge(width = dodge_width)
  
  p <- ggplot2::ggplot(
    d,
    ggplot2::aes(
      x     = estimate,
      y     = subset_label,
      color = point_color,
      shape = window_label,
      group = window_label
    )
  ) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
    ggplot2::geom_pointrange(
      ggplot2::aes(xmin = ci_lower, xmax = ci_upper),
      position  = pos,
      linewidth = 0.5,
      size      = 0.5
    ) +
    ggplot2::scale_color_identity(
      guide  = "legend",
      name   = glue::glue("Window ({sig_col} < {sig_threshold})"),
      labels = c(unname(window_label_vec), "Non-significant"),
      breaks = c(unname(window_colors[window_order]), "grey60")
    ) +
    ggplot2::scale_shape_manual(
      values = all_shapes,
      name   = "Window"
    ) +
    ggplot2::labs(
      x     = "Difference in ATT vs fire only (95% CI)",
      y     = NULL,
      title = title
    ) +
    ggplot2::facet_wrap(~ facet_label, ncol = 3) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position    = "bottom",
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor   = ggplot2::element_blank()
    )
  
  if (!is.null(xlim)) {
    p <- p + ggplot2::coord_cartesian(xlim = xlim)
  }
  
  p
}


# USE FUNCTIONS ----

nfg_subsets <- agg_es_tbl |>
  filter(grepl("nfg_", subset_id)) |>
  pull(subset_id) |>
  unique()

ecor_subsets <- agg_es_tbl |>
  filter(grepl("ecor_", subset_id)) |>
  pull(subset_id) |>
  unique()


## event_study ----

for(i in nfg_subsets) {
  plot_event_study(
    agg_es           = agg_es_tbl,
    subset_id_filter = i,
    outcome_filter   = "rap_tree",
    group_id_filter  = "b10_pdsisumn10",
    model_id_filter  = "sunab_twfe_glmatotopoclimnfg",
    vcov_id_filter   = "conley_75km_5km",
    support          = support,
    min_n_events     = 3,
    min_n_points     = 20,
    facet_by_dummy = FALSE,
    title = glue(i, " glmatotopoclimnfg, conley")
  ) + xlim(-10, 20)
  
  ggsave(here(dir_figs_nfg, glue("event_study_", i, ".png")))
}


for(i in ecor_subsets) {
  plot_event_study(
    agg_es           = agg_es_tbl,
    subset_id_filter = i,
    outcome_filter   = "rap_tree",
    group_id_filter  = "b10_pdsisumn10",
    model_id_filter  = "sunab_twfe_glmatotopoclimnfg",
    vcov_id_filter   = "conley_75km_5km",
    support          = support,
    min_n_events     = 3,
    min_n_points     = 20,
    facet_by_dummy = FALSE,
    title = glue(i, " glmatotopoclimnfg, conley")
  ) + xlim(-10, 20)
  
  ggsave(here(dir_figs_ecor, glue("event_study_", i, ".png")))
}



## Plot comparisons to F ----


plot_att_comparisons_vs_f_windows(
  att_comps       = att_comps |>
    filter(grepl("nfg", subset_id)),
  outcome_filter  = "rap_tree",
  group_id_filter = "b10_pdsisumn10",
  model_id_filter = "sunab_twfe_unweighted",
  vcov_id_filter  = "conley_75km_5km",
  dodge_width     = 0.6,
  title = "Unweighted RAP, pt_cluster",
  sort = TRUE,
  xlim = c(-15, 15)
)
ggsave(here(dir_figs_nfg, glue("fire_comparison_summary_rap_unweighted_conleyclust.png")))


plot_att_comparisons_vs_f_windows(
  att_comps       = att_comps |> filter(grepl("nfg", subset_id)),
  outcome_filter  = "rap_tree",
  group_id_filter = "b10_pdsisumn10",
  model_id_filter = "sunab_twfe_glmatotopoclimnfg",
  vcov_id_filter  = "cluster_pt",
  dodge_width     = 0.6,
  title = "Weighted RAP - pt se",
  sort = TRUE,
  xlim = c(-15, 15)
)
ggsave(here(dir_figs_nfg, glue("fire_comparison_summary_rap_weighted_ptclust.png")))




plot_att_comparisons_vs_f_windows(
  att_comps       = att_comps |> filter(grepl("nfg", subset_id)),
  outcome_filter  = "rap_tree",
  group_id_filter = "b10_pdsisumn10",
  model_id_filter = "sunab_twfe_glmatotopoclimnfg",
  vcov_id_filter  = "cluster_h3",
  dodge_width     = 0.6,
  title = "Weighted RAP - h3 se",
  sort = FALSE,
  xlim = c(-15, 15)
)
ggsave(here(dir_figs_nfg, glue("fire_comparison_summary_rap_weighted_h3clust.png")))



plot_att_comparisons_vs_f_windows(
  att_comps       = att_comps |> filter(grepl("nfg", subset_id)),
  outcome_filter  = "rap_tree",
  group_id_filter = "b10_pdsisumn10",
  model_id_filter = "sunab_twfe_glmatotopoclimnfg",
  vcov_id_filter  = "conley_75km_5km",
  dodge_width     = 0.6,
  title = "Weighted RAP - conley se",
  sort = FALSE,
  xlim = c(-15, 15)
)
ggsave(here(dir_figs_nfg, glue("fire_comparison_summary_rap_weighted_conleyclust.png")))



plot_att_comparisons_vs_f_windows(
  att_comps       = att_comps |> filter(grepl("nfg", subset_id)),
  outcome_filter  = "vcf_tree",
  group_id_filter = "b10_pdsisumn10",
  model_id_filter = "sunab_twfe_glmatotopoclimnfg",
  vcov_id_filter  = "conley_75km_5km",
  dodge_width     = 0.6,
  title = "Weighted VCF - conley se",
  sort = TRUE,
  xlim = c(-15, 22)
)
ggsave(here(dir_figs_nfg, glue("fire_comparison_summary_vcf_weighted_conleyclust.png")))




# Examine support ----

for(i in nfg_subsets) {
  s <- support |> filter(subset_id == i & outcome == "rap_tree" & model_id == "sunab_twfe_unweighted")
  p <- ggplot(s) +
    geom_col(aes(x = event_time, y = n_ptids)) +
    facet_wrap(~dummy_group, scales = "free_y") +
    labs(title = glue("Support: {i}")) +
    theme_minimal() +
    xlim(-20, 20)
  
  ggsave(here(dir_figs_nfg, glue("support_{i}.png")))
}



for(i in nfg_subsets) {
  s <- support |>
    filter(
      subset_id == i,
      outcome == "rap_tree",
      model_id == "sunab_twfe_unweighted"
    )
  
  scale_factor <- max(s$n_ptids, na.rm = TRUE) / max(s$n_fireids, na.rm = TRUE)
  
  p <- ggplot(s, aes(x = event_time)) +
    geom_col(aes(y = n_ptids)) +
    geom_line(
      aes(y = n_fireids * scale_factor, group = dummy_group),
      linewidth = 1
    ) +
    facet_wrap(~ dummy_group, scales = "free_y") +
    scale_y_continuous(
      name = "Number of ptids (bars)",
      sec.axis = sec_axis(
        ~ . / scale_factor,
        name = "Number of fireids (line)"
      )
    ) +
    labs(title = glue("Support: {i}")) +
    theme_minimal() +
    coord_cartesian(xlim = c(-20, 20))
  
  ggsave(here(dir_figs_nfg, glue("support_{i}.png")), p)
}
