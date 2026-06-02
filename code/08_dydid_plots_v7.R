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
dir_results <- here::here("results", version)
dir_figs    <- here::here("figs",    version)
dir.create(dir_figs, recursive = TRUE, showWarnings = FALSE)


# ══════════════════════════════════════════════════════════════════════════════
# Load merged tables ----
# ══════════════════════════════════════════════════════════════════════════════

dir_all      <- file.path(dir_results, "tables", "all")
dir_inf_all  <- file.path(dir_results, "tables", "inference", "all")

agg_es_tbl  <- arrow::read_parquet(file.path(dir_all, "agg_eventstudy.parquet"))
registry    <- arrow::read_parquet(file.path(dir_all, "run_registry.parquet"))
dg_summary  <- arrow::read_parquet(file.path(dir_all, "dummy_group_summary.parquet"))

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
  "cd_f"   = "#E69F00",
  "cd_bf"  = "#56B4E9",
  "cd_df"  = "#009E73",
  "cd_bdf" = "#CC79A7"
)


# ══════════════════════════════════════════════════════════════════════════════
# Event study plot ----
# ══════════════════════════════════════════════════════════════════════════════

# TODO: port and generalize plot_event_study from v6.
# Key changes from v6:
#   - column is dummy_group (not subgroup)
#   - facet by subset_id; color/shape by dummy_group
#   - vcov_id filter for sensitivity ribbon comparison
#   - reference period line at ref.p (read from run_spec via agg_es_tbl metadata
#     or hardcode -6)
#   - event_time on x axis, estimate on y axis
#   - zero reference line, pre/post shading, N annotation optional

plot_event_study <- function(agg_es,
                             subset_id_filter,
                             outcome_filter,
                             group_id_filter,
                             model_id_filter,
                             vcov_id_filter     = NULL,
                             event_time_range   = c(-15, 15),
                             ref_period         = -6,
                             palette            = dummy_group_palette,
                             group_labels       = dummy_group_labels,
                             title              = NULL,
                             facet_by_dummy     = TRUE,
                             free_y             = FALSE) {
  
  d <- agg_es |>
    dplyr::filter(
      subset_id == subset_id_filter,
      outcome   == outcome_filter,
      group_id  == group_id_filter,
      model_id  == model_id_filter,
      event_time >= event_time_range[1],
      event_time <= event_time_range[2]
    )
  
  if (!is.null(vcov_id_filter)) d <- dplyr::filter(d, vcov_id == vcov_id_filter)
  
  if (nrow(d) == 0) {
    warning("No rows after filtering. Check subset/outcome/group/model filters.")
    return(NULL)
  }
  
  d <- d |>
    dplyr::mutate(
      dummy_group_label = factor(
        dplyr::recode(dummy_group, !!!group_labels),
        levels = group_labels[dummy_group_labels %in% names(group_labels)]
      )
    )
  
  p <- ggplot2::ggplot(d, ggplot2::aes(x = event_time, y = estimate,
                                       color = dummy_group, fill = dummy_group,
                                       group = dummy_group)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.15, color = NA) +
    ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::geom_point(size = 1.5) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    ggplot2::geom_vline(xintercept = ref_period, linetype = "dotted", color = "grey60") +
    ggplot2::annotate("rect", xmin = min(d$event_time), xmax = -0.5,
                      ymin = -Inf, ymax = Inf, alpha = 0.03, fill = "grey50") +
    ggplot2::scale_color_manual(values = palette, labels = group_labels, name = NULL) +
    ggplot2::scale_fill_manual(values = palette,  labels = group_labels, name = NULL) +
    ggplot2::labs(
      x     = "Event time (years since fire)",
      y     = glue::glue("ATT: {outcome_filter}"),
      title = title %||% glue::glue("{subset_id_filter} | {outcome_filter} | {group_id_filter} | {model_id_filter}")
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "bottom")
  
  if (facet_by_dummy) {
    scale_y <- if (free_y) "free_y" else "fixed"
    p <- p + ggplot2::facet_wrap(~ dummy_group_label, scales = scale_y)
  }
  
  p
}


# ══════════════════════════════════════════════════════════════════════════════
# ATT window dot plot ----
# ══════════════════════════════════════════════════════════════════════════════

# TODO: plot ATT window estimates per dummy_group across subset_ids.
# Intended use: forest-plot style comparison of ATTs across ecoregions,
# optionally faceted by window_id.

plot_att_windows <- function(att_windows,
                             outcome_filter,
                             model_id_filter,
                             vcov_id_filter,
                             window_id_filter = NULL,
                             palette          = dummy_group_palette,
                             group_labels     = dummy_group_labels) {
  stop("plot_att_windows() not yet implemented.")
}


# ══════════════════════════════════════════════════════════════════════════════
# Pairwise comparison dot plot ----
# ══════════════════════════════════════════════════════════════════════════════

# TODO: plot wald_compare_att results across subsets.
# Color by significance (p_bh < 0.05), facet by contrast.

plot_att_comparisons <- function(att_comps,
                                 outcome_filter,
                                 model_id_filter,
                                 vcov_id_filter,
                                 window_id_filter = NULL) {
  stop("plot_att_comparisons() not yet implemented.")
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