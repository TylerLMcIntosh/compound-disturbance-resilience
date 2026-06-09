# 07_dydid_aggregate_inference_v7.R
# Post-estimation aggregation and inference
# -------------------------------------------------------
# Inputs (per run, from 06_dydid_estimation_v7.R output):
#   tables/by_run/{run_stub}/
#     agg_{vcov_id}.rds    — aggregated coef + vcov (list) for each vcov spec
#     model.rds            — lean fixest model (primary vcov; for cohort-group agg)
#     run_spec.rds         — full run spec including dummy_cols, vcov_specs, etc.
#     registry.parquet     — one row per dummy_group
#
# Outputs (tables/by_run/{run_stub}/inference/):
#   pretrend_flat_{vcov_id}.parquet
#   pretrend_slope_{vcov_id}.parquet
#   att_windows_{vcov_id}.parquet
#   post_slopes_{vcov_id}.parquet
#   att_comparisons_{vcov_id}.parquet     (with BH + Holm corrections)
#   slope_comparisons_{vcov_id}.parquet
#   cohortgroup_att_windows_primary.parquet
#   cohortgroup_post_slopes_primary.parquet
#
# These are merged by rebuild_inference_tables() into tables/inference/all/.
# -------------------------------------------------------


Sys.setenv(LD_LIBRARY_PATH = paste("/opt/conda/lib", Sys.getenv("LD_LIBRARY_PATH"), sep = ":"))
Sys.setenv(PATH = paste("/usr/bin:/bin:/usr/local/bin", Sys.getenv("PATH"), sep = ":"))

rm(list = ls())

if (!requireNamespace("here", quietly = TRUE)) install.packages("here")
library(here)
here::i_am("code/07_dydid_aggregate_inference_v7.R")

required_pkgs <- c("dplyr", "readr", "purrr", "tibble", "stringr", "arrow", "glue", "here", "fixest")
missing_pkgs  <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) install.packages(missing_pkgs)

library(dplyr); library(readr); library(purrr); library(tibble)
library(stringr); library(arrow); library(glue); library(fixest)

source(here::here("code", "weight_dydid_pipeline_v7.R"))


# ══════════════════════════════════════════════════════════════════════════════
# 1. Configuration — point to estimation output directory ----
# ══════════════════════════════════════════════════════════════════════════════

run_name    <- "GEE_resilience_v7_operational_ss500_ts50000"
version     <- "v7"
dir_results <- here::here("results", version)
dir_infout  <- file.path(dir_results, "tables", "inference", "all")
dir.create(dir_infout, recursive = TRUE, showWarnings = FALSE)


# ══════════════════════════════════════════════════════════════════════════════
# 2. Inference specification ----
# ══════════════════════════════════════════════════════════════════════════════

# -- Pre-trend windows --------------------------------------------------------
# NOTE: ref.p from feols (here -6) should not be included in pre_years, as the
# reference period estimate is mechanically zero. Exclude it by design when
# setting the years vector.
pretrend_specs <- list(
  list(id = "pre_7_to_15", years = -15:-7)
)

# -- Post-treatment windows ---------------------------------------------------
posttreatment_window_specs <- list(
  list(id = "years_2_6",   years = 2:6),
  list(id = "years_2_15",  years = 2:15),
  list(id = "years_10_15", years = 10:15)
)

# -- Cohort group specs (for cohort-stratified aggregation from lean model) ---
# Only the primary vcov is available for cohort-group inference. Non-primary
# vcov cohort-group results will be point estimates only (no inference).
cohort_group_specs <- list(
  list(id = "early", min_cohort = 2000, max_cohort = 2010),
  list(id = "late",  min_cohort = 2011, max_cohort = 2020)
)

# -- Pairwise comparison pairs (reference dummy col first) --------------------
# Using dummy column names (cd_f, cd_bf, ...) directly.
comparison_pairs <- list(
  c("cd_f", "cd_bf"),
  c("cd_f", "cd_df"),
  c("cd_f", "cd_bdf")
)

# -- CI level (used for att_window point estimates) ---------------------------
ci_level <- 0.95


# ══════════════════════════════════════════════════════════════════════════════
# 3. Discover runs to process ----
# ══════════════════════════════════════════════════════════════════════════════

dir_by_run <- file.path(dir_results, "tables", "by_run")
run_spec_files <- list.files(dir_by_run, pattern = "^run_spec\\.rds$",
                             recursive = TRUE, full.names = TRUE)

if (length(run_spec_files) == 0) stop("No run_spec.rds files found under: ", dir_by_run)
message(glue::glue("Discovered {length(run_spec_files)} estimation runs to process."))


# ══════════════════════════════════════════════════════════════════════════════
# 4. Per-run inference helper ----
# ══════════════════════════════════════════════════════════════════════════════

run_one_inference <- function(run_spec_file,
                              pretrend_specs,
                              posttreatment_window_specs,
                              cohort_group_specs,
                              comparison_pairs,
                              ci_level,
                              skip_existing = TRUE) {
  
  run_spec  <- readRDS(run_spec_file)
  run_dir   <- dirname(run_spec_file)
  inf_dir   <- file.path(run_dir, "inference")
  dir.create(inf_dir, recursive = TRUE, showWarnings = FALSE)
  
  run_id      <- run_spec$run_id
  run_stub    <- run_spec$run_stub
  dummy_cols  <- run_spec$dummy_cols
  vcov_specs  <- run_spec$vcov_specs
  
  # run metadata carried into every output table
  meta <- list(
    run_id   = run_id, subset_id = run_spec$subset_id, outcome = run_spec$outcome,
    group_id = run_spec$group_id, model_id = run_spec$model_id
  )
  
  message(glue::glue("  Processing: {run_id}"))
  
  results <- list()
  
  # ── Per-vcov-spec inference ─────────────────────────────────────────────────
  for (vi in seq_len(nrow(vcov_specs))) {
    v_id <- vcov_specs$vcov_id[[vi]]
    agg_rds_file <- file.path(run_dir, glue::glue("agg_{safe_path_component(v_id)}.rds"))
    
    if (!file.exists(agg_rds_file)) {
      warning(glue::glue("[{run_id}] agg RDS not found for vcov_id={v_id}: {agg_rds_file}"))
      next
    }
    
    agg_obj <- readRDS(agg_rds_file)
    
    # pre-trend flatness tests
    pretrend_flat_rows <- purrr::map_dfr(pretrend_specs, function(ps) {
      purrr::map_dfr(dummy_cols, function(dc) {
        wald_pretrend_flat(agg_obj, dummy_group = dc, pre_years = ps$years) |>
          dplyr::mutate(pretrend_id = ps$id)
      })
    }) |>
      add_inference_meta(meta, v_id)
    
    # pre-trend slope tests
    pretrend_slope_rows <- purrr::map_dfr(pretrend_specs, function(ps) {
      purrr::map_dfr(dummy_cols, function(dc) {
        wald_pretrend_slope(agg_obj, dummy_group = dc, pre_years = ps$years) |>
          dplyr::mutate(pretrend_id = ps$id)
      })
    }) |>
      add_inference_meta(meta, v_id)
    
    # post-treatment ATT windows
    att_window_rows <- purrr::map_dfr(posttreatment_window_specs, function(ws) {
      purrr::map_dfr(dummy_cols, function(dc) {
        att_window(agg_obj, dummy_group = dc, years = ws$years) |>
          dplyr::mutate(window_id = ws$id)
      })
    }) |>
      add_inference_meta(meta, v_id)
    
    # post-treatment slopes
    post_slope_rows <- purrr::map_dfr(posttreatment_window_specs, function(ws) {
      purrr::map_dfr(dummy_cols, function(dc) {
        gls_slope_one_group(agg_obj, dummy_group = dc, years = ws$years) |>
          dplyr::mutate(window_id = ws$id)
      })
    }) |>
      add_inference_meta(meta, v_id)
    
    # pairwise ATT comparisons (with MTC)
    att_comp_rows_raw <- purrr::map_dfr(posttreatment_window_specs, function(ws) {
      purrr::map_dfr(comparison_pairs, function(pair) {
        wald_compare_att(agg_obj, group_a = pair[1], group_b = pair[2], years = ws$years) |>
          dplyr::mutate(window_id = ws$id)
      })
    }) |>
      add_inference_meta(meta, v_id)
    
    att_comp_rows <- apply_mtc(att_comp_rows_raw, p_col = "p", grouping = "window_id")
    
    # pairwise slope comparisons
    slope_comp_rows <- purrr::map_dfr(posttreatment_window_specs, function(ws) {
      purrr::map_dfr(comparison_pairs, function(pair) {
        compare_gls_slopes(agg_obj, group_a = pair[1], group_b = pair[2], years = ws$years) |>
          dplyr::mutate(window_id = ws$id)
      })
    }) |>
      add_inference_meta(meta, v_id) |>
      apply_mtc(p_col = "p", grouping = "window_id")
    
    # write per-vcov tables
    write_if_nonempty(pretrend_flat_rows,  file.path(inf_dir, glue::glue("pretrend_flat_{safe_path_component(v_id)}.parquet")))
    write_if_nonempty(pretrend_slope_rows, file.path(inf_dir, glue::glue("pretrend_slope_{safe_path_component(v_id)}.parquet")))
    write_if_nonempty(att_window_rows,     file.path(inf_dir, glue::glue("att_windows_{safe_path_component(v_id)}.parquet")))
    write_if_nonempty(post_slope_rows,     file.path(inf_dir, glue::glue("post_slopes_{safe_path_component(v_id)}.parquet")))
    write_if_nonempty(att_comp_rows,       file.path(inf_dir, glue::glue("att_comparisons_{safe_path_component(v_id)}.parquet")))
    write_if_nonempty(slope_comp_rows,     file.path(inf_dir, glue::glue("slope_comparisons_{safe_path_component(v_id)}.parquet")))
    
    results[[v_id]] <- list(
      pretrend_flat  = pretrend_flat_rows,
      pretrend_slope = pretrend_slope_rows,
      att_windows    = att_window_rows,
      post_slopes    = post_slope_rows,
      att_comparisons = att_comp_rows,
      slope_comparisons = slope_comp_rows
    )
  }
  
  # ── Cohort-group inference (primary vcov only) ──────────────────────────────
  # Loads lean model, rebuilds aggregated coef+vcov by cohort group using
  # the stored primary vcov. For Script 4 only: no MTC here since cohort
  # groups are descriptive stratifications, not simultaneous hypothesis tests.
  
  primary_vcov_id   <- run_spec$primary_vcov_id
  model_file        <- file.path(run_dir, "model.rds")
  
  cg_att_file    <- file.path(inf_dir, "cohortgroup_att_windows_primary.parquet")
  cg_slope_file  <- file.path(inf_dir, "cohortgroup_post_slopes_primary.parquet")
  
  if (!is.na(primary_vcov_id) && file.exists(model_file) && length(cohort_group_specs) > 0) {
    tryCatch({
      model <- readRDS(model_file)
      b     <- coef(model)
      V     <- stats::vcov(model)  # returns stored primary vcov from lean model
      
      cg_agg <- build_agg_eventstudy_cohortgroup_from_bV(
        b                = b,
        V                = V,
        dummy_cols       = dummy_cols,
        cohort_group_defs = cohort_group_specs
      )
      rm(model)
      
      if (length(cg_agg$coef) > 0) {
        z_val      <- stats::qnorm(1 - (1 - ci_level) / 2)
        cg_ses     <- sqrt(diag(cg_agg$vcov))
        
        cg_es_tbl <- cg_agg$groups |>
          dplyr::mutate(estimate = cg_agg$coef, se = cg_ses,
                        ci_lower = estimate - z_val * se,
                        ci_upper = estimate + z_val * se,
                        vcov_id = primary_vcov_id) |>
          add_inference_meta(meta, primary_vcov_id)
        
        # ATT windows by cohort group
        cg_att_rows <- purrr::map_dfr(posttreatment_window_specs, function(ws) {
          purrr::map_dfr(cohort_group_specs, function(cgs) {
            # subset agg_obj to this cohort_group only
            cg_sub <- subset_agg_by_cohort_group(cg_agg, cohort_group_id = cgs$id)
            purrr::map_dfr(dummy_cols, function(dc) {
              att_window(cg_sub, dummy_group = dc, years = ws$years) |>
                dplyr::mutate(window_id = ws$id, cohort_group = cgs$id)
            })
          })
        }) |>
          add_inference_meta(meta, primary_vcov_id)
        
        # Slopes by cohort group
        cg_slope_rows <- purrr::map_dfr(posttreatment_window_specs, function(ws) {
          purrr::map_dfr(cohort_group_specs, function(cgs) {
            cg_sub <- subset_agg_by_cohort_group(cg_agg, cohort_group_id = cgs$id)
            purrr::map_dfr(dummy_cols, function(dc) {
              gls_slope_one_group(cg_sub, dummy_group = dc, years = ws$years) |>
                dplyr::mutate(window_id = ws$id, cohort_group = cgs$id)
            })
          })
        }) |>
          add_inference_meta(meta, primary_vcov_id)
        
        write_if_nonempty(cg_att_rows,   cg_att_file)
        write_if_nonempty(cg_slope_rows, cg_slope_file)
        results$cohortgroup <- list(att_windows = cg_att_rows, post_slopes = cg_slope_rows)
      }
    }, error = function(e) {
      warning(glue::glue("[{run_id}] cohort-group inference failed: {e$message}"))
    })
  }
  
  invisible(results)
}


# ── Helper: slice agg_obj to a single cohort_group ──────────────────────────

subset_agg_by_cohort_group <- function(cg_agg, cohort_group_id) {
  # cg_agg has groups with cohort_group column; we remap to a standard agg_obj
  idx <- which(cg_agg$groups$cohort_group == cohort_group_id)
  list(
    coef   = cg_agg$coef[idx],
    vcov   = cg_agg$vcov[idx, idx, drop = FALSE],
    groups = cg_agg$groups[idx, c("event_time", "dummy_group")]  # drop cohort_group col for inference fns
  )
}


# ── Helper: attach run/vcov metadata to any output tibble ──────────────────

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


# ── Helper: apply Benjamini-Hochberg + Holm corrections ─────────────────────
# Corrections applied within each window (grouped by window_id),
# because tests within a window share the same multiple-comparison family.

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


# ── Helper: only write non-empty tibbles ────────────────────────────────────

write_if_nonempty <- function(tbl, path) {
  if (is.null(tbl) || nrow(tbl) == 0) return(invisible(NULL))
  arrow::write_parquet(tbl, path)
  invisible(path)
}


# ══════════════════════════════════════════════════════════════════════════════
# 5. Process all runs ----
# ══════════════════════════════════════════════════════════════════════════════

message(glue::glue("[{Sys.time()}] Running inference for {length(run_spec_files)} runs..."))

inference_results <- purrr::map(run_spec_files, function(f) {
  result <- tryCatch(
    run_one_inference(
      run_spec_file             = f,
      pretrend_specs            = pretrend_specs,
      posttreatment_window_specs = posttreatment_window_specs,
      cohort_group_specs        = cohort_group_specs,
      comparison_pairs          = comparison_pairs,
      ci_level                  = ci_level,
      skip_existing             = TRUE
    ),
    error = function(e) {
      message(glue::glue("  [ERROR] {basename(dirname(f))}: {e$message}"))
      NULL
    }
  )
  result
})

message(glue::glue("[{Sys.time()}] Done."))


# ══════════════════════════════════════════════════════════════════════════════
# 6. Rebuild merged inference tables ----
# ══════════════════════════════════════════════════════════════════════════════

rebuild_inference_tables <- function(dir_results, write_csv = TRUE) {
  dir_by_run <- file.path(dir_results, "tables", "by_run")
  dir_inf    <- file.path(dir_results, "tables", "inference", "all")
  dir.create(dir_inf, recursive = TRUE, showWarnings = FALSE)
  
  # find all inference parquets — filenames encode vcov_id and test type
  inf_files <- list.files(file.path(dir_by_run), pattern = "\\.parquet$",
                          recursive = TRUE, full.names = TRUE) |>
    grep("inference/", x = _, value = TRUE, fixed = TRUE)
  
  if (length(inf_files) == 0) {
    warning("No inference parquets found.")
    return(invisible(NULL))
  }
  
  # group by filename stem (strips vcov_id suffix)
  file_types <- c(
    "pretrend_flat",             # pretrend_flat_{vcov_id}.parquet
    "pretrend_slope",            # pretrend_slope_{vcov_id}.parquet
    "att_windows",               # att_windows_{vcov_id}.parquet
    "post_slopes",               # post_slopes_{vcov_id}.parquet
    "att_comparisons",           # att_comparisons_{vcov_id}.parquet
    "slope_comparisons",         # slope_comparisons_{vcov_id}.parquet
    "cohortgroup_att_windows_primary",   # flat name, primary vcov only
    "cohortgroup_post_slopes_primary"    # flat name, primary vcov only
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
    purrr::map(file_types, function(ft) {
      ft_files <- inf_files[stringr::str_detect(basename(inf_files), paste0("^", ft))]
      read_and_merge(ft_files)
    }),
    file_types
  )
  
  files <- purrr::imap(merged, function(tbl, nm) write_pair(tbl, nm))
  
  message(glue::glue("Inference tables merged:"))
  purrr::iwalk(merged, function(tbl, nm) {
    if (!is.null(tbl)) message(glue::glue("  {nm}: {nrow(tbl)} rows, {dplyr::n_distinct(tbl$run_id)} runs"))
  })
  
  invisible(list(tables = merged, files = files))
}


all_inference <- rebuild_inference_tables(dir_results = dir_results, write_csv = TRUE)