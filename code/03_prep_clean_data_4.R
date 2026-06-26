

# 5.6 hours in series locally


rm(list = ls())

if(!requireNamespace("here", quietly = TRUE)) {
  install.packages("here")
}
library(here)

source(here::here("utils", "functions.R"))

install_load_packages(c(
  "tidyverse",
  "terra",
  "jsonlite",
  "sf",
  "janitor",
  "mblm",
  "arrow",
  "forcats",
  "tictoc",
  "glue",
  "rlang",
  "tibble",
  "future",
  "furrr",
  "readr",
  "fs",
  "fect",
  "mgcv",
  "h3jsr",
  "data.table"))

set.seed(1234)
gdrive_folder <- "GEE_resilience_v6_operational_ss500_ts50000"

dir_figs <- here::here('figs')
dir_derived <- here::here('data', 'derived', gdrive_folder)
dir_raw <- here::here('data', 'raw')
dir_manual <- here::here('data', 'manual')

dir_ensure(c(dir_figs,
             dir_derived,
             dir_manual,
             dir_raw))



forested_ecoregions <- tibble(
  na_l3name = c(
    "Blue Mountains",
    "Cascades",
    "Coast Range",
    "Eastern Cascades Slopes and Foothills",
    "Klamath Mountains",
    "North Cascades",
    "Straight of Georgia/Puget Lowland",
    "Willamette Valley",
    "California Coastal Sage, Chaparral, and Oak Woodlands",
    "Sierra Nevada",
    "Southern and Baja California Pine-Oak Mountains",
    "Canadian Rockies",
    "Idaho Batholith",
    "Middle Rockies",
    "Columbia Mountains/Northern Rockies",
    "Southern Rockies",
    "Wasatch and Uinta Mountains",
    "Arizona/New Mexico Mountains",
    "Colorado Plateaus"
  ),
  na_l3code = c(
    "6.2.9", "6.2.7", "7.1.8", "6.2.8", "6.2.11", "6.2.5", "7.1.7",
    "7.1.9", "11.1.1", "6.2.12", "11.1.3", "6.2.4", "6.2.15",
    "6.2.10", "6.2.3", "6.2.14", "6.2.13", "13.1.1", "10.1.6"
  ),
  short_name = c(
    "Blue Mtns", "Cascades", "Coast Range", "Eastern Cascades",
    "Klamath Mtns", "North Cascades", "Puget Lowland",
    "Willamette Valley", "Central California Mtns", "Sierra Nevada",
    "Southern California Mtns", "Canadian Rockies", "Idaho Batholith",
    "Middle Rockies", "Northern Rockies", "Southern Rockies",
    "Wasatch Uinta Mtns", "AZ/NM Mtns", "Colorado Plateaus"
  ),
  region = c(
    "Pacific Northwest", "Pacific Northwest", "Pacific Northwest",
    "Pacific Northwest", "Pacific Northwest", "Pacific Northwest",
    "Pacific Northwest", "Pacific Northwest",
    "California", "California", "California",
    "Upper Rockies", "Upper Rockies", "Upper Rockies", "Upper Rockies",
    "Lower Rockies", "Lower Rockies",
    "Southwest", "Southwest"
  ),
  code_name = c(
    "bluemtns", "cascades", "coastrange", "eastcascades",
    "klamathmtns", "northcascades", "pugetlowland",
    "willamettevalley", "centralcaliforniamtns", "sierranevada",
    "southerncaliforniamtns", "canadianrockies", "idahobatholith",
    "middlerockies", "northernrockies", "southernrockies",
    "wasatchuintamtns", "aznmmtns", "coloradoplateaus"
  ),
  us_l3code = c(
    "11", "4", "1", "9", "78", "77", "2", "3", "6", "5", "8",
    "41", "16", "17", "15", "21", "19", "23", "20"
  )
)


# NFG metadata: https://data.fs.usda.gov/geodata/rastergateway/forest_type/conus_forest_type_group_metadata.php
nfg_lookup <- tribble(
  ~nfg, ~nfg_factor, ~nfg_factor_clean, ~nfg_broad_type,
  100, "White/Red/Jack Pine Group",       "white_red_jack_pine",       "conifer",
  120, "Spruce/Fir Group",                "spruce_fir",                "conifer",
  140, "Longleaf/Slash Pine Group",       "longleaf_slash_pine",       "conifer",
  160, "Loblolly/Shortleaf Pine Group",   "loblolly_shortleaf_pine",   "conifer",
  180, "Pinyon/Juniper Group",            "pinyon_juniper",            "conifer",
  200, "Douglas-fir Group",               "douglas_fir",               "conifer",
  220, "Ponderosa Pine Group",            "ponderosa_pine",            "conifer",
  240, "Western White Pine Group",        "western_white_pine",        "conifer",
  260, "Fir/Spruce/Mountain Hemlock Group","fir_spruce_mountain_hemlock","conifer",
  280, "Lodgepole Pine Group",            "lodgepole_pine",            "conifer",
  300, "Hemlock/Sitka Spruce Group",      "hemlock_sitka_spruce",      "conifer",
  320, "Western Larch Group",             "western_larch",             "conifer",
  340, "Redwood Group",                   "redwood",                   "conifer",
  360, "Other Western Softwood Group",    "other_western_softwood",    "conifer",
  370, "California Mixed Conifer Group",  "california_mixed_conifer",  "conifer",
  380, "Exotic Softwoods Group",          "exotic_softwoods",          "conifer",
  400, "Oak/Pine Group",                  "oak_pine",                  "mixed",
  500, "Oak/Hickory Group",               "oak_hickory",               "broadleaf",
  600, "Oak/Gum/Cypress Group",           "oak_gum_cypress",           "broadleaf",
  700, "Elm/Ash/Cottonwood Group",        "elm_ash_cottonwood",        "broadleaf",
  800, "Maple/Beech/Birch Group",         "maple_beech_birch",         "broadleaf",
  900, "Aspen/Birch Group",               "aspen_birch",               "broadleaf",
  910, "Alder/Maple Group",               "alder_maple",               "broadleaf",
  920, "Western Oak Group",               "western_oak",               "broadleaf",
  940, "Tanoak/Laurel Group",             "tanoak_laurel",             "broadleaf",
  950, "Other Western Hardwoods Group",   "other_western_hardwoods",   "broadleaf",
  980, "Tropical Hardwoods Group",        "tropical_hardwoods",        "broadleaf",
  990, "Exotic Hardwoods Group",          "exotic_hardwoods",          "broadleaf"
)



# Remove ecoregions with fewer than 100 fire samples
forested_ecoregions <- forested_ecoregions |> filter(code_name != "coastrange" &
                                code_name != "pugetlowland" & 
                                code_name != "willamettevalley" &
                                code_name != "coloradoplateaus" &
                                code_name != "centralcaliforniamtns")


process_ecoregion <- function(l3, parquet_path, out_dir, utils_path, clean_prep_path) {
  
  source(utils_path)
  source(clean_prep_path)
  
  print(glue("Operating on {l3}"))
  
  d <- arrow::open_dataset(parquet_path, format = "parquet") |>
    dplyr::filter(ecoregion_code_name == l3) |>
    dplyr::collect()
  print("prepping")
  d <- basic_prep_fire(d)
  print("biotic_drought")
  d <- biotic_drought_process(d) #thresholds are set inside function
  
  print("gams")
  d <- compute_gam_stats(
    df = d,
    column_prefix = "rap_tree_",
    reference_time_col = "mock_burn_year",
    baseline = -6,
    gam_method = "REML",
    k = 8,
    min_n = 15,
    parallel = FALSE,
    debug = FALSE
  )
  

  d <- compute_gam_stats(
    df = d,
    column_prefix = "vcf_tree_",
    reference_time_col = "mock_burn_year",
    baseline = -6,
    gam_method = "REML",
    k = 8,
    min_n = 15,
    parallel = FALSE,
    debug = FALSE
  )
  
  print("raw")
  d <- d |>
    compute_avgs(prefix = "rap_tree_") |>
    compute_avgs(prefix = "vcf_tree_") |>
    compute_response_differences_recovery(prefix = "rap_tree_")  |>
    compute_response_differences_recovery(prefix = "vcf_tree_") |>
    compute_response_slopes(prefix = "rap_tree_") |>
    compute_response_slopes(prefix = "vcf_tree_")
  
  print("generating")
  # Final clean
  d <- d |>
    remove_columns_with_prefix(prefixes = c("cbi",
                                            "distance_forest_after"))
  
  # Generate and write datasets
  
  # Short data
  short_fln <- here(out_dir, "parquet_short", glue("dats_short_{l3}.parquet"))
  d_short <- d |>
    remove_columns_with_prefix(prefixes = c("^pdsi_summer_\\d{4}$",
                                            "^pdsi_annual_\\d{4}$",
                                            "^hd_fingerprint_\\d{4}$",
                                            "^biotic_relaxedforestnorm\\d{4}$",
                                            "^biotic6roll_relaxedforestnorm_",
                                            "rap_tree",
                                            "vcf_tree"))   
  
  arrow::write_parquet(x = d_short, sink = short_fln)
  
  # # GPKG
  # gpkg_fln <- here(out_dir, glue("dats_gpkg_{l3}.gpkg"))
  # gpkg <- d_short |>
  #   select(pt_id, lat, long, ecoregion_code_name) |>
  #   st_as_sf(
  #     coords = c("long", "lat"),
  #     crs = 4326,
  #     remove = FALSE
  #   )
  
  # sf::st_write(gpkg, gpkg_fln, append = FALSE)
  
  # Long data
  long_fln <- here(out_dir, "parquet_long", glue("dats_long_{l3}.parquet"))
  
  # Set up for dydid
  d_long <- d |>
    pivot_longer(
      cols = matches("_(19|20)\\d{2}$"),        # ends with a 4-digit year
      names_to = c(".value", "year"),           # .value = "pdsi_summer", "vpd_fall", ...
      names_pattern = "(.*)_(\\d{4})"           # capture "pdsi_summer" and the year
    ) |>
    mutate(year = as.integer(year)) |>
    mutate(year_from_fire_index = year - mock_burn_year,
           fire_timeframe = ifelse(year_from_fire_index < 0, "before_burn", ifelse(year_from_fire_index > 0, "after_burn", "year_of_burn")),
           treated = ifelse(year_from_fire_index >= 0 & fire == 1, 1, 0)) |>
    fect::get.cohort(D = "treated",
                     index = c("pt_id", "year"),
                     start0 = TRUE) |>
    select(-starts_with('gam_'), -starts_with('raw')) |>
    filter(year >= 1992)
  
  d_long[which(is.na(d_long$FirstTreat)),"FirstTreat"] <- 1000  
  
  arrow::write_parquet(x = d_long, sink = long_fln)
  
  flnm_list <- list("short_fln" = short_fln,
                    "long_fln" = long_fln#,
                    # "gpkg_fln" = gpkg_fln
  )
  
  return(flnm_list)
}





# Use functions ----

# Started 5:24pm
parquet_path <- file.path(dir_derived, "merged_dats_raw.parquet")

dir_ensure(c(
  file.path(dir_derived, "parquet_short"),
  file.path(dir_derived, "parquet_long")
))



future::plan(future::multisession, workers = 4)
opts <- furrr::furrr_options(
  seed = TRUE,
  globals = list(
    nfg_lookup = nfg_lookup),
  packages = c(
    "arrow",
    "dplyr",
    "tidyr",
    "tibble",
    "rlang",
    "glue",
    "data.table",
    "mgcv",
    "mblm",
    "fect",
    "here"
  )
)


dats_by_l3 <- furrr::future_map(
  forested_ecoregions$code_name,
  process_ecoregion,
  parquet_path = parquet_path,
  out_dir = dir_derived,
  utils_path = here::here("utils", "functions.R"),
  clean_prep_path = here::here("code", "clean_prep_functions_v4.R"),
  .options = opts
)


# tic()
# dats_by_l3 <- purrr::map(
#   forested_ecoregions$code_name[7],
#   process_ecoregion,
#   parquet_path = parquet_path,
#   out_dir = dir_derived,
#   utils_path = here::here("utils", "functions.R")
# )
# toc()






# Use full dataset to assign analysis groups ----

dir_long  <- here(dir_derived, "parquet_long")
dir_short <- here(dir_derived, "parquet_short")

long_files <- list.files(dir_long, pattern = "\\.parquet$", full.names = TRUE)
short_files <- list.files(dir_short, pattern = "\\.parquet$", full.names = TRUE)

pwalk(
  list(lf = long_files, sf = short_files),
  function(lf, sf) {

    short_dat <- arrow::read_parquet(sf)

    short_dat$h3jsr_5 <- h3jsr::point_to_cell(
      input = as.matrix(short_dat[, c("long", "lat")]),
      res   = 5
    )


    # assign analysis_subset 1-20 within spatial groups for controls,
    # within fire event groups for treated units
    controls <- short_dat |>
      dplyr::filter(fire != 1) |>
      dplyr::group_by(h3jsr_5) |>
      dplyr::mutate(
        analysis_subset = sample(rep(1:20, length.out = dplyr::n()))
      ) |>
      dplyr::ungroup()

    treated <- short_dat |>
      dplyr::filter(fire == 1) |>
      dplyr::group_by(fireid) |>
      dplyr::mutate(
        analysis_subset = sample(rep(1:20, length.out = dplyr::n()))
      ) |>
      dplyr::ungroup()

    short_dat <- dplyr::bind_rows(controls, treated)
    rm(controls, treated)

    forest_groups_over10k <- short_dat |>
      group_by(nfg_factor) |>
      summarize(n = n()) |>
      arrange(desc(n)) |>
      filter(n >= 10000) |>
      pull(nfg_factor)

    short_dat <- short_dat |>
      mutate(nfg_10k = nfg_factor %in% forest_groups_over10k)


    # Join data to long
    subset_lookup <- short_dat |>
      dplyr::select(pt_id, analysis_subset) |>
      dplyr::distinct()

    h3_lookup <- short_dat |>
      dplyr::select(pt_id, h3jsr_5) |>
      dplyr::distinct()

    nfg10k_lookup <- short_dat |>
      dplyr::select(pt_id, nfg_10k) |>
      dplyr::distinct()

    long_dat <- arrow::read_parquet(lf) |>
      dplyr::left_join(h3_lookup,     by = "pt_id") |>
      dplyr::left_join(subset_lookup, by = "pt_id") |>
      dplyr::left_join(nfg10k_lookup, by = "pt_id")

    arrow::write_parquet(long_dat,  lf)
    arrow::write_parquet(short_dat, sf)
  }
)














# 
# # Deprecated gam function ----
# 
# compute_gam_stats <- function(
#     df,
#     column_prefix,
#     reference_time_col,
#     k = 10,
#     min_n = 8,
#     post_min_search_max_year = Inf,
#     gam_method = "REML",
#     parallel = FALSE,
#     n_cores = max(1L, parallel::detectCores() - 1L),
#     debug = FALSE
# ) {
#   stopifnot(is.data.frame(df))
#   if (!requireNamespace("data.table", quietly = TRUE)) stop("Please install 'data.table'.")
#   if (!requireNamespace("mgcv", quietly = TRUE)) stop("Please install 'mgcv'.")
#   if (!requireNamespace("tibble", quietly = TRUE)) stop("Please install 'tibble'.")
#   
#   # ---- validate required column ----
#   if (!reference_time_col %in% names(df)) {
#     stop("reference_time_col not found in df: ", reference_time_col)
#   }
#   
#   # ---- identify time-series cols ----
#   ts_cols <- grep(paste0("^", column_prefix, "\\d{4}$"), names(df), value = TRUE)
#   
#   prefix_tag <- sub("_$", "", column_prefix)
#   
#   nm_pre6_fit <- paste0("gam_", prefix_tag, "_pre6_fit")
#   nm_post_min <- paste0("gam_", prefix_tag, "_post_min")
#   nm_post_min_year <- paste0("gam_", prefix_tag, "_post_min_year")
#   nm_yrs_ref_to_post_min <- paste0("gam_", prefix_tag, "_yrs_ref_to_post_min")
#   nm_diff_pre6_to_post_min <- paste0("gam_", prefix_tag, "_diff_pre6_to_post_min")
#   
#   nm_fit_ref_p10 <- paste0("gam_", prefix_tag, "_fit_ref_p10")
#   nm_fit_ref_p15 <- paste0("gam_", prefix_tag, "_fit_ref_p15")
#   nm_fit_ref_p20 <- paste0("gam_", prefix_tag, "_fit_ref_p20")
#   
#   nm_diff_post_min_to_p10 <- paste0("gam_", prefix_tag, "_diff_post_min_to_p10")
#   nm_diff_post_min_to_p15 <- paste0("gam_", prefix_tag, "_diff_post_min_to_p15")
#   nm_diff_post_min_to_p20 <- paste0("gam_", prefix_tag, "_diff_post_min_to_p20")
#   
#   nm_slope_post_min_to_p10 <- paste0("gam_", prefix_tag, "_slope_post_min_to_p10")
#   nm_slope_post_min_to_p15 <- paste0("gam_", prefix_tag, "_slope_post_min_to_p15")
#   nm_slope_post_min_to_p20 <- paste0("gam_", prefix_tag, "_slope_post_min_to_p20")
#   
#   nm_min_at_end <- paste0("gam_", prefix_tag, "_min_at_end")
#   
#   nm_perc_recov_p10 <- paste0("gam_", prefix_tag, "_perc_recov_p10")
#   nm_perc_recov_p15 <- paste0("gam_", prefix_tag, "_perc_recov_p15")
#   nm_perc_recov_p20 <- paste0("gam_", prefix_tag, "_perc_recov_p20")
#   
#   metric_names_num <- c(
#     nm_pre6_fit, nm_post_min, nm_post_min_year, nm_yrs_ref_to_post_min, nm_diff_pre6_to_post_min,
#     nm_fit_ref_p10, nm_fit_ref_p15, nm_fit_ref_p20,
#     nm_diff_post_min_to_p10, nm_diff_post_min_to_p15, nm_diff_post_min_to_p20,
#     nm_slope_post_min_to_p10, nm_slope_post_min_to_p15, nm_slope_post_min_to_p20,
#     nm_perc_recov_p10, nm_perc_recov_p15, nm_perc_recov_p20
#   )
#   metric_name_logical <- nm_min_at_end
#   
#   # ---- convert once ----
#   dt <- data.table::as.data.table(df)
#   dt[, id__temp := .I]
#   
#   make_stats_na <- function(ids_dt) {
#     out <- data.table::data.table(id__temp = ids_dt$id__temp)
#     out[, (metric_names_num) := NA_real_]
#     out[, (metric_name_logical) := NA]
#     out
#   }
#   
#   # ---- Robust exit: no TS columns found ----
#   # Fill NA metric columns and return (do not error).
#   if (length(ts_cols) == 0) {
#     stats_dt <- make_stats_na(dt[, .(id__temp)])
#     out <- merge(dt, stats_dt, by = "id__temp", all.x = TRUE)
#     out[, id__temp := NULL]
#     return(tibble::as_tibble(out))
#   }
#   
#   # ---- melt ----
#   long_dt <- data.table::melt(
#     dt,
#     id.vars = c("id__temp", reference_time_col),
#     measure.vars = ts_cols,
#     variable.name = "year_col",
#     value.name = "value"
#   )
#   
#   long_dt[, year := as.integer(sub(column_prefix, "", year_col))]
#   long_dt[, ref_year := suppressWarnings(as.integer(get(reference_time_col)))]
#   
#   # Keep only rows with a usable reference year
#   long_dt <- long_dt[!is.na(ref_year)]
#   
#   # ---- Robust exit: no rows survive ref-year filter ----
#   if (nrow(long_dt) == 0L) {
#     stats_dt <- make_stats_na(dt[, .(id__temp)])
#     out <- merge(dt, stats_dt, by = "id__temp", all.x = TRUE)
#     out[, id__temp := NULL]
#     return(tibble::as_tibble(out))
#   }
#   
#   dbg <- new.env(parent = emptyenv())
#   dbg$fail_bad_ref <- 0L
#   dbg$fail_min_n <- 0L
#   dbg$fail_years_lt3 <- 0L
#   dbg$fail_gam <- 0L
#   dbg$ok <- 0L
#   dbg$first_gam_error <- NULL
#   
#   na_row <- function(id) {
#     out <- data.table::data.table(id__temp = id)
#     out[, (metric_names_num) := NA_real_]
#     out[, (metric_name_logical) := NA]
#     out
#   }
#   
#   fit_one <- function(d) {
#     id <- unique(d$id__temp)[1]
#     
#     ref <- unique(d$ref_year)
#     if (length(ref) != 1L || is.na(ref)) {
#       if (isTRUE(debug)) dbg$fail_bad_ref <- dbg$fail_bad_ref + 1L
#       return(na_row(id))
#     }
#     
#     d_fit <- d[!is.na(value) & !is.na(year)]
#     if (nrow(d_fit) < min_n) {
#       if (isTRUE(debug)) dbg$fail_min_n <- dbg$fail_min_n + 1L
#       return(na_row(id))
#     }
#     
#     years_obs <- sort(unique(d_fit$year))
#     if (length(years_obs) < 3) {
#       if (isTRUE(debug)) dbg$fail_years_lt3 <- dbg$fail_years_lt3 + 1L
#       return(na_row(id))
#     }
#     
#     d_fit[, t := year - ref]
#     t_obs <- sort(unique(d_fit$t))
#     
#     k_eff <- min(as.integer(k), max(3L, length(t_obs) - 1L))
#     
#     g <- tryCatch(
#       mgcv::gam(value ~ s(t, k = k_eff), data = d_fit, method = gam_method),
#       error = function(e) {
#         if (isTRUE(debug) && is.null(dbg$first_gam_error)) dbg$first_gam_error <- conditionMessage(e)
#         NULL
#       }
#     )
#     if (is.null(g)) {
#       if (isTRUE(debug)) dbg$fail_gam <- dbg$fail_gam + 1L
#       return(na_row(id))
#     }
#     
#     t_min <- min(t_obs)
#     t_max <- max(t_obs)
#     
#     t_pre6 <- -6L
#     pre6_fit <- if (t_pre6 >= t_min && t_pre6 <= t_max) {
#       as.numeric(stats::predict(g, newdata = data.frame(t = t_pre6), type = "response"))
#     } else {
#       NA_real_
#     }
#     
#     t_cap <- if (is.finite(post_min_search_max_year)) post_min_search_max_year - ref else Inf
#     post_ts <- t_obs[t_obs > 0 & t_obs <= t_cap]
#     
#     post_min_val <- NA_real_
#     post_min_t <- NA_real_
#     
#     if (length(post_ts) > 0) {
#       post_fits <- as.numeric(stats::predict(g, newdata = data.frame(t = post_ts), type = "response"))
#       i_min <- which.min(post_fits)
#       post_min_val <- post_fits[i_min]
#       post_min_t <- post_ts[i_min]
#     }
#     
#     post_min_year <- if (!is.na(post_min_t)) ref + post_min_t else NA_real_
#     yrs_ref_to_min <- if (!is.na(post_min_t)) post_min_t else NA_real_
#     
#     diff_pre6_to_min <- if (!is.na(pre6_fit) && !is.na(post_min_val)) post_min_val - pre6_fit else NA_real_
#     
#     t_p10 <- 10L
#     t_p15 <- 15L
#     t_p20 <- 20L
#     
#     fit_p10 <- if (t_p10 >= t_min && t_p10 <= t_max) {
#       as.numeric(stats::predict(g, newdata = data.frame(t = t_p10), type = "response"))
#     } else {
#       NA_real_
#     }
#     
#     fit_p15 <- if (t_p15 >= t_min && t_p15 <= t_max) {
#       as.numeric(stats::predict(g, newdata = data.frame(t = t_p15), type = "response"))
#     } else {
#       NA_real_
#     }
#     
#     fit_p20 <- if (t_p20 >= t_min && t_p20 <= t_max) {
#       as.numeric(stats::predict(g, newdata = data.frame(t = t_p20), type = "response"))
#     } else {
#       NA_real_
#     }
#     
#     diff_min_to_p10 <- if (!is.na(post_min_val) && !is.na(fit_p10)) fit_p10 - post_min_val else NA_real_
#     diff_min_to_p15 <- if (!is.na(post_min_val) && !is.na(fit_p15)) fit_p15 - post_min_val else NA_real_
#     diff_min_to_p20 <- if (!is.na(post_min_val) && !is.na(fit_p20)) fit_p20 - post_min_val else NA_real_
#     
#     slope_to_p10 <- if (!is.na(post_min_t) && !is.na(post_min_val) && !is.na(fit_p10) && (t_p10 != post_min_t)) {
#       (fit_p10 - post_min_val) / (t_p10 - post_min_t)
#     } else {
#       NA_real_
#     }
#     
#     slope_to_p15 <- if (!is.na(post_min_t) && !is.na(post_min_val) && !is.na(fit_p15) && (t_p15 != post_min_t)) {
#       (fit_p15 - post_min_val) / (t_p15 - post_min_t)
#     } else {
#       NA_real_
#     }
#     
#     slope_to_p20 <- if (!is.na(post_min_t) && !is.na(post_min_val) && !is.na(fit_p20) && (t_p20 != post_min_t)) {
#       (fit_p20 - post_min_val) / (t_p20 - post_min_t)
#     } else {
#       NA_real_
#     }
#     
#     min_at_end <- if (!is.na(post_min_t) && length(post_ts) > 0) {
#       isTRUE(post_min_t == max(post_ts))
#     } else {
#       NA
#     }
#     
#     denom <- abs(diff_pre6_to_min)
#     perc_recov_p10 <- if (!is.na(diff_min_to_p10) && !is.na(denom) && denom > 0) (diff_min_to_p10 / denom) * 100 else NA_real_
#     perc_recov_p15 <- if (!is.na(diff_min_to_p15) && !is.na(denom) && denom > 0) (diff_min_to_p15 / denom) * 100 else NA_real_
#     perc_recov_p20 <- if (!is.na(diff_min_to_p20) && !is.na(denom) && denom > 0) (diff_min_to_p20 / denom) * 100 else NA_real_
#     
#     if (isTRUE(debug)) dbg$ok <- dbg$ok + 1L
#     
#     out <- data.table::data.table(
#       id__temp = id,
#       pre6_fit = pre6_fit,
#       post_min_val = post_min_val,
#       post_min_year = post_min_year,
#       yrs_ref_to_min = yrs_ref_to_min,
#       diff_pre6_to_min = diff_pre6_to_min,
#       fit_p10 = fit_p10,
#       fit_p15 = fit_p15,
#       fit_p20 = fit_p20,
#       diff_min_to_p10 = diff_min_to_p10,
#       diff_min_to_p15 = diff_min_to_p15,
#       diff_min_to_p20 = diff_min_to_p20,
#       slope_to_p10 = slope_to_p10,
#       slope_to_p15 = slope_to_p15,
#       slope_to_p20 = slope_to_p20,
#       min_at_end = min_at_end,
#       perc_recov_p10 = perc_recov_p10,
#       perc_recov_p15 = perc_recov_p15,
#       perc_recov_p20 = perc_recov_p20
#     )
#     
#     data.table::setnames(
#       out,
#       old = c(
#         "pre6_fit",
#         "post_min_val",
#         "post_min_year",
#         "yrs_ref_to_min",
#         "diff_pre6_to_min",
#         "fit_p10",
#         "fit_p15",
#         "fit_p20",
#         "diff_min_to_p10",
#         "diff_min_to_p15",
#         "diff_min_to_p20",
#         "slope_to_p10",
#         "slope_to_p15",
#         "slope_to_p20",
#         "min_at_end",
#         "perc_recov_p10",
#         "perc_recov_p15",
#         "perc_recov_p20"
#       ),
#       new = c(
#         nm_pre6_fit,
#         nm_post_min,
#         nm_post_min_year,
#         nm_yrs_ref_to_post_min,
#         nm_diff_pre6_to_post_min,
#         nm_fit_ref_p10,
#         nm_fit_ref_p15,
#         nm_fit_ref_p20,
#         nm_diff_post_min_to_p10,
#         nm_diff_post_min_to_p15,
#         nm_diff_post_min_to_p20,
#         nm_slope_post_min_to_p10,
#         nm_slope_post_min_to_p15,
#         nm_slope_post_min_to_p20,
#         nm_min_at_end,
#         nm_perc_recov_p10,
#         nm_perc_recov_p15,
#         nm_perc_recov_p20
#       )
#     )
#     
#     out
#   }
#   
#   # Split by id__temp. This can still be empty in edge cases; handle it.
#   split_list <- split(long_dt, long_dt$id__temp)
#   
#   if (length(split_list) == 0L) {
#     stats_dt <- make_stats_na(dt[, .(id__temp)])
#     out <- merge(dt, stats_dt, by = "id__temp", all.x = TRUE)
#     out[, id__temp := NULL]
#     return(tibble::as_tibble(out))
#   }
#   
#   if (isTRUE(parallel)) {
#     if (!requireNamespace("future.apply", quietly = TRUE)) stop("Please install 'future.apply' for parallel=TRUE.")
#     if (!requireNamespace("future", quietly = TRUE)) stop("Please install 'future' for parallel=TRUE.")
#     old_plan <- future::plan()
#     on.exit(future::plan(old_plan), add = TRUE)
#     
#     future::plan(future::multisession, workers = n_cores)
#     res_list <- future.apply::future_lapply(split_list, fit_one)
#   } else {
#     res_list <- lapply(split_list, fit_one)
#   }
#   
#   if (length(res_list) == 0L) {
#     stats_dt <- make_stats_na(dt[, .(id__temp)])
#   } else {
#     stats_dt <- data.table::rbindlist(res_list, use.names = TRUE, fill = TRUE)
#     
#     # If something went very wrong, force a safe NA stats_dt
#     if (!("id__temp" %in% names(stats_dt))) {
#       stats_dt <- make_stats_na(dt[, .(id__temp)])
#     }
#   }
#   
#   out <- merge(dt, stats_dt, by = "id__temp", all.x = TRUE)
#   out[, id__temp := NULL]
#   
#   if (isTRUE(debug)) {
#     message(
#       "compute_gam_stats debug counts:\n",
#       "ok: ", dbg$ok, "\n",
#       "fail_bad_ref: ", dbg$fail_bad_ref, "\n",
#       "fail_min_n: ", dbg$fail_min_n, "\n",
#       "fail_years_lt3: ", dbg$fail_years_lt3, "\n",
#       "fail_gam: ", dbg$fail_gam
#     )
#     if (!is.null(dbg$first_gam_error)) message("First GAM error: ", dbg$first_gam_error)
#   }
#   
#   tibble::as_tibble(out)
# }



# # Deprecated drought/biotic processing ----
# 
# # Get biotic & drought sums and means in time before/after burn
# 
# # helper: robust row/column extraction for data.table OR data.frame
# .row_vals <- function(df, i, cols) {
#   if (data.table::is.data.table(df)) {
#     unlist(df[i, ..cols], use.names = FALSE)
#   } else {
#     unlist(df[i, cols, drop = FALSE], use.names = FALSE)
#   }
# }
# 
# transform_annual_to_priorafter <- function(
#     df,
#     column_pattern,
#     summary_fn_name,
#     nyears,
#     before_inclusive_year_of = FALSE
# ) {
#   # Validate the summary function
#   if (!exists(summary_fn_name, mode = "function")) {
#     stop("Summary function '", summary_fn_name, "' is not a valid function.")
#   }
#   summary_fn <- match.fun(summary_fn_name)
#   
#   matched_cols <- grep(column_pattern, names(df), value = TRUE)
#   if (length(matched_cols) == 0) return(df)
#   
#   base_name <- sub("^(.*?)_\\d{4}$", "\\1", matched_cols[1])
#   
#   col_years <- as.integer(sub(".*_(\\d{4})$", "\\1", matched_cols))
#   col_year_map <- setNames(matched_cols, col_years)
#   
#   n <- nrow(df)
#   prior_vals <- rep(NA_real_, n)
#   after_vals <- rep(NA_real_, n)
#   
#   # window + suffix logic
#   if (isTRUE(before_inclusive_year_of)) {
#     prior_suffix <- "yot"  # inclusive year-of
#     prior_years_fn <- function(byr) (byr - (nyears - 1)):byr
#   } else {
#     prior_suffix <- "yof"  # as now (exclusive year-of)
#     prior_years_fn <- function(byr) (byr - nyears):(byr - 1)
#   }
#   
#   for (i in seq_len(n)) {
#     byr <- df$mock_burn_year[i]
#     if (!is.na(byr)) {
#       yrs_prior <- prior_years_fn(byr)
#       if (all(yrs_prior %in% col_years)) {
#         cols <- col_year_map[as.character(yrs_prior)]
#         vals <- .row_vals(df, i, cols)
#         if (!any(is.na(vals))) prior_vals[i] <- summary_fn(vals, na.rm = TRUE)
#       }
#       
#       yrs_after <- (byr + 1):(byr + nyears)
#       if (all(yrs_after %in% col_years)) {
#         cols <- col_year_map[as.character(yrs_after)]
#         vals <- .row_vals(df, i, cols)
#         if (!any(is.na(vals))) after_vals[i] <- summary_fn(vals, na.rm = TRUE)
#       }
#     }
#   }
#   
#   fn_suffix <- tolower(summary_fn_name)
#   prior_col <- paste0(base_name, "_", nyears, "_yrs_prior_", fn_suffix, "_", prior_suffix)
#   after_col <- paste0(base_name, "_", nyears, "_yrs_after_", fn_suffix)
#   
#   df[[prior_col]] <- prior_vals
#   df[[after_col]] <- after_vals
#   
#   df
# }
# 
# 
# n_beyond_threshold_x_yrs <- function(
#     df,
#     column_pattern,
#     nyears,
#     threshold,
#     threshold_nm,
#     before_inclusive_year_of = FALSE
# ) {
#   threshold_fn <- tryCatch(
#     eval(parse(text = paste0("function(x) x", threshold))),
#     error = function(e) stop("Invalid threshold expression: ", threshold)
#   )
#   
#   matched_cols <- grep(column_pattern, names(df), value = TRUE)
#   if (length(matched_cols) == 0) return(df)
#   
#   base_name <- sub("^(.*?)_\\d{4}$", "\\1", matched_cols[1])
#   col_years <- as.integer(sub(".*_(\\d{4})$", "\\1", matched_cols))
#   col_year_map <- setNames(matched_cols, col_years)
#   
#   n <- nrow(df)
#   prior_counts <- integer(n)
#   after_counts <- integer(n)
#   
#   # window + suffix logic
#   if (isTRUE(before_inclusive_year_of)) {
#     prior_suffix <- "yot"
#     prior_years_fn <- function(byr) (byr - (nyears - 1)):byr
#   } else {
#     prior_suffix <- "yof"
#     prior_years_fn <- function(byr) (byr - nyears):(byr - 1)
#   }
#   
#   for (i in seq_len(n)) {
#     byr <- df$mock_burn_year[i]
#     if (!is.na(byr)) {
#       yrs_prior <- prior_years_fn(byr)
#       if (all(yrs_prior %in% col_years)) {
#         cols <- col_year_map[as.character(yrs_prior)]
#         vals <- .row_vals(df, i, cols)
#         prior_counts[i] <- sum(threshold_fn(vals), na.rm = TRUE)
#       } else {
#         prior_counts[i] <- NA_integer_
#       }
#       
#       yrs_after <- (byr + 1):(byr + nyears)
#       if (all(yrs_after %in% col_years)) {
#         cols <- col_year_map[as.character(yrs_after)]
#         vals <- .row_vals(df, i, cols)
#         after_counts[i] <- sum(threshold_fn(vals), na.rm = TRUE)
#       } else {
#         after_counts[i] <- NA_integer_
#       }
#     } else {
#       prior_counts[i] <- NA_integer_
#       after_counts[i] <- NA_integer_
#     }
#   }
#   
#   prior_col <- paste0(base_name, "_", nyears, "_yrs_prior_threshold_", threshold_nm, "_", prior_suffix)
#   after_col <- paste0(base_name, "_", nyears, "_yrs_after_threshold_", threshold_nm)
#   
#   df[[prior_col]] <- prior_counts
#   df[[after_col]] <- after_counts
#   
#   df
# }






# 
# # Compare old version to new version
# 
# dir_short <- here(dir_derived, "parquet_short")
# 
# old <- read_parquet(here(dir_short, "deprecated", "dats_short_southerncaliforniamtns.parquet"))
# new <- read_parquet(here(dir_short, "dats_short_southerncaliforniamtns.parquet"))
# glimpse(old)
# glimpse(new)
# 
# 
# diff <- old |> filter(! pt_id %in% new$pt_id)
# hist(diff$mock_burn_year)
# 
# hist(old$mock_burn_year)
# hist(new$mock_burn_year)
# 
# 
# nrow(old |> filter(is.na(raw_post_fire_rap_tree_difference_n6_3_yr)))
# nrow(new |> filter(is.na(raw_post_fire_rap_tree_difference_n6_3_yr)))
# 
# nrow(old |> filter(is.na(gam_rap_tree_diff_pre6_to_post_min)))
# nrow(new |> filter(is.na(gam_rap_tree_diff_pre6_to_post_min)))
# 
# min_na <- old[which.min(rowSums(is.na(old))), ] |> pull(pt_id)
# 
# old_test_unit <- old |> filter(pt_id == min_na)
# new_test_unit <- new |> filter(pt_id == min_na)
# 
# 
# 
# raw <- arrow::open_dataset(parquet_path, format = "parquet") |>
#   dplyr::filter(ecoregion_code_name == "southerncaliforniamtns") |>
#   dplyr::collect()
# 
# raw_test_unit <- raw |> filter(pt_id == min_na)
# 
# 
# 
# -271 + -464 + -169
# 
# 217 + -307 + -381 + -271 + -464 + -169
# 
# 217 + -307 + -381 + -271 + -464
# 
# -307 + -381 + -271 + -464 + -169
