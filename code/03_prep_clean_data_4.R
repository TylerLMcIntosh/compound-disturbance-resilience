
# 5.6 hours in series locally


rm(list = ls())

cyverse = FALSE

if(!requireNamespace("here", quietly = TRUE)) {
  install.packages("here")
}
library(here)

source(here::here("utils", "functions.R"))

install_load_packages(c(
  "tidyverse",
  "terra",
  "mapview",
  "jsonlite",
  "sf",
  "janitor",
  "googledrive",
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
  "data.table"))

set.seed(1234)
gdrive_folder <- "GEE_resilience_v6_operational_ss500_ts50000"
#run <- "v6_ss1000_ts100000"
service_account_file <- here('config', 'secrets', 'tymc5571-utils-project-692f27c034dd.json')

if(cyverse) {
  dir_figs <- file.path("~/data-store/data/iplant/home/shared/earthlab/macrosystems/tlm", run, "figs")
  dir_derived <- file.path("~/data-store/data/iplant/home/shared/earthlab/macrosystems/tlm", run, "data/derived")
  dir_parquet <- here::here(dir_derived, 'parquet')
  dir_manual <- file.path("~/data-store/data/iplant/home/shared/earthlab/macrosystems/tlm", run, "data/manual")
  dir_raw <- file.path("~/data-store/data/iplant/home/shared/earthlab/macrosystems/tlm", run, "data/raw")
} else {
  dir_figs <- here::here('figs')
  dir_derived <- here::here('data', 'derived', gdrive_folder)
  dir_parquet <- here::here(dir_derived, 'parquet')
  dir_raw <- here::here('data', 'raw')
  dir_manual <- here::here('data', 'manual')
}

dir_ensure(c(dir_figs,
             dir_derived,
             dir_parquet,
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

# Remove ecoregions with fewer than 100 fire samples
forested_ecoregions <- forested_ecoregions |> filter(code_name != "coastrange" &
                                code_name != "pugetlowland" & 
                                code_name != "willamettevalley" &
                                code_name != "coloradoplateaus" &
                                code_name != "centralcaliforniamtns")


process_ecoregion <- function(l3, parquet_path, out_dir, utils_path) {
  
  source(utils_path)
  
  print(glue("Operating on {l3}"))
  
  d <- arrow::open_dataset(parquet_path, format = "parquet") |>
    dplyr::filter(ecoregion_code_name == l3) |>
    dplyr::collect()
  print("prepping")
  d <- basic_prep_fire(d)
  print("biotic_drought")
  d <- biotic_drought_process(d,
                              nyears = 5,
                              before_inclusive_year_of = TRUE) #thresholds are set inside function
  
  print("gams")
  d <- compute_gam_stats(
    df = d,
    column_prefix = "rap_tree_",
    reference_time_col = "burn_year",
    k = 8,
    min_n = 15,
    gam_method = "REML",
    parallel = FALSE,
    debug = FALSE
  )
  
  d <- compute_gam_stats(
    df = d,
    column_prefix = "vcf_tree_",
    reference_time_col = "burn_year",
    k = 8,
    min_n = 15,
    gam_method = "REML",
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
  short_fln <- here(out_dir, glue("dats_short_{l3}.parquet"))
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
  long_fln <- here(out_dir, glue("dats_long_{l3}.parquet"))
  
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
    select(-starts_with('gam_'), -starts_with('raw'))
  
  d_long[which(is.na(d_long$FirstTreat)),"FirstTreat"] <- 1000  
  
  arrow::write_parquet(x = d_long, sink = long_fln)
  
  flnm_list <- list("short_fln" = short_fln,
                    "long_fln" = long_fln#,
                    # "gpkg_fln" = gpkg_fln
  )
  
  return(flnm_list)
}



basic_prep_fire <- function(dats) {
  
  names(dats) <-
    gsub("(_yr_|_year_)(?!year$)", "_", names(dats), perl = TRUE)
  
  na_check_cols <- setdiff(names(dats), c('fireid', 'fire_year'))
  dats <- dats |>
    #select(-rand) |>
    tidyr::drop_na(dplyr::all_of(na_check_cols))
  
  if (nrow(dats) == 0) {
    dats$burn_year <- integer()
    dats$mock_burn_year <- integer()
    dats$burn_year_cbi <- numeric()
    dats$burn_sev <- factor(character(), levels = c("Never","Negligible","Low/Moderate","High"))
    return(dats)
  }
  
  # remove any rows with values that are - but not -1 in CBI; shouldn't be any
  dats <- dats |>
    filter(
      !if_any(
        starts_with("cbi_"),
        ~ .x < 0 & .x != -1
      )
    )
  
  cbi_cols <- grep("^cbi_\\d{4}$", names(dats), value = TRUE)
  
  # build the column name each row should use
  target_col <- paste0("cbi_", dats$fire_year)
  
  # extract by row + matched column index
  m <- match(target_col, cbi_cols)
  
  dats$burn_cbi <- NA_real_
  ok <- !is.na(m)
  
  cbi_mat <- as.matrix(dats[, ..cbi_cols, drop = FALSE])
  dats$burn_cbi[ok] <- cbi_mat[cbind(which(ok), m[ok])]
  
  dats <- dats |> rename(burn_year = fire_year)
  
  # NOT USING AT THE MOMENT
  # 
  # # GET THE OTHER CBI CONTEXT DATA
  # 
  # # Function to extract burn-year-aligned values from year-tagged columns
  # extract_burn_year_value <- function(df, burn_year, prefix) {
  #   # Get matching columns like 'cbi_1kmmwmedian_2008'
  #   pattern <- paste0("^", prefix, "_\\d{4}$")
  #   cols <- grep(pattern, names(df), value = TRUE)
  #   
  #   # Extract years from column names
  #   years <- as.integer(sub(".*_(\\d{4})$", "\\1", cols))
  #   
  #   # Create value matrix (rows = obs, cols = years)
  #   value_mat <- as.data.frame(df[cols], stringsAsFactors = FALSE)
  #   
  #   # Create lookup index: for each row, match burn_year to column year
  #   col_index <- match(burn_year, years)
  #   
  #   # Use cbind(row, col) indexing to extract values
  #   result <- mapply(function(row, col) {
  #     if (!is.na(col)) value_mat[row, col] else NA
  #   }, row = seq_along(burn_year), col = col_index)
  #   
  #   return(result)
  # }
  # 
  # # Apply the function to create the new columns
  # dats_no_reburn$burn_cbi_500mmwmean <- extract_burn_year_value(
  #   dats_no_reburn,
  #   burn_year = dats_no_reburn$burn_year,
  #   prefix = "cbi_500mmwmean"
  # )
  # 
  # dats_no_reburn$burn_cbi_500mmwstd <- extract_burn_year_value(
  #   dats_no_reburn,
  #   burn_year = dats_no_reburn$burn_year,
  #   prefix = "cbi_500mmwstd"
  # )
  # 
  # dats_no_reburn$burn_distance_forest_after <- extract_burn_year_value(
  #   dats_no_reburn,
  #   burn_year = dats_no_reburn$burn_year,
  #   prefix = "distance_forest_after"
  # )
  
  # Generate mock_burn_year
  dats$mock_burn_year <- dats$burn_year
  
  na_idx <- which(is.na(dats$mock_burn_year))
  dats$mock_burn_year[na_idx] <- sample(1984:2020, length(na_idx), replace = TRUE)
  
  
  dats <- dats |>
    mutate(
      across(
        c(starts_with("vcf_tree_SD"), starts_with("pdsi_")),
        ~ .x * 0.01
      )
    )
  
  # NFG metadata: https://data.fs.usda.gov/geodata/rastergateway/forest_type/conus_forest_type_group_metadata.php
  nfg_lookup <- c(
    "100" = "White/Red/Jack Pine Group",
    "120" = "Spruce/Fir Group",
    "140" = "Longleaf/Slash Pine Group",
    "160" = "Loblolly/Shortleaf Pine Group",
    "180" = "Pinyon/Juniper Group",
    "200" = "Douglas-fir Group",
    "220" = "Ponderosa Pine Group",
    "240" = "Western White Pine Group",
    "260" = "Fir/Spruce/Mountain Hemlock Group",
    "280" = "Lodgepole Pine Group",
    "300" = "Hemlock/Sitka Spruce Group",
    "320" = "Western Larch Group",
    "340" = "Redwood Group",
    "360" = "Other Western Softwood Group",
    "370" = "California Mixed Conifer Group",
    "380" = "Exotic Softwoods Group",
    "400" = "Oak/Pine Group",
    "500" = "Oak/Hickory Group",
    "600" = "Oak/Gum/Cypress Group",
    "700" = "Elm/Ash/Cottonwood Group",
    "800" = "Maple/Beech/Birch Group",
    "900" = "Aspen/Birch Group",
    "910" = "Alder/Maple Group",
    "920" = "Western Oak Group",
    "940" = "Tanoak/Laurel Group",
    "950" = "Other Western Hardwoods Group",
    "980" = "Tropical Hardwoods Group",
    "990" = "Exotic Hardwoods Group"
  )
  
  dats <- dats |>
    mutate(
      nfg_factor = factor(nfg_lookup[as.character(nfg)])
    )
  
  dats <- dats %>%
    mutate(
      burn_sev = case_when(
        fire == 1 & burn_cbi >= 2.25 ~ "High",
        fire == 1 & burn_cbi < 0.1 ~ "Negligible",
        #fire == 1 & burn_cbi >= 0.1 & burn_cbi < 2.25 ~ "Low/Moderate",
        fire == 1 & burn_cbi >= 0.1 & burn_cbi < 1.25  ~ "Low",
        fire == 1 & burn_cbi >= 1.25 & burn_cbi < 2.25 ~ "Moderate",
        TRUE ~ "Never"
      )
    )
  
  #Clean up factors & years as integers 
  dats <- dats |>
    mutate(frg_reclass_desc = case_when(frg_reclass == 1 ~ "<35 return, low/mixed",
                                        frg_reclass == 2 ~ "<35 return, replacement",
                                        frg_reclass == 3 ~ "35-200 return, low/mixed",
                                        frg_reclass == 4 ~ "35-200 return replacement",
                                        frg_reclass == 5 ~ ">200 return any"),
           frg_reclass_desc = as.factor(frg_reclass_desc),
           burn_sev = factor(
             burn_sev,
             levels = c("Never", "Negligible", "Low", "Moderate", "High")
           ),
           burn_year = as.integer(burn_year),
           mock_burn_year = as.integer(mock_burn_year)
    )
  
  return(dats)
}



# Get biotic & drought sums and means in time before/after burn

# helper: robust row/column extraction for data.table OR data.frame
.row_vals <- function(df, i, cols) {
  if (data.table::is.data.table(df)) {
    unlist(df[i, ..cols], use.names = FALSE)
  } else {
    unlist(df[i, cols, drop = FALSE], use.names = FALSE)
  }
}

transform_annual_to_priorafter <- function(
    df,
    column_pattern,
    summary_fn_name,
    nyears,
    before_inclusive_year_of = FALSE
) {
  # Validate the summary function
  if (!exists(summary_fn_name, mode = "function")) {
    stop("Summary function '", summary_fn_name, "' is not a valid function.")
  }
  summary_fn <- match.fun(summary_fn_name)
  
  matched_cols <- grep(column_pattern, names(df), value = TRUE)
  if (length(matched_cols) == 0) return(df)
  
  base_name <- sub("^(.*?)_\\d{4}$", "\\1", matched_cols[1])
  
  col_years <- as.integer(sub(".*_(\\d{4})$", "\\1", matched_cols))
  col_year_map <- setNames(matched_cols, col_years)
  
  n <- nrow(df)
  prior_vals <- rep(NA_real_, n)
  after_vals <- rep(NA_real_, n)
  
  # window + suffix logic
  if (isTRUE(before_inclusive_year_of)) {
    prior_suffix <- "yot"  # inclusive year-of
    prior_years_fn <- function(byr) (byr - (nyears - 1)):byr
  } else {
    prior_suffix <- "yof"  # as now (exclusive year-of)
    prior_years_fn <- function(byr) (byr - nyears):(byr - 1)
  }
  
  for (i in seq_len(n)) {
    byr <- df$mock_burn_year[i]
    if (!is.na(byr)) {
      yrs_prior <- prior_years_fn(byr)
      if (all(yrs_prior %in% col_years)) {
        cols <- col_year_map[as.character(yrs_prior)]
        vals <- .row_vals(df, i, cols)
        if (!any(is.na(vals))) prior_vals[i] <- summary_fn(vals, na.rm = TRUE)
      }
      
      yrs_after <- (byr + 1):(byr + nyears)
      if (all(yrs_after %in% col_years)) {
        cols <- col_year_map[as.character(yrs_after)]
        vals <- .row_vals(df, i, cols)
        if (!any(is.na(vals))) after_vals[i] <- summary_fn(vals, na.rm = TRUE)
      }
    }
  }
  
  fn_suffix <- tolower(summary_fn_name)
  prior_col <- paste0(base_name, "_", nyears, "_yrs_prior_", fn_suffix, "_", prior_suffix)
  after_col <- paste0(base_name, "_", nyears, "_yrs_after_", fn_suffix)
  
  df[[prior_col]] <- prior_vals
  df[[after_col]] <- after_vals
  
  df
}


n_beyond_threshold_x_yrs <- function(
    df,
    column_pattern,
    nyears,
    threshold,
    threshold_nm,
    before_inclusive_year_of = FALSE
) {
  threshold_fn <- tryCatch(
    eval(parse(text = paste0("function(x) x", threshold))),
    error = function(e) stop("Invalid threshold expression: ", threshold)
  )
  
  matched_cols <- grep(column_pattern, names(df), value = TRUE)
  if (length(matched_cols) == 0) return(df)
  
  base_name <- sub("^(.*?)_\\d{4}$", "\\1", matched_cols[1])
  col_years <- as.integer(sub(".*_(\\d{4})$", "\\1", matched_cols))
  col_year_map <- setNames(matched_cols, col_years)
  
  n <- nrow(df)
  prior_counts <- integer(n)
  after_counts <- integer(n)
  
  # window + suffix logic
  if (isTRUE(before_inclusive_year_of)) {
    prior_suffix <- "yot"
    prior_years_fn <- function(byr) (byr - (nyears - 1)):byr
  } else {
    prior_suffix <- "yof"
    prior_years_fn <- function(byr) (byr - nyears):(byr - 1)
  }
  
  for (i in seq_len(n)) {
    byr <- df$mock_burn_year[i]
    if (!is.na(byr)) {
      yrs_prior <- prior_years_fn(byr)
      if (all(yrs_prior %in% col_years)) {
        cols <- col_year_map[as.character(yrs_prior)]
        vals <- .row_vals(df, i, cols)
        prior_counts[i] <- sum(threshold_fn(vals), na.rm = TRUE)
      } else {
        prior_counts[i] <- NA_integer_
      }
      
      yrs_after <- (byr + 1):(byr + nyears)
      if (all(yrs_after %in% col_years)) {
        cols <- col_year_map[as.character(yrs_after)]
        vals <- .row_vals(df, i, cols)
        after_counts[i] <- sum(threshold_fn(vals), na.rm = TRUE)
      } else {
        after_counts[i] <- NA_integer_
      }
    } else {
      prior_counts[i] <- NA_integer_
      after_counts[i] <- NA_integer_
    }
  }
  
  prior_col <- paste0(base_name, "_", nyears, "_yrs_prior_threshold_", threshold_nm, "_", prior_suffix)
  after_col <- paste0(base_name, "_", nyears, "_yrs_after_threshold_", threshold_nm)
  
  df[[prior_col]] <- prior_counts
  df[[after_col]] <- after_counts
  
  df
}


biotic_drought_process <- function(dats,
                                   nyears,
                                   before_inclusive_year_of = FALSE) {
  
  dats <- transform_annual_to_priorafter(
    df = dats,
    column_pattern = "^biotic_relaxedforestnorm_\\d{4}$",
    summary_fn_name = "sum",
    nyears = nyears,
    before_inclusive_year_of = before_inclusive_year_of
  )
  
  dats <- transform_annual_to_priorafter(
    df = dats,
    column_pattern = "^pdsi_annual_\\d{4}$",
    summary_fn_name = "sum",
    nyears = nyears,
    before_inclusive_year_of = before_inclusive_year_of
  )
  
  dats <- transform_annual_to_priorafter(
    df = dats,
    column_pattern = "^pdsi_annual_\\d{4}$",
    summary_fn_name = "mean",
    nyears = nyears,
    before_inclusive_year_of = before_inclusive_year_of
  )
  
  dats <- n_beyond_threshold_x_yrs(
    df = dats,
    column_pattern = "^pdsi_annual_\\d{4}$",
    nyears = nyears,
    threshold = "<=-4",
    threshold_nm = "n4",
    before_inclusive_year_of = before_inclusive_year_of
  )
  
  dats <- n_beyond_threshold_x_yrs(
    df = dats,
    column_pattern = "^pdsi_annual_\\d{4}$",
    nyears = nyears,
    threshold = "<=-3",
    threshold_nm = "n3",
    before_inclusive_year_of = before_inclusive_year_of
  )
  
  dats <- n_beyond_threshold_x_yrs(
    df = dats,
    column_pattern = "^hd_fingerprint_\\d{4}$",
    nyears = nyears,
    threshold = ">=6",
    threshold_nm = "6",
    before_inclusive_year_of = before_inclusive_year_of
  )
  
  dats <- dats |>
    extract_year_of_value(column_pattern = "^pdsi_annual_\\d{4}$",
                          year_col = "burn_year") |>
    extract_year_of_value(column_pattern = "^hd_fingerprint_\\d{4}$",
                          year_col = "burn_year")
  
  dats
}



# PROCESS GAMS

compute_gam_stats <- function(
    df,
    column_prefix,
    reference_time_col,
    k = 10,
    min_n = 8,
    post_min_search_max_year = Inf,
    gam_method = "REML",
    parallel = FALSE,
    n_cores = max(1L, parallel::detectCores() - 1L),
    debug = FALSE
) {
  stopifnot(is.data.frame(df))
  if (!requireNamespace("data.table", quietly = TRUE)) stop("Please install 'data.table'.")
  if (!requireNamespace("mgcv", quietly = TRUE)) stop("Please install 'mgcv'.")
  if (!requireNamespace("tibble", quietly = TRUE)) stop("Please install 'tibble'.")
  
  # ---- validate required column ----
  if (!reference_time_col %in% names(df)) {
    stop("reference_time_col not found in df: ", reference_time_col)
  }
  
  # ---- identify time-series cols ----
  ts_cols <- grep(paste0("^", column_prefix, "\\d{4}$"), names(df), value = TRUE)
  
  prefix_tag <- sub("_$", "", column_prefix)
  
  nm_pre6_fit <- paste0("gam_", prefix_tag, "_pre6_fit")
  nm_post_min <- paste0("gam_", prefix_tag, "_post_min")
  nm_post_min_year <- paste0("gam_", prefix_tag, "_post_min_year")
  nm_yrs_ref_to_post_min <- paste0("gam_", prefix_tag, "_yrs_ref_to_post_min")
  nm_diff_pre6_to_post_min <- paste0("gam_", prefix_tag, "_diff_pre6_to_post_min")
  
  nm_fit_ref_p10 <- paste0("gam_", prefix_tag, "_fit_ref_p10")
  nm_fit_ref_p15 <- paste0("gam_", prefix_tag, "_fit_ref_p15")
  nm_fit_ref_p20 <- paste0("gam_", prefix_tag, "_fit_ref_p20")
  
  nm_diff_post_min_to_p10 <- paste0("gam_", prefix_tag, "_diff_post_min_to_p10")
  nm_diff_post_min_to_p15 <- paste0("gam_", prefix_tag, "_diff_post_min_to_p15")
  nm_diff_post_min_to_p20 <- paste0("gam_", prefix_tag, "_diff_post_min_to_p20")
  
  nm_slope_post_min_to_p10 <- paste0("gam_", prefix_tag, "_slope_post_min_to_p10")
  nm_slope_post_min_to_p15 <- paste0("gam_", prefix_tag, "_slope_post_min_to_p15")
  nm_slope_post_min_to_p20 <- paste0("gam_", prefix_tag, "_slope_post_min_to_p20")
  
  nm_min_at_end <- paste0("gam_", prefix_tag, "_min_at_end")
  
  nm_perc_recov_p10 <- paste0("gam_", prefix_tag, "_perc_recov_p10")
  nm_perc_recov_p15 <- paste0("gam_", prefix_tag, "_perc_recov_p15")
  nm_perc_recov_p20 <- paste0("gam_", prefix_tag, "_perc_recov_p20")
  
  metric_names_num <- c(
    nm_pre6_fit, nm_post_min, nm_post_min_year, nm_yrs_ref_to_post_min, nm_diff_pre6_to_post_min,
    nm_fit_ref_p10, nm_fit_ref_p15, nm_fit_ref_p20,
    nm_diff_post_min_to_p10, nm_diff_post_min_to_p15, nm_diff_post_min_to_p20,
    nm_slope_post_min_to_p10, nm_slope_post_min_to_p15, nm_slope_post_min_to_p20,
    nm_perc_recov_p10, nm_perc_recov_p15, nm_perc_recov_p20
  )
  metric_name_logical <- nm_min_at_end
  
  # ---- convert once ----
  dt <- data.table::as.data.table(df)
  dt[, id__temp := .I]
  
  make_stats_na <- function(ids_dt) {
    out <- data.table::data.table(id__temp = ids_dt$id__temp)
    out[, (metric_names_num) := NA_real_]
    out[, (metric_name_logical) := NA]
    out
  }
  
  # ---- Robust exit: no TS columns found ----
  # Fill NA metric columns and return (do not error).
  if (length(ts_cols) == 0) {
    stats_dt <- make_stats_na(dt[, .(id__temp)])
    out <- merge(dt, stats_dt, by = "id__temp", all.x = TRUE)
    out[, id__temp := NULL]
    return(tibble::as_tibble(out))
  }
  
  # ---- melt ----
  long_dt <- data.table::melt(
    dt,
    id.vars = c("id__temp", reference_time_col),
    measure.vars = ts_cols,
    variable.name = "year_col",
    value.name = "value"
  )
  
  long_dt[, year := as.integer(sub(column_prefix, "", year_col))]
  long_dt[, ref_year := suppressWarnings(as.integer(get(reference_time_col)))]
  
  # Keep only rows with a usable reference year
  long_dt <- long_dt[!is.na(ref_year)]
  
  # ---- Robust exit: no rows survive ref-year filter ----
  if (nrow(long_dt) == 0L) {
    stats_dt <- make_stats_na(dt[, .(id__temp)])
    out <- merge(dt, stats_dt, by = "id__temp", all.x = TRUE)
    out[, id__temp := NULL]
    return(tibble::as_tibble(out))
  }
  
  dbg <- new.env(parent = emptyenv())
  dbg$fail_bad_ref <- 0L
  dbg$fail_min_n <- 0L
  dbg$fail_years_lt3 <- 0L
  dbg$fail_gam <- 0L
  dbg$ok <- 0L
  dbg$first_gam_error <- NULL
  
  na_row <- function(id) {
    out <- data.table::data.table(id__temp = id)
    out[, (metric_names_num) := NA_real_]
    out[, (metric_name_logical) := NA]
    out
  }
  
  fit_one <- function(d) {
    id <- unique(d$id__temp)[1]
    
    ref <- unique(d$ref_year)
    if (length(ref) != 1L || is.na(ref)) {
      if (isTRUE(debug)) dbg$fail_bad_ref <- dbg$fail_bad_ref + 1L
      return(na_row(id))
    }
    
    d_fit <- d[!is.na(value) & !is.na(year)]
    if (nrow(d_fit) < min_n) {
      if (isTRUE(debug)) dbg$fail_min_n <- dbg$fail_min_n + 1L
      return(na_row(id))
    }
    
    years_obs <- sort(unique(d_fit$year))
    if (length(years_obs) < 3) {
      if (isTRUE(debug)) dbg$fail_years_lt3 <- dbg$fail_years_lt3 + 1L
      return(na_row(id))
    }
    
    d_fit[, t := year - ref]
    t_obs <- sort(unique(d_fit$t))
    
    k_eff <- min(as.integer(k), max(3L, length(t_obs) - 1L))
    
    g <- tryCatch(
      mgcv::gam(value ~ s(t, k = k_eff), data = d_fit, method = gam_method),
      error = function(e) {
        if (isTRUE(debug) && is.null(dbg$first_gam_error)) dbg$first_gam_error <- conditionMessage(e)
        NULL
      }
    )
    if (is.null(g)) {
      if (isTRUE(debug)) dbg$fail_gam <- dbg$fail_gam + 1L
      return(na_row(id))
    }
    
    t_min <- min(t_obs)
    t_max <- max(t_obs)
    
    t_pre6 <- -6L
    pre6_fit <- if (t_pre6 >= t_min && t_pre6 <= t_max) {
      as.numeric(stats::predict(g, newdata = data.frame(t = t_pre6), type = "response"))
    } else {
      NA_real_
    }
    
    t_cap <- if (is.finite(post_min_search_max_year)) post_min_search_max_year - ref else Inf
    post_ts <- t_obs[t_obs > 0 & t_obs <= t_cap]
    
    post_min_val <- NA_real_
    post_min_t <- NA_real_
    
    if (length(post_ts) > 0) {
      post_fits <- as.numeric(stats::predict(g, newdata = data.frame(t = post_ts), type = "response"))
      i_min <- which.min(post_fits)
      post_min_val <- post_fits[i_min]
      post_min_t <- post_ts[i_min]
    }
    
    post_min_year <- if (!is.na(post_min_t)) ref + post_min_t else NA_real_
    yrs_ref_to_min <- if (!is.na(post_min_t)) post_min_t else NA_real_
    
    diff_pre6_to_min <- if (!is.na(pre6_fit) && !is.na(post_min_val)) post_min_val - pre6_fit else NA_real_
    
    t_p10 <- 10L
    t_p15 <- 15L
    t_p20 <- 20L
    
    fit_p10 <- if (t_p10 >= t_min && t_p10 <= t_max) {
      as.numeric(stats::predict(g, newdata = data.frame(t = t_p10), type = "response"))
    } else {
      NA_real_
    }
    
    fit_p15 <- if (t_p15 >= t_min && t_p15 <= t_max) {
      as.numeric(stats::predict(g, newdata = data.frame(t = t_p15), type = "response"))
    } else {
      NA_real_
    }
    
    fit_p20 <- if (t_p20 >= t_min && t_p20 <= t_max) {
      as.numeric(stats::predict(g, newdata = data.frame(t = t_p20), type = "response"))
    } else {
      NA_real_
    }
    
    diff_min_to_p10 <- if (!is.na(post_min_val) && !is.na(fit_p10)) fit_p10 - post_min_val else NA_real_
    diff_min_to_p15 <- if (!is.na(post_min_val) && !is.na(fit_p15)) fit_p15 - post_min_val else NA_real_
    diff_min_to_p20 <- if (!is.na(post_min_val) && !is.na(fit_p20)) fit_p20 - post_min_val else NA_real_
    
    slope_to_p10 <- if (!is.na(post_min_t) && !is.na(post_min_val) && !is.na(fit_p10) && (t_p10 != post_min_t)) {
      (fit_p10 - post_min_val) / (t_p10 - post_min_t)
    } else {
      NA_real_
    }
    
    slope_to_p15 <- if (!is.na(post_min_t) && !is.na(post_min_val) && !is.na(fit_p15) && (t_p15 != post_min_t)) {
      (fit_p15 - post_min_val) / (t_p15 - post_min_t)
    } else {
      NA_real_
    }
    
    slope_to_p20 <- if (!is.na(post_min_t) && !is.na(post_min_val) && !is.na(fit_p20) && (t_p20 != post_min_t)) {
      (fit_p20 - post_min_val) / (t_p20 - post_min_t)
    } else {
      NA_real_
    }
    
    min_at_end <- if (!is.na(post_min_t) && length(post_ts) > 0) {
      isTRUE(post_min_t == max(post_ts))
    } else {
      NA
    }
    
    denom <- abs(diff_pre6_to_min)
    perc_recov_p10 <- if (!is.na(diff_min_to_p10) && !is.na(denom) && denom > 0) (diff_min_to_p10 / denom) * 100 else NA_real_
    perc_recov_p15 <- if (!is.na(diff_min_to_p15) && !is.na(denom) && denom > 0) (diff_min_to_p15 / denom) * 100 else NA_real_
    perc_recov_p20 <- if (!is.na(diff_min_to_p20) && !is.na(denom) && denom > 0) (diff_min_to_p20 / denom) * 100 else NA_real_
    
    if (isTRUE(debug)) dbg$ok <- dbg$ok + 1L
    
    out <- data.table::data.table(
      id__temp = id,
      pre6_fit = pre6_fit,
      post_min_val = post_min_val,
      post_min_year = post_min_year,
      yrs_ref_to_min = yrs_ref_to_min,
      diff_pre6_to_min = diff_pre6_to_min,
      fit_p10 = fit_p10,
      fit_p15 = fit_p15,
      fit_p20 = fit_p20,
      diff_min_to_p10 = diff_min_to_p10,
      diff_min_to_p15 = diff_min_to_p15,
      diff_min_to_p20 = diff_min_to_p20,
      slope_to_p10 = slope_to_p10,
      slope_to_p15 = slope_to_p15,
      slope_to_p20 = slope_to_p20,
      min_at_end = min_at_end,
      perc_recov_p10 = perc_recov_p10,
      perc_recov_p15 = perc_recov_p15,
      perc_recov_p20 = perc_recov_p20
    )
    
    data.table::setnames(
      out,
      old = c(
        "pre6_fit",
        "post_min_val",
        "post_min_year",
        "yrs_ref_to_min",
        "diff_pre6_to_min",
        "fit_p10",
        "fit_p15",
        "fit_p20",
        "diff_min_to_p10",
        "diff_min_to_p15",
        "diff_min_to_p20",
        "slope_to_p10",
        "slope_to_p15",
        "slope_to_p20",
        "min_at_end",
        "perc_recov_p10",
        "perc_recov_p15",
        "perc_recov_p20"
      ),
      new = c(
        nm_pre6_fit,
        nm_post_min,
        nm_post_min_year,
        nm_yrs_ref_to_post_min,
        nm_diff_pre6_to_post_min,
        nm_fit_ref_p10,
        nm_fit_ref_p15,
        nm_fit_ref_p20,
        nm_diff_post_min_to_p10,
        nm_diff_post_min_to_p15,
        nm_diff_post_min_to_p20,
        nm_slope_post_min_to_p10,
        nm_slope_post_min_to_p15,
        nm_slope_post_min_to_p20,
        nm_min_at_end,
        nm_perc_recov_p10,
        nm_perc_recov_p15,
        nm_perc_recov_p20
      )
    )
    
    out
  }
  
  # Split by id__temp. This can still be empty in edge cases; handle it.
  split_list <- split(long_dt, long_dt$id__temp)
  
  if (length(split_list) == 0L) {
    stats_dt <- make_stats_na(dt[, .(id__temp)])
    out <- merge(dt, stats_dt, by = "id__temp", all.x = TRUE)
    out[, id__temp := NULL]
    return(tibble::as_tibble(out))
  }
  
  if (isTRUE(parallel)) {
    if (!requireNamespace("future.apply", quietly = TRUE)) stop("Please install 'future.apply' for parallel=TRUE.")
    if (!requireNamespace("future", quietly = TRUE)) stop("Please install 'future' for parallel=TRUE.")
    old_plan <- future::plan()
    on.exit(future::plan(old_plan), add = TRUE)
    
    future::plan(future::multisession, workers = n_cores)
    res_list <- future.apply::future_lapply(split_list, fit_one)
  } else {
    res_list <- lapply(split_list, fit_one)
  }
  
  if (length(res_list) == 0L) {
    stats_dt <- make_stats_na(dt[, .(id__temp)])
  } else {
    stats_dt <- data.table::rbindlist(res_list, use.names = TRUE, fill = TRUE)
    
    # If something went very wrong, force a safe NA stats_dt
    if (!("id__temp" %in% names(stats_dt))) {
      stats_dt <- make_stats_na(dt[, .(id__temp)])
    }
  }
  
  out <- merge(dt, stats_dt, by = "id__temp", all.x = TRUE)
  out[, id__temp := NULL]
  
  if (isTRUE(debug)) {
    message(
      "compute_gam_stats debug counts:\n",
      "ok: ", dbg$ok, "\n",
      "fail_bad_ref: ", dbg$fail_bad_ref, "\n",
      "fail_min_n: ", dbg$fail_min_n, "\n",
      "fail_years_lt3: ", dbg$fail_years_lt3, "\n",
      "fail_gam: ", dbg$fail_gam
    )
    if (!is.null(dbg$first_gam_error)) message("First GAM error: ", dbg$first_gam_error)
  }
  
  tibble::as_tibble(out)
}


compute_forest_response_offset <- function(dt, offsets, new_col_name, prefix = "rap_TRE_") {
  # Ensure dt is a data.table
  dt <- as.data.table(copy(dt))
  dt[, id__temp := .I]  # Create unique ID for merge
  
  # Melt relevant columns
  melt_cols <- grep(paste0("^", prefix, "\\d{4}$"), names(dt), value = TRUE)
  
  long_dt <- melt(
    dt,
    id.vars = c("id__temp", "burn_year"),
    measure.vars = melt_cols,
    variable.name = "year_col",
    value.name = "value"
  )
  
  # Extract year from column name
  long_dt[, year := as.integer(sub(prefix, "", year_col))]
  
  # Calculate offset from burn year
  long_dt[, offset := year - burn_year]
  
  # Only keep rows with desired offsets
  long_dt <- long_dt[offset %in% offsets]
  
  # Now, summarize carefully
  summary_dt <- long_dt[, .(
    n_available = sum(!is.na(value)),
    n_expected = length(offsets),
    avg = if (sum(!is.na(value)) == length(offsets)) mean(value, na.rm = TRUE) else NA_real_
  ), by = id__temp]
  
  # Keep only id and avg
  summary_dt <- summary_dt[, .(id__temp, avg)]
  
  # Merge back
  dt <- merge(dt, summary_dt, by = "id__temp", all.x = TRUE)
  setnames(dt, "avg", new_col_name)
  dt[, id__temp := NULL]
  
  return(as_tibble(dt))
}


#' Compute Post-Fire Forest Trend Slope from forest cover Time Series
#'
#' This function calculates the slope and R-squared value of a thiel-sen regression
#' fitted to post-fire remote sensing time series data (e.g., RAP TRE) over
#' specified year offsets. The output includes the computed slope and R² for each
#' observation where sufficient data is available.
#'
#' @param dt A `data.table` or `data.frame` containing RAP data and a `burn_year` column.
#' @param offsets Integer vector indicating year offsets from `burn_year` to include in the slope calculation. Defaults to `2:10`.
#' @param req_yrs Integer indicating the minimum number of non-NA observations required to fit a regression. Defaults to `5`.
#' @param new_col_name Name of the column to store computed slope values. Defaults to `"post_fire_rap_slope"`.
#' @param r2_col_name Name of the column to store R-squared values. Defaults to `"post_fire_rap_r2"`.
#' @param prefix Prefix of column names containing yearly response values. These should be named like `prefixYYYY` (e.g., `"rap_TRE_2012"`).
#'
#' @return A `tibble` with the original data and two additional columns containing slope and R-squared values.
#'
#' @details
#' The function internally melts the  time series into long format, filters values
#' by offset from the `burn_year`, and fits a thiel-sen model (`value ~ offset`) for each
#' row. If the number of non-NA years is below `req_yrs`, NA values are returned.
#'
#' @importFrom data.table as.data.table copy melt setnames
#' @importFrom stats lm coef summary
#' @importFrom tibble as_tibble
#'
#' @examples
#' \dontrun{
#' dt <- data.table(
#'   burn_year = c(2010, 2012),
#'   rap_TRE_2011 = c(0.5, 0.6),
#'   rap_TRE_2012 = c(0.6, 0.65),
#'   rap_TRE_2013 = c(0.7, NA),
#'   rap_TRE_2014 = c(0.8, 0.75)
#' )
#' compute_post_fire_forest_slope_dt(dt, offsets = 1:4, req_yrs = 3)
#' }
#'
#' @export
compute_post_fire_forest_slope_dt <- function(
    dt,
    offsets = 2:10,
    req_yrs = 5,
    new_col_name = "post_fire_rap_slope",
    r2_col_name = "post_fire_rap_r2",
    prefix = "rap_TRE_"
) {
  dt <- as.data.table(copy(dt))
  dt[, id__temp := .I]  # Unique ID
  
  # Melt RAP TRE columns
  melt_cols <- grep(paste0("^", prefix, "\\d{4}$"), names(dt), value = TRUE)
  
  long_dt <- melt(
    dt,
    id.vars = c("id__temp", "burn_year"),
    measure.vars = melt_cols,
    variable.name = "year_col",
    value.name = "value"
  )
  
  # Extract year
  long_dt[, year := as.integer(sub(prefix, "", year_col))]
  
  # Calculate offset from burn year
  long_dt[, offset := year - burn_year]
  
  # Keep only desired offsets
  long_dt <- long_dt[offset %in% offsets]
  
  # Group and compute slope + R2 if enough data
  slope_dt <- long_dt[
    !is.na(value),
    {
      if (.N >= req_yrs) {
        # model <- lm(value ~ offset, data = .SD) ####################SWITCH TO USE THEIL-SEN REGRESSION
        # list(
        #   slope = coef(model)[["offset"]],
        #   r2 = summary(model)$r.squared
        # )
        model <- mblm::mblm(value ~ offset, data = .SD, repeated = FALSE)
        preds <- predict(model)
        r2 <- 1 - sum((.SD$value - preds)^2, na.rm = TRUE) / sum((.SD$value - mean(.SD$value, na.rm = TRUE))^2, na.rm = TRUE)
        
        list(
          slope = coef(model)[["offset"]],
          r2 = r2
        )
      } else {
        list(
          slope = NA_real_,
          r2 = NA_real_
        )
      }
    },
    by = id__temp
  ]
  
  # Merge results back
  dt <- merge(dt, slope_dt, by = "id__temp", all.x = TRUE)
  setnames(dt, c("slope", "r2"), c(new_col_name, r2_col_name))
  dt[, id__temp := NULL]
  
  return(as_tibble(dt))
}


extract_year_of_value <- function(df, column_pattern, year_col = "mock_burn_year") {
  
  matched_cols <- grep(column_pattern, names(df), value = TRUE)
  if (length(matched_cols) == 0) return(df)
  
  base_name <- sub("^(.*?)_\\d{4}$", "\\1", matched_cols[1])
  
  col_years <- as.integer(sub(".*_(\\d{4})$", "\\1", matched_cols))
  col_year_map <- setNames(matched_cols, as.character(col_years))
  
  n <- nrow(df)
  out_vals <- rep(NA_real_, n)
  
  yr_vec <- df[[year_col]]
  
  for (i in seq_len(n)) {
    
    yr <- yr_vec[i]
    if (is.na(yr)) next
    
    yr_i <- as.integer(yr)
    
    col <- col_year_map[as.character(yr_i)]  # returns NA if not present
    if (is.na(col)) next
    
    out_vals[i] <- df[[col]][i]
  }
  
  new_col <- paste0(base_name, "_year_of")
  df[[new_col]] <- out_vals
  
  df
}


compute_avgs <- function(dats, prefix) {
  dats <- dats |>
    compute_forest_response_offset(offsets = 9:11,
                                   new_col_name = glue("raw_avg_{prefix}10_yr_post_fire"),
                                   prefix = prefix) |>
    compute_forest_response_offset(offsets = 14:16,
                                   new_col_name = glue("raw_avg_{prefix}15_yr_post_fire"),
                                   prefix = prefix) |>
    compute_forest_response_offset(offsets = 19:21,
                                   new_col_name = glue("raw_avg_{prefix}20_yr_post_fire"),
                                   prefix = prefix) |>
    compute_forest_response_offset(offsets = -7:-5,
                                   new_col_name = glue("raw_avg_{prefix}6_yr_pre_fire"),
                                   prefix = prefix) |>
    compute_forest_response_offset(offsets = 2:4, 
                                   new_col_name = glue("raw_avg_{prefix}3_yr_post_fire"),  # USE 3 YRS POST FIRE SO THAT WINDOW CAPTURES 2-4, missing most delayed mortality
                                   prefix = prefix)
  return(dats)
}

compute_response_differences_recovery <- function(dats, prefix) {
  dats <- dats |>
    mutate(
      # n6 -> 3 year difference
      !!sym(glue("raw_post_fire_{prefix}difference_n6_3_yr")) :=
        !!sym(glue("raw_avg_{prefix}3_yr_post_fire")) -
        !!sym(glue("raw_avg_{prefix}6_yr_pre_fire")),
      
      # increase from 3 -> 10 years post-fire
      !!sym(glue("raw_post_fire_{prefix}increase_3_10_yr")) :=
        !!sym(glue("raw_avg_{prefix}10_yr_post_fire")) -
        !!sym(glue("raw_avg_{prefix}3_yr_post_fire")),
      
      # increase from 3 -> 15 years post-fire
      !!sym(glue("raw_post_fire_{prefix}increase_3_15_yr")) :=
        !!sym(glue("raw_avg_{prefix}15_yr_post_fire")) -
        !!sym(glue("raw_avg_{prefix}3_yr_post_fire")),
      
      # increase from 3 -> 20 years post-fire
      !!sym(glue("raw_post_fire_{prefix}increase_3_20_yr")) :=
        !!sym(glue("raw_avg_{prefix}20_yr_post_fire")) -
        !!sym(glue("raw_avg_{prefix}3_yr_post_fire")),
      
      # percent recovery (3–10yr)
      !!sym(glue("raw_post_fire_{prefix}perc_recovery_n6_3_10_yr")) :=
        (!!sym(glue("raw_post_fire_{prefix}increase_3_10_yr")) /
           (-1 * !!sym(glue("raw_post_fire_{prefix}difference_n6_3_yr")))) * 100,
      
      # percent recovery (3–15yr)
      !!sym(glue("raw_post_fire_{prefix}perc_recovery_n6_3_15_yr")) :=
        (!!sym(glue("raw_post_fire_{prefix}increase_3_15_yr")) /
           (-1 * !!sym(glue("raw_post_fire_{prefix}difference_n6_3_yr")))) * 100,
      
      # percent recovery (3–20yr)
      !!sym(glue("raw_post_fire_{prefix}perc_recovery_n6_3_20_yr")) :=
        (!!sym(glue("raw_post_fire_{prefix}increase_3_20_yr")) /
           (-1 * !!sym(glue("raw_post_fire_{prefix}difference_n6_3_yr")))) * 100
    )
  return(dats)
}  

# Allow for up to 1 missing year for each slope computation
compute_response_slopes <- function(dats, prefix) {
  dats <- dats |>
    compute_post_fire_forest_slope_dt(
      offsets = 3:10,
      req_yrs = 7,
      new_col_name = glue("raw_post_fire_{prefix}slope_3_10_yrs"),
      r2_col_name = glue("raw_post_fire_{prefix}3_10_yrs_r2"),
      prefix = prefix) |>
    compute_post_fire_forest_slope_dt(
      offsets = 3:15,
      req_yrs = 12,
      new_col_name = glue("raw_post_fire_{prefix}slope_3_15_yrs"),
      r2_col_name = glue("raw_post_fire_{prefix}3_15_yrs_r2"),
      prefix = prefix) |>
    compute_post_fire_forest_slope_dt(
      offsets = 3:20,
      req_yrs = 17,
      new_col_name = glue("raw_post_fire_{prefix}slope_3_20_yrs"),
      r2_col_name = glue("raw_post_fire_{prefix}3_20_yrs_r2"),
      prefix = prefix)
  return(dats)
}







future::plan(future::multisession, workers = 4)
opts <- furrr::furrr_options(
  seed = TRUE,
  packages = c(
    "arrow",      # open_dataset / write_parquet
    "dplyr",      # filter/collect/mutate/select
    "tidyr",      # pivot_longer / drop_na
    "tibble",     # as_tibble
    "rlang",      # sym / tidy eval (your compute_* fns)
    "glue",       # glue() / glue strings
    "data.table", # melt / as.data.table / rbindlist etc
    "mgcv",       # gam
    "mblm",       # Theil–Sen regression
    "fect"        # get.cohort
  )
)


# Started 5:24pm
parquet_path <- file.path(dir_derived, "merged_dats_raw.parquet")

# dats_by_l3 <- furrr::future_map(
#   forested_ecoregions$code_name,
#   process_ecoregion,
#   parquet_path = parquet_path,
#   out_dir = dir_parquet,
#   utils_path = here::here("utils", "functions.R"),
#   .options = opts
# )


tic()
dats_by_l3 <- purrr::map(
  forested_ecoregions$code_name,
  process_ecoregion,
  parquet_path = parquet_path,
  out_dir = dir_parquet,
  utils_path = here::here("utils", "functions.R")
)
toc()





