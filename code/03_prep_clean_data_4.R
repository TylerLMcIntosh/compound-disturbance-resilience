

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



future::plan(future::multisession, workers = max(1L, parallel::detectCores() / 2))

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
