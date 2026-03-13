
# ~6 hours for all ecoregions with run_grid 120

# Setup ----

rm(list = ls())

cyverse = FALSE

if(!requireNamespace("here", quietly = TRUE)) {
  install.packages("here")
}
library(here)

source(here::here("utils", "functions.R"))


install_load_packages(c(
  "tidyverse",
  "tictoc",
  "glue",
  "fixest",
  "arrow"
))

source(here::here("code", "portable_sunab2.R"))

run <- "GEE_resilience_v6_operational_ss500_ts50000"


if(cyverse) {
  dir_figs <- file.path("~/data-store/data/iplant/home/shared/earthlab/macrosystems/tlm", run, "figs")
  dir_derived <- file.path("~/data-store/data/iplant/home/shared/earthlab/macrosystems/tlm", run, "data/derived")
  dir_parquet <- here::here(dir_derived, 'parquet')
  dir_manual <- file.path("~/data-store/data/iplant/home/shared/earthlab/macrosystems/tlm", run, "data/manual")
  dir_raw <- file.path("~/data-store/data/iplant/home/shared/earthlab/macrosystems/tlm", run, "data/raw")
  dir_results <- file.path("~/data-store/data/iplant/home/shared/earthlab/macrosystems/tlm", run, "results")
} else {
  dir_figs <- here::here('figs')
  dir_derived <- here::here('data', 'derived', run)
  dir_parquet <- here::here(dir_derived, 'parquet')
  dir_raw <- here::here('data', 'raw')
  dir_manual <- here::here('data', 'manual')
  dir_results <- here::here('results')
}

dir_ensure(c(dir_figs,
             dir_derived,
             dir_parquet,
             dir_manual,
             dir_raw,
             dir_results))


seed = 1234
set.seed(seed)



# Set up ecoregions ----
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

# Remove ecoregions with <30% USFS forested land coverage
forested_ecoregions <- forested_ecoregions |> filter(code_name != "coastrange" &
                                                       code_name != "pugetlowland" & 
                                                       code_name != "willamettevalley" &
                                                       code_name != "coloradoplateaus" &
                                                       code_name != "centralcaliforniamtns")



# group setting function ----
set_cd_groups <- function(df,
                          b_nm,
                          d_nm,
                          cd_nm,
                          b_threshold,
                          d_threshold) {
  
  if(d_threshold < 0) {
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
  
  
  return(df_new)
}

# Specs ----


make_analysis_specs_ecoregions <- function(forested_ecoregions,
                                           dir_data,
                                           parquet_pattern_prefix = "dats_long_") {
  
  if (!"code_name" %in% names(forested_ecoregions)) {
    stop("forested_ecoregions must contain a 'code_name' column.")
  }
  
  parquet_dir <- file.path(dir_data, "parquet")
  
  analysis_specs <- forested_ecoregions |>
    dplyr::distinct(code_name) |>
    dplyr::transmute(
      analysis_id = code_name,
      parquet_files = purrr::map(
        code_name,
        \(x) file.path(parquet_dir, glue::glue("{parquet_pattern_prefix}{x}.parquet"))
      )
    )
  
  missing_any <- analysis_specs |>
    dplyr::mutate(
      missing = purrr::map_lgl(parquet_files, \(x) !all(file.exists(x)))
    ) |>
    dplyr::filter(missing)
  
  if (nrow(missing_any) > 0) {
    warning(
      "Some parquet files listed in analysis_specs do not exist yet: ",
      paste(missing_any$analysis_id, collapse = ", ")
    )
  }
  
  analysis_specs
}


ecoregion_analysis_specs <- make_analysis_specs_ecoregions(
  forested_ecoregions = forested_ecoregions,
  dir_data = dir_derived
)

full_analysis_specs <- tibble::tibble(
  analysis_id = "zalleco",
  parquet_files = list(
    list.files(
      file.path(dir_derived, "parquet"),
      pattern = "^dats_long_.*\\.parquet$",
      full.names = TRUE
    )
  )
)

analysis_specs <- rbind(ecoregion_analysis_specs,
                        full_analysis_specs)


outcome_specs <- tibble::tibble(
  outcome = c("rap_tree", "vcf_tree")
)


group_specs <- tibble::tibble(
  group_id = c(
    "b10_pdsin4t1",
    "b10_pdsin3t1",
    "b10_pdsisumn10",
    "b25_pdsisumn10"
  ),
  group_fun = list(
    set_cd_groups,
    set_cd_groups,
    set_cd_groups,
    set_cd_groups
  ),
  group_args = list(
    list(
      b_nm = "biotic_relaxedforestnorm_5_yrs_prior_sum_yot",
      d_nm = "pdsi_annual_5_yrs_prior_threshold_n4_yot",
      cd_nm = "b10_pdsin4t1",
      b_threshold = 10,
      d_threshold = 1
    ),
    list(
      b_nm = "biotic_relaxedforestnorm_5_yrs_prior_sum_yot",
      d_nm = "pdsi_annual_5_yrs_prior_threshold_n3_yot",
      cd_nm = "b10_pdsin3t1",
      b_threshold = 10,
      d_threshold = 1
    ),
    list(
      b_nm = "biotic_relaxedforestnorm_5_yrs_prior_sum_yot",
      d_nm = "pdsi_annual_5_yrs_prior_sum_yot",
      cd_nm = "b10_pdsisumn10",
      b_threshold = 10,
      d_threshold = -10
    ),
    list(
      b_nm = "biotic_relaxedforestnorm_5_yrs_prior_sum_yot",
      d_nm = "pdsi_annual_5_yrs_prior_sum_yot",
      cd_nm = "b25_pdsisumn10",
      b_threshold = 25,
      d_threshold = -10
    )
  )
)


model_specs <- tibble::tibble(
  model_id = c("b_pdsi_sunab_twfe"),
  formula_template = c(
    "{outcome} ~ biotic_relaxedforestnorm + pdsi_annual + sunab(FirstTreat, year, ref.p = -6) | pt_id + year"
  ),
  term_pattern = c("^year::"),
  unit_id = c("pt_id"),
  event_id = c("fireid")
)


ecor_vcov_specs <- tibble::tibble(
  vcov_id = c(
    "iid",
    "cluster_pt",
    "cluster_eco4",
    # "conley_50km_20km",
    # "conley_75km_20km",
    "conley_125km_20km"#,
    # "conley_150km_20km",
    # "conley_50km_10km"
  ),
  vcov_label = c(
    "IID Classical",
    "Cluster pt_id",
    "Clustered SEs: eco4",
    # "Conley SEs: 50 km cutoff, 20 km pixel",
    # "Conley SEs: 75 km cutoff, 20 km pixel",
    "Conley SEs: 125 km cutoff, 20 km pixel"#,
    # "Conley SEs: 150 km cutoff, 20 km pixel",
    # "Conley SEs: 50 km cutoff, 10 km pixel"
  ),
  vcov = list(
    "iid",
    stats::as.formula("~pt_id"),
    stats::as.formula("~us_l4code"),
    # fixest::vcov_conley(
    #   lat = "lat",
    #   lon = "long",
    #   cutoff = 50,
    #   pixel = 20,
    #   distance = "triangular"
    # ),
    # fixest::vcov_conley(
    #   lat = "lat",
    #   lon = "long",
    #   cutoff = 75,
    #   pixel = 20,
    #   distance = "triangular"
    # ),
    fixest::vcov_conley(
      lat = "lat",
      lon = "long",
      cutoff = 125,
      pixel = 20,
      distance = "triangular"
    )#,
    # fixest::vcov_conley(
    #   lat = "lat",
    #   lon = "long",
    #   cutoff = 150,
    #   pixel = 20,
    #   distance = "triangular"
    # ),
    # fixest::vcov_conley(
    #   lat = "lat",
    #   lon = "long",
    #   cutoff = 50,
    #   pixel = 10,
    #   distance = "triangular"
    # )
  ),
  vcov_vars = list(
    character(0),
    c("pt_id"),
    c("us_l4code"),
    # c("lat", "long"),
    # c("lat", "long"),
    c("lat", "long")#,
    # c("lat", "long"),
    # c("lat", "long")
  )
)


# Run analysis ----

run_grid <- tidyr::crossing(
  analysis_specs,
  outcome_specs,
  group_specs,
  model_specs
)

run_grid |> 
  dplyr::select(analysis_id, outcome, group_id, model_id)



tic()
dydid_results <- run_dydid_experiment(
  analysis_specs = (analysis_specs |> arrange(analysis_id))[1,],
  #analysis_specs = ecoregion_analysis_specs,
  outcome_specs = outcome_specs[1,],
  group_specs = group_specs[2,],
  model_specs = model_specs,
  vcov_specs = ecor_vcov_specs[1,],
  dir_out = dir_results,
  trt_col = "fire",
  ci_level = 0.95,
  data_filter = ~ year >= 1997,
  skip_existing = TRUE,
  .progress = TRUE
)
toc()



all_tables <- rebuild_dydid_tables_from_by_run(
  dir_out = dir_results,
  write_csv = FALSE
)





# Descriptive analysis ----

descriptive_results <- run_descriptive_experiment_simple(
  analysis_specs = analysis_specs,
  outcome_specs = outcome_specs,
  group_specs = group_specs,
  model_specs = model_specs,
  dir_out = dir_results,
  trt_col = "fire",
  data_filter = ~ year >= 1997,
  treated_year_var = "burn_year",
  control_year_var = "mock_burn_year",
  time_var = "year",
  skip_existing = TRUE,
  .progress = TRUE
)

descriptive_tables <- rebuild_descriptive_tables_from_by_run_simple(
  dir_out = dir_results,
  write_csv = TRUE
)

