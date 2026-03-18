
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
  "arrow",
  "brms",
  "rstan",
  "tidybayes",
  "remotes",
  "loo",
  "posterior",
  "bayesplot"
))

# Set up stan and use cmdstanr
remotes::install_github("stan-dev/cmdstanr", force = TRUE)
library(cmdstanr)
cmdstanr::install_cmdstan()
options(brms.backend = "cmdstanr")



source(here::here("code", "portable_sunab2.R"))

run <- "GEE_resilience_v6_operational_ss500_ts50000"


if(cyverse) {
  dir_figs <- file.path("~/data-store/data/iplant/home/shared/earthlab/macrosystems/tlm", run, "figs")
  dir_derived <- file.path("~/data-store/data/iplant/home/shared/earthlab/macrosystems/tlm", run, "data/derived")
  dir_parquet <- here::here(dir_derived, 'parquet')
  dir_manual <- file.path("~/data-store/data/iplant/home/shared/earthlab/macrosystems/tlm", run, "data/manual")
  dir_raw <- file.path("~/data-store/data/iplant/home/shared/earthlab/macrosystems/tlm", run, "data/raw")
  dir_results <- file.path("~/data-store/data/iplant/home/shared/earthlab/macrosystems/tlm", run, "results")
  dir_bayes <- here::here(dir_results, "bayes")
} else {
  dir_figs <- here::here('figs')
  dir_derived <- here::here('data', 'derived', run)
  dir_parquet <- here::here(dir_derived, 'parquet')
  dir_raw <- here::here('data', 'raw')
  dir_manual <- here::here('data', 'manual')
  dir_results <- here::here('results')
  dir_bayes <- here::here(dir_results, "bayes")
}

dir_ensure(c(dir_figs,
             dir_derived,
             dir_parquet,
             dir_manual,
             dir_raw,
             dir_results,
             dir_bayes))


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


d_s_f <- list.files(
  here::here(dir_derived, "parquet"),
  pattern = "dats_short.*\\.parquet$",
  full.names = TRUE
)

d_s <- open_dataset(d_s_f, format = "parquet") |>
  collect()

d_s_fire <- d_s |> filter(fire == 1)


d_s_fire_test <- d_s_fire |>
  filter(ecoregion_code_name == "sierranevada") |>
  slice_sample(n = 10000)

d_s_fire_test <- set_cd_groups(df = d_s_fire_test,
                     b_nm = "biotic_relaxedforestnorm_5_yrs_prior_sum_yot",
                     d_nm = "pdsi_annual_5_yrs_prior_threshold_n4_yot",
                     cd_nm = "cd_group_b10_pdsin3t1",
                     b_threshold = 10,
                     d_threshold = 1) |>
  mutate(cd_group_b10_pdsin3t1 =
           fct_relevel(cd_group_b10_pdsin3t1, "f"))



# Null model
priors_null <- c(
  prior(normal(0.7, 0.5), class = "Intercept"),
  prior(exponential(2), class = "sd"),
  prior(exponential(2), class = "sigma")
)

# Main models with scaled predictors
priors_scaled <- c(
  prior(normal(0.7, 0.5), class = "Intercept"),
  prior(normal(0, 0.3), class = "b"),
  prior(exponential(2), class = "sd"),
  prior(exponential(2), class = "sigma")
)

# Random slopes (if used)
priors_scaled_rs <- c(
  prior(normal(0.7, 0.5), class = "Intercept"),
  prior(normal(0, 0.3), class = "b"),
  prior(exponential(2), class = "sd"),
  prior(student_t(3, 0, 1), class = "sigma")
)

hist(d_s_fire_test$raw_post_fire_rap_tree_difference_n6_3_yr)

tic("m1_null")
m1_null <- brm(
  formula = raw_post_fire_rap_tree_difference_n6_3_yr ~ 1 + 
    (1 | fireid) +
    (1 | us_l4code),
  data = d_s_fire_test |> filter(!is.na(cd_group_b10_pdsin3t1)),
  family = student(),
  prior = priors_null,
  chains = 4, 
  cores = 4,
  iter = 1000, 
  warmup = 500,
  control = list(adapt_delta = 0.99, max_treedepth = 20),
  save_pars = save_pars(all = TRUE)
)
toc()
saveRDS(m1_null, file = here(dir_bayes, "m1_null.rds"))



m2 <- brm(
  formula = raw_post_fire_rap_tree_difference_n6_3_yr ~ cd_group_b10_pdsin3t1 + 
    (1 | fireid) +
    (1 | us_l4code),
  data = d_s_fire_test |> filter(!is.na(cd_group_b10_pdsin3t1)),
  family = student(),
  prior = priors_null,
  chains = 4, 
  cores = 4,
  iter = 1000, 
  warmup = 500,
  control = list(adapt_delta = 0.99, max_treedepth = 20),
  save_pars = save_pars(all = TRUE)
)

saveRDS(m2, file = here(dir_bayes, "m2.rds"))

m2_1 <- brm(
  formula = raw_post_fire_rap_tree_difference_n6_3_yr ~ biotic_relaxedforestnorm_5_yrs_prior_sum_yot + 
    pdsi_annual_5_yrs_prior_sum_yot +
    (1 | fireid) +
    (1 | us_l4code),
  data = d_s_fire_test |> filter(!is.na(cd_group_b10_pdsin3t1)),
  family = student(),
  prior = priors_null,
  chains = 4, 
  cores = 4,
  iter = 1000, 
  warmup = 500,
  control = list(adapt_delta = 0.99, max_treedepth = 20),
  save_pars = save_pars(all = TRUE)
)

saveRDS(m2_1, file = here(dir_bayes, "m2_1.rds"))

m2_2 <- brm(
  formula = raw_post_fire_rap_tree_difference_n6_3_yr ~ biotic_relaxedforestnorm_5_yrs_prior_sum_yot + 
    pdsi_annual_5_yrs_prior_threshold_n3_yot +
    (1 | fireid) +
    (1 | us_l4code),
  data = d_s_fire_test |> filter(!is.na(cd_group_b10_pdsin3t1)),
  family = student(),
  prior = priors_null,
  chains = 4, 
  cores = 4,
  iter = 1000, 
  warmup = 500,
  control = list(adapt_delta = 0.99, max_treedepth = 20),
  save_pars = save_pars(all = TRUE)
)

saveRDS(m2_2, file = here(dir_bayes, "m2_2.rds"))

m2_3 <- brm(
  formula = raw_post_fire_rap_tree_difference_n6_3_yr ~ biotic_relaxedforestnorm_5_yrs_prior_sum_yot + 
    hd_fingerprint_5_yrs_prior_threshold_6_yot +
    (1 | fireid) +
    (1 | us_l4code),
  data = d_s_fire_test |> filter(!is.na(cd_group_b10_pdsin3t1)),
  family = student(),
  prior = priors_null,
  chains = 4, 
  cores = 4,
  iter = 1000, 
  warmup = 500,
  control = list(adapt_delta = 0.99, max_treedepth = 20),
  save_pars = save_pars(all = TRUE)
)

saveRDS(m2_3, file = here(dir_bayes, "m2_3.rds"))



# add pre-fire cover
m3 <- brm(
  formula = raw_post_fire_rap_tree_difference_n6_3_yr ~ cd_group_b10_pdsin3t1 + gam_rap_tree_pre6_fit +
    (1 | fireid) +
    (1 | us_l4code),
  data = d_s_fire_test,
  family = student(),
  prior = priors_null,
  chains = 4, 
  cores = 4,
  iter = 1000, 
  warmup = 500,
  control = list(adapt_delta = 0.99, max_treedepth = 20),
  save_pars = save_pars(all = TRUE)
)

saveRDS(m3, file = here(dir_bayes, "m3.rds"))

# add topo
m4 <- brm(
  formula = raw_post_fire_rap_tree_difference_n6_3_yr ~ cd_group_b10_pdsin3t1 + gam_rap_tree_pre6_fit +
    chili + srtm + tpi +
    (1 | fireid) +
    (1 | us_l4code),
  data = d_s_fire_test,
  family = student(),
  prior = priors_null,
  chains = 4, 
  cores = 4,
  iter = 1000, 
  warmup = 500,
  control = list(adapt_delta = 0.99, max_treedepth = 20),
  save_pars = save_pars(all = TRUE)
)

saveRDS(m4, file = here(dir_bayes, "m4.rds"))

# add climate
m5 <- brm(
  formula = raw_post_fire_rap_tree_difference_n6_3_yr ~ cd_group_b10_pdsin3t1 + gam_rap_tree_pre6_fit +
    chili + srtm + tpi +
    tt_normal_aet + tt_normal_def +
    (1 | fireid) +
    (1 | us_l4code),
  data = d_s_fire_test,
  family = student(),
  prior = priors_null,
  chains = 4, 
  cores = 4,
  iter = 1000, 
  warmup = 500,
  control = list(adapt_delta = 0.99, max_treedepth = 20),
  save_pars = save_pars(all = TRUE)
)

saveRDS(m5, file = here(dir_bayes, "m5.rds"))

# disturbance information
# mean chain time with 10,000 samples: 1082 second (20 min)
m6 <- brm(
  formula = raw_post_fire_rap_tree_difference_n6_3_yr ~ cd_group_b10_pdsin3t1 + gam_rap_tree_pre6_fit +
    chili + srtm + tpi +
    tt_normal_aet + tt_normal_def +
    pdsi_annual_5_yrs_prior_sum_yot + biotic_relaxedforestnorm_5_yrs_prior_sum_yot + #pdsi_annual_5_yrs_after_sum +
    (1 | fireid) +
    (1 | us_l4code),
  data = d_s_fire_test,
  family = student(),
  prior = priors_null,
  chains = 4, 
  cores = 4,
  iter = 1000, 
  warmup = 500,
  control = list(adapt_delta = 0.99, max_treedepth = 20),
  save_pars = save_pars(all = TRUE)
)

saveRDS(m6, file = here(dir_bayes, "m6.rds"))



# Initial model comparison
loo_m1_null <- loo(m1_null)
saveRDS(loo_m1_null, here(dir_bayes, "loo_m1_null.rds"))

loo_m2 <- loo(m2)
loo_m2_1 <- loo(m2_1)
loo_m2_2 <- loo(m2_2)
loo_m2_3 <- loo(m2_3)
saveRDS(loo_m2, here(dir_bayes, "loo_m2.rds"))

loo_m3 <- loo(m3)
loo_m4 <- loo(m4)
loo_m5 <- loo(m5)
loo_m6 <- loo(m6)



# Compare
loo_compare(loo_m1_null, 
            loo_m2,
            loo_m2_1,
            loo_m2_2,
            loo_m2_3,
            loo_m3,
            loo_m4,
            loo_m5#,
            #loo_m6
            )



summary(m4)
bayes_R2(m4)
performance::r2_bayes(m4)
plot(m4)
conditional_effects(m4)


bayesplot::mcmc_intervals(as.array(m4), regex_pars = "^b_")

