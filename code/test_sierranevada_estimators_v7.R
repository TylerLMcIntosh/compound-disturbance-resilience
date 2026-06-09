
rm(list = ls())

cyverse = FALSE

if(!requireNamespace("here", quietly = TRUE)) {
  install.packages("here")
}
library(here)

source(here::here("utils", "functions.R"))
source(here::here("code", "sunab_aggregate_vcov.R"))


install_load_packages(c(
  "tidyverse",
  "tictoc",
  "glue",
  "fixest",
  "arrow",
  "did",
  "didimputation"
))

run <- "GEE_resilience_v6_operational_ss500_ts50000"


if(cyverse) {
  dir_figs <- file.path("~/data-store/data/iplant/home/shared/earthlab/macrosystems/tlm", run, "figs")
  dir_derived <- file.path("~/data-store/data/iplant/home/shared/earthlab/macrosystems/tlm", run, "data/derived")
  dir_parquet <- here::here(dir_derived, 'parquet')
  dir_manual <- file.path("~/data-store/data/iplant/home/shared/earthlab/macrosystems/tlm", run, "data/manual")
  dir_raw <- file.path("~/data-store/data/iplant/home/shared/earthlab/macrosystems/tlm", run, "data/raw")
  dir_results <- file.path("~/data-store/data/iplant/home/shared/earthlab/macrosystems/tlm", run, "results")
} else {
  dir_figs <- here::here('figs', run)
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

set_cd_groups <- function(df,
                          group_col,
                          b_nm,
                          d_nm,
                          b_threshold,
                          d_threshold,
                          dummy_cols      = c("cd_f", "cd_bf", "cd_df", "cd_bdf"),
                          include_control_group = FALSE) {
  
  if (d_threshold < 0) {
    df_new <- df |>
      dplyr::mutate(
        "{group_col}" := dplyr::case_when(
          fire == 1 & .data[[b_nm]] <  b_threshold & .data[[d_nm]] >  d_threshold ~ "f",
          fire == 1 & .data[[b_nm]] >= b_threshold & .data[[d_nm]] >  d_threshold ~ "bf",
          fire == 1 & .data[[b_nm]] <  b_threshold & .data[[d_nm]] <= d_threshold ~ "df",
          fire == 1 & .data[[b_nm]] >= b_threshold & .data[[d_nm]] <= d_threshold ~ "bdf",
          include_control_group & fire == 0                                               ~ "control",
          TRUE ~ NA_character_
        )
      )
  } else {
    df_new <- df |>
      dplyr::mutate(
        "{group_col}" := dplyr::case_when(
          fire == 1 & .data[[b_nm]] <  b_threshold & .data[[d_nm]] <  d_threshold ~ "f",
          fire == 1 & .data[[b_nm]] >= b_threshold & .data[[d_nm]] <  d_threshold ~ "bf",
          fire == 1 & .data[[b_nm]] <  b_threshold & .data[[d_nm]] >= d_threshold ~ "df",
          fire == 1 & .data[[b_nm]] >= b_threshold & .data[[d_nm]] >= d_threshold ~ "bdf",
          include_control_group & fire == 0                                               ~ "control",
          TRUE ~ NA_character_
        )
      )
  }
  
  # f as reference level for WeightIt factor
  df_new <- df_new |>
    dplyr::mutate("{group_col}" := relevel(factor(.data[[group_col]]), ref = "f"))
  
  # binary dummies for feols; NAs become 0 (control units get 0 for all dummies)
  level_names <- c("f", "bf", "df", "bdf")
  for (i in seq_along(dummy_cols)) {
    lv <- level_names[i]
    df_new[[dummy_cols[i]]] <- tidyr::replace_na(
      as.integer(!is.na(df_new[[group_col]]) & as.character(df_new[[group_col]]) == lv),
      0L
    )
  }
  
  df_new
}

# Set up data ----
set <- "sierranevada"
set_col <- "ecoregion_code_name"

small_test <- FALSE
n_test <- 1000

aggregation_code <- "(year::-?[0-9]+):cohort::[0-9]+:(cd.*)"


dats_long <- arrow::open_dataset(here(dir_derived, "parquet_long_filtered_test")) |> 
  filter(.data[[set_col]] == set) |>
#  filter(year >= 1997) |>
  collect()
dats_short <- arrow::open_dataset(here(dir_derived, "parquet_short_filtered_test")) |>
  filter(.data[[set_col]] == set) |>
  collect()

# create minimal set
if(small_test) {
  test_pts <- sample(unique(dats_long$pt_id), n_test)
  
  dats_long <- dats_long |>
    dplyr::filter(
      pt_id %in% test_pts
    )
  
  dats_short <- dats_short |>
    dplyr::filter(
      pt_id %in% test_pts
    )
}


# add categories

cd_nm <- "b10_pdsisumn10"

dats_long <- dats_long |>
  set_cd_groups(
    group_col = "b10_pdsisumn10",
    b_nm = "biotic_relaxedforestnorm_5_yrs_prior_sum_yot",
    d_nm = "pdsi_annual_5_yrs_prior_sum_yot",
    b_threshold = 10,
    d_threshold = -10,
    include_control_group = TRUE
  )




# Test dummy version ----

est_sunab_dummy <- feols(
  rap_tree ~
    sunab(FirstTreat, year, ref.p = -6, no_agg = TRUE):cd_f +
    sunab(FirstTreat, year, ref.p = -6, no_agg = TRUE):cd_bf +
    sunab(FirstTreat, year, ref.p = -6, no_agg = TRUE):cd_df +
    sunab(FirstTreat, year, ref.p = -6, no_agg = TRUE):cd_bdf |
    pt_id + year,
  data = dats_long,
  cluster = ~ pt_id
)




# see dummy coefficient names
cn_dummy <- names(coef(est_sunab_dummy, agg = FALSE))
cat(cn_dummy, sep = "\n")


# Aggregate using internal fixest function (does not expose full vcov)
agg_cd_event_dummy <- aggregate(
  est_sunab_dummy,
  agg = aggregation_code
)

glimpse(agg_cd_event_dummy)

# Double-check aggregation using new function (full test in test_aggregate_with_vcov_v7.R)

aggnew_cd_event_dummy <- sunab_aggregate_vcov(
  est_sunab_dummy,
  agg = aggregation_code,
  weight_method = "model_matrix"
)


# Test weighting ----

weight_dats <- dats_long |>
  distinct(pt_id, .keep_all = TRUE) |>
  filter(!is.na(b10_pdsisumn10)) |>
  left_join(dats_short |> select(pt_id, gam_rap_tree_pre6_fit)
)


# generalized overlap weighting via multinomial propensity scores
w_out_glm <- weightit(
  b10_pdsisumn10 ~ chili + def + aet + srtm + tpi + nfg_factor,
  data = weight_dats,
  method = "glm",
  estimand = "ATO"
)

summary(w_out_glm)

cobalt::love.plot(w_out_glm)

dats_long <- dats_long |>
  left_join(weight_dats |>
              mutate(weights = w_out_glm$weights) |>
              select(pt_id, weights), by = c("pt_id"))

est_sunab_dummy_weighted <- feols(
  rap_tree ~
    sunab(FirstTreat, year, ref.p = -6, no_agg = TRUE):cd_f +
    sunab(FirstTreat, year, ref.p = -6, no_agg = TRUE):cd_bf +
    sunab(FirstTreat, year, ref.p = -6, no_agg = TRUE):cd_df +
    sunab(FirstTreat, year, ref.p = -6, no_agg = TRUE):cd_bdf |
    pt_id + year,
  data = dats_long,
  cluster = ~ pt_id,
  weights = dats_long$weights
)



# Aggregate using internal fixest function (does not expose full vcov)
agg_cd_event_dummy_weighted <- aggregate(
  est_sunab_dummy_weighted,
  agg = aggregation_code
)

# Double-check aggregation using new function (full test in test_aggregate_with_vcov_v7.R)

aggnew_cd_event_dummy_weighted <- sunab_aggregate_vcov(
  est_sunab_dummy_weighted,
  agg = aggregation_code,
  weight_method = "model_matrix"
)


# PLOTS ----

create_plot <- function(model, version) {
  
  es_df_dummy <- as.data.frame(model) |>
    rownames_to_column("term") |>
    rename(
      estimate = Estimate,
      se = `Std. Error`
    ) |>
    mutate(
      event_time = as.integer(str_extract(term, "(?<=year::)-?[0-9]+")),
      cd = str_extract(term, "cd_.*$"),
      conf_low = estimate - 1.96 * se,
      conf_high = estimate + 1.96 * se
    )
  
  p1 <- ggplot(es_df_dummy, aes(x = event_time, y = estimate, color = cd)) +
    geom_hline(yintercept = 0, linewidth = 0.3) +
    geom_vline(xintercept = -1, linetype = "dashed", linewidth = 0.3) +
    geom_ribbon(
      aes(ymin = conf_low, ymax = conf_high, fill = cd),
      alpha = 0.15,
      color = NA
    ) +
    geom_line() +
    geom_point() +
    labs(
      x = "Years since treatment",
      y = "Estimated effect on rap_tree",
      color = "Group",
      fill = "Group",
      title = version
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  
  return(p1)  
  
}

create_plot(agg_cd_event_dummy, version = "raw dummy")
create_plot(aggnew_cd_event_dummy$feols_structure, version = "new aggregation dummy")

create_plot(agg_cd_event_dummy_weighted, version = "weighted raw dummy")
create_plot(aggnew_cd_event_dummy_weighted$feols_structure, version = "weighted new aggregation dummy")








