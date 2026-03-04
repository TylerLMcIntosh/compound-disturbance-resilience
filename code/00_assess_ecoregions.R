
rm(list = ls())

if(!requireNamespace("here", quietly = TRUE)) {
  install.packages("here")
}
library(here)

source(here::here("utils", "functions.R"))

install_load_packages(c(
  "tidyverse"
))

dir_figs <- here::here('figs')
dir_derived <- here::here('data', 'derived')
dir_raw <- here::here('data', 'raw')
dir_manual <- here::here('data', 'manual')
dir_ensure(c(dir_figs,
             dir_derived,
             dir_manual,
             dir_raw))

d <- read_csv(here(dir_derived, "L3_western_conservative_forest_usfs.csv"))

dd <- d |>
  group_by(na_l3name, us_l3name) |>
  summarize(area_km2_both = sum(area_km2_both),
            area_km2_conservative_forest = sum(area_km2_conservative_forest),
            area_km2_total_in_west = sum(area_km2_total_in_west)) |>
  ungroup() |>
  mutate(perc_both = (area_km2_both / area_km2_total_in_west) * 100,
         perc_forest = (area_km2_conservative_forest / area_km2_total_in_west) * 100,) |>
  arrange(desc(perc_both)) |>
  select(na_l3name, us_l3name, perc_both, perc_forest)


