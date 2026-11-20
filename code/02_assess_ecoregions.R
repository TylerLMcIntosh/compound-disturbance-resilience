# This script computes percentages of ecoregions covered by USFS-managed forest and outputs a csv and table for the supplement


rm(list = ls())

if(!requireNamespace("here", quietly = TRUE)) {
  install.packages("here")
}
library(here)

source(here::here("utils", "functions.R"))

install_load_packages(c(
  "tidyverse",
  "flextable"
))

dir_figs <- here::here('figs/v9')
dir_derived <- here::here('data', 'derived')
dir_raw <- here::here('data', 'raw')
dir_manual <- here::here('data', 'manual')
dir_ensure(c(dir_figs,
             dir_derived,
             dir_manual,
             dir_raw))

#d <- read_csv(here(dir_derived, "L3_western_conservative_forest_usfs.csv"))
d <- read_csv(here(dir_derived, "L3_western_wffrc_forest_usfs.csv"))

dd <- d |>
  group_by(na_l3name, us_l3name) |>
  summarize(area_km2_both = sum(area_km2_both),
            #area_km2_conservative_forest = sum(area_km2_conservative_forest),
            area_km2_wffrc_forest = sum(area_km2_wffrc_forest),
            area_km2_total_in_west = sum(area_km2_total_in_west)) |>
  ungroup() |>
  mutate(perc_both = (area_km2_both / area_km2_total_in_west) * 100,
         #perc_forest = (area_km2_conservative_forest / area_km2_total_in_west) * 100,) |>
         perc_forest = (area_km2_wffrc_forest / area_km2_total_in_west) * 100,) |>
  arrange(desc(perc_both)) |>
  select(na_l3name, us_l3name, perc_both, perc_forest)

write_csv(dd, here(dir_derived, "l3_perc_usfs_forest.csv"))

# SI fig

# keep only the columns we need, sorted by USFS-managed share (highest first)
supp_tbl <- dd |>
  arrange(desc(perc_both)) |>
  select(us_l3name, perc_forest, perc_both)

ft <- flextable(supp_tbl) |>
  set_header_labels(
    us_l3name   = "U.S. Level III Ecoregion",
    perc_forest = "Percent forest",
    perc_both   = "Percent USFS-managed forest"
  ) |>
  colformat_double(j = c("perc_forest", "perc_both"), digits = 1) |>
  theme_booktabs() |>       # clean top/mid/bottom rules, typical journal look
  align(j = 2:3, align = "right", part = "all") |>
  font(fontname = "Times New Roman", part = "all") |>
  fontsize(size = 10, part = "all") |>
  autofit()

# grey out ecoregions where less than 30% of forest is USFS-managed
ft <- color(ft, i = ~ perc_both < 30, color = "gray50", part = "body")
ft <- bg(ft, bg = "white", part = "all")
# trim vertical padding (default is ~5 pt top and bottom) and tighten line height
ft <- ft |>
  padding(padding.top = 2, padding.bottom = 2, part = "all") |>
  line_spacing(space = 1, part = "all") |>
  autofit()  # re-fit so the new heights are used

ft  # preview in the viewer

# write out for the supplement
save_as_image(ft, path = here(dir_figs, "Table_ecoregion_forest.png"), res = 300)

