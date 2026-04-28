
rm(list = ls())

if (!requireNamespace("here", quietly = TRUE)) install.packages("here")
library(here)

source(here::here("utils", "functions.R"))
source(here::here("code", "portable_sunab3.R"))


install_load_packages(c(
  "tidyverse",
  "glue",
  "arrow",
  "here",
  "sf"
))

# ── Directories ───────────────────────────────────────────────────────────────
# Point these at the same dir_results used in run_analysis.R.

run_name    <- "GEE_resilience_v6_operational_ss500_ts50000"
version     <- "v6"
dir_results <- here::here("results", version)   # shared dir_out from run_analysis.R
dir_figs    <- here::here("figs", run_name, version)

dir_ensure(c(dir_results, dir_figs))

seed <- 1234
set.seed(seed)


weights <- arrow::read_parquet(here::here(dir_results, "weights/by_run/all_ecoregions__b10_pdsin3t1__glm_ato/weights.parquet"))
hist(weights$glm_ato_topoclimnfg_weights)


dats_short <- arrow::open_dataset(here("data/derived/GEE_resilience_v6_operational_ss500_ts50000/parquet_short_filtered")) |> 
  collect() |> 
  select(lat, long, pt_id, fire)


dats_short <- dats_short |> left_join(weights) |> select(lat, long, glm_ato_topoclimnfg_weights, fire)
dats_short_controls <- dats_short |> filter(fire == 0)



weighting_geographic_figs <- function(dats, cd_set, version) {
  # total weight for normalizing cells to relative weight share
  total_weight <- sum(dats$glm_ato_topoclimnfg_weights)
  
  ggplot(dats, aes(x = long, y = lat, z = glm_ato_topoclimnfg_weights)) +
    stat_summary_hex(
      fun = function(w) sum(w) / total_weight,  # each hex = its share of total weight mass
      bins = 60
    ) +
    scale_fill_viridis_c(
      name = "Weight share",
      labels = scales::label_percent(accuracy = 0.01),
      limits = c(0, 0.01)
    ) +
    coord_fixed(ratio = 1.3) +  # roughly corrects for lat/lon distortion at ~32°N
    labs(
      title = glue("Effective post-weighting spatial domain: {cd_set}, {version}"),
      x = "Longitude", y = "Latitude"
    ) +
    theme_minimal()
  ggsave(here(dir_figs, glue("post_glm_weightmass_{cd_set}_{version}.png")))
  
  ggplot(dats, aes(x = long, y = lat)) +
    stat_bin_hex(
      bins = 60,
      aes(fill = after_stat(count / sum(count)))  # each hex = share of total points
    ) +
    scale_fill_viridis_c(
      name = "Point share",
      labels = scales::label_percent(accuracy = 0.01),
      limits = c(0, 0.01)
    ) +
    coord_fixed(ratio = 1.3) +
    labs(
      title = glue("Effective pre-weighting spatial domain: {cd_set}, {version}"),
      x = "Longitude", y = "Latitude"
    ) +
    theme_minimal()
  ggsave(here(dir_figs, glue("pre_glm_weightmass_{cd_set}_{version}.png")))
  
  
  
  
  total_weight <- sum(dats$glm_ato_topoclimnfg_weights)
  total_pts <- length(dats$glm_ato_topoclimnfg_weights)
  
  ggplot(dats, aes(x = long, y = lat, z = glm_ato_topoclimnfg_weights)) +
    stat_summary_hex(
      fun = function(w) (sum(w) / total_weight) - (length(w) / total_pts),  # each hex = its share of total weight mass
      bins = 60
    ) +
    scale_fill_gradientn(
      values = scales::rescale(c(-1, -0.6, -0.3, -0.1, -0.05, 0, 0.05, 0.1, 0.3, 0.6, 1)),
      colors = c(
        "#2166ac", "#4393c3", "#92c5de", "#d1e5f0", "#ddeeff",  # blue side
        "white",
        "#ffffcc", "#ffeda0", "#feb24c", "#f03b20", "#bd0026"   # yellow -> orange -> red
      ),
      limits = \(x) c(-1, 1) * max(abs(x)),
      name   = "Change\n in\n weight\n share",
      labels = scales::label_percent(accuracy = 0.01)
    ) +
    coord_fixed(ratio = 1.3) +  # roughly corrects for lat/lon distortion at ~32°N
    labs(
      title = glue("Shift in spatial domain (pre-to-post-weighted): {cd_set}, {version}"),
      x = "Longitude", y = "Latitude"
    ) +
    theme_minimal() +
    theme(panel.background = element_rect(fill = "gray90", color = NA))
  ggsave(here(dir_figs, glue("shift_glm_weightmass_{cd_set}_{version}.png")))
  
  
  
  
  ggplot(dats, aes(x = long, y = lat, z = glm_ato_topoclimnfg_weights)) +
    stat_summary_hex(
      fun = function(w) ((sum(w) / total_weight) - (length(w) / total_pts)) / (length(w) / total_pts),
      bins = 60
    ) +
    scale_fill_gradientn(
      colors = c(
        "#2166ac", "#4393c3", "#92c5de", "#d1e5f0", "#ddeeff",
        "white",
        "#ffffcc", "#ffeda0", "#feb24c", "#f03b20", "#bd0026"
      ),
      values = scales::rescale(c(-1, -0.6, -0.3, -0.1, -0.05, 0, 0.5, 1, 2.0, 4.0, 8),
                               from = c(-1, 8)),
      limits = c(-0.9, 8),
      name   = "Change\n in\n weight\n share\n normalized",
      labels = scales::label_percent(accuracy = 1)
    ) +
    coord_fixed(ratio = 1.3) +  # roughly corrects for lat/lon distortion at ~32°N
    labs(
      title = glue("Shift in spatial domain (pre-to-post-weighted): {cd_set}, {version}"),
      x = "Longitude", y = "Latitude"
    ) +
    theme_minimal() +
    theme(panel.background = element_rect(fill = "gray90", color = NA))
  ggsave(here(dir_figs, glue("shift_normalized_glm_weightmass_{cd_set}_{version}.png")))
  
}


weighting_geographic_figs(dats = dats_short,
                          cd_set = "b10_pdsi3t1",
                          version = "all data")


weighting_geographic_figs(dats = dats_short_controls,
                          cd_set = "b10_pdsi3t1",
                          version = "controls_only")


