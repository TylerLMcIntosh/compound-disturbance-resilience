
rm(list = ls())

if (!requireNamespace("here", quietly = TRUE)) install.packages("here")
library(here)

source(here::here("utils", "functions.R"))
#source(here::here("code", "portable_sunab3.R"))


install_load_packages(c(
  "tidyverse",
  "glue",
  "arrow",
  "here",
  "sf",
  "h3jsr"
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




## H3JSR version ----

weighting_geographic_figs_h3 <- function(
    dats,
    cd_set,
    version,
    h3_res = 5,
    weight_col = "glm_ato_topoclimnfg_weights",
    lon_col = "long",
    lat_col = "lat",
    plot_crs = 4326
) {
  
  # Required packages:
  # library(dplyr)
  # library(sf)
  # library(ggplot2)
  # library(h3jsr)
  # library(scales)
  # library(glue)
  # library(here)
  
  # Prepare point data
  dats_sf <- dats |>
    dplyr::filter(
      !is.na(.data[[lon_col]]),
      !is.na(.data[[lat_col]]),
      !is.na(.data[[weight_col]])
    ) |>
    sf::st_as_sf(
      coords = c(lon_col, lat_col),
      crs = 4326,
      remove = FALSE
    )
  
  total_weight <- sum(dats_sf[[weight_col]], na.rm = TRUE)
  total_pts <- nrow(dats_sf)
  
  # Assign H3 cell IDs
  dats_h3 <- dats_sf |>
    dplyr::mutate(
      h3_cell = h3jsr::point_to_cell(geometry, res = h3_res)
    )
  
  # Aggregate to H3 cells
  h3_summary <- dats_h3 |>
    sf::st_drop_geometry() |>
    dplyr::group_by(h3_cell) |>
    dplyr::summarise(
      weight_sum = sum(.data[[weight_col]], na.rm = TRUE),
      n_points = dplyr::n(),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      weight_share = weight_sum / total_weight,
      point_share = n_points / total_pts,
      weight_share_shift = weight_share - point_share,
      weight_share_shift_norm = weight_share_shift / point_share
    )
  
  # Convert H3 cells to polygons
  h3_summary_sf <- h3_summary |>
    dplyr::mutate(
      geometry = h3jsr::cell_to_polygon(h3_cell)
    ) |>
    sf::st_as_sf(crs = 4326)
  
  # Optional display projection
  h3_plot_sf <- sf::st_transform(h3_summary_sf, plot_crs)
  
  # ---------------------------------------------------------------------------
  # 1. Post-weighting spatial domain: weight share
  # ---------------------------------------------------------------------------
  
  p_post <- ggplot(h3_plot_sf) +
    geom_sf(aes(fill = weight_share), color = NA) +
    scale_fill_viridis_c(
      name = "Weight share",
      labels = scales::label_percent(accuracy = 0.01),
      limits = c(0, 0.01)
    ) +
    coord_sf() +
    labs(
      title = glue::glue(
        "Effective post-weighting spatial domain: {cd_set}, {version}"
      ),
      x = "Longitude",
      y = "Latitude"
    ) +
    theme_minimal()
  
  ggsave(
    filename = here::here(
      dir_figs,
      glue::glue("post_glm_weightmass_{cd_set}_{version}.png")
    ),
    plot = p_post
  )
  
  # ---------------------------------------------------------------------------
  # 2. Pre-weighting spatial domain: point share
  # ---------------------------------------------------------------------------
  
  p_pre <- ggplot(h3_plot_sf) +
    geom_sf(aes(fill = point_share), color = NA) +
    scale_fill_viridis_c(
      name = "Point share",
      labels = scales::label_percent(accuracy = 0.01),
      limits = c(0, 0.01)
    ) +
    coord_sf() +
    labs(
      title = glue::glue(
        "Effective pre-weighting spatial domain: {cd_set}, {version}"
      ),
      x = "Longitude",
      y = "Latitude"
    ) +
    theme_minimal()
  
  ggsave(
    filename = here::here(
      dir_figs,
      glue::glue("pre_glm_weightmass_{cd_set}_{version}.png")
    ),
    plot = p_pre
  )
  
  # ---------------------------------------------------------------------------
  # 3. Absolute shift in spatial domain: post share - pre share
  # ---------------------------------------------------------------------------
  
  p_shift <- ggplot(h3_plot_sf) +
    geom_sf(aes(fill = weight_share_shift), color = NA) +
    scale_fill_gradientn(
      values = scales::rescale(
        c(-1, -0.6, -0.3, -0.1, -0.05, 0, 0.05, 0.1, 0.3, 0.6, 1)
      ),
      colors = c(
        "#2166ac", "#4393c3", "#92c5de", "#d1e5f0", "#ddeeff",
        "white",
        "#ffffcc", "#ffeda0", "#feb24c", "#f03b20", "#bd0026"
      ),
      limits = function(x) {
        c(-1, 1) * max(abs(x), na.rm = TRUE)
      },
      name = "Change\n in\n weight\n share",
      labels = scales::label_percent(accuracy = 0.01)
    ) +
    coord_sf() +
    labs(
      title = glue::glue(
        "Shift in spatial domain (pre-to-post-weighted): {cd_set}, {version}"
      ),
      x = "Longitude",
      y = "Latitude"
    ) +
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "gray90", color = NA)
    )
  
  ggsave(
    filename = here::here(
      dir_figs,
      glue::glue("shift_glm_weightmass_{cd_set}_{version}.png")
    ),
    plot = p_shift
  )
  
  # ---------------------------------------------------------------------------
  # 4. Normalized shift: (post share - pre share) / pre share
  # ---------------------------------------------------------------------------
  
  p_shift_norm <- ggplot(h3_plot_sf) +
    geom_sf(aes(fill = weight_share_shift_norm), color = NA) +
    scale_fill_gradientn(
      colors = c(
        "#2166ac", "#4393c3", "#92c5de", "#d1e5f0", "#ddeeff",
        "white",
        "#ffffcc", "#ffeda0", "#feb24c", "#f03b20", "#bd0026"
      ),
      values = scales::rescale(
        c(-1, -0.6, -0.3, -0.1, -0.05, 0, 0.5, 1, 2.0, 4.0, 8),
        from = c(-1, 8)
      ),
      limits = c(-0.9, 8),
      name = "Change\n in\n weight\n share\n normalized",
      labels = scales::label_percent(accuracy = 1)
    ) +
    coord_sf() +
    labs(
      title = glue::glue(
        "Shift in spatial domain (pre-to-post-weighted): {cd_set}, {version}"
      ),
      x = "Longitude",
      y = "Latitude"
    ) +
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "gray90", color = NA)
    )
  
  ggsave(
    filename = here::here(
      dir_figs,
      glue::glue("shift_normalized_glm_weightmass_{cd_set}_{version}.png")
    ),
    plot = p_shift_norm
  )
  
  invisible(list(
    h3_summary = h3_summary,
    h3_summary_sf = h3_summary_sf,
    post_plot = p_post,
    pre_plot = p_pre,
    shift_plot = p_shift,
    shift_norm_plot = p_shift_norm
  ))
}






library(sf)
library(dplyr)
library(ggplot2)
library(tigris)
library(h3jsr)
library(glue)

options(tigris_use_cache = TRUE)

# Western states
west_states <- tigris::states(cb = TRUE, year = 2023) |>
  sf::st_as_sf() |>
  dplyr::filter(STUSPS %in% c(
    "WA", "OR", "CA", "NV", "ID", "MT", "WY",
    "UT", "AZ", "CO", "NM"
  )) |>
  sf::st_transform(4326)

# Dissolve to one western U.S. polygon
west_union <- west_states |>
  dplyr::summarise(geometry = sf::st_union(geometry)) |>
  sf::st_make_valid()

# H3 resolution
h3_res <- 5

# Get H3 cells whose centers fall inside the western U.S. polygon
west_h3_ids <- h3jsr::polygon_to_cells(
  geometry = west_union,
  res = h3_res,
  simple = TRUE
) |>
  unlist(use.names = FALSE)

# Convert H3 cells to polygons
west_h3_sf <- tibble::tibble(h3_address = west_h3_ids) |>
  dplyr::mutate(
    geometry = h3jsr::cell_to_polygon(input = h3_address)
  ) |>
  sf::st_as_sf(crs = 4326)

# Optional: drop any NA cells, just in case
west_h3_sf <- west_h3_sf |>
  dplyr::filter(!is.na(h3_address))

# Plot in Albers for western U.S. display
west_states_5070 <- sf::st_transform(west_states, 5070)
west_h3_5070 <- sf::st_transform(west_h3_sf, 5070)

ggplot() +
  geom_sf(data = west_states_5070, fill = "gray95", color = "gray60") +
  geom_sf(data = west_h3_5070, fill = NA, color = "black", linewidth = 0.05) +
  coord_sf() +
  labs(
    title = "H3 resolution 5 hexes over the western U.S.",
    subtitle = glue::glue("{nrow(west_h3_sf)} H3 cells"),
    x = NULL,
    y = NULL
  ) +
  theme_minimal()
