# Download USFS


rm(list = ls())

if(!requireNamespace("here", quietly = TRUE)) {
  install.packages("here")
}
library(here)

source(here::here("utils", "functions.R"))

install_and_load_packages(c(
  "tidyverse",
  "terra",
  "mapview",
  "sf",
  "janitor",
  "tictoc"))

dir_figs <- here::here('figs')
dir_derived <- here::here('data', 'derived')
dir_raw <- here::here('data', 'raw')
dir_manual <- here::here('data', 'manual')
dir_ensure(c(dir_figs,
             dir_derived,
             dir_manual,
             dir_raw))

options(timeout = 600)

download_unzip_file(url = "https://blm-egis.maps.arcgis.com/sharing/rest/content/items/6bf2e737c59d4111be92420ee5ab0b46/data",
                    extract_to = dir_raw,
                    keep_zip = FALSE)

sma_fl <- here(dir_raw, "SMA_WM.gdb")

st_layers(sma_fl)

# Filter to western USFS only 
usfs <- sf::st_read(sma_fl, layer = "SurfaceMgtAgy_USFS")

states <- tigris::states()
state_list <- c("WA", "OR", "CA", "NV", "ID", "MT", "WY", "CO", "UT", "AZ", "NM")
west <- states |>
  filter(STUSPS %in% state_list) |>
  sf::st_transform(st_crs(usfs))

west_usfs <- usfs |>
  st_cast("POLYGON", group_or_split = TRUE) |>
  sf::st_filter(west, .predicate = st_intersects)

mapview(west_usfs)


# st_write_shp(shp = west_usfs,
#              location = here::here("data/derived"),
#              filename = "sma_usfs",
#              zip_only = TRUE,
#              overwrite = TRUE)

#Rasterize

r_template <- rast(ext(west_usfs), resolution = 30)  # Adjust resolution
crs(r_template) <- st_crs(west_usfs)$wkt  # Match CRS

tic() #393 seconds
poly_vect <- vect(west_usfs)
r_poly <- rasterize(poly_vect, r_template, field = 1, background = 0)
toc()


plot(r_poly, col = c("white", "forestgreen"), legend = FALSE)

terra::writeRaster(r_poly,
                   filename = here(dir_derived, "usfs_binary.tif"),
                   overwrite = TRUE,
                   gdal = c("COMPRESS=DEFLATE")
                   )




