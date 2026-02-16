
rm(list = ls())

if(!requireNamespace("here", quietly = TRUE)) {
  install.packages("here")
}
library(here)

source(here::here("utils", "functions.R"))

install_load_packages(c(
  "tidyverse",
  "terra"))

dir_figs <- here::here('figs')
dir_derived <- here::here('data', 'derived')
dir_raw <- here::here('data', 'raw')
dir_manual <- here::here('data', 'manual')
dir_ensure(c(dir_figs,
             dir_derived,
             dir_manual,
             dir_raw))

seed = 1234
set.seed(seed)

# --- load biotic as before ---
forest_dir <- dir_derived  # where forest_masks_XXXX.tif live
biotic <- terra::rast(here(dir_derived, "biotic_test_toDrive.tif"))

years_biotic <- 1997:2023
tmpl <- biotic[[1]]

names(biotic) <- paste0("biotic_", years_biotic)

stopifnot(all(paste0("biotic_", years_biotic) %in% names(biotic)))

# Helper: read the 'relaxed' band for a given year, project to biotic CRS if needed
read_relaxed_forest <- function(year, tmpl, forest_dir) {
  fp <- here(forest_dir, sprintf("forest_masks_%d.tif", year))
  if (!file.exists(fp)) stop("Missing forest mask file: ", fp)
  
  f <- terra::rast(fp)
  
  if (!("relaxed" %in% names(f))) {
    stop("Band 'relaxed' not found in: ", fp, "\nBands are: ", paste(names(f), collapse = ", "))
  }
  
  f_relaxed <- f[["relaxed"]]
  
  # ensure CRS matches tmpl before resampling
  if (!terra::compareGeom(f_relaxed, tmpl, crs=TRUE, ext=FALSE, rowcol=FALSE, res=FALSE, stopOnError=FALSE)) {
    f_relaxed <- terra::project(f_relaxed, terra::crs(tmpl), method = "near")
  }
  
  f_relaxed
}

# Function to get forest % on the biotic grid for a given forest layer (0/1 + NA)
forest_pct_on_biotic_grid <- function(forest_layer, tmpl) {
  f01 <- terra::ifel(is.na(forest_layer), 0, forest_layer)
  p <- terra::resample(f01, tmpl, method = "average")  # proportion in [0,1]
  p * 100
}

# Build normalized biotic stack
norm_list <- vector("list", length(years_biotic))

for (i in seq_along(years_biotic)) {
  y <- years_biotic[i]
  
  b <- biotic[[paste0("biotic_", y)]]
  
  # previous year forest mask file + relaxed band
  f_prev <- read_relaxed_forest(y - 1, tmpl = tmpl, forest_dir = forest_dir)
  
  pct <- forest_pct_on_biotic_grid(f_prev, tmpl)  # 0..100
  
  denom <- pct / 100
  norm <- b / denom
  norm <- terra::ifel(denom <= 0, NA, norm)
  norm <- terra::clamp(norm, upper = 100, values = TRUE)
  
  names(norm) <- paste0("biotic_norm_", y)
  norm_list[[i]] <- norm
}

biotic_norm <- terra::rast(norm_list)

terra::writeRaster(
  biotic_norm,
  filename = here(dir_derived, "biotic_norm_relaxedforest.tif"),
  overwrite = TRUE
)



x <- rast(here(dir_derived, "biotic_norm_relaxedforest.tif"))
names(x)


