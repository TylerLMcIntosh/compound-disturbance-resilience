cyverse = FALSE

run_name <- "GEE_resilience_v6_operational_ss500_ts50000"


if(cyverse) {
  Sys.setenv(LD_LIBRARY_PATH = paste(
    "/opt/conda/lib",
    Sys.getenv("LD_LIBRARY_PATH"),
    sep = ":"
  ))
  
  Sys.setenv(PATH = paste("/usr/bin:/bin:/usr/local/bin", Sys.getenv("PATH"), sep = ":"))
  Sys.setenv(PKG_CONFIG_PATH = "/usr/lib/x86_64-linux-gnu/pkgconfig")
  
  dir_base    <- file.path(
    "~/data-store/data/iplant/home/shared/earthlab/macrosystems/tlm", run_name
  )
  dir_derived    <- file.path(dir_base, "data", "derived")
} else {
  dir_derived <- here("data/derived", run_name)
}



required_script_pkgs <- c(
  "here", "arrow", "dplyr", "purrr"
)
missing_script_pkgs <- required_script_pkgs[
  !vapply(required_script_pkgs, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_script_pkgs) > 0) install.packages(missing_script_pkgs)


library(here)
library(arrow)
library(dplyr)
library(purrr)



# run to fix files
dir_long  <- here(dir_derived, "parquet_long")
dir_short <- here(dir_derived, "parquet_short")

dir_long_filtered  <- here(dir_derived, "parquet_long_filtered")
dir_short_filtered <- here(dir_derived, "parquet_short_filtered")

dir.create(dir_long_filtered, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_short_filtered, recursive = TRUE, showWarnings = FALSE)

long_files <- list.files(dir_long, pattern = "\\.parquet$", full.names = TRUE)
short_files <- list.files(dir_short, pattern = "\\.parquet$", full.names = TRUE)



# pwalk(list(long_files, short_files), function(lf, sf) {
#   long_out_file <- file.path(dir_long_filtered, basename(lf))
#   short_out_file <- file.path(dir_short_filtered, basename(sf))
#   
#   read_parquet(lf) |>
#     filter((year >= 1992 & burn_year >= 2002) | (year > 1992 & is.na(burn_year))) |>
#     write_parquet(long_out_file)
#   
#   read_parquet(sf) |>
#     filter(burn_year >= 2002 | is.na(burn_year)) |>
#     write_parquet(short_out_file)
# })


long_files <- sort(long_files)
short_files <- sort(short_files)

stopifnot(length(long_files) == length(short_files))

pwalk(
  list(lf = long_files, sf = short_files),
  function(lf, sf) {
    
    long_out_file <- file.path(dir_long_filtered, basename(lf))
    short_out_file <- file.path(dir_short_filtered, basename(sf))
    
    short_dat <- arrow::read_parquet(sf) |>
      dplyr::filter(burn_year >= 2002 | is.na(burn_year))
    
    short_dat$h3jsr_5 <- h3jsr::point_to_cell(
      input = as.matrix(short_dat[, c("long", "lat")]),
      res = 5
    )
    
    h3_lookup <- short_dat |>
      dplyr::select(pt_id, h3jsr_5) |>
      dplyr::distinct()

    
    long_dat <- arrow::read_parquet(lf) |>
      dplyr::filter(
        (year >= 1992 & burn_year >= 2002) |
          (year > 1992 & is.na(burn_year))
      ) |>
      dplyr::left_join(h3_lookup, by = "pt_id")
    
    arrow::write_parquet(long_dat, long_out_file)
    arrow::write_parquet(short_dat, short_out_file)
  }
)

