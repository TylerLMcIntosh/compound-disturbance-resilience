
run_name <- "GEE_resilience_v6_operational_ss500_ts50000"

dir_base    <- file.path(
    "~/data-store/data/iplant/home/shared/earthlab/macrosystems/tlm", run_name
  )
dir_derived    <- file.path(dir_base, "data", "derived")

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
# 
# # check that filtering ends up with same thing
# dir_long  <- here(dir_derived, "parquet_long")
# dir_short <- here(dir_derived, "parquet_short")
# 
# long_file  <- list.files(dir_long, pattern = "\\.parquet$", full.names = TRUE)[1]
# short_file <- list.files(dir_short, pattern = "\\.parquet$", full.names = TRUE)[1]
# 
# long_ids <- read_parquet(long_file) |>
#   filter((year >= 1997 & burn_year >= 2002) | (year > 1997 & is.na(burn_year))) |>
#   distinct(pt_id) |>
#   arrange(pt_id)
# 
# short_ids <- read_parquet(short_file) |>
#   filter(burn_year >= 2002 | is.na(burn_year)) |>
#   distinct(pt_id) |>
#   arrange(pt_id)
# 
# n_long  <- nrow(long_ids)
# n_short <- nrow(short_ids)
# 
# ids_only_long  <- anti_join(long_ids, short_ids, by = "pt_id")
# ids_only_short <- anti_join(short_ids, long_ids, by = "pt_id")
# 
# list(
#   long_file        = basename(long_file),
#   short_file       = basename(short_file),
#   n_long_ptid      = n_long,
#   n_short_ptid     = n_short,
#   same_ids         = identical(long_ids$pt_id, short_ids$pt_id),
#   only_in_long_n   = nrow(ids_only_long),
#   only_in_short_n  = nrow(ids_only_short),
#   only_in_long     = head(ids_only_long, 20),
#   only_in_short    = head(ids_only_short, 20)
# )
# 
# 
# xxx <- read_parquet(short_file) |> filter(pt_id %in% ids_only_long$pt_id) |> select(pt_id, fire, burn_year)




# run to fix files
dir_long  <- here(dir_derived, "parquet_long")
dir_short <- here(dir_derived, "parquet_short")

dir_long_filtered  <- here(dir_derived, "parquet_long_filtered")
dir_short_filtered <- here(dir_derived, "parquet_short_filtered")

dir.create(dir_long_filtered, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_short_filtered, recursive = TRUE, showWarnings = FALSE)

long_files <- list.files(dir_long, pattern = "\\.parquet$", full.names = TRUE)

walk(long_files, function(f) {
  out_file <- file.path(dir_long_filtered, basename(f))
  
  read_parquet(f) |>
    filter((year >= 1997 & burn_year >= 2002) | (year > 1997 & is.na(burn_year))) |>
    write_parquet(out_file)
})

short_files <- list.files(dir_short, pattern = "\\.parquet$", full.names = TRUE)

walk(short_files, function(f) {
  out_file <- file.path(dir_short_filtered, basename(f))
  
  read_parquet(f) |>
    filter(burn_year >= 2002 | is.na(burn_year)) |>
    write_parquet(out_file)
})