
Sys.setenv(LD_LIBRARY_PATH = paste(
  "/opt/conda/lib",
  Sys.getenv("LD_LIBRARY_PATH"),
  sep = ":"
))

Sys.setenv(PATH = paste("/usr/bin:/bin:/usr/local/bin", Sys.getenv("PATH"), sep = ":"))
Sys.setenv(PKG_CONFIG_PATH = "/usr/lib/x86_64-linux-gnu/pkgconfig")


# ══════════════════════════════════════════════════════════════════════════════
# 1. Environment setup
# ══════════════════════════════════════════════════════════════════════════════

rm(list = ls())

if (!requireNamespace("here", quietly = TRUE)) install.packages("here")
library(here)
here::i_am("code/06_dydid_v6.ipynb")

required_script_pkgs <- c(
  "dplyr", "ggplot2", "tidyr", "arrow"
)


missing_script_pkgs <- required_script_pkgs[
  !vapply(required_script_pkgs, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_script_pkgs) > 0) install.packages(missing_script_pkgs)

library(dplyr)
library(ggplot2)
library(tidyr)
library(arrow)



source(here::here("code", "portable_sunab3.R"))

seed <- 1234
set.seed(seed)


# ── Directory layout ──────────────────────────────────────────────────────────
# One shared dir_results per project. Multiple run_experiment() calls over time
# all write into the same directory; rebuild_*() merges them later.

run_name <- "GEE_resilience_v6_operational_ss500_ts50000"
version <- "v6"

cyverse <- TRUE

if (cyverse) {
  dir_base    <- file.path(
    "~/data-store/data/iplant/home/shared/earthlab/macrosystems/tlm", run_name
  )
  dir_data    <- file.path(dir_base, "data", "derived")
  dir_raw     <- file.path(dir_base, "data", "raw")
  dir_manual  <- file.path(dir_base, "data", "manual")
  dir_results <- file.path(dir_base, "results", version)
  dir_figs    <- file.path(dir_base, "figs", version)
} else {
  dir_data    <- here::here("data", "derived", run_name)
  dir_raw     <- here::here("data", "raw")
  dir_manual  <- here::here("data", "manual")
  dir_results <- here::here("results", version)
  dir_figs    <- here::here("figs", version)
}

# Long (panel) data — one row per unit × year, used for estimation
dir_parquet_long <- file.path(dir_data, "parquet_long_filtered")

# Short (cross-sectional) data — one row per unit, used for weighting.
# Set to NULL and omit short_data_source in subset specs when not weighting.
dir_parquet_short <- file.path(dir_data, "parquet_short_filtered")

dir_ensure_local(c(dir_data, dir_parquet_long, dir_raw, dir_manual, dir_results, dir_figs))



dats <- arrow::open_dataset(dir_parquet_long) |> collect()

hist(dats$year)

nb <- dats |> filter(is.na(biotic_relaxedforestnorm))
nb_f <- nb |> filter(fire == 1)
length(unique(nb$pt_id))

hist(nb$year)


nd <- dats |> filter(is.na(pdsi_annual_5_yrs_prior_threshold_n4_yot))
nd_f <- nd |> filter(fire == 1)
length(unique(nd$pt_id))

hist(nd$year)





