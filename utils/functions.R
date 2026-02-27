

# Operational ----


# Function to merge biotic data
merge_biotic_year <- function(tifs) {
  if (length(tifs) == 0) stop("No raster files provided.")
  
  year_totals <- terra::rast(tifs[1])
  values(year_totals) <- 0
  
  for(tif in tifs) {
    r <- terra::rast(tif)
    
    if (!terra::compareGeom(year_totals, r, stopOnError = FALSE)) {
      stop("Rasters do not have matching dimensions/resolution/projection.")
    }
    
    year_totals <- year_totals + r
  } 
  
  return(year_totals)
}


# Data ----


#' Download and Unzip a File
#'
#' Downloads a ZIP file from a specified URL and extracts its contents to a specified directory.
#' Optionally, the ZIP file can be retained after extraction.
#'
#' @param url Character. The URL of the ZIP file to download.
#' @param extract_to Character. The directory where the contents should be extracted.
#' @param keep_zip Logical. If `TRUE`, retains the ZIP file after extraction. Defaults to `FALSE`.
#'
#' @return Invisible `NULL`. The function is used for its side effects of downloading and extracting files.
#'
#' @details The function downloads a ZIP file from a URL and extracts its contents to a specified directory.
#' If `keep_zip` is set to `FALSE`, the ZIP file will be deleted after extraction.
#'
#' @importFrom utils download.file unzip
#' @export
#'
#' @examples
#' \dontrun{
#' download_unzip_file("https://example.com/data.zip", "path/to/extract", keep_zip = TRUE)
#' }
download_unzip_file <- function(url, extract_to, keep_zip = FALSE) {
  # Validate URL and extraction path
  if (!is.character(url) || length(url) != 1) stop("`url` must be a single character string.")
  if (!is.character(extract_to) || length(extract_to) != 1) stop("`extract_to` must be a single character string.")
  if (!is.logical(keep_zip) || length(keep_zip) != 1) stop("`keep_zip` must be a single logical value.")
  
  # Ensure the extraction directory exists
  if (!dir.exists(extract_to)) dir.create(extract_to, recursive = TRUE)
  
  # Determine the path to save the ZIP file
  zip_path <- if (keep_zip) {
    # Save the ZIP file to the specified extraction directory
    file.path(extract_to, basename(url))
  } else {
    # Use a temporary file path for the ZIP file
    tempfile(fileext = ".zip")
  }
  
  # Ensure temporary file cleanup if there's an error and keep_zip is FALSE
  on.exit({
    if (!keep_zip && file.exists(zip_path)) {
      unlink(zip_path)
    }
  }, add = TRUE)
  
  # Attempt to download the ZIP file
  tryCatch({
    download.file(url, zip_path, mode = "wb")
  }, error = function(e) {
    stop("Failed to download the file from the specified URL: ", e$message)
  })
  
  # Attempt to unzip the file to the specified extraction directory
  tryCatch({
    unzip(zip_path, exdir = extract_to)
  }, error = function(e) {
    stop("Failed to unzip the file: ", e$message)
  })
  
  # Delete the ZIP file if 'keep_zip' is FALSE
  if (!keep_zip) {
    unlink(zip_path)
  }
  
  gc()
  
  invisible(NULL)
}



#' Access EPA Level III Ecoregions Data via VSI
#'
#' This function retrieves the U.S. EPA Level III ecoregions shapefile from a remote server via VSI (Virtual Spatial Infrastructure).
#' The shapefile is stored in a ZIP file, and the function accesses it without downloading the file locally.
#'
#' @return A `sf` (simple features) object containing the EPA Level III ecoregions shapefile data.
#' 
#' @details
#' The function accesses the EPA Level III ecoregions shapefile directly from the EPA's data commons, utilizing the `/vsizip/` 
#' and `/vsicurl/` mechanisms to stream the shapefile from the zipped file. The file is accessed via a URL without the need to 
#' download it locally. This method allows efficient access to the shapefile data using the `sf` package.
#'
#' @source
#' U.S. EPA Ecoregions Data: \url{https://gaftp.epa.gov/EPADataCommons/ORD/Ecoregions/us/}
#' 
#' @references
#' U.S. EPA Ecoregions Information: \url{https://www.epa.gov/eco-research/ecoregions-north-america}
#'
#' @importFrom sf st_read
#' @export
#' @examples
#' # Example usage
#' epa_ecoregions <- access_data_epa_l3_ecoregions_vsi()
#'
access_data_epa_l3_ecoregions_vsi <- function() {
  epa_l3 <- paste0(
    "/vsizip/vsicurl/",
    "https://dmap-prod-oms-edc.s3.us-east-1.amazonaws.com/ORD/Ecoregions/us/us_eco_l3.zip",
    "/us_eco_l3.shp"
  ) |>
    sf::st_read()
  
  return(epa_l3)
}

#' Read and Combine CSVs from Google Drive with Safe Column Binding
#'
#' Searches for CSV files in a Google Drive folder that match a pattern in their filename,
#' reads them using `read_csv_from_gdrive_v2()`, and combines them with `bind_rows()` (fills in NAs).
#'
#' @param drive_folder A character string. Google Drive folder name or ID.
#' @param pattern A character string. Substring or regex pattern to match in filenames.
#'
#' @return A data frame combining all matching files. Columns not present in all files are filled with NA.
#' @export
#'
#' @importFrom googledrive drive_ls
#' @importFrom dplyr filter bind_rows
#'
#' @examples
#' \dontrun{
#' df <- read_matching_csvs_from_gdrive("GEE_Exports", "2023")
#' }
# 
# read_matching_csvs_from_gdrive <- function(drive_folder, pattern, force_cols = NA, force_types = NA) {
#   # List all files in the folder
#   folder_contents <- googledrive::drive_ls(path = drive_folder)
#   
#   # Filter files by name pattern
#   matching_files <- folder_contents |>
#     dplyr::filter(grepl(pattern, name))
#   
#   if (nrow(matching_files) == 0) {
#     stop("No files found in folder '", drive_folder, "' matching pattern: ", pattern)
#   } else {
#     message("Found ", nrow(matching_files), " matching file(s).")
#   }
#   
#   # Safely read each file, collect into a list
#   data_list <- lapply(matching_files$name, function(fname) {
#     tryCatch({
#       read_csv_from_gdrive_v2(drive_folder, fname)
#     }, error = function(e) {
#       warning("Failed to read file: ", fname, " - ", e$message)
#       NULL
#     })
#   })
#   
#   # Filter out any NULLs from failed reads
#   data_list <- Filter(Negate(is.null), data_list)
#   
#   if(!is.na(force_cols)) {
#     
#     
#   }
#   
#   # Combine safely: dplyr::bind_rows fills missing columns with NA
#   combined_df <- dplyr::bind_rows(data_list)
#   
#   return(combined_df)
# }

read_matching_csvs_from_gdrive <- function(drive_folder, pattern, force_cols = NA, force_types = NA) {
  # List all files in the folder
  folder_contents <- googledrive::drive_ls(path = drive_folder)
  
  # Filter files by name pattern
  matching_files <- folder_contents |>
    dplyr::filter(grepl(pattern, name))
  
  if (nrow(matching_files) == 0) {
    stop("No files found in folder '", drive_folder, "' matching pattern: ", pattern)
  } else {
    message("Found ", nrow(matching_files), " matching file(s).")
  }
  
  # Safely read each file, collect into a list
  data_list <- lapply(matching_files$name, function(fname) {
    tryCatch({
      read_csv_from_gdrive_v2(drive_folder, fname)
    }, error = function(e) {
      warning("Failed to read file: ", fname, " - ", e$message)
      NULL
    })
  })
  
  # Filter out any NULLs from failed reads
  data_list <- Filter(Negate(is.null), data_list)
  
  if (!all(is.na(force_cols))) {
    if (all(is.na(force_types))) {
      stop("If 'force_cols' is provided, you must also provide 'force_types'.")
    }
    if (!is.character(force_cols) || !is.character(force_types)) {
      stop("'force_cols' and 'force_types' must both be character vectors.")
    }
    if (length(force_cols) != length(force_types)) {
      stop("'force_cols' and 'force_types' must be the same length.")
    }
    
    coerce_one <- function(df, col, typ) {
      # ensure column exists so bind_rows doesn't surprise you later
      if (!col %in% names(df)) df[[col]] <- NA
      
      typ <- tolower(typ)
      
      if (typ == "character") {
        df[[col]] <- as.character(df[[col]])
      } else if (typ %in% c("double", "numeric")) {
        df[[col]] <- suppressWarnings(as.numeric(df[[col]]))
      } else if (typ == "integer") {
        df[[col]] <- suppressWarnings(as.integer(df[[col]]))
      } else if (typ == "logical") {
        df[[col]] <- as.logical(df[[col]])
      } else if (typ == "date") {
        df[[col]] <- as.Date(df[[col]])
      } else if (typ == "posixct") {
        df[[col]] <- as.POSIXct(df[[col]])
      } else {
        stop("Unsupported force_types value: '", typ, "' for column '", col, "'.")
      }
      
      df
    }
    
    data_list <- lapply(data_list, function(df) {
      for (i in seq_along(force_cols)) {
        df <- coerce_one(df, force_cols[[i]], force_types[[i]])
      }
      df
    })
  }
  
  # Combine safely: dplyr::bind_rows fills missing columns with NA
  combined_df <- dplyr::bind_rows(data_list)
  
  return(combined_df)
}



#' Read a CSV File from Google Drive
#'
#' Downloads and reads a CSV file from a specified Google Drive folder using the file name.
#' This function is SIGNIFICANTLY faster than the original function
#'
#' This function locates a file by name within a Google Drive folder, retrieves its file ID,
#' reads its content as a string, and converts it into a data frame.
#'
#' @param drive_folder A character string. The path to the Google Drive folder (can be a name or ID).
#' @param file_name A character string. The exact name of the file to be read from the Drive folder.
#'
#' @return A data frame containing the contents of the CSV file.
#' @export
#'
#' @importFrom googledrive drive_ls drive_get drive_read_string as_id
#' @importFrom dplyr filter pull
#'
#' @examples
#' \dontrun{
#' df <- read_csv_from_gdrive_v2("my-folder", "my-file.csv")
#' }
read_csv_from_gdrive_v2 <- function(drive_folder, file_name) {
  # Get file ID
  folder_contents <- googledrive::drive_ls(path = drive_folder)
  id <- folder_contents |>
    dplyr::filter(name == file_name) |>
    dplyr::pull(id)
  
  if (length(id) == 0) {
    stop("File not found in specified Drive folder.")
  } else {
    message("File found: ", file_name)
  }
  
  # Get file metadata
  f <- googledrive::drive_get(googledrive::as_id(id))
  
  # Read the content of the file as a string and convert to data frame
  csv <- f |>
    googledrive::drive_read_string() %>%
    read.csv(text = .)
  
  return(csv)
}


# Utility ----

#' Generate a Timestamp
#'
#' This function generates a timestamp in various formats representing the current date and time along with the time zone.
#'
#' @param type Either 'human_read' or 'for_file'
#' @return A character string representing the current date and time in the selected
#' @examples
#' # Generate a timestamp
#' timestamp()
#'
#' @export
timestamp <- function(type = 'human_read'){
  if(type == 'human_read') {
    t <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")    
  }
  if(type == 'for_file') {
    t <- format(Sys.time(), "%Y%m%d%H%M%S")
  }
  
  
  return(t)
}


#' Write Shapefile to a New Directory and Create a Zipped Version
#'
#' This function writes an `sf` object to a shapefile in a new, file-specific directory and optionally creates a zipped version of the shapefile.
#' It also allows for the removal of the original unzipped files and handles overwriting existing files.
#'
#' @param shp An `sf` object to write as a shapefile.
#' @param location A character string specifying the path of the directory to create the new file-specific subdirectory in.
#' @param filename A character string specifying the name of the file without the `.shp` extension.
#' @param zip_only A logical value indicating whether the original (unzipped) files should be removed after zipping. Defaults to `FALSE`.
#' @param overwrite A logical value indicating whether existing files should be overwritten. Defaults to `FALSE`.
#' @return No return value. The function writes a shapefile to a specified directory, optionally zips the files, and manages file cleanup based on user input.
#' @examples
#' \dontrun{
#' # Example usage
#' st_write_shp(shp = prepped_for_parks_etal,
#'              location = here::here("data/derived"),
#'              filename = "career_lba_for_parks_v1",
#'              zip_only = TRUE,
#'              overwrite = TRUE)
#' }
#' @importFrom sf st_write
#' @importFrom zip zip
#' @export
st_write_shp <- function(shp, location, filename, zip_only = FALSE, overwrite = FALSE) {
  
  # Define paths
  out_dir <- file.path(location, filename)
  zip_file <- file.path(out_dir, paste0(filename, ".zip"))
  zip_file_dest <- file.path(location, paste0(filename, ".zip"))
  
  # Manage overwriting and directory creation
  if (dir.exists(out_dir)) {
    if (overwrite) {
      unlink(out_dir, recursive = TRUE)
    } else {
      stop("Directory '", out_dir, "' already exists and overwrite is set to FALSE.")
    }
  }
  
  if (file.exists(zip_file_dest) && zip_only) {
    if (overwrite) {
      unlink(zip_file_dest)
    } else {
      stop("Zip file '", zip_file_dest, "' already exists and overwrite is set to FALSE.")
    }
  }
  
  # Create the directory if not there
  dir_ensure(out_dir)
  
  # Write the shapefile
  shapefile_path <- file.path(out_dir, paste0(filename, ".shp"))
  sf::st_write(shp, shapefile_path, append = FALSE)
  
  # Get all shapefile components
  all_shp_files <- list.files(out_dir, pattern = paste0(filename, ".*"), full.names = TRUE)
  
  # Create zip file
  zip::zip(zipfile = zip_file, files = all_shp_files, mode = "cherry-pick")
  
  # Remove raw files if zip_only is TRUE
  if (zip_only) {
    file.copy(zip_file, zip_file_dest)
    unlink(out_dir, recursive = TRUE)
  }
}

#' Remove Columns Based on One or More Prefixes
#'
#' Removes columns from a `data.frame` or `data.table` whose names begin with 
#' any of the specified prefix strings.
#'
#' @param df A `data.frame` or `data.table`.
#' @param prefixes A character vector of prefix strings. Columns whose names start 
#' with any of these prefixes will be removed.
#'
#' @return The input data with matching columns removed, preserving the original class.
#'
#' @examples
#' df <- data.frame(
#'   a_1 = 1:3,
#'   a_2 = 4:6,
#'   b_1 = 7:9,
#'   keep_col = 10:12
#' )
#' remove_columns_with_prefix(df, prefixes = c("a_", "b_"))
#'
#' @export
remove_columns_with_prefix <- function(df, prefixes) {
  pattern <- paste0("^(", paste0(prefixes, collapse = "|"), ")")
  df <- df[, !grepl(pattern, names(df)), drop = FALSE]
  return(df)
}

#' Ensure Directories Exist
#'
#' This function checks if one or more directories exist at the specified paths,
#' and creates any that do not exist.
#'
#' @param path A character string or a vector of strings specifying directory paths.
#' @return A character vector of all directory paths that were checked/created.
#' @examples
#' # Ensure a single directory
#' dir_ensure("data")
#'
#' # Ensure multiple directories
#' dir_ensure(c("data", "output", "logs"))
#'
#' @export
dir_ensure <- function(path) {
  if (!is.character(path)) {
    stop("`path` must be a character string or a vector of character strings.")
  }
  
  created_paths <- character()
  
  for (p in path) {
    if (!dir.exists(p)) {
      tryCatch({
        dir.create(p, recursive = TRUE)
        message("Directory created: ", p)
        created_paths <- c(created_paths, p)
      }, error = function(e) {
        warning("Failed to create directory: ", p, " — ", conditionMessage(e))
      })
    } else {
      message("Directory already exists: ", p)
    }
  }
  
  return(invisible(path))
}

#' Install and Load R Packages with Optional Version Control
#'
#' Ensures that specified R packages are installed and loaded, optionally using
#' \pkg{pak} for efficient installation or \pkg{groundhog} for reproducible
#' versioned installations. The function handles missing packages, installs
#' dependencies, sets repository sources, and provides informative messages
#' about installation and loading status.
#'
#' @param pkgs Character vector of package names to install and load.
#' @param date Optional character string in the format \code{"YYYY-MM-DD"} specifying
#'   a CRAN snapshot date. If provided, the function uses the corresponding Posit Package
#'   Manager repository (\url{https://packagemanager.posit.co/}) for reproducibility.
#'   Required if \code{groundhog = TRUE}.
#' @param groundhog Logical; if \code{TRUE}, uses the \pkg{groundhog} package to
#'   load packages at the specific version available on the given date. Defaults to \code{FALSE}.
#' @param pak_quiet Logical; if \code{TRUE} (default), suppresses verbose output
#'   from \pkg{pak} during installation.
#' @param try_pak Logical; if \code{TRUE} (default), attempts to use \pkg{pak}
#'   (installing it if needed) for faster installs. If \code{FALSE}, skips
#'   all \pkg{pak} logic and uses the base fallback installer only.
#'
#' @importFrom utils install.packages packageVersion installed.packages
#' @importFrom base requireNamespace library
#' @export
install_load_packages <- function(
    pkgs,
    date = NULL,
    groundhog = FALSE,
    pak_quiet = TRUE,
    try_pak = TRUE
) {
  
  # --- Helper: quiet install with base R ---
  safe_install <- function(pkg, repos = "https://cloud.r-project.org") {
    tryCatch(
      suppressWarnings(install.packages(pkg, repos = repos, dependencies = TRUE)),
      error = function(e) message("Could not install ", pkg, ": ", e$message)
    )
  }
  
  # --- Check which packages are missing ---
  not_installed <- vapply(pkgs, function(p) !requireNamespace(p, quietly = TRUE), logical(1))
  missing_pkgs <- pkgs[not_installed]
  
  if (length(missing_pkgs) == 0) {
    message("All requested packages are already installed.")
  }
  
  # --- Optionally ensure groundhog ---
  if (groundhog) {
    if (!requireNamespace("groundhog", quietly = TRUE)) {
      message("Installing 'groundhog'...")
      safe_install("groundhog")
    }
    if (is.null(date)) {
      stop("groundhog = TRUE requires a non-null 'date' argument (YYYY-MM-DD).")
    }
  }
  
  # --- Repository selection ---
  repo <- if (!is.null(date)) {
    sprintf("https://packagemanager.posit.co/cran/%s", date)
  } else {
    "https://cloud.r-project.org"
  }
  message("Using CRAN repository: ", repo)
  options(repos = c(CRAN = repo))
  
  # --- pak availability/bootstrapping (optional) ---
  has_pak <- FALSE
  if (isTRUE(try_pak) && length(missing_pkgs) > 0) {
    has_pak <- requireNamespace("pak", quietly = TRUE)
    
    # Only try to install pak if missing and we're allowed to try pak
    if (!has_pak) {
      message("Some packages are missing; installing 'pak'...")
      
      pak_install_success <- FALSE
      try({
        suppressWarnings(
          install.packages("pak", repos = "https://cloud.r-project.org", dependencies = TRUE)
        )
        pak_install_success <- requireNamespace("pak", quietly = TRUE)
      }, silent = TRUE)
      
      if (!pak_install_success) {
        message("Standard install failed; trying pak bootstrap installer...")
        try({
          source("https://pak.r-lib.org/install.R")
          pak_install_success <- requireNamespace("pak", quietly = TRUE)
        }, silent = TRUE)
      }
      
      if (!pak_install_success) {
        warning("Failed to install 'pak' by any method; will fall back to base installers only.")
      }
      
      has_pak <- requireNamespace("pak", quietly = TRUE)
    }
  }
  
  # --- Install missing packages ---
  installed_or_updated <- FALSE
  
  if (length(missing_pkgs) > 0) {
    message("Missing packages detected: ", paste(missing_pkgs, collapse = ", "))
    
    if (isTRUE(try_pak) && has_pak) {
      tryCatch({
        if (pak_quiet) {
          message("Attempting install with pak (quietly)...")
          suppressMessages(suppressWarnings(
            pak::pkg_install(missing_pkgs, ask = FALSE, upgrade = FALSE)
          ))
        } else {
          message("Attempting install with pak...")
          pak::pkg_install(missing_pkgs, ask = FALSE, upgrade = FALSE)
        }
        installed_or_updated <- TRUE
      }, error = function(e) {
        message("pak installation failed: ", e$message)
        message("Falling back to install.packages()...")
        for (p in missing_pkgs) safe_install(p, repos = repo)
        installed_or_updated <- TRUE
      })
    } else {
      if (!isTRUE(try_pak)) {
        message("try_pak = FALSE; installing missing packages with install.packages()...")
      } else {
        message("pak unavailable; installing missing packages with install.packages()...")
      }
      for (p in missing_pkgs) safe_install(p, repos = repo)
      installed_or_updated <- TRUE
    }
  }
  
  # --- Load packages ---
  failed_to_load <- character()
  
  if (groundhog) {
    message("Loading packages with groundhog (date = ", date, ")...")
    tryCatch({
      groundhog::groundhog.library(pkgs, date = date)
    }, error = function(e) {
      message("groundhog loading error: ", e$message)
      failed_to_load <<- pkgs
    })
  } else {
    message("Loading packages...")
    for (p in pkgs) {
      ok <- tryCatch({
        library(p, character.only = TRUE, quietly = TRUE)
        TRUE
      }, error = function(e) {
        message("Failed to load ", p, ": ", e$message)
        FALSE
      })
      if (!ok) failed_to_load <- c(failed_to_load, p)
    }
  }
  
  # --- Restart message if needed ---
  if (installed_or_updated || length(failed_to_load) > 0) {
    message("\nSome packages were newly installed, updated, or failed to load.\n",
            "This may be due to updated dependencies already loaded in memory.\n",
            "Please restart R and re-run this function to ensure all packages load correctly.\n")
  }
  
  # --- Report loaded versions ---
  loaded_versions <- sapply(pkgs, function(p) {
    if (requireNamespace(p, quietly = TRUE)) as.character(packageVersion(p)) else NA_character_
  })
  message("Packages loaded:\n",
          paste(names(loaded_versions), loaded_versions, collapse = "\n"))
  invisible(loaded_versions)
}



#' Safe Unzip a File (with Optional Recursive Unzipping and ZIP Cleanup)
#'
#' Safely unzips a ZIP file to a specified directory, skipping if all expected contents already exist.
#' Optionally removes the original and/or nested ZIP files after extraction.
#'
#' @param zip_path Character. Path to the local ZIP file.
#' @param extract_to Character. Directory where the contents should be extracted. Defaults to the ZIP's directory.
#' @param recursive Logical. If TRUE, recursively unzip nested ZIP files. Defaults to FALSE.
#' @param keep_zip Logical. If FALSE, deletes the original ZIP and any nested ZIPs after unzipping. Defaults to TRUE.
#'
#' @return A character vector of full paths of the extracted files (excluding directories).
#'
#' @importFrom utils unzip
#' @export
#'
#' @examples
#' \dontrun{
#' files <- safe_unzip("data/archive.zip", recursive = TRUE, keep_zip = FALSE)
#' print(files)  # Only unzipped files, not folders
#' }
safe_unzip <- function(zip_path,
                       extract_to = dirname(zip_path),
                       recursive = FALSE,
                       keep_zip = TRUE) {
  # Validate inputs
  if (!file.exists(zip_path)) stop("ZIP file does not exist: ", zip_path)
  if (!is.character(extract_to) || length(extract_to) != 1) stop("`extract_to` must be a single character string.")
  if (!is.logical(recursive) || length(recursive) != 1) stop("`recursive` must be a single logical value.")
  if (!is.logical(keep_zip) || length(keep_zip) != 1) stop("`keep_zip` must be a single logical value.")
  
  # List expected files from the archive
  zip_listing <- unzip(zip_path, list = TRUE)
  expected_paths <- file.path(extract_to, zip_listing$Name)
  
  # Skip if already fully extracted
  if (all(file.exists(expected_paths))) {
    message("Skipping unzip: All expected files already exist in ", extract_to)
    
    # Get all unzipped files (excluding directories)
    all_files <- list.files(extract_to, recursive = TRUE, full.names = TRUE)
    file_paths <- all_files[file.info(all_files)$isdir == FALSE]
    
    return(normalizePath(file_paths, mustWork = FALSE))
    
  } else {
    if (!dir.exists(extract_to)) dir.create(extract_to, recursive = TRUE)
    tryCatch({
      unzip(zip_path, exdir = extract_to)
    }, error = function(e) {
      stop("Failed to unzip: ", e$message)
    })
    
    # Recursive unzip of nested ZIPs
    if (recursive) {
      nested_zips <- list.files(extract_to, pattern = "\\.zip$", recursive = TRUE, full.names = TRUE)
      for (nz in nested_zips) {
        unzip(nz, exdir = dirname(nz))
        if (!keep_zip) unlink(nz)
      }
    }
    
    # Optionally remove original zip
    if (!keep_zip) unlink(zip_path)
    
    # Get all unzipped files (excluding directories)
    all_files <- list.files(extract_to, recursive = TRUE, full.names = TRUE)
    file_paths <- all_files[file.info(all_files)$isdir == FALSE]
    
    return(normalizePath(file_paths, mustWork = FALSE))
  }
}



#' Safely Download a File to a Directory
#'
#' Downloads a file from a URL to a specified directory, only if it doesn't already exist there.
#'
#' @param url Character. The URL to download from.
#' @param dest_dir Character. The directory where the file should be saved.
#' @param mode Character. Mode passed to `download.file()`. Default is "wb" (write binary).
#' @param timeout Integer. Optional timeout in seconds. Will be reset afterward.
#'
#' @return A character string with the full path to the downloaded file.
#'
#' @importFrom utils download.file
#' @export
#'
#' @examples
#' \dontrun{
#' path <- safe_download("https://example.com/data.zip", "data/")
#' }
safe_download <- function(url,
                          dest_dir,
                          mode = "wb",
                          timeout = NA) {
  # Validate input
  if (!is.character(url) || length(url) != 1) stop("`url` must be a single character string.")
  if (!is.character(dest_dir) || length(dest_dir) != 1) stop("`dest_dir` must be a single character string.")
  
  # Ensure destination directory exists
  if (!dir.exists(dest_dir)) dir.create(dest_dir, recursive = TRUE)
  
  # Derive destination file path from URL and directory
  filename <- basename(url)
  destfile <- file.path(dest_dir, filename)
  
  # Skip download if file already exists
  if (file.exists(destfile)) {
    message("Skipping download: File already exists at ", destfile)
    return(normalizePath(destfile, mustWork = FALSE))
  }
  
  # Handle optional timeout
  original_timeout <- getOption("timeout")
  if (!is.na(timeout) && timeout > original_timeout) {
    options(timeout = timeout)
    on.exit(options(timeout = original_timeout), add = TRUE)
  }
  
  # Attempt to download
  tryCatch({
    download.file(url, destfile, mode = mode)
    message("Downloaded: ", destfile)
  }, error = function(e) {
    stop("Failed to download file from URL: ", e$message)
  })
  
  return(normalizePath(destfile, mustWork = FALSE))
}





#Function to clip a raster to a vector, ensuring in same projection
#Returns raster in original projection, but clipped to vector
#Returns raster in the same form that it came in
# PARAMETERS
# raster : a SpatRaster, PackedSpatRaster, RasterLayer, RasterStack, or RasterBrick object
# vector : a SpatVector, PackedSpatVector or SF object
# mask : TRUE or FALSE; whether terra::clip should mask the raster as well
crop_careful_universal <- function(raster, vector, mask = FALSE, verbose = FALSE) {
  pack <- FALSE
  
  #Unpack if parallelized inputs
  if(class(raster)[1] == "PackedSpatRaster") {
    raster <- terra::unwrap(raster)
    pack <- TRUE
  }
  if(class(vector)[1] == "PackedSpatVector") {
    vector <- sf::st_as_sf(terra::unwrap(vector))
  }
  
  #Handle unpacked spatVector
  if(class(vector)[1] == "SpatVector") {
    vector <- sf::st_as_sf(vector)
  }
  
  #If using raster package
  if(class(raster)[1] == "RasterLayer" | class(raster)[1] == "RasterStack" | class(raster)[1] == "RasterBrick") {
    
    #Perform operation
    if (raster::crs(vector) != raster::crs(raster)) { #if raster and vector aren't in same projection, change vector to match
      if(verbose) {print("Projecting vector")}
      vector <- sf::st_transform(vector, raster::crs(raster)) 
    } else {
      if(verbose) {print("Vector already in raster CRS")}
    }
    if(verbose) {print("Clipping")}
    r <- raster::crop(raster,
                      vector)
    if(mask) {
      r <- r |> raster::mask(vector)
    }
    
    return(r)
    
  } else { #terra package
    
    #Perform operation
    if (terra::crs(vector) != terra::crs(raster)) { #if raster and vector aren't in same projection, change vector to match
      if(verbose) {print("Projecting vector")}
      vector <- sf::st_transform(vector, terra::crs(raster)) 
    } else {
      if(verbose) {print("Vector already in raster CRS")}
    }
    if(verbose) {print("Clipping")}
    r <- terra::crop(raster,
                     vector,
                     mask = mask) #crop & mask
    
    #Repack if was packed coming in (i.e. parallelized)
    if(pack) {
      r <- terra::wrap(r)
    }
    return(r)
    
  }
}
