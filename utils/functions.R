

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

read_matching_csvs_from_gdrive <- function(drive_folder, pattern) {
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


install_and_load_packages <- function(package_list, auto_install = "n") {
  # Ensure pak is available
  if (!requireNamespace("pak", quietly = TRUE)) {
    cat("The 'pak' package is required for fast installation of packages, installing now.\n")
    install.packages("pak")
  }
  
  # Helper: Extract base name of a package for require()
  parse_pkg_name <- function(pkg) {
    if (grepl("/", pkg)) {
      sub("^.+/(.+?)(@.+)?$", "\\1", pkg)  # GitHub: extract repo name
    } else {
      sub("@.*$", "", pkg)  # CRAN: remove @version if present
    }
  }
  
  # Classify and separate packages
  missing_pkgs <- c()
  for (pkg in package_list) {
    pkg_name <- parse_pkg_name(pkg)
    if (!requireNamespace(pkg_name, quietly = TRUE)) {
      missing_pkgs <- c(missing_pkgs, pkg)
    }
  }
  
  # Install missing ones (CRAN or GitHub), with version support
  if (length(missing_pkgs) > 0) {
    pak::pkg_install(missing_pkgs, upgrade = TRUE)
  }
  
  # Load all packages
  for (pkg in package_list) {
    pkg_name <- parse_pkg_name(pkg)
    success <- require(pkg_name, character.only = TRUE, quietly = TRUE)
    if (!success) cat("Failed to load package:", pkg_name, "\n")
  }
  
  cat("All specified packages installed and loaded.\n")
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
