# rename_shorten_paths.R
# Shorten file/directory names and parquet/csv content under a results directory.
#
# Replacement order is fixed — extended_strict must precede extended to avoid
# partial substitution (extended_strict → ext_strict instead of extstr).
# ─────────────────────────────────────────────────────────────────────────────

# ── Configuration ─────────────────────────────────────────────────────────────

base_dir <- "C:/Users/tymc5571/dev/compound-disturbance-resilience/results-exo/v8"
dry_run  <- FALSE   # flip to FALSE to apply changes

required_pkgs <- c("arrow", "readr", "tools")
invisible(lapply(required_pkgs, function(p) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
  library(p, character.only = TRUE)
}))


# ── Substitution rules (ORDER MATTERS) ────────────────────────────────────────

rules <- list(
  c("sunab_twfe_",               "sunab_"),
  c("extended_strict",           "extstr"),   # must precede "extended"
  c("extended",                  "ext"),
  c("b90global_pdsisum90global", "b90d90")
)

apply_rules <- function(x) {
  for (r in rules) x <- gsub(r[1], r[2], x, fixed = TRUE)
  x
}

apply_rules_vec <- function(v) {
  # handles character vectors; returns unchanged if not character
  if (!is.character(v)) return(v)
  apply_rules(v)
}


# ── Logging ───────────────────────────────────────────────────────────────────

change_log <- list()

log_change <- function(type, old, new, detail = NULL) {
  tag  <- if (dry_run) "[DRY RUN]" else "[DONE]   "
  desc <- if (is.null(detail)) "" else paste0("  (", detail, ")")
  message(sprintf("%s %-8s  %s\n%s             --> %s%s",
                  tag, type, old, strrep(" ", nchar(tag) + 10), new, desc))
  change_log[[length(change_log) + 1L]] <<- list(type = type, old = old, new = new)
}


# ── Step 1: Rename directories (deepest first) ────────────────────────────────

message("\n══ Step 1: Directory renames ══════════════════════════════════════")

all_dirs <- list.dirs(base_dir, full.names = TRUE, recursive = TRUE)
all_dirs <- setdiff(all_dirs, base_dir)

# sort deepest first by counting path separators
depth <- nchar(gsub("[^/\\\\]", "", all_dirs))
all_dirs <- all_dirs[order(-depth)]

for (d in all_dirs) {
  old_nm <- basename(d)
  new_nm <- apply_rules(old_nm)
  if (old_nm == new_nm) next
  
  new_path <- file.path(dirname(d), new_nm)
  log_change("DIR", d, new_path)
  if (!dry_run) {
    ok <- file.rename(d, new_path)
    if (!ok) warning("  Rename failed: ", d)
  }
}


# ── Step 2: Rename files ──────────────────────────────────────────────────────

message("\n══ Step 2: File renames ═══════════════════════════════════════════")

# re-list after dir renames (no-op in dry_run since dirs weren't moved)
all_files <- list.files(base_dir, full.names = TRUE, recursive = TRUE)

for (f in all_files) {
  old_nm <- basename(f)
  new_nm <- apply_rules(old_nm)
  if (old_nm == new_nm) next
  
  new_path <- file.path(dirname(f), new_nm)
  log_change("FILE", f, new_path)
  if (!dry_run) {
    ok <- file.rename(f, new_path)
    if (!ok) warning("  Rename failed: ", f)
  }
}


# ── Step 3: Replace string content in .parquet and .csv files ─────────────────

message("\n══ Step 3: File content ═══════════════════════════════════════════")

# re-list after file renames
all_files  <- list.files(base_dir, full.names = TRUE, recursive = TRUE)
data_files <- all_files[grepl("\\.(parquet|csv)$", all_files, ignore.case = TRUE)]

for (f in data_files) {
  ext <- tolower(tools::file_ext(f))
  
  tryCatch({
    
    # read
    df <- if (ext == "parquet") arrow::read_parquet(f) else readr::read_csv(f, show_col_types = FALSE)
    
    # identify string/factor columns and apply rules
    changed_cols <- character(0)
    df_new <- df
    
    for (col in names(df)) {
      v <- df[[col]]
      
      if (is.character(v)) {
        v_new <- apply_rules(v)
        if (!identical(v, v_new)) {
          df_new[[col]] <- v_new
          n <- sum(v != v_new, na.rm = TRUE)
          changed_cols <- c(changed_cols, sprintf("%s (%d cells)", col, n))
        }
        
      } else if (is.factor(v)) {
        # replace in levels only, then reconstruct
        old_lvls <- levels(v)
        new_lvls <- apply_rules(old_lvls)
        if (!identical(old_lvls, new_lvls)) {
          levels(df_new[[col]]) <- new_lvls
          n <- sum(old_lvls != new_lvls)
          changed_cols <- c(changed_cols, sprintf("%s (%d levels)", col, n))
        }
      }
    }
    
    if (length(changed_cols) == 0L) next
    
    log_change("CONTENT", f, f,
               detail = paste(changed_cols, collapse = ", "))
    
    if (!dry_run) {
      if (ext == "parquet") {
        # Arrow memory-maps the source file on read; writing back to the same
        # path while the map is open triggers Windows error 1224.
        # Write to a temp file first, then atomically replace the original.
        tmp <- tempfile(fileext = ".parquet")
        arrow::write_parquet(df_new, tmp)
        file.rename(tmp, f)
      } else {
        readr::write_csv(df_new, f)
      }
    }
    
  }, error = function(e) {
    message(sprintf("[ERROR]   %s\n          %s", basename(f), e$message))
  })
}


# ── Step 4: Patch string fields in run_spec.rds files ────────────────────────

message("\n══ Step 4: run_spec.rds fields ════════════════════════════════════")

# fields to patch — all expected to be single character strings
rds_fields <- c("run_id", "run_stub", "subset_id", "group_id", "model_id", "group_col")

all_files <- list.files(base_dir, full.names = TRUE, recursive = TRUE)
rds_files <- all_files[basename(all_files) == "run_spec.rds"]

for (f in rds_files) {
  tryCatch({
    spec <- readRDS(f)
    
    # find which fields exist, have a value, and would change
    changed_fields <- character(0)
    spec_new <- spec
    
    for (field in rds_fields) {
      if (!field %in% names(spec)) next
      v <- spec[[field]]
      if (!is.character(v) || length(v) == 0) next
      v_new <- apply_rules(v)
      if (!identical(v, v_new)) {
        spec_new[[field]] <- v_new
        changed_fields <- c(changed_fields,
                            sprintf("%s: '%s' --> '%s'", field, v, v_new))
      }
    }
    
    if (length(changed_fields) == 0L) next
    
    log_change("RDS", f, f,
               detail = paste(changed_fields, collapse = "; "))
    
    if (!dry_run) saveRDS(spec_new, f)
    
  }, error = function(e) {
    message(sprintf("[ERROR]   %s\n          %s", basename(f), e$message))
  })
}


# ── Summary ───────────────────────────────────────────────────────────────────

n_dir  <- sum(sapply(change_log, `[[`, "type") == "DIR")
n_file <- sum(sapply(change_log, `[[`, "type") == "FILE")
n_cont <- sum(sapply(change_log, `[[`, "type") == "CONTENT")
n_rds  <- sum(sapply(change_log, `[[`, "type") == "RDS")

message(sprintf(
  "\n══ Summary (%s) ════════════════════════════════════════\n  Dirs renamed:              %d\n  Files renamed:             %d\n  Files with content changes:%d\n  run_spec.rds patched:      %d\n  Total changes:             %d",
  if (dry_run) "DRY RUN — no changes written" else "CHANGES APPLIED",
  n_dir, n_file, n_cont, n_rds, n_dir + n_file + n_cont + n_rds
))