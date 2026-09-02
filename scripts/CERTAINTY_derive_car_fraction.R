# name: CERTAINTY_derive_car_fraction.R
# ==============================================================================
# Derive CAR-positive T-cell fraction per specimen, OUTSIDE OMOP.
# ==============================================================================
#
# One row per PBMC specimen (135) with the derived feature, provenance and QC.
# Emits NO OMOP rows; the ETL step (CERTAINTY_etl_scrna_features.R) lifts the
# result into MEASUREMENT.
#
# Feature: car_pos_fraction
#   numerator   = T cells with CAR-BCMA count > 0   (CAR_BY_EXPRS == TRUE)
#   denominator = T cells                          (celltype_short_3 ~ ^CD4|^CD8|^gd)
#   value       = numerator / denominator           (dimensionless fraction in [0,1])
#   min_cells   = 50                              (below this -> MISSING, never zero)
#
# Inputs:
#   - split manifest:    $VTI_DATA_ROOT/h5ad/splits/_split_manifest.csv  (135 specimens, sha256)
#   - per-specimen h5ad: $VTI_DATA_ROOT/h5ad/splits/<orig.ident>.h5ad  (obs only: CAR_BY_EXPRS, celltype_short_3)
#
# Outputs ($VTI_FEATURE_DIR/):
#   - car_fraction_derivation.csv         (135 rows, one per specimen)
#   - car_fraction_provenance.csv         (single-row provenance summary)
#   - car_fraction_validation_checks.csv  (one row per check)
#
# Provenance pin requires the Rade_et_al_CAR_2025 source checkout at
# $VTI_RADE_REPO if set; otherwise git revision is recorded as NA.
#
# Run inside the RStudio container (bind your data root to the path
# configured by VTI_DATA_ROOT):
#   apptainer exec --bind $PWD/data:/data singularity-rstudio-4-3-2.sif \
#     Rscript scripts/CERTAINTY_derive_car_fraction.R \
#     [VTI_DATA_ROOT=/data VTI_FEATURE_DIR=/data/output/feature_derivations]

suppressMessages({
  library(rhdf5)
  library(dplyr)
})

# ==============================================================================
# PATHS + CONSTANTS  (env-var driven; defaults are repo-relative)
# ==============================================================================

DATA_ROOT    <- Sys.getenv("VTI_DATA_ROOT", "./data")
FEATURE_DIR  <- Sys.getenv("VTI_FEATURE_DIR", "./output/feature_derivations")
SPLIT_DIR    <- file.path(DATA_ROOT, "h5ad/splits")
MANIFEST_CSV <- file.path(SPLIT_DIR, "_split_manifest.csv")
RADE_REPO    <- Sys.getenv("VTI_RADE_REPO", "")

FEATURE_VERSION <- "car_pos_fraction"
MIN_CELLS       <- 0L
TOL             <- 1e-9
DERIV_RUN_ID    <- format(Sys.time(), "carfrac_%Y%m%d_%H%M%S")

dir.create(FEATURE_DIR, showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# HELPERS
# ==============================================================================

sha256_file <- function(path) {
  if (!file.exists(path)) return(NA_character_)
  v <- system(paste0("sha256sum \"", path, "\" 2>/dev/null | cut -d' ' -f1"),
              intern = TRUE)
  if (length(v) == 0 || is.na(v) || v == "") NA_character_ else v
}

pkg_version <- function(pkg) as.character(packageVersion(pkg))

# Per-specimen read of the two obs datasets required for the fraction.
# Returns data.frame(numerator, denominator) for one h5ad file.
derive_one <- function(h5_path) {
  cbe <- h5read(h5_path, "obs/CAR_BY_EXPRS")
  ct3 <- h5read(h5_path, "obs/celltype_short_3")
  H5close()
  is_t      <- grepl("^CD4|^CD8", ct3)
  is_carpos <- (cbe == "TRUE")
  data.frame(
    numerator   = sum(is_t & is_carpos, na.rm = TRUE),
    denominator = sum(is_t, na.rm = TRUE)
  )
}

# ==============================================================================
# PROVENANCE: pin sources + environment
# ==============================================================================

manifest_sha <- sha256_file(MANIFEST_CSV)
r_version    <- R.version.string
dplyr_v      <- pkg_version("dplyr")
code_commit <- if (nzchar(RADE_REPO)) {
  tryCatch(
    system(paste0("git -C \"", RADE_REPO, "\" rev-parse HEAD 2>/dev/null"), intern = TRUE),
    error = function(e) NA_character_, warning = function(w) NA_character_
  )
} else NA_character_
if (length(code_commit) == 0) code_commit <- NA_character_

cat("=== Provenance ===\n")
cat("split_manifest sha256:", manifest_sha, "\n")
cat("R:", r_version, "\n")
cat("rhdf5:", pkg_version("rhdf5"), " | dplyr:", dplyr_v, "\n")
cat("extraction-code commit:", code_commit, "\n")
cat("feature version:", FEATURE_VERSION, " | min_cells:", MIN_CELLS, "\n")
cat("derivation_run_id:", DERIV_RUN_ID, "\n\n")

# ==============================================================================
# LOAD MANIFEST (specimen existence frame, 135 rows)
# ==============================================================================

manifest <- read.csv(MANIFEST_CSV, stringsAsFactors = FALSE)
manifest <- manifest[order(manifest$PATIENT_ID, manifest$TIMEPOINT), ]
stopifnot(nrow(manifest) == 135L)
cat("Specimens (manifest):", nrow(manifest), "\n")

# ==============================================================================
# PER-SPECIMEN DERIVATION (h5ad obs -> one row per specimen)
# ==============================================================================

rows <- vector("list", nrow(manifest))
for (i in seq_len(nrow(manifest))) {
  sp    <- manifest$orig.ident[i]
  fpath <- file.path(SPLIT_DIR, paste0(sp, ".h5ad"))
  if (!file.exists(fpath)) {
    cat("WARNING: h5ad missing for", sp, "-> skipped\n")
    rows[[i]] <- data.frame(orig.ident = sp, numerator = NA_integer_,
                            denominator = NA_integer_,
                            source_file_present = FALSE)
    next
  }
  d <- derive_one(fpath)
  rows[[i]] <- data.frame(orig.ident = sp, numerator = d$numerator,
                          denominator = d$denominator,
                          source_file_present = TRUE)
  if (i %% 25 == 0) cat("  derived", i, "of", nrow(manifest), "specimens\n")
}
deriv <- do.call(rbind, rows)
cat("Per-specimen derivation complete.\n")

# ==============================================================================
# MERGE: 135 manifest rows x derivation (left join)
# ==============================================================================
# QC status distinguishes four non-ok cases that were previously conflated:
#   missing_source_file    : the .h5ad file is absent (workflow failure)
#   missing_below_min_cells: T cells counted but < MIN_CELLS (biological)
#   zero_tcells            : T cells counted == 0 (biological)
#   invalid_out_of_range   : fraction outside [0,1] (should never happen)
# Only qc_status == "ok" rows keep a numeric value; below-min cells and the
# zero-T cases keep value = NA so a downstream consumer cannot mistake a
# workflow/low-cell outcome for a real fraction. A true biological zero
# (numerator == 0 with denominator >= MIN_CELLS) legitimately keeps value = 0.
full <- manifest %>%
  left_join(deriv, by = "orig.ident") %>%
  mutate(
    value       = ifelse(!is.na(denominator) & denominator >= MIN_CELLS,
                         numerator / denominator, NA_real_),
    qc_status   = case_when(
      !source_file_present               ~ "missing_source_file",
      is.na(denominator) | denominator == 0 ~ "zero_tcells",
      denominator < MIN_CELLS            ~ "missing_below_min_cells",
      value < 0 | value > 1              ~ "invalid_out_of_range",
      TRUE                               ~ "ok"
    ),
    exclusion_reason = case_when(
      qc_status == "missing_source_file"      ~ "source .h5ad file not found",
      qc_status == "zero_tcells"              ~ "0 eligible T cells in specimen",
      qc_status == "missing_below_min_cells"  ~
        paste0("denominator < MIN_CELLS (", MIN_CELLS, ")"),
      qc_status == "invalid_out_of_range"     ~ "fraction outside [0,1]",
      TRUE                                     ~ NA_character_
    )
  )

# Attach per-row provenance. source_object_sha = each specimen's own h5ad sha256
# (carried from the manifest); all_cell_object_sha = the split manifest sha256
# (existence frame); t_cell_object_sha = NA (no separate T-cell object is read).
full$feature_version        <- FEATURE_VERSION
full$source_object_sha      <- full$sha256
full$all_cell_object_sha    <- manifest_sha
full$t_cell_object_sha      <- NA_character_
full$extraction_code_commit <- code_commit
full$r_version              <- r_version
full$seurat_version         <- manifest$seurat_version   # from the split step
full$dplyr_version          <- dplyr_v
full$min_cells              <- MIN_CELLS
full$derivation_run_id      <- DERIV_RUN_ID

# ==============================================================================
# VALIDATION CHECKS
# ==============================================================================

cat("\n=== Validation checks ===\n")
checks <- list()

bad_nd <- full %>% filter(!is.na(value), numerator > denominator)
checks$numerator_le_denominator <- nrow(bad_nd) == 0
cat("1. numerator <= denominator:", checks$numerator_le_denominator,
    "(violations:", nrow(bad_nd), ")\n")

val_check <- full %>%
  filter(!is.na(value), denominator > 0) %>%
  mutate(expected = numerator / denominator, diff = abs(value - expected))
checks$value_equals_num_over_den <- all(val_check$diff <= TOL)
cat("2. value == num/den (tol", TOL, "):", checks$value_equals_num_over_den,
    "(max diff:", max(val_check$diff, 0), ")\n")

zero_den <- full %>% filter(denominator == 0)
checks$no_emit_for_zero_den <- all(is.na(zero_den$value))
cat("3. no value for zero denominator:", checks$no_emit_for_zero_den,
    "(zero-den rows:", nrow(zero_den), ")\n")

# Stronger invariant: value is present iff qc_status == "ok". Catches the old
# failure mode where NA-coerced-to-0 produced a value 0/0=0 for a missing
# specimen that was then quietly emitted to MEASUREMENT.
checks$value_present_iff_ok <-
  all((!is.na(full$value)) == (full$qc_status == "ok"))
cat("4. value present iff qc_status==ok:", checks$value_present_iff_ok, "\n")

range_check <- full %>% filter(!is.na(value))
checks$range_0_1 <- all(range_check$value >= 0 & range_check$value <= 1)
cat("5. fractions in [0,1]:", checks$range_0_1, "\n")

dup <- full %>% group_by(orig.ident, feature_version) %>% tally() %>% filter(n > 1)
checks$no_duplicate_specimen_feature <- nrow(dup) == 0
cat("6. no duplicate specimen-feature-version:", checks$no_duplicate_specimen_feature, "\n")

bij1 <- full %>% group_by(orig.ident) %>% summarise(np = n_distinct(PATIENT_ID)) %>% filter(np > 1)
bij2 <- full %>% group_by(orig.ident) %>% summarise(nt = n_distinct(TIMEPOINT)) %>% filter(nt > 1)
checks$identifiers_map_uniquely <- nrow(bij1) == 0 && nrow(bij2) == 0
cat("7. identifiers map uniquely:", checks$identifiers_map_uniquely, "\n")

all_pass <- all(unlist(checks))
cat("\nALL CHECKS PASS:", all_pass, "\n")

# ==============================================================================
# WRITE OUTPUTS
# ==============================================================================

out_csv <- file.path(FEATURE_DIR, "car_fraction_derivation.csv")
write.csv(full, out_csv, row.names = FALSE)
cat("Wrote:", out_csv, "(", nrow(full), "rows)\n")

checks_df <- data.frame(check = names(checks), passed = unlist(checks), row.names = NULL)
write.csv(checks_df,
          file.path(FEATURE_DIR, "car_fraction_validation_checks.csv"),
          row.names = FALSE)
cat("Wrote: car_fraction_validation_checks.csv\n")

prov <- data.frame(
  key = c("feature_version", "min_cells", "tolerance", "derivation_run_id",
          "split_manifest", "split_manifest_sha256",
          "extraction_code_commit", "r_version", "rhdf5_version", "dplyr_version",
          "n_specimens", "n_with_value", "n_missing", "all_checks_pass"),
  value = c(FEATURE_VERSION, MIN_CELLS, TOL, DERIV_RUN_ID,
            MANIFEST_CSV, manifest_sha,
            code_commit, r_version, pkg_version("rhdf5"), dplyr_v,
            nrow(full), sum(!is.na(full$value)), sum(is.na(full$value)), all_pass)
)
write.csv(prov, file.path(FEATURE_DIR, "car_fraction_provenance.csv"), row.names = FALSE)
cat("Wrote: car_fraction_provenance.csv\n")

# ==============================================================================
# SUMMARY
# ==============================================================================

cat("\n=== Summary ===\n")
cat("Specimens:                       ", nrow(full), "\n")
cat("With value (ok):                 ", sum(full$qc_status == "ok"), "\n")
cat("Missing (source file):          ", sum(full$qc_status == "missing_source_file"), "\n")
cat("Missing (below min cells):      ", sum(full$qc_status == "missing_below_min_cells"), "\n")
cat("Missing (zero T cells):         ", sum(full$qc_status == "zero_tcells"), "\n")
cat("Invalid:                         ", sum(full$qc_status == "invalid_out_of_range"), "\n")

if (!all_pass) quit(status = 1)
