# name: CERTAINTY_etl_scrna_features.R
# ==============================================================================
# ETL: scRNA-derived CAR-positive T-cell fraction -> OMOP MEASUREMENT (CDM 5.4)
# ==============================================================================
#
# Lifts the per-specimen car_pos_fraction derived outside OMOP into the
# MEASUREMENT table. One MEASUREMENT row per specimen with qc_status == "ok".
# Specimens with a missing/invalid fraction get NO row.
#
# Schema: 13 columns (CDM 5.4 MEASUREMENT minimal set + specimen link):
#   measurement_id, person_id, measurement_concept_id, measurement_date,
#   measurement_type_concept_id, value_as_number, value_as_concept_id,
#   unit_concept_id, measurement_source_value,
#   measurement_event_id, meas_event_field_concept_id,
#   range_low, range_high
# The specimen link (measurement_event_id = specimen_id,
# meas_event_field_concept_id = 1147822 = specimen.specimen_id) lets a consumer
# join MEASUREMENT back to SPECIMEN through the standard *_event_id /
# *_event_field_concept_id pair instead of parsing measurement_source_value.
# Clinical rows (from CERTAINTY_etl_clinical.R) are 9-column; bind_rows grows
# them with NA for the 4 event/range columns, matching the CDM 5.4 DDL.
#
# Inputs:
#   - derivation:  $VTI_FEATURE_DIR/car_fraction_derivation.csv  (from CERTAINTY_derive_car_fraction.R)
#   - specimen:    $VTI_OUTPUT_DIR/specimen.csv     (from CERTAINTY_etl_scrna.R)
#   - measurement: $VTI_OUTPUT_DIR/measurement.csv (from CERTAINTY_etl_clinical.R)
#
# Outputs ($VTI_OUTPUT_DIR/):
#   - measurement.csv   OVERWRITTEN with clinical rows + new feature rows (13-col)
#   - _aux/measurement_scrna_feature_manifest.csv  per-row provenance map
#
# ID range: feature rows start at 5001, disjoint from clinical (1001..1364).
#
# Run inside the RStudio container (bind your data root to the path
# configured by VTI_DATA_ROOT):
#   apptainer exec --bind $PWD/data:/data singularity-rstudio-4-3-2.sif \
#     Rscript scripts/CERTAINTY_etl_scrna_features.R \
#     [VTI_DATA_ROOT=/data VTI_OUTPUT_DIR=/data/output]

suppressMessages({
  library(readr)    # read_csv, write_csv
  library(dplyr)    # mutate, transmute, filter, n, n_distinct, case_when
})

# ==============================================================================
# PATHS + CONSTANTS  (env-var driven; defaults are repo-relative)
# ==============================================================================

DATA_ROOT   <- Sys.getenv("VTI_DATA_ROOT", "./data")
OUTPUT_DIR  <- Sys.getenv("VTI_OUTPUT_DIR", "./output")
FEATURE_DIR <- Sys.getenv("VTI_FEATURE_DIR", "./output/feature_derivations")
AUX_DIR     <- file.path(OUTPUT_DIR, "_aux")

DERIV_CSV    <- file.path(FEATURE_DIR, "car_fraction_derivation.csv")
SPECIMEN_CSV <- file.path(OUTPUT_DIR,   "specimen.csv")
MEAS_CSV     <- file.path(OUTPUT_DIR,   "measurement.csv")
MANIFEST_CSV <- file.path(AUX_DIR, "measurement_scrna_feature_manifest.csv")

MEAS_OFFSET         <- 5000L                       # disjoint from clinical 1001..1364
TYPE_CRF            <- 32809L                       # Case Report Form (matches clinical ETL)
FEATURE_CONCEPT_ID  <- 2000010003L                  # CAR-positive T-cell fraction by scRNA-seq
MEAS_EVENT_FIELD    <- 1147822L                     # specimen.specimen_id (CDM 5.4 field concept)
FEATURE_VERSION     <- "car_pos_fraction"

# ==============================================================================
# LOAD INPUTS
# ==============================================================================

dir.create(AUX_DIR, showWarnings = FALSE, recursive = TRUE)

stopifnot("run CERTAINTY_derive_car_fraction.R first (car_fraction_derivation.csv missing)" =
            file.exists(DERIV_CSV))
stopifnot("run CERTAINTY_etl_scrna.R first (specimen.csv missing)" =
            file.exists(SPECIMEN_CSV))
stopifnot("run CERTAINTY_etl_clinical.R first (measurement.csv missing)" =
            file.exists(MEAS_CSV))

deriv    <- read_csv(DERIV_CSV,    show_col_types = FALSE)
specimen <- read_csv(SPECIMEN_CSV, show_col_types = FALSE)
existing <- read_csv(MEAS_CSV,    show_col_types = FALSE) %>%
  filter(measurement_concept_id != FEATURE_CONCEPT_ID)  # drop prior feature rows

cat(sprintf("Loaded derivation: %d rows (%d ok, %d missing/invalid)\n",
            nrow(deriv),
            sum(deriv$qc_status == "ok"),
            sum(deriv$qc_status != "ok")))
cat(sprintf("Loaded specimen:   %d rows\n", nrow(specimen)))
cat(sprintf("Loaded measurement: %d clinical rows\n", nrow(existing)))

# ==============================================================================
# JOIN derivation -> specimen (on orig.ident = specimen_source_id)
# ==============================================================================
# Left join from specimen guarantees one row per specimen.

joined <- specimen %>%
  left_join(deriv, by = join_by(specimen_source_id == orig.ident))

lost <- sum(is.na(joined$qc_status))
if (lost > 0) {
  cat(sprintf("  WARNING: %d specimens with no derivation row (dropped)\n", lost))
  joined <- joined %>% filter(!is.na(qc_status))
}

cat(sprintf("  joined: %d specimens (one per specimen.csv row)\n", nrow(joined)))

# ==============================================================================
# BUILD MEASUREMENT ROWS (qc_status == "ok" only)
# ==============================================================================

ok_rows <- joined %>% filter(qc_status == "ok")

if (nrow(ok_rows) == 0) {
  stop("No specimens with qc_status == 'ok'. Aborting.")
}

ok_rows <- ok_rows %>% arrange(specimen_id) %>%
  mutate(measurement_id = row_number() + MEAS_OFFSET)

meas_new <- ok_rows %>%
  transmute(
    measurement_id              = measurement_id,
    person_id                   = person_id,
    measurement_concept_id      = FEATURE_CONCEPT_ID,
    measurement_date            = specimen_date,
    measurement_type_concept_id = TYPE_CRF,
    value_as_number             = value,
    value_as_concept_id         = 0L,
    unit_concept_id             = NA_integer_,           # dimensionless, no unit concept
    measurement_source_value    = paste0(FEATURE_VERSION, "|", specimen_source_id),
    measurement_event_id        = specimen_id,           # FK to SPECIMEN
    meas_event_field_concept_id = MEAS_EVENT_FIELD,      # specimen.specimen_id
    range_low                   = 0,
    range_high                  = 1
  )

cat(sprintf("\nNew MEASUREMENT rows: %d (qc_status == 'ok')\n", nrow(meas_new)))
cat(sprintf("  id range: %d to %d\n",
            min(meas_new$measurement_id), max(meas_new$measurement_id)))

# ==============================================================================
# VALIDATION CHECKS
# ==============================================================================

cat("\n=== Validation checks ===\n")
checks <- list()

checks$id_disjoint_from_clinical <-
  !any(meas_new$measurement_id %in% existing$measurement_id)
cat("1. new measurement_id disjoint from clinical:",
    checks$id_disjoint_from_clinical, "\n")

combined <- bind_rows(existing, meas_new)
checks$ids_unique_globally <- nrow(combined) == n_distinct(combined$measurement_id)
cat("2. all measurement_id unique (clinical + new):",
    checks$ids_unique_globally, "\n")

checks$all_join_to_specimen <-
  all(ok_rows$specimen_id %in% specimen$specimen_id)
cat("3. every feature row has a linked SPECIMEN:",
    checks$all_join_to_specimen, "\n")

spec_sub <- specimen %>% select(specimen_id, person_id, specimen_date)
date_check <- ok_rows %>%
  select(specimen_id, person_id, measurement_date = specimen_date) %>%
  left_join(spec_sub, by = "specimen_id", suffix = c("_meas", "_spec"))
checks$person_date_agree <-
  all(date_check$person_id_meas == date_check$person_id_spec, na.rm = TRUE) &&
  all(date_check$measurement_date == date_check$specimen_date, na.rm = TRUE)
cat("4. person_id + date agree with linked SPECIMEN:",
    checks$person_date_agree, "\n")

checks$value_in_0_1 <- all(meas_new$value_as_number >= 0 &
                           meas_new$value_as_number <= 1, na.rm = TRUE)
cat("5. value_as_number in [0,1]:", checks$value_in_0_1, "\n")

checks$no_missing_as_zero <- !any(joined$value == 0 & joined$qc_status != "ok"
                                  & !is.na(joined$value))
cat("6. no missing feature encoded as zero:",
    checks$no_missing_as_zero, "\n")

checks$source_value_reconstructs_specimen <-
  n_distinct(meas_new$measurement_source_value) == nrow(meas_new)
cat("7. measurement_source_value unique per row:",
    checks$source_value_reconstructs_specimen, "\n")

checks$meas_event_id_in_specimen <-
  all(meas_new$measurement_event_id %in% specimen$specimen_id)
cat("8. measurement_event_id -> SPECIMEN.specimen_id:",
    checks$meas_event_id_in_specimen, "\n")

checks$meas_event_field_correct <-
  all(meas_new$meas_event_field_concept_id == MEAS_EVENT_FIELD)
cat("9. meas_event_field_concept_id == 1147822:",
    checks$meas_event_field_correct, "\n")

checks$range_low_high_0_1 <-
  all(meas_new$range_low == 0 & meas_new$range_high == 1, na.rm = TRUE)
cat("10. range_low=0, range_high=1:", checks$range_low_high_0_1, "\n")

all_pass <- all(unlist(checks))
cat("\nALL CHECKS PASS:", all_pass, "\n")

# ==============================================================================
# WRITE COMBINED MEASUREMENT.CSV + PROVENANCE MANIFEST
# ==============================================================================

write_csv(combined, MEAS_CSV)
cat(sprintf("\nWrote: measurement.csv (%d clinical + %d feature = %d total rows)\n",
            nrow(existing), nrow(meas_new), nrow(combined)))

# Manifest: one row per feature MEASUREMENT, with the specimen link +
# provenance carried through from the derivation (sha256, git commit, run id,
# num/den). Stays in _aux/ so the OMOP MEASUREMENT row keeps the standard 13-col
# link while the heavier provenance lives in the sidecar.
manifest <- ok_rows %>%
  transmute(
    measurement_id              = measurement_id,
    specimen_id                 = specimen_id,
    person_id                   = person_id,
    PATIENT_ID                  = PATIENT_ID,
    TIMEPOINT                   = TIMEPOINT,
    orig.ident                  = specimen_source_id,
    feature_concept_id          = FEATURE_CONCEPT_ID,
    feature_version             = FEATURE_VERSION,
    measurement_date            = specimen_date,
    value_as_number             = value,
    numerator                   = numerator,
    denominator                 = denominator,
    qc_status                   = qc_status,
    min_cells                   = min_cells,
    source_object_sha           = source_object_sha,
    all_cell_object_sha         = all_cell_object_sha,
    t_cell_object_sha          = t_cell_object_sha,
    extraction_code_commit      = extraction_code_commit,
    r_version                   = r_version,
    seurat_version              = seurat_version,
    dplyr_version               = dplyr_version,
    derivation_run_id           = derivation_run_id,
    measurement_event_id        = specimen_id,
    meas_event_field_concept_id = MEAS_EVENT_FIELD
  )

write_csv(manifest, MANIFEST_CSV)
cat(sprintf("Wrote: _aux/measurement_scrna_feature_manifest.csv (%d rows)\n",
            nrow(manifest)))

# ==============================================================================
# SUMMARY
# ==============================================================================

cat("\n")
cat(strrep("=", 60), "\n")
cat("scRNA-feature ETL COMPLETE — Output directory:", OUTPUT_DIR, "\n")
cat(strrep("=", 60), "\n")
cat(sprintf("  measurement.csv (combined):   %4d rows\n", nrow(combined)))
cat(sprintf("    clinical (unchanged):         %4d rows\n", nrow(existing)))
cat(sprintf("    feature (concept %d):       %4d rows\n",
            FEATURE_CONCEPT_ID, nrow(meas_new)))
cat(sprintf("  ----------------------------------------\n"))
cat(sprintf("  specimens with value:          %4d / %d\n",
            nrow(meas_new), nrow(joined)))
cat(sprintf("  specimens skipped (no row):   %4d   (qc_status != 'ok')\n",
            nrow(joined) - nrow(meas_new)))
cat(sprintf("  ----------------------------------------\n"))
cat(sprintf("  provenance manifest:           %4d rows (_aux/)\n", nrow(manifest)))
cat(strrep("=", 60), "\n")

skipped <- joined %>% filter(qc_status != "ok")
if (nrow(skipped) > 0) {
  cat("\nSkipped specimens by reason:\n")
  for (r in names(table(skipped$qc_status))) {
    cat(sprintf("  %-40s %d\n", r, sum(skipped$qc_status == r)))
  }
}

if (!all_pass) quit(status = 1)
