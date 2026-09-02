# name: CERTAINTY_etl_scrna.R
# ==============================================================================
# ETL: Rade et al. CAR-T scRNA-seq linkage to OMOP CDM v5.4 (scRNA domain)
# ==============================================================================
#
# Inputs:
#   - person_lookup:   $VTI_OUTPUT_DIR/_aux/person_lookup.csv  (from CERTAINTY_etl_clinical.R)
#   - split manifest:  $VTI_DATA_ROOT/h5ad/splits/_split_manifest.csv
#   - specimen days:   $VTI_DATA_ROOT/clinical-data/clinical_data_for_seurat.Rds
#
# Outputs ($VTI_OUTPUT_DIR/):
#   - specimen.csv       (135 rows)
#   - observation.csv    (135 rows)
#
# Note: scRNA-derived feature encoding in MEASUREMENT is handled by
#   CERTAINTY_etl_scrna_features.R (CAR+ T-cell fraction, concept 2000010003,
#   linked to SPECIMEN via meas_event_field_concept_id = 1147822).

library(readr)    # read_csv, write_csv
library(dplyr)    # mutate, transmute, case_when, filter, select, n, n_distinct
library(stringr)  # strrep (also base::strrep)
library(hms)      # time parsing
library(lubridate)  # date arithmetic (NA_Date_ is base, lubridate for parity)

# ==============================================================================
# PATHS + CONSTANTS  (env-var driven; defaults are repo-relative)
# ==============================================================================

OUTPUT_DIR <- Sys.getenv("VTI_OUTPUT_DIR", "./output")
DATA_ROOT  <- Sys.getenv("VTI_DATA_ROOT", "./data")
AUX_DIR    <- file.path(OUTPUT_DIR, "_aux")
dir.create(AUX_DIR, showWarnings = FALSE, recursive = TRUE)

ID_OFFSET <- 1000L  # must match CERTAINTY_etl_clinical.R so specimen/measurement
                    # ids share the 1001+ space with clinical rows

TYPE_CRF              <- 32809L    # Case Report Form (matches clinical ETL)
SPECIMEN_PBMC         <- 4047495L  # Peripheral blood specimen
ANATOMIC_SITE_PBMC    <- 4229415L  # Peripheral blood mononuclear cell
SPECIMEN_MATERIAL     <- "PBMC"    # material string used in source_value rendering
OBS_SCRNA_DATASET     <- 2000010001L  # Custom: Single-cell RNA-seq dataset reference
OBS_EVENT_FIELD_SPEC  <- 1147822L  # specimen.specimen_id (CDM 5.4 field concept)
SPLIT_DIR           <- file.path(DATA_ROOT, "h5ad/splits")
SCRNA_CONTAINER_URI <- paste0("file://", normalizePath(SPLIT_DIR, mustWork = FALSE))

# ==============================================================================
# LOAD person_lookup (bridge from CERTAINTY_etl_clinical.R)
# ==============================================================================

person_lookup_path <- file.path(AUX_DIR, "person_lookup.csv")
if (!file.exists(person_lookup_path)) {
  stop(
    "_aux/person_lookup.csv not found. Run CERTAINTY_etl_clinical.R first"
  )
}
person_lookup <- read_csv(person_lookup_path, show_col_types = FALSE)
cat(sprintf("Loaded person_lookup: %d patients\n", nrow(person_lookup)))

# ==============================================================================
# TABLE: SPECIMEN + OBSERVATION (scRNA-seq linkage)
# ==============================================================================
# Sources:
#   - split manifest: one row per external scRNA-seq object (.h5ad file)
#   - specimen days: per-timepoint day offsets from the publication's frame
#   - person_lookup: maps PATIENT_ID -> person_id + anchored dates
#
# Design:
#   - SPECIMEN: one row per physical PBMC blood draw (patient x timepoint)
#   - OBSERVATION: one row per scRNA-seq dataset reference, 1:1 with SPECIMEN
#   - observation_event_id -> specimen_id (FK)
#   - obs_event_field_concept_id = 1147822 (specimen.specimen_id, CDM 5.4)
#
# Two columns serialize orig.ident for distinct roles (intentional, not a bug):
#   - SPECIMEN.specimen_source_id   = orig.ident  (source-system specimen id;
#     the lab sample tube id, also the .h5ad filename stem below)
#   - OBSERVATION.value_as_string    = orig.ident  (digital-asset key for the
#     scRNA dataset reference; OBSERVATION is about the *dataset*, SPECIMEN is
#     about the *biospecimen*; the dataset happens to be named after the biospecimen
#     it was measured on). The 1:1 join observation_event_id -> specimen_id is the
#     canonical link; the string match between the source ids is by construction,
#     NOT a substitute for the FK.
#   - SPECIMEN.specimen_source_value = "PBMC|<TIMEPOINT>"  (human-readable,
#     material + study-defined timepoint; WHAT was sampled)

# --- Load split manifest (ground truth: one row per external scRNA-seq object) ---
# The split manifest is content-addressable: one row per .h5ad file with
# orig.ident (durable object id), (PATIENT_ID, TIMEPOINT), n_cells, file, sha256.
# It replaces seu_t@meta.data as the SPECIMEN driver: 135 rows vs the 130 that
# the old distinct(PATIENT_ID, TIMEPOINT, TIME_*_CAR) produced ??? the 5
# "full-only" specimens that lack timing columns are now retained.
SPLIT_MANIFEST_PATH <- file.path(DATA_ROOT, "h5ad/splits/_split_manifest.csv")
SPLIT_MANIFEST <- read_csv(SPLIT_MANIFEST_PATH, show_col_types = FALSE)
cat(sprintf("Loaded split manifest: %d specimens (%d patients)\n",
            nrow(SPLIT_MANIFEST), n_distinct(SPLIT_MANIFEST$PATIENT_ID)))

# Per-specimen day offsets from the publication's per-timepoint frame. Keyed
# by orig.ident (unique in both frames)
SPECIMEN_DAYS_PATH <- file.path(DATA_ROOT, "clinical-data/clinical_data_for_seurat.Rds")
specimen_days <- readRDS(SPECIMEN_DAYS_PATH)
if (is.list(specimen_days) && !is.data.frame(specimen_days)) {
  specimen_days <- specimen_days[[1]]
}
specimen_days <- specimen_days %>%
  select(orig.ident, TIME_CAR_DAY_30, TIME_CAR_DAY_100)
cat(sprintf("Loaded specimen day offsets: %d rows keyed by orig.ident\n",
            nrow(specimen_days)))

# --- Join manifest to person_lookup + specimen_days and
# compute specimen_id + specimen_date once ---
# specimen_date: per-patient offsets from pdata.seurat joined by orig.ident.
# DDL allows NULL specimen_date ??? NA passthrough only, do not fabricate.
specimen_src <- SPLIT_MANIFEST %>%
  left_join(person_lookup, by = "PATIENT_ID") %>%
  left_join(specimen_days, by = "orig.ident")

# Hard-fail on missing person_id mapping instead of silently dropping rows.
# A dropped specimen means the manifest claims a biospecimen for a patient the
# clinical ETL never registered: DQD's fkClass would flag every downstream
# MEASUREMENT/OBSERVATION on that specimen. Fail loudly here so the operator
# fixes the person_lookup bridge instead of producing a silently truncated OMOP.
stopifnot("scRNA specimens missing person_id mapping: fix person_lookup bridge" =
            !anyNA(specimen_src$person_id))

specimen_src <- specimen_src %>%
  mutate(
    specimen_id = row_number() + ID_OFFSET,
    specimen_date = case_when(
      TIMEPOINT == "LP"         ~ apheresis_date,
      TIMEPOINT == "Late"       ~ car_infusion_date + TIME_CAR_DAY_30,
      TIMEPOINT == "Very Late"  ~ car_infusion_date + TIME_CAR_DAY_100,
      TRUE                      ~ NA_Date_
    )
  )

# specimen_id determinism: row_number() is taken over the manifest's row
# order. The split manifest is sorted PATIENT_ID, TIMEPOINT by the writer
# (CERTAINTY_split_all_cell_by_orig_ident.R), so the same inputs reproduce
# the same ids. If you reorder the manifest, specimen_id values will shift
# ORPHAN every prior OMOP load's *_event_id FKs. Reorder = re-run full ETL.
stopifnot("Late specimen missing TIME_CAR_DAY_30" =
            all(!is.na(specimen_src$TIME_CAR_DAY_30[specimen_src$TIMEPOINT == "Late"])))
stopifnot("Very Late specimen missing TIME_CAR_DAY_100" =
            all(!is.na(specimen_src$TIME_CAR_DAY_100[specimen_src$TIMEPOINT == "Very Late"])))
cat(sprintf("  specimens after person_lookup join: %d (unmatched: %d, all hard-failed above)\n",
            nrow(specimen_src), nrow(SPLIT_MANIFEST) - nrow(specimen_src)))

# --- Build SPECIMEN table (135 rows, driven by the manifest) ---
specimen <- specimen_src %>%
  transmute(
    specimen_id               = specimen_id,
    person_id                 = person_id,
    specimen_concept_id       = SPECIMEN_PBMC,
    specimen_type_concept_id  = TYPE_CRF,
    specimen_date             = specimen_date,
    anatomic_site_concept_id  = ANATOMIC_SITE_PBMC,
    specimen_source_value     = paste0(SPECIMEN_MATERIAL, "|", TIMEPOINT),
    specimen_source_id        = orig.ident      # lab sample tube id (= .h5ad stem)
  )

# Hard-fail rather than silently emit NULL specimen_date rows (the DDL allows
# NULL, but the features ETL joins on specimen_date for measurement_date, and a
# NULL there would produce a MEASUREMENT with no observation day). Only the
# "Full" timepoint has no day offset; if any non-Full specimen loses its date
# here, the join upstream is broken.
nonfull_missing_date <- specimen %>%
  filter(is.na(specimen_date) & specimen_source_value != "PBMC|Full")
stopifnot("non-Full specimen missing specimen_date" = nrow(nonfull_missing_date) == 0)

write_csv(specimen, file.path(OUTPUT_DIR, "specimen.csv"))
cat(sprintf("SPECIMEN: %d rows\n", nrow(specimen)))
cat(sprintf("  PBMC|LP=%d, PBMC|Late=%d, PBMC|Very Late=%d, PBMC|Full=%d\n",
            sum(specimen$specimen_source_value == "PBMC|LP"),
            sum(specimen$specimen_source_value == "PBMC|Late"),
            sum(specimen$specimen_source_value == "PBMC|Very Late"),
            sum(specimen$specimen_source_value == "PBMC|Full")))
cat(sprintf("  specimen_date NA (full-only): %d\n", sum(is.na(specimen$specimen_date))))

# --- Build OBSERVATION table (scRNA-seq linkage, 1:1 with SPECIMEN) ---
# value_as_string stores a file URI using the path visible inside the
# scrna-features container. SPECIMEN.specimen_source_id retains the original
# orig.ident (the durable .h5ad filename stem). This is the *digital-asset* role
# of the URI; the
# SPECIMEN.specimen_source_id above is the *biospecimen* role. Same string,
# different FK semantics (see the design block). The SHA-256 stays in the
# governed source manifest (_split_manifest.csv), not in this OMOP column.
# The modified PostgreSQL DDL declares value_as_string as TEXT so the
# container-visible URI is not constrained by the standard varchar(60) limit.
observation_scrna <- specimen_src %>%
  transmute(
    observation_id                = row_number() + ID_OFFSET,
    person_id                     = person_id,
    observation_concept_id        = OBS_SCRNA_DATASET,
    observation_date              = specimen_date,
    observation_type_concept_id   = TYPE_CRF,
    value_as_string               = paste0(SCRNA_CONTAINER_URI, "/", orig.ident, ".h5ad"),
    observation_event_id          = specimen_id,
    obs_event_field_concept_id    = OBS_EVENT_FIELD_SPEC
  )

write_csv(observation_scrna, file.path(OUTPUT_DIR, "observation.csv"))
cat(sprintf("OBSERVATION (scRNA-seq): %d rows\n", nrow(observation_scrna)))

# --- Verification ---
cat(sprintf("\n  Linkage check:\n"))
cat(sprintf("    specimen_id range: %d to %d\n",
            min(specimen$specimen_id), max(specimen$specimen_id)))
cat(sprintf("    observation_event_id range: %d to %d\n",
            min(observation_scrna$observation_event_id), max(observation_scrna$observation_event_id)))
cat(sprintf("    All observation_event_id in specimen_id: %s\n",
            all(observation_scrna$observation_event_id %in% specimen$specimen_id)))
cat(sprintf("    1:1 specimen <-> observation: %s\n",
            nrow(specimen) == nrow(observation_scrna)))

# --- Show example value_as_string ---
cat("\n  Example value_as_string (first 3):\n")
cat(paste0("    ", head(observation_scrna$value_as_string, 3), collapse = "\n"), "\n")

# ==============================================================================
# SUMMARY
# ==============================================================================

cat("\n")
cat(strrep("=", 60), "\n")
cat("scRNA ETL COMPLETE")
cat("Output directory:", OUTPUT_DIR, "\n")
cat(strrep("=", 60), "\n")
cat(sprintf("  specimen.csv:               %4d rows\n", nrow(specimen)))
cat(sprintf("  observation.csv:            %4d rows\n", nrow(observation_scrna)))
cat(sprintf("  ----------------------------------------\n"))
cat(sprintf("  TOTAL:                      %4d rows\n",
            nrow(specimen) + nrow(observation_scrna)))
cat(strrep("=", 60), "\n")

