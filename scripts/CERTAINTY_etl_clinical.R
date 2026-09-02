# name: CERTAINTY_etl_clinical.R
# ==============================================================================
# ETL: Rade et al. CAR-T Clinical Data to OMOP CDM v5.4 (clinical domain only)
# ==============================================================================
# Inputs:
#   - Source data:  $VTI_DATA_ROOT/clinical-data/clinical_table_DF_2024_10_28.Rds
#   - Dictionary:   $VTI_VOCAB_DIR/source_to_concept_map.csv  (SOURCE_TO_CONCEPT_MAP)
#
# Outputs ($VTI_OUTPUT_DIR/):
#   - person.csv, observation_period.csv, condition_occurrence.csv,
#     procedure_occurrence.csv, drug_exposure.csv, death.csv,
#     measurement.csv
#   - _aux/person_lookup.csv (PATIENT_ID -> person_id + anchored dates

library(readr)    # read_csv, write_csv
library(dplyr)    # mutate, transmute, case_when, filter, select, n, n_distinct
library(tidyr)    # pivot_longer
library(purrr)    # list helpers
library(stringr)  # strrep (also base::strrep)
library(hms)      # time parsing
library(lubridate)  # date arithmetic (NA_Date_ is base, lubridate for parity)

# ==============================================================================
# PATHS  (env-var driven; defaults are repo-relative)
# ==============================================================================

VOCAB_DIR  <- Sys.getenv("VTI_VOCAB_DIR", "./vocabulary")
DATA_ROOT  <- Sys.getenv("VTI_DATA_ROOT", "./data")
OUTPUT_DIR <- Sys.getenv("VTI_OUTPUT_DIR", "./output")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
file.path(OUTPUT_DIR, "_aux") |> dir.create(showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# LOAD SOURCE-TO-CONCEPT MAP
# ==============================================================================

s2cm <- read_csv(file.path(VOCAB_DIR, "source_to_concept_map.csv"), show_col_types = FALSE)
cat(sprintf("Loaded source_to_concept_map.csv: %d mappings across %d source vocabularies\n",
            nrow(s2cm),
            n_distinct(s2cm$source_vocabulary_id)))

# Helper: look up target_concept_id by source_code (optionally filtered by vocabulary)
lookup_concept <- function(code, vocab = NULL) {
  matches <- s2cm %>% filter(source_code == code)
  if (!is.null(vocab)) {
    matches <- matches %>% filter(source_vocabulary_id == vocab)
  }
  if (nrow(matches) == 0) {
    warning(paste("No s2cm match for:", code))
    return(0L)
  }
  if (nrow(matches) > 1) {
    warning(paste("Multiple s2cm matches for:", code, "- using first"))
  }
  as.integer(matches$target_concept_id[1])
}

# ==============================================================================
# LOAD SOURCE DATA
# ==============================================================================

clinical_data <- readRDS(
  file.path(DATA_ROOT, "clinical-data/clinical_table_DF_2024_10_28.Rds")
)

pdata_clin  <- clinical_data$pdata.clin
pdata_crs   <- clinical_data$pdata.crs
pdata_imm_s <- clinical_data$pdata.imm.s
pdata_t_s   <- clinical_data$pdata.t.s
pdata_elisa <- clinical_data$pdata.elisa

cat(sprintf("  pdata_clin:  %d patients\n", nrow(pdata_clin)))
cat(sprintf("  pdata_crs:   %d rows\n", nrow(pdata_crs)))
cat(sprintf("  pdata_imm_s: %d rows\n", nrow(pdata_imm_s)))
cat(sprintf("  pdata_t_s:   %d rows\n", nrow(pdata_t_s)))
cat(sprintf("  pdata_elisa: %d rows\n", nrow(pdata_elisa)))

# ==============================================================================
# RESOLVE CONCEPTS INTO NAMED CONSTANTS
# ==============================================================================

cat("Resolving concepts from s2cm...\n")

# --- Infrastructure (not source-data mappings → hardcoded per OMOP convention) ---
TYPE_CRF <- 32809L   # Case Report Form (§4.4 codebook)
TYPE_EHR <- 32817L   # EHR

# --- Gender ---
MALE   <- lookup_concept("m", "ULEI_CLIN")
FEMALE <- lookup_concept("w", "ULEI_CLIN")

# --- Conditions ---
CONDITION_MM    <- lookup_concept("Multiple myeloma", "ULEI_CLIN")
CONDITION_CRS   <- lookup_concept("CRS", "ULEI_CRS")
CONDITION_ICANS <- lookup_concept("ICANS", "ULEI_CRS")

# --- Procedures ---
PROC_LEUKAPHERESIS  <- lookup_concept("Leukapheresis", "ULEI_CLIN")
PROC_CAR_T_INFUSION <- lookup_concept("CAR-T infusion", "ULEI_CLIN")

# --- Drugs ---
DRUG_IDE_CEL   <- lookup_concept("ide", "ULEI_CLIN")
DRUG_CILTA_CEL <- lookup_concept("cilta", "ULEI_CLIN")
DRUG_TOCI      <- lookup_concept("TOCI", "ULEI_CRS")

# --- Measurements ---
MEAS_RISS         <- lookup_concept("R_ISS_measurement", "ULEI_STAGING")
MEAS_DURIE_SALMON <- lookup_concept("Durie_Salmon_measurement", "ULEI_STAGING")
MEAS_CRP          <- lookup_concept("CRP", "ULEI_CRS")

# --- Measurement values (R-ISS stages) ---
RISS_1 <- lookup_concept("R_ISS_1", "ULEI_STAGING")
RISS_2 <- lookup_concept("R_ISS_2", "ULEI_STAGING")
RISS_3 <- lookup_concept("R_ISS_3", "ULEI_STAGING")

# --- Treatment response values ---
RESPONSE_CR <- lookup_concept("CR", "ULEI_RESPONSE")
RESPONSE_PR <- lookup_concept("PR", "ULEI_RESPONSE")
RESPONSE_SD <- lookup_concept("SD", "ULEI_RESPONSE")
RESPONSE_PD <- lookup_concept("PD", "ULEI_RESPONSE")

# --- Units ---
UNIT_MG_L     <- lookup_concept("mg/L", "ULEI_UNIT")
UNIT_PG_ML    <- lookup_concept("pg/mL", "ULEI_UNIT")
UNIT_NG_ML    <- lookup_concept("ng/mL", "ULEI_UNIT")
UNIT_CELLS_UL <- lookup_concept("cells/uL", "ULEI_UNIT")
UNIT_PERCENT  <- lookup_concept("%", "ULEI_UNIT")

# --- Verification ---
cat(sprintf("
  Gender:     MALE=%d, FEMALE=%d
  Conditions: MM=%d, CRS=%d, ICANS=%d
  Procedures: LEUK=%d, CAR_T=%d
  Drugs:      IDE=%d, CILTA=%d, TOCI=%d
  Meas:       RISS=%d, CRP=%d, DURIE_SALMON=%d
  R-ISS vals: I=%d, II=%d, III=%d
  Responses:  CR=%d, PR=%d, SD=%d, PD=%d
  Units:      mg/L=%d, pg/mL=%d, ng/mL=%d, cells/uL=%d, %%=%d
",
            MALE, FEMALE,
            CONDITION_MM, CONDITION_CRS, CONDITION_ICANS,
            PROC_LEUKAPHERESIS, PROC_CAR_T_INFUSION,
            DRUG_IDE_CEL, DRUG_CILTA_CEL, DRUG_TOCI,
            MEAS_RISS, MEAS_CRP, MEAS_DURIE_SALMON,
            RISS_1, RISS_2, RISS_3,
            RESPONSE_CR, RESPONSE_PR, RESPONSE_SD, RESPONSE_PD,
            UNIT_MG_L, UNIT_PG_ML, UNIT_NG_ML, UNIT_CELLS_UL, UNIT_PERCENT
))

ID_OFFSET <- 1000

# ==============================================================================
# SYNTHETIC DATE ANCHORING
# ==============================================================================
# Assumptions:
#   a) Synthetic anchor date = 2023-01-01 (same for all patients)
#   b) No inter-patient temporal ordering is implied
#   c) All OMOP dates = anchor +/- offset (intra-patient intervals exact)
#   d) No output date is a true calendar date
#   e) year_of_birth = 2023 - AGE_AT_CAR

ANCHOR_DATE <- as.Date("2023-01-01")
DAYS_PER_MONTH <- 365.25 / 12

pdata_clin <- pdata_clin %>%
  arrange(PATIENT_ID) %>%
  mutate(
    person_id         = row_number() + ID_OFFSET,
    car_infusion_date = ANCHOR_DATE,
    apheresis_date    = ANCHOR_DATE - TIME_APHERESIS_CAR,
    diagnosis_date    = ANCHOR_DATE - TIME_DIAGNOSE_CAR * DAYS_PER_MONTH,
    lfu_date          = ANCHOR_DATE + LFU,
    death_date        = if_else(DEATH == 1, ANCHOR_DATE + OS, NA_Date_)
  )

# Person lookup for downstream joins
person_lookup <- pdata_clin %>%
  select(PATIENT_ID, SAMPLE_ID, person_id,
         car_infusion_date, apheresis_date, diagnosis_date, lfu_date)

cat(sprintf("Anchored %d patients to %s\n", nrow(person_lookup), ANCHOR_DATE))
cat(sprintf("  Apheresis range:  %s to %s\n",
            min(pdata_clin$apheresis_date, na.rm = TRUE),
            max(pdata_clin$apheresis_date, na.rm = TRUE)))
cat(sprintf("  Diagnosis range:  %s to %s\n",
            min(pdata_clin$diagnosis_date, na.rm = TRUE),
            max(pdata_clin$diagnosis_date, na.rm = TRUE)))
cat(sprintf("  LFU range:        %s to %s\n",
            min(pdata_clin$lfu_date, na.rm = TRUE),
            max(pdata_clin$lfu_date, na.rm = TRUE)))

# ==============================================================================
# TABLE: PERSON
# ==============================================================================

person <- pdata_clin %>%
  transmute(
    person_id            = person_id,
    gender_concept_id    = case_when(
      SEX == "m" ~ MALE,
      SEX == "w" ~ FEMALE,
      TRUE       ~ 0L
    ),
    year_of_birth        = as.integer(2023 - AGE_AT_CAR),
    month_of_birth       = NA_integer_,
    day_of_birth         = NA_integer_,
    race_concept_id      = 0L,
    ethnicity_concept_id = 0L,
    person_source_value  = PATIENT_ID,
    gender_source_value  = SEX
  )

write_csv(person, file.path(OUTPUT_DIR, "person.csv"))
cat(sprintf("PERSON: %d rows written\n", nrow(person)))
cat(sprintf("  Gender: M=%d, F=%d, Unknown=%d\n",
            sum(person$gender_concept_id == MALE),
            sum(person$gender_concept_id == FEMALE),
            sum(person$gender_concept_id == 0L)))
cat(sprintf("  year_of_birth range: %d to %d\n",
            min(person$year_of_birth, na.rm = TRUE),
            max(person$year_of_birth, na.rm = TRUE)))

# ==============================================================================
# TABLE: OBSERVATION_PERIOD
# ==============================================================================
# observation_period_start_date = diagnosis_date (earliest recorded event)
# observation_period_end_date   = lfu_date (last follow-up = censoring)

observation_period <- pdata_clin %>%
  transmute(
    observation_period_id         = person_id,
    person_id                     = person_id,
    observation_period_start_date = diagnosis_date,
    observation_period_end_date   = lfu_date,
    period_type_concept_id        = TYPE_CRF
  )

write_csv(observation_period, file.path(OUTPUT_DIR, "observation_period.csv"))
cat(sprintf("OBSERVATION_PERIOD: %d rows\n", nrow(observation_period)))
cat(sprintf("  Start range: %s to %s\n",
            min(observation_period$observation_period_start_date, na.rm = TRUE),
            max(observation_period$observation_period_start_date, na.rm = TRUE)))
cat(sprintf("  End range:   %s to %s\n",
            min(observation_period$observation_period_end_date, na.rm = TRUE),
            max(observation_period$observation_period_end_date, na.rm = TRUE)))

# ==============================================================================
# TABLE: CONDITION_OCCURRENCE
# ==============================================================================
# Sources:
#   - MM diagnosis: all patients (from pdata_clin)
#   - CRS: patients with CRS == "yes" (from pdata_crs)
#   - ICANS: patients with ICANS > 0 (from pdata_crs)

cond_id_counter <- ID_OFFSET  # Start at ID_OFFSET to avoid collision with demo data

# --- MM: all patients ---
cond_mm <- pdata_clin %>%
  mutate(cond_id = cond_id_counter + row_number()) %>%
  transmute(
    condition_occurrence_id   = cond_id,
    person_id                 = person_id,
    condition_concept_id      = CONDITION_MM,
    condition_start_date      = diagnosis_date,
    condition_end_date        = NA_Date_,
    condition_type_concept_id = TYPE_CRF,
    condition_source_value    = "Multiple myeloma"
  )
cond_id_counter <- cond_id_counter + nrow(cond_mm)

# --- CRS ---
crs_patients <- pdata_crs %>%
  filter(CRS == "yes") %>%
  left_join(person_lookup, by = "PATIENT_ID")

stopifnot("CRS patients missing person_id" = !anyNA(crs_patients$person_id))

cond_crs <- crs_patients %>%
  mutate(cond_id = cond_id_counter + row_number()) %>%
  transmute(
    condition_occurrence_id   = cond_id,
    person_id                 = person_id,
    condition_concept_id      = CONDITION_CRS,
    condition_start_date      = car_infusion_date,
    condition_end_date        = NA_Date_,
    condition_type_concept_id = TYPE_CRF,
    condition_source_value    = paste0("CRS Grade ", CRS_GRADE)
  )
cond_id_counter <- cond_id_counter + nrow(cond_crs)

# --- ICANS ---
icans_patients <- pdata_crs %>%
  filter(ICANS > 0) %>%
  left_join(person_lookup, by = "PATIENT_ID")

stopifnot("ICANS patients missing person_id" = !anyNA(icans_patients$person_id))

cond_icans <- icans_patients %>%
  mutate(cond_id = cond_id_counter + row_number()) %>%
  transmute(
    condition_occurrence_id   = cond_id,
    person_id                 = person_id,
    condition_concept_id      = CONDITION_ICANS,
    condition_start_date      = car_infusion_date,
    condition_end_date        = NA_Date_,
    condition_type_concept_id = TYPE_CRF,
    condition_source_value    = paste0("ICANS Grade ", ICANS)
  )
cond_id_counter <- cond_id_counter + nrow(cond_icans)

# --- Combine ---
condition_occurrence <- bind_rows(cond_mm, cond_crs, cond_icans)
write_csv(condition_occurrence, file.path(OUTPUT_DIR, "condition_occurrence.csv"))
cat(sprintf("CONDITION_OCCURRENCE: %d rows (MM=%d, CRS=%d, ICANS=%d)\n",
            nrow(condition_occurrence), nrow(cond_mm), nrow(cond_crs), nrow(cond_icans)))

# ==============================================================================
# TABLE: PROCEDURE_OCCURRENCE
# ==============================================================================
# Sources:
#   - Leukapheresis: all patients (date = apheresis_date)
#   - CAR-T infusion: all patients (date = car_infusion_date)

proc_id_counter <- ID_OFFSET  # Start at ID_OFFSET to avoid collision with demo data

# --- Leukapheresis ---
proc_aph <- pdata_clin %>%
  mutate(proc_id = proc_id_counter + row_number()) %>%
  transmute(
    procedure_occurrence_id   = proc_id,
    person_id                 = person_id,
    procedure_concept_id      = PROC_LEUKAPHERESIS,
    procedure_date            = apheresis_date,
    procedure_type_concept_id = TYPE_CRF,
    procedure_source_value    = "Leukapheresis"
  )
proc_id_counter <- proc_id_counter + nrow(proc_aph)

# --- CAR-T infusion ---
proc_inf <- pdata_clin %>%
  mutate(proc_id = proc_id_counter + row_number()) %>%
  transmute(
    procedure_occurrence_id   = proc_id,
    person_id                 = person_id,
    procedure_concept_id      = PROC_CAR_T_INFUSION,
    procedure_date            = car_infusion_date,
    procedure_type_concept_id = TYPE_CRF,
    procedure_source_value    = paste0("CAR-T infusion (", PRODUCT, ")")
  )
proc_id_counter <- proc_id_counter + nrow(proc_inf)

# --- Combine ---
procedure_occurrence <- bind_rows(proc_aph, proc_inf)
write_csv(procedure_occurrence, file.path(OUTPUT_DIR, "procedure_occurrence.csv"))
cat(sprintf("PROCEDURE_OCCURRENCE: %d rows (leukapheresis=%d, CAR-T=%d)\n",
            nrow(procedure_occurrence), nrow(proc_aph), nrow(proc_inf)))

# ==============================================================================
# TABLE: DRUG_EXPOSURE
# ==============================================================================
# Sources:
#   - CAR-T product (ide-cel / cilta-cel): from pdata_clin$PRODUCT
#   - Tocilizumab: from pdata_crs$TOCI (CRS management)

drug_id_counter <- ID_OFFSET  # Start at ID_OFFSET to avoid collision with demo data

# --- CAR-T products ---
drug_cart <- pdata_clin %>%
  mutate(drug_id = drug_id_counter + row_number()) %>%
  transmute(
    drug_exposure_id         = drug_id,
    person_id                = person_id,
    drug_concept_id          = case_when(
      PRODUCT == "ide"   ~ DRUG_IDE_CEL,
      PRODUCT == "cilta" ~ DRUG_CILTA_CEL,
      TRUE               ~ 0L
    ),
    drug_exposure_start_date = car_infusion_date,
    drug_exposure_end_date   = car_infusion_date,
    drug_type_concept_id     = TYPE_CRF,
    drug_source_value        = PRODUCT
  )
drug_id_counter <- drug_id_counter + nrow(drug_cart)

# --- Tocilizumab (for CRS management) ---
toci_patients <- pdata_crs %>%
  filter(TOCI == "yes") %>%
  left_join(person_lookup, by = "PATIENT_ID")

stopifnot("Tocilizumab patients missing person_id" = !anyNA(toci_patients$person_id))

drug_toci <- toci_patients %>%
  mutate(drug_id = drug_id_counter + row_number()) %>%
  transmute(
    drug_exposure_id         = drug_id,
    person_id                = person_id,
    drug_concept_id          = DRUG_TOCI,
    drug_exposure_start_date = car_infusion_date,
    drug_exposure_end_date   = car_infusion_date,
    drug_type_concept_id     = TYPE_CRF,
    drug_source_value        = paste0("Tocilizumab (", TOCI, " doses)")
  )
drug_id_counter <- drug_id_counter + nrow(drug_toci)

# --- Combine ---
drug_exposure <- bind_rows(drug_cart, drug_toci)
write_csv(drug_exposure, file.path(OUTPUT_DIR, "drug_exposure.csv"))
cat(sprintf("DRUG_EXPOSURE: %d rows (ide=%d, cilta=%d, toci=%d)\n",
            nrow(drug_exposure),
            sum(drug_cart$drug_concept_id == DRUG_IDE_CEL),
            sum(drug_cart$drug_concept_id == DRUG_CILTA_CEL),
            nrow(drug_toci)))

# ==============================================================================
# TABLE: DEATH
# ==============================================================================
# Source: pdata_clin$DEATH (1 = died), pdata_clin$OS (days to death/censoring)
# Per codebook §3.2.4: death_date used for OS in survival analysis

death <- pdata_clin %>%
  filter(DEATH == 1) %>%
  transmute(
    person_id             = person_id,
    death_date            = death_date,
    death_type_concept_id = TYPE_CRF,
    cause_concept_id      = 0L,
    cause_source_value    = NA_character_
  )

write_csv(death, file.path(OUTPUT_DIR, "death.csv"))
cat(sprintf("DEATH: %d rows (%.0f%% of cohort)\n",
            nrow(death),
            100 * nrow(death) / nrow(pdata_clin)))

# ==============================================================================
# TABLE: MEASUREMENT
# ==============================================================================
# Sources:
#   - R-ISS staging (categorical): pdata_clin$R_ISS
#   - CRP baseline (continuous): pdata_crs$CRP_BASE
#   - CRP max (continuous): pdata_crs$CRP_MAX
#   - Treatment responses at 1/3/6/12 months: pdata_clin$RESPONSE_*_M
#
# Design:
#   - R-ISS: measurement_concept_id = stage (Standard), value_as_concept_id = 0
#   - CRP: measurement_concept_id = 3020460, value_as_number, unit = mg/L
#   - Responses: per [ro3] the Cancer Modifier concept IS the measurement
#     (standalone, no value_as_concept_id needed)

meas_id_counter <- ID_OFFSET  # Start at ID_OFFSET to avoid collision with demo data
all_measurements <- list()

# --- R-ISS staging ---
riss_meas <- pdata_clin %>%
  filter(!is.na(R_ISS)) %>%
  mutate(meas_id = meas_id_counter + row_number()) %>%
  transmute(
    measurement_id              = meas_id,
    person_id                   = person_id,
    measurement_concept_id      = case_when(
      R_ISS == 1 ~ RISS_1,
      R_ISS == 2 ~ RISS_2,
      R_ISS == 3 ~ RISS_3,
      TRUE       ~ 0L
    ),
    measurement_date            = diagnosis_date,
    measurement_type_concept_id = TYPE_CRF,
    value_as_number             = NA_real_,
    value_as_concept_id         = 0L,
    unit_concept_id             = NA_integer_,  # categorical staging, dimensionless
    measurement_source_value    = paste0("R-ISS Stage ", R_ISS)
  )
meas_id_counter <- meas_id_counter + nrow(riss_meas)
all_measurements[["riss"]] <- riss_meas

# --- CRP baseline (pre-infusion, Day -1) ---
crp_base_data <- pdata_crs %>%
  filter(!is.na(CRP_BASE)) %>%
  left_join(person_lookup, by = "PATIENT_ID")

crp_base_meas <- crp_base_data %>%
  mutate(meas_id = meas_id_counter + row_number()) %>%
  transmute(
    measurement_id              = meas_id,
    person_id                   = person_id,
    measurement_concept_id      = MEAS_CRP,
    measurement_date            = car_infusion_date - 1L,
    measurement_type_concept_id = TYPE_CRF,
    value_as_number             = CRP_BASE,
    value_as_concept_id         = 0L,
    unit_concept_id             = UNIT_MG_L,
    measurement_source_value    = "CRP_BASE"
  )
meas_id_counter <- meas_id_counter + nrow(crp_base_meas)
all_measurements[["crp_base"]] <- crp_base_meas

# --- CRP max (peak post-infusion) ---
crp_max_data <- pdata_crs %>%
  filter(!is.na(CRP_MAX)) %>%
  left_join(person_lookup, by = "PATIENT_ID")

crp_max_meas <- crp_max_data %>%
  mutate(meas_id = meas_id_counter + row_number()) %>%
  transmute(
    measurement_id              = meas_id,
    person_id                   = person_id,
    measurement_concept_id      = MEAS_CRP,
    measurement_date            = car_infusion_date + 7L,
    measurement_type_concept_id = TYPE_CRF,
    value_as_number             = CRP_MAX,
    value_as_concept_id         = 0L,
    unit_concept_id             = UNIT_MG_L,
    measurement_source_value    = "CRP_MAX"
  )
meas_id_counter <- meas_id_counter + nrow(crp_max_meas)
all_measurements[["crp_max"]] <- crp_max_meas

# --- Treatment responses (1/3/6/12 months) ---
# Per [ro3]: the Cancer Modifier concept IS the measurement (standalone)
response_to_concept <- function(resp) {
  case_when(
    resp %in% c("CR", "sCR", "nCR") ~ RESPONSE_CR,
    resp %in% c("PR", "VGPR")       ~ RESPONSE_PR,
    resp == "SD"                      ~ RESPONSE_SD,
    resp == "PD"                      ~ RESPONSE_PD,
    TRUE                              ~ NA_integer_
  )
}

response_long <- pdata_clin %>%
  select(person_id, car_infusion_date,
         RESPONSE_1_M, RESPONSE_3_M, RESPONSE_6_M, RESPONSE_12_M) %>%
  pivot_longer(
    cols = starts_with("RESPONSE_"),
    names_to = "timepoint",
    values_to = "response"
  ) %>%
  filter(!is.na(response)) %>%
  mutate(
    days_offset = case_when(
      timepoint == "RESPONSE_1_M"  ~ 30L,
      timepoint == "RESPONSE_3_M"  ~ 90L,
      timepoint == "RESPONSE_6_M"  ~ 180L,
      timepoint == "RESPONSE_12_M" ~ 365L
    ),
    measurement_concept_id = response_to_concept(response)
  ) %>%
  filter(!is.na(measurement_concept_id))

response_meas <- response_long %>%
  mutate(meas_id = meas_id_counter + row_number()) %>%
  transmute(
    measurement_id              = meas_id,
    person_id                   = person_id,
    measurement_concept_id      = measurement_concept_id,
    measurement_date            = car_infusion_date + days_offset,
    measurement_type_concept_id = TYPE_CRF,
    value_as_number             = NA_real_,
    value_as_concept_id         = 0L,
    unit_concept_id             = NA_integer_,  # categorical response, dimensionless
    measurement_source_value    = paste0(timepoint, "=", response)
  )
meas_id_counter <- meas_id_counter + nrow(response_meas)
all_measurements[["responses"]] <- response_meas

# --- Write clinical MEASUREMENT (clinical rows only) ---
measurement <- bind_rows(all_measurements)
write_csv(measurement, file.path(OUTPUT_DIR, "measurement.csv"))
cat(sprintf("MEASUREMENT (clinical): %d rows written -> measurement.csv\n", nrow(measurement)))
cat(sprintf("  R-ISS:     %d\n", nrow(riss_meas)))
cat(sprintf("  CRP_BASE:  %d\n", nrow(crp_base_meas)))
cat(sprintf("  CRP_MAX:   %d\n", nrow(crp_max_meas)))
cat(sprintf("  Responses: %d\n", nrow(response_meas)))

AUX_DIR <- file.path(OUTPUT_DIR, "_aux")
dir.create(AUX_DIR, showWarnings = FALSE, recursive = TRUE)
write_csv(person_lookup, file.path(AUX_DIR, "person_lookup.csv"))
cat("  + _aux/person_lookup.csv (consumed by CERTAINTY_etl_scrna.R)\n")

# ==============================================================================

# ==============================================================================
# SUMMARY
# ==============================================================================

cat("\n")
cat(strrep("=", 60), "\n")
cat("CLINICAL ETL COMPLETE Output directory:", OUTPUT_DIR, "\n")
cat(strrep("=", 60), "\n")
cat(sprintf("  person.csv:                 %4d rows\n", nrow(person)))
cat(sprintf("  observation_period.csv:     %4d rows\n", nrow(observation_period)))
cat(sprintf("  condition_occurrence.csv:   %4d rows\n", nrow(condition_occurrence)))
cat(sprintf("  procedure_occurrence.csv:   %4d rows\n", nrow(procedure_occurrence)))
cat(sprintf("  drug_exposure.csv:          %4d rows\n", nrow(drug_exposure)))
cat(sprintf("  measurement.csv:            %4d rows (clinical only)\n", nrow(measurement)))
cat(sprintf("  death.csv:                  %4d rows\n", nrow(death)))
cat(sprintf("  ----------------------------------------\n"))
cat(sprintf("  TOTAL:                      %4d rows\n",
            nrow(person) + nrow(observation_period) + nrow(condition_occurrence) +
              nrow(procedure_occurrence) + nrow(drug_exposure) + nrow(measurement) +
              nrow(death)))
cat(strrep("=", 60), "\n")

