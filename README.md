# OMOP ETL — Rade et al. 2026 (CAR-T scRNA-seq)

OMOP CDM v5.4 ETL for the published single-cell CAR-T multiple myeloma cohort of
Rade, Fandrei, Kreuz et al. (2026). Produces 9 CDM v5.4 CSV tables from the
published research dataset.

## Related resources

- Paper: https://www.cell.com/cancer-cell/fulltext/S1535-6108(25)00490-8
- Source code (pre-processing, figures): https://github.com/fraunhofer-izi/Rade_et_al_CAR_2025
- Source data (Seurat object, Singularity image, R packages): https://zenodo.org/records/14732727  (restricted access, request via Zenodo)

## Inputs (not included; obtain from the sources above)

| File | Source | Used by |
|---|---|---|
| `05_seurat_harmony_all_new.Rds` | Zenodo 14732727 | `CERTAINTY_split_all_cell_by_orig_ident.R` |
| `data/clinical_data_for_seurat.Rds` | GitHub Rade_et_al_CAR_2025 `data/` | `CERTAINTY_etl_scrna.R` (per-specimen day offsets) |
| `publication/clinicial_data/clinical_table_DF_2024_10_28.Rds` | GitHub Rade_et_al_CAR_2025 `publication/clinicial_data/` | `CERTAINTY_etl_clinical.R` (cohort source, 61 patients) |

The Seurat object requires a restricted-access request to Zenodo (GDPR).
The two `.Rds` files are git-tracked on GitHub and freely downloadable.

Place the inputs under a single data root (default `./data`), keep the
included vocabulary under a vocabulary dir (default `./vocabulary`), and
point the ETL at an output dir (default `./output`). Override with env vars
(see **Paths** below).

## Vocabulary (included)

| File | Rows | Role |
|---|---|---|
| `vocabulary/source_to_concept_map.csv` | 30 | `SOURCE_TO_CONCEPT_MAP` (9-col, official CDM schema). Source codes (`m`, `ide`, `cilta`, `R_ISS_1..3`, `CRS`, `ICANS`, `TOCI`, `CRP`, `mg/L`, …) → standard OMOP `concept_id` across 8 `ULEI_*` source vocabularies. |
| `vocabulary/vocabulary.csv` | 8 | Project subset of `VOCABULARY` (2-col: `source_vocabulary_id,Covers`). Not the official 5-col `vocabulary` table schema. |
| `vocabulary/concept.csv` | 2 | `CONCEPT` (10-col, official CDM schema). Custom concepts in `CERTAINTY` vocabulary: `2000010001` (scRNA-seq dataset reference, Observation) and `2000010003` (CAR-positive T-cell fraction by scRNA-seq, Measurement). Must be inserted into `cdm.concept` for full DQD compliance. |

## Scripts (included)

Run order is data-dependency-driven, not filename order.

| # | Script | Inputs -> Outputs |
|---|---|---|
| 1 | `CERTAINTY_split_all_cell_by_orig_ident.R` | Seurat object -> `data/h5ad/splits/_split_manifest.csv` (135 rows) + `data/h5ad/splits/<orig.ident>.h5ad` (135 files) |
| 2 | `CERTAINTY_etl_clinical.R` | cohort RDS + `source_to_concept_map.csv` -> `person`, `observation_period`, `condition_occurrence`, `procedure_occurrence`, `drug_exposure`, `death`, `measurement` (clinical rows) + `_aux/person_lookup.csv` |
| 3 | `CERTAINTY_derive_car_fraction.R` | split manifest + 135 `.h5ad` -> `car_fraction_derivation.csv`, `car_fraction_provenance.csv`, `car_fraction_validation_checks.csv` |
| 4 | `CERTAINTY_etl_scrna.R` | split manifest + `person_lookup.csv` + specimen-days RDS -> `specimen.csv` (135), `observation.csv` (135) |
| 5 | `CERTAINTY_etl_scrna_features.R` | `measurement.csv` (from step 2) + `specimen.csv` (from step 4) + `car_fraction_derivation.csv` (from step 3) -> final `measurement.csv` (clinical + 130 CAR+ fraction rows, `concept_id=2000010003`) + `_aux/measurement_scrna_feature_manifest.csv` |

Steps 2 and 3 are independent and may run in parallel after step 1.
Step 4 requires step 1 and the `person_lookup.csv` from step 2.
Step 5 requires steps 2, 3, 4.

## Output

9 CDM v5.4 tables (CSV), one row per schema row, header row matches the CDM
column order:

```
person.csv
observation_period.csv
condition_occurrence.csv
drug_exposure.csv
procedure_occurrence.csv
observation.csv
measurement.csv
specimen.csv
death.csv
```

Plus two auxiliary provenance files:

```
_aux/person_lookup.csv
_aux/measurement_scrna_feature_manifest.csv
```

Loading the CSVs into a typed CDM v5.4 schema (PostgreSQL DDL from
https://github.com/OHDSI/CommonDataModel at the `v5.4` tag) and inserting the
two custom concepts into `cdm.concept` is out of scope for this repo.

## Paths

All paths are env-var driven with repo-relative defaults. Set them once for the
session; every script reads the same vars.

| Env var | Default | Used for |
|---|---|---|
| `VTI_DATA_ROOT` | `./data` | Input data root: Seurat object at `$VTI_DATA_ROOT/meta/integration/05_seurat_harmony_all_new.Rds`, h5ad splits at `$VTI_DATA_ROOT/h5ad/splits/`, clinical RDS at `$VTI_DATA_ROOT/clinical-data/clinical_table_DF_2024_10_28.Rds`, specimen-days RDS at `$VTI_DATA_ROOT/clinical-data/clinical_data_for_seurat.Rds` |
| `VTI_VOCAB_DIR` | `./vocabulary` | `source_to_concept_map.csv` (read by `CERTAINTY_etl_clinical.R`) |
| `VTI_OUTPUT_DIR` | `./output` | CDM CSV outputs + `_aux/person_lookup.csv` (written by clinical/scRNA ETL) |
| `VTI_FEATURE_DIR` | `./output/feature_derivations` | `car_fraction_derivation.csv` (written by `CERTAINTY_derive_car_fraction.R`, read by `CERTAINTY_etl_scrna_features.R`) |
| `VTI_RADE_REPO` | unset | Optional path to a `Rade_et_al_CAR_2025` checkout; if set, provenance pins `git rev-parse HEAD`, else records `NA` |

Inputs are not bundled. Minimal layout to run end-to-end:

```
$VTI_DATA_ROOT/
  meta/integration/05_seurat_harmony_all_new.Rds
  h5ad/splits/
  clinical-data/
    clinical_table_DF_2024_10_28.Rds
    clinical_data_for_seurat.Rds
$VTI_VOCAB_DIR/
  source_to_concept_map.csv
$VTI_OUTPUT_DIR/
$VTI_FEATURE_DIR/
```

Override at the shell, e.g.:

```bash
export VTI_DATA_ROOT=/data/rade-2026
export VTI_OUTPUT_DIR=$VTI_DATA_ROOT/omop-output
Rscript scripts/CERTAINTY_etl_clinical.R
```

## Environment

Singularity image on Zenodo 14732727.

## Contact

For access requests, integration questions, or to report issues, contact the
corresponding author of the publication or contact below:

- christoph.kaempf@izi.fraunhofer.de
- georg.popp@izi.fraunhofer.de
