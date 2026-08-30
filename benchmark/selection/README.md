# Benchmark selection tables

This directory contains the intermediate and final tables generated during construction of the protein–ligand benchmark.

## Initial RCSB PDB search

The RCSB PDB was queried on **6 August 2026**. Entries were required to satisfy the following conditions:

- initial release date between **10 January 2025 and 6 August 2026**, inclusive;
- experimental method: **X-ray diffraction**;
- resolution ≤ **2.5 Å**;
- organism annotated as *Homo sapiens*;
- at least one protein polymer entity;
- at least one non-polymer entity.

The query used `rcsb_accession_info.initial_release_date`; therefore, the date range refers to the **initial PDB release date**, not the deposition date.

Results were ordered by initial release date in descending order. The query returned **3,647 PDB entries**. The first 2,000 recent entries were processed through the automated screening pipeline.

The complete machine-readable RCSB query is provided in `stage0_query.json`.

## Selection workflow

| Stage | Description | Number retained as INCLUDE |
|---|---|---:|
| Stage 0 | RCSB search satisfying the initial database-level criteria | 3,647 candidates |
| Stage 1 | First 2,000 candidates in descending initial-release-date order | 2,000 screened |
| Stage 2 | Header-level chemical and structural screening | 293 |
| Stage 3 | Target-diversity filtering, with at most two entries per target signature | 121 |
| Stage 4 | Coordinate-level screening of the experimental complexes | 98 |
| Stage 5 | Experimental ligand validation screening | 89 |
| Stage 6 | Automatic cohort selected for prediction | 70 |
| Final benchmark | Automatic cohort plus 30 previously manually curated systems | 100 |

### Stage 2 classifications

The 2,000 structure headers produced:

| Decision | Number |
|---|---:|
| INCLUDE | 293 |
| RESERVE | 786 |
| MANUAL_REVIEW | 243 |
| EXCLUDE | 678 |

Header screening considered the availability of an organic ligand candidate, ligand size and composition, experimental resolution, R-free where available, possible covalent attachment, the number of protein chains, metals, cofactors, mutations and multiple candidate ligands.

The exact implementation and thresholds are recorded in:

- `scripts/selection/classify_candidate_headers.py`

### Stage 3 target diversity

Only Stage 2 `INCLUDE` entries were admitted to the clean diversity stage.

Protein targets were assigned a target signature, preferentially using UniProt identifiers and otherwise using a sequence-derived signature. A maximum of two structures per target signature was retained.

This produced:

- 293 eligible header-level entries;
- 121 selected structures;
- 172 target-redundancy exclusions;
- 78 distinct target signatures;
- maximum of two structures per target signature.

### Stage 4 coordinate screening

Full coordinate mmCIF files were inspected to identify actual ligand instances and the protein chains forming their experimental binding environments.

The coordinate screen produced:

| Decision | Number |
|---|---:|
| INCLUDE | 98 |
| RESERVE | 19 |
| MANUAL_REVIEW | 0 |
| EXCLUDE | 4 |

This stage checked, among other properties:

- whether the nominated ligand was present in the coordinates;
- ligand occupancy and alternate conformations;
- protein chains contacting the ligand;
- the number and location of ligand instances;
- nearby metals, cofactors and other non-polymer components;
- whether the observed complex could be represented consistently in the model input and experimental reference.

The exact coordinate-based definitions are implemented in:

- `scripts/selection/classify_candidate_coordinates.py`

### Stage 5 ligand validation

RCSB PDB validation reports were used to assess the experimental evidence for each selected ligand instance. The screen considered ligand RSCC, RSR, occupancy and whether the required validation values were available.

The validation stage produced:

| Decision | Number |
|---|---:|
| INCLUDE | 89 |
| RESERVE | 25 |
| MANUAL_REVIEW | 1 |
| EXCLUDE | 6 |

The 89 `INCLUDE` entries passed the final automated ligand-quality criteria. The precise thresholds are implemented in:

- `scripts/selection/classify_ligand_validation.py`

### Stage 6 and final benchmark

`stage6_final_clean_candidates.csv` contains the 89 candidates that passed the complete automated selection workflow.

`stage6_core_70.csv` contains the 70 systems selected from this pool for Boltz-1 and Boltz-2 prediction. These systems formed seven prediction batches of ten complexes.

The final benchmark contained **100 complexes**:

- 70 complexes from the automated selection workflow;
- 30 complexes from the preceding manually curated cohort.

## Files

| File | Description |
|---|---|
| `stage0_query.json` | Exact RCSB Search API query |
| `stage0_candidate_ids.csv` | All 3,647 PDB identifiers returned by the query |
| `stage2_header_classification_1_2000.csv` | Header-level classifications for the first 2,000 candidates |
| `stage3_clean_diverse_1_2000.csv` | Target-diverse clean candidate pool |
| `stage4_clean_coordinates_1_2000.csv` | Coordinate-level classifications |
| `stage5_ligand_validation_1_2000.csv` | Ligand-validation classifications |
| `stage6_final_clean_candidates.csv` | 89 candidates passing all automated filters |
| `stage6_core_70.csv` | 70 automatically selected benchmark systems |

## Data not included

Raw downloaded PDB header files, full coordinate mmCIF files and validation XML reports are not stored in this repository because they can be retrieved from the RCSB PDB using the included scripts and identifiers.

Intermediate download reports and redundancy-rejection logs are also omitted because they are not required to identify the final benchmark systems.
