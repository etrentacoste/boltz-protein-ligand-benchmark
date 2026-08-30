# Benchmark selection scripts

These scripts implement the automated selection pipeline used to identify recent, high-quality human protein–ligand crystal structures for the benchmark.

## Requirements

The selection scripts require Python 3 and the following packages:

```text
gemmi
numpy
```

Internet access is required for the RCSB Search API, coordinate downloads and validation-report downloads.

## Directory layout

The commands below assume the following working directories:

```bash
REPOSITORY_ROOT=/path/to/boltz-protein-ligand-benchmark
SELECTION_SCRIPTS="$REPOSITORY_ROOT/scripts/selection"
TABLE_DIR="$REPOSITORY_ROOT/benchmark/selection"
RAW_DIR=/path/to/untracked_selection_data

mkdir -p \
    "$RAW_DIR/headers" \
    "$RAW_DIR/full_cifs" \
    "$RAW_DIR/validation_reports"
```

Raw downloaded structures and validation reports are intentionally kept outside the version-controlled repository.

## 1. Search the RCSB PDB

```bash
python "$SELECTION_SCRIPTS/search_rcsb_candidates.py" \
    --date-from 2025-01-10 \
    --date-to 2026-08-06 \
    --output-dir "$TABLE_DIR"
```

This writes:

- `stage0_query.json`;
- `stage0_candidate_ids.csv`.

The query returned 3,647 entries ordered by initial release date from newest to oldest.

## 2. Select the first 2,000 recent entries

```bash
python - "$TABLE_DIR/stage0_candidate_ids.csv" \
    "$RAW_DIR/stage1_newest_2000.csv" <<'PY'
import csv
import sys
from pathlib import Path

source = Path(sys.argv[1])
output = Path(sys.argv[2])

with source.open(newline="") as handle:
    rows = list(csv.DictReader(handle))

selected = rows[:2000]

with output.open("w", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=rows[0].keys(),
    )
    writer.writeheader()
    writer.writerows(selected)

print("Selected:", len(selected))
print("Written:", output)
PY
```

The candidates were therefore not randomly sampled. They were the first 2,000 entries after sorting by decreasing initial release date.

## 3. Download PDB header mmCIF files

```bash
python "$SELECTION_SCRIPTS/download_candidate_headers.py" \
    --input "$RAW_DIR/stage1_newest_2000.csv" \
    --output-dir "$RAW_DIR/headers" \
    --workers 8
```

The number of workers controls only download parallelism and does not affect scientific results.

## 4. Perform header-level classification

```bash
python "$SELECTION_SCRIPTS/classify_candidate_headers.py" \
    --input "$RAW_DIR/stage1_newest_2000.csv" \
    --headers "$RAW_DIR/headers" \
    --output "$TABLE_DIR/stage2_header_classification_1_2000.csv"
```

This classifies candidates as `INCLUDE`, `RESERVE`, `MANUAL_REVIEW` or `EXCLUDE`.

The script evaluates the experimental metadata and chemical components, including ligand candidates, protein-chain multiplicity, metals, cofactors, mutations, possible covalent connectivity, resolution and R-free.

## 5. Apply target-diversity filtering

`make_clean_diverse_pool.py` dynamically loads the target-signature functions from `make_diverse_candidate_pool.py`. Both scripts must therefore remain in this directory.

```bash
python "$SELECTION_SCRIPTS/make_clean_diverse_pool.py" \
    --classification \
        "$TABLE_DIR/stage2_header_classification_1_2000.csv" \
    --headers "$RAW_DIR/headers" \
    --signature-script \
        "$SELECTION_SCRIPTS/make_diverse_candidate_pool.py" \
    --output \
        "$TABLE_DIR/stage3_clean_diverse_1_2000.csv" \
    --rejected-output \
        "$RAW_DIR/stage3_clean_redundancy_rejections.csv" \
    --max-per-target 2
```

Only header-level `INCLUDE` entries are admitted. At most two entries are retained for each target signature.

## 6. Download complete coordinate mmCIF files

```bash
python "$SELECTION_SCRIPTS/download_full_candidate_cifs.py" \
    --input "$TABLE_DIR/stage3_clean_diverse_1_2000.csv" \
    --output-dir "$RAW_DIR/full_cifs" \
    --workers 8 \
    --report "$RAW_DIR/full_cif_download_report.csv"
```

## 7. Perform coordinate-level screening

```bash
python "$SELECTION_SCRIPTS/classify_candidate_coordinates.py" \
    --input "$TABLE_DIR/stage3_clean_diverse_1_2000.csv" \
    --cif-dir "$RAW_DIR/full_cifs" \
    --output "$TABLE_DIR/stage4_clean_coordinates_1_2000.csv"
```

This stage inspects the experimentally observed ligand instances, occupancy, alternate conformations, protein-chain contacts and nearby cofactors or metals.

The resulting table records the selected experimental ligand instance and the protein chains required to represent its observed binding environment.

## 8. Download ligand-validation reports

```bash
python "$SELECTION_SCRIPTS/download_validation_reports.py" \
    --input "$TABLE_DIR/stage4_clean_coordinates_1_2000.csv" \
    --output-dir "$RAW_DIR/validation_reports" \
    --report "$RAW_DIR/validation_download_report.csv" \
    --workers 8
```

## 9. Classify experimental ligand quality

```bash
python "$SELECTION_SCRIPTS/classify_ligand_validation.py" \
    --input "$TABLE_DIR/stage4_clean_coordinates_1_2000.csv" \
    --report-dir "$RAW_DIR/validation_reports" \
    --output "$TABLE_DIR/stage5_ligand_validation_1_2000.csv"
```

The classifier uses the ligand instance selected during coordinate screening and evaluates its RCSB validation record, including RSCC, RSR and occupancy.

## 10. Extract candidates passing all automated filters

```bash
python - \
    "$TABLE_DIR/stage5_ligand_validation_1_2000.csv" \
    "$TABLE_DIR/stage6_final_clean_candidates.csv" <<'PY'
import csv
import sys
from pathlib import Path

source = Path(sys.argv[1])
output = Path(sys.argv[2])

with source.open(newline="") as handle:
    rows = list(csv.DictReader(handle))

selected = [
    row for row in rows
    if row["decision"].upper() == "INCLUDE"
]

with output.open("w", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=rows[0].keys(),
    )
    writer.writeheader()
    writer.writerows(selected)

print("Selected:", len(selected))
print("Written:", output)
PY
```

This produces the 89 candidates that passed the automated workflow.

The retained automatic benchmark cohort is recorded in `stage6_core_70.csv`. The remaining 30 systems in the 100-complex benchmark came from the earlier manually curated cohort.

## Manual inspection utilities

The `manual_review/` directory contains exploratory utilities used during manual curation:

- `inspect_pdb_batch.py`: summarises sequences, protein chains and non-polymer components;
- `inspect_complex_batch.py`: examines ligand instances, contacting chains and nearby non-polymer components.

These scripts supported decisions involving:

- repeated protein chains;
- multiple copies of the nominated ligand;
- binding sites formed by more than one chain;
- metals and cofactors;
- solvent or buffer molecules;
- modified or covalently attached ligands;
- discrepancies between structure metadata and the experimental coordinates.

Manual inspection did not replace the recorded automated classifications. It was used to resolve systems whose biological or chemical context could not be determined reliably from header metadata alone.
