# PLINDER reference sets

This directory contains the PLINDER system and PDB-entry manifests used as
training-set reference universes in the similarity analysis.

These manifests should not be interpreted as official training manifests
released by the Boltz developers. They are PLINDER-based proxy reference sets
defined using the reported temporal cutoffs associated with Boltz-1 and
Boltz-2. They were used to determine which experimentally released structures
could plausibly have been available before each model cutoff.

## Reference-set definitions

| Reference | Temporal rule | PLINDER systems | Unique PDB entries | Last included release date |
|---|---:|---:|---:|---:|
| Boltz-1 | release date `< 2021-09-30` | 167,997 | 53,487 | 2021-09-29 |
| Boltz-2 | release date `< 2023-06-01` | 200,678 | 62,659 | 2023-05-31 |

The cutoff operator is strictly `<`. Therefore, structures released on the
cutoff date itself were not included.

The underlying reference data were derived from PLINDER release
`2024-06/v2`.

## Files

### Boltz-1 reference

- `boltz1_reference_systems.txt.gz`: one PLINDER system identifier per line.
- `boltz1_reference_pdb_ids.txt.gz`: one unique PDB identifier per line.
- `boltz1_reference_summary.csv`: cutoff, system counts, PDB counts and
  release-date range.

### Boltz-2 reference

- `boltz2_reference_systems.txt.gz`: one PLINDER system identifier per line.
- `boltz2_reference_pdb_ids.txt.gz`: one unique PDB identifier per line.
- `boltz2_reference_summary.csv`: cutoff, system counts, PDB counts and
  release-date range.

A PLINDER system identifier represents a specific protein-ligand system.
Multiple PLINDER systems may originate from the same PDB entry; consequently,
the number of systems is larger than the number of unique PDB identifiers.

## Construction

The Boltz-1 reference was recovered from the distinct target systems in the
PLINDER-derived similarity table after applying:

```text
target_release_date < 2021-09-30
```

The Boltz-2 reference was constructed by applying:
```text
target_release_date < 2023-06-01
```

The Boltz-2 reference contains all 53,487 PDB entries represented in the
Boltz-1 reference plus 9,172 additional PDB entries released before the later
cutoff.

The compressed manifests are provided to make the reference populations
auditable without distributing the complete PLINDER dataset, downloaded
structures, Foldseek databases or MMseqs2 databases.

## Role in the similarity analysis

The manifests define the candidate reference populations. Similarity
calculation was performed using PLINDER-derived protein, pocket, interaction
and ligand metrics. Foldseek was used as a structural prefilter before more
detailed comparisons were calculated.

The detailed similarity workflow and processed result tables are documented
in the repository's similarity-analysis directories.
