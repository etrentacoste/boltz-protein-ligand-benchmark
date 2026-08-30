# Benchmark dataset

This directory contains the final 100-complex protein–ligand benchmark and the
intermediate tables generated during benchmark construction.

## Final benchmark

`benchmark_100_systems.csv` contains the 100 experimentally determined
protein–ligand complexes evaluated with Boltz-1 and Boltz-2.

The benchmark combines:

- 30 manually curated pilot systems;
- 70 systems obtained using the automated RCSB selection pipeline.

## Selection workflow

The complete selection history, intermediate classification tables and detailed
selection criteria are provided in [`selection/`](selection/).

The workflow included:

1. RCSB candidate retrieval;
2. metadata-based screening;
3. target-redundancy control;
4. coordinate-level inspection;
5. ligand validation assessment;
6. manual review of ambiguous systems;
7. construction of the final benchmark.

For the exact filtering criteria and counts remaining after each stage, see
[`selection/README.md`](selection/README.md).
