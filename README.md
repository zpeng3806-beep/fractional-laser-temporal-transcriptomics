# Fractional-laser temporal transcriptomic reproducibility

## Purpose

This repository provides the frozen code, metadata, biological-program definitions, derived results, and figure source data for the manuscript **“Reproducible temporal remodeling programs across independent longitudinal human fractional-laser transcriptomic cohorts.”** The study tests whether two independent longitudinal human skin cohorts show conserved participant-level temporal transcriptomic programs after fractional-laser exposure.

Repository URL: <https://github.com/zpeng3806-beep/fractional-laser-temporal-transcriptomics>

## Public cohorts and design

- [GSE168760](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE168760): longitudinal human fractional-laser cohort.
- [GSE206495](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE206495): independent longitudinal human fractional-laser cohort.
- Primary time points: D0, D1, D7, and D14.
- Statistical unit: participant. Probes, genes, biopsies, and time points are not treated as independent biological replicates.

Only deposited processed Series Matrix expression data and official GEO platform/metadata records are required. Raw CEL, FASTQ, and SRA files were not used for the reported analyses and are not distributed here.

## Frozen analytical program

1. Process each cohort separately and audit platform annotation.
2. Construct the 16,942-gene symbol-and-Entrez-agreeing shared universe.
3. Estimate participant-blocked D1, D7, and D14 effects relative to D0.
4. Score seven frozen biological programs.
5. Test cross-cohort program concordance.
6. Run the matched-random program-architecture null and leave-one-participant-out (LOPO) analysis.
7. Evaluate genome-wide concordance, gene-label permutation, and reciprocal-rank replication.
8. Run the prespecified transition and collagen-sensitivity analyses.
9. Generate locked figures from frozen tables.

The programs, participant definitions, time points, random seed (`20260823`), and stopping rules are frozen. No cross-platform matrix integration, pathway fishing, post hoc program replacement, or analysis-method switching is part of this workflow.

## Software environment

- R session records: `environment/phase1b_sessionInfo.txt` and `environment/phase2b_sessionInfo.txt`
- Python figure generation: `scripts/phase3b_generate_all.py`
- Key R packages: `data.table`, `limma`, `Matrix`, `matrixStats`, and `fgsea`
- Key Python packages: `pandas`, `numpy`, `matplotlib`, and `seaborn`

The recorded session files are the authoritative version snapshot. Exact package restoration may require an R/Python environment compatible with those records.

## Data download

The repository does not redistribute GEO expression matrices or the full GPL platform tables. Download them from NCBI GEO and place them in a local `data/` directory using these exact names:

```text
data/
├── GSE168760_series_matrix.txt.gz
├── GSE206495_series_matrix.txt.gz
├── GPL13667_platform.soft.tsv
└── GPL15207_platform.soft.tsv
```

Official accession pages and local roles are recorded in `metadata/INPUT_FILE_MAP.tsv`. Download commands and verification guidance are in `documentation/DATA_DOWNLOAD.md`.

## Reproduction

From the repository root:

```bash
bash scripts/run_reproduction.sh
```

The wrapper verifies required input filenames, creates `reproduction_output/`, and runs the formal scripts in order. To run manually:

```bash
Rscript scripts/phase1b_temporal_validation.R
Rscript scripts/phase2b_formal_dual_cohort.R
python3 scripts/phase3b_generate_all.py
```

Before interpreting Phase 2B, confirm that the reproduced D1/D7/D14 program-level Spearman correlations are `0.82142857`, `0.82142857`, and `0.85714286`. A mismatch beyond the prespecified tolerance is a stopping condition; do not alter programs, thresholds, or methods to rescue it.

`minimal_laser_pilot.R` is retained only to document the feasibility-stage chronology. It is not required to reproduce the frozen formal analysis because the final program definitions are provided in `program_definitions/`.

## Directory structure

- `environment/`: recorded software environments.
- `scripts/`: analysis, figure-production, and manuscript-QC scripts.
- `metadata/`: cohort metadata, input map, mapping audit, and shared gene universe.
- `program_definitions/`: frozen program constitution and genes.
- `derived_results/`: selected frozen result tables supporting the manuscript.
- `figure_source_data/`: source data behind Figures 1–5.
- `documentation/`: constitution, methods, download guide, reproducibility map, availability statement, and security audit.

Generated outputs are written to the ignored `reproduction_output/` directory. Public input files remain in the ignored `data/` directory.

## Output interpretation

- Program concordance: `derived_results/cross_cohort_conservation.tsv`
- Matched-random null: `derived_results/MATCHED_RANDOM_PROGRAM_NULL_SUMMARY.tsv`
- LOPO: `derived_results/PHASE2B_FULL_LOPO.tsv`
- Genome-wide concordance: `derived_results/GENOMEWIDE_CONCORDANCE.tsv`
- Reciprocal rank: `derived_results/RECIPROCAL_RANK_REPLICATION.tsv`
- Transition analysis: `derived_results/TRANSITION_SUMMARY.tsv`
- Collagen sensitivity: `derived_results/COLLAGEN_PRESPECIFIED_SENSITIVITY.tsv`
- Main-figure source data: `figure_source_data/Figure1_SourceData.tsv` through `Figure5_SourceData.tsv`

## Citation

Use `CITATION.cff`. If citing the biological findings, cite the associated manuscript once published. No DOI is assigned to this repository at release.

## License and data rights

The MIT License applies only to the authors’ original code in this repository. GSE168760, GSE206495, and their associated GEO/platform records remain governed by NCBI GEO and the original studies’ applicable terms. This repository does not claim ownership or copyright in GEO source data.

## Limitations

- Both cohorts are small, longitudinal human datasets.
- Cohort procedures differ, so the analysis addresses temporal program conservation rather than treatment equivalence.
- Transition concordance is descriptive/supportive because the matched-random null result is `P=0.286`.
- Proliferation is exploratory.
- GSE320017 and GSE315794 are excluded from positive manuscript inference.
- The analysis does not establish laser specificity, AFL-versus-NAFL causality, clinical efficacy, clinical prediction, or universal transportability.
