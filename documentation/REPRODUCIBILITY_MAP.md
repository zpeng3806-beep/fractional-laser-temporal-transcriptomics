# Reproducibility Map

## Frozen chronology

The feasibility stage encoded seven predefined biological categories. The exact public gene-set constitution was frozen after feasibility assessment and before formal Phase 1B/Phase 2B analyses. This chronology was not an external preregistration and did not create an unseen validation cohort.

## Script-to-output map

| Script | Role | Key output family |
|---|---|---|
| `minimal_laser_pilot.R` | Feasibility-stage program constitution and pilot | program definitions and pilot audit |
| `phase1b_temporal_validation.R` | Frozen temporal cross-cohort validation | program effects, shared universe, LOPO, null summaries |
| `phase2b_formal_dual_cohort.R` | Formal two-cohort inference | joint program, genome-wide, reciprocal rank, transition, collagen results |
| `phase3b_generate_all.py` | Locked figure production from frozen tables | Figures 1-5, S1-S8, Tables 1-2, figure source data |
| manuscript-QC scripts | Frozen V1/V2 manuscript checks | numerical, reference, claim, and overclaim audits |

## Required stopping checks

Before interpreting new outputs, confirm Phase 1B reproduction at D1/D7/D14: 0.82142857, 0.82142857, and 0.85714286. Do not change programs, thresholds, shared-gene rules, time points, or the participant-level statistical unit to rescue a mismatch.

## Evidence-to-file map

- Program concordance: `derived_results/cross_cohort_conservation.tsv` and `derived_results/PROGRAM_PRIMARY_JOINT_STATISTICS.tsv`.
- Matched-random null: `derived_results/MATCHED_RANDOM_PROGRAM_NULL_SUMMARY.tsv`.
- LOPO: `derived_results/PHASE2B_FULL_LOPO.tsv` and `derived_results/PHASE2B_LOPO_SUMMARY.tsv`.
- Genome-wide concordance: `derived_results/GENOMEWIDE_CONCORDANCE.tsv`.
- Reciprocal rank transfer: `derived_results/RECIPROCAL_RANK_REPLICATION.tsv`.
- Transition: `derived_results/TRANSITION_SUMMARY.tsv`.
- Collagen sensitivity: `derived_results/COLLAGEN_PRESPECIFIED_SENSITIVITY.tsv`.
- Main-figure source data: `figure_source_data/Figure1_SourceData.tsv` through `Figure5_SourceData.tsv`.

## Prohibited reinterpretations

Do not treat cells, genes, probes, or biopsies as independent participants. Do not infer AFL-versus-NAFL causality, laser specificity, clinical efficacy, clinical prediction, or scar transportability. Transition concordance remains descriptive/supportive because matched-random P=0.286.
