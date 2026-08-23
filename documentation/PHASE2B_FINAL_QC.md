# Phase 2B final QC

`FINAL_QC=PASS`

## Structural and numerical checks

- Phase 1B reproduction: 3/3 PASS; absolute differences 0/0/0.
- Genome-wide effects: 50,826 rows; no missing day/gene/AFL-effect/NAFL-effect fields.
- Genome-wide null: 30,000 rows; exactly 10,000 per primary day.
- Matched-random programs: 10,000 unique iterations.
- Reciprocal ranks: exactly 12 tests; 12/12 direction-correct and 12/12 BH-supported.
- Full LOPO: exactly 31 deletions; 31/31 retained positive D1/D7/D14 rho.
- Transition: seven programs; 6/7 same-direction transitions.
- Collagen sensitivity: 144 rows; all D14 effects positive.
- Evidence table: exact required seven-column schema; 14 evidence layers.
- Formal report: all 15 required answers; completion marker present.
- Figure candidates: exactly six non-empty PNGs, dimensions 1700–2200 pixels wide, visually inspected.
- Output manifest: 31 entries and includes the formal report, evidence table and originality matrix.

Independent QC script result: `PASS_18_OF_18`.

## Core SHA256

| File | SHA256 |
|---|---|
| `PHASE2B_FORMAL_DUAL_COHORT_ANALYSIS.md` | `c007da3f805f15cd633c79414d0187e8e8c3703782f09751eb85ef07aab601ab` |
| `PHASE2B_EVIDENCE_TABLE.tsv` | `5cc8df22d226ef654e6774b66abacf106e24e9635c58df24ca548e58e7c73da4` |
| `PHASE2B_ORIGINALITY_MATRIX.md` | `13c422e977e2e01cb3d9e33b3a99dd7948dd1c152373dea5b1feee83d3e0a1c8` |
| `phase2b_formal_dual_cohort.R` | `bc43e60c6b6d59fd421919c925b43bdaf4a75ee9c90d941fa8c2ef1ce6ec1df1` |
| `PHASE2B_OUTPUT_MANIFEST.tsv` | `921e541b785240812b05b6b7dcf6346cc1ea6e6247cd714201e81ef7eac8f3c0` |

## Boundary check

- GSE320017 reanalysis: NO.
- GSE315794 validation use: NO.
- New dataset analysis: NO.
- Manuscript/abstract/journal selection: NO.
- Submission-figure finalization: NO.
- Phase 3: NOT AUTHORIZED.

`FRACTIONAL_LASER_PHASE2B_COMPLETE`
