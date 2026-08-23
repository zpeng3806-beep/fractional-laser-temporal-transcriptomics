# Phase 2B constitution

Locked before new Phase 2B inference on 2026-08-23.

## Scientific question

Do independent longitudinal human fractional-laser cohorts reproduce a conserved temporal architecture from early injury response to delayed extracellular-matrix remodeling?

## Frozen analysis objects

- Cohort A: GSE168760, single-treatment arm; D0/D1/D7/D14 primary and D3/D21/D28 secondary AFL-only context.
- Cohort B: GSE206495, forearm only and first treatment cycle; D0/D1/D7/D14 primary and D29 secondary treatment context only.
- Biological unit: participant; paired within-participant changes from D0.
- Standardized effect: paired Hedges g computed from participant-level score differences.
- Shared universe: the 16,942 reliable, cross-platform symbol-and-Entrez-agreeing genes frozen in Phase 1B.
- Programs: acute injury, inflammation, epidermal repair, ECM remodeling, collagen organization, proliferation, and matrix maturation, frozen without additions or deletions.
- Primary days: D1, D7, D14.
- Primary cross-cohort statistic: Spearman correlation of the seven program effect vectors.
- Random seed: 20260823.

## Frozen inferential additions

- Reproduce the Phase 1B Day 1/7/14 rhos before new analysis; absolute tolerance 0.05 for each anchor.
- Retain the original 2,000 matched-random-program null as the Phase 1B reproduction and add a 10,000-iteration Phase 2B computational sensitivity using the same frozen architecture.
- Use 10,000 gene-label permutations for each genome-wide day.
- Reciprocal top-10% positive and negative rank replication uses one-sided preranked GSEA in both cohort directions and all three primary days; BH correction is across exactly 12 tests.
- Transition is the program effect at D14 minus its effect at D1, with a matched-random-program transition null.
- Bootstrap intervals resample participants independently within each cohort; 5,000 replicates.

## Input SHA256

| Input | SHA256 |
|---|---|
| Complete Phase 2B protocol | `9e6bbd63570f22106cfb16dc866861126e7b7d034fbce80318da479257994079` |
| Phase 1B script | `0f5cc585225ce76cfe4210ffee036785434ccc2e586ceeb83371c6fbd692efaf` |
| GSE168760 series matrix | `bf73540536c3619018260fbd0ebedafa8503ce7553870ec1cfad69d476b57c16` |
| GSE206495 series matrix | `eab43e952ee34f8136c9e81d7173986bd7782c35282e9c007b930819f08a3f9e` |
| GPL13667 platform table | `7781652c02f8e7065cf5a807ebb2d0e022ae13af698f247f5420836fc94d8278` |
| GPL15207 platform table | `708b367671a480bb06984c9c48db24a84acee66e917ae754d47c5ee25466bf96` |
| GSE168760 metadata | `cef8e7a46fe35177f82bd36558617580485da59c26125cfba9fdee46115206c2` |
| GSE206495 metadata | `bab9a620380b08939499c09fed23c60702829d226ddd8499fc2e61818a3848b5` |
| Frozen program constitution | `3585194e6bc4c7b259d0f215a970e2d1851cc93acaf7828f6fa91defcaafc66c` |
| Frozen shared universe | `151cec0ef4ad3a18601e963f63267c3995957ed5bdfee5d57ad7843c5e51a375` |
| Frozen Phase 1B conservation anchors | `42d77dff4cdea17b0ca351b33143fc1be665ebd1a781f7197911ab9d05bf054a` |

## Claim boundary

Allowed: reproducible dual-cohort temporal architecture with participant-level uncertainty and explicit negative/boundary evidence.

Forbidden: AFL-versus-NAFL causal superiority, validated mechanism, clinical prediction, cell-type attribution, GSE320017 or GSE315794 validation, new pathway discovery, manuscript-ready claims, or Phase 3 inference.

`PHASE2B_CONSTITUTION_LOCKED`
