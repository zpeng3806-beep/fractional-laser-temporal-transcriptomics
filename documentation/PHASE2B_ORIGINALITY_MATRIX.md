# Phase 2B originality collision matrix

Search frozen on 2026-08-23. The audit used accession-, DOI-, title-, and concept-level searches across PubMed-indexed records, GEO, publisher/PMC full text, and general web indexing. It was designed to detect a direct collision with the Phase 2B scientific question, not to prove that no related publication exists.

## Matrix

| Source | What is already published | Classification | Phase 2B relationship |
|---|---|---|---|
| [GSE168760 original paper, PMID 34843523](https://pubmed.ncbi.nlm.nih.gov/34843523/) | In 14 women, AFL caused a time-dependent human-skin response: early inflammatory/epidermal changes and later dermal/ECM/collagen remodeling; repeated-treatment effects and reversal of age-associated expression were also reported. | `ALREADY_PUBLISHED` | Laser-induced wound healing, inflammation, collagen/ECM change, temporal response, and rejuvenation are not new claims. This paper did not independently validate the architecture in GSE206495. |
| [GSE206495 original paper, PMID 36055399](https://pubmed.ncbi.nlm.nih.gov/36055399/) | In 17 participants, NAFL-associated MMP, collagen, extracellular-component, TGF-beta, dsRNA, and retinoic-acid responses were reported; fast responders showed stronger lipid/barrier programs. | `ALREADY_PUBLISHED` | NAFL collagen/ECM response and fast-responder lipid/barrier biology are not new. The present analysis does not repeat or claim the responder contrast. |
| [Known GSE168760 network reanalysis, PMID 37744012](https://pubmed.ncbi.nlm.nih.gov/37744012/) | Compared day 7 with day 1 and repeated-treatment contexts using DEG, PPI, GO, and KEGG analyses. | `ALREADY_PUBLISHED` | Time-dependent reanalysis of GSE168760 alone is not novel. It did not perform a dual-cohort participant-blocked shared-universe validation. |
| Targeted search for a GSE206495 secondary reanalysis | No separate indexed secondary reanalysis of GSE206495 was identified in the accession/title/DOI searches. The Garza et al. JID article is the original cohort report. | `PARTIAL_OVERLAP` search boundary | Absence is not proof. Any later record not indexed by these searches would require reassessment before submission. |
| [GSE234691 / guinea-pig fractional CO2 time-course](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE234691) | An animal RNA-seq study reported time-varying barrier, collagen, vascular, immune, and defense responses after fractional CO2 laser. | `PARTIAL_OVERLAP` | It overlaps the generic temporal biology but is not an independent human participant-level validation. The dataset was recorded only and was not analyzed, as required by the Phase 2B no-new-dataset boundary. |
| Current Phase 2B analysis | Common 16,942-gene universe; separate participant-blocked effects in two independent human cohorts; seven frozen cross-cohort programs; D1/D7/D14 concordance; full participant LOPO; 10,000 matched-random-program null; genome-wide effect concordance; 10,000 gene-label null; 12 reciprocal rank tests; D1-to-D14 transition-vector replication. | `OUR_INCREMENT` | The increment is the independent, participant-robust, multi-layer dual-cohort validation of temporal architecture. It is not the component biology, an AFL-versus-NAFL comparison, a mechanism, or a clinical predictor. |

## Collision decision

No directly identified publication completed the Phase 2B core combination of both GSE168760 and GSE206495 with the frozen shared-universe, participant-blocked, multi-time program architecture, full LOPO, matched-random-program null, genome-wide concordance, reciprocal rank replication, and transition-vector validation.

`DIRECT_ORIGINALITY_COLLISION=NO_IDENTIFIED_COLLISION_WITH_SEARCH_BOUNDARY`

`STOP_ORIGINALITY=NO`

## Frozen claim boundary

Allowed increment: **independent dual-cohort human replication of a conserved temporal architecture from early injury/inflammatory response to delayed ECM remodeling**.

Forbidden novelty claims: that fractional laser activates wound healing, inflammation, collagen, ECM, rejuvenation, or time-dependent repair; that AFL is superior to NAFL; that the computational programs establish mechanism; or that the analysis predicts clinical response.
