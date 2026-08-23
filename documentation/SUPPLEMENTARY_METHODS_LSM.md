# Supplementary Methods

## Deposited inputs and cohort-specific scaling

The analysis used the deposited, processed GEO Series Matrix files for GSE168760 and GSE206495. No FASTQ, CEL, or image-level raw data were reprocessed. For GSE168760, deposited intensity values were transformed as `log2(value+1)` exactly as encoded in the frozen analysis. GSE206495 values were analyzed on the deposited scale without an additional logarithmic transformation. All transformations, feature annotation, gene collapse, scaling, and model fitting were conducted separately by cohort.

Within GSE168760, each gene was standardized across retained baseline and single-treatment samples. Within GSE206495, standardization parameters were estimated from D0, D1, D7, and D14 forearm samples and applied within that cohort. Cohort-specific standardization supported effect-ordering comparisons but did not imply identical effect magnitudes or platform sensitivity.

## Platform annotation and shared-gene audit

GEO platform tables for GPL13667 and GPL15207 supplied feature annotation. A probe was retained only when gene symbol and Entrez identifier were both present and unambiguous. Symbols mapping to conflicting Entrez identifiers on a platform were removed. Multiple retained probes assigned to the same symbol were collapsed by mean expression with `avereps` from limma. Reliable platform-specific maps were joined by symbol, and genes were retained only when Entrez identifiers also agreed. The resulting frozen shared universe contained 16,942 genes. Expression matrices were never concatenated or batch-corrected across platforms.

## Program scores and standardized effects

For every retained sample, gene-wise standardized expression values were averaged over the shared genes in each frozen program. Repeated program scores were modeled within cohort using `program score ~ day + participant`, with D0 as reference. Model coefficients and model-based 95% intervals were used only for within-cohort displays.

Cross-cohort vectors used participant-paired Hedges g. For participant-specific treated-minus-D0 differences, the effect was calculated as the mean paired difference divided by the standard deviation of paired differences and multiplied by `1-3/(4n-5)`. The same definition was applied to program summaries and individual shared genes. All primary participants contributed D0, D1, D7, and D14 samples; no primary effect used an unpaired participant or imputation.

## Frozen program constitution

The feasibility protocol named seven categories before pilot execution. Exact MSigDB v2026.1 definitions were encoded for feasibility assessment, and the constitution was frozen without additions or deletions before formal Phase 1B and Phase 2B testing. The program set was not externally preregistered and did not supply an unseen validation cohort.

| Program | Public definition | Shared genes |
|---|---|---:|
| Acute injury | Hallmark TNFA signaling via NF-kappa B | 191 |
| Inflammation | Hallmark inflammatory response | 196 |
| Epidermal repair | GO positive regulation of epithelial cell proliferation involved in wound healing | 9 |
| Extracellular-matrix remodeling | Reactome extracellular matrix organization | 305 |
| Collagen organization | GO collagen fibril organization | 69 |
| Proliferation | Hallmark E2F targets | 183 |
| Matrix maturation | GO extracellular matrix assembly | 43 |

## Participant-cluster bootstrap and exact program-label test

For each of 5,000 bootstrap replicates, participants were sampled with replacement independently within each cohort. All seven program effects and the D1, D7, and D14 cross-cohort Spearman correlations were recalculated. Percentile 95% intervals were derived from valid bootstrap correlations. The exact program-label test evaluated each observed day-specific rho across all 7!=5,040 permutations and was two-sided with respect to absolute correlation.

## Matched-random-program architectures

Ten thousand random seven-program architectures were generated from the 16,942-gene universe. Each random set matched one frozen program's gene count; sampling was without replacement within a set, while genes could recur across different sets. Sample scores, participant-paired effects, and three day-specific correlations were recalculated for each architecture. Primary summaries were the mean and minimum of the D1/D7/D14 correlations. One-sided empirical P values used `(1 + count(null >= observed))/(10,000 + 1)`. Day-specific and transition distributions were secondary diagnostics. The null matched gene counts only and did not preserve gene expression, variance, dependence, functional overlap, or platform-probe properties.

## LOPO influence analysis

One participant was removed at a time from GSE168760 while GSE206495 remained intact, and then one participant was removed at a time from GSE206495 while GSE168760 remained intact. Program effects and all three primary correlations were recalculated after every deletion, producing 14 and 17 dependent influence analyses. The analysis tested whether a single participant reversed cross-cohort concordance; it did not create 31 independent replications.

## Genome-wide bootstrap and gene-label null

For every shared gene, D1-minus-D0, D7-minus-D0, and D14-minus-D0 standardized-expression differences were calculated per participant and summarized as paired Hedges g. Complete 16,942-gene vectors were correlated without significance or magnitude filtering. Genome-wide bootstrap intervals used 5,000 independent within-cohort participant resamples with recalculation of gene effects and rho.

For each day, GSE168760 gene effects were held fixed and GSE206495 gene labels were permuted 10,000 times. The one-sided plus-one empirical test asked whether same-gene alignment exceeded arbitrary gene pairing while preserving each cohort's marginal effect distribution. It did not preserve gene-gene dependence or pathway structure.

## Reciprocal rank transfer

At each day, the top 10% most positive and top 10% most negative genes in one cohort were treated as fixed sets and tested against the complete ranked effect vector in the other cohort. The negative tail was tested after reversing the target ranking so that the prespecified expected normalized enrichment score was positive. Analyses were performed for both tails, both transfer directions, and three days, giving 12 tests. `fgseaSimple` used 10,000 permutations per test. P values were adjusted together across exactly 12 tests with the Benjamini-Hochberg procedure.

## Transition and collagen sensitivity analyses

For each program, the transition effect was D14 paired Hedges g minus D1 paired Hedges g. The two seven-program vectors were compared by Spearman rho and sign concordance; the same 10,000 matched architectures supplied the one-sided transition null. Because observed transition rho did not exceed this null (P=0.286), transition findings remained descriptive/supportive.

The D14 collagen-organization effect was recalculated for the complete 69-gene program, after removing `COL1A1` and `COL1A2`, after removing every symbol beginning with `COL`, and after removing each of 69 genes in turn. The analysis tested major- and single-gene dependence only; it did not test new collagen induction, protein abundance, histological deposition, or clinical remodeling.

## Software and reproducibility

The frozen random seed was 20260823. Formal analyses used R v4.6.1 with data.table v1.18.4, limma v3.68.5, Matrix v1.7-5, matrixStats v1.5.0, fgsea v1.38.0, and msigdbr v26.1.1. The local repository-ready package contains scripts, input metadata, frozen program definitions, derived result tables, figure source data, environment records, and checksums. The package must not be published until Phase 4C authorization and author metadata/licensing decisions are complete.
