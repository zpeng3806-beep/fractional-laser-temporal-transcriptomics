# Fractional Laser Phase 2B — Formal Dual-Cohort Temporal Architecture Analysis

## Executive decision

**FINAL GATE = `STRONG_PASS_TO_MANUSCRIPT`**

GSE168760 与 GSE206495 两个独立纵向人皮肤队列，在预冻结的 7-program 层和 16,942-gene 全基因组层均支持跨队列时间架构一致性。Phase 1B 三个锚定 rho 完全复现；联合多时点 architecture 超过 10,000 次 gene-count-matched random-program null；31/31 个 participant LOPO 均保持 D1/D7/D14 正相关；全基因组三时点相关均为正，三次 gene-label permutation 均 P<0.0001；12/12 reciprocal rank tests 方向正确且 BH q<0.05；D1→D14 transition rho=0.857，6/7 programs 同方向，三个 delayed-remodeling programs 均在两个 cohort 中增强。

该 Gate 严格按预设 A–J 条件判定。它只授权下一阶段考虑 manuscript，不等于已经撰写、排版或投稿。本阶段未写 manuscript、未选择期刊、未生成正式投稿 Figure、未加入新数据，也未进入 Phase 3。

## Locked scope

- Cohort A：GSE168760，14 名 participants；single-treatment D0/D1/D7/D14 为 primary，D3/D21/D28 仅为 AFL secondary context。
- Cohort B：GSE206495，17 名 participants；forearm、first-cycle D0/D1/D7/D14 为 primary；face 排除，D29 仅作为 second-treatment context。
- Biological/statistical unit：participant；不把 probe、gene 或 sample 当作独立生物学重复。
- Standardized effect：participant-level day-minus-D0 differences 的 paired Hedges g。
- Shared universe：Phase 1B 冻结的 16,942 个 cross-platform symbol-and-Entrez-agreeing genes。
- Programs：7 个 Phase 1A/1B frozen programs，不增删、不换 pathway。
- Random seed：20260823。
- GSE320017：`PERMANENTLY_EXCLUDED_FROM_MANUSCRIPT_EVIDENCE_CHAIN`。
- GSE315794：`NOT_VALIDATION_EVIDENCE`。

## Methods summary

两个 cohort 分别完成 probe mapping、gene-level collapse 和 gene-wise standardization，表达矩阵从未直接合并。程序分数先在各 cohort 内计算，再以 participant 配对差异估计 Hedges g。D1/D7/D14 primary concordance 为两个 7-program effect vectors 的 Spearman rho；不确定性沿用 Phase 1B 的 5,000 次 participant-cluster bootstrap 与全部 7!=5,040 次 program-label exact permutation。

Phase 2B 增加：10,000 次 gene-count-matched random-program architectures、所有 31 名 participants 的双向 LOPO、16,942 个 gene-level paired Hedges g、每时点 10,000 次 gene-label permutation、预设 top-10% positive/negative reciprocal preranked GSEA（12 tests；每项 10,000 permutations；统一 BH）、以及 7-program D14-minus-D1 transition vectors。所有新增定义在查看 Phase 2B 结果前写入 constitution。

## Required 15 answers

### 1. Phase 1B 是否完全复现？

是。D1/D7/D14 的 reproduced rho 与 frozen anchors 完全相同，absolute difference 均为 0，小于预设容差 0.05。`PROGRAM_CONSTITUTION.tsv` 与 `shared_gene_universe.tsv` 的 SHA256 也保持不变。

### 2. D1/D7/D14 rho 最终值？

| Day | Spearman rho | Participant-bootstrap 95% CI | Exact 7! label P | Direction concordance |
|---:|---:|---:|---:|---:|
| D1 | 0.8214 | 0.5714–0.9643 | 0.0341 | 7/7 |
| D7 | 0.8214 | 0.4286–0.9643 | 0.0341 | 6/7 |
| D14 | 0.8571 | 0.2143–0.9286 | 0.0238 | 5/7 |

三项 rho 全部明显为正，且 bootstrap interval 下限均大于 0。

### 3. mean-rho / minimum-rho null？

在 10,000 次 frozen gene-count-matched random-program architectures 中：

- observed mean rho=0.8333；null mean=0.4820；null 95th percentile=0.7976；empirical one-sided P=0.0313；percentile=96.88。
- observed minimum rho=0.8214；null mean=0.1844；null 95th percentile=0.7143；P=0.0110；percentile=98.91。

结论是联合三时点 architecture 超过 random-program expectation。单独 D1/D7/D14 的 matched-random P=0.3634/0.1465/0.0745，不应声称每个独立时点都超过 random-program null。

### 4. LOPO 最差情况？

共 31 个 participant deletions（AFL 14 + NAFL 17）。31/31 次删除后，D1/D7/D14 rho 均保持正值。

- worst D1 rho=0.8214；
- worst D7 rho=0.8214；
- worst D14 rho=0.7500；
- worst mean rho=0.7976；
- worst minimum rho=0.7500。

不存在单一 participant 驱动 primary architecture 的证据。

### 5. Genome-wide D1/D7/D14 effect correlations？

| Day | Genes | Genome-wide rho | Participant-bootstrap 95% CI | All-gene sign concordance | Top-half sign concordance |
|---:|---:|---:|---:|---:|---:|
| D1 | 16,942 | 0.7143 | 0.5856–0.7294 | 73.8% | 84.0% |
| D7 | 16,942 | 0.5052 | 0.2799–0.4811 | 66.4% | 72.4% |
| D14 | 16,942 | 0.3730 | 0.1085–0.4230 | 63.2% | 67.3% |

三时点均为正，且 participant-bootstrap interval 下限均大于 0。相关强度随时间减弱，因此允许的表述是 genome-wide concordance persists but attenuates，不是 effect sizes identical。

### 6. Gene-label permutation？

每个时点固定 AFL gene effects、随机置换 NAFL gene labels 10,000 次。D1、D7、D14 的 empirical one-sided P 均为 1/10,001=0.00009999，observed rho 均位于 null 的 100th percentile（按 10,000 次有限抽样精度报告为 P<0.0001）。因此 3/3 timepoints 明确高于 random gene-label expectation。

### 7. Reciprocal rank replication 通过多少？

12/12 tests 的 oriented NES>0，且 12/12 tests 经全体 12 项 BH 校正后 q<0.05。NES 范围 2.788–4.339，最大 q=0.00117。该层同时覆盖：

- AFL top-10% positive → NAFL；
- AFL top-10% negative → NAFL；
- NAFL top-10% positive → AFL；
- NAFL top-10% negative → AFL；
- 三个 primary days。

这是预设 reciprocal rank replication，不是新的 pathway discovery。

### 8. D1→D14 transition rho？

`rho_transition=0.8571`，超过预设 STRONG threshold 0.60。

但 transition matched-random-program empirical P=0.2862，percentile=71.39。因此可主张的是 **frozen programs 的 transition ordering 在两 cohort 中高度一致**；不可主张这个 transition rho 本身显著超过任意 gene-count-matched transition architecture。

### 9. 多少 programs transition 同方向？

6/7 programs 的 D14-minus-D1 effect 在两个 cohort 中同方向。唯一不一致的是 9-gene shared epidermal-repair program。

- Early-weighted：acute injury 与 inflammation 均在两个 cohort 中从 D1 到 D14 下降（2/2）。
- Delayed-weighted：ECM remodeling、collagen organization、matrix maturation 均在两个 cohort 中从 D1 到 D14 增强（3/3）。

### 10. Collagen 是否仍稳健？

是。

- Full collagen-organization D14 g：AFL 1.2802；NAFL 1.0897。
- 删除 COL1A1+COL1A2：AFL 1.2399；NAFL 1.0465。
- 删除所有 `COL*` prefix genes：AFL 1.1594；NAFL 0.9656。
- 69 个 leave-one-gene-out analyses 全部保持正向；AFL g range 1.2391–1.3162，NAFL 1.0530–1.1209。

因此 collagen organization 是稳健 component，不依赖 COL1A1/COL1A2 或任一单基因。它不是新发现的 collagen induction，也没有 protein-level validation。

### 11. Proliferation 是否仍只是 exploratory？

是，状态保持：`EXPLORATORY_CROSS_COHORT_DIVERGENCE`。

两个 cohort 在 D1 和 D7 均为正；D14 AFL g=+0.3320，NAFL g=-0.8666。该 directional separation 可作为 hypothesis-generating secondary result，但 AFL D14 estimate 仍不确定，不能升级为 confirmatory mechanism、headline、title 或 abstract 主结论。

### 12. Originality collision 是否存在？

未识别到 direct collision，但结论带检索边界。

[GSE168760 original paper](https://pubmed.ncbi.nlm.nih.gov/34843523/) 已发表 early epidermal/inflammatory 与 delayed dermal/ECM/collagen response；[GSE206495 original paper](https://pubmed.ncbi.nlm.nih.gov/36055399/) 已发表 NAFL matrix/collagen/TGF-beta 与 responder lipid/barrier biology；[GSE168760 secondary network reanalysis](https://pubmed.ncbi.nlm.nih.gov/37744012/) 已做 day 7 versus day 1 的 DEG/PPI/GO 分析。检索还发现动物时序转录组 GSE234691，仅记录、未分析。

本次 targeted search 未识别到已同时完成两个人类 cohort、shared universe、participant-blocked effects、frozen multi-time programs、full LOPO、matched-random null、genome-wide concordance、reciprocal ranks 与 transition vectors 的文章。真正保留的 increment 是 **independent dual-cohort human replication of a conserved temporal architecture**，而不是 generic laser wound-healing biology。

### 13. 最终 Gate？

`STRONG_PASS_TO_MANUSCRIPT`

全部预设 A–J 条件通过。该判定不消除所有局限：两 cohort 在 laser modality、anatomical site、platform 和 treatment schedule 上不同；所有 participants 为女性；transition random-program null 不显著；D14 proliferation 分离不是 confirmatory；GSE320017 与 GSE315794 均不提供验证。

### 14. 推荐进入 manuscript 的证据层有哪些？

按当前 Gate，后续若获得 manuscript 阶段授权，正文可优先保留：

1. 双 cohort participant-level design 与严格 exclusions；
2. 7 frozen programs 的 D1/D7/D14 effect architecture；
3. 三时点 cross-cohort rho、bootstrap CI 与 exact label permutation；
4. 10,000 matched-random-program joint mean/min null；
5. full participant LOPO；
6. 16,942-gene genome-wide concordance 与 gene-label permutation；
7. 12 reciprocal rank replications；
8. D1→D14 transition ordering，并同时披露 transition random-null 不显著；
9. delayed ECM/collagen/matrix-maturation strengthening；
10. collagen major-gene-independence sensitivity；
11. originality与非因果 claim boundary。

### 15. 哪些分析必须留在 Supplement？

- 完整 input SHA256、platform mapping 和 shared-universe audit；
- participant program scores/differences 与所有 blocked-model diagnostics；
- 31 个 LOPO 明细和 program-level re-estimates；
- 10,000 matched-random iterations 和 30,000 gene-label null iterations；
- 50,826 行 genome-wide effect table；
- 12 项 reciprocal rank 的完整 ES/NES/P/q；
- transition program table与随机-null明细；
- collagen 69-gene LOO 与 COL-family removals；
- AFL D3/D21/D28 和 NAFL D29 secondary context；
- proliferation divergence；
- GSE315794 negative orthogonal result 与 GSE320017 failed clinical-scar extension；
- originality search matrix 与搜索边界；
- sessionInfo、脚本和 QC manifest。

## Formal A–J Gate audit

| Criterion | Prespecified requirement | Result | Status |
|---|---|---|---|
| A | Phase 1B reproduced; D1/D7/D14 rho all positive | Exact reproduction; 0.821/0.821/0.857 | PASS |
| B | Mean and minimum matched-random null supported | P=0.0313 and 0.0110 | PASS |
| C | All participant LOPO overall positive | 31/31; worst rho=0.750 | PASS |
| D | Genome-wide rho all >0 | 0.714/0.505/0.373 | PASS |
| E | At least 2/3 gene-label permutation positive | 3/3, each P<0.0001 | PASS |
| F | At least 8/12 correct; 6/12 BH q<0.05 | 12/12 and 12/12 | PASS |
| G | Transition rho>=0.60 and >=5/7 same direction | 0.857 and 6/7 | PASS |
| H | >=2 delayed programs strengthen in both | 3/3 | PASS |
| I | Not dependent on major collagen genes | All prespecified removals and LOO positive | PASS |
| J | No direct originality collision | None identified, with search boundary | PASS |

## Figure candidates

Six non-submission `FIGURE_CANDIDATE` PNGs were generated and visually inspected:

1. study design and participant structure;
2. dual-cohort 7-program effect heatmap;
3. D1/D7/D14 cross-cohort effect scatters with rho and CI;
4. 10,000 random-program null and full LOPO;
5. genome-wide concordance and reciprocal ranks;
6. transition architecture and participant-level delayed-program trajectories.

They are analytical candidates, not final submission layouts.

## Reproducibility and negative evidence

The executable script is `02_scripts/phase2b_formal_dual_cohort.R`. It writes all primary, null, genome-wide, reciprocal, transition, sensitivity, context, QC, and figure-candidate outputs from the frozen processed matrices. Raw inputs were not modified. No FASTQ/SRA or new validation dataset was downloaded or analyzed.

Negative/boundary results remain frozen:

- `GSE315794=NOT_VALIDATION_EVIDENCE`;
- `GSE320017=PERMANENTLY_EXCLUDED_FROM_MANUSCRIPT_EVIDENCE_CHAIN`;
- `PROLIFERATION=EXPLORATORY_CROSS_COHORT_DIVERGENCE`;
- `TRANSITION_RANDOM_PROGRAM_NULL=NOT_SIGNIFICANT`。

`FRACTIONAL_LASER_PHASE2B_COMPLETE`
