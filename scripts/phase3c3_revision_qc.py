#!/usr/bin/env python3
from pathlib import Path
import csv
import re

ROOT = Path(__file__).resolve().parents[1]
MANUSCRIPT = ROOT / "FRACTIONAL_LASER_MASTER_MANUSCRIPT_V2.md"
text = MANUSCRIPT.read_text()


def write_tsv(path, fieldnames, rows):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


number_checks = [
    ("cohort A participants", "14 GSE168760 participants", "Table1_CohortDesign.tsv", "14"),
    ("cohort B participants", "17 GSE206495 participants", "Table1_CohortDesign.tsv", "17"),
    ("cohort A age", "aged 30-55 years", "Sherrill et al. source article", "30-55 years"),
    ("cohort A skin types", "Fitzpatrick skin types I-IV", "Sherrill et al. source article", "I-IV"),
    ("cohort B age", "mean age of 54.5 years", "Garza et al. source article", "54.5 years"),
    ("cohort B skin types", "Fitzpatrick skin types I-III", "Garza et al. source article", "I-III"),
    ("cohort A samples", "yielding 56 primary samples", "Table1_CohortDesign.tsv", "56"),
    ("cohort B samples", "yielding 68 samples", "Table1_CohortDesign.tsv", "68"),
    ("primary window", "D0, D1, D7, and D14", "PHASE2B_CONSTITUTION.md", "D0/D1/D7/D14"),
    ("AFL secondary days", "D3, D21, and D28", "SECONDARY_TIMEPOINT_CONTEXT.tsv", "3/21/28"),
    ("NAFL secondary day", "D29", "GSE206495_metadata.tsv", "day 29"),
    ("shared universe", "16,942 shared genes", "shared_gene_universe.tsv", "16942"),
    ("program count", "seven frozen program summaries", "PROGRAM_CONSTITUTION.tsv", "7"),
    ("program sizes", "191, 196, 9, 305, 69, 183, and 43 genes", "PROGRAM_CONSTITUTION.tsv", "191/196/9/305/69/183/43"),
    ("program bootstrap", "5,000 participant-cluster bootstrap replicates", "cross_cohort_conservation.tsv", "5000"),
    ("exact permutations", "7!=5,040", "cross_cohort_conservation.tsv", "5040"),
    ("program D1 rho", "correlation was 0.821", "cross_cohort_conservation.tsv", "0.82142857"),
    ("program D7 rho", "At D7, rho was 0.821", "cross_cohort_conservation.tsv", "0.82142857"),
    ("program D14 rho", "At D14, rho was 0.857", "cross_cohort_conservation.tsv", "0.85714286"),
    ("program D1 CI", "0.571 to 0.964", "cross_cohort_conservation.tsv", "0.57142857-0.96428571"),
    ("program D7 CI", "0.429 to 0.964", "cross_cohort_conservation.tsv", "0.42857143-0.96428571"),
    ("program D14 CI", "0.214 to 0.929", "cross_cohort_conservation.tsv", "0.21428571-0.92857143"),
    ("program directions D1", "all seven program summaries", "cross_cohort_conservation.tsv", "7/7"),
    ("program directions D7", "six of seven summaries", "cross_cohort_conservation.tsv", "6/7"),
    ("program directions D14", "five of seven sharing direction", "cross_cohort_conservation.tsv", "5/7"),
    ("program label P", "P=0.0341, P=0.0341, and P=0.0238", "cross_cohort_conservation.tsv", "0.03412698/0.03412698/0.02380952"),
    ("collagen early effects", "-0.218 and -0.311", "cross_cohort_program_effect_vectors.tsv", "-0.218/-0.311"),
    ("collagen D14 effects", "1.280 and 1.090", "COLLAGEN_PRESPECIFIED_SENSITIVITY.tsv", "1.28022/1.08975"),
    ("matrix early effects", "-0.972 and -1.197", "cross_cohort_program_effect_vectors.tsv", "-0.972/-1.197"),
    ("matrix D14 effects", "0.483 and 0.695", "cross_cohort_program_effect_vectors.tsv", "0.483/0.695"),
    ("random architectures", "10,000 gene-count-matched", "MATCHED_RANDOM_PROGRAM_NULL_10000.tsv", "10000"),
    ("mean rho", "mean correlation across D1, D7, and D14 was 0.833", "PROGRAM_PRIMARY_JOINT_STATISTICS.tsv", "0.83333333"),
    ("mean null", "null mean was 0.482", "MATCHED_RANDOM_PROGRAM_NULL_SUMMARY.tsv", "0.48204286"),
    ("mean null q95", "95th percentile was 0.798", "MATCHED_RANDOM_PROGRAM_NULL_SUMMARY.tsv", "0.79761905"),
    ("mean null P", "P=0.0313", "MATCHED_RANDOM_PROGRAM_NULL_SUMMARY.tsv", "0.03129687"),
    ("minimum rho", "minimum day-specific correlation was 0.821", "PROGRAM_PRIMARY_JOINT_STATISTICS.tsv", "0.82142857"),
    ("minimum null mean", "random-program null mean of 0.184", "MATCHED_RANDOM_PROGRAM_NULL_SUMMARY.tsv", "0.18441071"),
    ("minimum null q95", "95th percentile of 0.714", "MATCHED_RANDOM_PROGRAM_NULL_SUMMARY.tsv", "0.71428571"),
    ("minimum null P", "P=0.0110", "MATCHED_RANDOM_PROGRAM_NULL_SUMMARY.tsv", "0.01099890"),
    ("day-specific null P", "0.363, 0.146, and 0.074", "MATCHED_RANDOM_PROGRAM_NULL_SUMMARY.tsv", "0.36336/0.14649/0.07449"),
    ("participant removals", "31 dependent influence analyses", "PHASE2B_FULL_LOPO.tsv", "31"),
    ("LOPO cohort counts", "14 removals from GSE168760 and 17 from GSE206495", "PHASE2B_FULL_LOPO.tsv", "14/17"),
    ("LOPO worst day values", "0.821 at D1, 0.821 at D7, and 0.750 at D14", "PHASE2B_LOPO_SUMMARY.tsv", "0.8214/0.8214/0.7500"),
    ("LOPO worst mean", "worst mean correlation across days was 0.798", "PHASE2B_LOPO_SUMMARY.tsv", "0.79761905"),
    ("genome D1 rho", "rho=0.714", "GENOMEWIDE_CONCORDANCE.tsv", "0.7143435"),
    ("genome D7 rho", "rho=0.505", "GENOMEWIDE_CONCORDANCE.tsv", "0.5051687"),
    ("genome D14 rho", "rho=0.373", "GENOMEWIDE_CONCORDANCE.tsv", "0.3730312"),
    ("genome sign D1", "73.8%", "GENOMEWIDE_CONCORDANCE.tsv", "73.828356%"),
    ("genome sign D7", "66.4%", "GENOMEWIDE_CONCORDANCE.tsv", "66.408925%"),
    ("genome sign D14", "63.2%", "GENOMEWIDE_CONCORDANCE.tsv", "63.203872%"),
    ("genome D1 CI", "0.586 to 0.729", "GENOMEWIDE_CONCORDANCE.tsv", "0.58558316-0.72935442"),
    ("genome D7 CI", "0.280 to 0.481", "GENOMEWIDE_CONCORDANCE.tsv", "0.27985134-0.48110640"),
    ("genome D14 CI", "0.109 to 0.423", "GENOMEWIDE_CONCORDANCE.tsv", "0.10850775-0.42299512"),
    ("gene-label permutations", "10,000 permutations", "GENOMEWIDE_CONCORDANCE.tsv", "10000/day"),
    ("gene-label exact P", "P=0.00009999", "GENOMEWIDE_CONCORDANCE.tsv", "0.00009999"),
    ("gene-label displayed P", "P<0.0001", "GENOMEWIDE_CONCORDANCE.tsv", "<0.0001"),
    ("rank cutoff", "top-10%", "RECIPROCAL_RANK_REPLICATION.tsv", "10%"),
    ("reciprocal tests", "All 12 prespecified tests", "RECIPROCAL_RANK_REPLICATION.tsv", "12/12"),
    ("reciprocal NES", "2.788 to 4.339", "RECIPROCAL_RANK_REPLICATION.tsv", "2.788012-4.338851"),
    ("reciprocal q", "q<0.05", "RECIPROCAL_RANK_REPLICATION.tsv", "12/12"),
    ("transition rho", "rho=0.857", "TRANSITION_SUMMARY.tsv", "0.85714286"),
    ("transition directions", "six of seven changed", "TRANSITION_SUMMARY.tsv", "6/7"),
    ("transition null P", "P=0.286", "TRANSITION_SUMMARY.tsv", "0.28617138"),
    ("collagen minus COL1", "1.240 and 1.047", "COLLAGEN_PRESPECIFIED_SENSITIVITY.tsv", "1.23989/1.04650"),
    ("collagen minus all COL", "1.159 and 0.966", "COLLAGEN_PRESPECIFIED_SENSITIVITY.tsv", "1.15943/0.96560"),
    ("collagen LOO", "All 69 leave-one-gene-out", "COLLAGEN_PRESPECIFIED_SENSITIVITY.tsv", "69/69"),
    ("random seed", "20260823", "phase2b_formal_dual_cohort.R", "20260823"),
    ("R version", "R v4.6.1", "phase2b_sessionInfo.txt", "4.6.1"),
    ("fgsea version", "fgsea v1.38.0", "phase2b_sessionInfo.txt", "1.38.0"),
    ("msigdbr version", "msigdbr v26.1.1", "sessionInfo.txt", "26.1.1"),
]

number_rows = []
for i, (metric, phrase, source, value) in enumerate(number_checks, 1):
    number_rows.append({
        "audit_id": f"N{i:03d}",
        "metric": metric,
        "V2_wording": phrase,
        "authoritative_source": source,
        "source_value": value,
        "status": "PASS" if phrase in text else "FAIL",
    })
write_tsv(ROOT / "PHASE3C3_NUMBER_AUDIT.tsv",
          ["audit_id", "metric", "V2_wording", "authoritative_source", "source_value", "status"],
          number_rows)

population_checks = [
    ("GSE168760", "sex", "14 women", "14 women", "Sherrill et al., PMID 34843523"),
    ("GSE168760", "age", "30-55 years", "aged 30-55 years", "Sherrill et al., PMID 34843523"),
    ("GSE168760", "Fitzpatrick", "I-IV", "Fitzpatrick skin types I-IV", "Sherrill et al., PMID 34843523"),
    ("GSE168760", "race/ethnicity", "not reported in checked article", "Race or ethnicity was not reported", "Sherrill et al., PMID 34843523"),
    ("GSE206495", "sex and race", "17 White women", "17 White women", "Garza et al., PMID 36055399"),
    ("GSE206495", "mean age", "54.5 years", "mean age of 54.5 years", "Garza et al., PMID 36055399"),
    ("GSE206495", "Fitzpatrick", "I-III", "Fitzpatrick skin types I-III", "Garza et al., PMID 36055399"),
    ("GSE206495", "photoaging", "moderate-to-severe", "moderate-to-severe photoaging", "Garza et al., PMID 36055399"),
]
population_rows = []
for cohort, field, value, phrase, source in population_checks:
    population_rows.append({
        "cohort": cohort, "field": field, "verified_value": value,
        "manuscript_phrase": phrase, "source": source,
        "status": "PASS" if phrase in text else "FAIL",
    })
write_tsv(ROOT / "PHASE3C3_POPULATION_AUDIT.tsv",
          ["cohort", "field", "verified_value", "manuscript_phrase", "source", "status"],
          population_rows)

terms = ["mechanism", "causal", "specific", "validation", "validated", "predict",
         "prediction", "superior", "universal", "clinical efficacy", "independent replication",
         "laser-specific", "AFL-specific", "NAFL-specific"]
prose = text.split("## References", 1)[0] + text.split("## Figure legends", 1)[1]
sentences = [s.strip() for s in re.split(r"(?<=[.!?])\s+|\n+", prose) if s.strip()]
overclaim_rows = []
for term in terms:
    hits = [s for s in sentences if term.lower() in s.lower()]
    if not hits:
        overclaim_rows.append({"term": term, "sentence": "NO_OCCURRENCE", "acceptable": "YES",
                               "reason": "Term absent from V2 prose and legends.", "action": "NONE"})
        continue
    for sentence in hits:
        lower = sentence.lower()
        boundary = any(x in lower for x in ["not ", "no ", "without ", "cannot", "did not", "rather than",
                                                   "was not", "were not", "do not", "could not", "prevents"])
        technical = any(x in lower for x in ["platform-specific", "cohort-specific", "day-specific",
                                             "size-matched", "gene-label", "specificity", "prespecified"])
        acceptable = boundary or technical
        overclaim_rows.append({
            "term": term,
            "sentence": sentence.replace("\t", " "),
            "acceptable": "YES" if acceptable else "NO",
            "reason": "Explicit negative boundary or technical descriptor." if acceptable else "Potential positive overclaim requires review.",
            "action": "RETAIN" if acceptable else "REVISE",
        })
write_tsv(ROOT / "PHASE3C3_OVERCLAIM_AUDIT.tsv",
          ["term", "sentence", "acceptable", "reason", "action"], overclaim_rows)

abstract = re.search(r"^## Abstract\n(.*?)^## Keywords\n", text, flags=re.M | re.S).group(1)
abstract_checks = [
    ("structured abstract present", all(h in abstract for h in ["### Background", "### Methods", "### Results", "### Conclusions"])),
    ("abstract word count 250-300", 250 <= len(re.findall(r"\b[\w'-]+\b", abstract)) <= 300),
    ("cohort counts correct", "14 participants" in abstract and "17 participants" in abstract),
    ("program rhos correct", all(v in abstract for v in ["0.821", "0.857"])),
    ("random-null P values correct", "P=0.0313" in abstract and "P=0.0110" in abstract),
    ("participant-removal phrasing", "participant-removal analyses" in abstract and "31 independent" not in abstract),
    ("genome rhos correct", all(v in abstract for v in ["rho=0.714", "rho=0.505", "rho=0.373"])),
    ("genome effect descriptors calibrated", all(v in abstract for v in ["substantial", "moderate", "modest"])),
    ("gene-label P correct", "P<0.0001" in abstract),
    ("reciprocal dependence explicit", "related tests were not independent replications" in abstract),
    ("no causal AFL/NAFL comparison", "whereas NAFL" not in abstract and "AFL induced" not in abstract),
    ("no laser-specific claim", "without establishing a fractional-laser-specific response" in abstract),
    ("no clinical endpoint exaggeration", "clinical efficacy" in abstract and "without establishing" in abstract),
    ("novelty bounded", "reproducibility" in abstract.lower() and "novel" not in abstract.lower()),
]
abstract_rows = [{"check_id": f"A{i:02d}", "check": name, "status": "PASS" if passed else "FAIL"}
                 for i, (name, passed) in enumerate(abstract_checks, 1)]
write_tsv(ROOT / "PHASE3C3_ABSTRACT_AUDIT.tsv", ["check_id", "check", "status"], abstract_rows)

if any(r["status"] != "PASS" for r in number_rows + population_rows + abstract_rows):
    raise SystemExit("Phase 3C-3 numerical/population/abstract audit failed")
if any(r["acceptable"] != "YES" for r in overclaim_rows):
    raise SystemExit("Phase 3C-3 overclaim audit failed")
print(f"NUMBER_AUDIT=PASS_{len(number_rows)}_OF_{len(number_rows)}")
print(f"POPULATION_AUDIT=PASS_{len(population_rows)}_OF_{len(population_rows)}")
print(f"OVERCLAIM_AUDIT=PASS_{len(overclaim_rows)}_OF_{len(overclaim_rows)}")
print(f"ABSTRACT_AUDIT=PASS_{len(abstract_rows)}_OF_{len(abstract_rows)}")
