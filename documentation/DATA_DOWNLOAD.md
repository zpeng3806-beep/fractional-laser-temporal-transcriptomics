# Data download and verification

This analysis uses only public, de-identified processed human transcriptomic data and official platform annotations from NCBI GEO.

| Accession | Required local filename | Official source |
|---|---|---|
| GSE168760 | `GSE168760_series_matrix.txt.gz` | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE168760 |
| GSE206495 | `GSE206495_series_matrix.txt.gz` | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE206495 |
| GPL13667 | `GPL13667_platform.soft.tsv` | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL13667 |
| GPL15207 | `GPL15207_platform.soft.tsv` | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL15207 |

Create `data/` at the repository root and save the files under the exact names above. GEO download interfaces may expose the Series Matrix and full platform table through accession-specific links; retrieve the deposited files without modifying their contents.

After download, record SHA-256 checksums locally:

```bash
shasum -a 256 data/*
```

The manuscript analysis did not use raw CEL, FASTQ, or SRA files. Do not substitute such files without defining a new, independently reviewed workflow.

