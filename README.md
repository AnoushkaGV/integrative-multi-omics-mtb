# Blood mRNA and small RNA remodeling during early tuberculosis treatment

This repository contains the analysis code for a whole-blood multi-omics study of the host response during early tuberculosis treatment. Paired samples collected at baseline and after four weeks of standard therapy were analyzed across three South African cohorts. The analyses include an mRNA meta-analysis with external validation, immune cell deconvolution, and miRNA and tRNA-derived fragment (tRF) profiling. A baseline gene-pair signature, EARLY5, was developed to predict sputum culture conversion at week 4.

## Scripts

| Script | Description |
|---|---|
| `01_meta_mrna.R` | mRNA meta-analysis of GSE89403 and the study cohort; validation in GSE31348 |
| `02_cibersort.R` | Immune cell deconvolution using CIBERSORTx and MCP-counter |
| `03_mirna.R` | miRNA differential expression and miRNA–mRNA target analysis |
| `04_trfs.R` | tRF differential expression, anticodon-level and seed-match analyses |
| `05_predictions.R` | EARLY5 development, comparison with published signatures, and clinical models |
| `run_all.R` | Runs all analyses in order |

## Usage

All scripts were run in R 4.5.0 from the repository root:

```r
source("scripts/run_all.R")
```

`01_meta_mrna.R` must be run before `03_mirna.R` and `04_trfs.R`, which use its output. Results are written to `output_revised/`.

## Data

Sequencing data from the study cohort have been deposited in GEO under accession **GSE341867**. GSE89403 and GSE31348 are publicly available from GEO. The scripts expect input files at the following locations:

```
data/mrna/gse89403/counts/              GSE89403 gene counts
data/mrna/gse31348/GSE31348/CEL/        GSE31348 CEL files
data/mrna/study_cohort/counts_paired/   Study cohort gene counts
data/mirna/counts_paired.csv            Study cohort miRNA counts
data/trf_counts.csv                     Study cohort tRF counts
data/trfs/                              Study cohort tRF read reports
metadata/gse89403/SraRunTable.csv       GSE89403 sample metadata
```

This repository also includes the CIBERSORTx input and output files for GSE89403, Human Protein Atlas immune cell expression data, the Ensembl gene symbol map used in the analysis (release 116), and the tRF cluster reference files.

## Contact

For questions regarding the code or data, please contact G V Anoushka Chinmayi (g.v.anoushka@gmail.com).

## License

This project is licensed under the MIT License.
