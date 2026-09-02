# Blood Based Multi-Omics during early TB Treatment Response

Analysis code for the whole blood multi-omics study conducted to understand host responses during TB infections and therapy. We have characterized mRNAs, miRNAs and tRFs changing after 4 weeks of standard anti-tuberculous therapy. We have proposed EARLY5 as a candidate gene-pair signature for predicting sputum culture conversion at the 4 week timepoint.

**Authors:** G V Anoushka Chinmayi, Darla Quijada, Artur T L Queiroz, Neil Martinson, Bruno B. Andrade, Marc K. Halushka, Joel S. Bader, Petros C. Karakousis

## Repository structure

```
.
├── scripts/
│   ├── 01_mirna.R           miRNA DESeq2 differential expression analysis, multiMiR target mapping (GSE341867)
│   ├── 02_cibersort.R       CIBERSORTx / MCPcounter immune deconvolution, paired baseline-vs-week4 comparison
│   ├── 03_meta_mrna.R       mRNA meta-analysis across GSE89403, GSE31348, and the study cohort (GSE341867)
│   ├── 04_trfs.R            tRF-level DESeq2 differential expression and annotation (GSE341867)
│   └── 05_predictions.R     k-TSP classifier for sputum culture conversion (LOOCV, clinical validation)
│   └── 06_predictions.R     EARLY5 signature score, ROC comparison with RISK6 and DIAG3
├── data/
│   ├── CIBERSORTx_tables/   CIBERSORTx inputs and outputs
│   ├── rna_immune_cell.tsv  Immune cell transcript expression from the Human Protein Atlas
└── README.md
```

## Requirements and setup

For raw sequencing-derived inputs (count matrices, CEL files, tRF reports) see **Data availability** below.

## CIBERSORTx data

The exact input and output files from CIBERSORTx runs for GSE89403 and the study cohort are included below.

| File | Dataset | Contents |
|---|---|---|
| `CIBERSORTx_GSE89403_results.csv` | GSE89403 | CIBERSORTx deconvolution output |
| `cibersortx_GSE89403_input.tsv` | GSE89403 | Input TPM mixture matrix|
| `CIBERSORTx_study_cohort_results.csv` | Study cohort | CIBERSORTx deconvolution output |
| `cibersortx_study_cohort_input.tsv` | Study cohort | Input TPM mixture matrix |

## Data availability

| Dataset | Type | Access |
|---|---|---|
| GSE89403 | mRNA, RNA-seq | Public — GEO (fetched in `03_meta_mrna.R` via `GEOquery`/SRA) |
| GSE31348 | mRNA, microarray | Public — GEO (fetched in `03_meta_mrna.R` via `GEOquery`) |
| GSE341867 | Study cohort — miRNA, mRNA, tRF-seq | Public — GEO |
| CIBERSORTx results | Immune cell deconvolution | Generated via the [CIBERSORTx portal](https://cibersortx.stanford.edu/); inputs/outputs provided in `data/CIBERSORTx_tables/` |

For additional information on data, metadata, or the analysis code, please contact g.v.anoushka@gmail.com.
