# Blood mRNA and small RNA remodeling during early tuberculosis treatment

Analysis code for:

> Chinmayi GVA, Quijada D, Queiroz ATL, Martinson N, Andrade BB, Halushka MK, Bader JS, Karakousis PC. *Blood mRNA and small RNA remodeling during early tuberculosis treatment and a baseline signature of culture conversion.* (manuscript under review)

We profiled paired whole-blood samples collected at baseline and after four weeks of standard anti-tuberculosis therapy across three South African cohorts. The analysis combines an mRNA meta-analysis (GSE89403 and a prospective cohort, with validation in GSE31348), immune cell deconvolution, and miRNA and tRNA-derived fragment (tRF) profiling of the prospective cohort. A baseline k-top scoring pairs signature, EARLY5, was developed in GSE89403 to predict sputum culture conversion by week 4.

## Repository structure

```
.
├── scripts/
│   ├── 01_meta_mrna.R      mRNA meta-analysis (GSE89403 + study cohort, Stouffer), GSE31348 validation
│   ├── 02_cibersort.R      CIBERSORTx / MCP-counter immune deconvolution (GSE89403)
│   ├── 03_mirna.R          miRNA differential expression (DESeq2), multiMiR targets among meta-DEGs
│   ├── 04_trfs.R           tRF differential expression, read lengths, anticodons, 7-mer seed matching
│   ├── 05_predictions.R    EARLY5 k-TSP (LOOCV), Sweeney3/RISK6 comparison, logistic and Cox models
│   └── run_all.R           runs 01-05 in order and saves a log
├── data/
│   ├── mrna/gse89403/cibersort/           CIBERSORTx input (TPM mixture) and output for GSE89403
│   ├── rna_immune_cell.tsv                Human Protein Atlas immune cell expression (nTPM)
│   ├── ensembl_symbol_map_biomart.csv     Ensembl -> HGNC symbol map used (Ensembl release 116)
│   └── ensembl_symbol_map_biomart_release.txt
├── trf-clusters.csv                       tRF cluster reference (Lu et al.)
├── trf_clusters_clean.csv                 tRF cluster reference with fragment type and origin
├── LICENSE
└── README.md
```

## Scripts and outputs

| Script | Manuscript items | Main methods |
|---|---|---|
| `01_meta_mrna.R` | Figure 2a–c, Supplementary Table S1 | edgeR TMM, limma-voom paired designs, Stouffer weighted Z (√n weights), BH; GSE31348 RMA + paired limma, binomial sign test |
| `02_cibersort.R` | Figure 2d, Supplementary Table S8 | CIBERSORTx (LM22, samples with p < 0.05), MCP-counter on TPM, paired Wilcoxon, BH |
| `03_mirna.R` | Figure 4a–c, Supplementary Tables S2–S3 | DESeq2 (~ subject + condition), apeglm, multiMiR (validated and predicted) |
| `04_trfs.R` | Figure 5a–f, Supplementary Table S4 | DESeq2, apeglm, fragment-type-restricted read analysis, RP100K anticodon log2FC, 7-mer seed matching (TxDb/BSgenome hg38) |
| `05_predictions.R` | Figure 3, Table 2, Supplementary Tables S5–S7, Supplementary Figures S1–S2 | switchBox k-TSP nested in LOOCV, pROC (DeLong), TBSignatureProfiler (ssGSEA), logistic regression (LRT, profile CIs, Firth), Cox (Efron and exact ties) |

## Running the analysis

Run from the repository root in R (4.5.0):

```r
source("scripts/run_all.R")
```

`01_meta_mrna.R` must run first: `03_mirna.R` and `04_trfs.R` use its meta-DEG list. All tables and figures are written to `output_revised/<analysis>/{tables,figures}`, and the console output is saved to `output_revised/run_log_*.txt` together with `sessionInfo()`.

Main R/Bioconductor packages: limma, edgeR, DESeq2, apeglm, multiMiR, biomaRt, GEOquery, oligo, org.Hs.eg.db, MCPcounter, switchBox (1.44.0), pROC (1.19.1), TBSignatureProfiler (1.20.0; GSVA 2.2.1), survival (3.8-12), logistf (1.26.1), TxDb.Hsapiens.UCSC.hg38.knownGene, BSgenome.Hsapiens.UCSC.hg38, ComplexHeatmap, ggplot2. CIBERSORTx was run on the [CIBERSORTx web portal](https://cibersortx.stanford.edu/); its input and output are included.

## Input data

Small reference files are included in this repository. Sequencing-derived inputs are not, and should be placed at the paths below.

| Path | Contents | Source |
|---|---|---|
| `data/mrna/gse89403/counts/` | Per-sample gene counts (STAR + featureCounts, GRCh38/GENCODE v47) | Generated from GSE89403 (SRA) |
| `metadata/gse89403/SraRunTable.csv` | GSE89403 run metadata | SRA Run Selector |
| `data/mrna/gse31348/GSE31348/CEL/` | Affymetrix HG-U133 Plus 2.0 CEL files | GEO GSE31348 |
| `data/mrna/study_cohort/counts_paired/` | Study cohort gene counts | GEO GSE341867 |
| `data/mirna/counts_paired.csv` | Study cohort miRNA counts (miRge3.0) | GEO GSE341867 |
| `data/trf_counts.csv` | Study cohort tRF cluster counts (miRge3.0) | GEO GSE341867 |
| `data/trfs/*.tRFs.report.tsv` | Study cohort tRF read-level reports (miRge3.0) | GEO GSE341867 |
| `data/mrna/gse89403/cibersort/` | CIBERSORTx input and output | Included |
| `data/rna_immune_cell.tsv` | HPA immune cell expression | Included (Human Protein Atlas) |
| `data/ensembl_symbol_map_biomart.csv` | Ensembl -> HGNC map (release 116) | Included |
| `trf-clusters.csv`, `trf_clusters_clean.csv` | tRF cluster reference | Included |

## Data availability

| Dataset | Data | Access |
|---|---|---|
| GSE341867 | Prospective study cohort: whole-blood RNA-seq and small RNA-seq (baseline and week 4) | GEO (public on publication) |
| GSE89403 | Whole-blood RNA-seq | GEO / SRA |
| GSE31348 | Whole-blood microarray | GEO |

## Contact

For questions about the data or code, contact g.v.anoushka@gmail.com.

## License

MIT (see `LICENSE`).
