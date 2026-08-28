# [Paper Title Here]

Analysis code accompanying: [Author list]. "[Paper title]." *[Journal]*, [year]. [DOI/link when available]

This repository contains the R scripts used to generate the miRNA, mRNA, and tRF
results and figures reported in the manuscript. Raw sequencing data are not
hosted here (see **Data availability** below) — this repo is scripts + small
derived summary tables only.

## Repository structure

```
.
├── scripts/
│   ├── 01_mirna.R          # miRNA DESeq2 DE analysis, heatmap, volcano, multiMiR target mapping
│   ├── 02_cibersort.R      # CIBERSORTx / MCPcounter immune deconvolution, paired baseline-vs-week4 comparison
│   ├── 03_meta_mrna.R      # mRNA meta-analysis across GSE89403, GSE31348, and the study cohort
│   ├── 04_trfs.R           # tRF-level DESeq2 differential expression and annotation
│   └── 05_predictions.R    # k-TSP / elastic-net classifier for sputum culture conversion, ROC validation
├── results/                # small derived tables (DEG lists, top tables) — see below
├── docs/
└── README.md
```

Each script expects a local `data/` directory (not tracked in this repo — see
`.gitignore`) laid out as:

```
data/
├── mirna/
│   └── counts_paired.csv
├── mrna/
│   ├── gse89403/
│   │   ├── counts/                  # per-sample featureCounts output
│   │   └── cibersort/
│   ├── gse31348/
│   │   └── GSE31348/CEL/            # raw CEL files, fetched via GEOquery
│   └── study_cohort/
│       ├── counts_paired/
│       └── cibersort/
├── trf_counts.csv
└── trfs/                            # per-sample tRF report .tsv files
```

and writes intermediate/output tables to `output/meta/tables/` (also not
tracked).

## Pipeline order

Run scripts in the numbered order above. `03_meta_mrna.R` and `01_mirna.R`
both read tables written by earlier steps (`output/meta/tables/*.tsv`,
`output/meta/tables/meta_degs.csv`) — see the header comments in each script
for its specific inputs.

## Data availability

- **GSE89403** (mRNA, RNA-seq) — publicly available on GEO, fetched directly
  in `03_meta_mrna.R` via `GEOquery`/SRA.
- **GSE31348** (mRNA, microarray) — publicly available on GEO, fetched
  directly in `03_meta_mrna.R` via `GEOquery`.
- **Study cohort** (miRNA, mRNA, tRF-seq) — deposited on GEO under accession
  [**GSE341867**](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE341867).
- **CIBERSORTx** results were generated via the CIBERSORTx web portal
  (https://cibersortx.stanford.edu/), which is not open-source and cannot be
  run as part of this pipeline; input mixture files and output fraction
  tables are described in the Methods.

## Software

Developed and run under R [fill in version, e.g. 4.3.x]. Key packages:
DESeq2, limma, edgeR, switchBox, ComplexHeatmap, multiMiR, MCPcounter, pROC,
glmnet, biomaRt, GEOquery, oligo, ggplot2. A full `sessionInfo()` dump can be
added here once finalized — consider committing one as `docs/sessionInfo.txt`
for exact reproducibility.

## Citation

If you use this code, please cite:

> [Full citation once available]

## License

Released under the MIT License — see `LICENSE`.
