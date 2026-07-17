# Liver Cancer Genomic Analysis (TCGA-LIHC)

An R-based analysis pipeline that integrates clinical, mutation, and RNA-seq data from TCGA liver hepatocellular carcinoma (LIHC) patients, groups patients by their mutation profiles, and runs survival analysis to identify genomic and clinical factors associated with patient outcomes.

---

## Overview

The project takes three separate TCGA-LIHC datasets, harmonizes them into a single patient cohort, and then runs three layers of analysis:

1. **Clinical exploration** — distributions and survival patterns across demographic and clinical variables.
2. **Mutational analysis** — identifying the most frequently mutated genes and clustering patients by their mutation profiles.
3. **Survival modeling** — Kaplan–Meier analysis with log-rank testing across clinical and molecular variables to find factors linked to prognosis.

---

## Dataset

TCGA Liver Hepatocellular Carcinoma (LIHC), three files:

| File | Contents |
|------|----------|
| `data_clinical_patient.txt` | Patient clinical data (age, sex, tumor stage, survival) |
| `data_mutations.txt` | Somatic mutation calls per patient |
| `RNAseq_LIHC.csv` | RNA-seq expression matrix |

Because the three files use slightly different TCGA barcode formats, custom regex-based functions extract and standardize patient IDs, and the analysis is restricted to the cohort of patients present in **all three** datasets.

---

## Pipeline

### 1. Data integration & cleaning
- Load all three datasets and extract standardized patient IDs from TCGA barcodes.
- Identify the intersection of patients common to clinical, mutation, and RNA-seq data.
- Filter each dataset down to that shared cohort.

### 2. Clinical exploratory analysis
- Age distributions, survival status parsing, and demographic breakdowns.
- Comparative plots (violin, box, bar) of age and survival across sex and diagnosis groups.
- Summary statistics such as mortality rate by subgroup.

### 3. Mutational analysis
- Filter mutations to non-synonymous variants (missense, nonsense, frameshift, splice-site, etc.).
- Rank and visualize the most frequently mutated genes.
- Build a binary gene-by-patient mutation matrix.
- Cluster patients into subgroups using **hierarchical clustering (Ward.D2)** and **k-means** (on mutation burden), visualized with clinically annotated heatmaps.

### 4. Survival analysis
Kaplan–Meier survival curves with log-rank tests across multiple variables:
- **Age group** (median split)
- **Sex**
- **Tumor stage** (Stage I–III)
- **T-stage** (T1–T4)
- **Driver-mutation burden** — using a literature-based LIHC driver-gene panel (TP53, CTNNB1, AXIN1, etc.), patients are grouped by number of mutated driver genes.

Chi-square tests are also used to check whether the mutation-based clusters align with clinical tumor staging.

---

## Key findings

- **Tumor T-stage** was strongly associated with overall survival (log-rank χ² ≈ 35, p ≈ 1×10⁻⁷) — later-stage tumors had markedly worse outcomes.
- **Overall tumor stage** was also significant (p ≈ 0.0001).
- **Age and sex** showed the expected clinical trends but were not statistically significant in this cohort.
- **Driver-mutation burden** was analyzed as a prognostic factor by splitting patients into low (0–1) vs. high (2+) driver-gene groups.

> Disease-specific survival (DSS) was excluded from analysis because the filtered cohort contained no recorded disease-specific death events, making the log-rank statistics undefined. Overall survival (OS) was used as the primary endpoint instead.

---

## Tech stack

- **R**
- **dplyr / tidyr** — data wrangling
- **ggplot2** — visualization
- **survival / survminer** — Kaplan–Meier estimation and log-rank testing
- **pheatmap** — annotated mutation heatmaps

---

## Repository structure

```
.
├── Initial_Analysis.Rmd       # Data integration + clinical exploratory analysis
├── Mutational_Analyisis.Rmd   # Mutation matrix, clustering, cluster survival
├── RNA_Seq_Analysis.Rmd       # Cohort integration + multi-variable survival analysis
└── data/                      # TCGA-LIHC files (not tracked — see below)
```

---

## Getting started

### Requirements
```r
install.packages(c("dplyr", "tidyr", "ggplot2", "survival", "survminer", "pheatmap", "scales"))
```

### Data
Place the three TCGA-LIHC files (`data_clinical_patient.txt`, `data_mutations.txt`, `RNAseq_LIHC.csv`) in the project directory. These are available through [cBioPortal](https://www.cbioportal.org/) / the [GDC Data Portal](https://portal.gdc.cancer.gov/). The data files are **not** included in this repo.

### Run
Open each `.Rmd` in RStudio and knit, or run the chunks top to bottom. Start with `Initial_Analysis.Rmd` (it builds the integrated cohort the other analyses rely on).

---

## Limitations & future work

- **RNA-seq expression is not yet analyzed.** The RNA-seq matrix is currently used only to define the shared patient cohort — no differential expression or dimensionality reduction is performed on it. A natural next step is adding **DESeq2** differential-expression analysis and **PCA** on the expression data to connect gene-expression patterns to the mutation-based subgroups.
- Some subgroups (e.g., Stage IV, T4) have small sample sizes, which limits statistical power.
- Clustering parameters (number of clusters, distance metric) were chosen manually and could be validated more rigorously.
