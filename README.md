# MAPK/ERK Pathway Analysis in Colorectal Cancer
### A Machine Learning Approach to Drug Target Identification

## Overview

This project uses machine learning to identify potential drug targets within the MAPK/ERK signalling pathway in colorectal cancer. Rather than the conventional approach of predicting tumour vs. normal from gene expression data, this project reframes the problem as **pathway membership prediction** — using the 84 Oncogenic MAPK signalling genes (Reactome) as the target class and the remaining ~60,000 genes in the transcriptome as predictor features.

Genes whose expression patterns are predictive of MAPK pathway activity are identified as candidate drug targets or pathway regulators.

## Dataset

- **Source:** TCGA-COAD (colon adenocarcinoma) + TCGA-READ (rectal adenocarcinoma) via [GDC Portal](https://portal.gdc.cancer.gov/)
- **Samples:** 698 (647 Primary Tumour, 51 Solid Tissue Normal)
- **Genes:** 60,660 raw → 27,384 after filtering (classification) / 20,830 (regression)
- **Pathway genes:** 84 Oncogenic MAPK signalling genes from [Reactome](https://reactome.org/)

Raw TCGA data is not included in this repository due to size. The R preprocessing script downloads it automatically via TCGAbiolinks.

## Pipeline

```
TCGA Data Download (TCGAbiolinks)
        ↓
Sample Filtering (remove Metastatic/Recurrent)
        ↓
Log2 Normalisation → Ensembl-to-Symbol Conversion → Low Expression Filtering
        ↓
┌───────────────────────────┐     ┌──────────────────────────────┐
│   Classification Approach │     │     Regression Approach      │
│                           │     │                              │
│ Transpose matrix          │     │ Compute mean MAPK expression │
│ Label: MAPK / Other       │     │ Variance-based top 500 genes │
│ InfoGain top 100 features │     │ ANN regression (MLP)         │
│ RF & NB ± ClassBalancer   │     │ 100/500/1000 epochs          │
│ 10-fold cross-validation  │     │ 10-fold cross-validation     │
└───────────────────────────┘     └──────────────────────────────┘
        ↓                                     ↓
        └──────────── Pearson Correlation ─────┘
                          ↓
              Ranked Candidate Gene List
```

## Key Results

### Classification (WEKA 3.9.6)

| Classifier | Balancing | Accuracy | MAPK Recall | AUC | Kappa |
|---|---|---|---|---|---|
| Random Forest | None | 99.70% | 0.000 | 0.681 | 0.000 |
| Random Forest | ClassBalancer | 50.61% | 0.012 | 0.746 | 0.012 |
| Naive Bayes | ClassBalancer | 73.80% | **0.866** | **0.776** | **0.476** |
| Naive Bayes | None | 62.13% | 0.866 | 0.767 | 0.008 |

### Regression (WEKA 3.9.6 — MultilayerPerceptron)

| Features | Epochs | Correlation | MAE | RMSE | RAE |
|---|---|---|---|---|---|
| 500 | 100 | 0.7143 | 0.4239 | 0.6351 | 61.78% |
| 500 | **500** | **0.7755** | **0.3805** | **0.5775** | **55.47%** |
| 500 | 1000 | 0.7660 | 0.3795 | 0.5838 | 55.32% |

### Top Drug Target Candidates (Pearson correlation > 0.5)

18 genes identified, with the most biologically credible being:
- **AGTR1** — Angiotensin II Receptor Type 1; known activator of MAPK/ERK signalling via the RAS-RAF-MEK-ERK cascade
- **HAND2-AS1** — long non-coding RNA that directly binds ERK and reduces its phosphorylation in CRC
- **MYH11** — smooth muscle myosin associated with CRC prognosis

## Tools and Packages

### R (preprocessing)
- `TCGAbiolinks` — TCGA data download and preparation
- `SummarizedExperiment` — genomics data container
- `org.Hs.eg.db` — Ensembl to gene symbol conversion
- `tidyverse` — data manipulation

### WEKA 3.9.6 (machine learning)
- InfoGainAttributeEval + Ranker (feature selection)
- ClassBalancer (class imbalance handling)
- NaiveBayes (classification)
- RandomForest (classification)
- MultilayerPerceptron (ANN regression)

## Repository Structure

```
├── README.md
├── scripts/
│   ├── 01_data_download.R          # TCGA data acquisition via TCGAbiolinks
│   ├── 02_preprocessing.R          # Filtering, normalisation, ID conversion
│   ├── 03_classification_prep.R    # Matrix reframing, MAPK/Other labelling
│   ├── 04_regression_prep.R        # Mean MAPK target variable, variance selection
│   └── 05_correlation_analysis.R   # Pearson correlation and candidate ranking
├── weka/
│   └── WEKA_configurations.md      # Full WEKA settings for reproducibility
├── data/
│   └── reactome_mapk_genes.tsv     # 84 Oncogenic MAPK pathway genes from Reactome
├── results/
│   ├── mapk_candidate_genes_ranked.csv
│   ├── mapk_strong_candidates.csv
│   └── mapk_negative_candidates.csv
└── figures/
    └── regression_scatter_plot.png
```

## Future Work

- Network inference (ARACNE / GENIE3) to determine directionality of regulatory relationships
- Graph theory analysis to classify candidates as bottlenecks, hubs, or master regulators
- Druggability and toxicity assessment via DGIdb and ChEMBL

## References

- Hall, M. et al. (2009) 'The WEKA Data Mining Software: An Update', *SIGKDD Explorations*, 11(1), pp. 10–18.
- Li, Q. et al. (2024) 'Signaling pathways involved in colorectal cancer', *Signal Transduction and Targeted Therapy*, 9(1), p. 266.
- Zelli, V. et al. (2022) 'Concurrent RAS and RAS/BRAF V600E Variants in Colorectal Cancer', *Frontiers in Oncology*, 12, p. 863639.

## Author

**Fathimath Nubla Latheef**

OmicsLogic AI-ML for Omics Data Analysis Programme — Capstone Project (2026)

Mentor: Prof Graham Ball
