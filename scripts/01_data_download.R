# =============================================================================
# 01_data_download.R
# MAPK/ERK Pathway Analysis in Colorectal Cancer
# Step 1: Download TCGA-COAD and TCGA-READ RNA-Seq data via TCGAbiolinks
# =============================================================================

# --- Install required packages (run once) ---
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("TCGAbiolinks")
BiocManager::install("SummarizedExperiment")
BiocManager::install("org.Hs.eg.db")
install.packages("tidyverse")

# --- Load libraries ---
library(TCGAbiolinks)
library(SummarizedExperiment)
library(tidyverse)

# --- Query TCGA-COAD (colon adenocarcinoma) ---
query.coad <- GDCquery(
  project = "TCGA-COAD",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

# --- Query TCGA-READ (rectal adenocarcinoma) ---
query.read <- GDCquery(
  project = "TCGA-READ",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

# --- Download data files ---
# Downloads 524 COAD + 177 READ files into GDCdata subfolder
# Files are cached locally so this only needs to run once
GDCdownload(query.coad)
GDCdownload(query.read)

# --- Prepare (merge individual files into SummarizedExperiment objects) ---
data.coad <- GDCprepare(query.coad)
data.read <- GDCprepare(query.read)

# --- Extract expression matrices and merge ---
matrix.coad <- assay(data.coad, "unstranded")
matrix.read <- assay(data.read, "unstranded")
combined.matrix <- cbind(matrix.coad, matrix.read)

# Verify: should be 60660 genes x 701 samples
dim(combined.matrix)

# --- Extract and merge sample metadata ---
coldata.coad <- as.data.frame(colData(data.coad))
coldata.read <- as.data.frame(colData(data.read))
combined.coldata <- rbind(
  coldata.coad[, c("sample_type", "project_id")],
  coldata.read[, c("sample_type", "project_id")]
)

# Check sample type distribution
table(combined.coldata$sample_type)
# Expected: Primary Tumor=647, Solid Tissue Normal=51, Metastatic=1, Recurrent Tumor=2
