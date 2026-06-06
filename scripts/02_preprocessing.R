# =============================================================================
# 02_preprocessing.R
# MAPK/ERK Pathway Analysis in Colorectal Cancer
# Step 2: Sample filtering, log normalisation, Ensembl-to-symbol conversion,
#         low expression filtering
# =============================================================================
# Requires: data.coad, data.read from 01_data_download.R

library(TCGAbiolinks)
library(SummarizedExperiment)
library(org.Hs.eg.db)
library(tidyverse)

# --- Filter samples ---
# Remove Metastatic (1) and Recurrent Tumour (2) as biologically distinct
keep.samples <- combined.coldata$sample_type %in%
  c("Primary Tumor", "Solid Tissue Normal")

filtered.matrix <- combined.matrix[, keep.samples]
filtered.coldata <- combined.coldata[keep.samples, ]

# Verify: should be 60660 x 698
dim(filtered.matrix)
table(filtered.coldata$sample_type)

# --- Log2 normalisation ---
# Add pseudocount of 1 to avoid log(0)
# Compresses the range and brings data closer to normal distribution
log.matrix <- log2(filtered.matrix + 1)

# --- Transpose so samples are rows, genes are columns (WEKA format) ---
expression.df <- as.data.frame(t(log.matrix))

# --- Convert Ensembl IDs to HGNC gene symbols ---
ensembl.ids <- gsub("\\..*", "", colnames(expression.df))
gene.symbols <- mapIds(org.Hs.eg.db,
                       keys = ensembl.ids,
                       column = "SYMBOL",
                       keytype = "ENSEMBL",
                       multiVals = "first")

# Only keep columns that successfully converted to gene symbols
# This removes non-coding RNAs and pseudogenes without standard annotations
converted <- !is.na(gene.symbols[seq_along(gene.symbols)])
expression.df <- expression.df[, c(which(converted), ncol(expression.df))]
colnames(expression.df) <- c(gene.symbols[converted], "class")
cat("Genes after removing unconverted Ensembl IDs:", sum(converted), "\n")

# --- Load Reactome MAPK pathway gene list ---
mapk.genes <- read.delim("Gene list.tsv", header = TRUE)
mapk.gene.list <- mapk.genes$Gene.name

# --- Handle duplicate column names ---
# Multiple Ensembl entries can map to the same gene symbol
sum(duplicated(colnames(expression.df)))
colnames(expression.df) <- make.unique(colnames(expression.df))
sum(duplicated(colnames(expression.df)))  # Should be 0

# --- Filter lowly expressed genes ---
# Genes with mean log2 expression below 1 are essentially silent
gene.means <- colMeans(expression.df[, colnames(expression.df) != "class"])
summary(gene.means)

keep.genes <- gene.means > 1
sum(keep.genes)  # Number of genes retained

expression.df.filtered <- expression.df[, c(keep.genes, TRUE)]
dim(expression.df.filtered)
