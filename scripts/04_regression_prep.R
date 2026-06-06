# =============================================================================
# 04_regression_prep.R
# MAPK/ERK Pathway Analysis in Colorectal Cancer
# Step 4: Prepare regression dataset
#         Target variable = mean MAPK pathway expression per patient
#         Features = top 500 genes by variance
# =============================================================================
# Requires: mapk_weka_input.csv from 03_classification_prep.R

# --- Load the classification dataset ---
expr <- read.csv("mapk_weka_input.csv", row.names = 1, check.names = FALSE)

# --- Remove class column and transpose ---
# Now patients are rows, genes are columns
X_reg <- expr[, colnames(expr) != "class"]
X_reg <- as.data.frame(t(X_reg))

# --- Construct target variable: mean MAPK expression per patient ---
# Identify the 82 MAPK pathway genes
mapk_gene_list <- rownames(expr)[expr[, "class"] == "MAPK"]
cat("Number of MAPK genes:", length(mapk_gene_list), "\n")

# Extract MAPK gene expression and compute column means
# This produces a single continuous value per patient (range ~8-12 on log2 scale)
mapk_expr <- expr[mapk_gene_list, colnames(expr) != "class"]
mean_mapk_per_patient <- colMeans(mapk_expr)
cat("Calculated mean MAPK expression across", length(mean_mapk_per_patient), "patients\n")

# Add as target column
X_reg$mean_mapk <- mean_mapk_per_patient

# --- Variance-based feature selection: top 500 genes ---
# Genes with highest variance across patients are most informative
gene_cols <- colnames(X_reg)[colnames(X_reg) != "mean_mapk"]
feature_vars <- apply(X_reg[, gene_cols], 2, var)

top500_genes <- names(sort(feature_vars, decreasing = TRUE)[1:500])
X_reg_500 <- X_reg[, c(top500_genes, "mean_mapk")]

# --- Export to CSV for WEKA ---
write.csv(X_reg_500, "mapk_regression_500.csv", row.names = FALSE)
cat("Done!", nrow(X_reg_500), "rows x", ncol(X_reg_500), "columns\n")
# Expected: 698 patients x 501 columns (500 gene features + mean_mapk target)
