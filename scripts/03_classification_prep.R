# =============================================================================
# 03_classification_prep.R
# MAPK/ERK Pathway Analysis in Colorectal Cancer
# Step 3: Reframe matrix for pathway membership classification
#         Genes become instances, patients become features
#         Each gene labelled as MAPK or Other
# =============================================================================
# Requires: expression.df.filtered, mapk.gene.list from 02_preprocessing.R

# --- Transpose: genes become rows, patient samples become features ---
# This reframes the problem from predicting sample phenotype
# to predicting gene pathway membership
gene.matrix <- as.data.frame(
  t(expression.df.filtered[, colnames(expression.df.filtered) != "class"])
)

# Verify dimensions
dim(gene.matrix)

# --- Add MAPK/Other class label ---
# Genes in the 84-gene Reactome Oncogenic MAPK Signalling list = "MAPK"
# All other genes = "Other"
gene.matrix$class <- ifelse(rownames(gene.matrix) %in% mapk.gene.list,
                            "MAPK", "Other")

# Check class distribution
table(gene.matrix$class)
# Expected: 82 MAPK, ~27302 Other
# Note: 2 of the original 84 MAPK genes were removed during low expression filtering

# --- Export to CSV for WEKA ---
write.csv(gene.matrix, "mapk_weka_input.csv", row.names = TRUE)
file.exists("mapk_weka_input.csv")

cat("Classification dataset exported successfully\n")
cat("Dimensions:", nrow(gene.matrix), "genes x", ncol(gene.matrix), "columns\n")
cat("(", ncol(gene.matrix) - 1, "patient features + 1 class label)\n")
