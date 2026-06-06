# =============================================================================
# 05_correlation_analysis.R
# MAPK/ERK Pathway Analysis in Colorectal Cancer
# Step 5: Calculate Pearson correlation between each gene feature and
#         mean MAPK expression to rank candidate drug targets
# =============================================================================
# Requires: X_reg_500 from 04_regression_prep.R
#           (or load mapk_regression_500.csv)

# --- Load regression dataset if not in memory ---
if (!exists("X_reg_500")) {
  X_reg_500 <- read.csv("mapk_regression_500.csv", check.names = FALSE)
  cat("Loaded regression dataset:", nrow(X_reg_500), "rows x", ncol(X_reg_500), "columns\n")
}

# --- Calculate Pearson correlation ---
# Correlates each of the 500 gene features with the mean_mapk target
# Produces a continuous association score per gene: -1 to +1
cat("Calculating gene correlations with mean MAPK expression...\n")
gene_correlations <- cor(X_reg_500[, colnames(X_reg_500) != "mean_mapk"],
                         X_reg_500$mean_mapk)
cat("Done!\n")

# --- Convert to data frame and sort by absolute correlation ---
gene_cor_df <- data.frame(
  gene = rownames(gene_correlations),
  correlation = gene_correlations[, 1],
  stringsAsFactors = FALSE
)
gene_cor_df <- gene_cor_df[order(abs(gene_cor_df$correlation), decreasing = TRUE), ]

# --- Display top 20 candidates ---
cat("Top 20 MAPK-associated genes:\n")
print(head(gene_cor_df, 20))

# --- Summary statistics ---
cat("\nCorrelation summary:\n")
print(summary(gene_cor_df$correlation))
cat("\nGenes with correlation > 0.5 (strong candidates):", sum(gene_cor_df$correlation > 0.5), "\n")
cat("Genes with correlation > 0.3 (moderate-strong):", sum(gene_cor_df$correlation > 0.3), "\n")
cat("Genes with negative correlation:", sum(gene_cor_df$correlation < 0), "\n")

# --- Export ranked candidate lists ---
# All 500 genes ranked
write.csv(gene_cor_df, "mapk_candidate_genes_ranked.csv", row.names = FALSE)

# Strong positive candidates (correlation > 0.5)
strong_candidates <- gene_cor_df[gene_cor_df$correlation > 0.5, ]
write.csv(strong_candidates, "mapk_strong_candidates.csv", row.names = FALSE)

# Negatively correlated genes (potential tumour suppressors)
negative_candidates <- gene_cor_df[gene_cor_df$correlation < 0, ]
write.csv(negative_candidates, "mapk_negative_candidates.csv", row.names = FALSE)

cat("\nSaved!\n")
cat("All genes ranked:", nrow(gene_cor_df), "\n")
cat("Strong candidates (r > 0.5):", nrow(strong_candidates), "\n")
cat("Negative candidates:", nrow(negative_candidates), "\n")
