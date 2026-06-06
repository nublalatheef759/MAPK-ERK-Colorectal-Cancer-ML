# WEKA Configuration Documentation
### MAPK/ERK Pathway Analysis in Colorectal Cancer

WEKA is a GUI-based machine learning platform. Since there is no exportable code, 
this document records the exact configurations used to reproduce the analysis.

**WEKA Version:** 3.9.6  
**Memory:** 4 GB (set via `maxheap=4096m` in RunWeka.ini)

---

## 1. Data Loading

- **Input file (Classification):** `mapk_weka_input.csv` (27,384 genes × 699 columns)
- **Input file (Regression):** `mapk_regression_500.csv` (698 patients × 501 columns)
- **Interface:** Explorer → Preprocess tab → Open file

---

## 2. Classification Preprocessing

### Remove gene identifier column
- Filter: `filters → unsupervised → attribute → Remove`
- Attribute index: 1 (removes the row name column from CSV export)

### Feature selection (InfoGain, top 100)
- Filter: `filters → supervised → attribute → AttributeSelection`
- Evaluator: `InfoGainAttributeEval`
- Search: `Ranker` with `numToSelect = 100`
- Result: 101 attributes (100 patient features + 1 class)

### Class balancing (applied for balanced runs only)
- Filter: `filters → supervised → instance → ClassBalancer`
- Applied before classification to equalise MAPK/Other class weights

---

## 3. Classification Models

All evaluated using **10-fold stratified cross-validation** with class set to `(Nom) class`.

### Naive Bayes
- Classifier: `classifiers → bayes → NaiveBayes`
- Settings: all defaults
  - batchSize: 100
  - useKernelEstimator: False
  - useSupervisedDiscretization: False
- Run configurations:
  - With ClassBalancer
  - Without ClassBalancer

### Random Forest
- Classifier: `classifiers → trees → RandomForest`
- Settings: all defaults
  - numIterations: 100
  - maxDepth: 0 (unlimited)
  - numFeatures: 0 (auto = sqrt of total features)
- Run configurations:
  - With ClassBalancer
  - Without ClassBalancer

---

## 4. Classification Results

| Classifier | Balancing | Accuracy | MAPK Recall | AUC | Kappa |
|---|---|---|---|---|---|
| Random Forest | None | 99.70% | 0.000 | 0.681 | 0.000 |
| Random Forest | ClassBalancer | 50.61% | 0.012 | 0.746 | 0.012 |
| Naive Bayes | ClassBalancer | 73.80% | **0.866** | **0.776** | **0.476** |
| Naive Bayes | None | 62.13% | 0.866 | 0.767 | 0.008 |

---

## 5. Regression Model

### MultilayerPerceptron (ANN Regression)
- Classifier: `classifiers → functions → MultilayerPerceptron`
- Evaluation: 10-fold cross-validation
- Input: `mapk_regression_500.csv` (no additional preprocessing in WEKA)

### Model settings (constant across all runs):
- Learning rate (`-L`): 0.3
- Momentum (`-M`): 0.2
- Hidden layers (`-H`): a (auto-calculated)
- normalizeNumericClass: True
- normalizeAttributes: True
- Validation set size (`-V`): 0

### Training time configurations tested:

| Epochs | Correlation | MAE | RMSE | RAE |
|---|---|---|---|---|
| 100 | 0.7143 | 0.4239 | 0.6351 | 61.78% |
| **500** | **0.7755** | **0.3805** | **0.5775** | **55.47%** |
| 1000 | 0.7660 | 0.3795 | 0.5838 | 55.32% |

**Optimal configuration:** 500 epochs (best correlation; 1000 epochs shows mild overfitting)

---

## 6. Visualising Results

### Predicted vs Actual scatter plot
- Right-click result in Result list → Visualize classifier errors
- X axis: mean_mapk (Num)
- Y axis: predictedmean_mapk (Num)

### ROC curve
- Right-click result in Result list → Visualize threshold curve → select MAPK class
