# DEDML

`DEDML` provides a donor-blocked double machine learning pipeline for
single-cell differential expression.

## Why this package

This package refactors the real-data script workflow into reusable functions so users can:

- choose treatment confounders based on their own dataset
- choose outcome confounders based on their own dataset
- switch nuisance learners without rewriting scripts
- run the full pipeline directly from package functions

## Core API

- `dedml_fit()` for matrix + metadata input
- `dedml_fit_seurat()` for Seurat input

Default learners:

- treatment model: `glm`
- outcome model: `lightgbm`

Supported learners currently:

- treatment: `glm`, `lightgbm`
- outcome: `lightgbm`, `glm`

## Minimal example

```r
library(DEDML)

fit <- dedml_fit(
  counts = counts_mat,
  meta = meta_df,
  donor_id_col = "HATIMID",
  sample_id_col = "sample_id",
  treatment_col = "D",
  cell_type_col = "CellType",
  treatment_confounders = c("Age", "Sex"),
  treatment_cell_summaries = c("log_nCount_RNA", "nFeature_RNA", "percent.mt"),
  outcome_confounders = c("Age", "Sex", "log_nCount_RNA", "nFeature_RNA", "percent.mt"),
  n_folds = 3,
  n_cores = 4,
  treatment_learner = "glm",
  outcome_learner = "lightgbm"
)

head(fit$results)
```

## Notes

- For `lightgbm` learners, install the `lightgbm` R package first.
- Legacy script used for refactoring is kept at:
  `inst/scripts/legacy_20260415_run_dedml_pb.R`.
