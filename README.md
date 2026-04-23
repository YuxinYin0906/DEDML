# DEDML

`DEDML` is an R package for donor-blocked double machine learning on single-cell data.

It refactors the real-data script workflow into reusable functions so users can run analyses directly and choose confounders based on their own dataset.

## Install

From source tarball:

```r
install.packages("DEDML_0.1.0.tar.gz", repos = NULL, type = "source")
```

## Main workflow

1. Build a confounder spec (donor-level + cell-level).
2. Run `dedml_fit()` (or `dedml_fit_seurat()`).
3. Read `fit$results`.

## Main function inputs (`dedml_fit`)

- `counts`: gene-by-cell count matrix (rows = genes, cols = cells)
- `meta`: metadata data.frame for cells
- `donor_id_col`: donor ID column in `meta`
- `sample_id_col`: sample ID column in `meta`
- `treatment_col`: binary treatment column (0/1) in `meta`
- `cell_type_col`: cell type column in `meta`
- `confounder_spec`: object from `dedml_make_confounder_spec()`
- `n_cores`: number of CPU cores used for gene-level parallelization
- `treatment_learner`: default `"glm"`
- `outcome_learner`: default `"lightgbm"`

## Confounder helper functions

- `dedml_make_confounder_spec()` creates a validated confounder specification object
- `dedml_validate_confounder_spec()` checks that required confounder columns exist in `meta`

## Example

```r
library(DEDML)

conf_spec <- dedml_make_confounder_spec(
  donor_confounders = c("Age", "Sex"),
  cell_confounders = c("log_nCount_RNA", "nFeature_RNA", "percent.mt")
)

fit <- dedml_fit(
  counts = counts_mat,
  meta = meta_df,
  donor_id_col = "HATIMID",
  sample_id_col = "sample_id",
  treatment_col = "D",
  cell_type_col = "CellType",
  confounder_spec = conf_spec,
  n_folds = 3,
  n_cores = 8,
  treatment_learner = "glm",
  outcome_learner = "lightgbm"
)

head(fit$results)
```

## Defaults

- treatment nuisance learner: `glm`
- outcome nuisance learner: `lightgbm`
- supported treatment learners: `glm`, `lightgbm`
- supported outcome learners: `lightgbm`, `glm`
