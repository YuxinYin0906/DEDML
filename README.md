# DEDML

`DEDML` is an R package for donor-blocked double machine learning on single-cell data.

## Install options for users

From GitHub (recommended once pushed):

```r
remotes::install_github("YuxinYin0906/DEDML")
```

From a release/source tarball:

```r
install.packages("DEDML_0.1.0.tar.gz", repos = NULL, type = "source")
```

## Main function

Use `dedml_fit()` with simplified input names:

- `donor_id`
- `sample_id`
- `treatment`
- `cell_type`
- `cell_types` (which cell types to test)
- `donor_model` and `outcome_model`
- `n_cores`
- `parallel_backend` (`"auto"`, `"fork"`, or `"psock"`)

Parallel gene-level fitting defaults to `parallel_backend = "auto"`. In this
mode, DEDML uses PSOCK workers when `n_cores > 1`, which is portable across
Unix-like systems, Windows, and LightGBM-backed runs. `parallel_backend =
"fork"` is still available as an explicit opt-in on platforms that support it.
If local socket worker creation is restricted by the environment, DEDML falls
back to serial execution instead of failing the analysis.

## Helper functions for confounders and metadata

- `dedml_prepare_metadata()` standardizes IDs/treatment/cell type fields and checks missingness.
- `dedml_make_confounder_spec()` builds the confounder input object.
- `dedml_validate_confounder_spec()` validates confounder columns in metadata.

## Example

```r
library(DEDML)

meta_ready <- dedml_prepare_metadata(
  meta = meta_df,
  donor_id = "HATIMID",
  sample_id = "sample_id",
  treatment = "D",
  cell_type = "CellType",
  donor_confounders = c("Age", "Sex"),
  cell_confounders = c("log_nCount_RNA", "nFeature_RNA", "percent.mt")
)

conf_spec <- dedml_make_confounder_spec(
  donor_confounders = c("Age", "Sex"),
  cell_confounders = c("log_nCount_RNA", "nFeature_RNA", "percent.mt")
)

fit <- dedml_fit(
  counts = counts_mat,
  meta = meta_ready,
  donor_id = "HATIMID",
  sample_id = "sample_id",
  treatment = "D",
  cell_type = "CellType",
  confounder_spec = conf_spec,
  cell_types = c("B cell Ig+", "CD4-EM like", "NK"),
  donor_model = "glm",
  outcome_model = "lightgbm",
  n_folds = 3,
  n_cores = 8,
  parallel_backend = "auto"
)

head(fit$results)
```
