# DEDML

`DEDML` is an R package for donor-blocked double machine learning on single-cell data.

## Install options for users

From GitHub (recommended once pushed):

```r
remotes::install_github("YuxinYin0906/DEDML")
```

From a release/source tarball:

```r
install.packages("DEDML_0.1.4.tar.gz", repos = NULL, type = "source")
```

## Main function

Use `dedml_fit()` with simplified input names:

- `donor_id`
- `sample_id`
- `treatment`
- `cell_type`
- `cell_types` (which cell types to test)
- `donor_model` and `outcome_model`
- `outcome_distribution` (`"poisson"` by default; `"gaussian"` and `"nb"` are available)
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
- `dedml_summarize_design()` reports donor/cell/fold balance and stops early
  when fatal overlap problems are present, such as batch/site/confounder levels
  observed in only treated or only control donors.

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

design_report <- dedml_summarize_design(
  meta = meta_ready,
  donor_id = "HATIMID",
  sample_id = "sample_id",
  treatment = "D",
  cell_type = "CellType",
  confounder_spec = conf_spec,
  n_folds = 3,
  stop_on_error = FALSE
)

design_report$diagnostics

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
  outcome_distribution = "poisson",
  n_folds = 3,
  n_cores = 8,
  parallel_backend = "auto"
)

head(fit$results)
```

DEDML runs preflight diagnostics before fitting and stores the full diagnostic
table in `fit$diagnostics`. Warnings are emitted for designs with weak
treated-control overlap, very small donor counts in cross-fitting folds,
LightGBM treatment settings that cannot split donor-level data, numeric
covariate non-overlap, and degenerate propensity estimates near 0/1 or nearly
constant across donors. DEDML stops before model fitting when fatal design
diagnostics are found. In particular, categorical batch/site/confounder levels
seen only in treated donors or only in controls are treated as positivity
violations and must be combined, filtered to overlapping levels, or removed by
changing the cohort definition.

`fit$results$estimate` is the DEDML effect on the residual scale used for
testing. The output also includes `estimate_logfc`, a practical post-hoc log
fold-change effect size targeting the untreated-baseline mean. The conversion
uses the out-of-fold nuisance mean and treatment propensity: `"gaussian"` and
`"poisson"` use the Poisson residual scale, while `"nb"` uses the
negative-binomial residual scale. The p-values, adjusted p-values, and test
statistics remain based on the residual-scale DEDML regression.
