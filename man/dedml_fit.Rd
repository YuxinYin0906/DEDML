\name{dedml_fit}
\alias{dedml_fit}
\title{Fit DEDML on Cell-Level Data}
\usage{
dedml_fit(
  counts,
  meta,
  donor_id = "HATIMID",
  sample_id = "sample_id",
  treatment = "D",
  cell_type = "CellType",
  confounder_spec = NULL,
  treatment_confounders = c("Age", "Sex"),
  treatment_cell_summaries = c("log_nCount_RNA", "nFeature_RNA", "percent.mt"),
  outcome_confounders = c("Age", "Sex", "log_nCount_RNA", "nFeature_RNA", "percent.mt"),
  cell_types = NULL,
  gene_subset = NULL,
  n_folds = 3L,
  n_cores = 1L,
  donor_model = "glm",
  outcome_model = "lightgbm",
  outcome_distribution = c("poisson", "gaussian", "nb"),
  treatment_params = list(),
  outcome_params = list(),
  min_cells_per_sample_ct = 10L,
  min_samples_per_celltype = 3L,
  nb_size_min = 0.1,
  nb_size_max = 1e+06,
  pvalue_calibration = c("none", "tail"),
  calibration_quantile = 0.9,
  parallel_backend = c("auto", "fork", "psock"),
  gene_chunk_size = NULL,
  seed = 123L,
  verbose = TRUE,
  donor_id_col = NULL,
  sample_id_col = NULL,
  treatment_col = NULL,
  cell_type_col = NULL,
  focus_celltypes = NULL,
  treatment_learner = NULL,
  outcome_learner = NULL
)
}
\arguments{
\item{counts}{Gene-by-cell count matrix. Rows are genes, columns are cells.}
\item{meta}{Cell metadata data.frame.}
\item{donor_id}{Donor id column name in \code{meta} (or vector of length \code{nrow(meta)}).}
\item{sample_id}{Sample id column name in \code{meta} (or vector). If \code{NULL}, donor id is used.}
\item{treatment}{Binary treatment column name in \code{meta} (0/1), or vector.}
\item{cell_type}{Cell type column name in \code{meta}, or vector.}
\item{confounder_spec}{Optional list from \code{dedml_make_confounder_spec()}; overrides confounder vectors below.}
\item{treatment_confounders}{Donor-level confounders for the donor treatment model.}
\item{treatment_cell_summaries}{Cell-level covariates summarized by donor for the donor treatment model.}
\item{outcome_confounders}{Confounders for the outcome model.}
\item{cell_types}{Optional vector of cell types to test.}
\item{gene_subset}{Optional vector of genes to test.}
\item{n_folds}{Number of donor-level cross-fitting folds.}
\item{n_cores}{Number of CPU cores for gene-level parallelization.}
\item{donor_model}{Donor treatment model: \code{"glm"} or \code{"lightgbm"}.}
\item{outcome_model}{Outcome model: \code{"lightgbm"} or \code{"glm"}.}
\item{outcome_distribution}{Outcome likelihood/working distribution:
\code{"poisson"} by default, with \code{"gaussian"} and \code{"nb"} also
available. For \code{outcome_model = "glm"}, Poisson and Gaussian use \code{stats::glm.fit()}
and NB uses \code{mgcv::gam()} with \code{mgcv::nb()}. For
\code{outcome_model = "lightgbm"}, these map to LightGBM L2 regression,
Poisson, and free-dispersion NB objectives. Gaussian/L2 changes only the
nuisance mean learner; its final residual conversion still uses Poisson
variance.}
\item{treatment_params}{Optional donor model hyperparameters.}
\item{outcome_params}{Optional outcome learner hyperparameters. The outcome
objective/family is controlled by \code{outcome_distribution}, so redundant
\code{objective}, \code{metric}, and \code{family} entries are ignored. For
\code{outcome_model = "glm"} and \code{outcome_distribution = "poisson"},
\code{batched_irls = TRUE} enables an exact shared-design IRLS implementation;
it uses \code{n_cores} and \code{parallel_backend} to split genes into
independent blocks. It is opt-in because the regular parallel scalar
\code{glm.fit()} path can be faster for some backend/data-size combinations.}
\item{min_cells_per_sample_ct}{Minimum cells per sample-celltype pseudobulk unit.}
\item{min_samples_per_celltype}{Minimum pseudobulk samples required per cell type.}
\item{nb_size_min}{Minimum NB size estimate used for residual conversion.}
\item{nb_size_max}{Maximum NB size estimate used for residual conversion; large
values approximate Poisson residuals.}
\item{pvalue_calibration}{Optional p-value calibration. \code{"tail"} estimates
cell-type-specific upper-tail statistic inflation before BH adjustment.}
\item{calibration_quantile}{Chi-square quantile used when
\code{pvalue_calibration = "tail"}.}
\item{parallel_backend}{Gene-level parallel backend. With \code{n_cores > 1},
\code{"auto"} uses PSOCK workers for portability; \code{"fork"} remains
available as an explicit opt-in on platforms that support it.}
\item{gene_chunk_size}{Number of genes per parallel task. If \code{NULL},
DEDML chooses a backend-aware value. This also controls the maximum number of
genes materialized together for batched residual regression.}
\item{seed}{Random seed.}
\item{verbose}{Whether to print progress messages.}
\item{donor_id_col}{Deprecated alias for \code{donor_id}.}
\item{sample_id_col}{Deprecated alias for \code{sample_id}.}
\item{treatment_col}{Deprecated alias for \code{treatment}.}
\item{cell_type_col}{Deprecated alias for \code{cell_type}.}
\item{focus_celltypes}{Deprecated alias for \code{cell_types}.}
\item{treatment_learner}{Deprecated alias for \code{donor_model}.}
\item{outcome_learner}{Deprecated alias for \code{outcome_model}.}
}
\value{
A list with \code{results}, \code{donor_meta}, \code{settings},
\code{diagnostics}, and \code{errors}.
The \code{results} table reports the residual-scale DEDML effect in
\code{estimate}. It also includes \code{estimate_logfc}, a practical post-hoc
log fold-change effect size targeting the untreated-baseline mean. The
conversion uses the out-of-fold nuisance mean and treatment propensity:
\code{"gaussian"} and \code{"poisson"} use the Poisson residual scale, while
\code{"nb"} uses the negative-binomial residual scale. DEDML inference uses the
residual-scale statistic and p-value. The \code{diagnostics} table records
design and treatment model issues such as weak fold overlap, numeric
non-overlap, and degenerate propensity estimates. Fatal design diagnostics,
including separated categorical batch/site/confounder levels observed in only
one treatment group, stop the fit before model training.
}
\description{
Run donor-blocked cross-fitted double machine learning with configurable
confounders, selectable cell types, model choices, and parallel cores.
}
