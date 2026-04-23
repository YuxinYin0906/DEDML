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
  treatment_params = list(),
  outcome_params = list(),
  min_cells_per_sample_ct = 10L,
  min_samples_per_celltype = 3L,
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
\item{treatment_params}{Optional donor model hyperparameters.}
\item{outcome_params}{Optional outcome model hyperparameters.}
\item{min_cells_per_sample_ct}{Minimum cells per sample-celltype pseudobulk unit.}
\item{min_samples_per_celltype}{Minimum pseudobulk samples required per cell type.}
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
A list with \code{results}, \code{donor_meta}, \code{settings}, and \code{errors}.
}
\description{
Run donor-blocked cross-fitted double machine learning with configurable
confounders, selectable cell types, model choices, and parallel cores.
}
