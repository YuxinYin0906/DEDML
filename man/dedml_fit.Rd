\name{dedml_fit}
\alias{dedml_fit}
\title{Fit DEDML on Cell-Level Data}
\usage{
dedml_fit(
  counts,
  meta,
  donor_id_col = "HATIMID",
  sample_id_col = "sample_id",
  treatment_col = "D",
  cell_type_col = "CellType",
  treatment_confounders = c("Age", "Sex"),
  treatment_cell_summaries = c("log_nCount_RNA", "nFeature_RNA", "percent.mt"),
  outcome_confounders = c("Age", "Sex", "log_nCount_RNA", "nFeature_RNA", "percent.mt"),
  focus_celltypes = NULL,
  gene_subset = NULL,
  n_folds = 3L,
  n_cores = 1L,
  treatment_learner = "glm",
  outcome_learner = "lightgbm",
  treatment_params = list(),
  outcome_params = list(),
  min_cells_per_sample_ct = 10L,
  min_samples_per_celltype = 3L,
  seed = 123L,
  verbose = TRUE
)
}
\arguments{
\item{counts}{Gene-by-cell count matrix.}
\item{meta}{Cell metadata data.frame.}
\item{donor_id_col}{Donor id column in metadata.}
\item{sample_id_col}{Sample id column in metadata.}
\item{treatment_col}{Binary treatment column (0/1).}
\item{cell_type_col}{Cell type column.}
\item{treatment_confounders}{Donor-level confounders for treatment model.}
\item{treatment_cell_summaries}{Cell-level covariates summarized by donor for treatment model.}
\item{outcome_confounders}{Cell-level confounders for outcome model.}
\item{focus_celltypes}{Optional cell types to keep.}
\item{gene_subset}{Optional gene subset to keep.}
\item{n_folds}{Number of donor folds.}
\item{n_cores}{Number of cores.}
\item{treatment_learner}{Treatment nuisance learner.}
\item{outcome_learner}{Outcome nuisance learner.}
\item{treatment_params}{Learner parameters for treatment model.}
\item{outcome_params}{Learner parameters for outcome model.}
\item{min_cells_per_sample_ct}{Minimum cells in sample-celltype pseudobulk.}
\item{min_samples_per_celltype}{Minimum samples required per cell type.}
\item{seed}{Random seed.}
\item{verbose}{Print progress.}
}
\value{
A list with \code{results}, \code{donor_meta}, \code{settings}, and \code{errors}.
}
\description{
Run donor-blocked cross-fitted double machine learning with configurable
confounders and learners.
}
