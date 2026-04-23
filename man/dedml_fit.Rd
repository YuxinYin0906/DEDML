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
  confounder_spec = NULL,
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
\item{counts}{Gene-by-cell count matrix. Rows are genes and columns are cells.}
\item{meta}{Cell metadata data.frame. Row order should match columns of \code{counts}, or row names should match column names of \code{counts}.}
\item{donor_id_col}{Column in \code{meta} containing donor IDs.}
\item{sample_id_col}{Column in \code{meta} containing sample IDs.}
\item{treatment_col}{Binary treatment column (0/1) in \code{meta}; must be constant within donor.}
\item{cell_type_col}{Column in \code{meta} containing cell type labels.}
\item{confounder_spec}{Optional list from \code{dedml_make_confounder_spec()}. If supplied, it overrides \code{treatment_confounders}, \code{treatment_cell_summaries}, and \code{outcome_confounders}.}
\item{treatment_confounders}{Donor-level confounders for treatment nuisance model.}
\item{treatment_cell_summaries}{Cell-level confounders to summarize by donor for treatment nuisance model.}
\item{outcome_confounders}{Cell-level confounders for outcome nuisance model (can include donor-level columns too).}
\item{focus_celltypes}{Optional vector of cell types to keep.}
\item{gene_subset}{Optional vector of genes to run.}
\item{n_folds}{Number of donor-level cross-fitting folds.}
\item{n_cores}{Number of CPU cores for gene-level parallelization.}
\item{treatment_learner}{Treatment nuisance learner (\code{"glm"} or \code{"lightgbm"}).}
\item{outcome_learner}{Outcome nuisance learner (\code{"lightgbm"} or \code{"glm"}).}
\item{treatment_params}{Optional hyperparameters for treatment learner.}
\item{outcome_params}{Optional hyperparameters for outcome learner.}
\item{min_cells_per_sample_ct}{Minimum number of cells per sample-celltype pseudobulk unit.}
\item{min_samples_per_celltype}{Minimum number of pseudobulk samples required for a cell type.}
\item{seed}{Random seed.}
\item{verbose}{Whether to print progress messages.}
}
\value{
A list with \code{results}, \code{donor_meta}, \code{settings}, and \code{errors}.
}
\description{
Run donor-blocked cross-fitted double machine learning with configurable
confounders and learners.
}
