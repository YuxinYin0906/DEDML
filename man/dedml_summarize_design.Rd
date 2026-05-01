\name{dedml_summarize_design}
\alias{dedml_summarize_design}
\title{Summarize and validate a DEDML design}
\usage{
dedml_summarize_design(
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
  n_folds = 3L,
  donor_model = c("glm", "lightgbm"),
  treatment_params = list(),
  seed = 123L,
  drop_missing = TRUE,
  stop_on_error = TRUE
)
}
\arguments{
\item{meta}{Metadata data.frame.}
\item{donor_id}{Donor id column name or vector.}
\item{sample_id}{Sample id column name or vector. If \code{NULL}, donor id is used.}
\item{treatment}{Treatment column name or vector (binary 0/1).}
\item{cell_type}{Cell type column name or vector.}
\item{confounder_spec}{Optional list from \code{dedml_make_confounder_spec()}.}
\item{treatment_confounders}{Donor-level confounders for the treatment model.}
\item{treatment_cell_summaries}{Cell-level covariates summarized by donor for the treatment model.}
\item{outcome_confounders}{Covariates for the outcome model.}
\item{cell_types}{Optional cell types to keep before summarizing.}
\item{n_folds}{Number of donor-blocked cross-fitting folds.}
\item{donor_model}{Treatment model learner, either \code{"glm"} or \code{"lightgbm"}.}
\item{treatment_params}{Optional treatment learner hyperparameters.}
\item{seed}{Random seed for fold construction.}
\item{drop_missing}{Whether to drop incomplete rows for selected confounders, matching \code{dedml_fit()} behavior.}
\item{stop_on_error}{Whether to stop when fatal design diagnostics are found.}
}
\value{
A \code{dedml_design_summary} list with \code{summary}, \code{treatment},
\code{folds}, \code{cell_types}, \code{categorical}, \code{numeric}, and
\code{diagnostics} tables.
}
\description{
Checks donor balance, cross-fitting folds, cell-type representation, and
treatment overlap for categorical and numeric confounders before running the
gene-level DEDML fit.
}
\details{
Categorical confounder levels observed in only one treatment group are treated
as fatal errors because the treatment effect is not identifiable for those
levels. This catches common failures such as disease status being completely
separated by batch, site, chemistry, or another design variable.
}
