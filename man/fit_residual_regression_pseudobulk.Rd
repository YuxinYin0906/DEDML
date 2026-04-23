\name{fit_residual_regression_pseudobulk}
\alias{fit_residual_regression_pseudobulk}
\title{Pseudobulk residual regression}
\usage{
fit_residual_regression_pseudobulk(
  y_resid,
  d_resid,
  cell_type,
  sample_id,
  donor_id,
  min_cells_per_sample_ct = 10L,
  min_samples = 3L,
  cluster_robust = TRUE
)
}
\arguments{
\item{y_resid}{Outcome residuals at cell level.}
\item{d_resid}{Treatment residuals at cell level.}
\item{cell_type}{Cell type labels.}
\item{sample_id}{Sample identifiers.}
\item{donor_id}{Donor identifiers.}
\item{min_cells_per_sample_ct}{Minimum cell count per sample-celltype unit.}
\item{min_samples}{Minimum number of pseudobulk samples required per cell type.}
\item{cluster_robust}{Use donor-cluster robust standard errors when possible.}
}
\description{
Fits weighted pseudobulk residual regression with treatment-by-celltype
interactions and returns celltype-specific estimates.
}
\value{
A data.frame with effect estimates, standard errors, and p-values per cell type.
}
