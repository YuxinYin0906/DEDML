\name{fit_residual_regression_pseudobulk_matrix}
\alias{fit_residual_regression_pseudobulk_matrix}
\title{Batched pseudobulk residual regression}
\usage{
fit_residual_regression_pseudobulk_matrix(
  y_resid_mat,
  d_resid,
  cell_type,
  sample_id,
  donor_id,
  gene_names = NULL,
  min_cells_per_sample_ct = 10L,
  min_samples = 3L,
  cluster_robust = TRUE,
  fallback_scalar = TRUE
)
}
\arguments{
\item{y_resid_mat}{Gene-by-cell matrix of outcome residuals.}
\item{d_resid}{Treatment residuals at cell level.}
\item{cell_type}{Cell type labels.}
\item{sample_id}{Sample identifiers.}
\item{donor_id}{Donor identifiers.}
\item{gene_names}{Optional gene names for rows of \code{y_resid_mat}.}
\item{min_cells_per_sample_ct}{Minimum cell count per sample-celltype unit.}
\item{min_samples}{Minimum number of pseudobulk samples required per cell type.}
\item{cluster_robust}{Use donor-cluster robust standard errors when possible.}
\item{fallback_scalar}{Use the scalar implementation for genes with non-finite
pseudobulk residuals or rank-deficient batch design.}
}
\description{
Fits the same weighted pseudobulk residual regression as
\code{fit_residual_regression_pseudobulk()}, but solves all genes jointly when
the design matrix is shared across genes.
}
\value{
A data.frame with effect estimates, standard errors, and p-values per
gene-celltype pair.
}
