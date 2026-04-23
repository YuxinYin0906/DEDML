\name{dedml_fit_seurat}
\alias{dedml_fit_seurat}
\title{Fit DEDML from Seurat}
\usage{
dedml_fit_seurat(
  object,
  assay = "RNA",
  layer = "counts",
  donor_id_col = "HATIMID",
  sample_id_col = NULL,
  treatment_col = "D",
  cell_type_col = "celltype_manual",
  ...
)
}
\arguments{
\item{object}{Seurat object.}
\item{assay}{Assay name.}
\item{layer}{Layer name for counts.}
\item{donor_id_col}{Donor id column.}
\item{sample_id_col}{Sample id column; defaults to donor id.}
\item{treatment_col}{Binary treatment column.}
\item{cell_type_col}{Cell type column.}
\item{...}{Arguments passed to \code{dedml_fit}.}
}
\value{
A list with \code{results}, \code{donor_meta}, \code{settings}, and \code{errors}.
}
\description{
Convenience wrapper that extracts counts and metadata from Seurat and runs
\code{dedml_fit}.
}
