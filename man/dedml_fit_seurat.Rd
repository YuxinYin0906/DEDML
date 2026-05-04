\name{dedml_fit_seurat}
\alias{dedml_fit_seurat}
\title{Fit DEDML from Seurat}
\usage{
dedml_fit_seurat(
  object,
  assay = "RNA",
  layer = "counts",
  donor_id = "HATIMID",
  sample_id = NULL,
  treatment = "D",
  cell_type = "celltype_manual",
  ...
)
}
\arguments{
\item{object}{Seurat object.}
\item{assay}{Assay name.}
\item{layer}{Layer name for counts.}
\item{donor_id}{Donor id metadata column.}
\item{sample_id}{Sample id metadata column; defaults to donor id.}
\item{treatment}{Binary treatment metadata column.}
\item{cell_type}{Cell type metadata column.}
\item{...}{Arguments passed to \code{dedml_fit}.}
}
\value{
A list with \code{results}, \code{gene_results}, \code{donor_meta},
\code{settings}, \code{diagnostics}, and \code{errors}.
}
\description{
Convenience wrapper that extracts counts and metadata from Seurat and runs
\code{dedml_fit}.
}
