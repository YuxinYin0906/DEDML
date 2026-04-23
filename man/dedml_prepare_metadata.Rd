\name{dedml_prepare_metadata}
\alias{dedml_prepare_metadata}
\title{Prepare metadata for DEDML input}
\usage{
dedml_prepare_metadata(
  meta,
  donor_id = "HATIMID",
  sample_id = "sample_id",
  treatment = "D",
  cell_type = "CellType",
  donor_confounders = NULL,
  cell_confounders = NULL,
  drop_missing = TRUE
)
}
\arguments{
\item{meta}{Metadata data.frame.}
\item{donor_id}{Donor id column name or vector.}
\item{sample_id}{Sample id column name or vector. If \code{NULL}, donor id is used.}
\item{treatment}{Treatment column name or vector (binary 0/1).}
\item{cell_type}{Cell type column name or vector.}
\item{donor_confounders}{Optional donor-level confounders to validate.}
\item{cell_confounders}{Optional cell-level confounders to validate.}
\item{drop_missing}{Whether to remove incomplete rows in required fields/confounders.}
}
\value{
Metadata data.frame with standardized columns \code{.donor_id}, \code{.sample_id}, \code{.D}, and \code{.CellType}.
}
\description{
Standardizes and validates metadata fields used by \code{dedml_fit}.
}
