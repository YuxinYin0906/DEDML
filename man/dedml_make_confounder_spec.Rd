\name{dedml_make_confounder_spec}
\alias{dedml_make_confounder_spec}
\title{Build a confounder specification for DEDML}
\usage{
dedml_make_confounder_spec(
  donor_confounders,
  cell_confounders,
  treatment_cell_summaries = cell_confounders
)
}
\arguments{
\item{donor_confounders}{Donor-level confounder column names.}
\item{cell_confounders}{Cell-level confounder column names.}
\item{treatment_cell_summaries}{Cell-level confounders to summarize by donor in treatment nuisance model.}
}
\value{
A list with \code{treatment_confounders}, \code{treatment_cell_summaries}, and \code{outcome_confounders}.
}
\description{
Creates a reusable confounder object that can be passed to \code{dedml_fit(confounder_spec = ...)}.
}
