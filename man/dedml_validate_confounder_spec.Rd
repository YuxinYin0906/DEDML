\name{dedml_validate_confounder_spec}
\alias{dedml_validate_confounder_spec}
\title{Validate a DEDML confounder specification}
\usage{
dedml_validate_confounder_spec(meta, confounder_spec)
}
\arguments{
\item{meta}{Metadata data.frame.}
\item{confounder_spec}{List returned by \code{dedml_make_confounder_spec()}.}
}
\value{
Returns \code{TRUE} invisibly. Errors if required columns are missing.
}
\description{
Checks whether a confounder specification is valid for a given metadata table.
}
