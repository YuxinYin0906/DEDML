\name{make_stratified_donor_folds}
\alias{make_stratified_donor_folds}
\title{Create stratified donor folds}
\usage{
make_stratified_donor_folds(donor_id, treatment, n_folds = 3L, seed = 123L)
}
\arguments{
\item{donor_id}{Donor identifiers.}
\item{treatment}{Binary treatment vector at donor level (0/1).}
\item{n_folds}{Number of folds.}
\item{seed}{Random seed.}
}
\description{
Creates donor-level cross-fitting folds stratified by binary treatment.
}
\value{
A data.frame with \code{donor_id} and \code{fold}.
}
