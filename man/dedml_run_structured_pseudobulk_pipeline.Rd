\name{dedml_run_structured_pseudobulk_pipeline}
\alias{dedml_run_structured_pseudobulk_pipeline}
\title{Run structured sum-pseudobulk DEDML}
\usage{
dedml_run_structured_pseudobulk_pipeline(
  counts,
  meta,
  donor_id = "HATIMID",
  sample_id = "sample_id",
  treatment = "D",
  cell_type = "CellType",
  cell_types = NULL,
  families = c("nb", "poisson"),
  treatment_confounders = c("Age", "Sex"),
  treatment_cell_summaries = c("log_lib_size", "log_n_cells", "nFeature_RNA", "percent.mt"),
  outcome_confounders = c("Age", "Sex", "log_n_cells", "nFeature_RNA", "percent.mt"),
  outdir = NULL,
  prefix = "dedml_structured_sum_pseudobulk",
  tuning_n_genes = "ALL",
  tune_only = FALSE,
  tune_then_run = TRUE,
  run_full_without_tuning = FALSE,
  n_cores = 1L,
  parallel_backend = c("auto", "fork", "psock"),
  min_samples_per_celltype = 5L,
  pvalue_calibration = c("none", "tail"),
  save_fit = TRUE,
  verbose = TRUE
)
}
\arguments{
\item{counts}{Gene-by-cell count matrix.}
\item{meta}{Cell metadata data.frame.}
\item{donor_id}{Donor id column name or vector.}
\item{sample_id}{Sample id column name or vector.}
\item{treatment}{Binary treatment column name or vector.}
\item{cell_type}{Cell type column name or vector.}
\item{cell_types}{Optional cell types to keep.}
\item{families}{Outcome families to run.}
\item{treatment_confounders}{Donor-level treatment confounders.}
\item{treatment_cell_summaries}{Pseudobulk summaries for the treatment model.}
\item{outcome_confounders}{Outcome covariates before offset and alpha terms.}
\item{outdir}{Optional output directory.}
\item{prefix}{Output file prefix.}
\item{tuning_n_genes}{\code{ALL} or number of genes for tuning.}
\item{tune_only}{Run tuning only.}
\item{tune_then_run}{Run selected full model after tuning.}
\item{run_full_without_tuning}{Run the smallest grid directly.}
\item{n_cores}{Number of gene-parallel workers.}
\item{parallel_backend}{Gene-level backend.}
\item{min_samples_per_celltype}{Minimum pseudobulk samples per cell type.}
\item{pvalue_calibration}{P-value calibration mode.}
\item{save_fit}{Save fit RDS files.}
\item{verbose}{Print progress.}
}
\description{
Builds sum pseudobulk counts by sample and cell type, adds the log library-size
offset and cell-type alpha indicators, tunes the structured LightGBM grid, and
runs the selected DEDML model.
}
\value{
A list with pseudobulk metadata, design report, tuning grid, selected grid, and metrics.
}
