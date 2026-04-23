utils::globalVariables(c(
  ".donor_id", ".D", "CellType", "n_cells", "p_value", "estimate", "statistic"
))

.onAttach <- function(libname, pkgname) {
  packageStartupMessage("DEDML loaded. Use dedml_fit() or dedml_fit_seurat() to run the pipeline.")
}
