available_dedml_learners <- function() {
  list(
    donor = c("glm", "lightgbm"),
    outcome_model = c("cell", "pseudobulk"),
    outcome_learner = c("lightgbm", "glm"),
    outcome_distribution = c("poisson", "gaussian", "nb")
  )
}

if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    ".CellType", ".DEDML_POISSON_IRLS_STATE", ".DEDML_PSOCK_STATE",
    ".sample_id", ".unit_key", "donor_id", "fold", "level", "sample_id",
    "tail_lambda", "variable"
  ))
}

#' Build a Confounder Specification for DEDML
#'
#' Helper to define donor-level and cell-level confounders once, then pass
#' them into [dedml_fit()] through `confounder_spec`.
#'
#' @param donor_confounders Character vector of donor-level confounders.
#' @param cell_confounders Character vector of cell-level confounders.
#' @param treatment_cell_summaries Optional cell-level confounders summarized by donor
#' for treatment nuisance modeling. Defaults to `cell_confounders`.
#'
#' @return A named list with `treatment_confounders`, `treatment_cell_summaries`,
#' and `outcome_confounders`.
#' @export
#'
dedml_make_confounder_spec <- function(
    donor_confounders,
    cell_confounders,
    treatment_cell_summaries = cell_confounders) {
  donor_confounders <- unique(as.character(donor_confounders))
  cell_confounders <- unique(as.character(cell_confounders))
  treatment_cell_summaries <- unique(as.character(treatment_cell_summaries))

  donor_confounders <- donor_confounders[nzchar(donor_confounders)]
  cell_confounders <- cell_confounders[nzchar(cell_confounders)]
  treatment_cell_summaries <- treatment_cell_summaries[nzchar(treatment_cell_summaries)]

  if (length(donor_confounders) == 0L) {
    stop("donor_confounders must contain at least one column name.", call. = FALSE)
  }

  list(
    treatment_confounders = donor_confounders,
    treatment_cell_summaries = treatment_cell_summaries,
    outcome_confounders = unique(c(donor_confounders, cell_confounders))
  )
}

#' Validate a DEDML Confounder Specification Against Metadata
#'
#' @param meta Metadata data.frame.
#' @param confounder_spec List returned by [dedml_make_confounder_spec()].
#'
#' @return Invisibly returns `TRUE`; throws an error if invalid.
#' @export
#'
dedml_validate_confounder_spec <- function(meta, confounder_spec) {
  required_names <- c("treatment_confounders", "treatment_cell_summaries", "outcome_confounders")
  if (!is.list(confounder_spec) || !all(required_names %in% names(confounder_spec))) {
    stop(
      "confounder_spec must be a list with: ",
      paste(required_names, collapse = ", "),
      call. = FALSE
    )
  }

  all_vars <- unique(c(
    confounder_spec$treatment_confounders,
    confounder_spec$treatment_cell_summaries,
    confounder_spec$outcome_confounders
  ))
  all_vars <- as.character(all_vars)
  all_vars <- all_vars[nzchar(all_vars)]

  missing_vars <- setdiff(all_vars, colnames(meta))
  if (length(missing_vars) > 0L) {
    stop(
      "Missing confounder columns in meta: ",
      paste(missing_vars, collapse = ", "),
      call. = FALSE
    )
  }

  invisible(TRUE)
}

.dedml_default_empty_diagnostics <- function() {
  data.frame(
    severity = character(),
    stage = character(),
    check = character(),
    message = character(),
    recommendation = character(),
    stringsAsFactors = FALSE
  )
}

.dedml_format_diagnostic_bullets <- function(diagnostics, max_items = 8L) {
  if (!is.data.frame(diagnostics) || nrow(diagnostics) == 0L) {
    return(character())
  }
  show_df <- utils::head(diagnostics, max_items)
  bullets <- paste0(
    "- [", show_df$stage, "/", show_df$check, "] ",
    show_df$message,
    " Recommendation: ", show_df$recommendation
  )
  if (nrow(diagnostics) > nrow(show_df)) {
    bullets <- c(
      bullets,
      paste0("- ... ", nrow(diagnostics) - nrow(show_df), " more issue(s).")
    )
  }
  bullets
}

.dedml_stop_for_fatal_diagnostics <- function(diagnostics) {
  if (!is.data.frame(diagnostics) || nrow(diagnostics) == 0L) {
    return(invisible(FALSE))
  }
  fatal_diagnostics <- diagnostics[diagnostics$severity == "error", , drop = FALSE]
  if (nrow(fatal_diagnostics) == 0L) {
    return(invisible(FALSE))
  }
  stop(
    paste0(
      "DEDML cannot continue because fatal design diagnostics were detected:\n",
      paste(.dedml_format_diagnostic_bullets(fatal_diagnostics), collapse = "\n"),
      "\nRun dedml_summarize_design(..., stop_on_error = FALSE) to inspect the full design report."
    ),
    call. = FALSE
  )
}

.dedml_prepare_design_metadata <- function(
    meta,
    donor_id,
    sample_id,
    treatment,
    cell_type,
    confounder_spec,
    treatment_confounders,
    treatment_cell_summaries,
    outcome_confounders,
    cell_types = NULL,
    drop_missing = TRUE) {
  meta <- as.data.frame(meta, stringsAsFactors = FALSE)
  if (is.null(sample_id)) {
    sample_id <- donor_id
  }

  donor_res <- .resolve_input_column(meta, donor_id, "donor_id")
  sample_res <- .resolve_input_column(meta, sample_id, "sample_id")
  treatment_res <- .resolve_input_column(meta, treatment, "treatment")
  celltype_res <- .resolve_input_column(meta, cell_type, "cell_type")

  meta$.donor_id <- as.character(donor_res$values)
  meta$.sample_id <- as.character(sample_res$values)
  meta$.D <- as.integer(treatment_res$values)
  meta$.CellType <- as.character(celltype_res$values)

  if (!all(meta$.D %in% c(0L, 1L, NA_integer_))) {
    stop("Treatment column must be coded as 0/1.", call. = FALSE)
  }

  if (!is.null(confounder_spec)) {
    dedml_validate_confounder_spec(meta = meta, confounder_spec = confounder_spec)
    treatment_confounders <- confounder_spec$treatment_confounders
    treatment_cell_summaries <- confounder_spec$treatment_cell_summaries
    outcome_confounders <- confounder_spec$outcome_confounders
  }

  if (!is.null(cell_types)) {
    meta <- meta[meta$.CellType %in% cell_types, , drop = FALSE]
  }

  needed_cols <- unique(c(
    treatment_confounders,
    treatment_cell_summaries,
    outcome_confounders
  ))
  needed_cols <- setdiff(needed_cols, c("", NA))

  missing_needed <- setdiff(needed_cols, colnames(meta))
  if (length(missing_needed) > 0L) {
    stop("Missing confounder columns in meta: ", paste(missing_needed, collapse = ", "), call. = FALSE)
  }

  if (isTRUE(drop_missing) && length(needed_cols) > 0L) {
    keep_complete <- stats::complete.cases(meta[, needed_cols, drop = FALSE])
    meta <- meta[keep_complete, , drop = FALSE]
  }

  if (nrow(meta) == 0L) {
    stop("No metadata rows available after filtering and complete-case selection.", call. = FALSE)
  }

  donor_consistency <- tapply(meta$.D, meta$.donor_id, function(x) length(unique(x)))
  if (any(donor_consistency != 1L)) {
    bad <- names(donor_consistency)[donor_consistency != 1L]
    stop("Treatment must be constant within donor. Violations: ", paste(utils::head(bad, 10), collapse = ", "), call. = FALSE)
  }

  list(
    meta = meta,
    treatment_confounders = treatment_confounders,
    treatment_cell_summaries = treatment_cell_summaries,
    outcome_confounders = outcome_confounders,
    needed_cols = needed_cols
  )
}

.dedml_build_donor_design <- function(
    meta,
    treatment_confounders,
    treatment_cell_summaries,
    n_folds,
    seed) {
  if (length(treatment_confounders) > 0L) {
    donor_meta <- dplyr::summarise(
      dplyr::group_by(meta, .donor_id),
      .D = .D[1],
      dplyr::across(dplyr::all_of(treatment_confounders), ~ .x[1]),
      .groups = "drop"
    )
  } else {
    donor_meta <- dplyr::summarise(
      dplyr::group_by(meta, .donor_id),
      .D = .D[1],
      .groups = "drop"
    )
  }

  if (length(treatment_cell_summaries) > 0L) {
    donor_summ <- dplyr::summarise(
      dplyr::group_by(meta, .donor_id),
      dplyr::across(dplyr::all_of(treatment_cell_summaries), ~ mean(as.numeric(.x), na.rm = TRUE), .names = "mean_{.col}"),
      .groups = "drop"
    )
    donor_meta <- dplyr::left_join(donor_meta, donor_summ, by = ".donor_id")
    treatment_features <- c(treatment_confounders, paste0("mean_", treatment_cell_summaries))
  } else {
    treatment_features <- treatment_confounders
  }

  fold_df <- make_stratified_donor_folds(
    donor_id = donor_meta$.donor_id,
    treatment = donor_meta$.D,
    n_folds = n_folds,
    seed = seed
  )

  donor_meta <- dplyr::left_join(donor_meta, fold_df, by = c(".donor_id" = "donor_id"))
  meta <- dplyr::left_join(meta, fold_df, by = c(".donor_id" = "donor_id"))

  list(
    meta = meta,
    donor_meta = donor_meta,
    treatment_features = treatment_features,
    fold_df = fold_df
  )
}

#' Prepare Metadata for DEDML Input
#'
#' Standardizes ID and phenotype fields used by [dedml_fit()] and optionally
#' drops incomplete rows for selected confounders.
#'
#' @param meta Metadata data.frame.
#' @param donor_id Donor id column name or vector.
#' @param sample_id Sample id column name or vector. If `NULL`, donor id is used.
#' @param treatment Treatment column name or vector (binary 0/1).
#' @param cell_type Cell type column name or vector.
#' @param donor_confounders Optional donor-level confounders to validate.
#' @param cell_confounders Optional cell-level confounders to validate.
#' @param drop_missing Logical; if `TRUE`, removes rows with missing values in
#' required fields and selected confounders.
#'
#' @return Metadata data.frame with standardized columns:
#' `.donor_id`, `.sample_id`, `.D`, `.CellType`.
#' @export
#'
dedml_prepare_metadata <- function(
    meta,
    donor_id = "HATIMID",
    sample_id = "sample_id",
    treatment = "D",
    cell_type = "CellType",
    donor_confounders = NULL,
    cell_confounders = NULL,
    drop_missing = TRUE) {
  meta <- as.data.frame(meta, stringsAsFactors = FALSE)

  .resolve <- function(x, label) {
    if (length(x) == 1L && is.character(x) && x %in% colnames(meta)) {
      return(meta[[x]])
    }
    if (length(x) == nrow(meta)) {
      return(x)
    }
    stop(label, " must be a metadata column name or a vector with nrow(meta) elements.", call. = FALSE)
  }

  if (is.null(sample_id)) {
    sample_id <- donor_id
  }

  meta$.donor_id <- as.character(.resolve(donor_id, "donor_id"))
  meta$.sample_id <- as.character(.resolve(sample_id, "sample_id"))
  meta$.D <- as.integer(.resolve(treatment, "treatment"))
  meta$.CellType <- as.character(.resolve(cell_type, "cell_type"))

  if (!all(meta$.D %in% c(0L, 1L, NA_integer_))) {
    stop("treatment must be binary 0/1.", call. = FALSE)
  }

  needed <- unique(c(
    ".donor_id", ".sample_id", ".D", ".CellType",
    donor_confounders, cell_confounders
  ))
  needed <- needed[nzchar(needed)]

  missing_cols <- setdiff(c(donor_confounders, cell_confounders), colnames(meta))
  if (length(missing_cols) > 0L) {
    stop("Missing confounder columns in meta: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  if (isTRUE(drop_missing)) {
    keep <- stats::complete.cases(meta[, needed, drop = FALSE])
    meta <- meta[keep, , drop = FALSE]
  }

  meta
}

#' Summarize and Validate a DEDML Design
#'
#' Checks donor balance, cross-fitting folds, cell-type representation, and
#' treatment overlap for categorical and numeric confounders before running the
#' expensive gene-level DEDML fit. Categorical confounder levels observed in
#' only one treatment group are treated as fatal errors because the treatment
#' effect is not identifiable for those levels. This catches common failures
#' such as disease status being completely separated by batch, site, chemistry,
#' or another design variable.
#'
#' @param meta Metadata data.frame.
#' @param donor_id Donor id column name or vector.
#' @param sample_id Sample id column name or vector. If `NULL`, donor id is used.
#' @param treatment Treatment column name or vector (binary 0/1).
#' @param cell_type Cell type column name or vector.
#' @param confounder_spec Optional list from [dedml_make_confounder_spec()].
#' @param treatment_confounders Donor-level confounders for the treatment model.
#' @param treatment_cell_summaries Cell-level covariates summarized by donor for
#' the treatment model.
#' @param outcome_confounders Covariates for the outcome model.
#' @param cell_types Optional cell types to keep before summarizing.
#' @param n_folds Number of donor-blocked cross-fitting folds.
#' @param donor_model Treatment model learner, either `"glm"` or `"lightgbm"`.
#' @param treatment_params Optional treatment learner hyperparameters.
#' @param seed Random seed for fold construction.
#' @param drop_missing Whether to drop incomplete rows for selected confounders,
#' matching [dedml_fit()] behavior.
#' @param stop_on_error Whether to stop when fatal design diagnostics are found.
#'
#' @return A `dedml_design_summary` list with `summary`, `treatment`, `folds`,
#' `cell_types`, `categorical`, `numeric`, and `diagnostics` tables.
#' @export
#'
dedml_summarize_design <- function(
    meta,
    donor_id = "HATIMID",
    sample_id = "sample_id",
    treatment = "D",
    cell_type = "CellType",
    confounder_spec = NULL,
    treatment_confounders = c("Age", "Sex"),
    treatment_cell_summaries = c("log_nCount_RNA", "nFeature_RNA", "percent.mt"),
    outcome_confounders = c("Age", "Sex", "log_nCount_RNA", "nFeature_RNA", "percent.mt"),
    cell_types = NULL,
    n_folds = 3L,
    donor_model = c("glm", "lightgbm"),
    treatment_params = list(),
    seed = 123L,
    drop_missing = TRUE,
    stop_on_error = TRUE) {
  donor_model <- match.arg(donor_model)
  n_folds <- as.integer(n_folds)
  if (is.na(n_folds) || n_folds < 2L) {
    stop("n_folds must be at least 2.", call. = FALSE)
  }

  design_meta <- .dedml_prepare_design_metadata(
    meta = meta,
    donor_id = donor_id,
    sample_id = sample_id,
    treatment = treatment,
    cell_type = cell_type,
    confounder_spec = confounder_spec,
    treatment_confounders = treatment_confounders,
    treatment_cell_summaries = treatment_cell_summaries,
    outcome_confounders = outcome_confounders,
    cell_types = cell_types,
    drop_missing = drop_missing
  )

  donor_design <- .dedml_build_donor_design(
    meta = design_meta$meta,
    treatment_confounders = design_meta$treatment_confounders,
    treatment_cell_summaries = design_meta$treatment_cell_summaries,
    n_folds = n_folds,
    seed = seed
  )

  x_treat <- .build_design_matrix(donor_design$donor_meta, donor_design$treatment_features)
  diagnostics <- .dedml_preflight_diagnostics(
    meta = donor_design$meta,
    donor_meta = donor_design$donor_meta,
    treatment_confounders = design_meta$treatment_confounders,
    treatment_features = donor_design$treatment_features,
    outcome_confounders = design_meta$outcome_confounders,
    n_folds = n_folds,
    donor_model = donor_model,
    treatment_params = treatment_params,
    x_treat = x_treat,
    separated_categorical_severity = "error"
  )

  report <- .dedml_design_report(
    meta = donor_design$meta,
    donor_meta = donor_design$donor_meta,
    treatment_confounders = design_meta$treatment_confounders,
    treatment_cell_summaries = design_meta$treatment_cell_summaries,
    treatment_features = donor_design$treatment_features,
    outcome_confounders = design_meta$outcome_confounders,
    diagnostics = diagnostics
  )

  if (isTRUE(stop_on_error)) {
    .dedml_stop_for_fatal_diagnostics(report$diagnostics)
  }
  report
}

#' @export
print.dedml_design_summary <- function(x, ...) {
  cat("DEDML design summary\n")
  print(x$summary, row.names = FALSE)

  if (nrow(x$treatment) > 0L) {
    cat("\nTreatment balance\n")
    print(x$treatment, row.names = FALSE)
  }

  if (nrow(x$diagnostics) > 0L) {
    cat("\nDiagnostics\n")
    print(x$diagnostics, row.names = FALSE)
  } else {
    cat("\nDiagnostics: no issues detected\n")
  }

  separated <- x$categorical[x$categorical$separated %in% TRUE, , drop = FALSE]
  if (nrow(separated) > 0L) {
    cat("\nSeparated categorical levels\n")
    print(separated, row.names = FALSE)
  }

  invisible(x)
}

make_stratified_donor_folds <- function(donor_id, treatment, n_folds = 3L, seed = 123L) {
  stopifnot(length(donor_id) == length(treatment))
  donor_id <- as.character(donor_id)
  treatment <- as.integer(treatment)

  donor_tbl <- unique(data.frame(donor_id = donor_id, treatment = treatment, stringsAsFactors = FALSE))
  if (any(is.na(donor_tbl$treatment))) {
    stop("Treatment contains NA at donor level.", call. = FALSE)
  }

  if (!all(donor_tbl$treatment %in% c(0L, 1L))) {
    stop("Treatment must be binary 0/1 at donor level.", call. = FALSE)
  }

  if (n_folds < 2L) {
    stop("n_folds must be at least 2.", call. = FALSE)
  }

  set.seed(seed)
  donors0 <- donor_tbl$donor_id[donor_tbl$treatment == 0L]
  donors1 <- donor_tbl$donor_id[donor_tbl$treatment == 1L]

  donors0 <- sample(donors0)
  donors1 <- sample(donors1)

  fold0 <- rep(seq_len(n_folds), length.out = length(donors0))
  fold1 <- rep(seq_len(n_folds), length.out = length(donors1))

  out <- rbind(
    data.frame(donor_id = donors0, fold = fold0, stringsAsFactors = FALSE),
    data.frame(donor_id = donors1, fold = fold1, stringsAsFactors = FALSE)
  )
  out
}

.build_design_matrix <- function(df, cols) {
  if (length(cols) == 0L) {
    return(matrix(numeric(0), nrow = nrow(df), ncol = 0))
  }

  mm <- stats::model.matrix(
    stats::as.formula(paste("~", paste(cols, collapse = " + "))),
    data = df,
    na.action = stats::na.pass
  )
  mm <- mm[, colnames(mm) != "(Intercept)", drop = FALSE]
  storage.mode(mm) <- "double"
  mm
}

.resolve_input_column <- function(meta, value_or_name, label) {
  if (length(value_or_name) == 1L && is.character(value_or_name) && value_or_name %in% colnames(meta)) {
    return(list(values = meta[[value_or_name]], source = value_or_name))
  }
  if (length(value_or_name) == nrow(meta)) {
    return(list(values = value_or_name, source = label))
  }
  stop(label, " must be a metadata column name or a vector with nrow(meta) elements.", call. = FALSE)
}

.dedml_diag_row <- function(severity, stage, check, message, recommendation) {
  data.frame(
    severity = severity,
    stage = stage,
    check = check,
    message = message,
    recommendation = recommendation,
    stringsAsFactors = FALSE
  )
}

.dedml_emit_diagnostic_warnings <- function(diagnostics, max_items = 6L) {
  if (!is.data.frame(diagnostics) || nrow(diagnostics) == 0L) {
    return(invisible(FALSE))
  }
  warn_df <- diagnostics[diagnostics$severity %in% c("warning", "error"), , drop = FALSE]
  if (nrow(warn_df) == 0L) {
    return(invisible(FALSE))
  }
  show_df <- utils::head(warn_df, max_items)
  bullets <- .dedml_format_diagnostic_bullets(show_df, max_items = max_items)
  more <- if (nrow(warn_df) > nrow(show_df)) {
    paste0("\n- ... ", nrow(warn_df) - nrow(show_df), " more issue(s); inspect fit$diagnostics.")
  } else {
    "\nInspect fit$diagnostics for the full diagnostic table."
  }
  warning(
    paste0(
      "DEDML design diagnostics detected ", nrow(warn_df),
      " potential issue(s):\n",
      paste(bullets, collapse = "\n"),
      more
    ),
    call. = FALSE
  )
  invisible(TRUE)
}

.dedml_effective_lgb_min_data_in_leaf <- function(params) {
  candidates <- c("min_data_in_leaf", "min_data", "min_child_samples")
  for (nm in candidates) {
    if (!is.null(params[[nm]])) {
      val <- suppressWarnings(as.integer(params[[nm]])[1])
      if (is.finite(val) && val > 0L) {
        return(val)
      }
    }
  }
  50L
}

.dedml_donor_numeric_summary <- function(meta, vars) {
  vars <- vars[vars %in% colnames(meta)]
  vars <- vars[vapply(meta[vars], function(x) is.numeric(x) || is.integer(x), logical(1L))]
  if (length(vars) == 0L) {
    return(NULL)
  }
  out <- dplyr::summarise(
    dplyr::group_by(meta, .donor_id),
    .D = .D[1],
    dplyr::across(dplyr::all_of(vars), ~ mean(as.numeric(.x), na.rm = TRUE)),
    .groups = "drop"
  )
  as.data.frame(out)
}

.dedml_design_report <- function(
    meta,
    donor_meta,
    treatment_confounders,
    treatment_cell_summaries,
    treatment_features,
    outcome_confounders,
    diagnostics) {
  treatment_summary <- dplyr::summarise(
    dplyr::group_by(meta, .D),
    donors = dplyr::n_distinct(.donor_id),
    samples = dplyr::n_distinct(.sample_id),
    cells = dplyr::n(),
    .groups = "drop"
  )
  treatment_summary$treatment <- ifelse(treatment_summary$.D == 1L, "treated", "control")
  treatment_summary <- treatment_summary[, c("treatment", ".D", "donors", "samples", "cells"), drop = FALSE]

  fold_summary <- dplyr::summarise(
    dplyr::group_by(donor_meta, fold),
    donors = dplyr::n(),
    treated_donors = sum(.D == 1L),
    control_donors = sum(.D == 0L),
    .groups = "drop"
  )

  celltype_summary <- dplyr::summarise(
    dplyr::group_by(meta, .CellType),
    cells = dplyr::n(),
    donors = dplyr::n_distinct(.donor_id),
    samples = dplyr::n_distinct(.sample_id),
    treated_donors = dplyr::n_distinct(.donor_id[.D == 1L]),
    control_donors = dplyr::n_distinct(.donor_id[.D == 0L]),
    .groups = "drop"
  )
  names(celltype_summary)[names(celltype_summary) == ".CellType"] <- "cell_type"

  cat_vars <- unique(c(treatment_confounders, outcome_confounders))
  cat_vars <- cat_vars[cat_vars %in% colnames(meta)]
  cat_vars <- cat_vars[vapply(meta[cat_vars], function(x) {
    is.factor(x) || is.character(x) || is.logical(x)
  }, logical(1L))]
  categorical_summary <- lapply(cat_vars, function(var) {
    df <- data.frame(
      variable = var,
      .donor_id = meta$.donor_id,
      .D = meta$.D,
      level = as.character(meta[[var]]),
      stringsAsFactors = FALSE
    )
    df <- df[!is.na(df$level) & nzchar(df$level), , drop = FALSE]
    if (nrow(df) == 0L) {
      return(NULL)
    }
    donor_level <- unique(df[, c("variable", ".donor_id", ".D", "level"), drop = FALSE])
    out <- dplyr::summarise(
      dplyr::group_by(donor_level, variable, level),
      donors = dplyr::n_distinct(.donor_id),
      treated_donors = dplyr::n_distinct(.donor_id[.D == 1L]),
      control_donors = dplyr::n_distinct(.donor_id[.D == 0L]),
      .groups = "drop"
    )
    cell_counts <- dplyr::summarise(
      dplyr::group_by(df, variable, level),
      cells = dplyr::n(),
      .groups = "drop"
    )
    out <- dplyr::left_join(out, cell_counts, by = c("variable", "level"))
    out$role <- vapply(out$variable, function(v) {
      in_treat <- v %in% treatment_confounders
      in_outcome <- v %in% outcome_confounders
      if (in_treat && in_outcome) {
        "treatment,outcome"
      } else if (in_treat) {
        "treatment"
      } else {
        "outcome"
      }
    }, character(1L))
    out$separated <- out$treated_donors == 0L | out$control_donors == 0L
    out[, c("variable", "role", "level", "donors", "treated_donors", "control_donors", "cells", "separated"), drop = FALSE]
  })
  categorical_summary <- categorical_summary[!vapply(categorical_summary, is.null, logical(1L))]
  categorical_summary <- if (length(categorical_summary) > 0L) {
    as.data.frame(data.table::rbindlist(categorical_summary, fill = TRUE))
  } else {
    data.frame(
      variable = character(),
      role = character(),
      level = character(),
      donors = integer(),
      treated_donors = integer(),
      control_donors = integer(),
      cells = integer(),
      separated = logical(),
      stringsAsFactors = FALSE
    )
  }

  num_vars <- unique(c(treatment_confounders, treatment_cell_summaries, treatment_features, outcome_confounders))
  donor_num <- .dedml_donor_numeric_summary(meta, num_vars)
  numeric_summary <- if (!is.null(donor_num)) {
    out <- lapply(setdiff(colnames(donor_num), c(".donor_id", ".D")), function(var) {
      x0 <- donor_num[[var]][donor_num$.D == 0L]
      x1 <- donor_num[[var]][donor_num$.D == 1L]
      x0 <- x0[is.finite(x0)]
      x1 <- x1[is.finite(x1)]
      if (length(x0) == 0L || length(x1) == 0L) {
        return(NULL)
      }
      data.frame(
        variable = var,
        control_min = min(x0),
        control_max = max(x0),
        treated_min = min(x1),
        treated_max = max(x1),
        control_mean = mean(x0),
        treated_mean = mean(x1),
        overlap = !(max(x0) < min(x1) || max(x1) < min(x0)),
        stringsAsFactors = FALSE
      )
    })
    out <- out[!vapply(out, is.null, logical(1L))]
    if (length(out) > 0L) {
      as.data.frame(data.table::rbindlist(out, fill = TRUE))
    } else {
      NULL
    }
  } else {
    NULL
  }
  if (is.null(numeric_summary)) {
    numeric_summary <- data.frame(
      variable = character(),
      control_min = numeric(),
      control_max = numeric(),
      treated_min = numeric(),
      treated_max = numeric(),
      control_mean = numeric(),
      treated_mean = numeric(),
      overlap = logical(),
      stringsAsFactors = FALSE
    )
  }

  summary <- data.frame(
    n_cells = nrow(meta),
    n_donors = nrow(donor_meta),
    n_samples = dplyr::n_distinct(meta$.sample_id),
    n_celltypes = dplyr::n_distinct(meta$.CellType),
    n_treated_donors = sum(donor_meta$.D == 1L),
    n_control_donors = sum(donor_meta$.D == 0L),
    n_diagnostics = nrow(diagnostics),
    n_errors = sum(diagnostics$severity == "error"),
    n_warnings = sum(diagnostics$severity == "warning"),
    stringsAsFactors = FALSE
  )

  out <- list(
    summary = summary,
    treatment = as.data.frame(treatment_summary),
    folds = as.data.frame(fold_summary),
    cell_types = as.data.frame(celltype_summary),
    categorical = categorical_summary,
    numeric = numeric_summary,
    diagnostics = diagnostics
  )
  class(out) <- "dedml_design_summary"
  out
}

.dedml_preflight_diagnostics <- function(
    meta,
    donor_meta,
    treatment_confounders,
    treatment_features,
    outcome_confounders,
    n_folds,
    donor_model,
    treatment_params,
    x_treat = NULL,
    separated_categorical_severity = "error") {
  separated_categorical_severity <- match.arg(separated_categorical_severity, c("error", "warning"))
  diag_list <- list()
  add <- function(...) {
    diag_list[[length(diag_list) + 1L]] <<- .dedml_diag_row(...)
  }

  n_donors <- nrow(donor_meta)
  n_treated <- sum(donor_meta$.D == 1L, na.rm = TRUE)
  n_control <- sum(donor_meta$.D == 0L, na.rm = TRUE)
  min_group <- min(n_treated, n_control)

  if (n_treated == 0L || n_control == 0L) {
    add(
      "error", "treatment_model", "one_treatment_class",
      paste0("Only one treatment class is present at donor level: treated=", n_treated, ", control=", n_control, "."),
      "DEDML requires both treated and control donors; check treatment coding and cohort selection."
    )
  } else if (min_group < 5L) {
    add(
      "warning", "treatment_model", "small_treatment_group",
      paste0("Smallest donor treatment group has ", min_group, " donor(s) out of ", n_donors, "."),
      "Use fewer folds, simplify treatment confounders, or collect more donors before relying on treatment residualization."
    )
  }

  fold_tab <- table(donor_meta$fold, donor_meta$.D)
  for (val in c("0", "1")) {
    if (!val %in% colnames(fold_tab)) {
      fold_tab <- cbind(fold_tab, stats::setNames(rep(0L, nrow(fold_tab)), val))
    }
  }
  test_min_treated <- min(fold_tab[, "1"], na.rm = TRUE)
  test_min_control <- min(fold_tab[, "0"], na.rm = TRUE)
  train_counts <- lapply(sort(unique(donor_meta$fold)), function(k) {
    idx <- donor_meta$fold != k
    c(
      fold = k,
      n_train = sum(idx),
      train_treated = sum(donor_meta$.D[idx] == 1L),
      train_control = sum(donor_meta$.D[idx] == 0L)
    )
  })
  train_counts <- as.data.frame(do.call(rbind, train_counts))
  min_train_treated <- min(train_counts$train_treated, na.rm = TRUE)
  min_train_control <- min(train_counts$train_control, na.rm = TRUE)
  min_train_total <- min(train_counts$n_train, na.rm = TRUE)

  if (n_treated < n_folds || n_control < n_folds || test_min_treated == 0L || test_min_control == 0L) {
    add(
      "warning", "cross_fitting", "fold_missing_class",
      paste0(
        "At least one held-out fold lacks treated or control donors ",
        "(min held-out treated=", test_min_treated, ", min held-out control=", test_min_control, ")."
      ),
      "Reduce n_folds so every fold has both treatment groups, or use more donors."
    )
  }
  if (min_train_treated == 0L || min_train_control == 0L) {
    add(
      "error", "cross_fitting", "training_fold_missing_class",
      paste0(
        "At least one treatment-model training split lacks a class ",
        "(min training treated=", min_train_treated, ", min training control=", min_train_control, ")."
      ),
      "Use fewer folds or increase the smaller treatment group; treatment propensities cannot be learned from one class."
    )
  } else if (min(min_train_treated, min_train_control) < 5L) {
    add(
      "warning", "cross_fitting", "small_training_class",
      paste0(
        "Smallest treatment-model training split has treated=", min_train_treated,
        " and control=", min_train_control, " donor(s)."
      ),
      "Use fewer folds and a simpler treatment model; cross-fitted propensity estimates may be unstable."
    )
  }

  if (donor_model == "lightgbm") {
    min_leaf <- .dedml_effective_lgb_min_data_in_leaf(treatment_params)
    if (min_train_total < 2L * min_leaf) {
      add(
        "warning", "treatment_model", "lightgbm_cannot_split",
        paste0(
          "LightGBM treatment model has min_data_in_leaf=", min_leaf,
          " but the smallest training split has ", min_train_total,
          " donors; tree splits are unlikely or impossible."
        ),
        "Lower treatment_params$min_data_in_leaf, use donor_model='glm', use fewer folds, or increase donor count."
      )
    }
  }

  if (!is.null(x_treat)) {
    n_features <- ncol(x_treat)
    if (n_features >= max(1L, min_train_total - 2L)) {
      add(
        "warning", "treatment_model", "too_many_treatment_features",
        paste0("Treatment model uses ", n_features, " design columns with as few as ", min_train_total, " training donors."),
        "Reduce donor confounders/categorical levels or use regularized treatment modeling."
      )
    }
  }

  cat_vars <- unique(c(treatment_confounders, outcome_confounders))
  cat_vars <- cat_vars[cat_vars %in% colnames(meta)]
  cat_vars <- cat_vars[vapply(meta[cat_vars], function(x) {
    is.factor(x) || is.character(x) || is.logical(x)
  }, logical(1L))]
  for (var in cat_vars) {
    df <- unique(data.frame(
      .donor_id = meta$.donor_id,
      .D = meta$.D,
      level = as.character(meta[[var]]),
      stringsAsFactors = FALSE
    ))
    df <- df[!is.na(df$level) & nzchar(df$level), , drop = FALSE]
    if (nrow(df) == 0L) {
      next
    }
    lvl <- dplyr::summarise(
      dplyr::group_by(df, level),
      donors = dplyr::n_distinct(.donor_id),
      treated = dplyr::n_distinct(.donor_id[.D == 1L]),
      control = dplyr::n_distinct(.donor_id[.D == 0L]),
      .groups = "drop"
    )
    bad <- lvl[lvl$treated == 0L | lvl$control == 0L, , drop = FALSE]
    if (nrow(bad) > 0L) {
      examples <- paste0(
        utils::head(bad$level, 5L),
        "(treated=", utils::head(bad$treated, 5L),
        ", control=", utils::head(bad$control, 5L), ")"
      )
      stage <- if (var %in% treatment_confounders) "treatment_model" else "outcome_model"
      add(
        separated_categorical_severity, stage, "separated_categorical_level",
        paste0("Covariate '", var, "' has level(s) observed in only one treatment group: ", paste(examples, collapse = ", "), "."),
        "Combine sparse levels, restrict to overlapping batches/sites, or interpret disease effects as not identifiable from this covariate."
      )
    }
  }

  num_vars <- unique(c(treatment_confounders, treatment_features, outcome_confounders))
  donor_num <- .dedml_donor_numeric_summary(meta, num_vars)
  if (!is.null(donor_num)) {
    for (var in setdiff(colnames(donor_num), c(".donor_id", ".D"))) {
      x0 <- donor_num[[var]][donor_num$.D == 0L]
      x1 <- donor_num[[var]][donor_num$.D == 1L]
      x0 <- x0[is.finite(x0)]
      x1 <- x1[is.finite(x1)]
      if (length(x0) == 0L || length(x1) == 0L) {
        next
      }
      no_overlap <- max(x0) < min(x1) || max(x1) < min(x0)
      if (isTRUE(no_overlap)) {
        add(
          "warning", "confounding", "numeric_nonoverlap",
          paste0("Numeric covariate '", var, "' has non-overlapping donor-level ranges between treatment groups."),
          "Check positivity/overlap; DEDML cannot separate treatment from a confounder with no treated-control overlap."
        )
      }
    }
  }

  if (length(diag_list) == 0L) {
    return(.dedml_default_empty_diagnostics())
  }
  as.data.frame(data.table::rbindlist(diag_list, fill = TRUE))
}

.dedml_propensity_diagnostics <- function(donor_meta) {
  diag_list <- list()
  add <- function(...) {
    diag_list[[length(diag_list) + 1L]] <<- .dedml_diag_row(...)
  }
  e <- suppressWarnings(as.numeric(donor_meta$propensity_oof))
  if (length(e) == 0L || all(!is.finite(e))) {
    add(
      "warning", "treatment_model", "invalid_propensity",
      "Treatment model produced no finite cross-fitted propensity estimates.",
      "Check treatment coding, treatment model settings, and donor-level confounders."
    )
  } else {
    finite_e <- e[is.finite(e)]
    near_zero <- mean(finite_e <= 0.01)
    near_one <- mean(finite_e >= 0.99)
    sd_e <- stats::sd(finite_e)
    if (near_zero + near_one >= 0.5) {
      add(
        "warning", "treatment_model", "boundary_propensity",
        paste0(
          round(100 * (near_zero + near_one), 1),
          "% of donors have propensity <=0.01 or >=0.99."
        ),
        "This suggests separation, severe imbalance, or overfitting; simplify confounders, use fewer folds, or use a more stable treatment model."
      )
    }
    if (is.finite(sd_e) && sd_e < 0.01) {
      add(
        "warning", "treatment_model", "near_constant_propensity",
        paste0("Cross-fitted propensities are nearly constant (sd=", signif(sd_e, 3), ")."),
        "Treatment residualization is close to using raw treatment; check whether the treatment model is too constrained or underfit."
      )
    }
    e_treated <- e[donor_meta$.D == 1L]
    e_control <- e[donor_meta$.D == 0L]
    if (length(e_treated) > 0L && mean(e_treated <= 0.05, na.rm = TRUE) >= 0.5) {
      add(
        "warning", "treatment_model", "treated_low_propensity",
        "Most treated donors have fitted propensity <=0.05.",
        "Check treatment coding and positivity; treated donors look impossible under the fitted treatment model."
      )
    }
    if (length(e_control) > 0L && mean(e_control >= 0.95, na.rm = TRUE) >= 0.5) {
      add(
        "warning", "treatment_model", "control_high_propensity",
        "Most control donors have fitted propensity >=0.95.",
        "Check treatment coding and positivity; control donors look impossible under the fitted treatment model."
      )
    }
  }
  if (length(diag_list) == 0L) {
    return(.dedml_default_empty_diagnostics())
  }
  as.data.frame(data.table::rbindlist(diag_list, fill = TRUE))
}

.fit_predict_glm_binomial <- function(x_train, y_train, x_test) {
  x_train_i <- cbind("(Intercept)" = 1, x_train)
  fit <- suppressWarnings(stats::glm.fit(x = x_train_i, y = y_train, family = stats::binomial()))

  beta <- fit$coefficients
  beta[is.na(beta)] <- 0

  pred <- stats::plogis(cbind("(Intercept)" = 1, x_test) %*% beta)
  pmin(pmax(as.numeric(pred), 1e-4), 1 - 1e-4)
}

.fit_predict_glm_gaussian <- function(x_train, y_train, x_test) {
  x_train_i <- cbind("(Intercept)" = 1, x_train)
  fit <- suppressWarnings(stats::glm.fit(x = x_train_i, y = y_train, family = stats::gaussian()))

  beta <- fit$coefficients
  beta[is.na(beta)] <- 0

  pred <- cbind("(Intercept)" = 1, x_test) %*% beta
  pred <- as.numeric(pred)
  attr(pred, "nb_size") <- NA_real_
  pred
}

.fit_predict_glm_poisson <- function(x_train, y_train, x_test) {
  x_train_i <- cbind("(Intercept)" = 1, x_train)
  fit <- suppressWarnings(stats::glm.fit(x = x_train_i, y = y_train, family = stats::poisson()))

  beta <- fit$coefficients
  beta[is.na(beta)] <- 0

  pred <- exp(cbind("(Intercept)" = 1, x_test) %*% beta)
  pred <- pmax(as.numeric(pred), 1e-6)
  attr(pred, "nb_size") <- NA_real_
  pred
}

.fit_glm_poisson_coef_one <- function(x, y, maxit = 25L, epsilon = 1e-8) {
  y <- as.numeric(y)
  n <- length(y)
  p <- ncol(x)
  fam <- stats::poisson()
  dev_resids <- fam$dev.resids
  tol <- min(1e-7, epsilon / 1000)

  mustart <- y + 0.1
  eta <- log(mustart)
  mu <- exp(eta)
  devold <- sum(dev_resids(y, mu, rep.int(1, n)))
  coef <- rep.int(0, p)
  coefold <- NULL
  rank <- p
  pivot <- seq_len(p)

  for (iter in seq_len(maxit)) {
    good <- is.finite(mu) & mu > 0
    if (!all(good)) {
      stop("invalid Poisson IRLS mean", call. = FALSE)
    }

    z <- eta + (y - mu) / mu
    w <- sqrt(mu)
    fit <- .Call(stats:::C_Cdqrls, x * w, z * w, tol, FALSE)
    if (any(!is.finite(fit$coefficients))) {
      stop("non-finite Poisson IRLS coefficients", call. = FALSE)
    }

    start <- numeric(p)
    start[fit$pivot] <- fit$coefficients
    eta_new <- as.numeric(x %*% start)
    mu_new <- exp(eta_new)
    dev <- sum(dev_resids(y, mu_new, rep.int(1, n)))

    if (!is.finite(dev)) {
      if (is.null(coefold)) {
        stop("no valid Poisson IRLS coefficient set", call. = FALSE)
      }
      ii <- 1L
      while (!is.finite(dev)) {
        if (ii > maxit) {
          stop("Poisson IRLS step halving failed", call. = FALSE)
        }
        ii <- ii + 1L
        start <- (start + coefold) / 2
        eta_new <- as.numeric(x %*% start)
        mu_new <- exp(eta_new)
        dev <- sum(dev_resids(y, mu_new, rep.int(1, n)))
      }
    }

    rank <- fit$rank
    pivot <- fit$pivot
    coef <- start
    eta <- eta_new
    mu <- mu_new
    if (abs(dev - devold) / (0.1 + abs(dev)) < epsilon) {
      break
    }
    devold <- dev
    coefold <- coef
  }

  if (rank < p) {
    coef_p <- coef[pivot]
    coef_p[seq.int(rank + 1L, p)] <- NA_real_
    coef[pivot] <- coef_p
  }
  coef
}

.fit_predict_glm_poisson_matrix <- function(x_train, y_train_mat, x_test, maxit = 25L, epsilon = 1e-8) {
  x_train_i <- cbind("(Intercept)" = 1, x_train)
  x_test_i <- cbind("(Intercept)" = 1, x_test)
  y_train_mat <- as.matrix(y_train_mat)
  storage.mode(y_train_mat) <- "double"

  coef_mat <- vapply(
    seq_len(ncol(y_train_mat)),
    function(j) {
      coef <- tryCatch(
        .fit_glm_poisson_coef_one(
          x = x_train_i,
          y = y_train_mat[, j],
          maxit = maxit,
          epsilon = epsilon
        ),
        error = function(e) {
          fit <- suppressWarnings(stats::glm.fit(x = x_train_i, y = y_train_mat[, j], family = stats::poisson()))
          coef <- fit$coefficients
          coef[is.na(coef)] <- 0
          coef
        }
      )
      coef[is.na(coef)] <- 0
      coef
    },
    numeric(ncol(x_train_i))
  )
  pred <- x_test_i %*% coef_mat
  pmax(exp(pred), 1e-6)
}

.fit_predict_gam_nb <- function(x_train, y_train, x_test, params = list()) {
  if (!requireNamespace("mgcv", quietly = TRUE)) {
    stop("Package 'mgcv' is required for outcome_learner='glm' and outcome_distribution='nb'.", call. = FALSE)
  }

  nb_train_df_base <- params$.nb_train_df_base
  nb_test_df <- params$.nb_test_df
  nb_form <- params$.nb_form
  params$.nb_train_df_base <- NULL
  params$.nb_test_df <- NULL
  params$.nb_form <- NULL

  if (!is.null(nb_train_df_base) && !is.null(nb_test_df) && !is.null(nb_form)) {
    train_df <- nb_train_df_base
    train_df$.y <- as.numeric(y_train)
    test_df <- nb_test_df
    form <- nb_form
  } else {
    x_train <- as.data.frame(x_train, check.names = FALSE)
    x_test <- as.data.frame(x_test, check.names = FALSE)
    covar_names <- make.names(colnames(x_train), unique = TRUE)
    colnames(x_train) <- covar_names
    colnames(x_test) <- covar_names

    train_df <- data.frame(.y = as.numeric(y_train), x_train, check.names = FALSE)
    test_df <- data.frame(x_test, check.names = FALSE)
    form <- if (length(covar_names) == 0L) {
      stats::as.formula(".y ~ 1")
    } else {
      stats::as.formula(paste(".y ~", paste(covar_names, collapse = " + ")))
    }
  }

  gam_args <- utils::modifyList(
    list(formula = form, data = train_df, family = mgcv::nb(), method = "REML"),
    params
  )
  fit <- suppressWarnings(do.call(mgcv::gam, gam_args))
  pred <- as.numeric(stats::predict(fit, newdata = test_df, type = "response"))
  pred <- pmax(pred, 1e-6)
  nb_size <- tryCatch(as.numeric(fit$family$getTheta(TRUE)), error = function(e) NA_real_)
  attr(pred, "nb_size") <- nb_size
  pred
}

.fit_predict_lightgbm <- function(x_train, y_train, x_test, objective, params = list(), weight_train = NULL) {
  if (!requireNamespace("lightgbm", quietly = TRUE)) {
    stop("Package 'lightgbm' is required for learner='lightgbm'.", call. = FALSE)
  }

  default_params <- list(
    objective = objective,
    learning_rate = 0.05,
    num_leaves = 31,
    min_data_in_leaf = 50,
    feature_fraction = 0.8,
    bagging_fraction = 0.8,
    bagging_freq = 1,
    verbosity = -1,
    num_threads = 1
  )

  dtrain <- params$.lgb_dataset
  params$.lgb_dataset <- NULL
  gc_after_fit <- isTRUE(params$.gc_after_fit)
  params$.gc_after_fit <- NULL
  params$.reuse_lgb_dataset <- NULL

  merged_params <- utils::modifyList(default_params, params)
  effective_objective <- merged_params$objective
  if (is.null(merged_params$metric)) {
    merged_params$metric <- if (effective_objective == "binary") {
      "auc"
    } else if (effective_objective %in% c(
      "poisson", "negative_binomial", "nb", "negativebinomial",
      "negative_binomial_free", "nb_free", "negativebinomial_free"
    )) {
      effective_objective
    } else {
      "l2"
    }
  }
  nrounds <- if (!is.null(merged_params$nrounds)) as.integer(merged_params$nrounds) else 150L
  merged_params$nrounds <- NULL

  reuse_dataset <- !is.null(dtrain)
  if (isTRUE(reuse_dataset)) {
    dtrain$set_field("label", as.numeric(y_train))
    if (!is.null(weight_train)) {
      dtrain$set_field("weight", as.numeric(weight_train))
    }
  } else {
    dtrain <- lightgbm::lgb.Dataset(data = x_train, label = y_train, weight = weight_train)
  }
  fit <- lightgbm::lgb.train(
    params = merged_params,
    data = dtrain,
    nrounds = nrounds,
    verbose = -1
  )

  pred <- as.numeric(stats::predict(fit, x_test))
  nb_size <- NA_real_
  if (effective_objective %in% c("negative_binomial", "nb", "negativebinomial")) {
    nb_size <- if (!is.null(merged_params$negative_binomial_size)) {
      as.numeric(merged_params$negative_binomial_size)
    } else if (!is.null(merged_params$nb_size)) {
      as.numeric(merged_params$nb_size)
    } else {
      NA_real_
    }
  } else if (effective_objective %in% c("negative_binomial_free", "nb_free", "negativebinomial_free") &&
      "get_negative_binomial_size" %in% names(fit)) {
    nb_size <- tryCatch(as.numeric(fit$get_negative_binomial_size()), error = function(e) NA_real_)
  }
  if (effective_objective == "binary") {
    pred <- pmin(pmax(pred, 1e-4), 1 - 1e-4)
  } else {
    pred <- pmax(pred, 1e-6)
  }
  attr(pred, "nb_size") <- nb_size
  attr(pred, "lgb_nb_size") <- nb_size
  try(fit$.__enclos_env__$private$finalize(), silent = TRUE)
  if (!isTRUE(reuse_dataset)) {
    try(dtrain$.__enclos_env__$private$finalize(), silent = TRUE)
  }
  rm(fit)
  if (isTRUE(gc_after_fit)) {
    gc(FALSE)
  }
  pred
}

.fit_predict_treatment <- function(x_train, y_train, x_test, learner, params = list()) {
  if (learner == "glm") {
    return(.fit_predict_glm_binomial(x_train, y_train, x_test))
  }

  if (learner == "lightgbm") {
    return(.fit_predict_lightgbm(x_train, y_train, x_test, objective = "binary", params = params))
  }

  stop("Unsupported treatment learner: ", learner, call. = FALSE)
}

.fit_predict_outcome <- function(x_train, y_train, x_test, learner, params = list()) {
  distribution <- match.arg(params$outcome_distribution, c("poisson", "gaussian", "nb"))
  params$outcome_distribution <- NULL

  if (learner == "lightgbm") {
    objective <- switch(
      distribution,
      gaussian = "regression",
      poisson = "poisson",
      nb = "negative_binomial_free"
    )
    return(.fit_predict_lightgbm(x_train, y_train, x_test, objective = objective, params = params))
  }

  if (learner == "glm") {
    if (distribution == "gaussian") {
      return(.fit_predict_glm_gaussian(x_train, y_train, x_test))
    }
    if (distribution == "poisson") {
      return(.fit_predict_glm_poisson(x_train, y_train, x_test))
    }
    if (distribution == "nb") {
      return(.fit_predict_gam_nb(x_train, y_train, x_test, params = params))
    }
  }

  stop("Unsupported outcome learner: ", learner, call. = FALSE)
}

.add_lgb_datasets_to_fold_plan <- function(fold_plan, outcome_params) {
  if (!requireNamespace("lightgbm", quietly = TRUE)) {
    stop("Package 'lightgbm' is required for learner='lightgbm'.", call. = FALSE)
  }
  lapply(fold_plan, function(fold) {
    dtrain <- lightgbm::lgb.Dataset(
      data = fold$x_train,
      label = rep.int(0, nrow(fold$x_train))
    )
    fold$params <- utils::modifyList(
      outcome_params,
      list(.lgb_dataset = dtrain)
    )
    fold
  })
}

.finalize_lgb_datasets_in_fold_plan <- function(fold_plan) {
  for (fold in fold_plan) {
    dtrain <- fold$params$.lgb_dataset
    if (!is.null(dtrain)) {
      try(dtrain$.__enclos_env__$private$finalize(), silent = TRUE)
    }
  }
  invisible(TRUE)
}

.compute_outcome_residual <- function(
    y,
    mu,
    cell_type,
    outcome_distribution = c("poisson", "gaussian", "nb"),
    nb_size = NULL,
    nb_size_min = 0.1,
    nb_size_max = 1e6) {
  outcome_distribution <- match.arg(outcome_distribution)
  y <- as.numeric(y)
  mu <- as.numeric(mu)

  mu <- pmax(mu, 1e-6)

  if (outcome_distribution %in% c("gaussian", "poisson")) {
    resid <- (y - mu) / sqrt(mu)
    attr(resid, "nb_size") <- NA_real_
    return(resid)
  }

  if (outcome_distribution == "nb") {
    size_raw <- suppressWarnings(as.numeric(nb_size))
    size_valid <- size_raw[is.finite(size_raw) & size_raw > 0]
    if (length(size_valid) == 0L) {
      stop(
        "outcome_distribution='nb' requires an outcome model that returns a positive ",
        "negative-binomial dispersion/size parameter.",
        call. = FALSE
      )
    }
    size <- stats::median(size_valid)
    size <- pmin(pmax(size, nb_size_min), nb_size_max)
    var_y <- mu + mu^2 / size
    resid <- (y - mu) / sqrt(pmax(var_y, 1e-6))
    attr(resid, "nb_size") <- size
    return(resid)
  }
}

.dedml_logfc_from_residual <- function(estimate, factor) {
  x <- as.numeric(estimate) / pmax(as.numeric(factor), 1e-8)
  out <- rep(NA_real_, length(x))
  ok <- is.finite(x) & x > -0.999999
  out[ok] <- log1p(x[ok])
  out
}

.dedml_logfc_from_counterfactual_residual <- function(
    estimate,
    mu,
    propensity,
    outcome_distribution = c("poisson", "gaussian", "nb"),
    nb_size = NA_real_,
    weights = NULL) {
  outcome_distribution <- match.arg(outcome_distribution)
  estimate <- suppressWarnings(as.numeric(estimate))[1]
  mu <- pmax(suppressWarnings(as.numeric(mu)), 1e-8)
  propensity <- suppressWarnings(as.numeric(propensity))
  if (is.null(weights)) {
    weights <- rep.int(1, length(mu))
  } else {
    weights <- suppressWarnings(as.numeric(weights))
  }

  ok <- is.finite(mu) & is.finite(propensity) & is.finite(weights) & weights > 0
  mu <- mu[ok]
  propensity <- pmin(pmax(propensity[ok], 1e-6), 1 - 1e-6)
  weights <- weights[ok]
  if (!is.finite(estimate) || length(mu) == 0L) {
    return(NA_real_)
  }
  if (abs(estimate) < 1e-12) {
    return(0)
  }

  if (outcome_distribution == "nb") {
    nb_size <- suppressWarnings(as.numeric(nb_size))[1]
    if (!is.finite(nb_size) || nb_size <= 0) {
      return(NA_real_)
    }
    scale <- mu / sqrt(pmax(mu + mu^2 / nb_size, 1e-8))
  } else {
    scale <- sqrt(mu)
  }

  objective <- function(delta) {
    stats::weighted.mean(scale * delta / (1 + propensity * delta), w = weights, na.rm = TRUE) - estimate
  }

  lower <- -0.999999
  if (estimate < 0) {
    f_lower <- objective(lower)
    if (!is.finite(f_lower) || f_lower > 0) {
      return(NA_real_)
    }
    root <- stats::uniroot(objective, lower = lower, upper = 0, tol = 1e-8)$root
    return(log1p(root))
  }

  upper <- 1
  f_upper <- objective(upper)
  while (is.finite(f_upper) && f_upper < 0 && upper < 1e6) {
    upper <- upper * 2
    f_upper <- objective(upper)
  }
  if (!is.finite(f_upper) || f_upper < 0) {
    return(NA_real_)
  }
  root <- stats::uniroot(objective, lower = 0, upper = upper, tol = 1e-8)$root
  log1p(root)
}

.summarise_oof_mu_factors <- function(mu, cell_type, nb_size = NA_real_, propensity = NULL, weights = NULL) {
  mu <- pmax(as.numeric(mu), 1e-8)
  cell_type <- as.character(cell_type)
  nb_size <- suppressWarnings(as.numeric(nb_size))[1]
  has_propensity <- !is.null(propensity)
  if (is.null(propensity)) {
    propensity <- rep(NA_real_, length(mu))
  }
  propensity <- suppressWarnings(as.numeric(propensity))
  if (is.null(weights)) {
    weights <- rep.int(1, length(mu))
  }
  weights <- suppressWarnings(as.numeric(weights))
  cell_types <- sort(unique(cell_type))
  out <- lapply(cell_types, function(ct) {
    idx <- cell_type == ct
    mu_ct <- mu[idx]
    propensity_ct <- propensity[idx]
    weight_ct <- weights[idx]
    ok_w <- is.finite(weight_ct) & weight_ct > 0
    if (!any(ok_w)) {
      weight_ct <- rep.int(1, length(mu_ct))
      ok_w <- rep(TRUE, length(mu_ct))
    }
    factor_poisson <- stats::weighted.mean(sqrt(mu_ct[ok_w]), w = weight_ct[ok_w], na.rm = TRUE)
    factor_nb <- if (is.finite(nb_size) && nb_size > 0) {
      stats::weighted.mean(
        mu_ct[ok_w] / sqrt(mu_ct[ok_w] + mu_ct[ok_w]^2 / nb_size),
        w = weight_ct[ok_w],
        na.rm = TRUE
      )
    } else {
      NA_real_
    }
    ans <- data.table::data.table(
      CellType = ct,
      muhat_mean = stats::weighted.mean(mu_ct[ok_w], w = weight_ct[ok_w], na.rm = TRUE),
      muhat_factor_poisson = factor_poisson,
      muhat_factor_nb = factor_nb,
      muhat_nb_size = nb_size,
      muhat_propensity_mean = stats::weighted.mean(propensity_ct[ok_w], w = weight_ct[ok_w], na.rm = TRUE)
    )
    if (has_propensity) {
      ans$.muhat_mu <- list(mu_ct)
      ans$.muhat_propensity <- list(propensity_ct)
      ans$.muhat_weight <- list(weight_ct)
    }
    ans
  })
  data.table::rbindlist(out, fill = TRUE)
}

.summarise_oof_mu_factors_matrix <- function(
    mu_hat_mat,
    cell_type,
    nb_size = NA_real_,
    propensity = NULL,
    weights = NULL) {
  mu_hat_mat <- as.matrix(mu_hat_mat)
  gene_names <- rownames(mu_hat_mat)
  if (is.null(gene_names)) {
    gene_names <- paste0("gene_", seq_len(nrow(mu_hat_mat)))
  }
  out <- lapply(seq_len(nrow(mu_hat_mat)), function(i) {
    ans <- .summarise_oof_mu_factors(
      mu = mu_hat_mat[i, ],
      cell_type = cell_type,
      nb_size = if (length(nb_size) >= i) nb_size[i] else nb_size[1],
      propensity = propensity,
      weights = weights
    )
    ans$gene <- gene_names[i]
    ans
  })
  data.table::rbindlist(out, fill = TRUE)
}

.append_logfc_effect_sizes <- function(
    result_df,
    mu_factor_df,
    outcome_distribution = c("poisson", "gaussian", "nb")) {
  outcome_distribution <- match.arg(outcome_distribution)
  if (!is.data.frame(result_df) || nrow(result_df) == 0L ||
      is.null(mu_factor_df) || nrow(mu_factor_df) == 0L) {
    return(result_df)
  }
  result_df <- dplyr::left_join(
    result_df,
    as.data.frame(mu_factor_df),
    by = c("gene", "CellType")
  )
  result_df$muhat_factor <- if (outcome_distribution == "nb") {
    result_df$muhat_factor_nb
  } else {
    result_df$muhat_factor_poisson
  }
  has_counterfactual_inputs <- all(c(".muhat_mu", ".muhat_propensity") %in% names(result_df))
  if (has_counterfactual_inputs) {
    if (!(".muhat_weight" %in% names(result_df))) {
      result_df$.muhat_weight <- replicate(nrow(result_df), NULL, simplify = FALSE)
    }
    result_df$estimate_logfc <- mapply(
      FUN = .dedml_logfc_from_counterfactual_residual,
      estimate = result_df$estimate,
      mu = result_df$.muhat_mu,
      propensity = result_df$.muhat_propensity,
      nb_size = result_df$muhat_nb_size,
      weights = result_df$.muhat_weight,
      MoreArgs = list(outcome_distribution = outcome_distribution),
      SIMPLIFY = TRUE
    )
  } else {
    result_df$estimate_logfc <- .dedml_logfc_from_residual(
      result_df$estimate,
      result_df$muhat_factor
    )
  }
  result_df$muhat_factor_poisson <- NULL
  result_df$muhat_factor_nb <- NULL
  result_df$.muhat_mu <- NULL
  result_df$.muhat_propensity <- NULL
  result_df$.muhat_weight <- NULL
  result_df
}

fit_residual_regression_pseudobulk <- function(
    y_resid,
    d_resid,
    cell_type,
    sample_id,
    donor_id,
    min_cells_per_sample_ct = 10L,
    min_samples = 3L,
    cluster_robust = TRUE) {
  df <- data.frame(
    y_resid = as.numeric(y_resid),
    d_resid = as.numeric(d_resid),
    CellType = factor(as.character(cell_type)),
    sample_id = as.character(sample_id),
    donor_id = as.character(donor_id),
    stringsAsFactors = FALSE
  )

  sample_ct_counts <- dplyr::summarise(
    dplyr::group_by(df, sample_id, CellType),
    n_cells = dplyr::n(),
    .groups = "drop"
  )

  valid_sample_ct <- dplyr::filter(sample_ct_counts, n_cells >= min_cells_per_sample_ct)
  if (nrow(valid_sample_ct) == 0L) {
    return(NULL)
  }

  df <- dplyr::inner_join(df, valid_sample_ct, by = c("sample_id", "CellType"))
  if (nrow(df) == 0L) {
    return(NULL)
  }

  pb_df <- dplyr::summarise(
    dplyr::group_by(df, sample_id, CellType),
    y_resid_pb = mean(y_resid, na.rm = TRUE),
    d_resid = d_resid[1],
    donor_id = donor_id[1],
    n_cells = n_cells[1],
    .groups = "drop"
  )

  sample_counts_by_ct <- dplyr::summarise(
    dplyr::group_by(pb_df, CellType),
    n_samples = dplyr::n(),
    .groups = "drop"
  )

  keep_ct <- as.character(sample_counts_by_ct$CellType[sample_counts_by_ct$n_samples >= min_samples])
  pb_df <- dplyr::filter(pb_df, as.character(CellType) %in% keep_ct)

  if (nrow(pb_df) == 0L || length(unique(pb_df$d_resid)) < 2L) {
    return(NULL)
  }

  pb_df$CellType <- droplevels(pb_df$CellType)

  weights <- as.numeric(pb_df$n_cells)
  fit <- stats::lm(
    y_resid_pb ~ d_resid * CellType,
    data = pb_df,
    weights = weights
  )

  if (length(stats::coef(fit)) == 0L) {
    return(NULL)
  }

  repeated_donor_samples <- any(table(pb_df$donor_id) > 1L)
  n_unique_donors <- length(unique(pb_df$donor_id))
  use_cluster <- isTRUE(cluster_robust) && repeated_donor_samples && n_unique_donors >= 2L

  vc <- if (use_cluster) {
    sandwich::vcovCL(fit, cluster = pb_df$donor_id, type = "HC1")
  } else {
    stats::vcov(fit)
  }

  ct_levels <- levels(pb_df$CellType)
  nd_1 <- data.frame(
    d_resid = rep(1, length(ct_levels)),
    CellType = factor(ct_levels, levels = ct_levels)
  )
  nd_0 <- data.frame(
    d_resid = rep(0, length(ct_levels)),
    CellType = factor(ct_levels, levels = ct_levels)
  )

  trm_noy <- stats::delete.response(stats::terms(fit))
  x1 <- stats::model.matrix(trm_noy, data = nd_1)
  x0 <- stats::model.matrix(trm_noy, data = nd_0)

  contrast <- x1 - x0
  est <- as.numeric(contrast %*% stats::coef(fit))

  var_est <- pmax(diag(contrast %*% vc %*% t(contrast)), 0)
  std_error <- sqrt(var_est)
  statistic <- est / std_error
  p_value <- 2 * stats::pnorm(-abs(statistic))

  out <- data.frame(
    CellType = ct_levels,
    estimate = est,
    std_error = std_error,
    statistic = statistic,
    p_value = p_value,
    se_method = if (use_cluster) "donor_cluster_robust" else "model_based",
    repeated_donor_samples = repeated_donor_samples,
    stringsAsFactors = FALSE
  )

  ct_summary <- dplyr::summarise(
    dplyr::group_by(pb_df, CellType),
    n_pseudobulk = dplyr::n(),
    n_cells = sum(n_cells),
    n_samples = dplyr::n_distinct(sample_id),
    n_donors = dplyr::n_distinct(donor_id),
    .groups = "drop"
  )

  dplyr::left_join(out, ct_summary, by = "CellType")
}

fit_residual_regression_pseudobulk_matrix <- function(
    y_resid_mat,
    d_resid,
    cell_type,
    sample_id,
    donor_id,
    gene_names = NULL,
    min_cells_per_sample_ct = 10L,
    min_samples = 3L,
    cluster_robust = TRUE,
    fallback_scalar = TRUE) {
  if (is.null(dim(y_resid_mat))) {
    y_resid_mat <- matrix(as.numeric(y_resid_mat), nrow = 1L)
  }
  y_resid_mat <- as.matrix(y_resid_mat)
  storage.mode(y_resid_mat) <- "double"

  if (ncol(y_resid_mat) != length(d_resid)) {
    stop("y_resid_mat must have genes in rows and cells in columns.", call. = FALSE)
  }

  if (is.null(gene_names)) {
    gene_names <- rownames(y_resid_mat)
  }
  if (is.null(gene_names)) {
    gene_names <- paste0("gene_", seq_len(nrow(y_resid_mat)))
  }

  df <- data.frame(
    d_resid = as.numeric(d_resid),
    CellType = factor(as.character(cell_type)),
    sample_id = as.character(sample_id),
    donor_id = as.character(donor_id),
    stringsAsFactors = FALSE
  )
  df$.cell_index <- seq_len(nrow(df))

  sample_ct_counts <- dplyr::summarise(
    dplyr::group_by(df, sample_id, CellType),
    n_cells = dplyr::n(),
    .groups = "drop"
  )

  valid_sample_ct <- dplyr::filter(sample_ct_counts, n_cells >= min_cells_per_sample_ct)
  if (nrow(valid_sample_ct) == 0L) {
    return(NULL)
  }

  valid_sample_ct$.unit_key <- paste(valid_sample_ct$sample_id, valid_sample_ct$CellType, sep = "\r")
  df$.unit_key <- paste(df$sample_id, df$CellType, sep = "\r")
  keep_cell <- df$.unit_key %in% valid_sample_ct$.unit_key
  df <- df[keep_cell, , drop = FALSE]
  if (nrow(df) == 0L) {
    return(NULL)
  }

  first_idx <- match(valid_sample_ct$.unit_key, df$.unit_key)
  pb_df <- data.frame(
    sample_id = valid_sample_ct$sample_id,
    CellType = valid_sample_ct$CellType,
    d_resid = df$d_resid[first_idx],
    donor_id = df$donor_id[first_idx],
    n_cells = valid_sample_ct$n_cells,
    .unit_key = valid_sample_ct$.unit_key,
    stringsAsFactors = FALSE
  )

  sample_counts_by_ct <- dplyr::summarise(
    dplyr::group_by(pb_df, CellType),
    n_samples = dplyr::n(),
    .groups = "drop"
  )

  keep_ct <- as.character(sample_counts_by_ct$CellType[sample_counts_by_ct$n_samples >= min_samples])
  keep_unit <- as.character(pb_df$CellType) %in% keep_ct
  pb_df <- pb_df[keep_unit, , drop = FALSE]
  if (nrow(pb_df) == 0L || length(unique(pb_df$d_resid)) < 2L) {
    return(NULL)
  }

  unit_index <- match(df$.unit_key, pb_df$.unit_key)
  keep_cell <- !is.na(unit_index)
  unit_index <- unit_index[keep_cell]
  unit_factor <- factor(unit_index, levels = seq_len(nrow(pb_df)))
  cell_index <- df$.cell_index[keep_cell]

  y_resid_keep <- y_resid_mat[, cell_index, drop = FALSE]
  if (requireNamespace("Matrix", quietly = TRUE) && all(is.finite(y_resid_keep))) {
    unit_n <- tabulate(unit_index, nbins = nrow(pb_df))
    agg <- Matrix::sparseMatrix(
      i = unit_index,
      j = seq_along(unit_index),
      x = 1 / unit_n[unit_index],
      dims = c(nrow(pb_df), length(unit_index))
    )
    y_pb <- as.matrix(agg %*% t(y_resid_keep))
  } else {
    y_cell_gene <- t(y_resid_keep)
    finite_y <- is.finite(y_cell_gene)
    y_cell_gene[!finite_y] <- 0
    y_sum <- rowsum(y_cell_gene, group = unit_factor, reorder = TRUE)
    y_n <- rowsum(matrix(as.numeric(finite_y), nrow = nrow(finite_y)), group = unit_factor, reorder = TRUE)
    y_pb <- y_sum / y_n
  }
  colnames(y_pb) <- gene_names

  finite_gene <- colSums(!is.finite(y_pb)) == 0L
  batch_genes <- which(finite_gene)
  fallback_genes <- which(!finite_gene)
  out <- list()

  if (length(batch_genes) > 0L) {
    pb_df$CellType <- droplevels(pb_df$CellType)
    x <- stats::model.matrix(~ d_resid * CellType, data = pb_df)
    w <- as.numeric(pb_df$n_cells)
    y <- y_pb[, batch_genes, drop = FALSE]

    fit <- stats::lm.wfit(x = x, y = y, w = w)
    if (fit$rank < ncol(x)) {
      fallback_genes <- c(fallback_genes, batch_genes)
      batch_genes <- integer(0L)
    } else {
      coef_mat <- fit$coefficients
      if (is.null(dim(coef_mat))) {
        coef_mat <- matrix(coef_mat, ncol = 1L)
      }
      colnames(coef_mat) <- gene_names[batch_genes]

      a_inv <- chol2inv(qr.R(fit$qr))
      rownames(a_inv) <- colnames(a_inv) <- colnames(x)

      ct_levels <- levels(pb_df$CellType)
      nd_1 <- data.frame(
        d_resid = rep(1, length(ct_levels)),
        CellType = factor(ct_levels, levels = ct_levels)
      )
      nd_0 <- data.frame(
        d_resid = rep(0, length(ct_levels)),
        CellType = factor(ct_levels, levels = ct_levels)
      )
      contrast <- stats::model.matrix(~ d_resid * CellType, data = nd_1) -
        stats::model.matrix(~ d_resid * CellType, data = nd_0)
      colnames(contrast) <- colnames(x)

      est_mat <- contrast %*% coef_mat
      repeated_donor_samples <- any(table(pb_df$donor_id) > 1L)
      n_unique_donors <- length(unique(pb_df$donor_id))
      use_cluster <- isTRUE(cluster_robust) && repeated_donor_samples && n_unique_donors >= 2L

      if (use_cluster) {
        b_mat <- contrast %*% a_inv
        var_mat <- matrix(0, nrow = nrow(contrast), ncol = ncol(y))
        donor_levels <- unique(pb_df$donor_id)
        wx <- x * w
        for (donor in donor_levels) {
          idx <- which(pb_df$donor_id == donor)
          score_d <- crossprod(wx[idx, , drop = FALSE], fit$residuals[idx, , drop = FALSE])
          z_d <- b_mat %*% score_d
          var_mat <- var_mat + z_d^2
        }
        n_obs <- nrow(x)
        n_coef <- ncol(x)
        df_correction <- (n_unique_donors / (n_unique_donors - 1L)) * ((n_obs - 1L) / (n_obs - n_coef))
        var_mat <- df_correction * var_mat
        se_method <- "donor_cluster_robust"
      } else {
        sigma2 <- colSums(w * fit$residuals^2) / fit$df.residual
        contrast_var <- diag(contrast %*% a_inv %*% t(contrast))
        var_mat <- outer(contrast_var, sigma2)
        se_method <- "model_based"
      }

      std_error <- sqrt(pmax(var_mat, 0))
      statistic <- est_mat / std_error

      ct_summary <- dplyr::summarise(
        dplyr::group_by(pb_df, CellType),
        n_pseudobulk = dplyr::n(),
        n_cells = sum(n_cells),
        n_samples = dplyr::n_distinct(sample_id),
        n_donors = dplyr::n_distinct(donor_id),
        .groups = "drop"
      )
      p_value <- 2 * stats::pnorm(-abs(statistic))

      batch_out <- do.call(rbind, lapply(seq_along(batch_genes), function(j) {
        data.frame(
          CellType = ct_levels,
          estimate = as.numeric(est_mat[, j]),
          std_error = as.numeric(std_error[, j]),
          statistic = as.numeric(statistic[, j]),
          p_value = as.numeric(p_value[, j]),
          se_method = se_method,
          repeated_donor_samples = repeated_donor_samples,
          gene = gene_names[batch_genes[j]],
          stringsAsFactors = FALSE
        )
      }))
      out <- c(out, list(dplyr::left_join(batch_out, ct_summary, by = "CellType")))
    }
  }

  if (length(fallback_genes) > 0L && isTRUE(fallback_scalar)) {
    scalar_out <- lapply(fallback_genes, function(g_idx) {
      rr <- fit_residual_regression_pseudobulk(
        y_resid = y_resid_mat[g_idx, ],
        d_resid = d_resid,
        cell_type = cell_type,
        sample_id = sample_id,
        donor_id = donor_id,
        min_cells_per_sample_ct = min_cells_per_sample_ct,
        min_samples = min_samples,
        cluster_robust = cluster_robust
      )
      if (is.null(rr)) {
        return(NULL)
      }
      rr$gene <- gene_names[g_idx]
      rr
    })
    scalar_out <- Filter(Negate(is.null), scalar_out)
    out <- c(out, scalar_out)
  }

  if (length(out) == 0L) {
    return(NULL)
  }
  as.data.frame(data.table::rbindlist(out, fill = TRUE))
}

.dedml_make_pseudobulk_outcome_data <- function(
    counts,
    meta,
    outcome_confounders,
    min_cells_per_sample_ct = 10L) {
  df <- as.data.frame(meta, stringsAsFactors = FALSE)
  df$.cell_index <- seq_len(nrow(df))
  df$.unit_key <- paste(df$.sample_id, df$.CellType, sep = "\r")

  sample_ct_counts <- dplyr::summarise(
    dplyr::group_by(df, .unit_key, .sample_id, .CellType),
    n_cells = dplyr::n(),
    .groups = "drop"
  )
  valid_sample_ct <- dplyr::filter(sample_ct_counts, n_cells >= min_cells_per_sample_ct)
  if (nrow(valid_sample_ct) == 0L) {
    return(NULL)
  }

  keep_cell <- df$.unit_key %in% valid_sample_ct$.unit_key
  df_keep <- df[keep_cell, , drop = FALSE]
  if (nrow(df_keep) == 0L) {
    return(NULL)
  }

  valid_sample_ct <- valid_sample_ct[match(unique(df_keep$.unit_key), valid_sample_ct$.unit_key), , drop = FALSE]
  unit_index <- match(df_keep$.unit_key, valid_sample_ct$.unit_key)
  first_idx <- match(valid_sample_ct$.unit_key, df_keep$.unit_key)

  pb_df <- data.frame(
    sample_id = valid_sample_ct$.sample_id,
    CellType = as.character(valid_sample_ct$.CellType),
    d_resid = df_keep$D_resid[first_idx],
    donor_id = df_keep$.donor_id[first_idx],
    fold = df_keep$fold[first_idx],
    propensity_oof = df_keep$propensity_oof[first_idx],
    n_cells = valid_sample_ct$n_cells,
    .unit_key = valid_sample_ct$.unit_key,
    stringsAsFactors = FALSE
  )

  covars <- unique(outcome_confounders)
  covars <- covars[covars %in% colnames(df_keep)]
  for (covar in covars) {
    x <- df_keep[[covar]]
    if (is.numeric(x) || is.integer(x)) {
      vals <- rowsum(as.numeric(x), group = unit_index, reorder = FALSE)[, 1] / as.numeric(valid_sample_ct$n_cells)
    } else {
      vals <- vapply(seq_len(nrow(pb_df)), function(i) {
        z <- x[unit_index == i]
        z <- z[!is.na(z)]
        if (length(z) == 0L) {
          return(NA_character_)
        }
        as.character(z[1])
      }, character(1L))
      if (is.factor(x)) {
        vals <- factor(vals, levels = levels(x))
      }
    }
    pb_df[[covar]] <- vals
  }
  if ("CellType" %in% covars) {
    pb_df$CellType <- as.character(pb_df$CellType)
  }

  cell_index <- df_keep$.cell_index
  unit_n <- as.numeric(valid_sample_ct$n_cells)
  if (requireNamespace("Matrix", quietly = TRUE)) {
    agg <- Matrix::sparseMatrix(
      i = unit_index,
      j = seq_along(unit_index),
      x = 1 / unit_n[unit_index],
      dims = c(nrow(pb_df), length(unit_index))
    )
    y_pb_mat <- t(as.matrix(agg %*% t(counts[, cell_index, drop = FALSE])))
  } else {
    y_cell_gene <- t(as.matrix(counts[, cell_index, drop = FALSE]))
    y_pb_mat <- t(rowsum(y_cell_gene, group = unit_index, reorder = FALSE) / unit_n)
  }
  rownames(y_pb_mat) <- rownames(counts)
  colnames(y_pb_mat) <- pb_df$.unit_key

  list(
    pb_df = pb_df,
    y_pb_mat = y_pb_mat,
    x_pb = .build_design_matrix(pb_df, outcome_confounders)
  )
}

.fit_predict_pseudobulk_outcome <- function(
    x_train,
    y_train,
    x_test,
    weights_train,
    learner = c("lightgbm", "glm"),
    distribution = c("poisson", "gaussian", "nb"),
    params = list()) {
  learner <- match.arg(learner)
  distribution <- match.arg(distribution)
  params$outcome_distribution <- NULL
  params$.reuse_lgb_dataset <- NULL
  params$.gc_after_fit <- NULL
  weights_train <- suppressWarnings(as.numeric(weights_train))
  weights_train[!is.finite(weights_train) | weights_train <= 0] <- 1

  if (learner == "lightgbm") {
    objective <- switch(
      distribution,
      gaussian = "regression",
      poisson = "poisson",
      nb = "negative_binomial_free"
    )
    return(.fit_predict_lightgbm(
      x_train = x_train,
      y_train = y_train,
      x_test = x_test,
      objective = objective,
      params = params,
      weight_train = weights_train
    ))
  }

  x_train_i <- cbind("(Intercept)" = 1, x_train)
  x_test_i <- cbind("(Intercept)" = 1, x_test)

  if (distribution == "gaussian") {
    fit <- suppressWarnings(stats::lm.wfit(x = x_train_i, y = as.numeric(y_train), w = weights_train))
    beta <- fit$coefficients
    beta[is.na(beta)] <- 0
    pred <- as.numeric(x_test_i %*% beta)
    pred <- pmax(pred, 1e-6)
    attr(pred, "nb_size") <- NA_real_
    return(pred)
  }

  if (distribution == "poisson") {
    fit <- suppressWarnings(stats::glm.fit(
      x = x_train_i,
      y = as.numeric(y_train),
      weights = weights_train,
      family = stats::poisson()
    ))
    beta <- fit$coefficients
    beta[is.na(beta)] <- 0
    pred <- pmax(as.numeric(exp(x_test_i %*% beta)), 1e-6)
    attr(pred, "nb_size") <- NA_real_
    return(pred)
  }

  if (distribution == "nb") {
    if (!requireNamespace("mgcv", quietly = TRUE)) {
      stop("Package 'mgcv' is required for outcome_model='pseudobulk' and outcome_distribution='nb'.", call. = FALSE)
    }
    x_train_df <- as.data.frame(x_train, check.names = FALSE)
    x_test_df <- as.data.frame(x_test, check.names = FALSE)
    covar_names <- make.names(colnames(x_train_df), unique = TRUE)
    colnames(x_train_df) <- covar_names
    colnames(x_test_df) <- covar_names
    train_df <- data.frame(.y = as.numeric(y_train), .w = weights_train, x_train_df, check.names = FALSE)
    test_df <- data.frame(x_test_df, check.names = FALSE)
    form <- if (length(covar_names) == 0L) {
      stats::as.formula(".y ~ 1")
    } else {
      stats::as.formula(paste(".y ~", paste(covar_names, collapse = " + ")))
    }
    gam_args <- utils::modifyList(
      list(formula = form, data = train_df, weights = train_df$.w, family = mgcv::nb(), method = "REML"),
      params
    )
    fit <- suppressWarnings(do.call(mgcv::gam, gam_args))
    pred <- pmax(as.numeric(stats::predict(fit, newdata = test_df, type = "response")), 1e-6)
    attr(pred, "nb_size") <- tryCatch(as.numeric(fit$family$getTheta(TRUE)), error = function(e) NA_real_)
    return(pred)
  }

  stop("Unsupported pseudo-bulk outcome distribution: ", distribution, call. = FALSE)
}

fit_residual_regression_pseudobulk_matrix_from_pb <- function(
    y_resid_pb_mat,
    pb_df,
    gene_names = NULL,
    min_samples = 3L,
    cluster_robust = TRUE) {
  if (is.null(dim(y_resid_pb_mat))) {
    y_resid_pb_mat <- matrix(as.numeric(y_resid_pb_mat), nrow = 1L)
  }
  y_resid_pb_mat <- as.matrix(y_resid_pb_mat)
  storage.mode(y_resid_pb_mat) <- "double"

  if (ncol(y_resid_pb_mat) != nrow(pb_df)) {
    stop("y_resid_pb_mat must have genes in rows and pseudo-bulk units in columns.", call. = FALSE)
  }
  if (is.null(gene_names)) {
    gene_names <- rownames(y_resid_pb_mat)
  }
  if (is.null(gene_names)) {
    gene_names <- paste0("gene_", seq_len(nrow(y_resid_pb_mat)))
  }

  pb_df <- as.data.frame(pb_df, stringsAsFactors = FALSE)
  pb_df$CellType <- factor(as.character(pb_df$CellType))
  sample_counts_by_ct <- dplyr::summarise(
    dplyr::group_by(pb_df, CellType),
    n_samples = dplyr::n(),
    .groups = "drop"
  )
  keep_ct <- as.character(sample_counts_by_ct$CellType[sample_counts_by_ct$n_samples >= min_samples])
  keep_unit <- as.character(pb_df$CellType) %in% keep_ct
  pb_df <- pb_df[keep_unit, , drop = FALSE]
  y_resid_pb_mat <- y_resid_pb_mat[, keep_unit, drop = FALSE]
  if (nrow(pb_df) == 0L || length(unique(pb_df$d_resid)) < 2L) {
    return(NULL)
  }

  finite_gene <- rowSums(!is.finite(y_resid_pb_mat)) == 0L
  if (!any(finite_gene)) {
    return(NULL)
  }
  y <- t(y_resid_pb_mat[finite_gene, , drop = FALSE])
  batch_gene_names <- gene_names[finite_gene]

  pb_df$CellType <- droplevels(pb_df$CellType)
  x <- stats::model.matrix(~ d_resid * CellType, data = pb_df)
  w <- as.numeric(pb_df$n_cells)
  fit <- stats::lm.wfit(x = x, y = y, w = w)
  if (fit$rank < ncol(x)) {
    return(NULL)
  }
  coef_mat <- fit$coefficients
  if (is.null(dim(coef_mat))) {
    coef_mat <- matrix(coef_mat, ncol = 1L)
  }
  colnames(coef_mat) <- batch_gene_names

  a_inv <- chol2inv(qr.R(fit$qr))
  rownames(a_inv) <- colnames(a_inv) <- colnames(x)

  ct_levels <- levels(pb_df$CellType)
  nd_1 <- data.frame(d_resid = rep(1, length(ct_levels)), CellType = factor(ct_levels, levels = ct_levels))
  nd_0 <- data.frame(d_resid = rep(0, length(ct_levels)), CellType = factor(ct_levels, levels = ct_levels))
  contrast <- stats::model.matrix(~ d_resid * CellType, data = nd_1) -
    stats::model.matrix(~ d_resid * CellType, data = nd_0)
  colnames(contrast) <- colnames(x)

  est_mat <- contrast %*% coef_mat
  repeated_donor_samples <- any(table(pb_df$donor_id) > 1L)
  n_unique_donors <- length(unique(pb_df$donor_id))
  use_cluster <- isTRUE(cluster_robust) && repeated_donor_samples && n_unique_donors >= 2L
  if (use_cluster) {
    b_mat <- contrast %*% a_inv
    var_mat <- matrix(0, nrow = nrow(contrast), ncol = ncol(y))
    donor_levels <- unique(pb_df$donor_id)
    wx <- x * w
    for (donor in donor_levels) {
      idx <- which(pb_df$donor_id == donor)
      score_d <- crossprod(wx[idx, , drop = FALSE], fit$residuals[idx, , drop = FALSE])
      z_d <- b_mat %*% score_d
      var_mat <- var_mat + z_d^2
    }
    n_obs <- nrow(x)
    n_coef <- ncol(x)
    df_correction <- (n_unique_donors / (n_unique_donors - 1L)) * ((n_obs - 1L) / (n_obs - n_coef))
    var_mat <- df_correction * var_mat
    se_method <- "donor_cluster_robust"
  } else {
    sigma2 <- colSums(w * fit$residuals^2) / fit$df.residual
    contrast_var <- diag(contrast %*% a_inv %*% t(contrast))
    var_mat <- outer(contrast_var, sigma2)
    se_method <- "model_based"
  }

  std_error <- sqrt(pmax(var_mat, 0))
  statistic <- est_mat / std_error
  p_value <- 2 * stats::pnorm(-abs(statistic))

  ct_summary <- dplyr::summarise(
    dplyr::group_by(pb_df, CellType),
    n_pseudobulk = dplyr::n(),
    n_cells = sum(n_cells),
    n_samples = dplyr::n_distinct(sample_id),
    n_donors = dplyr::n_distinct(donor_id),
    .groups = "drop"
  )

  out <- do.call(rbind, lapply(seq_along(batch_gene_names), function(j) {
    data.frame(
      CellType = ct_levels,
      estimate = as.numeric(est_mat[, j]),
      std_error = as.numeric(std_error[, j]),
      statistic = as.numeric(statistic[, j]),
      p_value = as.numeric(p_value[, j]),
      se_method = se_method,
      repeated_donor_samples = repeated_donor_samples,
      gene = batch_gene_names[j],
      stringsAsFactors = FALSE
    )
  }))
  as.data.frame(dplyr::left_join(out, ct_summary, by = "CellType"))
}

.finalize_dedml_output <- function(
    result_list,
    error_list,
    donor_meta,
    settings,
    diagnostics = NULL,
    pvalue_calibration,
    calibration_quantile) {
  if (length(result_list) == 0L) {
    msg <- "No valid DEDML gene-level results were produced."
    if (length(error_list) > 0L) {
      err_df <- data.table::rbindlist(error_list, fill = TRUE)
      msg <- paste0(
        msg,
        " First errors: ",
        paste(utils::head(unique(err_df$message), 5L), collapse = " | ")
      )
    }
    stop(msg, call. = FALSE)
  }

  result_df <- data.table::rbindlist(result_list, fill = TRUE)
  if (pvalue_calibration == "tail") {
    result_df$p_value_uncalibrated <- result_df$p_value
    result_df$statistic_uncalibrated <- result_df$statistic
    result_df <- dplyr::group_by(result_df, CellType)
    result_df <- dplyr::mutate(
      result_df,
      tail_lambda = {
        lambda <- stats::quantile(statistic^2, probs = calibration_quantile, na.rm = TRUE, names = FALSE) /
          stats::qchisq(calibration_quantile, df = 1)
        pmax(lambda, 1)
      },
      statistic = statistic / sqrt(tail_lambda),
      p_value = 2 * stats::pnorm(-abs(statistic))
    )
    result_df <- dplyr::ungroup(result_df)
  }
  result_df <- dplyr::group_by(result_df, CellType)
  result_df <- dplyr::mutate(
    result_df,
    p_adj = stats::p.adjust(p_value, method = "BH"),
    rank_score = sign(estimate) * abs(statistic)
  )
  result_df <- dplyr::ungroup(result_df)

  list(
    results = as.data.frame(result_df),
    donor_meta = donor_meta,
    settings = settings,
    diagnostics = diagnostics,
    errors = if (length(error_list) > 0L) data.table::rbindlist(error_list, fill = TRUE) else NULL
  )
}

.dedml_psock_eval_block <- function(indices) {
  get(".DEDML_PSOCK_STATE", envir = .GlobalEnv)$run_gene_block(indices)
}

.dedml_fork_available <- function() {
  .Platform$OS.type != "windows"
}

.dedml_make_psock_cluster <- function(n_cores, lib_paths = .libPaths()) {
  tryCatch({
    cl <- parallel::makeCluster(n_cores)
    tryCatch(
      parallel::clusterCall(cl, function(paths) {
        .libPaths(paths)
        invisible(.libPaths())
      }, lib_paths),
      error = function(e) {
        parallel::stopCluster(cl)
        stop(e)
      }
    )
    cl
  }, error = function(e) e)
}

.dedml_poisson_irls_psock_block <- function(indices) {
  state <- get(".DEDML_POISSON_IRLS_STATE", envir = .GlobalEnv)
  y_resid_mat <- matrix(NA_real_, nrow = length(indices), ncol = ncol(state$counts))
  rownames(y_resid_mat) <- state$gene_names[indices]
  for (fold in state$fold_plan) {
    y_train_mat <- t(as.matrix(state$counts[indices, fold$train_idx, drop = FALSE]))
    pred <- .fit_predict_glm_poisson_matrix(
      x_train = fold$x_train,
      y_train_mat = y_train_mat,
      x_test = fold$x_test,
      maxit = state$glm_maxit,
      epsilon = state$glm_epsilon
    )
    y_resid_mat[, fold$test_idx] <- t(pred)
  }
  mu_hat_mat <- y_resid_mat
  mu_factor_df <- .summarise_oof_mu_factors_matrix(
    mu_hat_mat = mu_hat_mat,
    cell_type = state$cell_type,
    nb_size = NA_real_,
    propensity = state$propensity
  )
  y_obs_mat <- as.matrix(state$counts[indices, , drop = FALSE])
  y_resid_mat <- (y_obs_mat - pmax(y_resid_mat, 1e-6)) / sqrt(pmax(y_resid_mat, 1e-6))
  rownames(y_resid_mat) <- state$gene_names[indices]
  list(y_resid_mat = y_resid_mat, mu_factor_df = mu_factor_df)
}

#' Fit DEDML on Cell-Level Data
#'
#' Runs donor-blocked cross-fitted DEDML with configurable confounders,
#' treatment learner, outcome model level, and outcome learner. With
#' `outcome_model = "cell"`, the outcome learner estimates nuisance means per
#' cell before the final pseudo-bulk residual regression. With
#' `outcome_model = "pseudobulk"`, the outcome learner estimates nuisance means
#' directly on sample-celltype pseudo-bulk means with cell-count weights, then
#' runs the final residual regression on those same pseudo-bulk units.
#'
#' @param counts Gene-by-cell count matrix (rows are genes, columns are cells).
#' @param meta Cell metadata data.frame.
#' @param donor_id Donor identifier column name (or vector).
#' @param sample_id Sample identifier column name (or vector). If `NULL`, donor id is used.
#' @param treatment Binary treatment indicator column name (or vector), coded 0/1.
#' @param cell_type Cell type column name (or vector).
#' @param confounder_spec Optional list from [dedml_make_confounder_spec()]. If provided,
#' this overrides `treatment_confounders`, `treatment_cell_summaries`, and `outcome_confounders`.
#' @param treatment_confounders Donor-level confounders for treatment model.
#' @param treatment_cell_summaries Cell-level covariates to summarize by donor for treatment model.
#' @param outcome_confounders Cell-level confounders for outcome model.
#' @param cell_types Optional vector of cell types to keep.
#' @param gene_subset Optional vector of gene names to keep.
#' @param n_folds Number of donor folds.
#' @param n_cores Number of cores for gene-level parallelization.
#' @param donor_model Donor treatment model (`"glm"` or `"lightgbm"`).
#' @param outcome_model Outcome nuisance model level (`"cell"` or
#' `"pseudobulk"`). Legacy values `"lightgbm"` and `"glm"` are accepted and
#' interpreted as `outcome_model = "cell"` with the corresponding
#' `outcome_learner`.
#' @param outcome_learner Outcome learner (`"lightgbm"` or `"glm"`). For
#' `outcome_model = "pseudobulk"`, LightGBM is fit to pseudo-bulk means with
#' cell-count weights, while `"glm"` uses weighted GLM/GAM fits.
#' @param outcome_distribution Outcome likelihood/working distribution:
#' `"poisson"` by default, with `"gaussian"` and `"nb"` also available.
#' For `outcome_learner = "glm"`, Poisson and Gaussian use `stats::glm.fit()`
#' and NB uses `mgcv::gam()`
#' with `mgcv::nb()`. For `outcome_learner = "lightgbm"`, these map to
#' LightGBM L2 regression, Poisson, and free-dispersion NB objectives.
#' Gaussian/L2 changes only the nuisance mean learner; its final residual
#' conversion still uses Poisson variance.
#' @param treatment_params Optional list of learner hyperparameters.
#' @param outcome_params Optional outcome learner hyperparameters. The
#' outcome objective/family is controlled by `outcome_distribution`, so
#' redundant `objective`, `metric`, and `family` entries are ignored. For
#' `outcome_model = "cell"`, `outcome_learner = "glm"`, and
#' `outcome_distribution = "poisson"`,
#' `batched_irls = TRUE` enables an exact shared-design IRLS implementation;
#' it uses `n_cores` and `parallel_backend` to split genes into independent
#' blocks. It is opt-in because the regular parallel scalar `glm.fit()` path can
#' be faster for some backend/data-size combinations.
#' @param min_cells_per_sample_ct Minimum cells in each sample-celltype pseudobulk.
#' @param min_samples_per_celltype Minimum samples required per cell type.
#' @param nb_size_min Minimum NB size estimate for residual conversion.
#' @param nb_size_max Maximum NB size estimate for residual conversion; values
#' near this upper bound approximate Poisson residuals.
#' @param pvalue_calibration Optional tail calibration for residual-regression p-values.
#' @param calibration_quantile Chi-square quantile used when `pvalue_calibration = "tail"`.
#' @param parallel_backend Gene-level parallel backend (`"auto"`, `"fork"`, or `"psock"`).
#' With `n_cores > 1`, `"auto"` uses PSOCK workers for portability; `"fork"`
#' remains available as an explicit opt-in on platforms that support it.
#' @param gene_chunk_size Number of genes per parallel task. If `NULL`, DEDML
#' chooses a backend-aware value. This also controls the maximum number of
#' genes materialized together for batched residual regression.
#' @param seed Random seed.
#' @param verbose Print progress.
#' @param donor_id_col Deprecated alias for `donor_id`.
#' @param sample_id_col Deprecated alias for `sample_id`.
#' @param treatment_col Deprecated alias for `treatment`.
#' @param cell_type_col Deprecated alias for `cell_type`.
#' @param focus_celltypes Deprecated alias for `cell_types`.
#' @param treatment_learner Deprecated alias for `donor_model`.
#' @param outcome_learner Outcome nuisance learner (`"lightgbm"` or `"glm"`).
#'
#' @return A list with `results`, `donor_meta`, `settings`, `diagnostics`,
#' and `errors`.
#' The `results` table reports the residual-scale DEDML effect in `estimate`.
#' It also includes `estimate_logfc`, a practical post-hoc log fold-change
#' effect size targeting the untreated-baseline mean. The conversion uses the
#' out-of-fold nuisance mean and treatment propensity; `"gaussian"` and
#' `"poisson"` use the Poisson residual scale, while `"nb"` uses the
#' negative-binomial residual scale. DEDML inference uses the residual-scale
#' statistic and p-value. The `diagnostics` table records design and treatment
#' model issues such as weak fold overlap, numeric non-overlap, and degenerate
#' propensity estimates. Fatal design diagnostics, including separated
#' categorical batch/site/confounder levels observed in only one treatment
#' group, stop the fit before model training.
#' @export
#'
dedml_fit <- function(
    counts,
    meta,
    donor_id = "HATIMID",
    sample_id = "sample_id",
    treatment = "D",
    cell_type = "CellType",
    confounder_spec = NULL,
    treatment_confounders = c("Age", "Sex"),
    treatment_cell_summaries = c("log_nCount_RNA", "nFeature_RNA", "percent.mt"),
    outcome_confounders = c("Age", "Sex", "log_nCount_RNA", "nFeature_RNA", "percent.mt"),
    cell_types = NULL,
    gene_subset = NULL,
    n_folds = 3L,
    n_cores = 1L,
    donor_model = "glm",
    outcome_model = c("cell", "pseudobulk"),
    outcome_learner = c("lightgbm", "glm"),
    outcome_distribution = c("poisson", "gaussian", "nb"),
    treatment_params = list(),
    outcome_params = list(),
    min_cells_per_sample_ct = 10L,
    min_samples_per_celltype = 3L,
    nb_size_min = 0.1,
    nb_size_max = 1e6,
    pvalue_calibration = c("none", "tail"),
    calibration_quantile = 0.9,
    parallel_backend = c("auto", "fork", "psock"),
    gene_chunk_size = NULL,
    seed = 123L,
    verbose = TRUE,
    donor_id_col = NULL,
    sample_id_col = NULL,
    treatment_col = NULL,
    cell_type_col = NULL,
    focus_celltypes = NULL,
    treatment_learner = NULL) {
  outcome_distribution <- match.arg(outcome_distribution)
  pvalue_calibration <- match.arg(pvalue_calibration)
  parallel_backend <- match.arg(parallel_backend)
  calibration_quantile <- as.numeric(calibration_quantile)
  if (!is.finite(calibration_quantile) || calibration_quantile <= 0.5 || calibration_quantile >= 1) {
    stop("calibration_quantile must be a finite number between 0.5 and 1.", call. = FALSE)
  }

  if (!is.matrix(counts) && !inherits(counts, "Matrix")) {
    stop("counts must be a matrix or Matrix object.", call. = FALSE)
  }

  counts <- as.matrix(counts)

  if (ncol(counts) != nrow(meta)) {
    if (!is.null(colnames(counts)) && !is.null(rownames(meta))) {
      shared <- intersect(colnames(counts), rownames(meta))
      if (length(shared) == 0L) {
        stop("Could not align counts columns and meta rows.", call. = FALSE)
      }
      counts <- counts[, shared, drop = FALSE]
      meta <- meta[shared, , drop = FALSE]
    } else {
      stop("counts columns must align with meta rows.", call. = FALSE)
    }
  }

  if (!is.null(donor_id_col)) donor_id <- donor_id_col
  if (!is.null(sample_id_col)) sample_id <- sample_id_col
  if (!is.null(treatment_col)) treatment <- treatment_col
  if (!is.null(cell_type_col)) cell_type <- cell_type_col
  if (!is.null(focus_celltypes)) cell_types <- focus_celltypes
  if (!is.null(treatment_learner)) donor_model <- treatment_learner
  donor_model <- match.arg(donor_model, c("glm", "lightgbm"))

  if (length(outcome_model) > 1L) {
    outcome_model <- outcome_model[[1L]]
  }
  if (length(outcome_learner) > 1L) {
    outcome_learner <- outcome_learner[[1L]]
  }
  if (outcome_model %in% c("lightgbm", "glm")) {
    if (isTRUE(verbose)) {
      message(
        "Interpreting legacy outcome_model='", outcome_model,
        "' as outcome_model='cell', outcome_learner='", outcome_model, "'."
      )
    }
    outcome_learner <- outcome_model
    outcome_model <- "cell"
  }
  outcome_model <- match.arg(outcome_model, c("cell", "pseudobulk"))
  outcome_learner <- match.arg(outcome_learner, c("lightgbm", "glm"))

  outcome_redundant <- intersect(names(outcome_params), c("objective", "metric", "family"))
  if (length(outcome_redundant) > 0L) {
    if (isTRUE(verbose)) {
      message(
        "Ignoring redundant outcome_params entries controlled by outcome_distribution: ",
        paste(outcome_redundant, collapse = ", ")
      )
    }
    outcome_params[outcome_redundant] <- NULL
  }
  outcome_batched_irls <- isTRUE(outcome_params$batched_irls)
  outcome_params$batched_irls <- NULL
  if (isTRUE(outcome_batched_irls) &&
      !(outcome_model == "cell" && outcome_learner == "glm" && outcome_distribution == "poisson") &&
      isTRUE(verbose)) {
    message(
      "Ignoring outcome_params$batched_irls; it is only used for ",
      "outcome_model='cell', outcome_learner='glm', outcome_distribution='poisson'."
    )
    outcome_batched_irls <- FALSE
  }
  outcome_params$outcome_distribution <- outcome_distribution

  meta <- as.data.frame(meta, stringsAsFactors = FALSE)
  if (is.null(sample_id)) {
    sample_id <- donor_id
  }

  donor_res <- .resolve_input_column(meta, donor_id, "donor_id")
  sample_res <- .resolve_input_column(meta, sample_id, "sample_id")
  treatment_res <- .resolve_input_column(meta, treatment, "treatment")
  celltype_res <- .resolve_input_column(meta, cell_type, "cell_type")

  meta$.donor_id <- as.character(donor_res$values)
  meta$.sample_id <- as.character(sample_res$values)
  meta$.D <- as.integer(treatment_res$values)
  meta$.CellType <- as.character(celltype_res$values)

  if (!all(meta$.D %in% c(0L, 1L))) {
    stop("Treatment column must be coded as 0/1.", call. = FALSE)
  }

  if (!is.null(confounder_spec)) {
    dedml_validate_confounder_spec(meta = meta, confounder_spec = confounder_spec)
    treatment_confounders <- confounder_spec$treatment_confounders
    treatment_cell_summaries <- confounder_spec$treatment_cell_summaries
    outcome_confounders <- confounder_spec$outcome_confounders
  }

  n_cores <- as.integer(n_cores)
  if (is.na(n_cores) || n_cores < 1L) {
    stop("n_cores must be a positive integer.", call. = FALSE)
  }
  n_cores <- min(n_cores, parallel::detectCores(logical = TRUE))

  if (!is.null(cell_types)) {
    keep_cells <- meta$.CellType %in% cell_types
    counts <- counts[, keep_cells, drop = FALSE]
    meta <- meta[keep_cells, , drop = FALSE]
  }

  if (!is.null(gene_subset)) {
    missing_genes <- setdiff(gene_subset, rownames(counts))
    if (length(missing_genes) > 0L && isTRUE(verbose)) {
      message("Ignoring missing genes: ", paste(utils::head(missing_genes, 10), collapse = ", "))
    }
    keep_genes <- intersect(gene_subset, rownames(counts))
    counts <- counts[keep_genes, , drop = FALSE]
  }

  if (nrow(meta) == 0L || nrow(counts) == 0L) {
    stop("No data left after filtering.", call. = FALSE)
  }

  needed_cols <- unique(c(
    treatment_confounders,
    treatment_cell_summaries,
    outcome_confounders
  ))
  needed_cols <- setdiff(needed_cols, c("", NA))

  missing_needed <- setdiff(needed_cols, colnames(meta))
  if (length(missing_needed) > 0L) {
    stop("Missing confounder columns in meta: ", paste(missing_needed, collapse = ", "), call. = FALSE)
  }

  keep_complete <- if (length(needed_cols) > 0L) {
    stats::complete.cases(meta[, needed_cols, drop = FALSE])
  } else {
    rep(TRUE, nrow(meta))
  }
  counts <- counts[, keep_complete, drop = FALSE]
  meta <- meta[keep_complete, , drop = FALSE]

  if (ncol(counts) == 0L) {
    stop("No complete-case cells available for analysis.", call. = FALSE)
  }

  donor_consistency <- tapply(meta$.D, meta$.donor_id, function(x) length(unique(x)))
  if (any(donor_consistency != 1L)) {
    bad <- names(donor_consistency)[donor_consistency != 1L]
    stop("Treatment must be constant within donor. Violations: ", paste(utils::head(bad, 10), collapse = ", "), call. = FALSE)
  }

  donor_meta <- dplyr::summarise(
    dplyr::group_by(meta, .donor_id),
    .D = .D[1],
    dplyr::across(dplyr::all_of(treatment_confounders), ~ .x[1]),
    .groups = "drop"
  )

  if (length(treatment_cell_summaries) > 0L) {
    donor_summ <- dplyr::summarise(
      dplyr::group_by(meta, .donor_id),
      dplyr::across(dplyr::all_of(treatment_cell_summaries), ~ mean(as.numeric(.x), na.rm = TRUE), .names = "mean_{.col}"),
      .groups = "drop"
    )
    donor_meta <- dplyr::left_join(donor_meta, donor_summ, by = ".donor_id")
    treatment_features <- c(treatment_confounders, paste0("mean_", treatment_cell_summaries))
  } else {
    treatment_features <- treatment_confounders
  }

  fold_df <- make_stratified_donor_folds(
    donor_id = donor_meta$.donor_id,
    treatment = donor_meta$.D,
    n_folds = n_folds,
    seed = seed
  )

  donor_meta <- dplyr::left_join(donor_meta, fold_df, by = c(".donor_id" = "donor_id"))
  meta <- dplyr::left_join(meta, fold_df, by = c(".donor_id" = "donor_id"))

  x_treat <- .build_design_matrix(donor_meta, treatment_features)
  diagnostics <- .dedml_preflight_diagnostics(
    meta = meta,
    donor_meta = donor_meta,
    treatment_confounders = treatment_confounders,
    treatment_features = treatment_features,
    outcome_confounders = outcome_confounders,
    n_folds = n_folds,
    donor_model = donor_model,
    treatment_params = treatment_params,
    x_treat = x_treat
  )
  .dedml_emit_diagnostic_warnings(diagnostics)
  .dedml_stop_for_fatal_diagnostics(diagnostics)

  donor_meta$propensity_oof <- NA_real_

  for (k in seq_len(n_folds)) {
    train_idx <- which(donor_meta$fold != k)
    test_idx <- which(donor_meta$fold == k)

    donor_meta$propensity_oof[test_idx] <- .fit_predict_treatment(
      x_train = x_treat[train_idx, , drop = FALSE],
      y_train = donor_meta$.D[train_idx],
      x_test = x_treat[test_idx, , drop = FALSE],
      learner = donor_model,
      params = treatment_params
    )
  }

  donor_meta$D_resid <- donor_meta$.D - donor_meta$propensity_oof
  propensity_diagnostics <- .dedml_propensity_diagnostics(donor_meta)
  if (nrow(propensity_diagnostics) > 0L) {
    diagnostics <- as.data.frame(data.table::rbindlist(
      list(diagnostics, propensity_diagnostics),
      fill = TRUE
    ))
    .dedml_emit_diagnostic_warnings(propensity_diagnostics)
  }

  meta <- dplyr::left_join(
    meta,
    donor_meta[, c(".donor_id", "D_resid", "propensity_oof"), drop = FALSE],
    by = ".donor_id"
  )

  if (donor_model == "lightgbm" && outcome_learner == "lightgbm" && n_cores > 1L &&
      "lightgbm" %in% loadedNamespaces()) {
    gc()
    try(unloadNamespace("lightgbm"), silent = TRUE)
    gc()
  }

  x_outcome <- .build_design_matrix(meta, outcome_confounders)
  gene_names <- rownames(counts)

  if (is.null(gene_names)) {
    gene_names <- paste0("gene_", seq_len(nrow(counts)))
    rownames(counts) <- gene_names
  }

  fold_plan <- lapply(seq_len(n_folds), function(k) {
    train_idx <- which(meta$fold != k)
    test_idx <- which(meta$fold == k)
    list(
      train_idx = train_idx,
      test_idx = test_idx,
      x_train = x_outcome[train_idx, , drop = FALSE],
      x_test = x_outcome[test_idx, , drop = FALSE]
    )
  })
  if (outcome_model == "cell" && outcome_learner == "glm" && outcome_distribution == "nb") {
    fold_plan <- lapply(fold_plan, function(fold) {
      x_train_df <- as.data.frame(fold$x_train, check.names = FALSE)
      x_test_df <- as.data.frame(fold$x_test, check.names = FALSE)
      covar_names <- make.names(colnames(x_train_df), unique = TRUE)
      colnames(x_train_df) <- covar_names
      colnames(x_test_df) <- covar_names
      form <- if (length(covar_names) == 0L) {
        stats::as.formula(".y ~ 1")
      } else {
        stats::as.formula(paste(".y ~", paste(covar_names, collapse = " + ")))
      }
      fold$nb_params <- utils::modifyList(
        outcome_params,
        list(
          .nb_train_df_base = data.frame(x_train_df, check.names = FALSE),
          .nb_test_df = data.frame(x_test_df, check.names = FALSE),
          .nb_form = form
        )
      )
      fold
    })
  }

  run_one_gene <- function(g_idx, gene_fold_plan = fold_plan) {
    y_obs <- as.numeric(counts[g_idx, ])
    mu_hat <- rep(NA_real_, length(y_obs))
    nb_sizes <- numeric(0)
    for (fold in gene_fold_plan) {
      fold_pred <- .fit_predict_outcome(
        x_train = fold$x_train,
        y_train = y_obs[fold$train_idx],
        x_test = fold$x_test,
        learner = outcome_learner,
        params = if (!is.null(fold$params)) {
          fold$params
        } else if (!is.null(fold$nb_params)) {
          fold$nb_params
        } else {
          outcome_params
        }
      )
      mu_hat[fold$test_idx] <- fold_pred
      nb_sizes <- c(nb_sizes, attr(fold_pred, "nb_size"))
    }
    lgb_nb_size <- stats::median(nb_sizes[is.finite(nb_sizes) & nb_sizes > 0], na.rm = TRUE)
    if (!is.finite(lgb_nb_size)) {
      lgb_nb_size <- NA_real_
    }
    mu_summary <- .summarise_oof_mu_factors(
      mu = mu_hat,
      cell_type = meta$.CellType,
      nb_size = lgb_nb_size,
      propensity = meta$propensity_oof
    )
    mu_summary$gene <- gene_names[g_idx]

    y_resid <- .compute_outcome_residual(
      y = y_obs,
      mu = mu_hat,
      cell_type = meta$.CellType,
      outcome_distribution = outcome_distribution,
      nb_size = lgb_nb_size,
      nb_size_min = nb_size_min,
      nb_size_max = nb_size_max
    )
    list(
      .resid = TRUE,
      gene_idx = g_idx,
      gene = gene_names[g_idx],
      y_resid = y_resid,
      nb_size = attr(y_resid, "nb_size"),
      mu_summary = mu_summary
    )
  }

  run_residual_batch <- function(resid_list) {
    if (length(resid_list) == 0L) {
      return(NULL)
    }

    ord <- order(vapply(resid_list, function(x) x$gene_idx, integer(1L)))
    resid_list <- resid_list[ord]
    resid_chunks <- split(
      seq_along(resid_list),
      ceiling(seq_along(resid_list) / gene_chunk_size)
    )
    result_chunks <- lapply(resid_chunks, function(idx) {
      block <- resid_list[idx]
      y_resid_mat <- do.call(rbind, lapply(block, function(x) x$y_resid))
      batch_gene_names <- vapply(block, function(x) x$gene, character(1L))
      rownames(y_resid_mat) <- batch_gene_names

      fit_residual_regression_pseudobulk_matrix(
        y_resid_mat = y_resid_mat,
        d_resid = meta$D_resid,
        cell_type = meta$.CellType,
        sample_id = meta$.sample_id,
        donor_id = meta$.donor_id,
        gene_names = batch_gene_names,
        min_cells_per_sample_ct = min_cells_per_sample_ct,
        min_samples = min_samples_per_celltype,
        cluster_robust = TRUE
      )
    })
    result_chunks <- Filter(Negate(is.null), result_chunks)
    if (length(result_chunks) == 0L) {
      return(NULL)
    }
    as.data.frame(data.table::rbindlist(result_chunks, fill = TRUE))
  }

  if (isTRUE(verbose)) {
    message("Running DEDML on ", length(gene_names), " genes with ", n_cores, " core(s)")
  }

  requested_parallel_backend <- parallel_backend
  if (parallel_backend == "auto") {
    parallel_backend <- if (n_cores > 1L) "psock" else "fork"
  }
  if (parallel_backend == "fork" && !.dedml_fork_available()) {
    if (isTRUE(verbose)) {
      message("parallel_backend='fork' is unavailable on this platform; using 'psock'.")
    }
    parallel_backend <- "psock"
  }

  if (is.null(gene_chunk_size)) {
    gene_chunk_size <- if (n_cores > 1L) {
      tasks_per_worker <- if (parallel_backend == "psock") 1L else 4L
      max(1L, ceiling(length(gene_names) / (n_cores * tasks_per_worker)))
    } else {
      length(gene_names)
    }
  }
  gene_chunk_size <- as.integer(gene_chunk_size)
  if (is.na(gene_chunk_size) || gene_chunk_size < 1L) {
    stop("gene_chunk_size must be NULL or a positive integer.", call. = FALSE)
  }
  gene_chunks <- split(
    seq_along(gene_names),
    ceiling(seq_along(gene_names) / gene_chunk_size)
  )

  if (outcome_model == "pseudobulk") {
    if (isTRUE(verbose)) {
      message("Using pseudo-bulk outcome nuisance fit")
    }
    pb_data <- .dedml_make_pseudobulk_outcome_data(
      counts = counts,
      meta = meta,
      outcome_confounders = outcome_confounders,
      min_cells_per_sample_ct = min_cells_per_sample_ct
    )
    if (is.null(pb_data)) {
      stop("No valid pseudo-bulk units are available for outcome_model='pseudobulk'.", call. = FALSE)
    }
    pb_df <- pb_data$pb_df
    y_pb_mat <- pb_data$y_pb_mat
    x_pb <- pb_data$x_pb
    pb_fold_plan <- lapply(seq_len(n_folds), function(k) {
      train_idx <- which(pb_df$fold != k)
      test_idx <- which(pb_df$fold == k)
      list(
        train_idx = train_idx,
        test_idx = test_idx,
        x_train = x_pb[train_idx, , drop = FALSE],
        x_test = x_pb[test_idx, , drop = FALSE],
        w_train = pb_df$n_cells[train_idx]
      )
    })

    run_pseudobulk_block <- function(indices) {
      tryCatch({
        block_gene_names <- gene_names[indices]
        mu_hat_mat <- matrix(NA_real_, nrow = length(indices), ncol = ncol(y_pb_mat))
        y_resid_pb_mat <- matrix(NA_real_, nrow = length(indices), ncol = ncol(y_pb_mat))
        rownames(mu_hat_mat) <- block_gene_names
        rownames(y_resid_pb_mat) <- block_gene_names
        nb_size_vec <- rep(NA_real_, length(indices))

        for (j in seq_along(indices)) {
          g_idx <- indices[j]
          y_obs <- as.numeric(y_pb_mat[g_idx, ])
          nb_sizes <- numeric(0)
          for (fold in pb_fold_plan) {
            fold_pred <- .fit_predict_pseudobulk_outcome(
              x_train = fold$x_train,
              y_train = y_obs[fold$train_idx],
              x_test = fold$x_test,
              weights_train = fold$w_train,
              learner = outcome_learner,
              distribution = outcome_distribution,
              params = outcome_params
            )
            mu_hat_mat[j, fold$test_idx] <- fold_pred
            nb_sizes <- c(nb_sizes, attr(fold_pred, "nb_size"))
          }
          nb_size <- stats::median(nb_sizes[is.finite(nb_sizes) & nb_sizes > 0], na.rm = TRUE)
          if (!is.finite(nb_size)) {
            nb_size <- NA_real_
          }
          y_resid <- .compute_outcome_residual(
            y = y_obs,
            mu = mu_hat_mat[j, ],
            cell_type = pb_df$CellType,
            outcome_distribution = outcome_distribution,
            nb_size = nb_size,
            nb_size_min = nb_size_min,
            nb_size_max = nb_size_max
          )
          y_resid_pb_mat[j, ] <- y_resid
          nb_size_vec[j] <- attr(y_resid, "nb_size")
        }

        rr <- fit_residual_regression_pseudobulk_matrix_from_pb(
          y_resid_pb_mat = y_resid_pb_mat,
          pb_df = pb_df,
          gene_names = block_gene_names,
          min_samples = min_samples_per_celltype,
          cluster_robust = TRUE
        )
        if (is.null(rr) || nrow(rr) == 0L) {
          return(NULL)
        }
        nb_size_df <- data.frame(gene = block_gene_names, nb_size = nb_size_vec, stringsAsFactors = FALSE)
        rr <- dplyr::left_join(rr, nb_size_df, by = "gene")
        mu_factor_df <- .summarise_oof_mu_factors_matrix(
          mu_hat_mat = mu_hat_mat,
          cell_type = pb_df$CellType,
          nb_size = nb_size_vec,
          propensity = pb_df$propensity_oof,
          weights = pb_df$n_cells
        )
        .append_logfc_effect_sizes(rr, mu_factor_df, outcome_distribution = outcome_distribution)
      }, error = function(e) {
        list(.error = TRUE, gene = paste(gene_names[indices], collapse = ","), message = conditionMessage(e))
      })
    }

    helper_names <- c(
      ".fit_predict_pseudobulk_outcome", ".fit_predict_lightgbm", ".compute_outcome_residual",
      ".summarise_oof_mu_factors", ".summarise_oof_mu_factors_matrix",
      ".dedml_logfc_from_residual", ".dedml_logfc_from_counterfactual_residual",
      ".append_logfc_effect_sizes", "fit_residual_regression_pseudobulk_matrix_from_pb"
    )
    helper_env <- environment(.fit_predict_pseudobulk_outcome)
    for (nm in helper_names[helper_names %in% ls(envir = helper_env, all.names = TRUE)]) {
      obj <- get(nm, envir = helper_env)
      if (is.function(obj)) {
        environment(obj) <- environment(run_pseudobulk_block)
      }
      assign(nm, obj, envir = environment(run_pseudobulk_block))
    }

    if (n_cores > 1L) {
      if (parallel_backend == "psock") {
        cl <- .dedml_make_psock_cluster(n_cores)
        if (inherits(cl, "error")) {
          if (isTRUE(verbose)) {
            message("PSOCK backend failed (", conditionMessage(cl), "); running pseudo-bulk outcome serially.")
          }
          raw_result_chunks <- lapply(gene_chunks, run_pseudobulk_block)
        } else {
          on.exit(parallel::stopCluster(cl), add = TRUE)
          psock_env <- new.env(parent = emptyenv())
          psock_env$.DEDML_PSOCK_STATE <- list(run_gene_block = run_pseudobulk_block)
          psock_env$.dedml_psock_eval_block <- .dedml_psock_eval_block
          parallel::clusterExport(
            cl,
            varlist = c(".DEDML_PSOCK_STATE", ".dedml_psock_eval_block"),
            envir = psock_env
          )
          raw_result_chunks <- parallel::parLapply(cl, gene_chunks, .dedml_psock_eval_block)
          parallel::clusterEvalQ(cl, rm(.DEDML_PSOCK_STATE, .dedml_psock_eval_block, envir = .GlobalEnv))
        }
      } else {
        raw_result_chunks <- parallel::mclapply(gene_chunks, run_pseudobulk_block, mc.cores = n_cores)
      }
    } else {
      raw_result_chunks <- lapply(gene_chunks, run_pseudobulk_block)
    }

    error_list <- Filter(function(x) is.list(x) && !is.data.frame(x) && isTRUE(x$.error), raw_result_chunks)
    result_chunks <- Filter(function(x) is.data.frame(x) && nrow(x) > 0L, raw_result_chunks)
    settings <- list(
      donor_id = donor_res$source,
      sample_id = sample_res$source,
      treatment = treatment_res$source,
      cell_type = celltype_res$source,
      treatment_confounders = treatment_confounders,
      treatment_cell_summaries = treatment_cell_summaries,
      outcome_confounders = outcome_confounders,
      n_folds = n_folds,
      donor_model = donor_model,
      outcome_model = outcome_model,
      outcome_learner = outcome_learner,
      outcome_distribution = outcome_distribution,
      nb_size_min = nb_size_min,
      nb_size_max = nb_size_max,
      pvalue_calibration = pvalue_calibration,
      calibration_quantile = calibration_quantile,
      outcome_batched_irls = outcome_batched_irls,
      requested_parallel_backend = requested_parallel_backend,
      parallel_backend = parallel_backend,
      gene_chunk_size = gene_chunk_size,
      n_gene_chunks = length(gene_chunks),
      outcome_fast_path = "pseudobulk_outcome_nuisance",
      residual_regression_memory = "pseudobulk",
      min_cells_per_sample_ct = min_cells_per_sample_ct,
      min_samples_per_celltype = min_samples_per_celltype,
      seed = seed
    )
    return(.finalize_dedml_output(
      result_list = result_chunks,
      error_list = error_list,
      donor_meta = donor_meta,
      settings = settings,
      diagnostics = diagnostics,
      pvalue_calibration = pvalue_calibration,
      calibration_quantile = calibration_quantile
    ))
  }

  if (outcome_model == "cell" && outcome_learner == "glm" && outcome_distribution == "gaussian") {
    if (isTRUE(verbose)) {
      message("Using batched Gaussian GLM outcome fit")
    }
    gaussian_fold_plan <- lapply(fold_plan, function(fold) {
      fold$x_test_i <- cbind("(Intercept)" = 1, fold$x_test)
      fold$fit_qr <- qr(cbind("(Intercept)" = 1, fold$x_train))
      fold
    })
    result_chunks <- lapply(gene_chunks, function(indices) {
      block_gene_names <- gene_names[indices]
      mu_hat_mat <- matrix(NA_real_, nrow = length(indices), ncol = ncol(counts))
      rownames(mu_hat_mat) <- block_gene_names
      for (fold in gaussian_fold_plan) {
        y_train_mat <- t(as.matrix(counts[indices, fold$train_idx, drop = FALSE]))
        beta <- qr.coef(fold$fit_qr, y_train_mat)
        beta[is.na(beta)] <- 0
        mu_hat_mat[, fold$test_idx] <- t(fold$x_test_i %*% beta)
      }
      y_obs_mat <- as.matrix(counts[indices, , drop = FALSE])
      y_resid_mat <- (y_obs_mat - pmax(mu_hat_mat, 1e-6)) / sqrt(pmax(mu_hat_mat, 1e-6))
      rownames(y_resid_mat) <- block_gene_names
      mu_factor_df <- .summarise_oof_mu_factors_matrix(
        mu_hat_mat = mu_hat_mat,
        cell_type = meta$.CellType,
        nb_size = NA_real_,
        propensity = meta$propensity_oof
      )
      rr <- fit_residual_regression_pseudobulk_matrix(
        y_resid_mat = y_resid_mat,
        d_resid = meta$D_resid,
        cell_type = meta$.CellType,
        sample_id = meta$.sample_id,
        donor_id = meta$.donor_id,
        gene_names = block_gene_names,
        min_cells_per_sample_ct = min_cells_per_sample_ct,
        min_samples = min_samples_per_celltype,
        cluster_robust = TRUE
      )
      .append_logfc_effect_sizes(rr, mu_factor_df, outcome_distribution = outcome_distribution)
    })
    result_chunks <- Filter(Negate(is.null), result_chunks)
    result_df_batch <- if (length(result_chunks) > 0L) {
      as.data.frame(data.table::rbindlist(result_chunks, fill = TRUE))
    } else {
      NULL
    }
    if (is.data.frame(result_df_batch) && nrow(result_df_batch) > 0L) {
      result_df_batch$nb_size <- NA_real_
    }
    settings <- list(
      donor_id = donor_res$source,
      sample_id = sample_res$source,
      treatment = treatment_res$source,
      cell_type = celltype_res$source,
      treatment_confounders = treatment_confounders,
      treatment_cell_summaries = treatment_cell_summaries,
      outcome_confounders = outcome_confounders,
      n_folds = n_folds,
      donor_model = donor_model,
      outcome_model = outcome_model,
      outcome_learner = outcome_learner,
      outcome_distribution = outcome_distribution,
      nb_size_min = nb_size_min,
      nb_size_max = nb_size_max,
      pvalue_calibration = pvalue_calibration,
      calibration_quantile = calibration_quantile,
      outcome_batched_irls = outcome_batched_irls,
      requested_parallel_backend = requested_parallel_backend,
      parallel_backend = parallel_backend,
      gene_chunk_size = gene_chunk_size,
      n_gene_chunks = length(gene_chunks),
      outcome_fast_path = "batched_gaussian_glm",
      residual_regression_memory = "chunked",
      min_cells_per_sample_ct = min_cells_per_sample_ct,
      min_samples_per_celltype = min_samples_per_celltype,
      seed = seed
    )
    return(.finalize_dedml_output(
      result_list = if (is.data.frame(result_df_batch) && nrow(result_df_batch) > 0L) list(result_df_batch) else list(),
      error_list = list(),
      donor_meta = donor_meta,
      settings = settings,
      diagnostics = diagnostics,
      pvalue_calibration = pvalue_calibration,
      calibration_quantile = calibration_quantile
    ))
  }

  if (outcome_model == "cell" && outcome_learner == "glm" && outcome_distribution == "poisson" && isTRUE(outcome_batched_irls)) {
    if (isTRUE(verbose)) {
      message("Using batched Poisson GLM outcome IRLS")
    }
    glm_maxit <- if (!is.null(outcome_params$maxit)) as.integer(outcome_params$maxit) else 25L
    glm_epsilon <- if (!is.null(outcome_params$epsilon)) as.numeric(outcome_params$epsilon) else 1e-8
    poisson_chunks <- gene_chunks
    poisson_resid_to_results <- function(resid_chunks) {
      result_chunks <- lapply(resid_chunks, function(chunk) {
        y_resid_mat <- if (is.list(chunk) && !is.null(chunk$y_resid_mat)) chunk$y_resid_mat else chunk
        fit_residual_regression_pseudobulk_matrix(
          y_resid_mat = y_resid_mat,
          d_resid = meta$D_resid,
          cell_type = meta$.CellType,
          sample_id = meta$.sample_id,
          donor_id = meta$.donor_id,
          gene_names = rownames(y_resid_mat),
          min_cells_per_sample_ct = min_cells_per_sample_ct,
          min_samples = min_samples_per_celltype,
          cluster_robust = TRUE
        )
      })
      result_chunks <- Filter(Negate(is.null), result_chunks)
      if (length(result_chunks) == 0L) {
        return(NULL)
      }
      result_df <- as.data.frame(data.table::rbindlist(result_chunks, fill = TRUE))
      mu_factor_df <- data.table::rbindlist(lapply(resid_chunks, function(chunk) {
        if (is.list(chunk) && !is.null(chunk$mu_factor_df)) chunk$mu_factor_df else NULL
      }), fill = TRUE)
      .append_logfc_effect_sizes(result_df, mu_factor_df, outcome_distribution = outcome_distribution)
    }
    if (n_cores > 1L) {
      poisson_state <- list(
        counts = counts,
        gene_names = gene_names,
        fold_plan = fold_plan,
        cell_type = meta$.CellType,
        propensity = meta$propensity_oof,
        glm_maxit = glm_maxit,
        glm_epsilon = glm_epsilon
      )
      if (parallel_backend == "psock") {
        cl <- .dedml_make_psock_cluster(n_cores)
        if (inherits(cl, "error")) {
          if (isTRUE(verbose)) {
            message(
              "PSOCK backend failed (", conditionMessage(cl),
              "); running batched Poisson IRLS serially."
            )
          }
          assign(".DEDML_POISSON_IRLS_STATE", poisson_state, envir = .GlobalEnv)
          on.exit(rm(.DEDML_POISSON_IRLS_STATE, envir = .GlobalEnv), add = TRUE)
          resid_chunks <- lapply(poisson_chunks, .dedml_poisson_irls_psock_block)
        } else {
          on.exit(parallel::stopCluster(cl), add = TRUE)
          psock_env <- new.env(parent = emptyenv())
          psock_env$.DEDML_POISSON_IRLS_STATE <- poisson_state
          psock_env$.fit_glm_poisson_coef_one <- .fit_glm_poisson_coef_one
          psock_env$.fit_predict_glm_poisson <- .fit_predict_glm_poisson
          psock_env$.fit_predict_glm_poisson_matrix <- .fit_predict_glm_poisson_matrix
          psock_env$.summarise_oof_mu_factors <- .summarise_oof_mu_factors
          psock_env$.summarise_oof_mu_factors_matrix <- .summarise_oof_mu_factors_matrix
          psock_env$.dedml_poisson_irls_psock_block <- .dedml_poisson_irls_psock_block
          parallel::clusterExport(
            cl,
            varlist = c(
              ".DEDML_POISSON_IRLS_STATE",
              ".fit_glm_poisson_coef_one",
              ".fit_predict_glm_poisson",
              ".fit_predict_glm_poisson_matrix",
              ".summarise_oof_mu_factors",
              ".summarise_oof_mu_factors_matrix",
              ".dedml_poisson_irls_psock_block"
            ),
            envir = psock_env
          )
          resid_chunks <- parallel::parLapply(cl, poisson_chunks, .dedml_poisson_irls_psock_block)
          parallel::clusterEvalQ(
            cl,
            rm(
              .DEDML_POISSON_IRLS_STATE,
              .fit_glm_poisson_coef_one,
              .fit_predict_glm_poisson,
              .fit_predict_glm_poisson_matrix,
              .summarise_oof_mu_factors,
              .summarise_oof_mu_factors_matrix,
              .dedml_poisson_irls_psock_block,
              envir = .GlobalEnv
            )
          )
        }
      } else {
        assign(".DEDML_POISSON_IRLS_STATE", poisson_state, envir = .GlobalEnv)
        on.exit(rm(.DEDML_POISSON_IRLS_STATE, envir = .GlobalEnv), add = TRUE)
        resid_chunks <- parallel::mclapply(
          poisson_chunks,
          .dedml_poisson_irls_psock_block,
          mc.cores = n_cores
        )
      }
      result_df_batch <- poisson_resid_to_results(resid_chunks)
    } else {
      resid_chunks <- lapply(poisson_chunks, function(indices) {
        block_gene_names <- gene_names[indices]
        mu_hat_mat <- matrix(NA_real_, nrow = length(indices), ncol = ncol(counts))
        rownames(mu_hat_mat) <- block_gene_names
        for (fold in fold_plan) {
          y_train_mat <- t(as.matrix(counts[indices, fold$train_idx, drop = FALSE]))
          pred <- .fit_predict_glm_poisson_matrix(
            x_train = fold$x_train,
            y_train_mat = y_train_mat,
            x_test = fold$x_test,
            maxit = glm_maxit,
            epsilon = glm_epsilon
          )
          mu_hat_mat[, fold$test_idx] <- t(pred)
        }
        y_obs_mat <- as.matrix(counts[indices, , drop = FALSE])
        y_resid_mat <- (y_obs_mat - pmax(mu_hat_mat, 1e-6)) / sqrt(pmax(mu_hat_mat, 1e-6))
        rownames(y_resid_mat) <- block_gene_names
        mu_factor_df <- .summarise_oof_mu_factors_matrix(
          mu_hat_mat = mu_hat_mat,
          cell_type = meta$.CellType,
          nb_size = NA_real_,
          propensity = meta$propensity_oof
        )
        list(y_resid_mat = y_resid_mat, mu_factor_df = mu_factor_df)
      })
      result_df_batch <- poisson_resid_to_results(resid_chunks)
    }
    if (is.data.frame(result_df_batch) && nrow(result_df_batch) > 0L) {
      result_df_batch$nb_size <- NA_real_
    }
    settings <- list(
      donor_id = donor_res$source,
      sample_id = sample_res$source,
      treatment = treatment_res$source,
      cell_type = celltype_res$source,
      treatment_confounders = treatment_confounders,
      treatment_cell_summaries = treatment_cell_summaries,
      outcome_confounders = outcome_confounders,
      n_folds = n_folds,
      donor_model = donor_model,
      outcome_model = outcome_model,
      outcome_learner = outcome_learner,
      outcome_distribution = outcome_distribution,
      nb_size_min = nb_size_min,
      nb_size_max = nb_size_max,
      pvalue_calibration = pvalue_calibration,
      calibration_quantile = calibration_quantile,
      outcome_batched_irls = outcome_batched_irls,
      requested_parallel_backend = requested_parallel_backend,
      parallel_backend = parallel_backend,
      gene_chunk_size = gene_chunk_size,
      n_gene_chunks = length(gene_chunks),
      outcome_fast_path = "batched_poisson_glm_irls",
      residual_regression_memory = "chunked",
      min_cells_per_sample_ct = min_cells_per_sample_ct,
      min_samples_per_celltype = min_samples_per_celltype,
      seed = seed
    )
    return(.finalize_dedml_output(
      result_list = if (is.data.frame(result_df_batch) && nrow(result_df_batch) > 0L) list(result_df_batch) else list(),
      error_list = list(),
      donor_meta = donor_meta,
      settings = settings,
      diagnostics = diagnostics,
      pvalue_calibration = pvalue_calibration,
      calibration_quantile = calibration_quantile
    ))
  }

  worker_fun <- function(i, gene_fold_plan = fold_plan) {
    tryCatch(
      run_one_gene(i, gene_fold_plan = gene_fold_plan),
      error = function(e) list(.error = TRUE, gene = gene_names[i], message = conditionMessage(e))
    )
  }

  run_gene_block <- function(indices) {
    block_fold_plan <- fold_plan
    if (outcome_model == "cell" && outcome_learner == "lightgbm" && !isFALSE(outcome_params$.reuse_lgb_dataset)) {
      block_fold_plan <- .add_lgb_datasets_to_fold_plan(fold_plan, outcome_params)
      on.exit({
        .finalize_lgb_datasets_in_fold_plan(block_fold_plan)
        gc(FALSE)
      }, add = TRUE)
    }
    lapply(indices, worker_fun, gene_fold_plan = block_fold_plan)
  }

  helper_names <- c(
    ".fit_predict_outcome", ".fit_predict_lightgbm",
    ".fit_predict_glm_gaussian", ".fit_predict_glm_poisson", ".fit_predict_gam_nb",
    ".add_lgb_datasets_to_fold_plan", ".finalize_lgb_datasets_in_fold_plan",
    ".compute_outcome_residual",
    ".summarise_oof_mu_factors", ".summarise_oof_mu_factors_matrix",
    ".dedml_logfc_from_residual", ".dedml_logfc_from_counterfactual_residual",
    ".append_logfc_effect_sizes",
    "fit_residual_regression_pseudobulk",
    "fit_residual_regression_pseudobulk_matrix"
  )
  helper_env <- environment(.fit_predict_outcome)
  for (nm in helper_names[helper_names %in% ls(envir = helper_env, all.names = TRUE)]) {
    obj <- get(nm, envir = helper_env)
    if (is.function(obj)) {
      environment(obj) <- environment(run_one_gene)
    }
    assign(nm, obj, envir = environment(run_one_gene))
  }

  if (n_cores > 1L) {
    if (parallel_backend == "psock") {
      cl <- .dedml_make_psock_cluster(n_cores)
      if (inherits(cl, "error")) {
        if (isTRUE(verbose)) {
          message("PSOCK backend failed (", conditionMessage(cl), "); running serially.")
        }
        raw_result_chunks <- lapply(gene_chunks, run_gene_block)
      } else {
        on.exit(parallel::stopCluster(cl), add = TRUE)
        psock_env <- new.env(parent = emptyenv())
        psock_env$.DEDML_PSOCK_STATE <- list(run_gene_block = run_gene_block)
        psock_env$.dedml_psock_eval_block <- .dedml_psock_eval_block
        parallel::clusterExport(
          cl,
          varlist = c(".DEDML_PSOCK_STATE", ".dedml_psock_eval_block"),
          envir = psock_env
        )
        raw_result_chunks <- parallel::parLapply(cl, gene_chunks, .dedml_psock_eval_block)
        parallel::clusterEvalQ(cl, rm(.DEDML_PSOCK_STATE, .dedml_psock_eval_block, envir = .GlobalEnv))
      }
    } else {
      raw_result_chunks <- parallel::mclapply(
        X = gene_chunks,
        FUN = run_gene_block,
        mc.cores = n_cores
      )
    }
  } else {
    raw_result_chunks <- lapply(gene_chunks, run_gene_block)
  }
  raw_results <- unlist(raw_result_chunks, recursive = FALSE, use.names = FALSE)

  error_list <- Filter(function(x) is.list(x) && !is.data.frame(x) && isTRUE(x$.error), raw_results)
  resid_list <- Filter(function(x) is.list(x) && !isTRUE(x$.error) && isTRUE(x$.resid), raw_results)
  result_df_batch <- run_residual_batch(resid_list)
  if (is.data.frame(result_df_batch) && nrow(result_df_batch) > 0L) {
    nb_size_df <- data.table::rbindlist(lapply(resid_list, function(x) {
      if (is.null(x$nb_size) || all(is.na(x$nb_size))) {
        return(data.frame(gene = x$gene, nb_size = NA_real_, stringsAsFactors = FALSE))
      }
      data.frame(gene = x$gene, nb_size = as.numeric(x$nb_size), stringsAsFactors = FALSE)
    }), fill = TRUE)
    if (nrow(nb_size_df) > 0L) {
      result_df_batch <- dplyr::left_join(result_df_batch, nb_size_df, by = "gene")
    }
    mu_factor_df <- data.table::rbindlist(lapply(resid_list, function(x) x$mu_summary), fill = TRUE)
    if (nrow(mu_factor_df) > 0L) {
      result_df_batch <- .append_logfc_effect_sizes(
        result_df_batch,
        mu_factor_df,
        outcome_distribution = outcome_distribution
      )
    }
  }
  result_list <- if (is.data.frame(result_df_batch) && nrow(result_df_batch) > 0L) {
    list(result_df_batch)
  } else {
    list()
  }

  settings <- list(
    donor_id = donor_res$source,
    sample_id = sample_res$source,
    treatment = treatment_res$source,
    cell_type = celltype_res$source,
    treatment_confounders = treatment_confounders,
    treatment_cell_summaries = treatment_cell_summaries,
    outcome_confounders = outcome_confounders,
    n_folds = n_folds,
    donor_model = donor_model,
    outcome_model = outcome_model,
    outcome_learner = outcome_learner,
    outcome_distribution = outcome_distribution,
    nb_size_min = nb_size_min,
    nb_size_max = nb_size_max,
    pvalue_calibration = pvalue_calibration,
    calibration_quantile = calibration_quantile,
    outcome_batched_irls = outcome_batched_irls,
    requested_parallel_backend = requested_parallel_backend,
    parallel_backend = parallel_backend,
    gene_chunk_size = gene_chunk_size,
    n_gene_chunks = length(gene_chunks),
    residual_regression_memory = "chunked",
    min_cells_per_sample_ct = min_cells_per_sample_ct,
    min_samples_per_celltype = min_samples_per_celltype,
    seed = seed
  )

  .finalize_dedml_output(
    result_list = result_list,
    error_list = error_list,
    donor_meta = donor_meta,
    settings = settings,
    diagnostics = diagnostics,
    pvalue_calibration = pvalue_calibration,
    calibration_quantile = calibration_quantile
  )
}

#' Fit DEDML Directly from a Seurat Object
#'
#' Convenience wrapper that extracts counts and metadata from a Seurat object
#' and calls [dedml_fit()].
#'
#' @param object Seurat object.
#' @param assay Assay to use.
#' @param layer Layer to use for counts.
#' @param donor_id Donor identifier column in `object@meta.data`.
#' @param sample_id Sample identifier column. If `NULL`, donor id is used.
#' @param treatment Binary treatment indicator column (0/1).
#' @param cell_type Cell type column.
#' @param ... Additional arguments passed to [dedml_fit()].
#'
#' @return A list with `results`, `donor_meta`, `settings`, `diagnostics`,
#' and `errors`.
#' @export
#'
dedml_fit_seurat <- function(
    object,
    assay = "RNA",
    layer = "counts",
    donor_id = "HATIMID",
    sample_id = NULL,
    treatment = "D",
    cell_type = "celltype_manual",
    ...) {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Package 'Seurat' is required for dedml_fit_seurat().", call. = FALSE)
  }

  meta <- object@meta.data

  if (!(donor_id %in% colnames(meta))) {
    stop("donor_id not found in Seurat metadata: ", donor_id, call. = FALSE)
  }

  if (is.null(sample_id)) {
    sample_id <- donor_id
  }

  if (!(sample_id %in% colnames(meta))) {
    stop("sample_id not found in Seurat metadata: ", sample_id, call. = FALSE)
  }

  required <- c(treatment, cell_type)
  missing_required <- setdiff(required, colnames(meta))
  if (length(missing_required) > 0L) {
    stop("Missing required Seurat metadata columns: ", paste(missing_required, collapse = ", "), call. = FALSE)
  }

  counts <- Seurat::GetAssayData(object, assay = assay, layer = layer)

  dedml_fit(
    counts = counts,
    meta = meta,
    donor_id = donor_id,
    sample_id = sample_id,
    treatment = treatment,
    cell_type = cell_type,
    ...
  )
}
