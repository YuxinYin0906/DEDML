available_dedml_learners <- function() {
  list(
    treatment = c("glm", "lightgbm"),
    outcome = c("lightgbm", "glm")
  )
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

.fit_predict_glm_binomial <- function(x_train, y_train, x_test) {
  x_train_i <- cbind("(Intercept)" = 1, x_train)
  fit <- suppressWarnings(stats::glm.fit(x = x_train_i, y = y_train, family = stats::binomial()))

  beta <- fit$coefficients
  beta[is.na(beta)] <- 0

  pred <- stats::plogis(cbind("(Intercept)" = 1, x_test) %*% beta)
  pmin(pmax(as.numeric(pred), 1e-4), 1 - 1e-4)
}

.fit_predict_glm_poisson <- function(x_train, y_train, x_test) {
  x_train_i <- cbind("(Intercept)" = 1, x_train)
  fit <- suppressWarnings(stats::glm.fit(x = x_train_i, y = y_train, family = stats::poisson()))

  beta <- fit$coefficients
  beta[is.na(beta)] <- 0

  pred <- exp(cbind("(Intercept)" = 1, x_test) %*% beta)
  pmax(as.numeric(pred), 1e-6)
}

.fit_predict_lightgbm <- function(x_train, y_train, x_test, objective, params = list()) {
  if (!requireNamespace("lightgbm", quietly = TRUE)) {
    stop("Package 'lightgbm' is required for learner='lightgbm'.", call. = FALSE)
  }

  default_params <- list(
    objective = objective,
    metric = if (objective == "binary") "auc" else if (objective == "poisson") "poisson" else "l2",
    learning_rate = 0.05,
    num_leaves = 31,
    min_data_in_leaf = 50,
    feature_fraction = 0.8,
    bagging_fraction = 0.8,
    bagging_freq = 1,
    verbosity = -1,
    num_threads = 1
  )

  merged_params <- utils::modifyList(default_params, params)
  nrounds <- if (!is.null(merged_params$nrounds)) as.integer(merged_params$nrounds) else 150L
  merged_params$nrounds <- NULL

  dtrain <- lightgbm::lgb.Dataset(data = x_train, label = y_train)
  fit <- lightgbm::lgb.train(
    params = merged_params,
    data = dtrain,
    nrounds = nrounds,
    verbose = -1
  )

  pred <- as.numeric(stats::predict(fit, x_test))
  if (objective == "binary") {
    pred <- pmin(pmax(pred, 1e-4), 1 - 1e-4)
  } else {
    pred <- pmax(pred, 1e-6)
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
  if (learner == "lightgbm") {
    return(.fit_predict_lightgbm(x_train, y_train, x_test, objective = "regression", params = params))
  }

  if (learner == "glm") {
    return(.fit_predict_glm_poisson(x_train, y_train, x_test))
  }

  stop("Unsupported outcome learner: ", learner, call. = FALSE)
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

  fit <- stats::lm(
    y_resid_pb ~ d_resid * CellType,
    data = pb_df,
    weights = n_cells
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

#' Fit DEDML on Cell-Level Data
#'
#' Runs donor-blocked cross-fitted DEDML with configurable confounders,
#' treatment learner, and outcome learner.
#'
#' @param counts Gene-by-cell count matrix (rows are genes, columns are cells).
#' @param meta Cell metadata data.frame.
#' @param donor_id_col Donor identifier column in `meta`.
#' @param sample_id_col Sample identifier column in `meta`.
#' @param treatment_col Binary treatment indicator column in `meta` (0/1).
#' @param cell_type_col Cell type column in `meta`.
#' @param confounder_spec Optional list from [dedml_make_confounder_spec()]. If provided,
#' this overrides `treatment_confounders`, `treatment_cell_summaries`, and `outcome_confounders`.
#' @param treatment_confounders Donor-level confounders for treatment model.
#' @param treatment_cell_summaries Cell-level covariates to summarize by donor for treatment model.
#' @param outcome_confounders Cell-level confounders for outcome model.
#' @param focus_celltypes Optional vector of cell types to keep.
#' @param gene_subset Optional vector of gene names to keep.
#' @param n_folds Number of donor folds.
#' @param n_cores Number of cores for gene-level parallelization.
#' @param treatment_learner Treatment learner (`"glm"` or `"lightgbm"`).
#' @param outcome_learner Outcome learner (`"lightgbm"` or `"glm"`).
#' @param treatment_params Optional list of learner hyperparameters.
#' @param outcome_params Optional list of learner hyperparameters.
#' @param min_cells_per_sample_ct Minimum cells in each sample-celltype pseudobulk.
#' @param min_samples_per_celltype Minimum samples required per cell type.
#' @param seed Random seed.
#' @param verbose Print progress.
#'
#' @return A list with `results`, `donor_meta`, `settings`, and `errors`.
#' @export
#'
dedml_fit <- function(
    counts,
    meta,
    donor_id_col = "HATIMID",
    sample_id_col = "sample_id",
    treatment_col = "D",
    cell_type_col = "CellType",
    confounder_spec = NULL,
    treatment_confounders = c("Age", "Sex"),
    treatment_cell_summaries = c("log_nCount_RNA", "nFeature_RNA", "percent.mt"),
    outcome_confounders = c("Age", "Sex", "log_nCount_RNA", "nFeature_RNA", "percent.mt"),
    focus_celltypes = NULL,
    gene_subset = NULL,
    n_folds = 3L,
    n_cores = 1L,
    treatment_learner = "glm",
    outcome_learner = "lightgbm",
    treatment_params = list(),
    outcome_params = list(),
    min_cells_per_sample_ct = 10L,
    min_samples_per_celltype = 3L,
    seed = 123L,
    verbose = TRUE) {

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

  required_cols <- c(donor_id_col, sample_id_col, treatment_col, cell_type_col)
  missing_required <- setdiff(required_cols, colnames(meta))
  if (length(missing_required) > 0L) {
    stop("Missing required columns in meta: ", paste(missing_required, collapse = ", "), call. = FALSE)
  }

  meta <- as.data.frame(meta, stringsAsFactors = FALSE)
  meta$.donor_id <- as.character(meta[[donor_id_col]])
  meta$.sample_id <- as.character(meta[[sample_id_col]])
  meta$.D <- as.integer(meta[[treatment_col]])
  meta$.CellType <- as.character(meta[[cell_type_col]])

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

  if (!is.null(focus_celltypes)) {
    keep_cells <- meta$.CellType %in% focus_celltypes
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
    donor_id_col,
    sample_id_col,
    treatment_col,
    cell_type_col,
    treatment_confounders,
    treatment_cell_summaries,
    outcome_confounders
  ))
  needed_cols <- setdiff(needed_cols, c("", NA))

  missing_needed <- setdiff(needed_cols, colnames(meta))
  if (length(missing_needed) > 0L) {
    stop("Missing confounder columns in meta: ", paste(missing_needed, collapse = ", "), call. = FALSE)
  }

  keep_complete <- stats::complete.cases(meta[, needed_cols, drop = FALSE])
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
  donor_meta$propensity_oof <- NA_real_

  for (k in seq_len(n_folds)) {
    train_idx <- which(donor_meta$fold != k)
    test_idx <- which(donor_meta$fold == k)

    donor_meta$propensity_oof[test_idx] <- .fit_predict_treatment(
      x_train = x_treat[train_idx, , drop = FALSE],
      y_train = donor_meta$.D[train_idx],
      x_test = x_treat[test_idx, , drop = FALSE],
      learner = treatment_learner,
      params = treatment_params
    )
  }

  donor_meta$D_resid <- donor_meta$.D - donor_meta$propensity_oof

  meta <- dplyr::left_join(
    meta,
    donor_meta[, c(".donor_id", "D_resid"), drop = FALSE],
    by = ".donor_id"
  )

  x_outcome <- .build_design_matrix(meta, outcome_confounders)
  gene_names <- rownames(counts)

  if (is.null(gene_names)) {
    gene_names <- paste0("gene_", seq_len(nrow(counts)))
    rownames(counts) <- gene_names
  }

  run_one_gene <- function(g_idx) {
    y_obs <- as.numeric(counts[g_idx, ])
    mu_hat <- rep(NA_real_, length(y_obs))

    for (k in seq_len(n_folds)) {
      train_idx <- which(meta$fold != k)
      test_idx <- which(meta$fold == k)

      mu_hat[test_idx] <- .fit_predict_outcome(
        x_train = x_outcome[train_idx, , drop = FALSE],
        y_train = y_obs[train_idx],
        x_test = x_outcome[test_idx, , drop = FALSE],
        learner = outcome_learner,
        params = outcome_params
      )
    }

    y_resid <- (y_obs - mu_hat) / sqrt(pmax(mu_hat, 1e-6))

    rr <- fit_residual_regression_pseudobulk(
      y_resid = y_resid,
      d_resid = meta$D_resid,
      cell_type = meta$.CellType,
      sample_id = meta$.sample_id,
      donor_id = meta$.donor_id,
      min_cells_per_sample_ct = min_cells_per_sample_ct,
      min_samples = min_samples_per_celltype,
      cluster_robust = TRUE
    )

    if (is.null(rr)) {
      return(NULL)
    }

    rr$gene <- gene_names[g_idx]
    rr
  }

  if (isTRUE(verbose)) {
    message("Running DEDML on ", length(gene_names), " genes with ", n_cores, " core(s)")
  }

  if (n_cores > 1L) {
    raw_results <- parallel::mclapply(
      X = seq_along(gene_names),
      FUN = function(i) {
        tryCatch(run_one_gene(i), error = function(e) list(.error = TRUE, gene = gene_names[i], message = conditionMessage(e)))
      },
      mc.cores = n_cores
    )
  } else {
    raw_results <- lapply(
      X = seq_along(gene_names),
      FUN = function(i) {
        tryCatch(run_one_gene(i), error = function(e) list(.error = TRUE, gene = gene_names[i], message = conditionMessage(e)))
      }
    )
  }

  error_list <- Filter(function(x) is.list(x) && !is.data.frame(x) && isTRUE(x$.error), raw_results)
  result_list <- Filter(function(x) is.data.frame(x) && nrow(x) > 0L, raw_results)

  if (length(result_list) == 0L) {
    stop("No valid DEDML gene-level results were produced.", call. = FALSE)
  }

  result_df <- data.table::rbindlist(result_list, fill = TRUE)
  result_df <- dplyr::group_by(result_df, CellType)
  result_df <- dplyr::mutate(
    result_df,
    p_adj = stats::p.adjust(p_value, method = "BH"),
    rank_score = sign(estimate) * abs(statistic)
  )
  result_df <- dplyr::ungroup(result_df)

  settings <- list(
    donor_id_col = donor_id_col,
    sample_id_col = sample_id_col,
    treatment_col = treatment_col,
    cell_type_col = cell_type_col,
    treatment_confounders = treatment_confounders,
    treatment_cell_summaries = treatment_cell_summaries,
    outcome_confounders = outcome_confounders,
    n_folds = n_folds,
    treatment_learner = treatment_learner,
    outcome_learner = outcome_learner,
    min_cells_per_sample_ct = min_cells_per_sample_ct,
    min_samples_per_celltype = min_samples_per_celltype,
    seed = seed
  )

  list(
    results = as.data.frame(result_df),
    donor_meta = donor_meta,
    settings = settings,
    errors = if (length(error_list) > 0L) data.table::rbindlist(error_list, fill = TRUE) else NULL
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
#' @param donor_id_col Donor identifier column in `object@meta.data`.
#' @param sample_id_col Sample identifier column. If `NULL`, uses donor id.
#' @param treatment_col Binary treatment indicator column (0/1).
#' @param cell_type_col Cell type column.
#' @param ... Additional arguments passed to [dedml_fit()].
#'
#' @return A list with `results`, `donor_meta`, `settings`, and `errors`.
#' @export
#'
dedml_fit_seurat <- function(
    object,
    assay = "RNA",
    layer = "counts",
    donor_id_col = "HATIMID",
    sample_id_col = NULL,
    treatment_col = "D",
    cell_type_col = "celltype_manual",
    ...) {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Package 'Seurat' is required for dedml_fit_seurat().", call. = FALSE)
  }

  meta <- object@meta.data

  if (!(donor_id_col %in% colnames(meta))) {
    stop("donor_id_col not found in Seurat metadata: ", donor_id_col, call. = FALSE)
  }

  if (is.null(sample_id_col)) {
    sample_id_col <- donor_id_col
  }

  if (!(sample_id_col %in% colnames(meta))) {
    stop("sample_id_col not found in Seurat metadata: ", sample_id_col, call. = FALSE)
  }

  required <- c(treatment_col, cell_type_col)
  missing_required <- setdiff(required, colnames(meta))
  if (length(missing_required) > 0L) {
    stop("Missing required Seurat metadata columns: ", paste(missing_required, collapse = ", "), call. = FALSE)
  }

  counts <- Seurat::GetAssayData(object, assay = assay, layer = layer)

  dedml_fit(
    counts = counts,
    meta = meta,
    donor_id_col = donor_id_col,
    sample_id_col = sample_id_col,
    treatment_col = treatment_col,
    cell_type_col = cell_type_col,
    ...
  )
}
