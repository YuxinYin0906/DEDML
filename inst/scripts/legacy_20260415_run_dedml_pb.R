#!/usr/bin/env Rscript

Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1",
  BLIS_NUM_THREADS = "1"
)

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(data.table)
  library(lightgbm)
  library(parallel)
  library(sandwich)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript 20260423_run_dedml_pb_pipeline.R <ncores> <obj_rds> <outdir>")
}

ncores <- as.integer(args[1])
obj_rds <- args[2]
outdir <- args[3]

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

focus_celltypes <- c(
  "B cell Ig+",
  "B cell Ig-",
  "CD4-non-EM like",
  "CD4-EM like",
  "NK",
  "Monocyte-non-classical",
  "Monocyte-classical",
  "CD8-non-EM like",
  "CD8-EM like",
  "Monocyte-intermediate",
  "DC"
)

K_folds <- 3

coerce_age_sex_d <- function(df, sex_levels = NULL) {
  if ("Age" %in% colnames(df)) {
    df$Age <- suppressWarnings(as.numeric(as.character(df$Age)))
  }

  if ("Sex" %in% colnames(df)) {
    if (is.null(sex_levels)) {
      df$Sex <- as.factor(trimws(as.character(df$Sex)))
    } else {
      df$Sex <- factor(trimws(as.character(df$Sex)), levels = sex_levels)
    }
  }

  if ("D" %in% colnames(df)) {
    df$D <- as.integer(as.character(df$D))
  }

  df
}

make_stratified_donor_folds <- function(donor_meta, K = 3, seed = 1) {
  set.seed(seed)

  donors0 <- sample(donor_meta$HATIMID[donor_meta$D == 0])
  donors1 <- sample(donor_meta$HATIMID[donor_meta$D == 1])

  bind_rows(
    data.frame(HATIMID = donors0, fold = rep(seq_len(K), length.out = length(donors0))),
    data.frame(HATIMID = donors1, fold = rep(seq_len(K), length.out = length(donors1)))
  )
}

fit_treatment_glm <- function(train_df, test_df, covars) {
  train_df <- coerce_age_sex_d(train_df)
  test_df <- coerce_age_sex_d(test_df, sex_levels = levels(train_df$Sex))

  form <- as.formula(paste("D ~", paste(covars, collapse = " + ")))
  fit <- glm(form, data = train_df, family = binomial())
  pred <- predict(fit, newdata = test_df, type = "response")
  pmin(pmax(pred, 1e-4), 1 - 1e-4)
}

build_model_matrix <- function(df, covars) {
  df <- coerce_age_sex_d(df)
  mm <- model.matrix(as.formula(paste("~", paste(covars, collapse = " + "))), data = df)
  mm[, colnames(mm) != "(Intercept)", drop = FALSE]
}

fit_predict_lgb_gene <- function(
    x_train,
    y_train,
    x_test,
    seed = 1,
    nrounds = 150,
    learning_rate = 0.05,
    num_leaves = 31,
    min_data_in_leaf = 50,
    feature_fraction = 0.8,
    bagging_fraction = 0.8,
    bagging_freq = 1) {
  dtrain <- lgb.Dataset(data = x_train, label = y_train)

  params <- list(
    objective = "regression",
    metric = "l2",
    learning_rate = learning_rate,
    num_leaves = num_leaves,
    min_data_in_leaf = min_data_in_leaf,
    feature_fraction = feature_fraction,
    bagging_fraction = bagging_fraction,
    bagging_freq = bagging_freq,
    verbosity = -1,
    seed = seed,
    num_threads = 1
  )

  fit <- lgb.train(params = params, data = dtrain, nrounds = nrounds)
  pmax(predict(fit, x_test), 1e-6)
}

fit_residual_regression_interaction_pseudobulk <- function(
    y_resid,
    d_resid,
    celltype,
    sample_id,
    donor_id,
    min_cells_per_sample_ct = 10,
    min_samples = 3) {
  df <- data.frame(
    y_resid = y_resid,
    d_resid = d_resid,
    CellType = factor(celltype),
    sample_id = sample_id,
    donor_id = donor_id,
    stringsAsFactors = FALSE
  )

  sample_ct_counts <- df %>%
    count(sample_id, CellType, name = "n_cells")

  valid_sample_ct <- sample_ct_counts %>%
    filter(n_cells >= min_cells_per_sample_ct)

  df <- df %>%
    inner_join(valid_sample_ct %>% select(sample_id, CellType, n_cells), by = c("sample_id", "CellType"))

  if (nrow(df) == 0) return(NULL)

  pb_df <- df %>%
    group_by(sample_id, CellType) %>%
    summarise(
      y_resid_pb = mean(y_resid, na.rm = TRUE),
      d_resid = first(d_resid),
      donor_id = first(donor_id),
      n_cells = first(n_cells),
      .groups = "drop"
    )

  keep_ct <- pb_df %>%
    count(CellType, name = "n_samples") %>%
    filter(n_samples >= min_samples) %>%
    pull(CellType) %>%
    as.character()

  pb_df <- pb_df %>% filter(as.character(CellType) %in% keep_ct)

  if (nrow(pb_df) == 0 || length(unique(pb_df$d_resid)) < 2) return(NULL)

  pb_df$CellType <- droplevels(pb_df$CellType)

  fit <- lm(y_resid_pb ~ d_resid * CellType, data = pb_df, weights = n_cells)
  if (length(coef(fit)) == 0) return(NULL)

  repeated_donor_samples <- any(table(pb_df$donor_id) > 1)
  n_unique_donors <- length(unique(pb_df$donor_id))
  use_cluster_se <- repeated_donor_samples && n_unique_donors >= 2

  vc <- if (use_cluster_se) {
    sandwich::vcovCL(fit, cluster = pb_df$donor_id, type = "HC1")
  } else {
    vcov(fit)
  }

  ct_levels <- levels(pb_df$CellType)
  nd_1 <- data.frame(d_resid = rep(1, length(ct_levels)), CellType = factor(ct_levels, levels = ct_levels))
  nd_0 <- data.frame(d_resid = rep(0, length(ct_levels)), CellType = factor(ct_levels, levels = ct_levels))

  trm_noy <- delete.response(terms(fit))
  x1 <- model.matrix(trm_noy, data = nd_1)
  x0 <- model.matrix(trm_noy, data = nd_0)

  L <- x1 - x0
  est <- as.numeric(L %*% coef(fit))
  var_est <- pmax(diag(L %*% vc %*% t(L)), 0)
  std_error <- sqrt(var_est)
  statistic <- est / std_error
  p_value <- 2 * pnorm(-abs(statistic))

  out <- data.frame(
    CellType = ct_levels,
    estimate = est,
    std_error = std_error,
    statistic = statistic,
    p_value = p_value,
    se_method = ifelse(use_cluster_se, "donor_cluster_robust", "model_based"),
    repeated_donor_samples = repeated_donor_samples,
    row.names = NULL
  )

  ct_cell_counts <- pb_df %>%
    group_by(CellType) %>%
    summarise(
      n_pseudobulk = n(),
      n_cells = sum(n_cells),
      median_cells_per_sample = median(n_cells),
      n_samples = n_distinct(sample_id),
      n_donors = n_distinct(donor_id),
      .groups = "drop"
    )

  out %>%
    left_join(ct_cell_counts, by = "CellType") %>%
    select(
      CellType, estimate, std_error, statistic, p_value,
      n_donors, n_samples, n_pseudobulk, n_cells, median_cells_per_sample,
      se_method, repeated_donor_samples
    )
}

cat("Running DEDML pseudobulk pipeline\n")
cat("ncores:", ncores, "\n")
cat("object:", obj_rds, "\n")
cat("outdir:", outdir, "\n")

obj <- readRDS(obj_rds)

required_cols <- c("celltype_manual", "D", "condition", "nCount_RNA", "nFeature_RNA", "percent.mt")


if (!("HATIMID" %in% colnames(obj@meta.data))) {
  if ("Patient_ID" %in% colnames(obj@meta.data)) {
    obj$HATIMID <- obj$Patient_ID
  } else {
    stop("Expected donor ID column 'HATIMID' or 'Patient_ID' was not found.")
  }
}

sample_id_candidates <- c("sample_id", "SampleID", "sample", "orig.ident", "library_id")
sample_id_col <- sample_id_candidates[sample_id_candidates %in% colnames(obj@meta.data)][1]
if (!is.na(sample_id_col)) {
  obj$sample_id <- as.character(obj@meta.data[[sample_id_col]])
} else {
  obj$sample_id <- as.character(obj$HATIMID)
  cat("No sample ID column found; using HATIMID as sample_id.\n")
}

obj$CellType <- as.character(obj$celltype_manual)
focus_celltypes <- intersect(focus_celltypes, unique(obj$CellType))
if (length(focus_celltypes) == 0) {
  stop("None of the requested focus cell types were found in the object.")
}

obj <- subset(obj, cells = rownames(obj@meta.data)[obj$CellType %in% focus_celltypes])

cts_all <- GetAssayData(obj, assay = "RNA", layer = "counts")
obj$nCount_RNA <- Matrix::colSums(cts_all)
obj$nFeature_RNA <- Matrix::colSums(cts_all > 0)
obj$log_nCount_RNA <- log1p(obj$nCount_RNA)

meta <- obj@meta.data
meta$cell_id <- rownames(meta)
meta <- coerce_age_sex_d(meta)

donor_covars_base <- c("Age", "Sex")
cell_covars <- c("log_nCount_RNA", "nFeature_RNA", "percent.mt")
needed_cols <- c("HATIMID", "sample_id", "D", "CellType", donor_covars_base, cell_covars)

missing_needed <- setdiff(needed_cols, colnames(meta))
if (length(missing_needed) > 0) {
  stop("Missing needed columns: ", paste(missing_needed, collapse = ", "))
}

keep <- complete.cases(meta[, needed_cols, drop = FALSE])
obj <- subset(obj, cells = rownames(meta)[keep])
meta <- obj@meta.data
meta$cell_id <- rownames(meta)
meta <- coerce_age_sex_d(meta)

if (length(unique(meta$D)) < 2) {
  stop("Only one treatment group remains after filtering.")
}

cat("Cells after filtering:", nrow(meta), "\n")
cat("Donors after filtering:", length(unique(meta$HATIMID)), "\n")
cat("Samples after filtering:", length(unique(meta$sample_id)), "\n")

donor_cell_summ <- meta %>%
  group_by(HATIMID) %>%
  summarise(
    mean_log_nCount_RNA = mean(log_nCount_RNA, na.rm = TRUE),
    mean_nFeature_RNA = mean(nFeature_RNA, na.rm = TRUE),
    mean_percent_mt = mean(percent.mt, na.rm = TRUE),
    .groups = "drop"
  )

donor_meta <- meta %>%
  group_by(HATIMID) %>%
  summarise(
    D = first(D),
    Age = first(Age),
    Sex = first(Sex),
    .groups = "drop"
  ) %>%
  left_join(donor_cell_summ, by = "HATIMID")

donor_meta <- coerce_age_sex_d(donor_meta)

donor_covars_treatment <- c(
  donor_covars_base,
  "mean_log_nCount_RNA",
  "mean_nFeature_RNA",
  "mean_percent_mt"
)
outcome_covars <- c(donor_covars_base, cell_covars)

fold_df <- make_stratified_donor_folds(donor_meta, K = K_folds, seed = 123)
donor_meta <- donor_meta %>% left_join(fold_df, by = "HATIMID")
meta <- meta %>% left_join(fold_df, by = "HATIMID")
rownames(meta) <- meta$cell_id
meta <- coerce_age_sex_d(meta, sex_levels = levels(donor_meta$Sex))

donor_meta$propensity_oof <- NA_real_
for (k in seq_len(K_folds)) {
  train_donors <- donor_meta %>% filter(fold != k)
  test_donors <- donor_meta %>% filter(fold == k)

  donor_meta$propensity_oof[donor_meta$fold == k] <-
    fit_treatment_glm(train_donors, test_donors, donor_covars_treatment)
}

donor_meta <- donor_meta %>% mutate(D_resid = D - propensity_oof)

# map to each cell 
meta <- meta %>%
  left_join(donor_meta %>% select(HATIMID, D_resid), by = "HATIMID")
rownames(meta) <- meta$cell_id

count_mat <- GetAssayData(obj, assay = "RNA", layer = "counts")
count_mat <- count_mat[, meta$cell_id, drop = FALSE]
expr <- as.matrix(count_mat)
X_all <- build_model_matrix(meta, outcome_covars)

gene_names <- rownames(expr)
sample_id_vec <- meta$sample_id
donor_id_vec <- meta$HATIMID
d_resid_vec <- meta$D_resid
celltype_vec <- meta$CellType

run_one_gene <- function(g_idx, min_cells_per_sample_ct = 10, min_samples = 3) {
  g <- gene_names[g_idx]
  y_obs <- as.numeric(expr[g_idx, ])
  mu_hat <- rep(NA_real_, length(y_obs))

  for (k in seq_len(K_folds)) {
    train_idx <- which(meta$fold != k)
    test_idx <- which(meta$fold == k)

    mu_hat[test_idx] <- fit_predict_lgb_gene(
      x_train = X_all[train_idx, , drop = FALSE],
      y_train = y_obs[train_idx],
      x_test = X_all[test_idx, , drop = FALSE],
      seed = 123 + g_idx + k
    )
  }

  y_resid <- (y_obs - mu_hat) / sqrt(pmax(mu_hat, 1e-6))

  rr <- fit_residual_regression_interaction_pseudobulk(
    y_resid = y_resid,
    d_resid = d_resid_vec,
    celltype = celltype_vec,
    sample_id = sample_id_vec,
    donor_id = donor_id_vec,
    min_cells_per_sample_ct = min_cells_per_sample_ct,
    min_samples = min_samples
  )

  if (is.null(rr)) return(NULL)
  rr$gene <- g
  rr
}

cat("Running genes in parallel\n")

n_chunks <- min(ncores, length(gene_names))
gene_chunks <- split(seq_along(gene_names), cut(seq_along(gene_names), breaks = n_chunks, labels = FALSE))

run_gene_chunk <- function(idx_vec) {
  lapply(idx_vec, function(i) {
    tryCatch(
      run_one_gene(i),
      error = function(e) {
        list(.error = TRUE, gene = gene_names[i], message = conditionMessage(e))
      }
    )
  })
}

res_chunks <- mclapply(gene_chunks, run_gene_chunk, mc.cores = n_chunks)
res_list_raw <- unlist(res_chunks, recursive = FALSE)

error_list <- Filter(function(x) is.list(x) && !is.data.frame(x) && isTRUE(x$.error), res_list_raw)
res_list <- Filter(function(x) is.data.frame(x) && nrow(x) > 0, res_list_raw)

if (length(error_list) > 0) {
  err_df <- bind_rows(lapply(error_list, as.data.frame))
  fwrite(err_df, file.path(outdir, "dedml_gene_errors.csv"))
  cat("Gene-level errors:", nrow(err_df), "\n")
}

if (length(res_list) == 0) {
  if (length(error_list) > 0) {
    unique_err <- unique(vapply(error_list, function(x) x$message, character(1)))
    stop(
      "No valid cell-type-specific regression results were produced.\n",
      "First errors:\n",
      paste(head(unique_err, 8), collapse = "\n")
    )
  }
  stop("No valid cell-type-specific regression results were produced.")
}

res_df <- bind_rows(res_list)
if (!("CellType" %in% colnames(res_df)) && ("celltype" %in% colnames(res_df))) {
  res_df <- res_df %>% rename(CellType = celltype)
}
if (!("CellType" %in% colnames(res_df)) || nrow(res_df) == 0) {
  stop("Result table does not contain valid CellType output.")
}

res_df <- res_df %>%
  group_by(CellType) %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    rank_score = sign(estimate) * abs(statistic)
  ) %>%
  ungroup() %>%
  arrange(CellType, p_value)

fwrite(res_df, file.path(outdir, "dedml_covid_flu_selected_cts.csv"))
#saveRDS(donor_meta, file.path(outdir, "dedml_covid_flu_donor_meta.rds"))

cat("Done\n")
