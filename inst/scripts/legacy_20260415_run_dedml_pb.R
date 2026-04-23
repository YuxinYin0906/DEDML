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
  library(lmtest)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Usage: Rscript run_dedml_covid_flu.R <ncores> <obj_rds> <outdir>")
}

ncores  <- as.integer(trimws(args[1]))
obj_rds <- trimws(args[2])
outdir  <- trimws(args[3])

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

cat("Mode: COVID vs Flu DEDML on selected immune cell types\n")
cat("Outcome uses raw count-like values directly\n")
cat("Outcome residual uses cross-fitted Pearson residuals\n")
cat("Treatment model: donor-level logistic GLM using donor + donor-summarized technical covariates\n")
cat("Outcome model: cell-level LightGBM using donor + cell-level covariates\n")
cat("Final stage: pseudobulk residual regression on sample x CellType (with intercept)\n")
cat("Cores:", ncores, "\n")
cat("Object:", obj_rds, "\n")
cat("Outdir:", outdir, "\n")

# CHANGED: updated to your requested cell types
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

# CHANGED: smaller donor count, use 3 folds instead of 5
K_folds <- 3

make_stratified_donor_folds <- function(donor_meta, K = 3, seed = 1) {
  set.seed(seed)
  
  donors0 <- donor_meta$HATIMID[donor_meta$D == 0]
  donors1 <- donor_meta$HATIMID[donor_meta$D == 1]
  
  donors0 <- sample(donors0)
  donors1 <- sample(donors1)
  
  fold_id0 <- rep(seq_len(K), length.out = length(donors0))
  fold_id1 <- rep(seq_len(K), length.out = length(donors1))
  
  bind_rows(
    data.frame(HATIMID = donors0, fold = fold_id0),
    data.frame(HATIMID = donors1, fold = fold_id1)
  )
}

fit_treatment_glm <- function(train_df, test_df, covars) {
  # CHANGED: force correct data types to avoid
  # "factor Age has new levels ..." during predict.glm
  if ("Age" %in% colnames(train_df)) {
    train_df$Age <- suppressWarnings(as.numeric(as.character(train_df$Age)))
  }
  if ("Age" %in% colnames(test_df)) {
    test_df$Age <- suppressWarnings(as.numeric(as.character(test_df$Age)))
  }
  
  if ("Sex" %in% colnames(train_df)) {
    train_df$Sex <- as.factor(trimws(as.character(train_df$Sex)))
  }
  if ("Sex" %in% colnames(test_df)) {
    test_df$Sex <- factor(
      trimws(as.character(test_df$Sex)),
      levels = levels(train_df$Sex)
    )
  }
  
  form <- as.formula(paste("D ~", paste(covars, collapse = " + ")))
  fit <- glm(form, data = train_df, family = binomial())
  pred <- predict(fit, newdata = test_df, type = "response")
  pmin(pmax(pred, 1e-4), 1 - 1e-4)
}

build_model_matrix <- function(df, covars) {
  # CHANGED: make sure Age is numeric and Sex is factor here too
  if ("Age" %in% colnames(df)) {
    df$Age <- suppressWarnings(as.numeric(as.character(df$Age)))
  }
  if ("Sex" %in% colnames(df)) {
    df$Sex <- as.factor(trimws(as.character(df$Sex)))
  }
  
  mm <- model.matrix(
    as.formula(paste("~", paste(covars, collapse = " + "))),
    data = df
  )
  mm[, colnames(mm) != "(Intercept)", drop = FALSE]
}

fit_predict_lgb_gene <- function(x_train, y_train, x_test,
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
  
  fit <- lgb.train(
    params = params,
    data = dtrain,
    nrounds = nrounds
  )
  
  pred <- predict(fit, x_test)
  pmax(pred, 1e-6)
}

fit_residual_regression_interaction_pseudobulk <- function(
    y_resid,
    d_resid,
    celltype,
    sample_id,
    donor_id,
    min_cells_per_sample_ct = 10,
    min_samples = 3
) {
  df <- data.frame(
    y_resid  = y_resid,
    d_resid  = d_resid,
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
    inner_join(
      valid_sample_ct %>% select(sample_id, CellType, n_cells),
      by = c("sample_id", "CellType")
    )
  
  if (nrow(df) == 0) return(NULL)
  
  # Pseudobulk unit is sample x celltype; donor_id is carried for clustered SE only.
  pb_df <- df %>%
    group_by(sample_id, CellType) %>%
    summarise(
      y_resid_pb = mean(y_resid, na.rm = TRUE),
      d_resid    = first(d_resid),
      donor_id   = first(donor_id),
      n_cells    = first(n_cells),
      .groups = "drop"
    )
  
  sample_counts_by_ct <- pb_df %>%
    count(CellType, name = "n_samples")
  
  keep_ct <- sample_counts_by_ct %>%
    filter(n_samples >= min_samples) %>%
    pull(CellType) %>%
    as.character()
  
  pb_df <- pb_df %>%
    filter(as.character(CellType) %in% keep_ct)
  
  if (nrow(pb_df) == 0) return(NULL)
  if (length(unique(pb_df$d_resid)) < 2) return(NULL)
  
  pb_df$CellType <- droplevels(pb_df$CellType)
  
  # Intercept-added model with treatment-by-celltype interactions.
  fit <- lm(
    y_resid_pb ~ d_resid * CellType,
    data = pb_df,
    weights = n_cells
  )
  
  if (length(coef(fit)) == 0) return(NULL)
  
  repeated_donor_samples <- any(table(pb_df$donor_id) > 1)
  n_unique_donors <- length(unique(pb_df$donor_id))
  use_donor_cluster_se <- repeated_donor_samples && n_unique_donors >= 2
  
  vc <- if (use_donor_cluster_se) {
    sandwich::vcovCL(fit, cluster = pb_df$donor_id, type = "HC1")
  } else {
    vcov(fit)
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
  
  # Build contrast rows from newdata (celltype levels), not from training frame.
  trm_noy <- delete.response(terms(fit))
  X1 <- model.matrix(trm_noy, data = nd_1)
  X0 <- model.matrix(trm_noy, data = nd_0)
  L <- X1 - X0
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
    se_method = ifelse(use_donor_cluster_se, "donor_cluster_robust", "model_based"),
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

obj <- readRDS(obj_rds)

# CHANGED: required columns adapted to your Seurat object
required_cols <- c(
  "celltype_manual", "D", "condition",
  "nCount_RNA", "nFeature_RNA", "percent.mt"
)

missing_required <- setdiff(required_cols, colnames(obj@meta.data))
if (length(missing_required) > 0) {
  stop("Missing required metadata columns: ", paste(missing_required, collapse = ", "))
}

# CHANGED: use HATIMID if present, otherwise create from Patient_ID
if (!("HATIMID" %in% colnames(obj@meta.data))) {
  if ("Patient_ID" %in% colnames(obj@meta.data)) {
    obj$HATIMID <- obj$Patient_ID
  } else {
    stop("Expected donor ID column 'HATIMID' or 'Patient_ID' was not found.")
  }
}

# CHANGED: derive a sample_id column for sample x CellType pseudobulk
sample_id_candidates <- c("sample_id", "SampleID", "sample", "orig.ident", "library_id")
sample_id_col <- sample_id_candidates[sample_id_candidates %in% colnames(obj@meta.data)][1]
if (!is.na(sample_id_col)) {
  obj$sample_id <- as.character(obj@meta.data[[sample_id_col]])
  cat("Using sample ID column:", sample_id_col, "\n")
} else {
  obj$sample_id <- as.character(obj$HATIMID)
  cat("Warning: no sample ID column found; using HATIMID as sample_id.\n")
}

# CHANGED: use your manual cell type column
obj$CellType <- as.character(obj$celltype_manual)

missing_ct <- setdiff(focus_celltypes, unique(obj$CellType))
if (length(missing_ct) > 0) {
  cat("Warning: these selected labels were not found and will be ignored:\n")
  print(missing_ct)
}
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

# CHANGED: force correct types early
if ("Age" %in% colnames(meta)) {
  meta$Age <- suppressWarnings(as.numeric(as.character(meta$Age)))
}
if ("Sex" %in% colnames(meta)) {
  meta$Sex <- as.factor(trimws(as.character(meta$Sex)))
}
if ("D" %in% colnames(meta)) {
  meta$D <- as.integer(as.character(meta$D))
}

cat("Selected cell types:\n")
print(sort(unique(meta$CellType)))
cat("Initial selected cells:", nrow(meta), "\n")
cat("Initial D table:\n")
print(table(meta$D, useNA = "ifany"))
cat("Initial condition x D table:\n")
print(table(meta$condition, meta$D, useNA = "ifany"))
cat("Initial CellType x D table:\n")
print(table(meta$CellType, meta$D, useNA = "ifany"))
cat("Unique donors:", length(unique(meta$HATIMID)), "\n")
cat("Unique samples:", length(unique(meta$sample_id)), "\n")

# CHANGED: COVID/flu covariates
donor_covars_base <- c("Age", "Sex")
cell_covars <- c("log_nCount_RNA", "nFeature_RNA", "percent.mt")
needed_cols <- c("HATIMID", "sample_id", "D", "CellType", donor_covars_base, cell_covars)

missing_needed <- setdiff(needed_cols, colnames(meta))
if (length(missing_needed) > 0) {
  stop("Missing needed columns after loading object: ", paste(missing_needed, collapse = ", "))
}

keep <- complete.cases(meta[, needed_cols, drop = FALSE])
obj <- subset(obj, cells = rownames(meta)[keep])
meta <- obj@meta.data
meta$cell_id <- rownames(meta)

# CHANGED: re-enforce types after filtering
meta$Age <- suppressWarnings(as.numeric(as.character(meta$Age)))
meta$Sex <- as.factor(trimws(as.character(meta$Sex)))
meta$D <- as.integer(as.character(meta$D))

cat("After complete-case filtering cells:", nrow(meta), "\n")
cat("After filtering D table:\n")
print(table(meta$D, useNA = "ifany"))
cat("After filtering CellType x D table:\n")
print(table(meta$CellType, meta$D, useNA = "ifany"))
cat("Unique donors after filtering:", length(unique(meta$HATIMID)), "\n")
cat("Unique samples after filtering:", length(unique(meta$sample_id)), "\n")
donor_sample_counts <- meta %>%
  distinct(HATIMID, sample_id) %>%
  count(HATIMID, name = "n_samples")
cat("Donors with >1 sample:", sum(donor_sample_counts$n_samples > 1), "\n")

if (length(unique(meta$D)) < 2) {
  stop("Only one treatment group remains after filtering.")
}

# CHANGED: summarize technical covariates that actually exist
donor_cell_summ <- meta %>%
  group_by(HATIMID) %>%
  summarise(
    mean_log_nCount_RNA = mean(log_nCount_RNA, na.rm = TRUE),
    mean_nFeature_RNA   = mean(nFeature_RNA, na.rm = TRUE),
    mean_percent_mt     = mean(percent.mt, na.rm = TRUE),
    .groups = "drop"
  )

# CHANGED: donor-level covariates are Age/Sex only; coerce types after summarize
donor_meta <- meta %>%
  group_by(HATIMID) %>%
  summarise(
    D   = first(D),
    Age = first(Age),
    Sex = first(Sex),
    .groups = "drop"
  ) %>%
  left_join(donor_cell_summ, by = "HATIMID")

donor_meta$Age <- suppressWarnings(as.numeric(as.character(donor_meta$Age)))
donor_meta$Sex <- as.factor(trimws(as.character(donor_meta$Sex)))
donor_meta$D <- as.integer(as.character(donor_meta$D))

cat("Donor-level D table:\n")
print(table(donor_meta$D, useNA = "ifany"))
cat("Number of donors in donor_meta:", nrow(donor_meta), "\n")
cat("Donor Age summary:\n")
print(summary(donor_meta$Age))
cat("Donor Sex table:\n")
print(table(donor_meta$Sex, useNA = "ifany"))

# CHANGED: treatment model uses donor summaries of technical covariates
donor_covars_treatment <- c(
  donor_covars_base,
  "mean_log_nCount_RNA",
  "mean_nFeature_RNA",
  "mean_percent_mt"
)

outcome_covars <- c(
  donor_covars_base,
  cell_covars
)

cat("Treatment model covariates:\n")
print(donor_covars_treatment)
cat("Outcome model covariates:\n")
print(outcome_covars)

fold_df <- make_stratified_donor_folds(donor_meta, K = K_folds, seed = 123)

donor_meta <- donor_meta %>% left_join(fold_df, by = "HATIMID")
meta <- meta %>% left_join(fold_df, by = "HATIMID")
rownames(meta) <- meta$cell_id

# CHANGED: re-enforce types again after joins
meta$Age <- suppressWarnings(as.numeric(as.character(meta$Age)))
meta$Sex <- factor(trimws(as.character(meta$Sex)), levels = levels(donor_meta$Sex))
meta$D <- as.integer(as.character(meta$D))

cat("Fold sizes:\n")
print(table(donor_meta$fold, donor_meta$D))

donor_meta$propensity_oof <- NA_real_
for (k in seq_len(K_folds)) {
  train_donors <- donor_meta %>% filter(fold != k)
  test_donors  <- donor_meta %>% filter(fold == k)
  
  donor_meta$propensity_oof[donor_meta$fold == k] <-
    fit_treatment_glm(train_donors, test_donors, donor_covars_treatment)
}

donor_meta <- donor_meta %>% mutate(D_resid = D - propensity_oof)

cat("Propensity summary:\n")
print(summary(donor_meta$propensity_oof))
cat("D residual summary:\n")
print(summary(donor_meta$D_resid))

meta <- meta %>%
  left_join(donor_meta %>% select(HATIMID, D_resid), by = "HATIMID")
rownames(meta) <- meta$cell_id

count_mat <- GetAssayData(obj, assay = "RNA", layer = "counts")
count_mat <- count_mat[, meta$cell_id, drop = FALSE]
expr <- as.matrix(count_mat)

cat("Genes retained:", nrow(expr), "\n")
cat("Using raw count-like values directly (no log transformation).\n")

X_all <- build_model_matrix(meta, outcome_covars)

gene_names   <- rownames(expr)
sample_id_vec <- meta$sample_id
donor_id_vec <- meta$HATIMID
d_resid_vec  <- meta$D_resid
celltype_vec <- meta$CellType

run_one_gene <- function(g_idx,
                         min_cells_per_sample_ct = 10,
                         min_samples = 3) {
  g <- gene_names[g_idx]
  y_obs <- as.numeric(expr[g_idx, ])
  mu_hat <- rep(NA_real_, length(y_obs))
  
  for (k in seq_len(K_folds)) {
    train_idx <- which(meta$fold != k)
    test_idx  <- which(meta$fold == k)
    
    pred_k <- fit_predict_lgb_gene(
      x_train = X_all[train_idx, , drop = FALSE],
      y_train = y_obs[train_idx],
      x_test  = X_all[test_idx, , drop = FALSE],
      seed    = 123 + g_idx + k
    )
    
    mu_hat[test_idx] <- pred_k
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

cat("Running genes in parallel by chunks...\n")

n_chunks <- min(ncores, length(gene_names))
gene_chunks <- split(
  seq_along(gene_names),
  cut(seq_along(gene_names), breaks = n_chunks, labels = FALSE)
)

run_gene_chunk <- function(idx_vec) {
  lapply(idx_vec, function(i) {
    tryCatch(
      run_one_gene(i),
      error = function(e) {
        list(
          .error = TRUE,
          gene = gene_names[i],
          message = conditionMessage(e)
        )
      }
    )
  })
}

res_chunks <- mclapply(
  gene_chunks,
  run_gene_chunk,
  mc.cores = n_chunks
)

res_list_raw <- unlist(res_chunks, recursive = FALSE)

error_list <- Filter(
  function(x) is.list(x) && !is.data.frame(x) && isTRUE(x$.error),
  res_list_raw
)

res_list <- Filter(
  function(x) is.data.frame(x) && nrow(x) > 0,
  res_list_raw
)

if (length(error_list) > 0) {
  err_df <- bind_rows(lapply(error_list, as.data.frame))
  fwrite(err_df, file.path(outdir, "dedml_gene_errors.csv"))
  cat("Gene-level errors detected:", nrow(err_df), "\n")
  cat("Wrote error details to:", file.path(outdir, "dedml_gene_errors.csv"), "\n")
}

if (length(res_list) == 0) {
  if (length(error_list) > 0) {
    unique_err <- unique(vapply(error_list, function(x) x$message, character(1)))
    stop(
      "No valid cell-type-specific regression results were produced.\n",
      "First errors:\n",
      paste(head(unique_err, 8), collapse = "\n")
    )
  } else {
    stop("No valid cell-type-specific regression results were produced. Check sample counts and cell counts per sample-celltype.")
  }
}

res_df <- bind_rows(res_list)

if (!("CellType" %in% colnames(res_df)) && ("celltype" %in% colnames(res_df))) {
  res_df <- res_df %>% rename(CellType = celltype)
}

if (!("CellType" %in% colnames(res_df))) {
  stop("Internal error: result table does not contain CellType after gene-level aggregation.")
}

if (nrow(res_df) == 0) {
  stop("No valid cell-type-specific regression results were produced. Check sample counts and cell counts per sample-celltype.")
}

res_df <- res_df %>%
  group_by(CellType) %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    rank_score = sign(estimate) * abs(statistic)
  ) %>%
  ungroup() %>%
  arrange(CellType, p_value)

fwrite(
  res_df,
  file.path(outdir, "dedml_covid_flu_selected_cts.csv")
)

saveRDS(
  donor_meta,
  file.path(outdir, "dedml_covid_flu_donor_meta.rds")
)

cat("Done: donor-blocked DEDML for COVID vs Flu with cell-level nuisance models and pseudobulk residual regression.\n")
