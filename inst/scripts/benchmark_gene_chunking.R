#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

find_dedml_source <- function() {
  candidates <- c(
    file.path("R", "dedml.R"),
    file.path("DEDML", "R", "dedml.R"),
    file.path("..", "R", "dedml.R"),
    file.path("..", "DEDML", "R", "dedml.R")
  )
  hit <- candidates[file.exists(candidates)][1]
  if (is.na(hit)) stop("Could not find R/dedml.R", call. = FALSE)
  hit
}

simulate_dedml_fixture <- function(
    n_donors = 50L,
    n_genes = 1000L,
    n_celltypes = 3L,
    cells_per_donor_ct = 8L,
    seed = 20260429L) {
  set.seed(seed)
  donors <- sprintf("donor_%03d", seq_len(n_donors))
  celltypes <- paste0("CT", seq_len(n_celltypes))
  donor <- data.frame(
    donor_id = donors,
    sample_id = donors,
    D = rep(c(0L, 1L), length.out = n_donors),
    Age = round(rnorm(n_donors, 55, 11)),
    Sex = factor(rep(c("F", "M"), length.out = n_donors)),
    stringsAsFactors = FALSE
  )
  donor$Age_z <- as.numeric(scale(donor$Age))

  meta <- rbindlist(lapply(seq_len(n_donors), function(i) {
    d <- donor[i, ]
    rbindlist(lapply(celltypes, function(ct) {
      data.frame(
        cell_id = sprintf("%s_%s_%03d", d$donor_id, ct, seq_len(cells_per_donor_ct)),
        donor_id = d$donor_id,
        sample_id = d$sample_id,
        D = d$D,
        Age = d$Age,
        Sex = d$Sex,
        CellType = ct,
        log_nCount_RNA = rnorm(cells_per_donor_ct, 8 + 0.15 * d$Age_z, 0.25),
        nFeature_RNA = rnorm(cells_per_donor_ct, 900 + 30 * d$Age_z, 60),
        percent.mt = pmax(0.1, rnorm(cells_per_donor_ct, 5 + 0.2 * d$Age_z, 1)),
        stringsAsFactors = FALSE
      )
    }))
  }))
  rownames(meta) <- meta$cell_id

  ct_idx <- match(meta$CellType, celltypes)
  donor_idx <- match(meta$donor_id, donors)
  genes <- paste0("gene_", seq_len(n_genes))
  base_gene <- rnorm(n_genes, -1.2, 0.45)
  ct_eff <- matrix(rnorm(n_genes * n_celltypes, 0, 0.25), n_genes, n_celltypes)
  donor_re <- matrix(rnorm(n_genes * n_donors, 0, 0.12), n_genes, n_donors)
  de <- matrix(0, n_genes, n_celltypes)
  de[seq_len(min(50L, n_genes)), 1] <- 0.45

  counts <- matrix(0L, nrow = n_genes, ncol = nrow(meta), dimnames = list(genes, meta$cell_id))
  for (g in seq_len(n_genes)) {
    eta <- base_gene[g] +
      ct_eff[g, ct_idx] +
      donor_re[g, donor_idx] +
      0.07 * as.numeric(scale(meta$Age)) +
      0.10 * (meta$Sex == "M") +
      0.15 * (meta$log_nCount_RNA - mean(meta$log_nCount_RNA)) +
      de[g, ct_idx] * meta$D
    counts[g, ] <- rnbinom(ncol(counts), mu = pmax(exp(eta), 1e-4), size = 8)
  }

  list(counts = counts, meta = as.data.frame(meta), celltypes = celltypes)
}

time_fit <- function(label, dat, donor_model, outcome_model, n_cores, gene_chunk_size, params = list()) {
  conf <- dedml_make_confounder_spec(
    donor_confounders = c("Age", "Sex"),
    cell_confounders = c("log_nCount_RNA", "nFeature_RNA", "percent.mt")
  )
  gc()
  t0 <- proc.time()[["elapsed"]]
  fit <- dedml_fit(
    counts = dat$counts,
    meta = dat$meta,
    donor_id = "donor_id",
    sample_id = "sample_id",
    treatment = "D",
    cell_type = "CellType",
    confounder_spec = conf,
    n_folds = 3L,
    n_cores = n_cores,
    donor_model = donor_model,
    outcome_model = outcome_model,
    treatment_params = params,
    outcome_params = params,
    min_cells_per_sample_ct = 5L,
    min_samples_per_celltype = 20L,
    final_stage = "residual_lm",
    parallel_backend = "auto",
    gene_chunk_size = gene_chunk_size,
    seed = 11L,
    verbose = FALSE
  )
  elapsed <- proc.time()[["elapsed"]] - t0
  list(label = label, fit = fit, elapsed = elapsed)
}

compare_results <- function(a, b) {
  ma <- merge(a$fit$results, b$fit$results, by = c("gene", "CellType"), suffixes = c("_a", "_b"))
  data.frame(
    n_rows_a = nrow(a$fit$results),
    n_rows_b = nrow(b$fit$results),
    n_rows_merged = nrow(ma),
    max_abs_estimate_diff = max(abs(ma$estimate_a - ma$estimate_b)),
    max_abs_se_diff = max(abs(ma$std_error_a - ma$std_error_b)),
    max_abs_p_diff = max(abs(ma$p_value_a - ma$p_value_b)),
    stringsAsFactors = FALSE
  )
}

run_pair <- function(name, dat, donor_model, outcome_model, n_cores, params = list()) {
  old_style <- time_fit(
    label = paste0(name, "_chunk1"),
    dat = dat,
    donor_model = donor_model,
    outcome_model = outcome_model,
    n_cores = n_cores,
    gene_chunk_size = 1L,
    params = params
  )
  auto_chunk <- time_fit(
    label = paste0(name, "_auto_chunk"),
    dat = dat,
    donor_model = donor_model,
    outcome_model = outcome_model,
    n_cores = n_cores,
    gene_chunk_size = NULL,
    params = params
  )
  cmp <- compare_results(old_style, auto_chunk)
  data.frame(
    scenario = name,
    donor_model = donor_model,
    outcome_model = outcome_model,
    n_cores = n_cores,
    n_genes = nrow(dat$counts),
    n_cells = ncol(dat$counts),
    chunk1_seconds = old_style$elapsed,
    auto_chunk_seconds = auto_chunk$elapsed,
    speedup = old_style$elapsed / auto_chunk$elapsed,
    chunk1_used_chunk = old_style$fit$settings$gene_chunk_size,
    auto_used_chunk = auto_chunk$fit$settings$gene_chunk_size,
    chunk1_n_chunks = old_style$fit$settings$n_gene_chunks,
    auto_n_chunks = auto_chunk$fit$settings$n_gene_chunks,
    cmp,
    stringsAsFactors = FALSE
  )
}

main <- function() {
  source(find_dedml_source())
  n_cores <- min(20L, parallel::detectCores(logical = TRUE))

  glm_dat <- simulate_dedml_fixture(n_genes = 1000L, seed = 1L)
  glm_summary <- run_pair(
    name = "glm_glm_1000g",
    dat = glm_dat,
    donor_model = "glm",
    outcome_model = "glm",
    n_cores = n_cores
  )

  lgb_dat <- simulate_dedml_fixture(n_genes = 200L, seed = 2L)
  lgb_params <- list(
    nrounds = 3L,
    learning_rate = 0.1,
    num_leaves = 4L,
    min_data_in_leaf = 5L,
    feature_fraction = 1,
    bagging_fraction = 1,
    bagging_freq = 0,
    verbosity = -1L,
    num_threads = 1L,
    seed = 20260429L
  )
  lgb_summary <- run_pair(
    name = "lgb_lgb_200g",
    dat = lgb_dat,
    donor_model = "lightgbm",
    outcome_model = "lightgbm",
    n_cores = n_cores,
    params = lgb_params
  )

  summary <- rbind(glm_summary, lgb_summary)
  outdir <- file.path("DEDML", "inst", "benchmark_results")
  if (!dir.exists(outdir)) {
    outdir <- file.path("inst", "benchmark_results")
  }
  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
  }
  fwrite(summary, file.path(outdir, "gene_chunking_benchmark_summary.csv"))
  print(summary)
}

main()
