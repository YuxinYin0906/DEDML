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

simulate_parallel_fixture <- function(
    n_donors = 24L,
    n_genes = 18L,
    n_celltypes = 3L,
    cells_per_donor_ct = 10L,
    seed = 20260429L) {
  set.seed(seed)
  donor_id <- sprintf("donor_%02d", seq_len(n_donors))
  donor <- data.frame(
    donor_id = donor_id,
    sample_id = donor_id,
    D = rep(c(0L, 1L), length.out = n_donors),
    Age = round(rnorm(n_donors, 55, 12)),
    Sex = factor(rep(c("F", "M"), length.out = n_donors)),
    stringsAsFactors = FALSE
  )
  donor$Age_z <- as.numeric(scale(donor$Age))
  celltypes <- paste0("CT", seq_len(n_celltypes))
  genes <- paste0("gene_", seq_len(n_genes))

  meta <- rbindlist(lapply(seq_len(n_donors), function(i) {
    d <- donor[i, ]
    rbindlist(lapply(celltypes, function(ct) {
      data.frame(
        cell_id = sprintf("%s_%s_cell_%02d", d$donor_id, ct, seq_len(cells_per_donor_ct)),
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
  donor_idx <- match(meta$donor_id, donor$donor_id)
  base_gene <- rnorm(n_genes, -1.2, 0.4)
  ct_eff <- matrix(rnorm(n_genes * n_celltypes, 0, 0.25), n_genes, n_celltypes)
  donor_re <- matrix(rnorm(n_genes * n_donors, 0, 0.12), n_genes, n_donors)
  de <- matrix(0, n_genes, n_celltypes)
  de[seq_len(3), 1] <- 0.5
  counts <- matrix(0L, nrow = n_genes, ncol = nrow(meta), dimnames = list(genes, meta$cell_id))
  for (g in seq_len(n_genes)) {
    eta <- base_gene[g] +
      ct_eff[g, ct_idx] +
      donor_re[g, donor_idx] +
      0.08 * as.numeric(scale(meta$Age)) +
      0.10 * (meta$Sex == "M") +
      0.20 * (meta$log_nCount_RNA - mean(meta$log_nCount_RNA)) +
      de[g, ct_idx] * meta$D
    mu <- pmax(exp(eta), 1e-4)
    counts[g, ] <- rnbinom(ncol(counts), mu = mu, size = 8)
  }

  list(counts = counts, meta = as.data.frame(meta), celltypes = celltypes, genes = genes)
}

run_case <- function(dat, case) {
  conf <- dedml_make_confounder_spec(
    donor_confounders = c("Age", "Sex"),
    cell_confounders = c("log_nCount_RNA", "nFeature_RNA", "percent.mt")
  )
  fit <- dedml_fit(
    counts = dat$counts,
    meta = dat$meta,
    donor_id = "donor_id",
    sample_id = "sample_id",
    treatment = "D",
    cell_type = "CellType",
    confounder_spec = conf,
    n_folds = 3L,
    n_cores = case$n_cores,
    donor_model = case$donor_model,
    outcome_model = case$outcome_model,
    treatment_params = case$params,
    outcome_params = case$params,
    min_cells_per_sample_ct = 5L,
    min_samples_per_celltype = 8L,
    parallel_backend = case$parallel_backend,
    seed = 11L,
    verbose = FALSE
  )
  res <- fit$results
  if (!is.data.frame(res) || nrow(res) == 0L) {
    stop("No result rows", call. = FALSE)
  }
  if (!all(is.finite(res$p_value))) {
    stop("Non-finite p-values", call. = FALSE)
  }
  data.frame(
    case = case$name,
    donor_model = case$donor_model,
    outcome_model = case$outcome_model,
    requested_backend = case$parallel_backend,
    used_backend = fit$settings$parallel_backend,
    n_cores = case$n_cores,
    n_results = nrow(res),
    min_p = min(res$p_value),
    max_p = max(res$p_value),
    stringsAsFactors = FALSE
  )
}

main <- function() {
  source(find_dedml_source())
  dat <- simulate_parallel_fixture()
  lgb_params <- list(
    nrounds = 3L,
    learning_rate = 0.1,
    num_leaves = 4L,
    min_data_in_leaf = 5L,
    verbosity = -1L,
    num_threads = 1L
  )
  cases <- list(
    list(name = "serial_glm_glm", donor_model = "glm", outcome_model = "glm", n_cores = 1L, parallel_backend = "auto", params = list()),
    list(name = "fork_glm_glm", donor_model = "glm", outcome_model = "glm", n_cores = 2L, parallel_backend = "auto", params = list()),
    list(name = "fork_glm_lgb", donor_model = "glm", outcome_model = "lightgbm", n_cores = 2L, parallel_backend = "auto", params = lgb_params),
    list(name = "psock_lgb_lgb_auto", donor_model = "lightgbm", outcome_model = "lightgbm", n_cores = 2L, parallel_backend = "auto", params = lgb_params),
    list(name = "psock_lgb_nb_auto", donor_model = "lightgbm", outcome_model = "nb_glm", n_cores = 2L, parallel_backend = "auto", params = lgb_params)
  )
  if (identical(Sys.getenv("DEDML_PARALLEL_STRESS"), "1")) {
    cases <- c(cases, list(list(
      name = "psock_lgb_lgb_20core_stress",
      donor_model = "lightgbm",
      outcome_model = "lightgbm",
      n_cores = min(20L, parallel::detectCores(logical = TRUE)),
      parallel_backend = "auto",
      params = lgb_params
    )))
  }
  out <- rbindlist(lapply(cases, function(case) {
    message("parallel check: ", case$name)
    run_case(dat, case)
  }))
  print(out)
}

main()
