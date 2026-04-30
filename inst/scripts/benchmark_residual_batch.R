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

simulate_residual_fixture <- function(
    n_donors = 50L,
    n_genes = 2000L,
    n_celltypes = 6L,
    cells_per_donor_ct = 20L,
    seed = 20260429L) {
  set.seed(seed)
  donors <- sprintf("donor_%03d", seq_len(n_donors))
  celltypes <- paste0("CT", seq_len(n_celltypes))
  donor_resid <- rnorm(n_donors, sd = 0.45)
  meta <- rbindlist(lapply(seq_len(n_donors), function(i) {
    rbindlist(lapply(celltypes, function(ct) {
      data.frame(
        sample_id = donors[i],
        donor_id = donors[i],
        CellType = ct,
        d_resid = donor_resid[i],
        cell_index = seq_len(cells_per_donor_ct),
        stringsAsFactors = FALSE
      )
    }))
  }))

  gene_names <- paste0("gene_", seq_len(n_genes))
  ct_idx <- match(meta$CellType, celltypes)
  donor_idx <- match(meta$donor_id, donors)
  gene_scale <- runif(n_genes, 0.6, 1.4)
  ct_shift <- matrix(rnorm(n_genes * n_celltypes, sd = 0.15), n_genes, n_celltypes)
  donor_shift <- matrix(rnorm(n_genes * n_donors, sd = 0.05), n_genes, n_donors)
  signal <- matrix(0, n_genes, n_celltypes)
  signal[seq_len(min(100L, n_genes)), 1] <- 0.25
  signal[seq_len(min(80L, n_genes)), 3] <- -0.20

  y_resid <- matrix(0, nrow = n_genes, ncol = nrow(meta), dimnames = list(gene_names, NULL))
  for (g in seq_len(n_genes)) {
    mean_g <- ct_shift[g, ct_idx] +
      donor_shift[g, donor_idx] +
      signal[g, ct_idx] * meta$d_resid
    y_resid[g, ] <- rnorm(nrow(meta), mean = mean_g, sd = gene_scale[g])
  }

  list(y_resid = y_resid, meta = as.data.frame(meta), gene_names = gene_names)
}

time_expr <- function(expr) {
  gc()
  t0 <- proc.time()[["elapsed"]]
  value <- force(expr)
  elapsed <- proc.time()[["elapsed"]] - t0
  list(value = value, elapsed = elapsed)
}

main <- function() {
  source(find_dedml_source())
  dat <- simulate_residual_fixture()
  meta <- dat$meta
  y_resid <- dat$y_resid

  scalar <- time_expr({
    rbindlist(lapply(seq_len(nrow(y_resid)), function(g) {
      rr <- fit_residual_regression_pseudobulk(
        y_resid = y_resid[g, ],
        d_resid = meta$d_resid,
        cell_type = meta$CellType,
        sample_id = meta$sample_id,
        donor_id = meta$donor_id,
        min_cells_per_sample_ct = 10L,
        min_samples = 20L,
        cluster_robust = TRUE
      )
      if (is.null(rr)) return(NULL)
      rr$gene <- rownames(y_resid)[g]
      rr
    }), fill = TRUE)
  })

  batched <- time_expr({
    fit_residual_regression_pseudobulk_matrix(
      y_resid_mat = y_resid,
      d_resid = meta$d_resid,
      cell_type = meta$CellType,
      sample_id = meta$sample_id,
      donor_id = meta$donor_id,
      gene_names = rownames(y_resid),
      min_cells_per_sample_ct = 10L,
      min_samples = 20L,
      cluster_robust = TRUE
    )
  })

  cmp <- merge(
    as.data.frame(scalar$value),
    as.data.frame(batched$value),
    by = c("gene", "CellType"),
    suffixes = c("_scalar", "_batched")
  )

  summary <- data.frame(
    n_donors = length(unique(meta$donor_id)),
    n_genes = nrow(y_resid),
    n_celltypes = length(unique(meta$CellType)),
    n_cells = ncol(y_resid),
    n_result_rows_scalar = nrow(scalar$value),
    n_result_rows_batched = nrow(batched$value),
    scalar_seconds = scalar$elapsed,
    batched_seconds = batched$elapsed,
    speedup = scalar$elapsed / batched$elapsed,
    max_abs_estimate_diff = max(abs(cmp$estimate_scalar - cmp$estimate_batched)),
    max_abs_se_diff = max(abs(cmp$std_error_scalar - cmp$std_error_batched)),
    max_abs_p_diff = max(abs(cmp$p_value_scalar - cmp$p_value_batched)),
    stringsAsFactors = FALSE
  )

  outdir <- file.path("DEDML", "inst", "benchmark_results")
  if (!dir.exists(outdir)) {
    outdir <- file.path("inst", "benchmark_results")
  }
  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
  }
  fwrite(summary, file.path(outdir, "residual_batch_benchmark_summary.csv"))
  print(summary)
}

main()
