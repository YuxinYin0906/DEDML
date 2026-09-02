test_that("structured pseudobulk uses sums, offset, and alpha columns", {
  counts <- matrix(
    c(
      1, 2, 3, 4,
      5, 7, 11, 13,
      0, 1, 0, 1
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(c("g1", "g2", "MT-g3"), paste0("c", 1:4))
  )
  meta <- data.frame(
    donor = c("d1", "d1", "d2", "d2"),
    sample = c("s1", "s1", "s2", "s2"),
    D = c(0, 0, 1, 1),
    ct = c("T", "T", "B", "B"),
    Age = c(30, 30, 50, 50),
    Sex = c("F", "F", "M", "M"),
    stringsAsFactors = FALSE
  )

  pb <- DEDML:::.dedml_make_sum_pseudobulk(
    counts = counts,
    meta = meta,
    donor_id = "donor",
    sample_id = "sample",
    treatment = "D",
    cell_type = "ct",
    cell_types = NULL,
    treatment_confounders = c("Age", "Sex"),
    outcome_confounders = c("Age", "Sex", "log_n_cells", "nFeature_RNA", "percent.mt")
  )

  expect_equal(as.numeric(pb$counts["g1", "s1__T"]), 3)
  expect_equal(as.numeric(pb$counts["g2", "s2__B"]), 24)
  expect_equal(pb$meta$.dedml_pb_offset, log(pmax(DEDML:::.dedml_matrix_col_sums(pb$counts), 1)))
  expect_true(all(paste0(".dedml_alpha_group_", c("B", "T")) %in% names(pb$meta)))
})

test_that("structured tuning selection uses the current score", {
  grid <- DEDML:::.dedml_structured_pb_tuning_grid("nb")
  metrics <- data.table::copy(grid)
  metrics[, `:=`(
    n_tests = c(10L, 10L, 10L),
    median_abs_estimate = c(1, 1, 1),
    median_abs_estimate_logfc = c(1, 1, 1),
    elapsed_sec = c(30, 10, 20),
    n_gene_errors = c(0L, 1L, 0L)
  )]

  selected <- DEDML:::.dedml_select_structured_pb_grid(metrics, grid)
  expect_equal(selected$variant, "fold5_n300_lr002_leaves31_min20")
  expect_equal(selected$tuning_score, 20 / 3600)
})

test_that("structured pseudobulk pipeline runs a small poisson fit", {
  testthat::skip_if_not_installed("lightgbm")
  testthat::skip_if_not_installed("Matrix")

  set.seed(1)
  donors <- paste0("d", 1:10)
  meta <- do.call(rbind, lapply(seq_along(donors), function(i) {
    data.frame(
      cell = paste0(donors[i], "_", rep(c("T", "B"), each = 2), "_", rep(1:2, times = 2)),
      donor = donors[i],
      sample = donors[i],
      D = as.integer(i > 5),
      ct = rep(c("T", "B"), each = 2),
      Age = 35 + i,
      Sex = ifelse(i %% 2, "F", "M"),
      percent.mt = 1,
      stringsAsFactors = FALSE
    )
  }))
  rownames(meta) <- meta$cell
  counts <- matrix(rpois(10 * nrow(meta), lambda = 5), nrow = 10)
  rownames(counts) <- paste0("g", 1:10)
  colnames(counts) <- meta$cell

  out <- suppressWarnings(dedml_run_structured_pseudobulk_pipeline(
    counts = counts,
    meta = meta,
    donor_id = "donor",
    sample_id = "sample",
    treatment = "D",
    cell_type = "ct",
    families = "poisson",
    outdir = tempdir(),
    prefix = "small_structured_pb",
    run_full_without_tuning = TRUE,
    tune_then_run = FALSE,
    n_cores = 1,
    min_samples_per_celltype = 3,
    save_fit = FALSE,
    verbose = FALSE
  ))

  expect_equal(out$full_metrics$n_tests, 20)
  expect_equal(out$selected_grid$variant, "fold5_n150_lr003_leaves7_min5")
})
