library(testthat)
library(data.table)

# ------------------------------------------------------------------------------
# Fixtures
# ------------------------------------------------------------------------------

bam_fixture_path <- normalizePath(
  test_path("testdata", "NanoMB1Rep1_sub0-001.bam"),
  mustWork = TRUE
)
reference_rinfo_path <- normalizePath(
  test_path("testdata", "NanoMB1Rep1_sub0-001_rinfo.txt.gz"),
  mustWork = TRUE
)

.sorted_rows <- function(dt) {
  dt <- as.data.table(dt)
  setorderv(dt, names(dt))
  dt
}

# ------------------------------------------------------------------------------
# generate_read_info_from_bam
# ------------------------------------------------------------------------------

test_that("generate_read_info_from_bam reproduces the legacy rinfo table exactly", {
  skip_if_not(nzchar(Sys.which("sort")), "GNU sort not available")

  out_file <- withr::local_tempfile(fileext = ".txt.gz")

  res <- generate_read_info_from_bam(
    bam_file = bam_fixture_path,
    output = out_file,
    chunk_size = 500000,
    cores = 1,
    sort_mem = "1G",
    tmp_dir = tempdir()
  )

  expect_true(file.exists(res$output))
  expect_gt(res$n_bundles, 0)

  generated <- .sorted_rows(fread(res$output))
  expected  <- .sorted_rows(fread(reference_rinfo_path))

  expect_equal(nrow(generated), nrow(expected))
  expect_equal(generated, expected)
})

test_that("generate_read_info_from_bam is robust to chunk_size choice", {
  skip_if_not(nzchar(Sys.which("sort")), "GNU sort not available")

  out_small <- withr::local_tempfile(fileext = ".txt.gz")
  out_large <- withr::local_tempfile(fileext = ".txt.gz")

  res_small <- generate_read_info_from_bam(
    bam_file = bam_fixture_path,
    output = out_small,
    chunk_size = 5000,
    cores = 1,
    sort_mem = "1G",
    tmp_dir = tempdir()
  )
  res_large <- generate_read_info_from_bam(
    bam_file = bam_fixture_path,
    output = out_large,
    chunk_size = 5000000,
    cores = 1,
    sort_mem = "1G",
    tmp_dir = tempdir()
  )

  expect_equal(res_small$n_bundles, res_large$n_bundles)
  expect_equal(
    .sorted_rows(fread(out_small)),
    .sorted_rows(fread(out_large))
  )
})

test_that("generate_read_info_from_bam errors on a missing BAM file", {
  expect_error(
    generate_read_info_from_bam(
      bam_file = "does_not_exist.bam",
      tmp_dir = tempdir()
    )
  )
})

# ------------------------------------------------------------------------------
# process_data end-to-end with --input_format bam
# ------------------------------------------------------------------------------

test_that("process_data computes metrics directly from a BAM file", {
  skip_if_not(nzchar(Sys.which("sort")), "GNU sort not available")

  out_csv <- withr::local_tempfile(fileext = ".csv")

  res <- process_data(
    input = bam_fixture_path,
    output = out_csv,
    rlen = 151,
    skips = 5,
    metrics = "efficiency,drop_out_rate,family",
    input_format = "bam",
    sort_mem = "1G",
    bam_chunk_size = 500000
  )

  expect_true(isTRUE(res$success))
  out_tbl <- fread(out_csv)
  expect_true(all(c("efficiency", "drop_out_rate", "total_families") %in%
                    out_tbl$metric))

  # sample name should be derived from the BAM filename, not a temp path
  expect_equal(unique(out_tbl$sample), "NanoMB1Rep1_sub0-001")
})

test_that("process_data with bam input matches rinfo input on the same data", {
  skip_if_not(nzchar(Sys.which("sort")), "GNU sort not available")

  out_bam_csv <- withr::local_tempfile(fileext = ".csv")
  out_rinfo_csv <- withr::local_tempfile(fileext = ".csv")

  res_bam <- process_data(
    input = bam_fixture_path,
    output = out_bam_csv,
    rlen = 151,
    skips = 5,
    metrics = "efficiency,drop_out_rate,family",
    input_format = "bam",
    sort_mem = "1G",
    bam_chunk_size = 500000
  )
  res_rinfo <- process_data(
    input = reference_rinfo_path,
    output = out_rinfo_csv,
    sample = "NanoMB1Rep1_sub0-001",
    rlen = 151,
    skips = 5,
    metrics = "efficiency,drop_out_rate,family",
    input_format = "rinfo"
  )

  expect_true(isTRUE(res_bam$success))
  expect_true(isTRUE(res_rinfo$success))

  bam_tbl <- fread(out_bam_csv)
  rinfo_tbl <- fread(out_rinfo_csv)

  setkey(bam_tbl, metric)
  setkey(rinfo_tbl, metric)
  expect_equal(bam_tbl$value, rinfo_tbl$value)
})

test_that("process_data rejects multiple BAM inputs (phase-1 restriction)", {
  out_csv <- withr::local_tempfile(fileext = ".csv")

  res <- process_data(
    input = c(bam_fixture_path, bam_fixture_path),
    output = out_csv,
    input_format = "bam"
  )

  expect_false(isTRUE(res$success))
  expect_match(res$error, "exactly one", ignore.case = TRUE)
})
