library(testthat)
library(data.table)
library(Rsamtools)
library(GenomicRanges)
library(Biostrings)

# ------------------------------------------------------------------------------
# Load test file
# ------------------------------------------------------------------------------

rinfo <- fread(test_path("testdata", "NanoMB1Rep1_test_10k.txt"))

rinfo_empty <- rinfo[0,]
rlen  <- 100
skips <- 0
rinfo[, chrom := "NARG01000001.1"]

# ------------------------------------------------------------------------------
# calculate_singletons
# ------------------------------------------------------------------------------

test_that("calculate_singletons returns fraction in [0,1]", {
  val <- calculate_singletons(rinfo)
  expect_type(val, "double")
  expect_true(val >= 0 && val <= 1)
  expect_equal(val, 0.0216, tolerance = 1e-3)
})

test_that("calculate_singletons returns NA with empty input", {
  val <- calculate_singletons(rinfo_empty)
  expect_true(is.na(val))
})

# ------------------------------------------------------------------------------
# calculate_family_stats
# ------------------------------------------------------------------------------

test_that("calculate_family_stats returns named vector with expected names", {
  stats <- calculate_family_stats(rinfo)
  expect_type(stats, "double")
  expect_named(stats, c("total_families","family_mean","family_median","family_max",
      "families_gt1","single_families","paired_families","paired_and_gt1"))
})

test_that("calculate_family_stats returns correct known values", {
  stats <- calculate_family_stats(rinfo)
  expect_equal(stats[["total_families"]], 9999)
  expect_equal(stats[["family_mean"]], 9, , tolerance = 1e-3)
  expect_equal(stats[["family_max"]], 42)
  expect_equal(stats[["single_families"]], 1943)
  expect_equal(stats[["paired_families"]], 6008)
})

# ------------------------------------------------------------------------------
# calculate_efficiency
# ------------------------------------------------------------------------------

test_that("calculate_efficiency returns valid output", {
  eff <- calculate_efficiency(rinfo, rlen = rlen, skips = skips)
  expect_type(eff, "double")
  expect_true(is.finite(eff))
  expect_true(eff >= 0 && eff <= 1)
  expect_equal(eff, 0.0567, tolerance = 1e-3)
})

test_that("calculate_efficiency errors on invalid rlen/skips", {
  expect_error(calculate_efficiency(rinfo, rlen = -1, skips = 0))
  expect_error(calculate_efficiency(rinfo, rlen = 10, skips = 10))
})

test_that("calculate_efficiency returns NA for zero reads", {
  eff <- calculate_efficiency(rinfo_empty, rlen = rlen, skips = skips)
  expect_true(is.na(eff))
})
# ------------------------------------------------------------------------------
# calculate_missed_fraction
# ------------------------------------------------------------------------------

test_that("calculate_missed_fraction returns numeric or NA", {
  val <- calculate_missed_fraction(rinfo)
  expect_type(val, "double")
  expect_equal(val,0.192, tolerance = 1e-3)
})

test_that("calculate_missed_fraction returns NA with empty input", {
  val <- calculate_missed_fraction(rinfo_empty)
  expect_true(is.na(val))
})

# ------------------------------------------------------------------------------
# calculate_gc
# ------------------------------------------------------------------------------

test_that("calculate_gc returns NA metrics if reference genome missing", {
  expect_warning(
    gc <- calculate_gc(
      rinfo,
      rlen = rlen,
      skips = skips,
      genomeFile = NULL,
      genome_max = NULL
    ),
    "genome_max is NULL or empty"
  )
  expect_named(gc, c("gc_single","gc_both","gc_deviation"))
  expect_true(all(is.na(gc)))
})


test_that("calculate_gc returns expected GC metrics when reference is provided", {  
 ref <- withr::local_tempfile(fileext = ".gz")
  download.file(
    "https://sra-download.ncbi.nlm.nih.gov/traces/wgs03/wgs_aux/NA/RG/NARG01/NARG01.1.fsa_nt.gz",
    destfile = ref, mode = "wb")
 fa <- R.utils::gunzip(ref, remove = FALSE)
 indexFa(fa)
 genomeFile <- FaFile(fa)
 genome_max <- seqlengths(genomeFile)

 gc <- suppressWarnings(
   calculate_gc(rinfo,rlen = rlen, skips = skips, genomeFile = genomeFile, genome_max = genome_max)
 )
  expect_named(gc, c("gc_single", "gc_both", "gc_deviation"))
  expect_equal(gc[["gc_single"]], 0.397, tolerance = 1e-3)
  expect_equal(gc[["gc_both"]], 0.398, tolerance = 1e-3)
  expect_equal(gc[["gc_deviation"]], 0.001, tolerance = 1e-3)
})

# ------------------------------------------------------------------------------
# resolve_metric_selection
# ------------------------------------------------------------------------------

test_that("resolve_metric_selection default selects all", {
  sel <- resolve_metric_selection()
  expect_true(all(c("gc","family") %in% sel$groups))
  expect_true(all(c("frac_singletons","efficiency","drop_out_rate") %in% sel$individual))
})

test_that("resolve_metric_selection errors on unknown metric", {
  expect_error(
    resolve_metric_selection("not_a_metric"),
    "Unknown metric/group"
  )
})

# ------------------------------------------------------------------------------
# calculate_metrics_selected
# ------------------------------------------------------------------------------

test_that("calculate_metrics_selected returns 1-row data.frame", {
  res <- calculate_metrics_selected(
    rinfo,
    groups = "family",
    individual = "frac_singletons",
    rlen = rlen,
    skips = skips )
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 1)
  expect_true("frac_singletons" %in% colnames(res))
  expect_true("family_mean" %in% colnames(res))
})

test_that("calculate_metrics_selected errors if GC selected without refernce genome", {
  expect_error(
    calculate_metrics_selected(
      rinfo,
      groups = "gc",
      individual = character(0),
      rlen = rlen,
      skips = skips,
      genomeFile = NULL,
      genome_max = NULL
    ),
    "GC metrics requested"
  )
})
