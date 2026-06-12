library(testthat)
library(data.table)
library(Rsamtools)
library(GenomicRanges)
library(Biostrings)

# ------------------------------------------------------------------------------
# Load test file
# ------------------------------------------------------------------------------

rinfo <- fread(test_path("testdata", "NanoMB1Rep1_test_10k.txt"))
fgbio_fixture_path <- normalizePath(
  test_path(
    "..",
    "..",
    "data",
    "NanoMB1Rep1.duplex_seq_metrics.duplex_family_sizes.txt"
  ),
  mustWork = TRUE
)
fgbio_family_sizes <- fread(fgbio_fixture_path)

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

test_that("fgbio frac_singletons matches expanded family-size input", {
  fgbio_small <- data.frame(
    ab_size = c(1, 2, 2, 3, 4),
    ba_size = c(0, 0, 2, 1, 4),
    count = c(2, 3, 1, 4, 2)
  )
  expanded <- fgbio_small[rep(seq_len(nrow(fgbio_small)), fgbio_small$count), ]
  expanded <- data.frame(x = expanded$ab_size, y = expanded$ba_size)

  expect_equal(
    calculate_singletons_fgbio(fgbio_small),
    calculate_singletons(expanded)
  )
})

# ------------------------------------------------------------------------------
# calculate_family_stats
# ------------------------------------------------------------------------------

test_that("calculate_family_stats returns named vector with expected names", {
  stats <- calculate_family_stats(rinfo)
  expect_type(stats, "double")
  expect_named(stats, c("total_families", "family_mean", "family_median",
                        "family_max", "families_gt1", "single_families",
                        "paired_families", "paired_and_gt1"))
})

test_that("calculate_family_stats returns correct known values", {
  stats <- calculate_family_stats(rinfo)
  expect_equal(stats[["total_families"]], 9999)
  expect_equal(stats[["family_mean"]], 9, , tolerance = 1e-3)
  expect_equal(stats[["family_max"]], 42)
  expect_equal(stats[["single_families"]], 1943)
  expect_equal(stats[["paired_families"]], 6008)
})

test_that("fgbio family stats match expanded family-size input", {
  fgbio_small <- data.frame(
    ab_size = c(1, 2, 2, 3, 4),
    ba_size = c(0, 0, 2, 1, 4),
    count = c(2, 3, 1, 4, 2)
  )
  expanded <- fgbio_small[rep(seq_len(nrow(fgbio_small)), fgbio_small$count), ]
  expanded <- data.frame(x = expanded$ab_size, y = expanded$ba_size)

  expect_equal(
    calculate_family_stats_fgbio(fgbio_small),
    calculate_family_stats(expanded)
  )
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

test_that("calculate_efficiency avoids integer overflow on large counts", {
  large_rinfo <- data.frame(
    x = c(20000000L, 1L),
    y = c(20000000L, 0L)
  )

  expect_no_warning(
    eff <- calculate_efficiency(large_rinfo, rlen = 151L, skips = 5L)
  )
  expect_type(eff, "double")
  expect_true(is.finite(eff))
  expect_equal(eff, 292 / (40000001 * 151 * 2), tolerance = 1e-15)
})

test_that("fgbio efficiency matches expanded family-size input", {
  fgbio_small <- data.frame(
    ab_size = c(1, 2, 2, 3, 4),
    ba_size = c(0, 0, 2, 1, 4),
    count = c(2, 3, 1, 4, 2)
  )
  expanded <- fgbio_small[rep(seq_len(nrow(fgbio_small)), fgbio_small$count), ]
  expanded <- data.frame(x = expanded$ab_size, y = expanded$ba_size)

  expect_equal(
    calculate_efficiency_fgbio(fgbio_small, rlen = rlen, skips = skips),
    calculate_efficiency(expanded, rlen = rlen, skips = skips)
  )
})

test_that("calculate_efficiency_fgbio avoids integer overflow", {
  total_reads <- sum(
    (as.double(fgbio_family_sizes$ab_size) +
       as.double(fgbio_family_sizes$ba_size)) *
      as.double(fgbio_family_sizes$count)
  )
  duplex_families <- sum(
    as.double(fgbio_family_sizes$count[
      fgbio_family_sizes$ab_size > 1 &
        fgbio_family_sizes$ba_size > 1
    ])
  )
  expected <- (duplex_families * ((151 - 5) * 2)) / (total_reads * 151 * 2)

  expect_no_warning(
    eff <- calculate_efficiency_fgbio(
      fgbio_family_sizes,
      rlen = 151L,
      skips = 5L
    )
  )
  expect_type(eff, "double")
  expect_true(is.finite(eff))
  expect_equal(eff, expected, tolerance = 1e-15)
})
# ------------------------------------------------------------------------------
# calculate_missed_fraction
# ------------------------------------------------------------------------------

test_that("calculate_missed_fraction returns numeric or NA", {
  val <- calculate_missed_fraction(rinfo)
  expect_type(val, "double")
  expect_equal(val, 0.192, tolerance = 1e-3)
})

test_that("calculate_missed_fraction returns NA with empty input", {
  val <- calculate_missed_fraction(rinfo_empty)
  expect_true(is.na(val))
})

test_that("fgbio drop-out rate matches expanded family-size input", {
  fgbio_small <- data.frame(
    ab_size = c(1, 2, 2, 3, 4),
    ba_size = c(0, 0, 2, 1, 4),
    count = c(2, 3, 1, 4, 2)
  )
  expanded <- fgbio_small[rep(seq_len(nrow(fgbio_small)), fgbio_small$count), ]
  expanded <- data.frame(x = expanded$ab_size, y = expanded$ba_size)

  expect_equal(
    calculate_missed_frac_fgbio(fgbio_small),
    calculate_missed_fraction(expanded)
  )
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
      genome_file = NULL,
      genome_max = NULL
    ),
    "genome_max is NULL or empty"
  )
  expect_named(gc, c("gc_single", "gc_both", "gc_deviation"))
  expect_true(all(is.na(gc)))
})


test_that("calculate_gc returns expected GC metrics when reference provided", {
  ref <- withr::local_tempfile(fileext = ".gz")
  download.file(
    "https://sra-download.ncbi.nlm.nih.gov/traces/wgs03/wgs_aux/NA/RG/NARG01/NARG01.1.fsa_nt.gz",  # nolint
    destfile = ref, mode = "wb")
  fa <- R.utils::gunzip(ref, remove = FALSE)
  indexFa(fa)
  genome_file <- FaFile(fa)
  genome_max <- seqlengths(genome_file)

  gc <- suppressWarnings(
    calculate_gc(rinfo, rlen = rlen, skips = skips, genome_file = genome_file,
                 genome_max = genome_max)
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
  expect_true(all(c("gc", "family") %in% sel$groups))
  expect_true(all(c("frac_singletons", "efficiency", "drop_out_rate") %in%
                    sel$individual))
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
    skips = skips 
  )
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 1)
  expect_true("frac_singletons" %in% colnames(res))
  expect_true("family_mean" %in% colnames(res))
})

test_that("calculate_metrics_selected errors if GC selected w/o ref genome", {
  expect_error(
    calculate_metrics_selected(
      rinfo,
      groups = "gc",
      individual = character(0),
      rlen = rlen,
      skips = skips,
      genome_file = NULL,
      genome_max = NULL
    ),
    "GC metrics requested"
  )
})

test_that("calculate_metrics_selected supports fgbio duplex family size input", {  # nolint
  res <- calculate_metrics_selected(
    fgbio_family_sizes,
    groups = "family",
    individual = c("frac_singletons", "efficiency", "drop_out_rate"),
    rlen = rlen,
    skips = skips,
    input_format = "fgbio"
  )

  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 1)
  expect_equal(res$frac_singletons[[1]], 0.0100687768531685, tolerance = 1e-12)
  expect_equal(res$efficiency[[1]], 0.051592899711472, tolerance = 1e-12)
  expect_equal(res$drop_out_rate[[1]], 0.323238117148393, tolerance = 1e-12)
  expect_equal(res$total_families[[1]], 18297772)
  expect_equal(res$family_mean[[1]], 8.693238061989186, tolerance = 1e-12)
  expect_equal(res$family_median[[1]], 8)
  expect_equal(res$family_max[[1]], 52)
  expect_equal(res$families_gt1[[1]], 16426257)
  expect_equal(res$single_families[[1]], 1601609)
  expect_equal(res$paired_families[[1]], 9892573)
  expect_equal(res$paired_and_gt1[[1]], 8206722)
})

test_that("calculate_metrics_selected rejects unsupported fgbio GC metrics", {
  expect_error(
    calculate_metrics_selected(
      fgbio_family_sizes,
      groups = "gc",
      individual = character(0),
      rlen = rlen,
      skips = skips,
      input_format = "fgbio"
    ),
    "GC metrics are not supported"
  )
})

test_that("process_data includes frac_singletons for fgbio all metrics", {
  out <- withr::local_tempfile(fileext = ".csv")

  expect_message(
    res <- process_data(
      input = fgbio_fixture_path,
      output = out,
      rlen = rlen,
      skips = skips,
      metrics = "all",
      input_format = "fgbio"
    ),
    "skipped"
  )
  expect_true(isTRUE(res$success))

  out_tbl <- fread(out)
  expect_false(any(out_tbl$metric %in% c(
    "gc_single",
    "gc_both",
    "gc_deviation"
  )))
  expect_true(all(c(
    "frac_singletons",
    "efficiency",
    "drop_out_rate",
    "total_families",
    "family_mean",
    "family_median",
    "family_max",
    "families_gt1",
    "single_families",
    "paired_families",
    "paired_and_gt1"
  ) %in% out_tbl$metric))
})
