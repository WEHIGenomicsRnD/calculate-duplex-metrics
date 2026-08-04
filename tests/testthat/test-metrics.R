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
                        "paired_families", "unpaired_families",
                        "paired_and_gt1", "frac_singletons"))
  expect_equal(tail(names(stats), 1), "frac_singletons")
})

test_that("calculate_family_stats returns correct known values", {
  stats <- calculate_family_stats(rinfo)
  expect_equal(stats[["total_families"]], 9999)
  expect_equal(stats[["family_mean"]], 9, , tolerance = 1e-3)
  expect_equal(stats[["family_max"]], 42)
  expect_equal(stats[["single_families"]], 1943)
  expect_equal(stats[["paired_families"]], 6008)
  # unpaired >= single_families (unpaired allows any count > 0 on one strand)
  expect_true(stats[["unpaired_families"]] >= stats[["single_families"]])
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
# calculate_on_target_rate_raw / calculate_on_target_rate_duplex
# ------------------------------------------------------------------------------

test_that("calculate_on_target_rate_raw and duplex return expected fractions", {
  rbs_small <- data.frame(
    chrom = c("chr1", "chr1"),
    pos = c(100, 500),
    mpos = c(120, 520),
    x = c(2, 3),
    y = c(2, 1)
  )
  grx <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(start = 90, end = 200)
  )

  expect_equal(calculate_on_target_rate_raw(rbs_small, grx, rlen = 100), 0.5)
  expect_equal(
    calculate_on_target_rate_dup(rbs_small, grx, rlen = 100,
                                 min_reads = c(4, 2, 2)),
    1
  )
  expect_equal(
    calculate_on_target_rate_dup(rbs_small, grx, rlen = 100,
                                 min_reads = c(4, 2, 1)),
    0.5
  )
  expect_equal(
    calculate_on_target_rate_dup(rbs_small, grx, rlen = 100,
                                 min_reads = c(4, 1, 2)),
    0.5
  )
})

# ------------------------------------------------------------------------------
# calculate_on_target_coverage
# ------------------------------------------------------------------------------

test_that("calculate_on_target_coverage returns expected coverages", {
  rbs_small <- data.frame(
    chrom = c("chr1", "chr1"),
    pos = c(100, 500),
    mpos = c(120, 520),
    x = c(2, 3),
    y = c(2, 1)
  )
  grx <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(start = 90, end = 200)
  )
  # target width = 111 (90-200 inclusive); rlen=100, skips=0 -> 200 bases/bundle
  # on-target bundle: pos=100 only (pos=500 doesn't overlap [90,200])
  # raw: (x+y) = 2+2 = 4 read-pairs * 200 bases = 800; 800/111
  # duplex (4 2 2): x=2,y=2 -> total=4, x>=2, y>=2 -> PASS; 1 bundle * 200/111
  cov <- calculate_on_target_coverage(rbs_small, rlen = 100, skips = 0,
                                      grx = grx, min_reads = c(4, 2, 2))
  expect_named(cov, c("on_target_coverage_raw", "on_target_coverage_duplex",
                      "on_target_duplex_ratio"))
  expect_equal(as.numeric(cov["on_target_coverage_raw"]),  800 / 111, tolerance = 1e-8)
  expect_equal(as.numeric(cov["on_target_coverage_duplex"]), 200 / 111, tolerance = 1e-8)
  expect_equal(as.numeric(cov["on_target_duplex_ratio"]), 4.0, tolerance = 1e-8)
})

test_that("calculate_on_target_coverage: empty target returns NA", {
  rbs_small <- data.frame(
    chrom = c("chr1"),
    pos = c(100),
    mpos = c(120),
    x = c(3),
    y = c(3)
  )
  grx_empty <- GenomicRanges::GRanges()
  cov <- calculate_on_target_coverage(rbs_small, rlen = 100, skips = 0,
                                      grx = grx_empty, min_reads = c(4, 2, 2))
  expect_true(is.na(cov["on_target_coverage_raw"]))
  expect_true(is.na(cov["on_target_coverage_duplex"]))
})

test_that("calculate_on_target_coverage: no bundles pass duplex filter -> duplex coverage 0", {
  rbs_small <- data.frame(
    chrom = c("chr1"),
    pos = c(100),
    mpos = c(120),
    x = c(1),
    y = c(1)
  )
  grx <- GenomicRanges::GRanges("chr1", IRanges::IRanges(50, 300))
  # total=2 < min_reads[1]=4, so no bundles pass duplex filter
  cov <- calculate_on_target_coverage(rbs_small, rlen = 100, skips = 0,
                                      grx = grx, min_reads = c(4, 2, 2))
  expect_equal(as.numeric(cov["on_target_coverage_duplex"]), 0)
})

test_that("calculate_on_target_coverage: epos used when available", {
  rbs <- data.frame(
    chrom = c("chr1"),
    pos   = c(100),
    mpos  = c(200),
    epos  = c(250),
    x = c(3),
    y = c(3)
  )
  grx <- GenomicRanges::GRanges("chr1", IRanges::IRanges(200, 400))
  # With epos=250, the range [100,250] overlaps [200,400] -> on-target
  cov_with_epos <- calculate_on_target_coverage(rbs, rlen = 100, skips = 0,
                                                grx = grx, min_reads = c(4, 2, 2))
  # Without epos, fallback end = pos + 2*(rlen-skips) = 100 + 200 = 300 -> also overlaps
  rbs_no_epos <- rbs[, c("chrom","pos","mpos","x","y")]
  cov_no_epos <- calculate_on_target_coverage(rbs_no_epos, rlen = 100, skips = 0,
                                              grx = grx, min_reads = c(4, 2, 2))
  # both should find the bundle on-target
  expect_true(cov_with_epos["on_target_coverage_raw"] > 0)
  expect_true(cov_no_epos["on_target_coverage_raw"] > 0)
})

test_that("calculate_on_target_coverage errors on invalid rlen/skips", {
  rbs <- data.frame(chrom="chr1", pos=100, mpos=200, x=2, y=2)
  grx <- GenomicRanges::GRanges("chr1", IRanges::IRanges(50, 300))
  expect_error(calculate_on_target_coverage(rbs, rlen = -1, skips = 0, grx = grx, min_reads = c(4,2,2)))
  expect_error(calculate_on_target_coverage(rbs, rlen = 10,  skips = 10, grx = grx, min_reads = c(4,2,2)))
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
  expect_true(all(c("efficiency", "drop_out_rate") %in% sel$individual))
  expect_false("frac_singletons" %in% sel$individual)
})

test_that("resolve_metric_selection maps frac_singletons to family group", {
  sel <- resolve_metric_selection("frac_singletons")
  expect_equal(sel$groups, "family")
  expect_equal(sel$individual, character(0))
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
    individual = character(0),
    rlen = rlen,
    skips = skips 
  )
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 1)
  expect_true("frac_singletons" %in% colnames(res))
  expect_true("family_mean" %in% colnames(res))
  expect_equal(tail(colnames(res), 1), "frac_singletons")
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
    individual = c("efficiency", "drop_out_rate"),
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
  expect_true(!is.null(res$unpaired_families[[1]]))
  expect_true(res$unpaired_families[[1]] >= res$single_families[[1]])
  expect_equal(res$paired_and_gt1[[1]], 8206722)
  expect_equal(tail(colnames(res), 1), "frac_singletons")
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

test_that("process_data rejects target_bed for fgbio input", {
  out <- withr::local_tempfile(fileext = ".csv")
  bed <- withr::local_tempfile(fileext = ".bed")
  writeLines("chr1\t0\t10", bed)

  res <- process_data(
    input = fgbio_fixture_path,
    output = out,
    rlen = rlen,
    skips = skips,
    metrics = "all",
    target_bed = bed,
    input_format = "fgbio"
  )

  expect_false(isTRUE(res$success))
  expect_match(res$error, "on_target_rate_\\* is not supported")
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
    "unpaired_families",
    "paired_and_gt1"
  ) %in% out_tbl$metric))
  expect_equal(tail(out_tbl$metric, 1), "frac_singletons")
})

test_that("process_data adds on_target rates when --target_bed is provided", {
  in_file <- withr::local_tempfile(fileext = ".txt")
  out <- withr::local_tempfile(fileext = ".csv")
  bed <- withr::local_tempfile(fileext = ".bed")
  fwrite(rinfo, in_file, sep = "\t")
  writeLines("NARG01000001.1\t0\t1000000", bed)

  res <- process_data(
    input = in_file,
    output = out,
    rlen = rlen,
    skips = skips,
    metrics = "efficiency",
    target_bed = bed,
    input_format = "rinfo"
  )

  expect_true(isTRUE(res$success))
  out_tbl <- fread(out)
  expect_true(all(c(
    "on_target_rate_raw",
    "on_target_rate_duplex",
    "on_target_coverage_raw",
    "on_target_coverage_duplex",
    "on_target_duplex_ratio"
  ) %in% out_tbl$metric))
})
