# ------------------------------------------------------------------
# metric_functions.R
#
# Core implementations of duplex-specific metric calculations.
#
# This script defines functions for computing individual and grouped
# duplex sequencing metrics from summarised rinfo data.
# It is not intended to be run as a standalone script.
#
# The functions are sourced by calculate.R and called with explicit
# parameters (e.g. rlen, skips, ref_fasta). They do not rely on
# implicit global variables for configuration.
#
# Metric selection, input validation, and single-/multi-file input
# are handled by higher-level scripts (main.R, cli.R, calculate.R).
#
# Code obtained from:
# https://github.com/WEHIGenomicsRnD/G000204_duplex/blob/main/code/
# efficiency_nanoseq_functions.R
# ------------------------------------------------------------------

suppressPackageStartupMessages({
  library(magrittr)
  library(Rsamtools)
  library(Biostrings)
  library(GenomicRanges)
  library(IRanges)
  library(seqinr)
})
`%>%` <- magrittr::`%>%`


# functions below are adapted from
# R/efficiency_nanoseq.R and perl/efficiency_nanoseq.pl
# from https://github.com/cancerit/NanoSeq

# Fraction of total reads that are from singleton read bundles
calculate_singletons <- function(rbs) {
  total_reads <- sum(rbs$x, rbs$y)
  singletons <- sum(rbs$x == 1 & rbs$y == 0 | rbs$x == 0 & rbs$y == 1)
  singletons / total_reads
}

calculate_family_stats <- function(rbs) {
  rbs$size <- rbs$x + rbs$y
  c(total_families = nrow(rbs),
    family_mean = mean(rbs$size),
    family_median = median(rbs$size),
    family_max = max(rbs$size),
    families_gt1 = sum(rbs$x > 1 | rbs$y > 1),
    single_families = sum(rbs$x == 1 & rbs$y == 0 | rbs$x == 0 & rbs$y == 1),
    paired_families = sum(rbs$x > 0 & rbs$y > 0),
    unpaired_families = sum(rbs$x > 0 & rbs$y == 0 | rbs$x == 0 & rbs$y > 0),
    paired_and_gt1 = sum(rbs$x > 1 & rbs$y > 1),
    frac_singletons = calculate_singletons(rbs)
  )
}


# from cancerit/NanoSeq documentation:
# "This is the number of duplex bases divided by the number of sequenced bases."
calculate_efficiency <- function(rbs, rlen, skips) {
  if (is.na(rlen) || rlen <= 0) stop("rlen must be positive")
  if (is.na(skips) || skips < 0) stop("skips must be >= 0")
  if (skips >= rlen) stop("skips must be < rlen")

  usable_bases_per_bundle <- (as.double(rlen) - as.double(skips)) * 2
  bases_ok_rbs <- as.double(nrow(rbs[rbs$x > 1 & rbs$y > 1, ])) *
    usable_bases_per_bundle
  total_reads <- as.double(sum(c(rbs$x, rbs$y)))
  if (total_reads == 0) return(NA_real_)
  bases_sequenced <- total_reads * as.double(rlen) * 2
  bases_ok_rbs / bases_sequenced
}


# from cancerit/NanoSeq documentation:
# "This shows the fraction of read bundles missing one of the two
# original strands beyond what would be expected under random sampling
# (assuming a binomial process).
calculate_missed_fraction <- function(rbs) {
  rbs <- data.frame(rbs)
  rbs$size <- rbs$x + rbs$y
  rbs$size <- pmin(rbs$size, 10)
  total_missed <- 0
  for (size in c(4:10)) {
    exp_orphan <- (0.5 ** size) * 2
    total_this_size <- nrow(rbs[which(rbs$size == size), ])
    if (total_this_size > 0) {
      with_both_strands <- nrow(
        rbs[which(rbs$size == size & rbs$x > 0 & rbs$y > 0), ]
      )
      obs_orphan <- 1 - with_both_strands / total_this_size
      missed <- (obs_orphan - exp_orphan) * total_this_size
      total_missed <- total_missed + missed
    }
  }
  den <- nrow(rbs[rbs$size >= 4, , drop = FALSE])
  if (den == 0) return(NA_real_)
  total_missed / den
}

calculate_on_target_coverage <- function(rbs, rlen, skips, grx, min_reads) {
  # Validate rlen/skips
  if (is.na(rlen) || rlen <= 0) stop("rlen must be positive")
  if (is.na(skips) || skips < 0) stop("skips must be >= 0")
  if (skips >= rlen) stop("skips must be < rlen")

  # Validate rbs input
  required_cols <- c("chrom", "pos", "mpos", "x", "y")
  missing_cols <- setdiff(required_cols, names(rbs))
  if (length(missing_cols) > 0) {
    stop(
      "calculate_coverage requires columns: ",
      paste(required_cols, collapse = ", "),
      ". Missing: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  # Handle min_reads input
  if (is.null(min_reads)) {
    min_reads <- c(4L, 2L, 2L)
  }
  if (is.character(min_reads)) {
    tokens <- unlist(strsplit(gsub(",", " ", min_reads), "\\s+"))
    tokens <- tokens[nzchar(tokens)]
    min_reads <- as.integer(tokens)
  } else {
    min_reads <- as.integer(min_reads)
  }
  if (length(min_reads) != 3 || anyNA(min_reads) || any(min_reads < 0)) {
    stop("--min_reads must contain exactly three non-negative integers.")
  }

  # epos is used only for overlap detection; usable bases are always
  # 2*(rlen-skips) per bundle (one read pair = forward + reverse read)
  epos_available <- "epos" %in% names(rbs)
  epos <- if (epos_available) rbs$epos else rbs$pos + ((rlen - skips) * 2)

  usable_per_bundle <- as.double(rlen - skips) * 2.0

  # Total bases in the target capture region
  target_area <- as.double(sum(GenomicRanges::width(grx)))
  if (target_area == 0) {
    return(c(raw_coverage = NA_real_, duplex_coverage = NA_real_))
  }

  # Find on-target read bundles (any overlap with target regions)
  rbx <- GenomicRanges::GRanges(
    seqnames = rbs$chrom,
    ranges = IRanges::IRanges(start = rbs$pos, end = epos)
  )
  on_target_idx <- unique(
    S4Vectors::queryHits(GenomicRanges::findOverlaps(rbx, grx))
  )
  # Raw coverage: each read pair (x + y) contributes usable_per_bundle bases
  raw_bases_on_target <- usable_per_bundle *
    as.double(sum(rbs$x[on_target_idx] + rbs$y[on_target_idx]))
  raw_coverage <- raw_bases_on_target / target_area

  # Duplex coverage: apply min_reads filter, then each passing bundle
  # collapses x/y into one duplex fragment contributing usable_per_bundle bases
  ss1 <- rbs$x >= min_reads[[2]] & rbs$y >= min_reads[[3]]
  ss2 <- rbs$x >= min_reads[[3]] & rbs$y >= min_reads[[2]]
  duplex_pass <- rbs$x + rbs$y >= min_reads[[1]] & (ss1 | ss2)

  n_duplex_on_target <- as.double(
    length(intersect(which(duplex_pass), on_target_idx))
  )
  duplex_coverage <- (n_duplex_on_target * usable_per_bundle) / target_area

  c(
    on_target_coverage_raw = raw_coverage,
    on_target_coverage_duplex = duplex_coverage,
    on_target_duplex_ratio = raw_coverage / duplex_coverage
  )
}

# Calculate on-target rate from read info and a GRanges object
# corresponding to target regions. It is approximate because we
# don't know the actual end of the read from the rinfo file.
# Uses epos (actual end position) if available, otherwise falls
# back to mpos + rlen*2.
calculate_on_target_rate_raw <- function(rbs, grx, rlen) {
  required_cols <- c("chrom", "pos", "mpos", "x", "y")
  missing_cols <- setdiff(required_cols, names(rbs))
  if (length(missing_cols) > 0) {
    stop(
      "on_target_rate_raw requires columns: ",
      paste(required_cols, collapse = ", "),
      ". Missing: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  total_reads <- (sum(rbs$x) + sum(rbs$y)) * 2
  total_reads <- as.double(total_reads)
  if (total_reads == 0) return(NA_real_)

  # Use epos (actual end position) if available, otherwise estimate
  # from mpos + rlen*2
  end_pos <- if ("epos" %in% names(rbs)) {
    rbs$epos
  } else {
    rbs$mpos + (rlen * 2)
  }

  rbx <- GenomicRanges::GRanges(seqnames = rbs$chrom,
                                ranges = IRanges(start = rbs$pos,
                                                 end = end_pos))
  rbs_olaps <- GenomicRanges::findOverlaps(rbx, grx) |>
    S4Vectors::queryHits() |>
    unique()
  on_target_reads <- (
    sum(rbs[rbs_olaps, ]$x) + sum(rbs[rbs_olaps, ]$y)
  ) * 2
  as.double(on_target_reads) / total_reads
}

# Calculate on-target rate from read info using duplex min_read filtering.
# Approximate calculation only: this estimates duplex consensus  reads from
# bundle counts, but those reads may be further filtered downstream.
calculate_on_target_rate_dup <- function(rbs, grx, rlen, min_reads) {
  required_cols <- c("chrom", "pos", "mpos", "x", "y")
  missing_cols <- setdiff(required_cols, names(rbs))
  if (length(missing_cols) > 0) {
    stop(
      "on_target_rate_duplex requires columns: ",
      paste(required_cols, collapse = ", "),
      ". Missing: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  if (is.null(min_reads)) {
    min_reads <- c(4L, 2L, 2L)
  }

  if (is.character(min_reads)) {
    tokens <- unlist(strsplit(gsub(",", " ", min_reads), "\\s+"))
    tokens <- tokens[nzchar(tokens)]
    min_reads <- as.integer(tokens)
  } else {
    min_reads <- as.integer(min_reads)
  }
  if (length(min_reads) != 3 || anyNA(min_reads) || any(min_reads < 0)) {
    stop("--min_reads must contain exactly three non-negative integers.")
  }

  # Test single-strand min read criteria on either strand
  ss1 <- rbs$x >= min_reads[[2]] & rbs$y >= min_reads[[3]]
  ss2 <- rbs$x >= min_reads[[3]] & rbs$y >= min_reads[[2]]
  duplex_rbs <- rbs[rbs$x + rbs$y >= min_reads[[1]] & (ss1 | ss2), ,
                    drop = FALSE]

  calculate_on_target_rate_raw(duplex_rbs, grx, rlen)
}

# Backward-compatible alias for existing callers.
calculate_on_target_rate <- calculate_on_target_rate_raw

# from cancerit/NanoSeq documentation:
# The GC content of RBs with both strands and with just one strand.
# I return the difference between the two values.
calculate_gc <- function(
  rbs,
  rlen,
  skips,
  genome_file,
  genome_max,
  sample_n = 10000,
  max_gap = 100000
) {
  rbs <- data.frame(rbs)

  # Rename x/y columns to plus/minus for consistency
  # Handle both old positional naming and new explicit naming
  if ("x" %in% names(rbs) && "y" %in% names(rbs)) {
    names(rbs)[names(rbs) == "x"] <- "plus"
    names(rbs)[names(rbs) == "y"] <- "minus"
  } else if (length(names(rbs)) >= 6) {
    # Fallback to positional if x/y not found (for compatibility)
    names(rbs)[5:6] <- c("plus", "minus")
  }

  if (is.null(genome_max) || length(genome_max) == 0) {
    warning(paste0(
      "calculate_gc: genome_max is NULL or empty - ",
      "no chromosome lengths available. Returning NA."
    ))
    return(c(gc_single = NA_real_, gc_both = NA_real_, gc_deviation = NA_real_))
  }

  n_before <- nrow(rbs)
  # remove any chroms not in the sizes vector
  rbs <- rbs[rbs$chrom %in% names(genome_max), ]
  if (nrow(rbs) == 0) {
    warning(paste0(
      "calculate_gc: no records remain after filtering for ",
      "chromosomes present in ref_fasta. ",
      "Check that chromosome names in the input file ",
      "match those in the reference FASTA."
    ))
    return(c(gc_single = NA_real_, gc_both = NA_real_, gc_deviation = NA_real_))
  }
  if (nrow(rbs) < n_before) {
    warning(paste0(
      sprintf(
        "calculate_gc: %d/%d records removed - ",
        n_before - nrow(rbs), n_before
      ),
      "chromosome not found in ref_fasta. ",
      "Check chromosome name formatting ",
      "between input and reference."
    ))
  }

  # compute end and drop invalid ranges early
  # Use epos if available, otherwise compute from pos + rlen - skips
  if (!("epos" %in% names(rbs))) {
    rbs$epos <- rbs$pos + rlen - skips
  }

  n_before <- nrow(rbs)
  rbs <- rbs[!is.na(rbs$pos) &
               !is.na(rbs$epos) & rbs$pos > 0 & rbs$epos >= rbs$pos, ]
  if (nrow(rbs) == 0) {
    warning(paste0(
      "calculate_gc: no records remain after removing ",
      "invalid/NA positions. Returning NA."
    ))
    return(c(gc_single = NA_real_, gc_both = NA_real_, gc_deviation = NA_real_))
  }
  if (nrow(rbs) < n_before) {
    warning(paste0(
      sprintf(
        "calculate_gc: %d/%d records removed due to ",
        n_before - nrow(rbs), n_before
      ),
      "invalid or NA position values."
    ))
  }

  # remove records with mate positions exceeding genome max
  n_before <- nrow(rbs)
  if (length(genome_max) > 1) {
    for (i in seq_along(genome_max)) {
      chrom <- names(genome_max)[i]
      chrom_max <- genome_max[[i]]
      rbs <- rbs[!(rbs$chrom == chrom & rbs$epos > chrom_max), ]
    }
  } else {
    rbs <- rbs[rbs$epos <= genome_max, ]
  }
  if (nrow(rbs) == 0) {
    warning(paste0(
      "calculate_gc: no records remain after removing reads ",
      "that extend beyond chromosome ends. ",
      "Check that --rlen and --skips are correct ",
      "for this data."
    ))
    return(c(gc_single = NA_real_, gc_both = NA_real_, gc_deviation = NA_real_))
  }
  if (nrow(rbs) < n_before) {
    warning(paste0(
      sprintf(
        "calculate_gc: %d/%d records removed - ",
        n_before - nrow(rbs), n_before
      ),
      "read end exceeds chromosome length. ",
      "Consider checking --rlen and --skips."
    ))
  }

  # remove records with large distances between mates
  n_before <- nrow(rbs)
  rbs <- rbs[(rbs$epos - rbs$pos) < max_gap, ]
  if (nrow(rbs) == 0) {
    warning(paste0(
      sprintf(
        "calculate_gc: no records remain after removing ",
        "reads with gap >= %d. ",
        max_gap
      ),
      "Returning NA."
    ))
    return(c(gc_single = NA_real_, gc_both = NA_real_, gc_deviation = NA_real_))
  }
  if (nrow(rbs) < n_before) {
    warning(paste0(
      sprintf(
        "calculate_gc: %d/%d records removed - ",
        n_before - nrow(rbs), n_before
      ),
      sprintf(
        "gap between read start and end >= %d ",
        max_gap
      ),
      "(max_gap)."
    ))
  }

  # remove zero-records (redundant but safe)
  rbs <- rbs[rbs$pos != 0, ]

  # split RBs
  rbs_both <- rbs[which(rbs$minus + rbs$plus >= 4 &
                          rbs$minus >= 2 &
                          rbs$plus  >= 2), ]
  rbs_single <- rbs[which(rbs$minus + rbs$plus > 4 &
                            (rbs$minus == 0 | rbs$plus == 0)), ]

  # guard: not enough RBs
  if (nrow(rbs_both) == 0 || nrow(rbs_single) == 0) {
    warning(paste0(
      sprintf(
        "calculate_gc: not enough read bundles after splitting - ",
        nrow(rbs_both), nrow(rbs_single)
      ),
      sprintf(
        "%d duplex (both-strand) and %d single-strand ",
        nrow(rbs_both), nrow(rbs_single)
      ),
      "bundles found. ",
      "GC metrics require both duplex and ",
      "single-strand bundles. Returning NA."
    ))
    return(c(gc_single = NA_real_, gc_both = NA_real_, gc_deviation = NA_real_))
  }

  # sample
  rbs_both <- rbs_both[sample(seq_len(nrow(rbs_both)),
                              min(sample_n, nrow(rbs_both))), ]
  rbs_single <- rbs_single[sample(seq_len(nrow(rbs_single)),
                                  min(sample_n, nrow(rbs_single))), ]

  # extract sequences
  seqs_both <- scanFa(
    genome_file,
    GRanges(
      rbs_both$chrom,
      IRanges(start = rbs_both$pos, end = rbs_both$epos)
    )
  ) |> as.vector()

  seqs_single <- scanFa(
    genome_file,
    GRanges(
      rbs_single$chrom,
      IRanges(start = rbs_single$pos, end = rbs_single$epos)
    )
  ) |> as.vector()

  # drop NA sequences
  seqs_both   <- seqs_both[!is.na(seqs_both)]
  seqs_single <- seqs_single[!is.na(seqs_single)]

  if (length(seqs_both) == 0 || length(seqs_single) == 0) {
    warning(paste0(
      sprintf(
        "calculate_gc: no sequences extracted from ref_fasta ",
        "(%d duplex, %d single-strand). ",
        length(seqs_both), length(seqs_single)
      ),
      "Check that the reference FASTA is indexed ",
      "and chromosome names match."
    ))
    return(c(gc_single = NA_real_, gc_both = NA_real_, gc_deviation = NA_real_))
  }

  seqs_both_collapsed   <- paste(seqs_both, collapse = "")
  seqs_single_collapsed <- paste(seqs_single, collapse = "")

  if (is.na(seqs_both_collapsed) || is.na(seqs_single_collapsed) ||
        nchar(seqs_both_collapsed) == 0 || nchar(seqs_single_collapsed) == 0) {
    warning(paste0(
      "calculate_gc: collapsed sequence strings are empty ",
      "or NA after extraction. Returning NA."
    ))
    return(c(gc_single = NA_real_, gc_both = NA_real_, gc_deviation = NA_real_))
  }

  gc_both <- s2c(seqs_both_collapsed) |>
    GC()
  gc_single <- s2c(seqs_single_collapsed) |>
    GC()

  c(gc_single = gc_single,
    gc_both = gc_both,
    gc_deviation = abs(gc_single - gc_both))
}

.supported_input_formats <- c("rinfo", "fgbio", "bam")

normalise_input_format <- function(input_format = "rinfo") {
  input_format <- tolower(trimws(input_format))

  if (!nzchar(input_format) || !input_format %in% .supported_input_formats) {
    stop(
      "Unknown --input_format: ", input_format,
      "\nValid input formats: ",
      paste(.supported_input_formats, collapse = ", ")
    )
  }

  input_format
}

# --- fgbio metric functions ---------
validate_fgbio_family_sizes <- function(metrics_tbl) {
  required_cols <- c("ab_size", "ba_size", "count")
  missing_cols <- setdiff(required_cols, names(metrics_tbl))

  if (length(missing_cols) > 0) {
    stop(
      "fgbio duplex family size input is missing required column(s): ",
      paste(missing_cols, collapse = ", ")
    )
  }

  metrics_tbl <- data.frame(metrics_tbl)
  metrics_tbl <- metrics_tbl[, required_cols, drop = FALSE]

  for (col in required_cols) {
    if (anyNA(metrics_tbl[[col]])) {
      stop("fgbio duplex family size input contains NA values in ", col)
    }
    if (any(metrics_tbl[[col]] < 0)) {
      stop("fgbio duplex family size input contains negative values in ", col)
    }
  }

  metrics_tbl
}

calculate_weighted_median <- function(values, weights) {
  if (length(values) == 0 || length(weights) == 0 || sum(weights) == 0) {
    return(NA_real_)
  }

  ord <- order(values)
  values <- values[ord]
  weights <- weights[ord]
  cutoff <- sum(weights) / 2

  values[which(cumsum(weights) >= cutoff)[1]]
}

calculate_family_stats_fgbio <- function(metrics_tbl) {
  metrics_tbl <- validate_fgbio_family_sizes(metrics_tbl)

  family_size <- metrics_tbl$ab_size + metrics_tbl$ba_size
  total_families <- sum(metrics_tbl$count)

  c(
    total_families = total_families,
    family_mean = if (total_families == 0) NA_real_
    else sum(family_size * metrics_tbl$count) / total_families,
    family_median = calculate_weighted_median(family_size, metrics_tbl$count),
    family_max = if (total_families == 0) NA_real_ else max(family_size),
    families_gt1 = sum(metrics_tbl$count[
      metrics_tbl$ab_size > 1 | metrics_tbl$ba_size > 1
    ]),
    single_families = sum(metrics_tbl$count[
      (metrics_tbl$ab_size == 1 & metrics_tbl$ba_size == 0) |
        (metrics_tbl$ab_size == 0 & metrics_tbl$ba_size == 1)
    ]),
    paired_families = sum(metrics_tbl$count[
      metrics_tbl$ab_size > 0 & metrics_tbl$ba_size > 0
    ]),
    unpaired_families = sum(metrics_tbl$count[
      (metrics_tbl$ab_size > 0 & metrics_tbl$ba_size == 0) |
        (metrics_tbl$ab_size == 0 & metrics_tbl$ba_size > 0)
    ]),
    paired_and_gt1 = sum(metrics_tbl$count[
      metrics_tbl$ab_size > 1 & metrics_tbl$ba_size > 1
    ]),
    frac_singletons = calculate_singletons_fgbio(metrics_tbl)
  )
}

calculate_singletons_fgbio <- function(metrics_tbl) {
  metrics_tbl <- validate_fgbio_family_sizes(metrics_tbl)

  total_reads <- sum(
    (metrics_tbl$ab_size + metrics_tbl$ba_size) * metrics_tbl$count
  )
  total_reads <- as.double(total_reads)
  if (total_reads == 0) return(NA_real_)

  singletons <- sum(metrics_tbl$count[
    (metrics_tbl$ab_size == 1 & metrics_tbl$ba_size == 0) |
      (metrics_tbl$ab_size == 0 & metrics_tbl$ba_size == 1)
  ])
  as.double(singletons) / total_reads
}

calculate_efficiency_fgbio <- function(metrics_tbl, rlen, skips) {
  metrics_tbl <- validate_fgbio_family_sizes(metrics_tbl)

  if (is.na(rlen) || rlen <= 0) stop("rlen must be positive")
  if (is.na(skips) || skips < 0) stop("skips must be >= 0")
  if (skips >= rlen) stop("skips must be < rlen")

  usable_bases_per_bundle <- (as.double(rlen) - as.double(skips)) * 2
  total_reads <- sum(
    (metrics_tbl$ab_size + metrics_tbl$ba_size) * metrics_tbl$count
  )
  total_reads <- as.double(total_reads)
  if (total_reads == 0) return(NA_real_)

  duplex_families <- as.double(sum(metrics_tbl$count[
    metrics_tbl$ab_size > 1 & metrics_tbl$ba_size > 1
  ]))
  bases_ok_rbs <- duplex_families * usable_bases_per_bundle
  bases_sequenced <- total_reads * as.double(rlen) * 2

  bases_ok_rbs / bases_sequenced
}

calculate_missed_frac_fgbio <- function(metrics_tbl) {
  metrics_tbl <- validate_fgbio_family_sizes(metrics_tbl)

  family_size <- pmin(metrics_tbl$ab_size + metrics_tbl$ba_size, 10)
  total_missed <- 0

  for (size in c(4:10)) {
    exp_orphan <- (0.5 ** size) * 2
    total_this_size <- sum(metrics_tbl$count[family_size == size])
    if (total_this_size > 0) {
      with_both_strands <- sum(metrics_tbl$count[
        family_size == size &
          metrics_tbl$ab_size > 0 &
          metrics_tbl$ba_size > 0
      ])
      obs_orphan <- 1 - with_both_strands / total_this_size
      total_missed <- total_missed + (obs_orphan - exp_orphan) * total_this_size
    }
  }

  den <- sum(metrics_tbl$count[family_size >= 4])
  if (den == 0) return(NA_real_)

  total_missed / den
}

# --- Metric grouping / selection ---------
.individual_metrics <- c("efficiency", "drop_out_rate")

.metric_groups <- list(
  gc = c("gc_single", "gc_both", "gc_deviation"),
  family = c(
    "total_families", "family_mean", "family_median", "family_max",
    "families_gt1", "single_families", "paired_families", "unpaired_families",
    "paired_and_gt1", "frac_singletons"
  )
)


# Resolve --metrics into:
# - groups: grouped metrics to compute (gc/family)
# - individual: other metrics to compute individually (efficiency,
#   drop_out_rate)
#
# Rules:
# - empty / NULL -> compute all available metrics
# - token "gc" or "family" -> compute that whole group
# - token is a individual metric name -> compute only that metric
# - token is a metric inside gc/family -> compute the whole group
resolve_metric_selection <- function(metrics_arg = NULL) {

  metrics_norm <- if (is.null(metrics_arg)) ""
  else tolower(gsub("\\s+", "", metrics_arg))

  # default/all mode
  if (!nzchar(metrics_norm) || identical(metrics_norm, "all")) {
    return(
      list(groups = names(.metric_groups), individual = .individual_metrics)
    )
  }

  tokens <- unlist(strsplit(metrics_norm, ","))
  tokens <- tokens[nzchar(tokens)]

  groups <- character(0)
  individual <- character(0)

  for (tok in tokens) {
    if (tok %in% names(.metric_groups)) {
      groups <- union(groups, tok)
      next
    }
    if (tok %in% .individual_metrics) {
      individual <- union(individual, tok)
      next
    }

    # if metric name inside a group, map to its group
    hit <- names(Filter(function(v) tok %in% v, .metric_groups))
    if (length(hit) > 0) {
      groups <- union(groups, hit)
      next
    }

    stop("Unknown metric/group in --metrics: ", tok,
         "\nValid groups: ", paste(names(.metric_groups), collapse = ", "),
         "\nIndividual metrics: ", paste(.individual_metrics, collapse = ", "),
         "\nGrouped metrics: ",
         paste(unique(unlist(.metric_groups)), collapse = ", "))
  }

  list(groups = groups, individual = individual)
}


# Compute selected metrics (returns 1-row data.frame)
calculate_metrics_selected <- function(
  rbs,
  groups = c("gc", "family"),
  individual = character(0),
  rlen,
  skips,
  input_format = "rinfo",
  genome_file = NULL,
  genome_max = NULL,
  target_regions = NULL,
  min_reads = NULL
) {
  input_format <- normalise_input_format(input_format)
  metrics <- list()

  if (input_format == "fgbio") {
    if ("gc" %in% groups) {
      stop(
        "GC metrics are not supported for ",
        "--input_format fgbio."
      )
    }

    fgbio_tbl <- validate_fgbio_family_sizes(rbs)

    if ("efficiency" %in% individual) {
      metrics$efficiency <- calculate_efficiency_fgbio(
        fgbio_tbl,
        rlen = rlen,
        skips = skips
      )
    }
    if ("drop_out_rate" %in% individual) {
      metrics$drop_out_rate <- calculate_missed_frac_fgbio(fgbio_tbl)
    }

    if ("family" %in% groups) {
      fam_stats <- calculate_family_stats_fgbio(fgbio_tbl)
      metrics <- c(metrics, as.list(fam_stats))
    }

    if (!is.null(target_regions)) {
      stop("on_target_rate_* is not supported for --input_format fgbio.")
    }

    return(as.data.frame(metrics, check.names = FALSE))
  }

  if (is.null(min_reads)) {
    min_reads <- c(4L, 2L, 2L)
  }

  if ("efficiency" %in% individual) {
    metrics$efficiency <- calculate_efficiency(rbs, rlen = rlen, skips = skips)
  }
  if ("drop_out_rate" %in% individual) {
    metrics$drop_out_rate <- calculate_missed_fraction(rbs)
  }

  if ("gc" %in% groups) {
    if (is.null(genome_file) || is.null(genome_max)) {
      stop(
        "GC metrics requested but required genome objects were not provided. ",
        "Please supply --ref_fasta (or ensure GC is not selected)."
      )
    }
    gc_stats <- calculate_gc(
      rbs,
      rlen = rlen,
      skips = skips,
      genome_file = genome_file,
      genome_max = genome_max
    )
    metrics <- c(metrics, as.list(gc_stats))
  }

  if ("family" %in% groups) {
    fam_stats <- calculate_family_stats(rbs)
    metrics <- c(metrics, as.list(fam_stats))
  }

  if (!is.null(target_regions)) {
    metrics$on_target_rate_raw <- calculate_on_target_rate_raw(
      rbs,
      grx = target_regions,
      rlen = rlen
    )
    metrics$on_target_rate_duplex <- calculate_on_target_rate_dup(
      rbs,
      grx = target_regions,
      rlen = rlen,
      min_reads = min_reads
    )

    # on-target coverage metrics: raw (all reads) and duplex (after min_reads)
    covs <- calculate_on_target_coverage(
      rbs,
      rlen = rlen,
      skips = skips,
      grx = target_regions,
      min_reads = min_reads
    )
    metrics <- c(metrics, as.list(covs))
  }

  as.data.frame(metrics, check.names = FALSE)
}
