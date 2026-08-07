# ------------------------------------------------------------------
# calculate.R
#
# This script is the main orchestration layer for duplex sequencing
# metric calculation.
#
# It is invoked upstream via:
#   main.R -> cli.R -> calculate.R
#
# Responsibilities:
#   - Resolves requested metric groups and individual metrics based on
#     user input
#   - Determines whether GC metrics should be computed based on the
#     requested metrics and the availability of a reference genome
#   - Supports both single-sample and multi-sample inputs
#   - Coordinates serial or parallel execution of metric calculation
#   - Writes results to a tidy CSV with one row per sample–metric pair
#
# When the package is installed, all R/*.R files (including
# metric_functions.R are compiled into the package
# namespace — no source() calls are needed. When running locally from
# the project root via main.R (renv), main.R sources calculate.R
# directly (and also sources cli.R); metric_functions.R is
# sourced explicitly below for that local-only path.
#
# All helper functions in this file focus on coordinating metric
# computation. They do not perform command-line argument parsing;
# inputs are assumed to have been validated upstream (in cli.R).
#
# This script centralises metric selection and GC gating logic so that
# individual metric functions are only called when required.
# ------------------------------------------------------------------


# Import packages
suppressPackageStartupMessages({
  library(parallel)
  library(data.table)
})


# ------------------------------------------------------------------
# Helper functions: no CLI argument parsing or validation;
# inputs are assumed to be validated upstream.
# ------------------------------------------------------------------

parse_samples <- function(sample_arg, inputs) {
  n <- length(inputs)

  # default: derive sample names from filenames
  if (is.null(sample_arg) || length(sample_arg) == 0) {
    return(sub("\\.txt(\\.gz)?$", "", basename(inputs), ignore.case = TRUE))
  }

  # if cli.R already parsed it into a vector (e.g. c("S1","S2"))
  if (length(sample_arg) > 1) {
    samples <- trimws(sample_arg)
    samples <- samples[nzchar(samples)]
  } else {
    # otherwise treat it as a comma-separated string
    if (!nzchar(sample_arg)) {
      return(sub("\\.txt(\\.gz)?$", "", basename(inputs), ignore.case = TRUE))
    }
    samples <- trimws(unlist(strsplit(sample_arg, ",", fixed = TRUE)))
    samples <- samples[nzchar(samples)]
  }

  if (n == 1 && length(samples) == 1) return(samples)

  if (length(samples) != n) {
    stop("--sample must contain the same number of names as --input files.")
  }

  samples
}

read_target_bed <- function(target_bed) {
  target_tbl <- data.table::fread(target_bed, header = FALSE,
                                  data.table = FALSE)
  if (nrow(target_tbl) == 0) {
    stop("Target BED is empty: ", target_bed)
  }
  if (ncol(target_tbl) < 3) {
    stop("Target BED must have at least 3 columns: chrom, start, end")
  }

  target_tbl <- target_tbl[, 1:3, drop = FALSE]
  colnames(target_tbl) <- c("chrom", "start", "end")

  if (anyNA(target_tbl$chrom) || anyNA(target_tbl$start) ||
        anyNA(target_tbl$end)) {
    stop("Target BED contains missing values in required columns.")
  }
  if (any(target_tbl$start < 0) || any(target_tbl$end <= target_tbl$start)) {
    stop(
      "Target BED has invalid coordinates; require start >= 0 and end > start."
    )
  }

  GenomicRanges::reduce(
    GenomicRanges::GRanges(
      seqnames = target_tbl$chrom,
      ranges = IRanges::IRanges(
        start = as.integer(target_tbl$start) + 1L,
        end = as.integer(target_tbl$end)
      )
    )
  )
}

parse_min_reads <- function(min_reads) {
  min_reads_default <- "4 2 2"
  if (is.null(min_reads)) {
    min_reads <- min_reads_default
  } else if (is.character(min_reads) && !nzchar(min_reads)) {
    min_reads <- min_reads_default
  }

  if (is.character(min_reads)) {
    tokens <- unlist(strsplit(gsub(",", " ", min_reads), "\\s+"))
    tokens <- tokens[nzchar(tokens)]
    vals <- as.integer(tokens)
  } else {
    vals <- as.integer(min_reads)
  }
  if (length(vals) != 3 || anyNA(vals) || any(vals < 0)) {
    stop("--min_reads must contain exactly three non-negative integers.")
  }

  vals
}

calc_metrics_one_file_df <- function(
  input,
  sample,
  rlen,
  skips,
  input_format,
  groups_to_compute,
  individual_to_compute,
  genome_file = NULL,
  genome_max = NULL,
  target_regions = NULL,
  min_reads = NULL
) {
  in_file <- normalizePath(input, mustWork = TRUE)

  if (is.null(sample) || !nzchar(sample)) {
    sample <- sub("\\.txt(\\.gz)?$", "", basename(in_file), ignore.case = TRUE)
  }

  rbs <- tryCatch(fread(in_file), error = function(e) e)
  if (inherits(rbs, "error")) stop("Failed to read input: ", rbs$message)

  tbl <- tryCatch( # nolint
    calculate_metrics_selected(
      rbs,
      groups = groups_to_compute,
      individual = individual_to_compute,
      rlen = rlen,
      skips = skips,
      input_format = input_format,
      genome_file = genome_file,
      genome_max = genome_max,
      target_regions = target_regions,
      min_reads = min_reads
    ),
    error = function(e) e
  )
  if (inherits(tbl, "error")) stop("Metric calculation failed: ", tbl$message)
  if (is.null(tbl) || nrow(tbl) == 0) stop("Empty metrics table returned")

  metric_names <- names(tbl)
  data.frame(
    sample = sample,
    metric = metric_names,
    value  = as.numeric(tbl[1, metric_names, drop = TRUE]),
    check.names = FALSE
  )
}

calc_metrics_many_files_df <- function(
  inputs,
  samples,
  rlen,
  skips,
  input_format,
  groups_to_compute,
  individual_to_compute,
  cores = 1,
  genome_file = NULL,
  genome_max = NULL,
  target_regions = NULL,
  min_reads = NULL
) {
  inputs <- normalizePath(inputs, mustWork = TRUE)

  process_one_file <- function(i) {
    calc_metrics_one_file_df(
      input = inputs[i],
      sample = samples[i],
      rlen = rlen,
      skips = skips,
      input_format = input_format,
      groups_to_compute = groups_to_compute,
      individual_to_compute = individual_to_compute,
      genome_file = genome_file,
      genome_max = genome_max,
      target_regions = target_regions,
      min_reads = min_reads
    )
  }

  idx <- seq_along(inputs)

  out_list <- if (cores > 1 && .Platform$OS.type != "windows") {
    parallel::mclapply(idx, process_one_file, mc.cores = cores)
  } else {
    lapply(idx, process_one_file)
  }

  as.data.frame(data.table::rbindlist(out_list, use.names = TRUE, fill = TRUE))
}

# ------------------------------------------------------------------
# process_data resolves metric selection, gates GC logic, and
# dispatches to single- or multi-file helpers. Helpers are pure and
# only compute metrics.
# ------------------------------------------------------------------

process_data <- function(
  input, output,
  sample = NULL,
  rlen = 151,
  skips = 5,
  ref_fasta = "",
  target_bed = "",
  min_reads = "4 2 2",
  metrics = "all",
  cores = 1,
  input_format = "rinfo",
  sort_mem = "4G",
  bam_chunk_size = 2000000L,
  rinfo_output = NULL
) {
  # validate input
  if (is.null(input) || length(input) == 0) {
    return(list(success = FALSE, error = "No input files provided"))
  }

  missing <- input[!file.exists(input)]
  if (length(missing) > 0) {
    return(list(success = FALSE,
                error = paste0("Input file(s) not found:\n",
                               paste(missing, collapse = "\n"))))
  }

  if (is.na(cores) || cores < 1) {
    return(list(success = FALSE, error = "--cores must be >= 1"))
  }

  input_format <- tryCatch(
    normalise_input_format(input_format),
    error = function(e) e
  )
  if (inherits(input_format, "error")) {
    return(list(success = FALSE, error = input_format$message))
  }

  # BAM preprocessing: convert a single BAM into a temporary rinfo table,
  # then continue through the rest of the pipeline as if --input_format
  # rinfo had been used. Phase 1 supports exactly one BAM per invocation.
  bam_tmp_file <- NULL
  if (identical(input_format, "bam")) {
    if (length(input) != 1) {
      return(list(
        success = FALSE,
        error = paste("--input_format bam currently supports exactly one",
                      "BAM file per invocation.")
      ))
    }

    bam_res <- tryCatch(
      generate_read_info_from_bam(
        bam_file = input[1],
        output = if (!is.null(rinfo_output) && nzchar(rinfo_output))
                   rinfo_output else NULL,
        chunk_size = bam_chunk_size,
        cores = cores,
        sort_mem = sort_mem,
        tmp_dir = tempdir()
      ),
      error = function(e) e
    )
    if (inherits(bam_res, "error")) {
      return(list(success = FALSE,
                  error = paste0("BAM preprocessing failed: ",
                                 bam_res$message)))
    }

    bam_tmp_file <- bam_res$output
    # only auto-delete if the user did not supply an explicit rinfo output path
    user_supplied_rinfo <- !is.null(rinfo_output) && nzchar(rinfo_output)
    if (!user_supplied_rinfo) {
      on.exit(if (!is.null(bam_tmp_file)) unlink(bam_tmp_file), add = TRUE)
    }

    # derive a sample name from the original BAM filename before we
    # overwrite `input` with the temporary rinfo path below
    if (is.null(sample) || (is.character(sample) && !nzchar(sample))) {
      sample <- sub("\\.bam$", "", basename(input[1]), ignore.case = TRUE)
    }

    input <- bam_tmp_file
    input_format <- "rinfo"
  }

  odir <- dirname(output)
  if (!dir.exists(odir))
    dir.create(odir, recursive = TRUE, showWarnings = FALSE)

  # Resolve metric selection once
  sel <- resolve_metric_selection(metrics)
  groups_to_compute <- sel$groups
  individual_to_compute <- sel$individual


  if (nzchar(ref_fasta)) {
    ref_fasta <- normalizePath(ref_fasta, mustWork = FALSE)
  }
  if (nzchar(target_bed)) {
    target_bed <- normalizePath(target_bed, mustWork = FALSE)
  }
  min_reads_parsed <- NULL
  if (nzchar(target_bed)) {
    min_reads_parsed <- tryCatch(
      parse_min_reads(min_reads),
      error = function(e) e
    )
    if (inherits(min_reads_parsed, "error")) {
      return(list(success = FALSE, error = min_reads_parsed$message))
    }
  }

  # Consolidated GC gating
  gc_requested <- "gc" %in% groups_to_compute
  has_ref <- nzchar(ref_fasta) && file.exists(ref_fasta)

  metrics_norm <- if (is.null(metrics)) ""
  else tolower(gsub("\\s+", "", metrics))

  is_default_all <- identical(metrics_norm, "") ||
    identical(metrics_norm, "all")
  explicit_metrics <- !is_default_all

  if (gc_requested && identical(input_format, "fgbio")) {
    msg <- paste(
      "GC metrics are not supported for",
      "--input_format fgbio."
    )

    if (explicit_metrics) {
      return(list(success = FALSE, error = msg))
    }

    message("GC metrics skipped for fgbio duplex family size input.")
    groups_to_compute <- setdiff(groups_to_compute, "gc")
  }

  if (nzchar(target_bed) && identical(input_format, "fgbio")) {
    return(list(
      success = FALSE,
      error = "on_target_rate_* is not supported for --input_format fgbio."
    ))
  }

  # If GC is requested but we can't compute it:
  # - default mode (metrics empty): skip GC
  # - explicit mode (user asked for gc/gc_single etc.): fail with clear error
  if (gc_requested && !has_ref) {
    msg <- if (!nzchar(ref_fasta)) {
      paste("GC metrics requested but --ref_fasta not provided.",
            "Please provide --ref_fasta to compute GC.")
    } else {
      paste0("GC metrics requested but --ref_fasta not found at: ", ref_fasta,
             ". Please supply a valid FASTA path.")
    }

    if (explicit_metrics) {
      return(list(success = FALSE, error = msg))
    }

    # default/all mode: skip GC and continue
    message(paste("GC metrics skipped because --ref_fasta was not provided",
                  "or was not found."))
    groups_to_compute <- setdiff(groups_to_compute, "gc")
    ref_fasta <- ""
  }

  if (length(groups_to_compute) == 0 && length(individual_to_compute) == 0 &&
        !nzchar(target_bed)) {
    return(list(success = FALSE, error = "No metrics selected to compute."))
  }

  # Prepare genome objects once (only when GC requested)
  genome_file <- NULL
  genome_max <- NULL
  target_regions <- NULL

  if ("gc" %in% groups_to_compute) {
    genome_file <- Rsamtools::FaFile(ref_fasta)
    fa <- Biostrings::readDNAStringSet(ref_fasta)
    fa_names <- sub("\\s.*$", "", names(fa))
    genome_max <- setNames(as.integer(Biostrings::width(fa)), fa_names)
  }
  if (nzchar(target_bed)) {
    if (!file.exists(target_bed)) {
      return(list(
        success = FALSE,
        error = paste0("Target BED not found at: ", target_bed)
      ))
    }
    target_regions <- tryCatch(
      read_target_bed(target_bed),
      error = function(e) e
    )
    if (inherits(target_regions, "error")) {
      return(list(success = FALSE, error = target_regions$message))
    }
  }


  # Parse samples for single/multi input
  samples <- parse_samples(sample, input)


  out_df <- tryCatch({
    if (length(input) == 1) {
      calc_metrics_one_file_df(
        input = input,
        sample = samples[1],
        rlen = rlen,
        skips = skips,
        input_format = input_format,
        groups_to_compute = groups_to_compute,
        individual_to_compute = individual_to_compute,
        genome_file = genome_file,
        genome_max = genome_max,
        target_regions = target_regions,
        min_reads = min_reads_parsed
      )
    } else {
      calc_metrics_many_files_df(
        inputs = input,
        samples = samples,
        rlen = rlen,
        skips = skips,
        input_format = input_format,
        groups_to_compute = groups_to_compute,
        individual_to_compute = individual_to_compute,
        cores = cores,
        genome_file = genome_file,
        genome_max = genome_max,
        target_regions = target_regions,
        min_reads = min_reads_parsed
      )
    }
  }, error = function(e) e)


  if (inherits(out_df, "error"))
    return(list(success = FALSE, error = out_df$message))

  write_result <- tryCatch({
    write.csv(out_df, output, row.names = FALSE, quote = FALSE)
    TRUE
  }, error = function(e) e)

  if (inherits(write_result, "error"))
    return(list(success = FALSE,
                error = paste0("Write failed: ", write_result$message)))

  if (!is.null(genome_file)) {
    try(close(genome_file), silent = TRUE)
  }
  list(success = TRUE)
}
