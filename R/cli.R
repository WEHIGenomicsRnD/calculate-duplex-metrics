# ------------------------------------------------------------------
# cli.R
#
# Command-line interface for calculating duplex sequencing metrics.
#
# Defines the main CLI entrypoint, parses and validates command-line
# arguments, resolves input files, and delegates computation to
# `process_data()` in calculate.R.
#
# When the package is installed, this file is part of the package
# namespace. When running locally from the project root (renv), it
# is sourced by main.R.
# ------------------------------------------------------------------


# Main CLI function
#' @export
main <- function() {
  # Parse arguments
  args <- parse_arguments()

  # Run metric computation
  res <- process_data(
    input     = args$input,
    output    = args$output,
    sample    = args$sample,
    rlen      = args$rlen,
    skips     = args$skips,
    ref_fasta = args$ref_fasta,
    metrics   = args$metrics,
    cores     = args$cores,
    input_format = args$input_format
  )


  # Exit
  if (is.list(res) && isTRUE(res$success)) {
    cat("Process completed successfully\n")
    quit(status = 0)
  } else {
    cat(
      "Process failed:",
      if (is.list(res)) res$error else "unknown error",
      "\n"
    )
    quit(status = 1)
  }
}

# Parse command line arguments
parse_arguments <- function() {
  p <- argparse::ArgumentParser(description = "Calculate duplex QC metrics")

  p$add_argument(
    "-i", "--input", nargs = "+", required = TRUE,
    help = paste(
      "Input path(s): either a directory OR one/more metric input files.",
      "You can also pass a comma-separated list.",
      "Examples:",
      "--input data/",
      "--input a.txt b.txt.gz",
      "--input a.txt,b.txt.gz"
    )
  )
  p$add_argument(
    "--pattern", default = "\\.txt(\\.gz)?$",
    help = paste("Regex pattern used only when --input is a directory",
                 "(default: \\.txt(\\.gz)?$)")
  )

  p$add_argument("--cores", type = "integer", default = 1,
                 help = "Number of cores for parallel processing (default: 1)")

  p$add_argument("-o", "--output", required = TRUE,
                 help = "Output CSV path (long format: sample,metric,value)")
  p$add_argument("-s", "--sample", default = NULL,
                 help = paste("Optional sample name(s). For multiple inputs, ",
                              "provide comma-separated names matching the",
                              "number of files."))
  p$add_argument("--rlen", type = "integer", default = 151,
                 help = "Read length (default: 151)")
  p$add_argument("--skips", type = "integer", default = 5,
                 help = paste("Trimmed/ignored bases per read",
                              "(NanoSeq = 5, xGen = 8)"))
  p$add_argument("--ref_fasta", default = "",
                 help = paste("Optional reference genome object (.fasta file).",
                              "GC metrics are skipped if not provided."))
  p$add_argument(
    "--input_format",
    default = "rinfo",
    choices = .supported_input_formats,
    help = paste(
      "Input format (default: rinfo).",
      "Use fgbio for",
      "CollectDuplexSeqMetrics duplex_family_sizes tables."
    )
  )
  p$add_argument(
    "--metrics", default = "all",
    help = paste(
      "Comma-separated metric groups or metric names. Default: all.",
      "Groups: gc,family.",
      "Individual: frac_singletons,efficiency,drop_out_rate.",
      "You may also specify grouped metrics (e.g. gc_single or family_mean)",
      "and the full group will be computed"
    )
  )

  args <- p$parse_args()

  if (!is.null(args$metrics) && nzchar(args$metrics)) {
    args$metrics <- tolower(gsub("\\s+", "", args$metrics))
  }
  args$input_format <- normalize_input_format(args$input_format)

  # validate cores
  if (is.na(args$cores) || args$cores < 1) {
    stop("--cores must be >= 1")
  }

  # ---- consolidated input resolution ----
  # flatten comma-separated items and trim whitespace
  raw_inputs <- unlist(strsplit(paste(args$input, collapse = ","), ",",
                                fixed = TRUE))
  raw_inputs <- trimws(raw_inputs)
  raw_inputs <- raw_inputs[nzchar(raw_inputs)]

  if (length(raw_inputs) == 0) {
    stop("You must provide at least one --input (directory or file path).")
  }

  # detect directory mode
  dir_inputs <- raw_inputs[dir.exists(raw_inputs)]
  file_inputs <- setdiff(raw_inputs, dir_inputs)

  if (length(dir_inputs) > 1) {
    stop(paste("Please provide only one input directory, or provide explicit",
               "file paths (not multiple directories)."))
  }

  using_dir <- length(dir_inputs) == 1

  if (using_dir && length(file_inputs) > 0) {
    stop(paste("Please provide either a single input directory OR explicit",
               "file paths, not both."))
  }

  if (using_dir) {
    if (!is.null(args$sample) && nzchar(args$sample)) {
      stop(paste("Do not use --sample when --input is a directory.",
                 "File ordering from a directory may not match the provided",
                 "sample name order."))
    }

    in_dir <- dir_inputs[[1]]
    files <- sort(list.files(in_dir, pattern = args$pattern, full.names = TRUE))

    if (length(files) == 0) {
      stop("No input files found in directory: ", in_dir,
           " (pattern: ", args$pattern, ")")
    }

    args$input <- files
  } else {
    # if directory not used, treat as explicit files
    args$input <- file_inputs
  }

  # validate that all input files exist
  missing <- args$input[!file.exists(args$input)]
  if (length(missing) > 0) {
    stop("Input file(s) not found:\n", paste(missing, collapse = "\n"))
  }

  # normalise paths
  args$input <- normalizePath(args$input, mustWork = TRUE)

  # validate --sample only for explicit file input (dir case already blocked)
  if (!is.null(args$sample) && nzchar(args$sample)) {
    sample_names <- trimws(unlist(strsplit(args$sample, ",", fixed = TRUE)))
    sample_names <- sample_names[nzchar(sample_names)]

    if (length(sample_names) != length(args$input)) {
      stop("--sample must contain exactly ", length(args$input),
           " name(s) to match the number of input files.")
    }

    args$sample <- sample_names
  }

  # validate numeric CLI arguments early to avoid invalid metric computation
  if (is.na(args$rlen) || args$rlen <= 0)  stop("--rlen must be positive")
  if (is.na(args$skips) || args$skips < 0) stop("--skips must be >= 0")
  if (args$skips >= args$rlen) stop("--skips must be < --rlen")

  # output dir
  odir <- dirname(args$output)
  if (!dir.exists(odir))
    dir.create(odir, recursive = TRUE, showWarnings = FALSE)

  args
}
