# ------------------------------------------------------------------
# read_info_from_bam.R
#
# BAM -> rinfo preprocessing for direct BAM input support.
#
# Converts a coordinate-sorted BAM file into the same summarised
# read-bundle ("rinfo") table produced by the legacy
# scripts/generate_read_info.sh + get_read_info.py + format_read_info.py # nolint
# pipeline, so that existing metric functions (metric_functions.R) can
# be reused unmodified.
#
# Output schema: chrom, pos, mpos, epos, umi, x, y
#   - pos/mpos are the 0-based leftmost/rightmost mate start positions
#     (equivalent to pysam's reference_start for each mate)
#   - epos is the rightmost mate end position (equivalent to pysam's
#     reference_end, computed here from CIGAR reference-consumed width
#     via GenomicAlignments)
#   - x/y are the plus-/minus-strand read-bundle counts
#
# Design (phase 1, single BAM per invocation, pure R):
#   1. Stream the BAM in chunks (Rsamtools::BamFile yieldSize) to bound
#      memory for 100GB+ inputs.
#   2. Within each chunk, pair up mates using vectorised data.table
#      joins on qname (not a per-read loop). Reads whose mate has not
#      yet been seen are carried forward ("pending") to later chunks,
#      mirroring the unbounded mate-cache semantics of the original
#      get_read_info.py (a read whose mate is filtered out or never
#      appears is simply never emitted).
#   3. Extracted rows are appended to a temporary flat file, then
#      aggregated with the same external `sort | cut | uniq -c`
#      pipeline used previously, which keeps global aggregation
#      disk-backed and scalable rather than an in-memory R hash.
#   4. The sorted/counted output is pivoted wide (+/- strand -> x/y)
#      with data.table and written out as the final rinfo table.
#
# This module is invoked by process_data() (calculate.R) when
# --input_format bam is selected; it does not perform CLI parsing.
# ------------------------------------------------------------------

suppressPackageStartupMessages({
  library(Rsamtools)
  library(GenomicAlignments)
  library(data.table)
})

# Build the flag filter used to select usable read pairs, matching the
# semantics of the previous get_read_info.py implementation:
#   - paired, with both mates mapped
#   - not secondary / supplementary / QC-fail
.bam_read_info_flag_filter <- function() {
  Rsamtools::scanBamFlag(
    isPaired = TRUE,
    isUnmappedQuery = FALSE,
    hasUnmappedMate = FALSE,
    isSecondaryAlignment = FALSE,
    isSupplementaryAlignment = FALSE,
    isNotPassingQualityControls = FALSE
  )
}

# Read one chunk of alignments into a flat data.table with only the
# fields required for mate-pairing and UMI/strand extraction.
.bam_chunk_to_dt <- function(aln) {
  m <- S4Vectors::mcols(aln)

  chunk <- data.table::data.table(
    qname      = m$qname,
    chrom      = as.character(GenomeInfoDb::seqnames(aln)),
    mate_chrom = as.character(m$mrnm),
    start0     = BiocGenerics::start(aln) - 1L,
    end0       = BiocGenerics::end(aln),
    is_reverse = as.character(BiocGenerics::strand(aln)) == "-",
    is_read1   = bitwAnd(m$flag, 64L) != 0,
    rx         = m$RX,
    za         = m$ZA,
    zb         = m$ZB
  )

  # keep only mates aligned to the same chromosome
  chunk <- chunk[chrom == mate_chrom] # nolint
  chunk[, mate_chrom := NULL] # nolint
  chunk
}

# Pair up mates within `combined` (pending + new chunk) using a
# vectorised qname join; returns list(paired = <extracted rows>,
# pending = <unmatched rows carried to next chunk>).
.pair_mates <- function(combined) {
  combined[, n_seen := .N, by = qname] # nolint

  extra <- combined[n_seen > 2L] # nolint
  if (nrow(extra) > 0) {
    warning(sprintf(
      "generate_read_info_from_bam: dropping %d record(s) with >2 ",
      "occurrences of the same qname (unexpected after flag filtering).",
      nrow(extra)
    ))
  }

  paired  <- combined[n_seen == 2L] # nolint
  pending <- combined[n_seen == 1L] # nolint
  paired[, n_seen := NULL] # nolint
  pending[, n_seen := NULL] # nolint

  if (nrow(paired) == 0) {
    return(list(
      rows = data.table::data.table(
        qname = character(0), chrom = character(0),
        pos = integer(0), mpos = integer(0), epos = integer(0),
        strand = character(0), umi = character(0)
      ),
      pending = pending
    ))
  }

  r1 <- paired[is_read1 == TRUE] # nolint
  r2 <- paired[is_read1 == FALSE] # nolint
  non_qname <- setdiff(names(r1), "qname")
  data.table::setnames(r1, non_qname, paste0(non_qname, "_r1"))
  data.table::setnames(r2, non_qname, paste0(non_qname, "_r2"))
  merged <- merge(r1, r2, by = "qname")

  merged[, left_is_r1 := start0_r1 <= start0_r2] # nolint
  merged[, strand_out := data.table::fifelse( # nolint
    left_is_r1, # nolint
    data.table::fifelse(!is_reverse_r1 & is_reverse_r2, "+", ""), # nolint
    data.table::fifelse(is_reverse_r1 & !is_reverse_r2, "-", "") # nolint
  )]
  merged <- merged[nzchar(strand_out)] # nolint

  if (nrow(merged) > 0) {
    merged[, umi := data.table::fifelse( # nolint
      strand_out == "-", # nolint
      data.table::fifelse(
        !is.na(za_r1) & !is.na(zb_r1), paste0(zb_r1, "-", za_r1), NA_character_ # nolint
      ),
      data.table::fifelse(!is.na(rx_r1), rx_r1, NA_character_) # nolint
    )]
    merged <- merged[!is.na(umi)] # nolint
  }

  if (nrow(merged) > 0) {
    merged[, pos := pmin(start0_r1, start0_r2)] # nolint
    merged[, mpos := pmax(start0_r1, start0_r2)] # nolint
    merged[, epos := pmax(end0_r1, end0_r2)] # nolint
    rows <- merged[, .( # nolint
      qname, chrom = chrom_r1, pos, mpos, epos, strand = strand_out, umi # nolint
    )]
  } else {
    rows <- data.table::data.table(
      qname = character(0), chrom = character(0),
      pos = integer(0), mpos = integer(0), epos = integer(0),
      strand = character(0), umi = character(0)
    )
  }

  list(rows = rows, pending = pending)
}

# Stream `bam_file` and write extracted per-read-pair rows
# (qname, chrom, pos, mpos, epos, strand, umi) to `extract_out`.
# Returns the number of pending (never-matched) reads at EOF, for
# diagnostics.
.extract_read_pairs <- function(bam_file, extract_out, chunk_size) {
  bf <- Rsamtools::BamFile(bam_file, yieldSize = chunk_size)
  Rsamtools::open.BamFile(bf)
  on.exit(Rsamtools::close.BamFile(bf), add = TRUE)

  param <- Rsamtools::ScanBamParam(
    flag = .bam_read_info_flag_filter(),
    what = c("qname", "flag", "mrnm"),
    tag = c("RX", "ZA", "ZB")
  )

  pending <- data.table::data.table(
    qname = character(0), chrom = character(0),
    start0 = integer(0), end0 = integer(0),
    is_reverse = logical(0), is_read1 = logical(0),
    rx = character(0), za = character(0), zb = character(0)
  )

  n_written <- 0L
  first_write <- TRUE

  # progress logging: emit a line for every `progress_every` reads
  # (alignment records) streamed from the BAM file
  progress_every <- 1000000L
  n_read <- 0L
  next_progress <- progress_every

  message("Extracting read info from BAM...")

  repeat {
    aln <- GenomicAlignments::readGAlignments(bf, param = param)
    if (length(aln) == 0) break

    n_read <- n_read + length(aln)

    chunk <- .bam_chunk_to_dt(aln)
    combined <- data.table::rbindlist(list(pending, chunk), use.names = TRUE)
    res <- .pair_mates(combined)
    pending <- res$pending

    if (nrow(res$rows) > 0) {
      data.table::fwrite(
        res$rows, extract_out,
        sep = "\t", col.names = FALSE, append = !first_write
      )
      first_write <- FALSE
      n_written <- n_written + nrow(res$rows)
    }

    while (n_read >= next_progress) {
      message(sprintf("  ...%d reads extracted", next_progress))
      next_progress <- next_progress + progress_every
    }
  }

  message(sprintf("Finished extracting reads: %d reads processed, %d read pairs written", # nolint
                  n_read, n_written))

  if (first_write) {
    # nothing was ever written; create an empty file so downstream
    # steps have something to operate on
    file.create(extract_out)
  }

  list(n_written = n_written, n_pending = nrow(pending))
}

# Run the external sort/aggregate step:
#  sort -k 2,7 -k 1,1 -u --parallel=<cores> <extract> -T <tmp_dir> -S <sort_mem>
#     | cut -f 2-7 | uniq --count > <sorted_out>
# Kept as an external pipeline (rather than an in-memory R aggregation)
# so global grouping stays disk-backed and scalable for 100GB+ inputs.
.sort_and_aggregate <- function(extract_out, sorted_out, cores, sort_mem,
                                tmp_dir) {
  if (!nzchar(Sys.which("sort")) || !nzchar(Sys.which("cut")) ||
        !nzchar(Sys.which("uniq"))) {
    stop(
      "generate_read_info_from_bam requires the 'sort', 'cut' and 'uniq' ",
      "system utilities (GNU coreutils) to be available on PATH."
    )
  }

  message("Sorting and aggregating read bundles...")

  if (file.size(extract_out) == 0) {
    file.create(sorted_out)
    return(invisible(TRUE))
  }

  cmd <- sprintf(
    "sort -k 2,7 -k 1,1 -u --parallel=%d %s -T %s -S %s | cut -f 2-7 | uniq --count > %s", # nolint
    as.integer(cores),
    shQuote(extract_out),
    shQuote(tmp_dir),
    shQuote(sort_mem),
    shQuote(sorted_out)
  )
  status <- system(cmd)
  if (status != 0) {
    stop("generate_read_info_from_bam: sort/cut/uniq pipeline failed ",
         "(exit status ", status, ").")
  }
  invisible(TRUE)
}

# Pivot the sorted/counted (chrom,pos,mpos,epos,strand,umi) rows wide
# into the final rinfo schema (chrom,pos,mpos,epos,umi,x,y), where x/y
# are plus-/minus-strand bundle counts. The first field of each line
# is "<count> <chrom>" separated by a single space (from `uniq -c`),
# with the remaining fields tab-separated -- so it must be re-split.
.pivot_sorted_counts <- function(sorted_out) {
  message("Pivoting strand counts into read-bundle table...")

  empty <- data.table::data.table(
    chrom = character(0), pos = integer(0), mpos = integer(0),
    epos = integer(0), umi = character(0), x = integer(0), y = integer(0)
  )

  if (file.size(sorted_out) == 0) {
    return(empty)
  }

  dt <- data.table::fread(sorted_out, sep = "\t", header = FALSE)
  if (nrow(dt) == 0 || ncol(dt) != 6) {
    return(empty)
  }
  data.table::setnames(
    dt, c("V1", "V2", "V3", "V4", "V5", "V6"),
    c("count_chrom", "pos", "mpos", "epos", "strand", "umi")
  )

  split_cols <- data.table::tstrsplit(trimws(dt$count_chrom), "\\s+")
  dt[, count := as.integer(split_cols[[1]])]
  dt[, chrom := split_cols[[2]]] # nolint
  dt[, count_chrom := NULL] # nolint

  plus  <- dt[strand == "+", .(chrom, pos, mpos, epos, umi, x = count)] # nolint
  minus <- dt[strand == "-", .(chrom, pos, mpos, epos, umi, y = count)] # nolint

  merged <- merge(
    plus, minus, by = c("chrom", "pos", "mpos", "epos", "umi"), all = TRUE
  )
  merged[is.na(x), x := 0L] # nolint
  merged[is.na(y), y := 0L] # nolint
  data.table::setcolorder(
    merged, c("chrom", "pos", "mpos", "epos", "umi", "x", "y")
  )
  merged
}

# Validate that the final table conforms to the expected rinfo schema
# before it is written out: required columns present, no NA in key
# columns, and non-negative integer read-bundle counts.
.validate_read_info_schema <- function(dt) {
  required_cols <- c("chrom", "pos", "mpos", "epos", "umi", "x", "y")
  missing_cols <- setdiff(required_cols, names(dt))
  if (length(missing_cols) > 0) {
    stop("generate_read_info_from_bam: internal error - generated table is ",
         "missing required column(s): ", paste(missing_cols, collapse = ", "))
  }

  if (nrow(dt) == 0) {
    return(invisible(TRUE))
  }

  key_cols <- c("chrom", "pos", "mpos", "epos", "umi") # nolint
  na_counts <- vapply(
    dt[, ..key_cols], function(col) sum(is.na(col)), integer(1) # nolint  
  )
  if (any(na_counts > 0)) {
    stop("generate_read_info_from_bam: internal error - NA values found in ",
         "required column(s): ",
         paste(names(na_counts)[na_counts > 0], collapse = ", "))
  }

  if (any(dt$x < 0L) || any(dt$y < 0L)) {
    stop("generate_read_info_from_bam: internal error - negative read-bundle ",
         "counts detected in x/y.")
  }

  if (any(dt$pos < 0L) || any(dt$mpos < 0L) || any(dt$epos < 0L)) {
    stop("generate_read_info_from_bam: internal error - negative genomic ",
         "coordinates detected in pos/mpos/epos.")
  }

  invisible(TRUE)
}

#' Generate a summarised read-info ("rinfo") table directly from a BAM file.
#'
#' Streams a coordinate-sorted BAM file in chunks, pairs mates, extracts
#' strand-adjusted UMIs, and aggregates into the same read-bundle schema
#' (chrom, pos, mpos, epos, umi, x, y) produced by the legacy
#' scripts/generate_read_info.sh pipeline. Designed to bound memory use
#' for very large (100GB+) BAM files by streaming extraction and using an
#' external disk-backed sort for global aggregation.
#'
#' @param bam_file Path to an input BAM file (coordinate-sorted).
#' @param output Path to write the resulting gzipped rinfo table to. If
#'   `NULL`, a temporary file is created and its path returned.
#' @param chunk_size Number of BAM records to read into memory per chunk
#'   (Rsamtools `yieldSize`). Larger values increase throughput at the
#'   cost of memory; default 2e6.
#' @param cores Number of threads passed to the external `sort` step
#'   (`sort --parallel`). Defaults to 1.
#' @param sort_mem Memory buffer passed to `sort -S` (e.g. "4G").
#' @param tmp_dir Directory for `sort`'s temporary spill files and this
#'   function's own intermediate files. Defaults to `tempdir()`.
#'
#' @return Invisibly, a list with `output` (path to the written rinfo
#'   file) and `n_bundles` (number of read-bundle rows written).
#' @export
generate_read_info_from_bam <- function(
  bam_file,
  output = NULL,
  chunk_size = 2e6,
  cores = 1,
  sort_mem = "4G",
  tmp_dir = tempdir()
) {
  bam_file <- normalizePath(bam_file, mustWork = TRUE)

  if (is.null(output) || !nzchar(output)) {
    output <- tempfile(pattern = "rinfo_from_bam_", tmpdir = tmp_dir,
                       fileext = ".txt.gz")
  }
  odir <- dirname(output)
  if (!dir.exists(odir))
    dir.create(odir, recursive = TRUE, showWarnings = FALSE)

  if (!dir.exists(tmp_dir)) {
    dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
  }

  extract_out <- tempfile(pattern = "rinfo_extract_", tmpdir = tmp_dir,
                          fileext = ".txt")
  sorted_out  <- tempfile(pattern = "rinfo_sorted_", tmpdir = tmp_dir,
                          fileext = ".txt")
  on.exit(unlink(c(extract_out, sorted_out)), add = TRUE)

  extract_res <- .extract_read_pairs(bam_file, extract_out, chunk_size)

  if (extract_res$n_pending > 0) {
    message(sprintf(
      paste0(
        "generate_read_info_from_bam: %d read(s) never found a mate ",
        "and were dropped (mate missing, filtered, or absent from the ",
        "BAM stream)."
      ),
      extract_res$n_pending
    ))
  }

  .sort_and_aggregate(extract_out, sorted_out, cores = cores,
                      sort_mem = sort_mem, tmp_dir = tmp_dir)

  final_tbl <- .pivot_sorted_counts(sorted_out)

  .validate_read_info_schema(final_tbl)

  message(sprintf("Writing %d read bundle(s) to %s...",
                  nrow(final_tbl), output))
  data.table::fwrite(final_tbl, output, sep = "\t")
  message("Done writing rinfo table from BAM.")

  invisible(list(output = output, n_bundles = nrow(final_tbl)))
}
