# Copilot Instructions

## Project Overview

An R CLI tool that takes summarised duplex sequencing read data (`.rinfo` files) and outputs a set of duplex QC metrics in long-format CSV.

## Running Tests

```bash
# Run all tests (from project root)
Rscript tests/testthat.R

# Run a single test file
Rscript -e "devtools::load_all(); testthat::test_file('tests/testthat/test-metrics.R')"
```

Tests live in `tests/testthat/test-metrics.R`. The test fixture is at `tests/testthat/testdata/NanoMB1Rep1_test_10k.txt` and referenced via `testthat::test_path()`. One test (`calculate_gc` with reference) downloads a FASTA from NCBI and requires internet access.

## Architecture

The call chain is strictly linear:

```
main.R  →  R/cli.R (parse_arguments)  →  R/calculate.R (process_data)  →  R/calculate_nanoseq_functions.R
```

- **`main.R`** — local renv-only entrypoint; sources all three R files explicitly for use from the project root. Not included in the installed package. The installed CLI is `inst/exec/calc-duplex-metrics`.
- **`R/cli.R`** — argument parsing and validation only; delegates to `process_data()`
- **`R/calculate.R`** — orchestration: resolves metric selection, gates GC logic, handles single/multi-file dispatch, writes output CSV
- **`R/calculate_nanoseq_functions.R`** — pure calculation functions; no argument parsing; no `source()` calls

## Key Conventions

### Input format
`.rinfo` files are tab-separated, read with `data.table::fread()`. Columns include `x` (forward strand reads), `y` (reverse strand reads), `chrom`, and `pos`. Column indices 5–6 are renamed to `plus`/`minus` inside `calculate_gc()`.

### Metric selection
Metrics are resolved once in `process_data()` via `resolve_metric_selection()` before any file is processed. There are two kinds:
- **Individual**: `frac_singletons`, `efficiency`, `drop_out_rate`
- **Groups**: `gc` (→ `gc_single`, `gc_both`, `gc_deviation`), `family` (→ 8 stats)

Specifying a sub-metric (e.g. `gc_single`) causes the whole group to be computed.

### GC gating
- `--metrics all` (default): GC is silently skipped if `--ref_fasta` is absent.
- Explicit GC request without `--ref_fasta`: exits with error.
- The `genomeFile` (Rsamtools `FaFile`) and `genome_max` (named integer vector of chromosome lengths) are prepared once in `process_data()` and passed down — never constructed inside metric functions.

### Error handling
Top-level functions (`process_data`, `calc_duplex_metrics_one_file_df`) return `list(success = FALSE, error = "...")` on failure rather than throwing. Lower-level calculation functions throw directly with `stop()`.

### Output
Always long-format CSV with columns `sample,metric,value` (no row names, no quoting). The output directory is created automatically if missing.

### Parallelism
Multi-file jobs use `parallel::mclapply` on non-Windows; falls back to `lapply` on Windows or when `cores = 1`.

### Sample names
Derived from input filenames by stripping `.txt` / `.txt.gz` suffixes when `--sample` is not provided.

### Environment / installation

Three installation modes:
- **Installed package** (Docker / conda): `R CMD INSTALL .` installs the package; `inst/exec/calc-duplex-metrics` is symlinked to `$PATH`. CLI: `calc-duplex-metrics --input ...`
- **renv (local dev)**: run `main.R` from the project root via `Rscript main.R`. `main.R` sources all three R files explicitly — it is the only file in the repo that uses `source()`.
- **conda setup**: `conda env create -f env/environment.yaml` then `bash env/install.sh` (installs package + symlinks CLI).

When running under conda **without** installing the package, use `Rscript --no-init-file main.R` to skip the `.Rprofile` that activates renv.
