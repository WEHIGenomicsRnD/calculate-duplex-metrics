# Calculate Duplex Metrics

An R CLI tool to take in summarised read information and output duplex
metrics.

### Available metrics

Individual metrics (selectable individually)
- efficiency
- drop_out_rate

Grouped metrics

GC metrics 
- gc_single
- gc_both
- gc_deviation

Family stats
- total_families
- family_mean
- family_median
- family_max
- families_gt1
- single_families
- paired_families
- unpaired_families
- paired_and_gt1
- frac_singletons

Additional optional metrics
- on_target_rate_raw (computed when `--target_bed` is supplied)
- on_target_rate_duplex (computed when `--target_bed` is supplied)
- on_target_coverage_raw (computed when `--target_bed` is supplied)
- on_target_coverage_duplex (computed when `--target_bed` is supplied)
- on_target_duplex_ratio (computed when `--target_bed` is supplied; ratio of raw to duplex coverage)

### Supported input formats

- `rinfo` (default): supports all currently available metrics.
- `fgbio`: supports `frac_singletons`, `efficiency`,
  `drop_out_rate`, and family stats from fgbio
  `CollectDuplexSeqMetrics` `*.duplex_family_sizes.txt` output.
  GC metrics are not available for this format because the table does
  not contain genomic coordinates. `on_target_rate_*` is also unavailable
  because fgbio duplex family-size input does not include read positions.
- `bam`: computes metrics directly from a single coordinate-sorted BAM
  file by first extracting read info internally. BAM mode supports the
  same downstream metrics as `rinfo`, including GC and on-target metrics
  when the relevant inputs are supplied.

## Implementation overview

- The CLI entrypoint is `main.R`
- Argument parsing and validation are handled in `cli.R`
- Metric execution logic is in `calculate.R`
- Core metric implementations are defined in `R/metric_functions.R`

Metric selection is resolved **before computation**.  
Only the requested individual metrics and/or metric groups are evaluated.

### GC metric behaviour
- GC metrics are computed only when a reference genome object (.fasta) is provided.
- `--metrics` defaults to `all`.
  - If `--ref_fasta` is not provided, GC metrics are skipped and a message is printed to the console.
  - If `--ref_fasta` is provided, GC metrics are computed (may return NA if insufficient data).
- If GC metrics are explicitly requested (e.g. `--metrics gc`) but no reference FASTA is supplied, the program exits with an error.
- If `--input_format fgbio` is used, GC metrics are
  always skipped in default mode and rejected when requested explicitly.

### On-target rate calculation

The on-target rate calculations are approximate because:

- If the `read_info` format includes an `epos` column (actual end position),
  it will be used for accurate end coordinates. Otherwise, the end coordinate
  is estimated using the read length.
- The on-target duplex rate is calculated from the bundle counts, but those
  reads may be further filtered downstream.
- For BAM input, `epos` is derived from the aligned read end and used for
  on-target calculations.

## Installation and Usage

### Option A: Using Docker (Recommended)

This method packages the script and all its dependencies into a self-contained environment. It is the most reliable way to run the analysis, as it guarantees that the exact same software versions are used every time.

#### Requirements
- **Docker:** You must have Docker installed and the Docker daemon running. You can download it from the [Docker website](https://www.docker.com/products/docker-desktop/).

#### Installation Steps

1.  **Clone the repository:**
    ```bash
    git clone https://github.com/WEHIGenomicsRnD/calculate-duplex-metrics.git
    ```

2.  **Navigate to the project directory:**
    ```bash
    cd calculate-duplex-metrics
    ```

3.  **Build the Docker image:** This command reads the `Dockerfile` and builds a container image named `calculate-duplex-metrics`. This may take several minutes the first time you run it.
    ```bash
    docker build -t calculate-duplex-metrics .
    ```

4.  **(Optional) Verify the image:** You can check that the image was built successfully by listing your local Docker images.
    ```bash
    docker images
    ```
    You should see `calculate-duplex-metrics` in the list.

#### Default Usage Example

To run the tool, you use the `docker run` command. The `-v` flags are essential for allowing the Docker container to access files on your local machine.

```bash
docker run --rm \
  -v "$(pwd)/data:/app/data" \
  -v "$(pwd)/out:/app/out" \
  calculate-duplex-metrics \
  calc-duplex-metrics \
  --input data/test.rinfo \
  --output metrics.csv
```

  To process a BAM file directly:

  ```bash
  calc-duplex-metrics \
    --input data/test.bam \
    --input_format bam \
    --output bam_metrics.csv
  ```

-   `-v "$(pwd)/data:/app/data"`: This "mounts" your local `data` directory into the `/app/data` directory inside the container, so the script can find the input file.
-   `-v "$(pwd)/out:/app/out"`: This mounts your local `out` directory into the `/app/out` directory inside the container, so the script can write the output file back to your machine.

### Option B: Local Installation with `renv`

This method uses the `renv` package to recreate the exact development environment, using the specific package versions defined in the `renv.lock` file. This is the recommended approach for local development and contributions.

> **Note:** `renv` mode uses `main.R` as the entrypoint (run from the project root). `main.R` is a local-only launcher that `source()`s the package files directly. The installed CLI (`calc-duplex-metrics`) is not available in this mode.

#### Requirements
- **R:** R version **4.4.1** 

#### Installation Steps

1.  **Clone the repository and navigate into it.**

2.  **Open an R console** in the project's root directory.

3.  **Restore the environment:** This command will install `renv` if needed, then install all the packages listed in `renv.lock` with their exact versions.
    ```r
    if (!require("renv")) install.packages("renv", repos = "https://cloud.r-project.org")
    renv::restore()
    ```

#### Default Usage Example

After restoring the environment, you can run the script directly from your terminal within the cloned repository directory.

```bash
Rscript main.R \
  --input data/test.rinfo \
  --output metrics.csv
```

To process a BAM file directly:

```bash
Rscript main.R \
  --input data/test.bam \
  --input_format bam \
  --output bam_metrics.csv
```

### Option C: Local Installation with `devtools`

This method installs the required R packages directly onto your system. It is more flexible if you have a different version of R, but it is less reproducible as it will use the latest available package versions.

#### Requirements
- **R:** Any modern version of R.
- **`devtools` R package:** This is used to install packages from GitHub.

#### Installation Steps

1.  **Open an R console.**

2.  **Install `devtools` and Bioconductor dependencies:**
    ```r
    # Install devtools from CRAN
    install.packages("devtools")

    # Install BiocManager and required Bioconductor packages
    if (!require("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
    BiocManager::install(c("Rsamtools", "GenomicAlignments", "GenomicRanges", "IRanges", "Biostrings"))
    ```

3.  **Install the package from GitHub:**
    ```r
    devtools::install_github('WEHIGenomicsRnD/calculate-duplex-metrics')
    ```

#### Default Usage Example

After installing the dependencies, you can run the script directly from your terminal within the cloned repository directory.

```bash
Rscript main.R \
  --input data/test.rinfo \
  --output metrics.csv
```

To process a BAM file directly:

```bash
Rscript main.R \
  --input data/test.bam \
  --input_format bam \
  --output bam_metrics.csv
```

### Option D: Using Conda

This method uses Conda to manage the R environment and dependencies, then installs the package into the environment so the `calc-duplex-metrics` CLI is available on `PATH`.

#### Requirements
- **Conda or Miniconda:** You must have Conda installed. You can download it from the [Conda website](https://docs.conda.io/en/latest/miniconda.html).

#### Installation Steps

1.  **Clone the repository and navigate into it.**

2.  **Create the conda environment:**
    ```bash
    conda env create -f env/environment.yaml
    ```

3.  **Activate the environment:**
    ```bash
    conda activate calculate-duplex-metrics
    ```

4.  **Install the package and register the CLI:**
    ```bash
    bash env/install.sh
    ```
    This runs `R CMD INSTALL .` and symlinks `calc-duplex-metrics` into the conda environment's `bin/`.

#### Default Usage Example

Validate that the CLI can be called:

```bash
# Installed package (Docker / conda):
calc-duplex-metrics --help

# Local renv development:
Rscript main.R --help
```

Basic usage is as follows:

```bash
calc-duplex-metrics \
  --input data/test.rinfo \
  --output metrics.csv
```

For fgbio duplex family size input, set `--input_format` explicitly:

```bash
calc-duplex-metrics \
  --input data/NanoMB1Rep1.duplex_seq_metrics.duplex_family_sizes.txt \
  --input_format fgbio \
  --output fgbio_metrics.csv \
  --metrics frac_singletons,efficiency,drop_out_rate,family
```

To deactivate the conda environment when finished:
```bash
conda deactivate
```

## Additional Usage Examples

#### Example: default mode with GC enabled (requires reference genome)

Note: The reference genome FASTA is user-provided and not included
in this repository. Any compatible reference genome may be used.

``` bash
Rscript main.R \
  --input data/NanoMB1Rep1_HJK2GDSX3_CGGCTAAT-CTCGTTCT_L001.txt \
  --output default_with_gc.csv \
  --ref_fasta ref/Escherichia_coli_ATCC_10798.fasta
```

#### Example: select individual metrics only

``` bash
Rscript main.R \
  --input data/NanoMB1Rep1_HJK2GDSX3_CGGCTAAT-CTCGTTCT_L001.txt \
  --output test_selected_metrics.csv \
  --metrics efficiency,drop_out_rate
```
Note: when listing multiple metrics, either omit spaces (efficiency,drop_out_rate) or quote the argument ("efficiency, drop_out_rate").

#### Example: select metric groups

``` bash
Rscript main.R \
  --input data/NanoMB1Rep1_HJK2GDSX3_CGGCTAAT-CTCGTTCT_L001.txt \
  --output test_family_metrics.csv \
  --metrics family
```

#### Example: mixed selection (individual + group)

``` bash
Rscript main.R \
  --input data/NanoMB1Rep1_HJK2GDSX3_CGGCTAAT-CTCGTTCT_L001.txt \
  --output test_mixed_metrics.csv \
  --metrics efficiency,family
```

#### Example: multiple input files

``` bash
Rscript main.R \
  --input data/a.txt data/b.txt \
  --output all_samples_metrics.csv \
  --metrics family \
  --cores 2
```

#### Example: BAM input mode

``` bash
Rscript main.R \
  --input data/test.bam \
  --input_format bam \
  --output bam_metrics.csv \
  --cores 4 \
  --sort_mem 8G \
  --bam_chunk_size 2000000
```


### CLI flags

``` bash
Required:
  -i, --input        One or more metric input files OR a directory
                     containing metric input files (.txt or .txt.gz)
                     Note: when --input is a directory, the tool selects matching files using --pattern (default: \.txt(\.gz)?$);
                           when --input is a list of files, --pattern is ignored
  -o, --output       Output CSV path (long format)

Optional:
  -s, --sample       Optional sample name(s). For multiple input files, provide
                     comma-separated names matching the number of files.
                     Note: if --input is a directory, --sample is not allowed.

      --pattern      Regex pattern used to select files when --input is a directory
                     (default: \.txt(\.gz)?$)

      --rlen         Read length (default: 151)
      --skips        Trimmed / ignored bases per read (NanoSeq = 5, xGen = 8)

      --ref_fasta    Reference genome FASTA (required for GC metrics)
      --target_bed   Optional BED file of target regions.
                     If supplied, on_target_rate_raw and
                     on_target_rate_duplex are computed.
                     Note that these calculations are approximate.
                     Not supported for --input_format fgbio.
      --min_reads    Minimum reads used for on_target_rate_duplex.
                     Used with --target_bed.
                     Format: "total duplex1 duplex2" (default: 4 2 2).

      --input_format Input format (default: rinfo)
                     - rinfo
                     - fgbio
                     - bam
                     BAM mode converts a single coordinate-sorted BAM
                     into read-info internally before metric calculation.

      --metrics      Comma-separated list of metrics and/or metric groups
                     - Individual: efficiency, drop_out_rate
                     - Groups: gc, family
                     (default: all)

                     Note: for fgbio input, supported
                     metrics are frac_singletons, efficiency,
                     drop_out_rate, and family.

      --cores        Number of CPU cores for parallel processing (default: 1)

      --sort_mem     Memory buffer passed to the external sort step used in
                     BAM mode (default: 4G)

      --bam_chunk_size  Number of BAM records read per chunk in BAM mode
                        (default: 2000000)


```
Note: when listing multiple metrics, either omit spaces (efficiency,drop_out_rate) or quote the argument ("efficiency, drop_out_rate").


## Outputs

Output is written in long format:
```
sample,metric,value

``` 

#### Example:

```
sample,metric,value
NanoMB1Rep1_HJK2GDSX3_CGGCTAAT-CTCGTTCT_L001,efficiency,0.0490258329591602
NanoMB1Rep1_HJK2GDSX3_CGGCTAAT-CTCGTTCT_L001,drop_out_rate,0.320805646128878
NanoMB1Rep1_HJK2GDSX3_CGGCTAAT-CTCGTTCT_L001,total_families,23825702
NanoMB1Rep1_HJK2GDSX3_CGGCTAAT-CTCGTTCT_L001,family_mean,6.748161712309
NanoMB1Rep1_HJK2GDSX3_CGGCTAAT-CTCGTTCT_L001,family_median,5
NanoMB1Rep1_HJK2GDSX3_CGGCTAAT-CTCGTTCT_L001,family_max,50
NanoMB1Rep1_HJK2GDSX3_CGGCTAAT-CTCGTTCT_L001,families_gt1,16771629
NanoMB1Rep1_HJK2GDSX3_CGGCTAAT-CTCGTTCT_L001,single_families,6731955
NanoMB1Rep1_HJK2GDSX3_CGGCTAAT-CTCGTTCT_L001,paired_families,9994045
NanoMB1Rep1_HJK2GDSX3_CGGCTAAT-CTCGTTCT_L001,paired_and_gt1,8152302
NanoMB1Rep1_HJK2GDSX3_CGGCTAAT-CTCGTTCT_L001,frac_singletons,0.0418706803079419

```

#### Sanity check the CLI

## Testing
Test that functions return valid numeric values, correct handling of edge cases (NA, zero reads, invalid inputs) and presence of expected metrics names.

#### Requirements
Packages: `devtools`, `testthat`

#### To run all tests
From the project root run:
```bash
Rscript tests/testthat.R
```

#### To run a single test file
```bash
Rscript -e "devtools::load_all(); testthat::test_file('tests/testthat/test-metrics.R')"
```
