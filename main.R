#!/usr/bin/env Rscript

# ------------------------------------------------------------------
# main.R
#
# LOCAL DEVELOPMENT ENTRYPOINT (renv / project-root only).
#
# This script is for running the tool directly from the project root
# without installing the package (e.g. via renv). It sources R/*.R
# using relative paths that only resolve from the project root.
#
# When the package is installed (R CMD INSTALL, Docker, conda), use
# the installed CLI instead:
#   calc-duplex-metrics --input <file> --output <file>
#
# The installed entry point is at:
#   inst/exec/calc-duplex-metrics  (within the package source)
# ------------------------------------------------------------------


# Check required packages are available (install via renv::restore())
required_packages <- c("argparse")
missing_packages <-
  required_packages[!required_packages %in% installed.packages()[, "Package"]]

if (length(missing_packages) > 0) {
  stop(
    "Missing required package(s): ", paste(missing_packages, collapse = ", "), "\n",
    "Please run renv::restore() before running this script.",
    call. = FALSE
  )
}

# Load required libraries
suppressPackageStartupMessages({
  library(argparse)
})

# Source all package files explicitly for local renv use
# (In the installed package these are loaded automatically by R)
source("R/calculate_nanoseq_functions.R")
source("R/calculate.R")
source("R/cli.R")
main()
