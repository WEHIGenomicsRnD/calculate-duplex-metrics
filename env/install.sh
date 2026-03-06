#!/usr/bin/env bash
# Install the CalcDuplexMetrics package into the active conda environment
# and place the CLI on PATH.
#
# Usage (from the project root):
#   conda env create -f env/environment.yaml
#   conda activate calculate-duplex-metrics
#   bash env/install.sh
#
# After this, `calc-duplex-metrics --help` should work from within the env.

set -euo pipefail

# Bypass the project .Rprofile (which activates renv) for all R calls in this script
export R_PROFILE_USER=/dev/null

if [[ -z "${CONDA_PREFIX:-}" ]]; then
  echo "ERROR: No active conda environment detected. Activate the environment first:"
  echo "  conda activate calculate-duplex-metrics"
  exit 1
fi

echo "Installing CalcDuplexMetrics into: $CONDA_PREFIX"

R CMD INSTALL .

SCRIPT_PATH="$(Rscript --vanilla -e 'cat(system.file("exec", "calc-duplex-metrics", package="CalcDuplexMetrics"))')"

if [[ -z "$SCRIPT_PATH" ]]; then
  echo "ERROR: Could not locate installed exec script. Check that R CMD INSTALL succeeded."
  exit 1
fi

ln -sf "$SCRIPT_PATH" "$CONDA_PREFIX/bin/calc-duplex-metrics"

echo "Done. calc-duplex-metrics is now available at: $CONDA_PREFIX/bin/calc-duplex-metrics"
echo "Test with: calc-duplex-metrics --help"
