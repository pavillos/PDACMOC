#!/usr/bin/env bash
# Create the "pdacmoc" conda environment and install PDACMOC.
# Usage: bash PDACMOC/install/install.sh PDACMOC_<version>.tar.gz [env_name]
# Needs conda (Miniforge or Miniconda). Tested on Linux.
set -euo pipefail

TARBALL="$1"
ENV_NAME="${2:-pdacmoc}"
HERE="$(cd "$(dirname "$0")" && pwd)"
ORG_URL="https://bioconductor.org/packages/3.18/data/annotation/src/contrib/org.Hs.eg.db_3.18.0.tar.gz"

conda env create -n "$ENV_NAME" -f "$HERE/environment.yml"
PREFIX="$(conda env list | awk -v n="$ENV_NAME" '$1 == n {print $NF}')"
R_CMD="$PREFIX/bin/R"

WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

# Pinned annotation package
curl -fsSL -o "$WORK/org.Hs.eg.db_3.18.0.tar.gz" "$ORG_URL"
"$R_CMD" CMD INSTALL "$WORK/org.Hs.eg.db_3.18.0.tar.gz"

# ADVOCATE (shipped inside the PDACMOC tarball) and PDACMOC
tar -xzf "$TARBALL" -C "$WORK"
"$R_CMD" CMD INSTALL "$WORK"/PDACMOC/inst/packages/ADVOCATE_*.tar.gz
"$R_CMD" CMD INSTALL "$WORK/PDACMOC"

echo "Done. In R, point reticulate to: $PREFIX/bin/python"
