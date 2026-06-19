#!/usr/bin/env bash
# run_magmaker.sh
#
# Runs the complete MAGmaker pipeline in sequence:
#   1. Main pipeline  — QC → assembly → profiling → prototype selection
#   2. Binning config — auto-generates binning.txt from selected prototypes
#   3. Binning pipeline — mapping → binning → DAS_Tool → MAG QC → rename
#
# All arguments are passed through to each snakemake invocation, so
# --profile, --cores, --use-conda, -n (dry run), etc. all work as expected.
#
# Usage:
#   ./run_magmaker.sh --profile resources/profiles/demon
#   ./run_magmaker.sh --cores 8 --use-conda
#   ./run_magmaker.sh --profile resources/profiles/demon -n   # dry run all steps
#
# To stop after the MAG summary table (before renaming, so you can review/edit
# output/mag_qc/mag_summary.tsv first), run the three stages manually:
#   snakemake [options]
#   snakemake [options] generate_binning_config
#   snakemake --snakefile Snakefile-bin [options] make_mag_summary \
#     --config binning=output/config/auto_binning.txt
# Then edit the table and run rename_mags when ready.

set -euo pipefail

if [[ ! -f Snakefile || ! -f Snakefile-bin ]]; then
    echo "ERROR: run_magmaker.sh must be run from the MAGmaker repository root." >&2
    exit 1
fi

echo "========================================"
echo " MAGmaker  Step 1/3 — Main pipeline"
echo "========================================"
snakemake "$@"

echo ""
echo "========================================"
echo " MAGmaker  Step 2/3 — Binning config"
echo "========================================"
snakemake "$@" generate_binning_config

echo ""
echo "========================================"
echo " MAGmaker  Step 3/3 — Binning pipeline"
echo "========================================"
snakemake --snakefile Snakefile-bin "$@" rename_mags \
    --config binning=output/config/auto_binning.txt

echo ""
echo "========================================"
echo " MAGmaker  Complete"
echo " Output: output/mag_qc/renamed_mags/"
echo " Summary: output/mag_qc/mag_summary.tsv"
echo "========================================"
