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
# Note: step 3 requires output/config/auto_binning.txt to exist before
# Snakemake can parse Snakefile-bin. On a first dry run (-n), this file
# won't exist yet and step 3 will be skipped with a message. Run without
# -n to execute stages 1-2 first; step 3 will then run on the same call
# or any subsequent call once the file exists.
#
# Usage:
#   ./run_magmaker.sh --profile resources/profiles/demon
#   ./run_magmaker.sh --cores 8 --use-conda
#   ./run_magmaker.sh --profile resources/profiles/demon -n   # dry run
#
# To stop after the MAG summary table so you can review/edit
# output/mag_qc/mag_summary.tsv before renaming, run stages manually:
#   snakemake [options]
#   snakemake [options] generate_binning_config
#   snakemake --snakefile Snakefile-bin [options] make_mag_summary \
#     --config binning=output/config/auto_binning.txt

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
if [[ ! -f output/config/auto_binning.txt ]]; then
    echo ""
    echo " SKIPPED: output/config/auto_binning.txt does not exist yet."
    echo " This is expected on a dry run (-n) — stages 1 and 2 must"
    echo " run first to generate it. Re-run without -n to execute the"
    echo " full pipeline, or run the binning stage manually:"
    echo ""
    echo "   snakemake --snakefile Snakefile-bin [options] rename_mags \\"
    echo "     --config binning=output/config/auto_binning.txt"
    echo ""
else
    snakemake --snakefile Snakefile-bin "$@" rename_mags \
        --config binning=output/config/auto_binning.txt

    echo ""
    echo "========================================"
    echo " MAGmaker  Complete"
    echo " Output:  output/mag_qc/renamed_mags/"
    echo " Summary: output/mag_qc/mag_summary.tsv"
    echo "========================================"
fi
