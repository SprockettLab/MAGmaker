#!/usr/bin/env bash
# run_magmaker.sh
#
# Runs the complete MAGmaker pipeline in sequence:
#   1. Main pipeline  — QC → assembly → profiling → prototype selection
#   2. Binning config — auto-generates binning.txt from selected prototypes
#   3. Binning pipeline — mapping → binning → DAS_Tool → MAG QC → rename
#   4. Viral track — geNomad + CheckV viral/plasmid discovery
#                    (on by default; turn off with --skip-virus)
#
# All arguments except --skip-virus are passed through to each snakemake
# invocation, so --profile, --cores, --use-conda, -n (dry run), etc. all work.
#
# --skip-virus is consumed by this script and NOT forwarded to snakemake; it
# turns off stage 4. The viral track needs the geNomad + CheckV databases set in
# your config (params.genomad.db_path / params.checkv.db_path), so use
# --skip-virus if those are not configured.
#
# Note: steps 3 and 4 require output/config/auto_binning.txt to exist before
# Snakemake can parse Snakefile-bin. On a first dry run (-n), this file won't
# exist yet and they are skipped with a message. Run without -n to execute
# stages 1-2 first; steps 3-4 then run on the same call or any subsequent call.
#
# Usage:
#   ./run_magmaker.sh --profile resources/profiles/demon
#   ./run_magmaker.sh --profile resources/profiles/demon --skip-virus
#   ./run_magmaker.sh --cores 8 --use-conda
#   ./run_magmaker.sh --profile resources/profiles/demon -n   # dry run
#   ./run_magmaker.sh --profile resources/profiles/demon \
#       --configfile resources/config/<Project>_config.yaml   # project config
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

# Pull run_magmaker.sh's own flags out of the snakemake pass-through args.
# --skip-virus is ours; everything else is forwarded verbatim.
RUN_VIRUS=1
ARGS=()
for arg in "$@"; do
    if [[ "$arg" == "--skip-virus" ]]; then
        RUN_VIRUS=0
    else
        ARGS+=("$arg")
    fi
done
# Reset positional params to the filtered args (guarded for an empty array
# under `set -u`).
set -- "${ARGS[@]+"${ARGS[@]}"}"

TOTAL=4
[[ "$RUN_VIRUS" -eq 0 ]] && TOTAL=3

echo "========================================"
echo " MAGmaker  Step 1/$TOTAL — Main pipeline"
echo "========================================"
snakemake "$@"

echo ""
echo "========================================"
echo " MAGmaker  Step 2/$TOTAL — Binning config"
echo "========================================"
# Target goes before "$@" so a trailing variable-length option in the
# user's args (e.g. --configfile FILE, which accepts one or more values)
# cannot swallow the target name.
snakemake generate_binning_config "$@"

echo ""
echo "========================================"
echo " MAGmaker  Step 3/$TOTAL — Binning pipeline"
echo "========================================"
if [[ ! -f output/config/auto_binning.txt ]]; then
    echo ""
    echo " SKIPPED: output/config/auto_binning.txt does not exist yet."
    echo " This is expected on a dry run (-n) — stages 1 and 2 must"
    echo " run first to generate it. Re-run without -n to execute the"
    echo " full pipeline, or run the remaining stages manually:"
    echo ""
    echo "   snakemake --snakefile Snakefile-bin [options] rename_mags \\"
    echo "     --config binning=output/config/auto_binning.txt"
    echo ""
    echo "   snakemake --snakefile Snakefile-bin [options] virus_all \\"
    echo "     --config binning=output/config/auto_binning.txt"
    echo ""
else
    # Target before "$@"; keep the variable-length --config last so it
    # doesn't consume, and isn't consumed by, options in the user's args.
    snakemake --snakefile Snakefile-bin rename_mags "$@" \
        --config binning=output/config/auto_binning.txt

    if [[ "$RUN_VIRUS" -eq 1 ]]; then
        echo ""
        echo "========================================"
        echo " MAGmaker  Step 4/$TOTAL — Viral track (geNomad + CheckV)"
        echo "                          (disable with --skip-virus)"
        echo "========================================"
        # virus_all lives in Snakefile-bin, so it needs the same binning
        # override to parse; the viral rules themselves are per-sample and do
        # not use the binning groups.
        snakemake --snakefile Snakefile-bin virus_all "$@" \
            --config binning=output/config/auto_binning.txt
    fi

    echo ""
    echo "========================================"
    echo " MAGmaker  Complete"
    echo " Output:  output/mag_qc/renamed_mags/"
    echo " Summary: output/mag_qc/mag_summary.tsv"
    if [[ "$RUN_VIRUS" -eq 1 ]]; then
        echo " Viral:   output/virus/virus_summary.tsv"
    fi
    echo "========================================"
fi
