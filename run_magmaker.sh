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
#   ./run_magmaker.sh --profile resources/profiles/demon \
#       --configfile resources/config/<Project>_config.yaml   # project config
#   ./run_magmaker.sh --profile resources/profiles/demon \
#       --skip-kraken --skip-metaphlan          # no read-level profiling
#   ./run_magmaker.sh --profile resources/profiles/demon \
#       --skip-cmseq                            # no strain heterogeneity
#   ./run_magmaker.sh --profile resources/profiles/demon \
#       --binette                               # Binette instead of DAS_Tool
#   ./run_magmaker.sh --profile resources/profiles/demon \
#       --skip-gtdbtk-classify                  # MSA_Percent, no pplacer
#
# --binette / --das-tool choose which tool reconciles the three per-binner
# bin sets into one set of MAGs, overriding params.consolidation.tool for
# this run only. DAS_Tool scores candidates with 51 single-copy genes;
# Binette scores them with CheckM2, the same estimator mag_qc then uses to
# judge the result. Their outputs land in DAS_Tool_Fastas/ and
# Binette_Fastas/ respectively, so both can exist side by side in one run
# directory and be compared.
#
# --skip-kraken and --skip-metaphlan drop the corresponding read-level
# profiling targets for this run only, leaving the config file alone. Use
# them when taxonomy comes from GTDB-Tk on the MAGs rather than from reads;
# kraken2 in particular reads a ~200 GB index per sample and is the
# heaviest I/O in the pipeline for what it contributes.
#
# To stop after the MAG summary table so you can review/edit
# output/mag_qc/mag_summary.tsv before renaming, run stages manually:
#   snakemake [options]
#   snakemake [options] generate_binning_config
#   snakemake --snakefile Snakefile-bin [options] make_mag_summary \
#     --config binning=output/config/auto_binning.txt

set -euo pipefail

# --skip-kraken / --skip-metaphlan are consumed here and turned into config
# overrides; every other argument is passed through untouched. They are
# separate flags because the two tools are useful independently: an
# analysis taking its taxonomy from GTDB-Tk on the MAGs needs neither, but
# one that wants community profiles may still want metaphlan without
# paying for kraken2's ~200 GB index read per sample.
SKIP_CONFIG=()
PASSTHRU=()
for arg in "$@"; do
    case "${arg}" in
        --skip-kraken)    SKIP_CONFIG+=("skip_kraken=True") ;;
        --skip-metaphlan) SKIP_CONFIG+=("skip_metaphlan=True") ;;
        --skip-cmseq)     SKIP_CONFIG+=("skip_cmseq=True") ;;
        # --------------------------------------------------------------
        # Everything collected here becomes ONE --config on each snakemake
        # command line. Snakemake's --config takes nargs="*" with argparse's
        # default store action, so a second --config REPLACES the first
        # rather than adding to it. Step 3 used to append its own
        # `--config binning=...` after these, which silently discarded every
        # override above it -- so --skip-cmseq still ran CMSeq,
        # --skip-gtdbtk-classify still ran pplacer, and --binette still ran
        # DAS_Tool, all without a word of warning. Steps 1 and 2 were
        # unaffected, which is why --skip-kraken appeared to work: kraken
        # only runs in step 1.
        # --------------------------------------------------------------
        # GTDB-Tk runs identify, align and classify. Nearly all of the cost
        # is classify, which runs pplacer and is why this rule reserves
        # 128-320 GB. --skip-gtdbtk-classify stops after align, which still
        # produces MSA_Percent -- the share of the marker alignment a genome
        # fills, and what decides whether it can be placed in a tree -- but
        # leaves the taxonomy columns NA. --skip-gtdbtk drops both.
        --skip-gtdbtk)         SKIP_CONFIG+=("skip_gtdbtk=True") ;;
        --skip-gtdbtk-classify) SKIP_CONFIG+=("skip_gtdbtk_classify=True") ;;
        # Which tool reconciles the per-binner bin sets. Given here rather
        # than only in the config file so one arm can be run both ways
        # without editing a tracked file between the two.
        --binette)        SKIP_CONFIG+=("consolidation_tool=binette") ;;
        --das-tool)       SKIP_CONFIG+=("consolidation_tool=das_tool") ;;
        *)                PASSTHRU+=("${arg}") ;;
    esac
done
if [[ ${#SKIP_CONFIG[@]} -gt 0 ]]; then
    echo "config overrides: ${SKIP_CONFIG[*]}"
fi

# A --config of the user's own would be a second one on the final command
# line and would silently drop everything this script sets, which is the
# failure described above. Refused rather than merged: merging would mean
# parsing arbitrary KEY=VALUE lists, and the flags exist precisely so that
# is not necessary.
for arg in ${PASSTHRU[@]+"${PASSTHRU[@]}"}; do
    if [[ "${arg}" == "--config" || "${arg}" == "-C" ]]; then
        echo "ERROR: pass --config through run_magmaker.sh is not supported." >&2
        echo "Snakemake keeps only the LAST --config on a command line, so" >&2
        echo "yours would silently discard the ones this script sets for" >&2
        echo "--skip-* and --binette/--das-tool, or be discarded by them." >&2
        echo "" >&2
        echo "Use the dedicated flags, or call snakemake directly:" >&2
        echo "  snakemake --snakefile Snakefile-bin rename_mags [options] \\" >&2
        echo "    --config binning=output/config/auto_binning.txt KEY=VALUE" >&2
        exit 1
    fi
done

# Steps 1 and 2 take the overrides as given. Step 3 needs `binning=` in the
# SAME --config, which is why it is built separately rather than appended
# to a shared one at the call site.
STAGE12=(${PASSTHRU[@]+"${PASSTHRU[@]}"})
if [[ ${#SKIP_CONFIG[@]} -gt 0 ]]; then
    STAGE12+=("--config" "${SKIP_CONFIG[@]}")
fi
STAGE3=(${PASSTHRU[@]+"${PASSTHRU[@]}"} "--config"
        ${SKIP_CONFIG[@]+"${SKIP_CONFIG[@]}"}
        "binning=output/config/auto_binning.txt")

if [[ ! -f Snakefile || ! -f Snakefile-bin ]]; then
    echo "ERROR: run_magmaker.sh must be run from the MAGmaker repository root." >&2
    exit 1
fi

echo "========================================"
echo " MAGmaker  Step 1/3 — Main pipeline"
echo "========================================"
snakemake ${STAGE12[@]+"${STAGE12[@]}"}

echo ""
echo "========================================"
echo " MAGmaker  Step 2/3 — Binning config"
echo "========================================"
# Target goes before "$@" so a trailing variable-length option in the
# user's args (e.g. --configfile FILE, which accepts one or more values)
# cannot swallow the target name.
snakemake generate_binning_config ${STAGE12[@]+"${STAGE12[@]}"}

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
    # Target before "$@"; keep the variable-length --config last so it
    # doesn't consume, and isn't consumed by, options in the user's args.
    snakemake --snakefile Snakefile-bin rename_mags ${STAGE3[@]+"${STAGE3[@]}"}

    echo ""
    echo "========================================"
    echo " MAGmaker  Complete"
    echo " Output:  output/mag_qc/renamed_mags/"
    echo " Summary: output/mag_qc/mag_summary.tsv"
    echo "========================================"
fi
