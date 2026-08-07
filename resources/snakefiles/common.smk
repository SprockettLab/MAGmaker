# =============================================================================
# common.smk — sample metadata loading, shared by Snakefile and Snakefile-bin
#
# Both entry points need the same metadata, so the loader lives here rather
# than being duplicated. This file defines no rules; it only exposes:
#
#   metadata_table   pandas DataFrame indexed by (Sample, Sequencing_Run)
#   samples          unique sample names, in file order
#   seqruns          the (Sample, Sequencing_Run) MultiIndex
#   reads            read identifiers from config, normally ['R1', 'R2']
#   get_read()       (sample, seqrun, read) -> FASTQ path
#   seqruns_for()    sample -> list of its sequencing runs
#
# Columns beyond the four required ones are carried through untouched and can
# be read from metadata_table inside a rule, so study covariates live
# alongside the file paths instead of in a separate sheet.
# =============================================================================

import sys
import pandas as pd

REQUIRED_COLUMNS = ["Sample", "Sequencing_Run", "R1_fp", "R2_fp"]


def _fail(message):
    """Abort with a message readable without a Python traceback."""
    sys.exit("\nMAGmaker metadata error\n" + "-" * 60 + "\n" + message + "\n")


def load_metadata(config):
    """in : the snakemake config dict
       out: (metadata_table, samples)

    Validation happens up front so a malformed sheet fails immediately with a
    specific message, instead of surfacing later as a missing-input or
    KeyError somewhere deep in the DAG."""

    if "metadata" not in config:
        # Point at the migration rather than just reporting a missing key.
        if "samples" in config or "units" in config:
            _fail(
                "config.yaml still uses the old 'samples:' and 'units:' keys.\n"
                "They are replaced by one metadata table:\n\n"
                "    metadata: resources/config/metadata.txt\n\n"
                "The metadata file is units.txt with renamed columns:\n"
                "    Sample  Sequencing_Run  R1_fp  R2_fp  [covariates...]\n\n"
                "samples.txt is no longer needed; the sample list comes from\n"
                "the Sample column. To convert an existing units.txt:\n\n"
                "    sed '1s/.*/Sample\\tSequencing_Run\\tR1_fp\\tR2_fp/' \\\n"
                "        units.txt > metadata.txt\n\n"
                "See docs/configuration.md."
            )
        _fail("config.yaml must define 'metadata' (path to the metadata table).")

    metadata_fp = config["metadata"]

    try:
        table = pd.read_csv(metadata_fp, sep="\t", header=0, dtype=str)
    except FileNotFoundError:
        _fail("metadata file not found: %s" % metadata_fp)
    except Exception as exc:
        _fail("could not read metadata file %s\n  %s" % (metadata_fp, exc))

    # Tabs are easily lost when a sheet is pasted through a terminal or
    # reformatted by an editor; a single-column frame is the signature.
    if table.shape[1] == 1:
        _fail(
            "metadata file %s parsed as a single column.\n"
            "It must be TAB-separated. Check with:  cat -A %s | head -2\n"
            "(real tabs appear as ^I)" % (metadata_fp, metadata_fp)
        )

    table.columns = [str(c).strip() for c in table.columns]

    missing = [c for c in REQUIRED_COLUMNS if c not in table.columns]
    if missing:
        _fail(
            "metadata file %s is missing required column(s): %s\n"
            "Required: %s\n"
            "Found   : %s"
            % (metadata_fp, ", ".join(missing),
               ", ".join(REQUIRED_COLUMNS), ", ".join(table.columns))
        )

    # fillna before astype: in pandas >=3 a missing value survives .astype(str)
    # as a float NaN rather than the string "nan", so comparing to "" alone
    # would silently miss empty fields.
    for column in REQUIRED_COLUMNS:
        table[column] = table[column].fillna("").astype(str).str.strip()

    blank = table[
        table[REQUIRED_COLUMNS].isin(["", "nan", "NaN", "None", "NA"]).any(axis=1)
    ]
    if not blank.empty:
        # +2 converts a 0-based row index to the line number in the file.
        lines = ", ".join(str(i + 2) for i in blank.index)
        _fail(
            "metadata file %s has empty required fields on line(s): %s\n"
            "Every row needs Sample, Sequencing_Run, R1_fp and R2_fp."
            % (metadata_fp, lines)
        )

    # A repeated (Sample, Sequencing_Run) makes the index ambiguous; caught
    # here because pandas would otherwise raise with no hint at the cause.
    duplicated = table.duplicated(subset=["Sample", "Sequencing_Run"], keep=False)
    if duplicated.any():
        pairs = (table.loc[duplicated, ["Sample", "Sequencing_Run"]]
                 .drop_duplicates().itertuples(index=False, name=None))
        _fail(
            "metadata file %s has duplicate (Sample, Sequencing_Run) pairs:\n  %s\n"
            "Each row must describe one sequencing run of one sample."
            % (metadata_fp, "\n  ".join("%s / %s" % p for p in pairs))
        )

    metadata_table = table.set_index(["Sample", "Sequencing_Run"])
    metadata_table.sort_index(inplace=True)

    # Order-preserving unique, so run order follows the file rather than sort
    # order -- keeps logs and reports in the order the user wrote them.
    samples = list(dict.fromkeys(table["Sample"]))

    return metadata_table, samples


metadata_table, samples = load_metadata(config)
seqruns = metadata_table.index
reads = config["reads"]


def get_read(sample, seqrun, read):
    """in : sample name, sequencing run name, read id ('R1' or 'R2')
       out: path to that FASTQ

    `read` is the bare identifier so it can double as a filename wildcard;
    the metadata column it maps to is <read>_fp."""
    return metadata_table.loc[(sample, seqrun), "%s_fp" % read]


def mem_escalate(key, base_default=8000, cap_multiple=4):
    """in : a config['mem_mb'] key, plus fallbacks
       out: a snakemake resource callable of (wildcards, attempt)

    Returns config['mem_mb'][key] on the first attempt and doubles it on each
    retry, capped by config['mem_mb'][key + '_max'].

    mem_mb is a scheduler RESERVATION, not a ceiling: a job requesting 300 GB
    will only run on a node with 300 GB free, however much it actually uses.
    Sizing every request for the worst-observed case therefore confines the
    whole rule to the largest nodes and serialises it, to benefit the few
    samples that need the headroom. Escalating instead keeps the common case
    schedulable anywhere and pays the queueing cost only where it is earned.

    Rules using this need a `retries:` directive, or `attempt` never exceeds 1
    and the escalation is inert."""
    def _resolve(wildcards, attempt):
        base = config["mem_mb"].get(key, base_default)
        cap = config["mem_mb"].get("%s_max" % key, base * cap_multiple)
        return min(cap, base * attempt)
    return _resolve


def seqruns_for(sample):
    """in : sample name
       out: list of that sample's sequencing run names, in index order"""
    return list(metadata_table.loc[sample].index)
