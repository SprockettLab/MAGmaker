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
#   is_single_end()  sample -> True for a single-end sample
#   reads_for()      sample -> ['R1', 'R2'] or ['SE']
#   nonhost_reads()  sample -> host-filtered FASTQ path(s)
#   samples_for_assembler()  assembler -> samples it can assemble
#
# Columns beyond the four required ones are carried through untouched and can
# be read from metadata_table inside a rule, so study covariates live
# alongside the file paths instead of in a separate sheet.
#
# SINGLE-END INPUT
#
# A row declares single-end sequencing by putting the literal NA in R2_fp.
# Nothing is inferred: not from the filesystem, not from read headers, not
# from a missing column. An empty, absent, "nan" or "None" R2_fp stays a
# hard error exactly as before, so an incomplete or half-written sample
# sheet still fails loudly instead of quietly assembling single-end.
#
# Single-end samples carry the read identifier 'SE' instead of 'R1'/'R2'.
# Because every QC rule is already parameterized by read identifier, that
# one substitution carries most of the way through the pipeline, and it
# keeps single-end outputs from ever colliding with paired-end ones.
# =============================================================================

import sys
import pandas as pd

REQUIRED_COLUMNS = ["Sample", "Sequencing_Run", "R1_fp", "R2_fp"]

# The literal that declares single-end input. Case-sensitive and exact:
# "na", "n/a" and "None" are all still errors.
SINGLE_END_MARKER = "NA"

# Read identifier used for single-end samples throughout the pipeline.
SINGLE_END_READ = "SE"

# Values that mean "this field was left blank", and are always an error.
BLANK_VALUES = ["", "nan", "NaN", "None"]


# Sample and sequencing-run identifiers never contain a path separator, but
# Snakemake wildcards match across "/" unless they are constrained. Without
# this, output/qc/fastp/se/<sample>.<run>.fastp.json is matchable by fastp_pe
# by binding sample="se/<sample>", so the paired rule is selected for a
# single-end report and dies in its input function with a KeyError naming a
# sample that was never in the metadata. Constraining the two wildcards the
# se/ subdirectory can be absorbed into sends each report to its own rule.
wildcard_constraints:
    sample=r"[^/]+",
    seqrun=r"[^/]+"


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
        # keep_default_na=False is load-bearing: pandas otherwise converts the
        # literal NA to NaN, which is exactly the value used to declare
        # single-end input, and it would arrive here indistinguishable from a
        # genuinely empty cell. It also stops NA in any carried-through
        # covariate column from being silently rewritten.
        table = pd.read_csv(metadata_fp, sep="\t", header=0, dtype=str,
                            keep_default_na=False)
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

    # R2_fp is checked separately below, because the literal NA is meaningful
    # there and an error everywhere else.
    always_required = [c for c in REQUIRED_COLUMNS if c != "R2_fp"]

    blank = table[
        table[always_required].isin(BLANK_VALUES + [SINGLE_END_MARKER]).any(axis=1)
    ]
    if not blank.empty:
        # +2 converts a 0-based row index to the line number in the file.
        lines = ", ".join(str(i + 2) for i in blank.index)
        _fail(
            "metadata file %s has empty required fields on line(s): %s\n"
            "Every row needs Sample, Sequencing_Run and R1_fp."
            % (metadata_fp, lines)
        )

    # A blank R2_fp is still an error. Only the exact literal NA declares
    # single-end, so a truncated sheet cannot silently downgrade a paired
    # library to single-end and assemble half the data.
    blank_r2 = table[table["R2_fp"].isin(BLANK_VALUES)]
    if not blank_r2.empty:
        lines = ", ".join(str(i + 2) for i in blank_r2.index)
        _fail(
            "metadata file %s has a blank R2_fp on line(s): %s\n\n"
            "R2_fp is never optional. Give the path to the reverse reads,\n"
            "or write the literal\n\n"
            "    NA\n\n"
            "to declare that the run is single-end. Declaring it explicitly\n"
            "is required so that an incomplete sample sheet cannot be\n"
            "mistaken for single-end data." % (metadata_fp, lines)
        )

    table["Layout"] = table["R2_fp"].apply(
        lambda v: "SINGLE" if v == SINGLE_END_MARKER else "PAIRED"
    )

    # merge_seqruns concatenates per read identifier, so a sample whose runs
    # disagree about layout would produce a reverse file covering only some
    # of its runs. Caught here rather than as a short R2 much later.
    mixed = (table.groupby("Sample")["Layout"].nunique() > 1)
    if mixed.any():
        names = ", ".join(sorted(mixed[mixed].index))
        _fail(
            "metadata file %s mixes single-end and paired-end runs within\n"
            "the same sample: %s\n\n"
            "Every sequencing run of one sample must have the same layout."
            % (metadata_fp, names)
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

# sample -> "SINGLE" or "PAIRED". Validated in load_metadata to be constant
# within a sample, so the first run of each is representative.
layout_index = (metadata_table.reset_index()
                .groupby("Sample")["Layout"].first().to_dict())


def get_read(sample, seqrun, read):
    """in : sample name, sequencing run name, read id ('R1', 'R2' or 'SE')
       out: path to that FASTQ

    `read` is the bare identifier so it can double as a filename wildcard;
    the metadata column it maps to is <read>_fp. 'SE' has no column of its
    own -- a single-end run keeps its only FASTQ in R1_fp -- so it reads
    R1_fp while carrying a distinct identifier through output filenames."""
    column = "R1_fp" if read == SINGLE_END_READ else "%s_fp" % read
    return metadata_table.loc[(sample, seqrun), column]


def is_single_end(sample):
    """in : sample name
       out: True if the sample was sequenced single-end

    Layout is validated to be consistent across a sample's runs, so the
    first run answers for all of them."""
    return layout_index[sample] == "SINGLE"


def reads_for(sample):
    """in : sample name
       out: ['SE'] for single-end, otherwise the configured read ids

    The list every per-read expansion should iterate, in place of the
    global `reads`, so single-end samples never have an R2 target."""
    return [SINGLE_END_READ] if is_single_end(sample) else reads


def nonhost_reads(sample):
    """in : sample name
       out: list of host-filtered FASTQ paths, one entry or two

    Rules take this instead of referring to host_filter's outputs directly.
    The paired and single-end host filters are separate rules with separate
    outputs, so a fixed `rules.host_filter.output.nonhost_R2` cannot serve
    both."""
    return expand("output/qc/host_filter/nonhost/{sample}.{read}.fastq.gz",
                  sample=sample, read=reads_for(sample))


def samples_for_assembler(assembler):
    """in : assembler name
       out: the samples that assembler can run on

    metaSPAdes rejects a single-end-only library outright, so single-end
    samples are assembled by MEGAHIT alone. Filtering the target list keeps
    that a planning decision rather than a run-time crash partway through a
    long job."""
    if assembler == "metaspades":
        return [s for s in samples if not is_single_end(s)]
    return list(samples)


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
